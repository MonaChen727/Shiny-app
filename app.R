library(tidyverse)  # tidyverse includes dplyr, ggplot2, etc.
library(DT)         # For the datatable
library(conflicted)
conflicts_prefer(dplyr::filter)

# Set the app directory
# if (requireNamespace("rprojroot", quietly = TRUE)) {
#   app_dir <- rprojroot::find_rstudio_root_file()
#   print(app_dir)
# }
app_dir <- here::here()

# Define paths relative to the app directory
data_dir <- file.path(app_dir, "gtex_v10_shiny/data")
raw_data_dir <- file.path(data_dir, "raw_data")

# Load tissue names
tissue_file <- file.path(data_dir, "tissue_names.txt")
tissues <- if (file.exists(tissue_file)) {
  readLines(tissue_file) %>% gsub("_", " ", .)
} else {
  stop("Tissue names file not found!")
}

# Load gene names
gene_file <- file.path(data_dir, "gene_names.txt")
genes <- if (file.exists(gene_file)) {
  readLines(gene_file)
} else {
  stop("Gene names file not found!")
}

# Define UI
ui <- fluidPage(
  titlePanel("GTEx Gene Expression Analysis"),
  sidebarLayout(
    sidebarPanel(
      selectInput("gene", "Select Gene", choices = genes, selectize = TRUE),
      selectInput("tissue", "Select tissue:", choices = tissues),
      actionButton("plot", "Generate Plot"),
      checkboxInput("show_pvalue", "Show P-value Analysis", value = TRUE),  # Add checkbox for p-value
      checkboxInput("plot_type", "Show Violin & Box Plot (Uncheck for Regression Plot)", value = FALSE)
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("tpmPlot")),
        tabPanel("P-Value Table", 
                 DT::dataTableOutput("pvalueTable"),
                 downloadButton("downloadTable", "Download P-Value Table"))
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Function to read and preprocess data
  read_and_preprocess_data <- function(gene, tissue) {
    req(gene, tissue)  # Ensure both inputs are available
    
    # exp.path <- file.path(raw_data_dir, sprintf("gene_tpm_v10_%s.gct.gz", gsub(" ", "_", tissue)))
    # if (!file.exists(exp.path)) return(NULL)
    # exp <- read.table(gzfile(exp.path), sep = "\t", skip = 2, header = TRUE)
    
    # Dynamically generate the URL for the tissue
    tissue_url <- sprintf("https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_%s.gct.gz", gsub(" ", "_", tissue))
    # Try to read the data from the URL
    exp <- tryCatch({
      temp_file <- tempfile(fileext = ".gz")
      download.file(tissue_url, temp_file, mode = "wb", quiet = TRUE)
      read.table(gzfile(temp_file), sep = "\t", skip = 2, header = TRUE)
    }, error = function(e) {
      return(NULL)
    })
    
    # Check the result
    if (is.null(exp)) {
      print("Error: Unable to download or read the file.")
    } 
    
    
    metadata_path <- file.path(raw_data_dir, "Updated_GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
    if (!file.exists(metadata_path)) return(NULL)
    metadata <- read.table(metadata_path, sep = "\t", header = TRUE)
    
    X <- exp %>% filter(Description == gene)
    if (nrow(X) == 0) return(NULL)
    
    X <- X %>%
      select(-Name, -Description) %>%
      pivot_longer(cols = everything(), names_to = "sample", values_to = "TPM") %>%
      mutate(donor = sub("^([^.]+)\\.([^.]+).*", "\\1-\\2", sample)) %>% 
      mutate(logTPM = log2(TPM + 1)) %>% 
      select(donor, logTPM)
    
    mergedData <- left_join(X, metadata, by = "donor")
    if (nrow(mergedData) == 0) return(NULL)
    
    return(mergedData)
  }
  
  # Reactive expression to generate plot data
  plotData <- eventReactive(input$plot, {
    tryCatch({
      df <- read_and_preprocess_data(input$gene, input$tissue)
      validate(need(!is.null(df), "Error: No matching data for this gene-tissue pair."))
      df
    }, error = function(e) {
      return(data.frame(error_msg = e$message))
    })
  })
  
  # Render the plot
  output$tpmPlot <- renderPlot({
    df <- plotData()
    
    # Check for errors in the reactive data
    validate(need(!("error_msg" %in% colnames(df) && !is.null(df$error_msg)), df$error_msg))
    
    # P-value analysis
    if (input$show_pvalue) {
      
      # Check for single level in sex_plot
      if (length(unique(df$sex_plot)) == 1) {
        fit <- lm(logTPM ~ age_plot, data = df)
      } else {
        fit <- lm(logTPM ~ age_plot + sex_plot, data = df)
      }
      
      coefs <- summary(fit)$coefficients
      p_val <- if("age_plot" %in% rownames(coefs)) signif(coefs["age_plot", "Pr(>|t|)"], 4) else NA_real_
      
      df$pred_logTPM <- predict(fit, newdata = df)
      pval_label <- paste0("p = ", p_val)
      
      ggplot(df, aes(x = age_plot, y = logTPM, color = sex_plot)) +
        geom_point(alpha = 0.7, size = 2) +
        geom_line(aes(y = pred_logTPM), linewidth = 1) +
        scale_color_manual(name = "Sex", values = c("Male" = "steelblue", "Female" = "red")) +
        annotate("text", x = min(df$age_plot, na.rm = TRUE) + 2, 
                 y = max(df$logTPM, na.rm = TRUE),
                 label = pval_label, hjust = -4, vjust = -0.5, size = 5) +
        labs(title = sprintf("Expression of %s in %s", input$gene, input$tissue),
             x = "Age", y = "log2TPM") +
        theme_minimal()
    } else {
      # Plot without p-value analysis
      if (input$plot_type) {
        df$age_group <- cut(
          df$age_plot,
          breaks = c(0, 20, 30, 40, 50, 60, 70, Inf), 
          right = FALSE,
          labels = c("0-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70+")
        )
        
        df_summary <- df %>%
          group_by(age_group, sex_plot) %>%
          summarize(mean_logTPM = mean(logTPM, na.rm = TRUE), .groups = "drop")
        
        ggplot(df, aes(x = age_group, y = logTPM, fill = sex_plot, color = sex_plot)) +
          geom_boxplot(
            position = position_dodge(width = 0.8), 
            alpha = 0.4, 
            outlier.shape = NA
          ) +
          geom_line(
            data = df_summary,
            aes(x = age_group, y = mean_logTPM, group = sex_plot),
            linetype = "dashed",
            position = position_dodge(width = 0.8)
          ) +
          geom_point(
            data = df_summary,
            aes(x = age_group, y = mean_logTPM),
            position = position_dodge(width = 0.8),
            size = 2
          ) +
          scale_fill_manual(name = "Sex", values = c("Male" = "steelblue", "Female" = "red")) +
          scale_color_manual(name = "Sex", values = c("Male" = "steelblue", "Female" = "red")) +
          theme_bw() +
          labs(
            title = sprintf("Expression of %s in %s (Box Plot by Age Group)", input$gene, input$tissue),
            x = "Age Range",
            y = "log2(TPM)"
          )
      } else {
        ggplot(df, aes(x = age_plot, y = logTPM, color = sex_plot)) +
          geom_smooth(method = "lm", formula = y ~ x, fill = "lightgray", alpha = 0.3) +
          geom_point(alpha = 0.7, size = 2) +
          scale_color_manual(name = "Sex", values = c("Male" = "steelblue", "Female" = "red")) +
          ggtitle(sprintf("Expression of %s in %s", input$gene, input$tissue)) +
          xlab("Age") + ylab("log2(TPM)") + theme_minimal()
      }
    }
  })
  
  # Function to load the p-value table
  load_pvalue_table <- function(tissue) {
    output_path <- file.path(here::here("gtex_v10_shiny/data/p_value"), 
                             paste0(gsub(" ", "_", tissue), "_pvalue_results.csv"))
    if (!file.exists(output_path)) {
      return(NULL)
    }
    df <- read.csv(output_path)
    
    if ("p_value" %in% colnames(df) && length(df$p_value) > 0) {
      df$p_value_rounded <- signif(df$p_value, 4)
      df$adj_p_value <- signif(p.adjust(df$p_value, method = "BH"), 4)
      df <- df[, setdiff(colnames(df), "p_value")]
    } else {
      warning("p_value does not exist")
      df$p_value_rounded <- NA
      df$adj_p_value <- NA
    }
    
    return(df)
  }
  
  # Render the p-value table
  output$pvalueTable <- DT::renderDataTable({
    req(input$tissue)
    pvalue_data <- load_pvalue_table(input$tissue)
    
    validate(
      need(!is.null(pvalue_data), "Error: No matching p-value data for this tissue.")
    )
    
    DT::datatable(pvalue_data, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  # Download handler
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste0(input$tissue, "_pvalue_results.csv")
    },
    content = function(file) {
      pvalue_data <- load_pvalue_table(input$tissue)
      if (is.null(pvalue_data)) {
        stop("P-Value table not found.")
      }
      write.csv(pvalue_data, file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
