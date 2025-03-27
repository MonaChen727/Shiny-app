library(shiny)
library(tidyverse)
library(ggplot2)

# Set the app directory
if (requireNamespace("rprojroot", quietly = TRUE)) {
  app_dir <- rprojroot::find_rstudio_root_file()
  print(app_dir)
}

# Define paths relative to the app directory
data_dir <- file.path(app_dir, "gtex_v10_shiny/data")  # Go one level up to access gtex_v10_shiny
raw_data_dir <- file.path(data_dir, "raw_data")

# Load tissue names
tissue_file <- file.path(data_dir, "tissue_names.txt")
if (file.exists(tissue_file)) {
  tissues <- readLines(tissue_file)
  tissues <- gsub("_", " ", tissues)
} else {
  stop("Tissue names file not found!")
}

# Load gene names
gene_file <- file.path(data_dir, "gene_names.txt")
if (file.exists(gene_file)) {
  genes <- readLines(gene_file)
} else {
  stop("Gene names file not found!")
}

# Define UI
ui <- fluidPage(
  titlePanel("GTEx Gene Expression Analysis"),
  sidebarLayout(
    sidebarPanel(
      selectInput("gene", "Select gene:", 
                  choices = genes),
      selectInput("tissue", "Select tissue:", 
                  choices = gsub("_", " ", tissues)),
      actionButton("plot", "Generate Plot"),
      checkboxInput("log2_transform", "Apply log2 transformation to TPM", value = TRUE),  # Checkbox to toggle log transformation
      checkboxInput("log10_transform", "Apply log10 transformation to TPM", value = TRUE)  # Checkbox to toggle log transformation
      ),
    mainPanel(
      plotOutput("tpmPlot"),      
      textOutput("errorMsg")  # Error message output

    )
  )
)

# Define Server
server <- function(input, output, session) {
  # Function to read and preprocess data
  read_and_preprocess_data <- function(gene, tissue) {
    req(gene, tissue)  
    exp.path <- file.path(raw_data_dir, sprintf("gene_tpm_v10_%s.gct.gz", gsub(" ", "_", tissue)))
    metadata.path <- file.path(raw_data_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
    
    if (!file.exists(exp.path)) return(paste("Error: Expression data not found at", exp.path))
    if (!file.exists(metadata.path)) return(paste("Error: Metadata file missing at", metadata.path))
    
    exp <- read.table(gzfile(exp.path), sep = "\t", skip = 2, header = TRUE)
    metadata <- read.table(metadata.path, sep = "\t", header = TRUE)
    
    colnames(exp) <- gsub("\\.", "-", colnames(exp))
    metadata <- metadata %>%
      rename(donor = SUBJID, sex = SEX, age = AGE, death_type = DTHHRDY) %>%
      mutate(
        age_plot = as.numeric(sub("-.*", "", age)),  
        sex_plot = ifelse(sex == 1, "Male", "Female")
      )
    
    # datasets based on tissue and gene
    X <- exp %>% filter(Description == gene)
    if (nrow(X) == 0) return("Error: Gene not found.")
    X <- X %>%
      select(-Name, -Description) %>%
      pivot_longer(cols = everything(), names_to = "sample", values_to = "TPM") %>%
      mutate(donor = str_extract(sample, "^[^-]+-[^-]+"))
    mergedData <- left_join(X, metadata, by = "donor")
    if (nrow(mergedData) == 0) return("Error: No matching data for this gene-tissue pair.")
    
    return(mergedData)
  }
  
  # Reactive expression to generate plot data
  plotData <- eventReactive(input$plot, {
    res <- read_and_preprocess_data(input$gene, input$tissue)
    if (is.character(res)) return(data.frame(age_plot = NA, TPM = NA, sex_plot = NA))  # Return empty dataframe
    res
  })

  
  # Render plot
  output$tpmPlot <- renderPlot({
    df <- plotData()
    validate(
      need(!is.null(df), "No data available for the selected gene and tissue.")
    )

    # Apply log transformation if checkbox is selected
    if (input$log2_transform) {
      df$TPM <- log2(df$TPM)
    }
    if (input$log10_transform) {
      df$TPM <- log10(df$TPM)
    }
    
    
    ggplot(data = df, aes(x = age_plot, y = (TPM), colour = sex_plot)) +
      geom_smooth(method = "lm", formula = y ~ x, fill = "lightgray", alpha = 0.3) +
      geom_point(alpha = 0.7, size = 2) +
      scale_color_manual(name = "Sex", values = c("Male" = "steelblue", "Female" = "red")) +
      ggtitle(sprintf(" Expression of %s in %s", input$gene, input$tissue)) +
      xlab("Age") +
      ylab("log2_TPM") + 
      theme_minimal()
  })
  
  
  # Error message output
  output$errorMsg <- renderText({
    df <- plotData()
    if (is.null(df) || nrow(df) == 0) {
      return("Error: Unable to load data for the selected gene and tissue.")
    }
    NULL
  })
}

# Run the application
shinyApp(ui = ui, server = server)

