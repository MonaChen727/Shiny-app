### This script plots the expression of a gene in a tissue
### using the GTEx v10 data. The transcripts per million (TPM)
### of the gene is plotted against the age of the subject,
### and the data are separated by sex.

library(tidyverse)
here::i_am("gtex_v10_shiny/scripts/gtex_age_sex_plots.R")
# specify gene name and tissue to plot
gene.name <- "VPS18"
tissue <- "brain_cortex"

# specify directory to save plots
plot.dir <- "C:/Users/ronni/Desktop/tpm_plots"
if (!dir.exists(plot.dir)) { dir.create(plot.dir) }

# specify file paths for data
exp.path <- sprintf(here::here("gtex_v10_shiny/data/gene_tpm_v10_%s.gct.gz"), tissue)
metadata.path <- here::here("gtex_v10_shiny/data/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

# ===== actual script starts here =====
# read in and format expression file
exp <- read.table(gzfile(exp.path), sep='\t', skip=2, header=TRUE)
colnames(exp) <- gsub("\\.","-", colnames(exp)) # change "." to "-" separators

# read in and format metadata
metadata <- read.table(metadata.path, sep='\t', header=TRUE)
colnames(metadata) <- c("donor","sex","age","death_type")
metadata$age_plot <- sapply(metadata$age, FUN=function(a) {as.numeric(strsplit(a,'-')[[1]][1])})
metadata$sex_plot <- gsub(2,"Female", gsub(1,"Male", metadata$sex))

# merge expression and metadata
X <- subset(exp, Description==gene.name) %>% subset(., select=-c(Name, Description)) %>% t() %>% as.data.frame()
X$donor <- sapply(rownames(X), FUN=function(s) {paste0(strsplit(s,'-')[[1]][1], '-', strsplit(s,'-')[[1]][2])})
colnames(X) <- c("TPM","donor")
X <- merge(X, metadata, by='donor')

# plot by age and sex
p <- ggplot(data=X, mapping=aes(x=age_plot, y=TPM, colour=sex_plot)) + 
  geom_smooth(method='lm', formula=y~x, fill='lightgray', alpha=0.3) + geom_point(alpha=0.7, size=2) + 
  scale_color_manual(name='Sex', values=c("Male"='steelblue', "Female"='red')) +
  ggtitle(sprintf("%s expression\n%s", gene.name, tissue)) + xlab('Age') + theme_minimal()

ggsave(file.path(plot.dir, sprintf("%s_age_sex.png", gene.name)), plot=p, scale=1, width=6, height=6, dpi=400, bg='white')
