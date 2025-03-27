library(tidyverse)
here::i_am("scripts/update_metadata.R")
raw_data_dir <- "gtex_v10_shiny/data/raw_data/"
metadata_path <- file.path(raw_data_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
if (!file.exists(metadata_path)) stop("Error: Metadata file missing at ", metadata_path)

metadata <- read.table(metadata_path, sep = "\t", header = TRUE)
metadata <- metadata %>%
  rename(donor = SUBJID, sex = SEX, age = AGE, death_type = DTHHRDY) %>%
  mutate(age_plot = as.numeric(sub("-.*", "", age)), 
         sex_plot = ifelse(sex == 1, "Male", "Female"))
metadata$sex <- factor(metadata$sex)
save_path <- file.path(raw_data_dir, "Updated_GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
write.table(metadata, save_path, sep = "\t", row.names = FALSE)
