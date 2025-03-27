# obtain the tissues name
tissue_names <- sub("^.*/gene_tpm_v10_(.*?)\\.gct\\.gz$", "\\1", files)
tissue_names <- gsub("_", " ", tissue_names)
writeLines(tissue_names, here::here("gtex_v10_shiny/data/tissue_names.txt"))

# obtain the genes name
exp.path <- sprintf(here::here("gtex_v10_shiny/data/raw_data/gene_tpm_v10_%s.gct.gz"), "uterus")
exp <- read.table(gzfile(exp.path), sep = "\t", skip = 2, header = TRUE)
# Extract gene names from the second column
gene_names <- sort(unique(exp[, 2]))
writeLines(gene_names, here::here("gtex_v10_shiny/data/gene_names.txt"))

