library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(optparse)

option_list = list(
  make_option("--workdir", action= "store", type= 'character', default=NA),
  make_option("--coverage", action= "store", type= 'character', default=NA)
  )

p <- OptionParser(option_list=option_list)

args <- parse_args(p)

thr <- args$coverage

mosdepth_ls <- list.files(args$workdir,
                   recursive = TRUE,
                   full.names = TRUE,
                   pattern = "regions.bed.gz$") 

df_list <- map(mosdepth_ls, fread)
names(df_list) <- gsub("\\..*", "", basename(mosdepth_ls))

df_list <- df_list %>%
  bind_rows(.id = "sampleID")

gene_pass <- df_list %>%
  group_by(V4) %>%
  summarise(midcov = median(V5)) %>%
  filter(midcov >= thr) %>%
  select(V4) %>%
  distinct() %>%
  pull()

passed <- signif(length(gene_pass)/length(unique(df_list$V4)), 4) * 100

passed <- paste(passed, "% exons with > ", thr, "x coverage", sep = "")
path_save <- paste("saving to: regions_passed_",
            thr,
            "x.bed", sep = "")

fileConn <- file("mosdepth_filter.log")
writeLines(c(passed,
             path_save), 
             fileConn)
close(fileConn)

df_list %>%
  bind_rows(.id = "sampleID") %>%
  filter(V4 %in% gene_pass) %>%
  select(V1:V3) %>%
  distinct() %>%
  write.table(paste("regions_passed_",
                    thr,
                    "x.bed", sep = ""),
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
  
