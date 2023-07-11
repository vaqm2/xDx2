#!/usr/bin/env Rscript

library(dplyr)
library(gprofiler2)

args = commandArgs(trailingOnly = TRUE)
files = list.files(path = getwd(), pattern = "*_DevPEC.txt")
out_file = paste0(args[1], "_Enrichments.txt")
file.create(out_file)

write(paste("SOURCE", "TERM", "TERM_SIZE", "INTERSECTION_SIZE", "P_FDR", "TEST", 
            sep = "\t"), 
      out_file,
      append = TRUE)

for (file in files) {
    if(!(file.size(file) == 0)) {
        test = gsub(".txt$", "", file)
        print(paste0("Processing", " ", file, "..."))
        genes = read.table(file, header = F)
        colnames(genes) = ("GENE")
        gost_out = gost(query = genes$GENE,
                        organism = "hsapiens", 
                        ordered_query = TRUE, 
                        significant = TRUE,
                        user_threshold = 0.05, 
                        correction_method = "fdr",
                        evcodes = FALSE,
                        domain_scope = "annotated",
                        sources = c("GO"))
        if(!is.null(gost_out$result)) {
            gost_result = gost_out$result %>% 
                arrange(p_value) %>%
                as.data.frame() %>%
                filter(term_size >= 10 & term_size <= 500 & intersection_size >= 5) %>%
                select(source, term_name, term_size, intersection_size, p_value) %>% 
                mutate(query = test)
    
            write.table(gost_result,
                        out_file, 
                        sep = "\t", 
                        row.names = F,
                        col.names = F,
                        quote = F,
                        append = TRUE)
        }
    print("Done!")
    }
}