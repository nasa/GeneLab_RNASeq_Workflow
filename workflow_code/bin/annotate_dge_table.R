#!/usr/bin/env Rscript

# Annotate Differential Gene Expression Table
# 
# This script adds/updates annotation columns in an existing DGE table for the dge_table entry point.
# It downloads annotation data from a URL or reads from a local file, then merges the annotations
# with the DGE table, overwriting existing annotation columns if they exist.
#
# Usage: Rscript annotate_dge_table.R <dge_table_path> <annotation_file_path> <gene_id_type> <output_path>
#
# Arguments:
#   dge_table_path: Path to input DGE table CSV file
#   annotation_file_path: URL or local path to annotation file (TSV format)
#   gene_id_type: Column name for gene ID matching (e.g., "ENSEMBL", "TAIR")
#   output_path: Path for output annotated DGE table CSV file
#
# The script handles download failures by proceeding with the original DGE table.

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
dge_table_path <- args[1]
annotation_file_path <- args[2]
gene_id_type <- args[3]
output_path <- args[4]

cat("=== DGE_TABLE ENTRY POINT: ANNOTATE DGE TABLE ===\n")
cat("Input DGE table:", dge_table_path, "\n")

if (!file.exists(dge_table_path)) {
    stop(paste("DGE table file not found at:", dge_table_path))
}

dge_table <- read.csv(dge_table_path, header = TRUE, check.names = FALSE)
cat("Loaded DGE table with", nrow(dge_table), "genes and", ncol(dge_table), "columns\n")

original_columns <- colnames(dge_table)
cat("Original columns:", paste(head(original_columns, 10), collapse=", "))
if (length(original_columns) > 10) cat(" ... (", length(original_columns), "total)")
cat("\n\n")

if (!is.null(annotation_file_path) && annotation_file_path != "" && annotation_file_path != "null") {
    cat("=== LOADING ANNOTATIONS ===\n")
    cat("Annotation source:", annotation_file_path, "\n")
    cat("Gene ID column:", gene_id_type, "\n")
    
    tryCatch({
        if (startsWith(annotation_file_path, "http")) {
            local_file <- "gene_annotations.tsv"
            cat("Downloading annotation file...\n")
            download.file(url = annotation_file_path, 
                         destfile = local_file, 
                         method = "libcurl", 
                         headers = c("User-Agent" = "Mozilla/5.0"))
            annot <- read.table(local_file, sep = "\t", header = TRUE, quote = "", comment.char = "")
        } else {
            cat("Reading local annotation file...\n")
            annot <- read.table(annotation_file_path, sep = "\t", header = TRUE, quote = "", comment.char = "")
        }
        
        cat("Loaded", nrow(annot), "annotations with", ncol(annot), "columns\n")
        cat("Annotation columns:", paste(colnames(annot), collapse=", "), "\n\n")
        
        cat("=== MERGING ANNOTATIONS ===\n")
        
        annot_cols <- setdiff(colnames(annot), gene_id_type)
        new_columns <- setdiff(annot_cols, original_columns)
        overwritten_columns <- intersect(annot_cols, original_columns)
        
        merged_table <- annot %>%
            merge(dge_table,
                by = gene_id_type,
                all.y = TRUE,
                suffixes = c("", "_original")
            ) %>%
            select(all_of(gene_id_type), everything())
        
        cols_to_drop <- grep("_original$", colnames(merged_table), value = TRUE)
        if (length(cols_to_drop) > 0) {
            merged_table <- merged_table %>% select(-all_of(cols_to_drop))
        }
        
        cat("New annotation columns added:", length(new_columns), "\n")
        if (length(new_columns) > 0) {
            for (col in new_columns) {
                cat("  + ", col, "\n", sep = "")
            }
        }
        
        if (length(overwritten_columns) > 0) {
            cat("\nColumns overwritten by annotations:", length(overwritten_columns), "\n")
            for (col in overwritten_columns) {
                cat("  * ", col, "\n", sep = "")
            }
        }
        
        if (length(annot_cols) > 0) {
            first_annot_col <- annot_cols[1]
            genes_with_annot <- sum(!is.na(merged_table[[first_annot_col]]))
            cat("\nGenes matched with annotations:", genes_with_annot, "/", nrow(merged_table), 
                " (", round(100 * genes_with_annot / nrow(merged_table), 1), "%)\n", sep = "")
        } else {
            cat("\nNo annotation columns to check\n")
        }
        
        dge_table <- merged_table
        
    }, error = function(e) {
        warning(paste("Could not load or merge annotation file:", e$message))
        cat("\nERROR: Annotation merging failed:", e$message, "\n")
        cat("Proceeding with original DGE table (no annotations added)\n")
    })
} else {
    cat("=== NO ANNOTATIONS PROVIDED ===\n")
    cat("No annotation file specified - outputting DGE table as-is\n")
}

cat("\n=== SAVING OUTPUT ===\n")
write.csv(dge_table, file = output_path, row.names = FALSE)
cat("Saved annotated DGE table to:", output_path, "\n")
cat("Final table dimensions:", nrow(dge_table), "genes x", ncol(dge_table), "columns\n")
