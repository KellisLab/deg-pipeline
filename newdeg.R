#' @return 
#' @export
new_deseq2 <- function(se,
                       pathology,
                       case,
                       control,
                       sample.col,
                       covariates=NULL, independentFiltering=as.logical(Sys.getenv("DESEQ2_INDEPENDENT", "TRUE")=="TRUE"),
                       shrinkage=Sys.getenv("DESEQ2_SHRINKAGE", "ashr"),
                       prefix="DESeq2") {
    pb = calculate_qc_metrics(se_make_pseudobulk(se, sample.col), assay="counts", qc_vars=c("mt", "ribo", "pc", "chrX", "chrY"))
    X = SummarizedExperiment::assays(pb)$counts
    cd = as.data.frame(SummarizedExperiment::colData(pb))
    cd[[pathology]] = relevel(as.factor(cd[[pathology]]), ref=control)
    if (!is.null(covariates)) {
        available_covariates = covariates[covariates %in% colnames(cd)]
        covariates = available_covariates[available_covariates %in% colnames(deg.filter.design(cd[c(available_covariates)]))]
    }
    design = model.matrix(as.formula(paste0("~", c(pathology, covariates), collapse=" + ")), data=cd)
    design = deg.filter.design(design)
    dds = DESeq2::DESeqDataSetFromMatrix(X, colData=cd, design=design)
    out = DESeq2::DESeq(dds)
    df = DESeq2::results(out, list(make.names(paste0(pathology, case))), independentFiltering=independentFiltering)
    rd = as.data.frame(SummarizedExperiment::rowData(se))
    result_genes = rownames(df)
    new_cols = c("log2FC", "FDR", "stat")
    for (col in new_cols) {
        rd[[paste0(prefix, "_", col)]] = NA
    }
    rd[result_genes, paste0(prefix, "_log2FC")] = df$log2FoldChange
    rd[result_genes, paste0(prefix, "_FDR")] = df$padj
    rd[result_genes, paste0(prefix, "_stat")] = df$stat
    if (shrinkage %in% c("apeglm", "ashr")) {
        coef_names = DESeq2::resultsNames(out)
        target_coef = make.names(paste0(pathology, case))
        IDX = which(grepl(target_coef, coef_names))[1]
        if (is.na(IDX)) IDX = 1
        dfs = DESeq2::lfcShrink(out, res=df, coef=coef_names[IDX], type=shrinkage)
        rd[rownames(dfs), paste0(prefix, "_", shrinkage, "_log2FC")] = dfs$log2FoldChange
    }
    SummarizedExperiment::rowData(se) = rd
    if ("deg" %in% names(S4Vectors::metadata(se))) {
        pval = df$pvalue
        pval[pval < 1e-300] = 1e-300
        pval[is.na(df$padj)] = NA
        S4Vectors::metadata(se)$deg[[paste0(prefix, "_harmonic_mean_pvalue")]] = 1./mean(1./pval, na.rm=TRUE)
    }
    return(se)
}

# 7/29
# move over to new repo
# test file in there too
# make a new package
# keep the tests!
# duplicated columns: can use linear combos to help out with that (proper error messaging)
# 8/1 1:30 eastern time

#' @export
filter.design <- function(design, max.ncol = -1, remove_linear_combos = FALSE, remove_overspecified = FALSE, remove_duplicates = FALSE) {
    # determine relationship message between columns
    get_relationship_message <- function(to_remove, col_names) {
        if (length(to_remove) == 1 && length(col_names) >= 3) {
            remaining_cols <- col_names[-to_remove]
            components <- remaining_cols[1:min(2, length(remaining_cols))]
            return(paste0("components: ", paste(components, collapse = ", ")))
        } else {
            return("Unable to determine specific relationships")
        }
    }
    
    # Duplicate detection and removal
    # Find duplicate columns using both forward and backward duplicate detection
    duplicate_cols <- duplicated(t(design)) | duplicated(t(design), fromLast = TRUE)
    duplicate_indices <- which(duplicate_cols)
    
    if (length(duplicate_indices) > 0) {
        # Create unique groups of duplicate columns
        # This prevents creating multiple groups for the same set of duplicates
        processed_cols <- logical(ncol(design))
        duplicate_groups <- list()
        group_counter <- 0
        
        # Iterate through each duplicate column
        for (col_idx in duplicate_indices) {
            # Skip if this column has already been processed
            if (!processed_cols[col_idx]) {
                # Find all columns identical to the current column
                identical_cols <- which(apply(design, 2, function(x) all(x == design[, col_idx])))
                
                # Only create a group if there are multiple identical columns
                if (length(identical_cols) > 1) {
                    group_counter <- group_counter + 1
                    duplicate_groups[[group_counter]] <- identical_cols
                    
                    # Mark all columns in this group as processed to avoid duplicates
                    processed_cols[identical_cols] <- TRUE
                }
            }
        }
        
        # Handle duplicates based on user preference
        if (remove_duplicates) {
            # Remove duplicate columns, keeping the first occurrence in each group
            cols_to_remove <- integer(0)
            
            for (group in duplicate_groups) {
                # Keep the first column, remove the rest
                cols_to_remove <- c(cols_to_remove, group[-1])
            }
            
            # Remove the duplicate columns if any were found
            if (length(cols_to_remove) > 0) {
                duplicate_colnames <- colnames(design)[cols_to_remove]
                message("Removing duplicate columns: ", paste(duplicate_colnames, collapse = ", "))
                
                # Remove columns from design matrix
                design <- design[, -cols_to_remove, drop = FALSE]
                
                # Update assign attribute if it exists
                if (!is.null(attr(design, "assign"))) {
                    attr(design, "assign") <- attr(design, "assign")[-cols_to_remove]
                }
            }
        } else {
            # Create detailed error message showing duplicate groups
            err_msg <- "Duplicate columns detected in design matrix:\n"
            
            for (group_idx in seq_along(duplicate_groups)) {
                group <- duplicate_groups[[group_idx]]
                col_names <- colnames(design)[group]
                err_msg <- paste0(err_msg, "Group ", group_idx, ": ", paste(col_names, collapse = ", "), "\n")
            }
            
            err_msg <- paste0(err_msg, "\nPlease remove these duplicate columns or set remove_duplicates = TRUE.")
            stop(err_msg)
        }
    }
    
    # remove overspecified
    overspecified <- ncol(design) > nrow(design)
    if (overspecified) {
        n_remove <- ncol(design) - nrow(design)
        col_sds <- apply(design, 2, sd)
        remove_order <- order(col_sds, seq_along(col_sds))
        to_remove_idx <- remove_order[1:n_remove]
        overspecified_colnames <- colnames(design)[to_remove_idx]
        overspecified_sds <- col_sds[to_remove_idx]
        
        if (remove_overspecified) {
            message("Removing overspecified columns: ", paste(overspecified_colnames, collapse = ", "))
            message("Standard deviations of removed columns: ", paste(overspecified_sds, collapse = ", "))
            design <- design[, -to_remove_idx, drop = FALSE]
            if (!is.null(attr(design, "assign"))) {
                attr(design, "assign") <- attr(design, "assign")[-to_remove_idx]
            }
        } else {
            err_msg <- if (length(overspecified_colnames) > 0) {
                paste0(
                    "Overspecified columns detected in design matrix:\n",
                    paste(overspecified_colnames, collapse = "\n"),
                    "\nPlease remove these columns or set remove_overspecified = TRUE."
                )
            } else {
                paste0(
                    "Design matrix has more columns (", ncol(design), ") than rows (", nrow(design), ").\n",
                    "Please be mindful of the number of columns in your design matrix.\n",
                    "Consider removing some columns or set remove_overspecified = TRUE."
                )
            }
            stop(err_msg)
        }
    }
    
    # Near-zero variance removal
    near_zero_var = caret::nearZeroVar(design, saveMetrics=TRUE)
    design = design[, !near_zero_var$nzv | rownames(near_zero_var) == "(Intercept)", drop = FALSE]
    
    # Linear combination removal or error
    linear_combos = caret::findLinearCombos(design)
    to_remove = linear_combos$remove
    
    if (remove_linear_combos) {
        repeat {
            to_remove <- caret::findLinearCombos(design)$remove
            if (length(to_remove) == 0) break
            
            col_names <- colnames(design)
            relationship_msg <- get_relationship_message(to_remove, col_names)
            
            message("Removing linear combination columns: ", paste(col_names[to_remove], collapse = ", "))
            message("Linear combination relationships: ", relationship_msg)
            
            design <- design[, -to_remove, drop = FALSE]
        }
    } else if (length(to_remove) > 0) {
        col_names <- colnames(design)
        relationship_msg <- get_relationship_message(to_remove, col_names)
        
        err_msg <- paste0(
            "Linear combinations detected in design matrix at columns: ",
            paste(to_remove, collapse = ", "),
            "\nColumn names: ", paste(col_names[to_remove], collapse = ", "),
            "\nLinear combination relationships: ", relationship_msg,
            "\nPlease remove these columns or set remove_linear_combos = TRUE."
        )
        stop(err_msg)
    }
    
    # Return modified design if removal requested, otherwise TRUE
    if (remove_overspecified || remove_linear_combos || remove_duplicates) {
        return(design)
    } else {
        return(TRUE)
    }
}

#old
#' @return 
#' @export
old_deseq2 <- function(se,
                       pathology,
                       case,
                       control,
                       sample.col,
                       covariates=NULL, independentFiltering=as.logical(Sys.getenv("DESEQ2_INDEPENDENT", "TRUE")=="TRUE"),
                       shrinkage=Sys.getenv("DESEQ2_SHRINKAGE", "ashr"),
                       prefix="DESeq2") {
    pb = calculate_qc_metrics(se_make_pseudobulk(se, sample.col), assay="counts", qc_vars=c("mt", "ribo", "pc", "chrX", "chrY"))
    X = SummarizedExperiment::assays(pb)$counts
    cd = as.data.frame(SummarizedExperiment::colData(pb))
    cd[[pathology]] = relevel(as.factor(cd[[pathology]]), ref=control)
    covariates = covariates[covariates %in% colnames(cd)]
    covariates = covariates[covariates %in% colnames(deg.filter.design(cd[c(covariates)]))]
    design = model.matrix(as.formula(paste0("~", c(pathology, covariates), collapse=" + ")),
                          data=cd)
    design = deg.filter.design(design)
    dds = DESeq2::DESeqDataSetFromMatrix(
                      X,
                      colData=cd,
                      design=design)
    out = DESeq2::DESeq(dds)
    df = DESeq2::results(out, list(make.names(paste0(pathology, case))), independentFiltering=independentFiltering)
    rd = as.data.frame(SummarizedExperiment::rowData(se))
    rd[[paste0(prefix, "_log2FC")]] = NA
    rd[rownames(df), paste0(prefix, "_log2FC")] = df$log2FoldChange
    rd[[paste0(prefix, "_FDR")]] = NA
    rd[rownames(df), paste0(prefix, "_FDR")] = df$padj
    rd[[paste0(prefix, "_stat")]] = NA
    rd[rownames(df), paste0(prefix, "_stat")] = df$stat
    if (shrinkage %in% c("apeglm", "ashr")) {
        cat("From coefficients\n")
        cat("\t", paste0(DESeq2::resultsNames(out), collapse=","), "\n")
        if (any(grepl(make.names(paste0(pathology, case)), DESeq2::resultsNames(out)))) {
            IDX = grep(make.names(paste0(pathology, case)), DESeq2::resultsNames(out))[1]
        } else {
            IDX = 1
        }
        cat("Shrinking", DESeq2::resultsNames(out)[IDX], "\n")
        dfs <- DESeq2::lfcShrink(out, res=df, coef=DESeq2::resultsNames(out)[IDX], type=shrinkage)
        rd[rownames(dfs), paste0(prefix, "_", shrinkage, "_log2FC")] <- dfs$log2FoldChange
    }
    SummarizedExperiment::rowData(se) <- rd
    
    
    if ("deg" %in% names(S4Vectors::metadata(se))) {
        pval <- df$pvalue
        pval[pval < 1e-300] <- 1e-300
        pval[is.na(df$padj)] <- NA
        S4Vectors::metadata(se)$deg[[paste0(prefix, "_harmonic_mean_pvalue")]] <- 1./mean(1./pval, na.rm=TRUE)
    }
    return(se)
}