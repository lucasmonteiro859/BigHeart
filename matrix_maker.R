library(DESeq2)

load_matrix <- function(txt) {
  counts_matrix <- read.delim(txt, comment.char = "#", check.names = FALSE)
  rownames(counts_matrix) <- counts_matrix$Geneid
  counts_matrix <- counts_matrix[, 7:ncol(counts_matrix)]  
  return(counts_matrix)
}

rename_matrix_ss <- function(counts_matrix_ss, condition) {
  long_names <- colnames(counts_matrix_ss)
  sham_count <- 0
  band_count <- 0
  simple_names <- sapply(condition, function(x) {
    if (x == "Sham") {
      sham_count <<- sham_count + 1
      paste0("Sham_", sham_count)
    } else {
      band_count <<- band_count + 1
      paste0("Band_", band_count)
    }
  })
  colnames(counts_matrix_ss) <- simple_names
  names(condition) <- simple_names  
  return(counts_matrix_ss)
}

rename_matrix_hs <- function(counts_matrix_hs, condition) {
  long_names <- colnames(counts_matrix_hs)
  simple_names <- sapply(long_names, function(x) {
    if (grepl("Ctrl", x)) {
      num <- sub(".*Ctrl_0*([0-9]+).*", "\\1", x)
      paste0("Ctrl_", num)
    } else if (grepl("_HCM_", x)) {
      num <- sub(".*_HCM_0*([0-9]+).*", "\\1", x)
      paste0("HCM_", num)
    } else if (grepl("_EA_", x)) {
      num <- sub(".*_EA_0*([0-9]+).*", "\\1", x)
      paste0("EA_", num)
    } else {
      x  
    }
  })
  colnames(counts_matrix_hs) <- simple_names
  names(condition) <- simple_names  
  return(counts_matrix_hs)
}

create_dds <- function(counts_matrix, condition) {
  coldata <- data.frame(
    row.names = colnames(counts_matrix),
    condition = factor(condition)
  )
  dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = coldata,
    design = ~ condition
  )
  return(dds)
}

norm_matrix <- function(dds) {
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  norm_counts <- as.matrix(norm_counts)
  return(norm_counts)
}

vst_matrix <- function(dds) {
  vst_dds <- vst(dds, blind = TRUE) 
  vst_counts <- assay(vst_dds)
  return(vst_counts)
}

run_DE <- function(dds) {
  dds <- DESeq(dds)
  res <- results(dds)               
  res <- lfcShrink(dds, coef=2)     
  return(res)
}



counts_matrix_hs <- load_matrix("data/counts_noMT.txt")
condition_hs <- ifelse(grepl("Ctrl", colnames(counts_matrix_hs)), "Ctrl", ifelse(grepl("EA_", colnames(counts_matrix_hs)), "EA", "HCM"))
counts_matrix_hs <- rename_matrix_hs(counts_matrix_hs, condition_hs)

remove_idx_hs <- which(colnames(counts_matrix_hs) %in% c("HCM_5", "HCM_10"))
remove_idx_hs <- c(remove_idx_hs, which(colnames(counts_matrix_hs) %in% c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Ctrl_4", "Ctrl_5")))
counts_matrix_hs <- counts_matrix_hs[, -remove_idx_hs, drop = FALSE]
condition_hs <- condition_hs[-remove_idx_hs]

counts_matrix_hs <- counts_matrix_hs[rowSums(counts_matrix_hs  > 10) >= 10, ]
dds_hs <- create_dds(counts_matrix_hs, condition_hs)
counts_matrix_hs <- norm_matrix(dds_hs)
counts_matrix_hs <- log2(counts_matrix_hs + 1)
# counts_matrix_hs <- vst_matrix(counts_matrix_hs, dds_hs)
# View(counts_matrix_hs)


counts_matrix_ss <- load_matrix("data/pig_counts.txt")
condition_ss <- ifelse(grepl("Ctrl", colnames(counts_matrix_ss)), "Sham", "Band")
counts_matrix_ss <- rename_matrix_ss(counts_matrix_ss, condition_ss)

remove_idx_ss <- which(colnames(counts_matrix_ss) %in% c("Band_1"))
counts_matrix_ss <- counts_matrix_ss[, -remove_idx_ss, drop = FALSE]
condition_ss <- condition_ss[-remove_idx_ss]

counts_matrix_ss <- counts_matrix_ss[rowSums(counts_matrix_ss  > 10) >= 4, ]
dds_ss <- create_dds(counts_matrix_ss, condition_ss)
# counts_matrix_ss <- norm_matrix(dds_ss)
# counts_matrix_ss <- log2(counts_matrix_ss + 1)
# counts_matrix_ss <- vst_matrix(counts_matrix_ss, dds_ss)
# View(counts_matrix_ss)

