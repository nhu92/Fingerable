# Load required libraries
library(Biostrings)
library(msa)
library(ape)
library(dplyr)
library(tidyr)
library(stringr)
library(combinat)

# Extract the shortest circular interval between two genes or extract full sequence of one gene
extract_shortest_interval_between_genes <- function(gene1, gene2 = NULL, input_dir = ".", output_fasta) {
  files <- list.files(input_dir, pattern = "\\.gb$", full.names = TRUE)
  fasta_seqs <- list()
  fasta_names <- character()
  
  for (file in files) {
    lines <- readLines(file)
    
    # Extract genome sequence
    origin_start <- grep("^ORIGIN", lines)
    origin_seq <- lines[(origin_start + 1):length(lines)]
    origin_seq <- gsub("[0-9]", "", origin_seq)
    origin_seq <- gsub(" ", "", origin_seq)
    origin_seq <- tolower(paste(origin_seq, collapse = ""))
    genome_len <- nchar(origin_seq)
    
    # Parse features
    feature_lines <- grep("^[[:space:]]{5}(CDS|gene)[[:space:]]", lines)
    gene1_coords <- list()
    gene2_coords <- list()
    
    for (idx in feature_lines) {
      type_line <- lines[idx]
      location_text <- gsub(".*?([0-9,<>()\\.]+)", "\\1", type_line)
      
      qualifiers <- lines[(idx + 1):min(idx + 20, length(lines))]
      gene_val <- sub('.*?/gene="(.*?)".*', "\\1", grep("/gene=", qualifiers, value = TRUE)[1])
      product_val <- sub('.*?/product="(.*?)".*', "\\1", grep("/product=", qualifiers, value = TRUE)[1])
      gene_val[is.na(gene_val)] <- ""
      product_val[is.na(product_val)] <- ""
      
      matches <- gregexpr("[0-9]+\\.\\.[0-9]+", location_text, perl = TRUE)
      joined <- unlist(regmatches(location_text, matches))
      coords_list <- lapply(joined, function(j) as.numeric(strsplit(j, "..", fixed = TRUE)[[1]]))
      
      if (gene_val == gene1 || grepl(gene1, product_val)) {
        gene1_coords <- c(gene1_coords, list(coords_list))
      }
      
      if (!is.null(gene2) && (gene_val == gene2 || grepl(gene2, product_val))) {
        gene2_coords <- c(gene2_coords, coords_list)
      }
    }
    
    if (length(gene1_coords) == 0 || (!is.null(gene2) && length(gene2_coords) == 0)) {
      message(sprintf("⚠️ Missing gene(s) in %s", file))
      next
    }
    
    base_name <- tools::file_path_sans_ext(basename(file))
    
    if (is.null(gene2)) {
      # Merge overlapping or contained blocks
      merged_blocks <- list()
      for (block in gene1_coords) {
        new_range <- range(unlist(block))
        overlap <- FALSE
        for (i in seq_along(merged_blocks)) {
          existing_range <- merged_blocks[[i]]
          if ((new_range[1] >= existing_range[1] && new_range[1] <= existing_range[2]) ||
              (new_range[2] >= existing_range[1] && new_range[2] <= existing_range[2]) ||
              (new_range[1] <= existing_range[1] && new_range[2] >= existing_range[2])) {
            merged_blocks[[i]] <- c(min(existing_range[1], new_range[1]), max(existing_range[2], new_range[2]))
            overlap <- TRUE
            break
          }
        }
        if (!overlap) {
          merged_blocks <- c(merged_blocks, list(new_range))
        }
      }
      
      # Extract and stitch sequences
      seq_blocks <- list()
      for (block in merged_blocks) {
        seq <- substr(origin_seq, block[1], block[2])
        seq_blocks <- c(seq_blocks, list(seq))
      }
      
      final_seq <- paste(seq_blocks, collapse = paste(rep("N", 10), collapse = ""))
      
      if (nchar(final_seq) > 10000) {
        message(sprintf("⚠️ Sequence for %s in %s exceeds 10k bp, skipping.", gene1, file))
        next
      }
      
      header <- sprintf(">%s|gene=%s|copies=%d", base_name, gene1, length(seq_blocks))
      fasta_names <- c(fasta_names, header)
      fasta_seqs <- c(fasta_seqs, final_seq)
    } else {
      # Two-gene mode: find shortest circular interval
      gene1_positions <- sapply(gene1_coords, function(x) mean(x[[1]]))
      gene2_positions <- sapply(gene2_coords, function(x) mean(x))
      
      shortest <- Inf
      best_start <- NULL
      best_end <- NULL
      best_dir <- NULL
      for (p1 in gene1_positions) {
        for (p2 in gene2_positions) {
          d1 <- (p2 - p1 + genome_len) %% genome_len
          d2 <- (p1 - p2 + genome_len) %% genome_len
          if (d1 <= d2 && d1 < shortest) {
            shortest <- d1; best_start <- ceiling(p1); best_end <- ceiling(p2); best_dir <- "forward"
          } else if (d2 < d1 && d2 < shortest) {
            shortest <- d2; best_start <- ceiling(p2); best_end <- ceiling(p1); best_dir <- "reverse"
          }
        }
      }
      
      if (best_start <= best_end) {
        seq <- substr(origin_seq, best_start, best_end)
      } else {
        seq <- paste0(substr(origin_seq, best_start, genome_len), substr(origin_seq, 1, best_end))
      }
      
      if (nchar(seq) > 10000) {
        message(sprintf("⚠️ Sequence from %s to %s in %s exceeds 10k bp, skipping.", gene1, gene2, file))
        next
      }
      
      header <- sprintf(">%s|start=%d|end=%d|dir=%s|from=%s|to=%s", base_name, best_start, best_end, best_dir, gene1, gene2)
      fasta_names <- c(fasta_names, header)
      fasta_seqs <- c(fasta_seqs, seq)
    }
  }
  
  if (length(fasta_seqs) > 0) {
    fasta_lines <- unlist(mapply(function(h, s) c(h, s), fasta_names, fasta_seqs, SIMPLIFY = FALSE))
    writeLines(fasta_lines, con = output_fasta)
    message(sprintf("✅ Extracted sequences written to %s", output_fasta))
  } else {
    message("❌ No sequences extracted.")
  }
}

# Function to align fasta sequences and calculate pairwise distances
align_and_calc_distance <- function(input_fasta, output_distances) {
  # Read sequences
  seqs <- readDNAStringSet(input_fasta)
  
  # Multiple sequence alignment
  alignment <- msa(seqs, method = "ClustalW")
  
  # Convert alignment to ape alignment format
  aligned_seqs <- as.DNAbin(alignment)
  
  # Calculate pairwise distance matrix
  dist_matrix <- dist.dna(aligned_seqs, model = "K80")  # Kimura 2-parameter
  
  # Convert matrix to a long-form data frame
  dist_df <- as.data.frame(as.table(as.matrix(dist_matrix)))
  colnames(dist_df) <- c("Sequence1", "Sequence2", "Distance")
  
  # Remove self-comparisons and duplicate comparisons
  dist_df <- dist_df[as.character(dist_df$Sequence1) < as.character(dist_df$Sequence2), ]
  
  # Write distances to a TSV file
  write.table(dist_df, file = output_distances, sep = "\t", row.names = FALSE, quote = FALSE)
  
  message(sprintf("✅ Pairwise distances written to %s", output_distances))
}

classify_species_by_distance <- function(tsv_file, offhit_threshold = 0.6, offhit_ratio = 0.7, indist_threshold = 0.01) {
  df <- read.delim(tsv_file, stringsAsFactors = FALSE)
  
  # Step 1: identify off-hit species
  all_species <- unique(c(df$Sequence1, df$Sequence2))
  offhit_species <- c()
  
  for (sp in all_species) {
    related <- df[df$Sequence1 == sp | df$Sequence2 == sp, ]
    too_far <- sum(related$Distance > offhit_threshold)
    total <- nrow(related)
    if (total > 0 && too_far / total >= offhit_ratio) {
      offhit_species <- c(offhit_species, sp)
    }
  }
  
  # Step 2: filter out off-hit rows
  df_filtered <- df[!(df$Sequence1 %in% offhit_species | df$Sequence2 %in% offhit_species), ]
  kept_species <- unique(c(df_filtered$Sequence1, df_filtered$Sequence2))
  
  # Step 3: identify not-distinguishable species (those involved in any close pair)
  close_pairs <- df_filtered[df_filtered$Distance < indist_threshold, ]
  not_distinguishable <- unique(c(close_pairs$Sequence1, close_pairs$Sequence2))
  
  # Step 4: the rest are distinguishable
  distinguishable <- setdiff(kept_species, not_distinguishable)
  
  # Step 5: print counts and return table
  cat("🧬 Final unique species classification:\n")
  cat("  Off-hit:               ", length(unique(offhit_species)), "\n")
  cat("  Not distinguishable:   ", length(unique(not_distinguishable)), "\n")
  cat("  Distinguishable:       ", length(unique(distinguishable)), "\n")
  
  # Output full classification table
  final_status <- data.frame(
    Species = all_species,
    Status = ifelse(all_species %in% offhit_species, "Off-hit",
                    ifelse(all_species %in% not_distinguishable, "Not distinguishable", "Distinguishable"))
  )
  
  return(final_status)
}

combine_and_evaluate_classification <- function(folder_path) {
  files <- list.files(folder_path, pattern = "classification.tsv$", full.names = TRUE)
  
  # Step 1: 读取每个 gene 的分类表
  gene_tables <- list()
  for (file in files) {
    df <- read.delim(file, stringsAsFactors = FALSE)
    df$Species <- sub("\\|.*", "", df$Species)
    gene_name <- tools::file_path_sans_ext(basename(file))
    gene_tables[[gene_name]] <- df[, c("Species", "Status")]
  }
  
  # 所有基因名
  gene_names <- names(gene_tables)
  
  # Step 2: 构建完整物种 × 基因分类表
  merged_df <- Reduce(function(x, y) merge(x, y, by = "Species", all = TRUE), gene_tables)
  colnames(merged_df)[-1] <- gene_names
  
  # Step 3: 定义分类逻辑函数
  classify_combination <- function(status_matrix) {
    apply(status_matrix, 1, function(statuses) {
      if ("Distinguishable" %in% statuses) {
        return("Distinguishable")
      } else if ("Not distinguishable" %in% statuses) {
        return("Not distinguishable")
      } else {
        return("Off-hit")
      }
    })
  }
  
  # Step 4: 对所有组合进行分析
  all_combinations <- unlist(lapply(1:length(gene_names), function(k) combn(gene_names, k, simplify = FALSE)), recursive = FALSE)
  
  # 统计每组组合的结果
  results <- data.frame()
  for (combo in all_combinations) {
    combo_label <- paste(combo, collapse = " + ")
    subset_status <- merged_df[, combo, drop = FALSE]
    final_status <- classify_combination(subset_status)
    counts <- table(factor(final_status, levels = c("Distinguishable", "Not distinguishable", "Off-hit")))
    results <- rbind(results, data.frame(
      Combination = combo_label,
      Distinguishable = counts["Distinguishable"],
      Not_distinguishable = counts["Not distinguishable"],
      Off_hit = counts["Off-hit"]
    ))
  }
  
  # 显示总表
  cat("📊 Summary for all gene combinations:\n")
  print(results, row.names = FALSE)
  
  return(results)
}

# # -------------
# 
# g_name <- ""
# 
# extract_shortest_interval_between_genes(g_name, input_dir = "./", output_fasta = paste0(g_name, ".fasta"))
# 
# align_and_calc_distance(paste0(g_name,".fasta"), paste0(g_name, "_distance.tsv"))
# 
# result <- classify_species_by_distance(paste0(g_name, "_distance.tsv"))
# 
# write.table(result, paste0(g_name, "_classification.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
# 
# # -------------
# 
# g_name1 <- ""
# g_name2 <- ""
# 
# extract_shortest_interval_between_genes(g_name1, g_name2, input_dir = "./", output_fasta = paste0(g_name1, "_", g_name2, ".fasta"))
# 
# align_and_calc_distance(paste0(g_name1, "_", g_name2, ".fasta"), paste0(g_name1, "_", g_name2,  "_distance.tsv"))
# 
# result <- classify_species_by_distance(paste0(g_name1, "_", g_name2,  "_distance.tsv"))
# 
# write.table(result, paste0(g_name1, "_", g_name2,  "_classification.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
# 
# # -------------
# 
# results <- combine_and_evaluate_classification("./")
# write.table(results, file = "gene_combination_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# 
# 
