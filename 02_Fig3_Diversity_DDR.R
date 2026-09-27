# ==============================================================================
# SCRIPT 02: TRI-KINGDOM DIVERSITY, TAXONOMY DONUTS & DISTANCE-DECAY (DDR)
# Reproducible Pipeline for Figure 3 (A–J), Figure S3 (A–E), Table S6
# Top-Tier Publication Standard (Nature / ISME Grade)
# ==============================================================================

# ----------------- 0. 环境准备与依赖包检测 -----------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, readxl, data.table, vegan, patchwork, cowplot,
  multcomp, multcompView, betapart, RColorBrewer, picante, ape,
  ggsci, svglite, readr, ggrastr
)

# 智能加载 linkET (用于 Figure S3A 的 Mantel 互作网络)
if (!require("linkET")) {
  message(">>> Installing linkET package...")
  tryCatch(install.packages("linkET"), error = function(e) {
    if (!require("devtools")) install.packages("devtools")
    devtools::install_github("Hy4m/linkET")
  })
  library(linkET)
}

select <- dplyr::select
filter <- dplyr::filter

# 相对工程路径配置 (杜绝任何本地绝对路径)
data_dir   <- "./data"
output_dir <- "./results/Figure3_Diversity"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

safe_write_csv <- function(df, path) {
  tryCatch({
    readr::write_excel_csv(df, path)
    message(">>> [Successfully Exported Table] ", path)
  }, error = function(e) {
    warning("Cannot write to file: ", path, " : ", e$message)
  })
}

safe_read <- function(file_path) {
  tryCatch({
    readxl::read_excel(file_path)
  }, error = function(e) {
    tryCatch({
      readr::read_delim(file_path, delim = "\t", show_col_types = FALSE)
    }, error = function(e2) {
      readr::read_csv(file_path, show_col_types = FALSE)
    })
  })
}

# 智能文件查找器 (支持横杠与下划线命名兼容)
find_file_safe <- function(dir, patterns, default_name) {
  target <- file.path(dir, default_name)
  if (file.exists(target)) return(target)
  for (pat in patterns) {
    hits <- list.files(dir, pattern = pat, full.names = TRUE, ignore.case = TRUE)
    if (length(hits) > 0) return(hits[1])
  }
  return(target)
}

# ----------------- 1. 顶刊极简主题系统 (8.0 pt 出版标准) -----------------
font_anno <- 2.0

theme_nature_v44_grid <- function(base_size = 8.0, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(size = base_size, family = base_family),
      panel.grid.major = element_line(color = "#EAECEE", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.40),
      strip.background = element_rect(fill = "#F2F4F4", color = "black", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 6.8, color = "black", margin = margin(t = 1.5, b = 1.5)),
      plot.title = element_text(face = "bold", size = 8.0, hjust = 0.5, margin = margin(b = 2.5, t = 0)),
      axis.text = element_text(color = "black", size = 6.8),
      axis.title = element_text(face = "bold", size = 7.5),
      plot.tag = element_text(face = "bold", size = 11.0),
      legend.title = element_text(face = "bold", size = 6.8, margin = margin(t = 1.5, b = 3.0, unit = "pt")),
      legend.text = element_text(size = 6.0),
      legend.key.size = unit(0.55, "lines"),
      legend.spacing.y = unit(4.5, "pt"),
      legend.box.spacing = unit(5.0, "pt"),
      legend.margin = margin(2, 2, 2, 2, "pt"),
      panel.spacing = unit(0.25, "lines"),
      plot.margin = margin(t = 2.0, r = 3.0, b = 2.0, l = 3.0, "pt")
    )
}

theme_set(theme_nature_v44_grid())

# 顶刊色彩与标记配置
report_cols    <- c("High" = "#2E8B57", "Med-High" = "#FF3B6E", "Med" = "#55B9FF", "Low" = "#FFC533")
veg_cols       <- c("Herb" = "#4A7C59", "Shrub" = "#D07A30", "Tree" = "#1F3A52")
veg_shapes     <- c("Herb" = 16, "Shrub" = 17, "Tree" = 15)
veg_lines      <- c("Herb" = "solid", "Shrub" = "dashed", "Tree" = "dotdash")
kingdom_cols   <- c("Bacteria" = "#BC3C29", "Fungi" = "#0072B5", "Protist" = "#E18727")
component_cols <- c("RichDiff" = "#BC3C29", "Replacement" = "#0072B5")

donut_palette_base <- c(
  "#1F4E79", "#D9534F", "#2E8B57", "#D2691E", 
  "#6A5ACD", "#E67E22", "#4682B4", "#8B4513", 
  "#95A5A6", "#BDC3C7"
)

# ----------------- 2. 自动化静默读取数据 (无弹窗) -----------------
env_file_path  <- find_file_safe(data_dir, c("(理化|Enzyme|Physico|env).*\\.xlsx?$"), "Soil_Physicochemical_Enzymes_4Seasons.xlsx")
bac_file_path  <- find_file_safe(data_dir, c(".*(bacteria|bac).*\\.(xls|xlsx|txt|tsv|csv)$"), "ASV_Taxon_Table_bacteria.xls")
fun_file_path  <- find_file_safe(data_dir, c(".*(fungi|fun).*\\.(xls|xlsx|txt|tsv|csv)$"), "ASV_Taxon_Table_fungi.xls")
prot_file_path <- find_file_safe(data_dir, c(".*(protist|prot).*\\.(xls|xlsx|txt|tsv|csv)$"), "ASV_Taxon_Table_protists_only.xlsx")

if (!file.exists(env_file_path)) stop(paste0("Environmental data file not found in: ", data_dir))
if (!file.exists(bac_file_path)) stop(paste0("Bacterial ASV table not found in: ", data_dir))
if (!file.exists(fun_file_path)) stop(paste0("Fungal ASV table not found in: ", data_dir))
if (!file.exists(prot_file_path)) stop(paste0("Protistan ASV table not found in: ", data_dir))

sheet_names <- readxl::excel_sheets(env_file_path)
env_raw <- lapply(sheet_names, function(sheet) {
  readxl::read_excel(env_file_path, sheet = sheet)
}) %>% bind_rows()

colnames(env_raw) <- gsub("NH4\\+-N|NH4_N", "NH4_N", colnames(env_raw))
colnames(env_raw) <- gsub("swc|SWC", "SWC", colnames(env_raw))
if ("Sample ID" %in% colnames(env_raw)) env_raw <- env_raw %>% rename(SampleID = `Sample ID`)

numeric_cols <- c(
  "SOM", "NH4_N", "AP", "AK", "TN", "TP", "TK", "pH", "SWC", 
  "SUE", "SSC", "SALP", 
  "Pb", "Cr", "Cu", "Ni", "Zn", "Cd", "As"
)
env_raw <- env_raw %>% mutate(across(any_of(numeric_cols), ~ as.numeric(as.character(.x))))

bac_raw  <- safe_read(bac_file_path)
fun_raw  <- safe_read(fun_file_path)
prot_raw <- safe_read(prot_file_path)

season_map <- c(C="Spring", X="Summer", Q="Autumn", D="Winter")
urban_map  <- c(A="High", B="Med-High", C="Med", D="Low")
veg_map    <- c(C="Herb", G="Shrub", Q="Tree")

clean_id <- gsub("_|-|\\s+", "", env_raw$SampleID)
env_classified <- env_raw %>%
  mutate(
    Season_code = substr(clean_id, 1, 1),
    Urban_code  = substr(clean_id, 2, 2),
    Veg_code    = substr(clean_id, 3, 3),
    Season = season_map[Season_code],
    Urban  = urban_map[Urban_code],
    Veg    = veg_map[Veg_code]
  ) %>%
  mutate(
    Season = factor(Season, levels = c("Spring", "Summer", "Autumn", "Winter")),
    Urban  = factor(Urban,  levels = c("High", "Med-High", "Med", "Low")), 
    Veg    = factor(Veg,    levels = c("Herb", "Shrub", "Tree"))
  )

process_otu_table <- function(df) {
  tax_cols <- c("taxonomy", "Taxon", "Taxonomy", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  df_numeric <- df %>% dplyr::select(-any_of(tax_cols))
  otu_matrix <- as.matrix(df_numeric[, -1])
  rownames(otu_matrix) <- df[[1]]
  return(t(otu_matrix))
}

bac_counts  <- process_otu_table(bac_raw)
fun_counts  <- process_otu_table(fun_raw)
prot_counts <- process_otu_table(prot_raw)

target_env_vars <- intersect(numeric_cols, colnames(env_classified))
env_clean <- env_classified %>% 
  dplyr::select(SampleID, Season, Urban, Veg, all_of(target_env_vars)) %>% 
  na.omit() %>%
  mutate(SampleID = as.character(SampleID))

common_samples <- intersect(rownames(bac_counts), rownames(fun_counts)) %>% 
  intersect(rownames(prot_counts)) %>% 
  intersect(env_clean$SampleID) %>%
  sort()

message(paste(">>> Successfully aligned", length(common_samples), "samples across all kingdoms!"))

env_aligned  <- env_clean %>% 
  filter(SampleID %in% common_samples) %>% 
  arrange(SampleID) %>%
  mutate(SampleID = as.character(SampleID))

bac_sub_rel  <- decostand(bac_counts[common_samples, ], method = "total") * 100
fun_sub_rel  <- decostand(fun_counts[common_samples, ], method = "total") * 100
prot_sub_rel <- decostand(prot_counts[common_samples, ], method = "total") * 100

bac_sub_cnt  <- bac_counts[common_samples, ]
fun_sub_cnt  <- fun_counts[common_samples, ]
prot_sub_cnt <- prot_counts[common_samples, ]

# ----------------- 3. 分类学规范与圆环图引擎 -----------------
clean_taxon_name <- function(name) {
  if (is.null(name) || length(name) == 0) return("Unclassified")
  name <- gsub("^[kpcofgs]__", "", name)
  name <- gsub(".*__", "", name)
  name <- trimws(name)
  name[name == "" | is.na(name) | tolower(name) %in% c("norank", "unclassified", "uncultured", "unknown", "metagenome", "null", "na")] <- "Unclassified"
  return(name)
}

standardize_protist_phylum <- function(p_name) {
  p_clean <- clean_taxon_name(p_name)
  if (p_clean == "Unclassified") return("Unclassified")
  
  if (p_clean %in% c("Tubulinea", "Discosea", "Evosea", "Variosea", "Mycetozoa", 
                     "Archamoebae", "Dictyostelia", "Myxogastria", "Protostelia", "Amoebozoa")) {
    return("Amoebozoa")
  }
  if (p_clean %in% c("Endomyxa", "Filosa", "Cercomonadida", "Phytomyxea", 
                     "Thecofilosea", "Imbricatea", "Chlorarachniophyta", "Cercozoa")) {
    return("Cercozoa")
  }
  if (p_clean %in% c("Heterolobosea", "Percolomonadida", "Tetramitia", "Percolozoa")) {
    return("Percolozoa")
  }
  if (p_clean %in% c("Kinetoplastea", "Euglenida", "Diplonemea", "Symbiontida", "Euglenozoa")) {
    return("Euglenozoa")
  }
  if (p_clean %in% c("Conoidasida", "Aconoidasida", "Gregarinasina", "Coccidia", "Apicomplexa")) {
    return("Apicomplexa")
  }
  if (p_clean %in% c("Spirotrichea", "Oligohymenophorea", "Litostomatea", "Colpodea", "Ciliophora")) {
    return("Ciliophora")
  }
  return(p_clean)
}

build_taxonomic_tree <- function(raw_df) {
  tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  present_cols <- intersect(tax_cols, colnames(raw_df))
  
  if (length(present_cols) >= 2) {
    tax_df <- raw_df %>% select(ID = 1, all_of(present_cols))
  } else {
    tax_col <- intersect(c("taxonomy", "Taxon", "Taxonomy"), colnames(raw_df))
    if (length(tax_col) > 0) {
      tax_str <- raw_df[[tax_col[1]]]
      parsed_list <- lapply(tax_str, function(x) {
        parts <- trimws(strsplit(as.character(x), ";|,")[[1]])
        parts <- gsub("^[a-z]__", "", parts)
        if (length(parts) < 7) parts <- c(parts, rep("Unclassified", 7 - length(parts)))
        parts[1:7]
      })
      tax_df <- as.data.frame(do.call(rbind, parsed_list), stringsAsFactors = FALSE)
      colnames(tax_df) <- tax_cols
      tax_df$ID <- raw_df[[1]]
    } else {
      tax_df <- data.frame(ID = raw_df[[1]], Kingdom = "Unclassified", Phylum = "Unclassified", Class = "Unclassified", stringsAsFactors = FALSE)
    }
  }
  
  tax_df_clean <- tax_df %>%
    mutate(across(-ID, ~ {
      x <- as.character(.x)
      x[is.na(x) | x == "" | x == "NA" | tolower(x) %in% c("null", "unclassified", "norank", "uncultured")] <- "Unclassified"
      factor(x)
    }))
  
  nested_cols <- colnames(tax_df_clean)[colnames(tax_df_clean) != "ID"]
  formula_str <- paste0("~", paste(nested_cols, collapse = "/"), "/ID")
  frm <- as.formula(formula_str)
  
  tr <- ape::as.phylo(frm, data = tax_df_clean, collapse = FALSE)
  tr$edge.length <- rep(1, nrow(tr$edge))
  return(tr)
}

calc_faith_pd <- function(counts_matrix, raw_df) {
  tryCatch({
    tr <- build_taxonomic_tree(raw_df)
    common_tips <- intersect(tr$tip.label, counts_matrix %>% colnames())
    if (length(common_tips) < 2) return(rep(0, nrow(counts_matrix)))
    
    tr_pruned <- ape::keep.tip(tr, common_tips)
    counts_pruned <- counts_matrix[, tr_pruned$tip.label, drop = FALSE]
    
    pd_res <- picante::pd(counts_pruned, tr_pruned, include.root = FALSE)
    pd_vector <- pd_res$PD[match(rownames(counts_matrix), rownames(pd_res))]
    pd_vector[is.na(pd_vector)] <- 0
    return(pd_vector)
  }, error = function(e) {
    return(as.numeric(vegan::specnumber(counts_matrix)))
  })
}

get_phylum_donut <- function(raw_df, kingdom_name) {
  tax_cols <- c("taxonomy", "Taxon", "Taxonomy", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_col <- intersect(tax_cols, colnames(raw_df))
  phylum_col <- if ("Phylum" %in% tax_col) "Phylum" else if (length(tax_col) > 0) tax_col[1] else NULL
  
  df_tax <- raw_df %>% select(ID = 1, all_of(phylum_col))
  if (phylum_col %in% c("taxonomy", "Taxonomy", "Taxon")) {
    df_tax <- df_tax %>%
      mutate(Phylum = sapply(get(phylum_col), function(x) {
        parts <- trimws(strsplit(as.character(x), ";|,")[[1]]) 
        p_part <- grep("^[pP]__", parts, value = TRUE)
        if (length(p_part) > 0) gsub("^[pP]__", "", p_part[1]) else if (length(parts) >= 2) parts[2] else if (length(parts) >= 1) parts[1] else "Unclassified"
      }))
  } else {
    df_tax <- df_tax %>% rename(Phylum = all_of(phylum_col))
  }
  
  if (grepl("protist", kingdom_name, ignore.case = TRUE)) {
    df_tax <- df_tax %>% mutate(Phylum = sapply(Phylum, standardize_protist_phylum))
  } else {
    df_tax <- df_tax %>% mutate(Phylum = clean_taxon_name(Phylum))
  }
  
  df_numeric <- raw_df %>% dplyr::select(-any_of(tax_cols))
  samples_in_df <- intersect(colnames(df_numeric), common_samples)
  if (length(samples_in_df) == 0) samples_in_df <- colnames(df_numeric)[-1]
  
  df_joined <- df_tax %>%
    left_join(df_numeric %>% select(ID = 1, all_of(samples_in_df)), by = "ID") %>%
    pivot_longer(cols = all_of(samples_in_df), names_to = "SampleID", values_to = "Abundance") %>%
    mutate(Abundance = as.numeric(Abundance)) %>%
    filter(!is.na(Abundance)) %>%
    group_by(Phylum) %>%
    summarise(Total_Abund = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(Percentage = Total_Abund / sum(Total_Abund) * 100)
  
  df_real <- df_joined %>% filter(Phylum != "Unclassified") %>% arrange(desc(Percentage))
  df_unclass <- df_joined %>% filter(Phylum == "Unclassified")
  
  max_real_phyla <- 8
  if (nrow(df_real) > max_real_phyla) {
    top_real <- df_real %>% head(max_real_phyla)
    others_abund <- sum(df_real$Total_Abund[(max_real_phyla + 1):nrow(df_real)])
    others_perc <- sum(df_real$Percentage[(max_real_phyla + 1):nrow(df_real)])
    df_others <- data.frame(Phylum = "Others", Total_Abund = others_abund, Percentage = others_perc)
    df_plot <- bind_rows(top_real, df_others, df_unclass)
  } else {
    df_plot <- bind_rows(df_real, df_unclass)
  }
  
  final_levels <- df_plot$Phylum[df_plot$Phylum != "Others" & df_plot$Phylum != "Unclassified"]
  if ("Others" %in% df_plot$Phylum) final_levels <- c(final_levels, "Others")
  if ("Unclassified" %in% df_plot$Phylum) final_levels <- c(final_levels, "Unclassified")
  
  df_plot <- df_plot %>%
    filter(Percentage > 0) %>%
    mutate(Phylum = factor(Phylum, levels = final_levels)) %>%
    arrange(Phylum) %>%
    mutate(ymax = cumsum(Percentage), ymin = c(0, head(ymax, -1)))
  
  n_real <- sum(!df_plot$Phylum %in% c("Others", "Unclassified"))
  plot_palette <- donut_palette_base[1:n_real]
  names(plot_palette) <- levels(df_plot$Phylum)[1:n_real]
  if ("Others" %in% df_plot$Phylum) plot_palette["Others"] <- "#95A5A6"
  if ("Unclassified" %in% df_plot$Phylum) plot_palette["Unclassified"] <- "#BDC3C7"
  
  ggplot(df_plot, aes(ymax = ymax, ymin = ymin, xmax = 3.92, xmin = 2.90, fill = Phylum)) +
    geom_rect(color = "white", linewidth = 0.25) +
    coord_polar(theta = "y") + xlim(c(0, 4.0)) + 
    scale_fill_manual(values = plot_palette, name = "Phylum") +
    annotate("text", x = 0, y = 50, label = kingdom_name, size = 1.95, fontface = "bold", color = "#1C2833") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 6.4, margin = margin(b = 2.0, unit = "pt")),
      legend.text = element_text(size = 5.4),
      legend.key.size = unit(0.38, "lines"),
      legend.spacing.y = unit(2.0, "pt"),
      legend.margin = margin(l = -2, r = 1, unit = "pt"),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
    )
}

p_donut_bac  <- get_phylum_donut(bac_raw, "Bacteria")
p_donut_fun  <- get_phylum_donut(fun_raw, "Fungi")
p_donut_prot <- get_phylum_donut(prot_raw, "Protist")

# 安全多重比较字母提取算法 (免疫连字符冲突)
get_clean_tukey_letters <- function(data, value_col, group_col) {
  sub_df <- data %>% filter(!is.na(.data[[value_col]]), !is.na(.data[[group_col]]))
  orig_levels <- levels(factor(sub_df[[group_col]]))
  if (length(orig_levels) <= 1 || nrow(sub_df) < 3) {
    return(data.frame(Group = orig_levels, Letters = "", stringsAsFactors = FALSE))
  }
  
  safe_tokens <- paste0("TOKEN", seq_along(orig_levels))
  forward_map <- setNames(safe_tokens, orig_levels)
  reverse_map <- setNames(orig_levels, safe_tokens)
  
  sub_df$Safe_Group <- factor(forward_map[as.character(sub_df[[group_col]])], levels = safe_tokens)
  fit <- aov(sub_df[[value_col]] ~ Safe_Group, data = sub_df)
  tuk <- TukeyHSD(fit)$Safe_Group
  p_val <- tuk[, "p adj"]
  
  letters_obj <- tryCatch({
    multcompLetters(p_val)$Letters
  }, error = function(e) {
    setNames(rep("a", length(safe_tokens)), safe_tokens)
  })
  
  data.frame(Group = reverse_map[names(letters_obj)], Letters = as.character(letters_obj), stringsAsFactors = FALSE)
}

calc_alpha_all_metrics <- function(counts_matrix, raw_df, prefix) {
  counts_integer <- round(counts_matrix)
  est <- estimateR(counts_integer)
  chao1_vec <- if (is.matrix(est)) est[2, ] else est[2]
  shannon_vec <- vegan::diversity(counts_matrix, index = "shannon")
  rich_vec <- vegan::specnumber(counts_matrix)
  
  evenness_vec <- shannon_vec / log(rich_vec)
  evenness_vec[rich_vec <= 1] <- 0
  evenness_vec[is.na(evenness_vec) | is.infinite(evenness_vec)] <- 0
  
  pd_vec <- calc_faith_pd(counts_matrix, raw_df)
  
  df <- data.frame(
    SampleID = rownames(counts_matrix),
    Chao1 = as.numeric(chao1_vec),
    Pd_faith = as.numeric(pd_vec),
    Shannon = as.numeric(shannon_vec),
    Pielou_J = as.numeric(evenness_vec),
    stringsAsFactors = FALSE
  )
  colnames(df)[2:5] <- paste0(prefix, "_", colnames(df)[2:5])
  return(df)
}

bac_alpha_df  <- calc_alpha_all_metrics(bac_counts[common_samples, ], bac_raw, "Bac")
fun_alpha_df  <- calc_alpha_all_metrics(fun_counts[common_samples, ], fun_raw, "Fun")
prot_alpha_df <- calc_alpha_all_metrics(prot_counts[common_samples, ], prot_raw, "Prot")

dist_bac  <- vegdist(bac_sub_rel, method = "bray")
dist_fun  <- vegdist(fun_sub_rel, method = "bray")
dist_prot <- vegdist(prot_sub_rel, method = "bray")

disp_bac  <- betadisper(dist_bac, env_aligned$Urban)
disp_fun  <- betadisper(dist_fun, env_aligned$Urban)
disp_prot <- betadisper(dist_prot, env_aligned$Urban)

bac_disp_df  <- data.frame(SampleID = names(disp_bac$distances), Bac_BetaDisp = as.numeric(disp_bac$distances), stringsAsFactors = FALSE)
fun_disp_df  <- data.frame(SampleID = names(disp_fun$distances), Fun_BetaDisp = as.numeric(disp_fun$distances), stringsAsFactors = FALSE)
prot_disp_df <- data.frame(SampleID = names(disp_prot$distances), Prot_BetaDisp = as.numeric(disp_prot$distances), stringsAsFactors = FALSE)

env_aligned <- env_aligned %>%
  left_join(bac_alpha_df, by = "SampleID") %>%
  left_join(fun_alpha_df, by = "SampleID") %>%
  left_join(prot_alpha_df, by = "SampleID") %>%
  left_join(bac_disp_df, by = "SampleID") %>%
  left_join(fun_disp_df, by = "SampleID") %>%
  left_join(prot_disp_df, by = "SampleID")

# ----------------- 4. DDR 距离衰减分析 -----------------
calculate_advanced_ddr <- function(comm_matrix, env_df, env_vars, kingdom_name, method = "bray") {
  comm_dist <- as.matrix(vegdist(comm_matrix, method = method))
  comm_sim <- 1 - comm_dist[lower.tri(comm_dist)]
  env_scaled <- scale(env_df[, env_vars])
  env_dist_scaled <- as.matrix(dist(env_scaled, method = "euclidean"))
  env_sim_scaled <- -env_dist_scaled[lower.tri(env_dist_scaled)]
  
  data.frame(
    Kingdom = kingdom_name,
    Distance = ifelse(method == "bray", "Bray-Curtis", "Jaccard"),
    Env_Sim_Scaled = env_sim_scaled,
    Comm_Sim = comm_sim,
    stringsAsFactors = FALSE
  )
}

ddr_bac_bc   <- calculate_advanced_ddr(bac_sub_rel,  env_aligned, target_env_vars, "Bacteria", "bray")
ddr_bac_jac  <- calculate_advanced_ddr(bac_sub_rel,  env_aligned, target_env_vars, "Bacteria", "jaccard")
ddr_fun_bc   <- calculate_advanced_ddr(fun_sub_rel,  env_aligned, target_env_vars, "Fungi",    "bray")
ddr_fun_jac  <- calculate_advanced_ddr(fun_sub_rel,  env_aligned, target_env_vars, "Fungi",    "jaccard")
ddr_prot_bc  <- calculate_advanced_ddr(prot_sub_rel, env_aligned, target_env_vars, "Protist",  "bray")
ddr_prot_jac <- calculate_advanced_ddr(prot_sub_rel, env_aligned, target_env_vars, "Protist",  "jaccard")

ddr_all_modes <- bind_rows(ddr_bac_bc, ddr_bac_jac, ddr_fun_bc, ddr_fun_jac, ddr_prot_bc, ddr_prot_jac)

p_ddr_bc_standard <- ggplot(ddr_all_modes %>% filter(Distance == "Bray-Curtis"), 
                            aes(x = Env_Sim_Scaled, y = Comm_Sim, color = Kingdom)) +
  ggrastr::rasterise(geom_point(alpha = 0.15, size = 0.45), dpi = 600) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = Kingdom), linewidth = 0.65) +
  scale_color_manual(values = kingdom_cols) +
  scale_fill_manual(values = kingdom_cols) +
  labs(x = "Env. Similarity (Scaled)", y = "Community Similarity\n(1 - Bray-Curtis)", title = "Distance Decay Relationship") +
  theme_nature_v44_grid() +
  theme(legend.position = "right")

# ----------------- 5. Podani Beta 多样性分解 -----------------
calc_podani_jaccard_components <- function(comm_matrix) {
  x <- as.matrix(comm_matrix); x[x > 0] <- 1
  A <- x %*% t(x); R <- rowSums(x)
  pairs_idx <- which(lower.tri(A), arr.ind = TRUE)
  j_col <- pairs_idx[, 2]; k_col <- pairs_idx[, 1]
  a_vec <- A[lower.tri(A)]; b_vec <- R[j_col] - a_vec; c_vec <- R[k_col] - a_vec
  denom <- a_vec + b_vec + c_vec; denom[denom == 0] <- 1
  data.frame(Similarity = a_vec / denom, Replacement = (2 * pmin(b_vec, c_vec)) / denom, RichDiff = abs(b_vec - c_vec) / denom)
}

calc_podani_ruzicka_components <- function(raw_counts_matrix) {
  x <- as.matrix(raw_counts_matrix); n <- nrow(x)
  pairs_idx <- which(lower.tri(matrix(0, n, n)), arr.ind = TRUE)
  j_col <- pairs_idx[, 2]; k_col <- pairs_idx[, 1]
  A_vec <- sapply(1:nrow(pairs_idx), function(i) sum(pmin(x[pairs_idx[i, 2], ], x[pairs_idx[i, 1], ])))
  R_sums <- rowSums(x); B_vec <- R_sums[j_col] - A_vec; C_vec <- R_sums[k_col] - A_vec
  denom <- A_vec + B_vec + C_vec; denom[denom == 0] <- 1
  data.frame(Similarity = A_vec / denom, Replacement = (2 * pmin(B_vec, C_vec)) / denom, RichDiff = abs(B_vec - C_vec) / denom)
}

get_partition_df <- function(cnt_matrix, kingdom_name) {
  unweighted_comp <- calc_podani_jaccard_components(cnt_matrix)
  df_unweighted <- data.frame(Kingdom = kingdom_name, Weighting = "Unweighted", Similarity = unweighted_comp$Similarity, Replacement = unweighted_comp$Replacement, RichDiff = unweighted_comp$RichDiff)
  weighted_comp <- calc_podani_ruzicka_components(cnt_matrix)
  df_weighted <- data.frame(Kingdom = kingdom_name, Weighting = "Weighted", Similarity = weighted_comp$Similarity, Replacement = weighted_comp$Replacement, RichDiff = weighted_comp$RichDiff)
  bind_rows(df_unweighted, df_weighted)
}

part_bac  <- get_partition_df(bac_sub_cnt, "Bacteria")
part_fun  <- get_partition_df(fun_sub_cnt, "Fungi")
part_prot <- get_partition_df(prot_sub_cnt, "Protist")
part_all  <- bind_rows(part_bac, part_fun, part_prot)
part_box_df <- part_all %>% pivot_longer(cols = c(Replacement, RichDiff), names_to = "Component", values_to = "Value")

plot_ternary_kingdom <- function(df_kingdom, kingdom_name, color_val) {
  pts <- df_kingdom %>% filter(Weighting == "Unweighted") %>% mutate(x = Similarity + 0.5 * Replacement, y = 0.866 * Replacement)
  centroid <- pts %>% summarise(S_c = mean(Similarity), R_c = mean(Replacement), D_c = mean(RichDiff))
  cx <- centroid$S_c + 0.5 * centroid$R_c; cy = 0.866 * centroid$R_c
  triangle_lines <- data.frame(x = c(0, 1, 0.5, 0), y = c(0, 0, 0.866, 0))
  
  ggplot() +
    geom_path(data = triangle_lines, aes(x=x, y=y), color="#5D6D7E", linewidth=0.40) +
    ggrastr::rasterise(geom_point(data = pts, aes(x=x, y=y), color=color_val, size=0.4, alpha=0.30), dpi = 600) +
    geom_point(aes(x=cx, y=cy), fill="white", color="#2C3E50", shape=21, size=1.8, stroke=0.8) +
    annotate("text", x = 0.5, y = -0.16, label = "Similarity (%)", size = 1.9, fontface = "bold") +
    coord_fixed(ratio = 0.75, xlim = c(-0.20, 1.20), ylim = c(-0.22, 0.98)) +
    labs(title = kingdom_name) + theme_void() + theme(plot.title = element_text(face = "bold", size = 7.8, hjust = 0.5))
}

p_ternary_bac  <- plot_ternary_kingdom(part_bac, "Bacteria", kingdom_cols["Bacteria"])
p_ternary_fun  <- plot_ternary_kingdom(part_fun, "Fungi", kingdom_cols["Fungi"])
p_ternary_prot <- plot_ternary_kingdom(part_prot, "Protist", kingdom_cols["Protist"])

p_box_l1 <- ggplot(part_box_df, aes(x = Component, y = Value, fill = Component)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, alpha = 0.70, color = "black", linewidth = 0.35) +
  facet_grid(Weighting ~ Kingdom, scales = "free_y") +
  scale_fill_manual(values = component_cols) + labs(x = NULL, y = "Partitioning Beta Diversity") + theme_nature_v44_grid()

p_box_l2 <- ggplot(part_box_df, aes(x = Kingdom, y = Value, fill = Kingdom)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, alpha = 0.70, color = "black", linewidth = 0.35) +
  facet_grid(Weighting ~ Component, scales = "free_y") +
  scale_fill_manual(values = kingdom_cols) + labs(x = NULL, y = "Partitioning Beta Diversity") + theme_nature_v44_grid()

# ----------------- 6. Alpha 多样性网格 (Panel E) -----------------
alpha_treatment <- env_aligned %>%
  dplyr::select(SampleID, Urban, Veg, Season, Bac_Chao1, Fun_Chao1, Prot_Chao1) %>%
  pivot_longer(cols = c(Bac_Chao1, Fun_Chao1, Prot_Chao1), names_to = "Kingdom", values_to = "Chao1") %>%
  mutate(Kingdom = case_match(Kingdom, "Bac_Chao1" ~ "Bacteria", "Fun_Chao1" ~ "Fungi", "Prot_Chao1" ~ "Protist"))

alpha_grouped_trends <- alpha_treatment %>%
  group_by(Kingdom, Season, Urban, Veg) %>%
  summarise(Mean = mean(Chao1, na.rm = TRUE), SE = sd(Chao1, na.rm = TRUE)/sqrt(n()), .groups = 'drop')

p_alpha_grid <- ggplot(alpha_grouped_trends, aes(x = Urban, y = Mean, color = Veg, group = Veg, shape = Veg, linetype = Veg)) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.22, linewidth = 0.32, position = position_dodge(0.65)) +
  geom_line(linewidth = 0.45, position = position_dodge(0.65)) +
  geom_point(size = 1.3, stroke = 0.5, fill = "white", position = position_dodge(0.65)) + 
  facet_grid(Kingdom ~ Season, scales = "free_y") +
  scale_color_manual(values = veg_cols, name = "Vegetation") +
  scale_shape_manual(values = veg_shapes, name = "Vegetation") +
  labs(x = "Urbanization Level", y = "Chao1 Richness (Mean ± SE)") + theme_nature_v44_grid()

# ----------------- 7. 组装并导出主图 Figure 3 -----------------
row1_composite <- ((p_donut_bac + labs(tag = "A")) | (p_donut_fun + labs(tag = "B")) | (p_donut_prot + labs(tag = "C")) | (p_ddr_bc_standard + labs(tag = "D"))) + plot_layout(widths = c(1.05, 1.0, 1.05, 1.35))
row2_alpha <- p_alpha_grid + labs(tag = "E")
row3_ternary <- ((p_ternary_bac + labs(tag = "F")) | (p_ternary_fun + labs(tag = "G")) | (p_ternary_prot + labs(tag = "H"))) + plot_layout(ncol = 3)
row4_boxes <- ((p_box_l1 + labs(tag = "I")) | (p_box_l2 + labs(tag = "J"))) + plot_layout(ncol = 2, widths = c(1, 1))

fig3_composite <- (row1_composite / row2_alpha / row3_ternary / row4_boxes) + 
  plot_layout(heights = c(0.95, 1.30, 0.95, 0.85)) & 
  theme(plot.tag = element_text(face = "bold", size = 11.0))

ggsave(file.path(output_dir, "Figure3_Diversity_Composite.pdf"), plot = fig3_composite, width = 19, height = 20, units = "cm", device = cairo_pdf)
ggsave(file.path(output_dir, "Figure3_Diversity_Composite.png"), plot = fig3_composite, width = 19, height = 20, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "Figure3_Diversity_Composite.svg"), plot = fig3_composite, width = 19, height = 20, units = "cm", device = "svg")

# ----------------- 8. 附图 Figure S3 构建 (Mantel, PCoA, 3x4 Alpha) -----------------
filter_top_taxa <- function(otu_matrix, top_n = 500) {
  mean_abund <- colMeans(otu_matrix)
  keep_cols <- names(sort(mean_abund, decreasing = TRUE))[1:min(top_n, ncol(otu_matrix))]
  return(otu_matrix[, keep_cols, drop = FALSE])
}

bac_fast  <- filter_top_taxa(bac_sub_rel, 500)
fun_fast  <- filter_top_taxa(fun_sub_rel, 500)
prot_fast <- filter_top_taxa(prot_sub_rel, 500)

combined_communities_fast <- cbind(bac_fast, fun_fast, prot_fast)
ncol_bac_fast  <- ncol(bac_fast)
ncol_fun_fast  <- ncol(fun_fast)
ncol_prot_fast <- ncol(prot_fast)

spec_groups_fast <- list(
  `Bacterial community` = 1:ncol_bac_fast,
  `Fungal community`    = (ncol_bac_fast + 1):(ncol_bac_fast + ncol_fun_fast),
  `Protist community`   = (ncol_bac_fast + ncol_fun_fast + 1):(ncol_bac_fast + ncol_fun_fast + ncol_prot_fast)
)

env_vars_only <- env_aligned %>% select(all_of(target_env_vars))
colnames(env_vars_only) <- gsub("^NH4_N$", "NH\u2084\u207A-N", colnames(env_vars_only))
colnames(env_vars_only) <- gsub("^swc$", "SWC", colnames(env_vars_only))

mantel_raw_results <- mantel_test(
  spec = combined_communities_fast, env = env_vars_only,
  spec_select = spec_groups_fast, mantel_fun = "mantel"
)

mantel_formatted <- mantel_raw_results %>%
  mutate(
    rd = cut(r, breaks = c(-Inf, 0.4, Inf), labels = c("< 0.4", ">= 0.4")),
    pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
  )

cor_env <- correlate(env_vars_only)

p_mantel_network <- qcorrplot(cor_env, type = "upper", diag = FALSE, fixed = FALSE) + 
  geom_square() +
  geom_couple(
    aes(colour = pd, size = rd),
    data = mantel_formatted,
    curvature = nice_curvature(0.15)
  ) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdYlBu"), limits = c(-1, 1), name = "Pearson's r") +
  scale_size_manual(values = c(0.45, 1.2), name = "Mantel's r") +
  scale_colour_manual(values = c("#BC3C29", "#20854E", "#BDC3C7"), name = "Mantel's p") +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.2, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 6.2, face = "bold", colour = "black"),
    legend.position = "right"
  )

get_pcoa_plot <- function(dist_matrix, env_df, kingdom_name, color_palette) {
  pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  var_explained <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)
  
  perm_res <- vegan::adonis2(dist_matrix ~ Urban, data = env_df, permutations = 999)
  r2_val <- round(perm_res$R2[1], 3)
  p_val <- perm_res$`Pr(>F)`[1]
  p_lab <- if (p_val <= 0.001) "P < 0.001" else paste0("P = ", round(p_val, 3))
  stat_label <- paste0("PERMANOVA:\nR² = ", sprintf("%.3f", r2_val), "\n", p_lab)
  
  df_pcoa <- data.frame(
    SampleID = rownames(pcoa_res$points),
    PCoA1 = pcoa_res$points[, 1],
    PCoA2 = pcoa_res$points[, 2],
    stringsAsFactors = FALSE
  ) %>% left_join(env_df, by = "SampleID")
  
  ggplot(df_pcoa, aes(x = PCoA1, y = PCoA2, color = Urban, shape = Veg)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.25) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.25) +
    stat_ellipse(aes(group = Urban, color = Urban, fill = Urban), geom = "polygon", level = 0.95, alpha = 0.06, linewidth = 0.35) +
    geom_point(size = 1.0, alpha = 0.85, stroke = 0.35, fill = "white") +
    annotate("text", x = -Inf, y = Inf, label = stat_label, hjust = -0.08, vjust = 1.15,
             fontface = "bold.italic", size = 2.0, family = "sans", color = "black") +
    scale_color_manual(values = color_palette, name = "Urbanization") +
    scale_fill_manual(values = color_palette, name = "Urbanization") +
    scale_shape_manual(values = veg_shapes, name = "Vegetation") +
    labs(x = paste0("PCoA1 (", var_explained[1], "%)"), y = paste0("PCoA2 (", var_explained[2], "%)"), title = paste(kingdom_name, "communities")) +
    theme_nature_v44_grid()
}

p_pcoa_bac  <- get_pcoa_plot(dist_bac, env_aligned, "Bacterial", report_cols)
p_pcoa_fun  <- get_pcoa_plot(dist_fun, env_aligned, "Fungal", report_cols)
p_pcoa_prot <- get_pcoa_plot(dist_prot, env_aligned, "Protist", report_cols)

p_mantel_full  <- p_mantel_network + labs(tag = "A")
p_pcoa_row     <- ((p_pcoa_bac + labs(tag = "B")) | (p_pcoa_fun + labs(tag = "C")) | (p_pcoa_prot + labs(tag = "D"))) + plot_layout(ncol = 3, guides = "collect")

fig_s3_composite <- (p_mantel_full / p_pcoa_row) + plot_layout(heights = c(1.3, 1.0)) &
  theme(plot.tag = element_text(face = "bold", size = 11.0))

ggsave(file.path(output_dir, "FigureS3_Ecology_Supplementary_Atlas.pdf"), plot = fig_s3_composite, width = 19, height = 20, units = "cm", device = cairo_pdf)
ggsave(file.path(output_dir, "FigureS3_Ecology_Supplementary_Atlas.png"), plot = fig_s3_composite, width = 19, height = 20, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "FigureS3_Ecology_Supplementary_Atlas.svg"), plot = fig_s3_composite, width = 19, height = 20, units = "cm", device = "svg")

# ----------------- 9. 导出附表 Table S6 -----------------
table_s6_ddr <- ddr_all_modes %>%
  group_by(Kingdom, Distance) %>%
  do({
    fit <- lm(Comm_Sim ~ Env_Sim_Scaled, data = .)
    sum_fit <- summary(fit)
    data.frame(
      Slope     = sprintf("%.4f", coef(fit)[2]),
      Intercept = sprintf("%.4f", coef(fit)[1]),
      R_squared = sprintf("%.4f", sum_fit$r.squared),
      P_value   = if (sum_fit$coefficients[2, 4] < 0.001) "< 0.001" else sprintf("%.4f", sum_fit$coefficients[2, 4]),
      stringsAsFactors = FALSE
    )
  }) %>% ungroup() %>%
  rename(
    `Taxonomic Kingdom` = Kingdom,
    `Similarity Metric` = Distance,
    `R²` = R_squared,
    `P value` = P_value
  )

safe_write_csv(table_s6_ddr, file.path(output_dir, "TableS6_Distance_Decay_Results.csv"))

cat("\n=======================================================================\n")
cat(">>> [Script 02 Completed Successfully!] <<<\n")
cat(">>> Output Directory: ", output_dir, "\n")
cat(">>> 1. Figure 3: Figure3_Diversity_Composite (PDF/PNG/SVG, 19cm x 20cm)\n")
cat(">>> 2. Figure S3: FigureS3_Ecology_Supplementary_Atlas (PDF/PNG/SVG)\n")
cat(">>> 3. Table S6: TableS6_Distance_Decay_Results.csv\n")
cat("=======================================================================\n")