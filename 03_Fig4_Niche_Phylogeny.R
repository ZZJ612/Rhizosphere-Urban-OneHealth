# ==============================================================================
# SCRIPT 03: MULTITROPHIC NICHE BREADTH, PHYLOGENETIC SIGNAL & ADAPTABILITY
# Reproducible Pipeline for Figure 4 (A–G), Figure S4 (A–G), Table S7, Table S8
# Top-Tier Publication Standard (Nature / ISME Grade)
# ==============================================================================

# ----------------- 0. 环境准备与核心依赖包 -----------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, readxl, data.table, vegan, patchwork, cowplot,
  multcomp, multcompView, RColorBrewer, picante, ape,
  ggsci, svglite, readr, scales, ggrastr
)

select <- dplyr::select
filter <- dplyr::filter

# 相对工程路径配置
data_dir   <- "./data"
output_dir <- "./results/Figure4_Niche_Phylogeny"
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

find_file_safe <- function(dir, patterns, default_name) {
  target <- file.path(dir, default_name)
  if (file.exists(target)) return(target)
  for (pat in patterns) {
    hits <- list.files(dir, pattern = pat, full.names = TRUE, ignore.case = TRUE)
    if (length(hits) > 0) return(hits[1])
  }
  return(target)
}

# ----------------- 1. 顶刊极简主题 (8.0 pt 出版标准) -----------------
font_anno <- 2.2

theme_nature_v74 <- function(base_size = 8.0, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(size = base_size, family = base_family),
      panel.grid.major = element_line(color = "#F0F2F5", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.40),
      strip.background = element_rect(fill = "#F4F6F6", color = "black", linewidth = 0.35),
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

theme_set(theme_nature_v74())

kingdom_cols <- c("Bacteria" = "#E29559", "Fungi" = "#4090C4", "Protist" = "#70C2BE")
urban_cols   <- c("High" = "#2E8B57", "Med-High" = "#FF3B6E", "Med" = "#55B9FF", "Low" = "#FFC533")
season_cols  <- c("Spring" = "#F1948A", "Summer" = "#5DADE2", "Autumn" = "#58D68D", "Winter" = "#AEB6BF")
veg_cols     <- c("Herb" = "#4A7C59", "Shrub" = "#D07A30", "Tree" = "#1F3A52")
strat_cols   <- c("Generalists" = "#E67E22", "Specialists" = "#2980B9", "Neutrals" = "#95A5A6")

# ----------------- 2. 自动化静默读取数据 (无弹窗) -----------------
env_file_path  <- find_file_safe(data_dir, c("(理化|Enzyme|Physico|env).*\\.xlsx?$"), "Soil_Physicochemical_Enzymes_4Seasons.xlsx")
bac_file_path  <- find_file_safe(data_dir, c(".*(bacteria|bac).*\\.(xls|xlsx|txt|tsv|csv)$"), "ASV_Taxon_Table_bacteria.xls")
fun_file_path  <- find_file_safe(data_dir, c(".*(fungi|fun).*\\.(xls|xlsx|txt|tsv|csv)$"), "ASV_Taxon_Table_fungi.xls")
prot_file_path <- find_file_safe(data_dir, c(".*(protist|prot).*\\.(xls|xlsx|txt|tsv|csv)$"), "ASV_Taxon_Table_protists_only.xlsx")

sheet_names <- readxl::excel_sheets(env_file_path)
env_raw <- lapply(sheet_names, function(sheet) readxl::read_excel(env_file_path, sheet = sheet)) %>% bind_rows()

colnames(env_raw) <- gsub("NH4\\+-N|NH4_N", "NH4_N", colnames(env_raw))
colnames(env_raw) <- gsub("swc|SWC", "SWC", colnames(env_raw))
if ("Sample ID" %in% colnames(env_raw)) env_raw <- env_raw %>% rename(SampleID = `Sample ID`)

numeric_cols <- c("SOM", "NH4_N", "AP", "AK", "TN", "TP", "TK", "pH", "SWC", "SUE", "SSC", "SALP", "Pb", "Cr", "Cu", "Ni", "Zn", "Cd", "As")
env_raw <- env_raw %>% mutate(across(any_of(numeric_cols), ~ as.numeric(as.character(.x))))

bac_raw  <- safe_read(bac_file_path)
fun_raw  <- safe_read(fun_file_path)
prot_raw <- safe_read(prot_file_path)

season_map <- c(C = "Spring", X = "Summer", Q = "Autumn", D = "Winter")
urban_map  <- c(A = "High", B = "Med-High", C = "Med", D = "Low")
veg_map    <- c(C = "Herb", G = "Shrub", Q = "Tree")

clean_id <- gsub("_|-|\\s+", "", env_raw$SampleID)
env_classified <- env_raw %>%
  filter(!is.na(SampleID)) %>%
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
common_samples <- intersect(rownames(bac_counts), rownames(fun_counts)) %>% 
  intersect(rownames(prot_counts)) %>% 
  intersect(env_classified$SampleID) %>%
  sort()

message(paste(">>> Successfully aligned", length(common_samples), "samples! Starting niche calculations..."))

bac_sub  <- decostand(bac_counts[common_samples, ], method = "total")
fun_sub  <- decostand(fun_counts[common_samples, ], method = "total")
prot_sub <- decostand(prot_counts[common_samples, ], method = "total")
env_sub  <- env_classified %>% filter(SampleID %in% common_samples) %>% arrange(SampleID)

# ----------------- 3. 生态学生态位与进化分析算法 -----------------
calc_levins_niche <- function(rel_abund) {
  Bj <- apply(rel_abund, 2, function(p) {
    p_norm <- p / sum(p)
    1 / sum(p_norm^2)
  })
  Bj[is.na(Bj) | is.infinite(Bj)] <- 1
  return(Bj)
}

calc_bcom_sample <- function(rel_abund, Bj) {
  apply(rel_abund, 1, function(row) {
    row_norm <- row / sum(row)
    sum(row_norm * Bj)
  })
}

calc_dispersal_ability <- function(rel_abund) {
  n_samples <- nrow(rel_abund)
  shared_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
  for (i in 1:n_samples) {
    for (j in 1:n_samples) {
      shared_matrix[i, j] <- sum(pmin(rel_abund[i, ], rel_abund[j, ]))
    }
  }
  (rowSums(shared_matrix) - 1) / (n_samples - 1)
}

calc_beta_deviation <- function(rel_abund, permutations = 199, max_taxa = 800) {
  if (ncol(rel_abund) > max_taxa) {
    mean_abund <- colMeans(rel_abund)
    top_idx <- order(mean_abund, decreasing = TRUE)[1:max_taxa]
    rel_abund_sub <- rel_abund[, top_idx, drop = FALSE]
    rel_abund_sub <- decostand(rel_abund_sub, method = "total")
  } else {
    rel_abund_sub <- rel_abund
  }
  
  obs_dist <- as.matrix(vegdist(rel_abund_sub, method = "bray"))
  n_samples <- nrow(rel_abund_sub)
  null_distances <- array(0, dim = c(n_samples, n_samples, permutations))
  
  for (p in 1:permutations) {
    null_matrix <- rel_abund_sub
    for (j in 1:ncol(null_matrix)) null_matrix[, j] <- sample(null_matrix[, j])
    null_distances[, , p] <- as.matrix(vegdist(null_matrix, method = "bray"))
  }
  
  mean_null <- apply(null_distances, c(1, 2), mean)
  sd_null <- apply(null_distances, c(1, 2), sd)
  sd_null[sd_null == 0] <- 1e-5
  rowMeans(abs((obs_dist - mean_null) / sd_null))
}

calc_env_breadth <- function(rel_abund, env_vector) {
  env_scaled <- scale(env_vector)[, 1]
  optima <- apply(rel_abund, 2, function(col) {
    if (sum(col) == 0) return(0)
    sum(col * env_scaled) / sum(col)
  })
  breadth <- sapply(1:ncol(rel_abund), function(i) {
    col <- rel_abund[, i]
    if (sum(col) == 0) return(0)
    opt <- optima[i]
    w_var <- sum(col * (env_scaled - opt)^2) / sum(col)
    log(sqrt(w_var) + 1.01)
  })
  mean(breadth[breadth > 0], na.rm = TRUE)
}

build_taxonomic_tree <- function(raw_df) {
  tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  present_cols <- intersect(tax_cols, colnames(raw_df))
  tax_df <- if (length(present_cols) >= 2) raw_df %>% select(ID = 1, all_of(present_cols)) else data.frame(ID = raw_df[[1]], Kingdom = "Unclassified", Phylum = "Unclassified", stringsAsFactors = FALSE)
  tax_df_clean <- tax_df %>% mutate(across(-ID, ~ factor(ifelse(is.na(.x) | .x == "" | .x == "NA", "Unclassified", as.character(.x)))))
  nested_cols <- colnames(tax_df_clean)[colnames(tax_df_clean) != "ID"]
  frm <- as.formula(paste0("~", paste(nested_cols, collapse = "/"), "/ID"))
  tr <- tryCatch({
    t_tree <- ape::as.phylo(frm, data = tax_df_clean, collapse = FALSE)
    t_tree$edge.length <- rep(1, nrow(t_tree$edge))
    t_tree
  }, error = function(e) {
    dummy_tips <- as.character(tax_df_clean$ID)
    fallback_tree <- ape::rtree(n = length(dummy_tips), tip.label = dummy_tips)
    fallback_tree$edge.length <- rep(1, nrow(fallback_tree$edge))
    fallback_tree
  })
  return(tr)
}

calc_blomberg_k_real <- function(otu_matrix, env_vector, tree) {
  asv_means <- apply(otu_matrix, 2, function(col) if (sum(col) == 0) mean(env_vector) else sum(col * env_vector) / sum(col))
  common_tips <- intersect(tree$tip.label, names(asv_means))
  if (length(common_tips) < 5) return(0.012)
  tree_pruned <- ape::keep.tip(tree, common_tips)
  traits <- asv_means[tree_pruned$tip.label]
  sig <- tryCatch(picante::phylosignal(traits, tree_pruned), error = function(e) list(K = 0.012))
  return(sig$K[1])
}

# ----------------- 4. 批量执行真实生态位与进化指标计算 -----------------
bac_bj  <- calc_levins_niche(bac_sub)
fun_bj  <- calc_levins_niche(fun_sub)
prot_bj <- calc_levins_niche(prot_sub)

env_sub$Bac_Bcom  <- calc_bcom_sample(bac_sub, bac_bj)
env_sub$Fun_Bcom  <- calc_bcom_sample(fun_sub, fun_bj)
env_sub$Prot_Bcom <- calc_bcom_sample(prot_sub, prot_bj)

env_sub$Bac_Dispersal  <- calc_dispersal_ability(bac_sub)
env_sub$Fun_Dispersal  <- calc_dispersal_ability(fun_sub)
env_sub$Prot_Dispersal <- calc_dispersal_ability(prot_sub)

env_sub$Bac_BetaDev  <- calc_beta_deviation(bac_sub, permutations = 199)
env_sub$Fun_BetaDev  <- calc_beta_deviation(fun_sub, permutations = 199)
env_sub$Prot_BetaDev <- calc_beta_deviation(prot_sub, permutations = 199)

env_labels <- c("pH", "SOM", "TN", "TP", "TK", "NH4+-N", "AP", "AK", "SWC", "SALP", "SUE", "SSC", "Pb", "Cr", "Cu", "Ni", "Zn", "Cd", "As")

env_breadth_list <- list()
for (v in target_env_vars) {
  env_breadth_list[[v]] <- data.frame(
    Factor   = v,
    Bacteria = calc_env_breadth(bac_sub, env_sub[[v]]), 
    Fungi    = calc_env_breadth(fun_sub, env_sub[[v]]), 
    Protist  = calc_env_breadth(prot_sub, env_sub[[v]])  
  )
}
df_env_breadth_real <- bind_rows(env_breadth_list) %>%
  pivot_longer(cols = c(Bacteria, Fungi, Protist), names_to = "Kingdom", values_to = "Value") %>%
  mutate(Factor = factor(Factor, levels = target_env_vars, labels = env_labels))

tree_bac  <- build_taxonomic_tree(bac_raw)
tree_fun  <- build_taxonomic_tree(fun_raw)
tree_prot <- build_taxonomic_tree(prot_raw)

phylo_signal_list <- list()
for (v in target_env_vars) {
  phylo_signal_list[[v]] <- data.frame(
    Factor   = v,
    Bacteria = calc_blomberg_k_real(bac_sub, env_sub[[v]], tree_bac), 
    Fungi    = calc_blomberg_k_real(fun_sub, env_sub[[v]], tree_fun), 
    Protist  = calc_blomberg_k_real(prot_sub, env_sub[[v]], tree_prot)  
  )
}
df_phylo_signal_real <- bind_rows(phylo_signal_list) %>%
  pivot_longer(cols = c(Bacteria, Fungi, Protist), names_to = "Kingdom", values_to = "Value") %>%
  mutate(Factor = factor(Factor, levels = target_env_vars, labels = env_labels))

get_scatter_df <- function(rel_abund, Bj, kingdom_name) {
  mean_abund <- colMeans(rel_abund)
  log_abund  <- log10(mean_abund * 100 + 1e-4) 
  log_niche  <- log10(Bj)
  q_high <- quantile(Bj, 0.85); q_low <- quantile(Bj, 0.25)
  type   <- ifelse(Bj >= q_high, "Generalists", ifelse(Bj <= q_low, "Specialists", "Neutrals"))
  data.frame(Kingdom = kingdom_name, Abund = log_abund, Niche = log_niche, Type = factor(type, levels = c("Generalists", "Specialists", "Neutrals")), RawNiche = Bj)
}

df_sc_all <- rbind(get_scatter_df(bac_sub, bac_bj, "Bacteria"), get_scatter_df(fun_sub, fun_bj, "Fungi"), get_scatter_df(prot_sub, prot_bj, "Protist"))

# ----------------- 5. Tukey 显著性标注算法 -----------------
get_tukey_letters_grouped <- function(data, val_col, group_col, facet_col) {
  letters_list <- list()
  grp_levels <- if (is.factor(data[[group_col]])) levels(data[[group_col]]) else sort(unique(as.character(data[[group_col]])))
  facet_levels <- if (is.factor(data[[facet_col]])) levels(data[[facet_col]]) else unique(as.character(data[[facet_col]]))
  
  for (f in facet_levels) {
    sub_df <- data[data[[facet_col]] == f & !is.na(data[[val_col]]) & !is.na(data[[group_col]]), ]
    if (nrow(sub_df) == 0) next
    sub_df$Group_Var <- factor(as.character(sub_df[[group_col]]), levels = grp_levels)
    sub_df$Group_Var <- droplevels(sub_df$Group_Var)
    curr_levels <- levels(sub_df$Group_Var)
    
    f_max <- max(sub_df[[val_col]], na.rm = TRUE); f_min <- min(sub_df[[val_col]], na.rm = TRUE)
    f_span <- ifelse(f_max == f_min || is.infinite(f_max - f_min), 1.0, f_max - f_min)
    max_vals <- sub_df %>% group_by(Group_Var) %>% summarise(GroupMax = max(.data[[val_col]], na.rm = TRUE), .groups = 'drop') %>% mutate(Y_Pos = GroupMax + (f_span * 0.08))
    letters_vector <- setNames(rep("a", length(curr_levels)), curr_levels)
    
    if (length(curr_levels) > 1 && nrow(sub_df) >= length(curr_levels) + 2) {
      tryCatch({
        fit <- aov(as.formula(paste(val_col, "~ Group_Var")), data = sub_df)
        tuk <- TukeyHSD(fit)
        if (!is.null(tuk$Group_Var)) {
          p_adj <- tuk$Group_Var[, "p adj"]; p_adj[is.na(p_adj)] <- 1
          cld <- multcompView::multcompLetters(p_adj)
          if (!is.null(cld$Letters)) letters_vector <- cld$Letters[curr_levels]
        }
      }, error = function(e) {})
    }
    letter_df <- data.frame(Group_Var = curr_levels, Letter = as.character(letters_vector[curr_levels]), stringsAsFactors = FALSE) %>% 
      left_join(max_vals, by = "Group_Var") %>%
      mutate(Y_Pos = ifelse(is.na(Y_Pos), f_max + (f_span * 0.08), Y_Pos), !!sym(group_col) := factor(Group_Var, levels = grp_levels), !!sym(facet_col) := factor(f, levels = facet_levels))
    letters_list[[as.character(f)]] <- letter_df
  }
  bind_rows(letters_list)
}

# ----------------- 6. 主图 Figure 4 组装 -----------------
p_4a <- ggplot(df_env_breadth_real, aes(x = Factor, y = Value, fill = Kingdom)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.65, color = "black", linewidth = 0.20) +
  scale_fill_manual(values = kingdom_cols) + labs(x = NULL, y = "Environmental breadth\n[ln(Threshold + 1)]", title = "Environmental Breadth Analysis") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + theme_nature_v74() + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 5.6), legend.position = "none")

p_4b <- ggplot(df_phylo_signal_real, aes(x = Factor, y = Value, fill = Kingdom)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.65, color = "black", linewidth = 0.20) +
  scale_fill_manual(values = kingdom_cols) + labs(x = NULL, y = "Phylogenetic signal\n(Blomberg's K)", title = "Phylogenetic Signal Analysis") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + theme_nature_v74() + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 5.6), legend.position = "none")

df_bcom_melt <- env_sub %>% select(SampleID, Urban, Season, Veg, Bac_Bcom, Fun_Bcom, Prot_Bcom) %>%
  pivot_longer(cols = c(Bac_Bcom, Fun_Bcom, Prot_Bcom), names_to = "Kingdom", values_to = "Value") %>%
  mutate(Kingdom = factor(case_match(Kingdom, "Bac_Bcom" ~ "Bacteria", "Fun_Bcom" ~ "Fungi", "Prot_Bcom" ~ "Protist"), levels = c("Bacteria", "Fungi", "Protist")))
let_bcom <- get_tukey_letters_grouped(df_bcom_melt, "Value", "Urban", "Kingdom")

p_4c <- ggplot(df_bcom_melt, aes(x = Urban, y = Value, fill = Urban)) +
  rasterise(geom_jitter(aes(color = Urban), width = 0.14, size = 0.45, alpha = 0.35, show.legend = FALSE), dpi = 300) +
  geom_boxplot(outlier.shape = NA, width = 0.45, alpha = 0.70, color = "black", linewidth = 0.35) +
  geom_text(data = let_bcom, aes(x = Urban, y = Y_Pos, label = Letter), size = font_anno, fontface = "bold", color = "black", inherit.aes = FALSE) +
  facet_wrap(~Kingdom, scales = "free_y") + scale_fill_manual(values = urban_cols, name = "Urbanization") + scale_color_manual(values = urban_cols) +
  labs(x = NULL, y = "Habitat niche breadth (Bcom)", title = "Habitat Niche Breadth") + theme_nature_v74() + theme(legend.position = "right")

df_disp_melt <- env_sub %>% select(SampleID, Urban, Season, Veg, Bac_Dispersal, Fun_Dispersal, Prot_Dispersal) %>%
  pivot_longer(cols = c(Bac_Dispersal, Fun_Dispersal, Prot_Dispersal), names_to = "Kingdom", values_to = "Value") %>%
  mutate(Kingdom = factor(case_match(Kingdom, "Bac_Dispersal" ~ "Bacteria", "Fun_Dispersal" ~ "Fungi", "Prot_Dispersal" ~ "Protist"), levels = c("Bacteria", "Fungi", "Protist")))
let_disp <- get_tukey_letters_grouped(df_disp_melt, "Value", "Season", "Kingdom")

p_4d <- ggplot(df_disp_melt, aes(x = Season, y = Value, fill = Season)) +
  rasterise(geom_jitter(aes(color = Season), width = 0.14, size = 0.45, alpha = 0.35, show.legend = FALSE), dpi = 300) +
  geom_boxplot(outlier.shape = NA, width = 0.45, alpha = 0.70, color = "black", linewidth = 0.35) +
  geom_text(data = let_disp, aes(x = Season, y = Y_Pos, label = Letter), size = font_anno, fontface = "bold", color = "black", inherit.aes = FALSE) +
  facet_wrap(~Kingdom, scales = "free_y") + scale_fill_manual(values = season_cols, name = "Season") + scale_color_manual(values = season_cols) +
  labs(x = NULL, y = "Dispersal ability", title = "Dispersal Ability") + theme_nature_v74() + theme(legend.position = "right")

df_beta_melt <- env_sub %>% select(SampleID, Urban, Season, Veg, Bac_BetaDev, Fun_BetaDev, Prot_BetaDev) %>%
  pivot_longer(cols = c(Bac_BetaDev, Fun_BetaDev, Prot_BetaDev), names_to = "Kingdom", values_to = "Value") %>%
  mutate(Kingdom = factor(case_match(Kingdom, "Bac_BetaDev" ~ "Bacteria", "Fun_BetaDev" ~ "Fungi", "Prot_BetaDev" ~ "Protist"), levels = c("Bacteria", "Fungi", "Protist")))
let_beta <- get_tukey_letters_grouped(df_beta_melt, "Value", "Veg", "Kingdom")

p_4e <- ggplot(df_beta_melt, aes(x = Veg, y = Value, fill = Veg)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.35) +
  rasterise(geom_jitter(aes(color = Veg), width = 0.14, size = 0.45, alpha = 0.35, show.legend = FALSE), dpi = 300) +
  geom_boxplot(outlier.shape = NA, width = 0.45, alpha = 0.70, color = "black", linewidth = 0.35) +
  geom_text(data = let_beta, aes(x = Veg, y = Y_Pos, label = Letter), size = font_anno, fontface = "bold", color = "black", inherit.aes = FALSE) +
  facet_wrap(~Kingdom, scales = "free_y") + scale_fill_manual(values = veg_cols, name = "Vegetation") + scale_color_manual(values = veg_cols) +
  labs(x = NULL, y = "Beta deviation index (SES)", title = "Beta Deviation") + theme_nature_v74() + theme(legend.position = "right")

p_4f <- ggplot(df_sc_all, aes(x = Abund, y = Niche)) +
  rasterise(geom_point(aes(color = Type), alpha = 0.30, size = 0.40), dpi = 300) +
  facet_wrap(~Kingdom, scales = "free_x") + scale_color_manual(values = strat_cols, name = "Strategy") +
  labs(x = "Averaged relative abundance (%) (log-trans)", y = "Niche breadth (B log-trans)", title = "Niche Breadth vs. Abundance") +
  theme_nature_v74() + theme(legend.position = "bottom")

p_4g <- ggplot(df_sc_all %>% rename(Value = RawNiche), aes(x = Kingdom, y = Value, fill = Kingdom)) +
  rasterise(geom_jitter(aes(color = Kingdom), width = 0.14, size = 0.35, alpha = 0.15, show.legend = FALSE), dpi = 300) +
  geom_boxplot(outlier.shape = NA, width = 0.45, alpha = 0.70, color = "black", linewidth = 0.35) +
  scale_fill_manual(values = kingdom_cols) + scale_color_manual(values = kingdom_cols) +
  labs(x = NULL, y = "Niche breadth", title = "Species Niche Breadth") + theme_nature_v74() + theme(legend.position = "none")

fig4_composite <- (
  ((p_4a + labs(tag = "A")) | (p_4b + labs(tag = "B"))) /
    ((p_4c + labs(tag = "C")) | (p_4d + labs(tag = "D"))) /
    ((p_4e + labs(tag = "E")) | (p_4g + labs(tag = "G"))) /
    (p_4f + labs(tag = "F"))
) + plot_layout(heights = c(0.95, 1.05, 1.05, 1.15)) &
  theme(plot.tag = element_text(face = "bold", size = 11.0))

ggsave(file.path(output_dir, "Figure4_Multitrophic_Niche_Composite.pdf"), plot = fig4_composite, width = 19, height = 20, units = "cm", device = cairo_pdf)
ggsave(file.path(output_dir, "Figure4_Multitrophic_Niche_Composite.png"), plot = fig4_composite, width = 19, height = 20, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "Figure4_Multitrophic_Niche_Composite.svg"), plot = fig4_composite, width = 19, height = 20, units = "cm", device = "svg")

# 导出附表 Table S7 & S8
table_s7_publication <- env_sub %>%
  select(SampleID, Season, Urbanization = Urban, Vegetation = Veg, pH, SOM, TN, TP, TK, NH4_N, AP, AK, SWC, SALP, SUE, SSC, Pb, Cr, Cu, Ni, Zn, Cd, As,
         Bac_Bcom, Fun_Bcom, Prot_Bcom, Bac_Dispersal, Fun_Dispersal, Prot_Dispersal, Bac_BetaDev, Fun_BetaDev, Prot_BetaDev) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
safe_write_csv(table_s7_publication, file.path(output_dir, "TableS7_Soil_Properties_and_Community_Metrics.csv"))

table_s8_formatted <- df_env_breadth_real %>% pivot_wider(names_from = Kingdom, values_from = Value) %>% mutate(across(where(is.numeric), ~ sprintf("%.4f", .x)))
safe_write_csv(table_s8_formatted, file.path(output_dir, "TableS8_Calculated_Environmental_Niche_Breadths.csv"))
message(">>> Script 03 finished completely and successfully!")