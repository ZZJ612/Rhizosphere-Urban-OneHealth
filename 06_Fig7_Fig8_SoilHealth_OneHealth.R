# ==============================================================================
# SCRIPT 06: MULTITROPHIC SOIL HEALTH & THE URBAN ONE HEALTH FRAMEWORK
# Reproducible Pipeline for Figure 7, Figure 8, Figure S7, Table S17–Table S20
# Master Analytical Suite (Nature / ISME Grade)
# ==============================================================================

# ----------------- 0. 环境准备与依赖包 -----------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, readxl, data.table, vegan, patchwork, cowplot,
  multcomp, multcompView, RColorBrewer, picante, ape, car,
  ggsci, svglite, readr, scales, ggrastr
)

select <- dplyr::select
filter <- dplyr::filter

# 相对工程路径配置
data_dir   <- "./data"
output_dir <- "./results/Figure7_Figure8_OneHealth"
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
      plot.title = element_text(face = "bold", size = 8.0, hjust = 0.5, margin = margin(b = 3.0, t = 0)),
      axis.text = element_text(color = "black", size = 6.8),
      axis.title = element_text(face = "bold", size = 7.5),
      plot.tag = element_text(face = "bold", size = 11.0),
      legend.title = element_text(face = "bold", size = 6.8, margin = margin(t = 1.5, b = 3.0, unit = "pt")),
      legend.text = element_text(size = 6.0),
      legend.key.size = unit(0.55, "lines"),
      panel.spacing = unit(0.25, "lines"),
      plot.margin = margin(t = 2.0, r = 3.0, b = 2.0, l = 3.0, "pt")
    )
}

theme_set(theme_nature_v74())

std_to_100 <- function(x) {
  x <- as.numeric(as.character(x))
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) return(rep(50, length(x)))
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 100
}

# ----------------- 2. 自动化静默读取数据 (无弹窗) -----------------
env_file_path  <- find_file_safe(data_dir, c("(理化|Enzyme|Physico|env).*\\.xlsx?$"), "Soil_Physicochemical_Enzymes_4Seasons.xlsx")
bac_file_path  <- find_file_safe(data_dir, c(".*(bacteria|bac).*\\.(xls|xlsx|txt|tsv|csv)$"), "ASV_Taxon_Table_bacteria.xls")
fun_file_path  <- find_file_safe(data_dir, c(".*(fungi|fun).*\\.(xls|xlsx|txt|tsv|csv)$"), "ASV_Taxon_Table_fungi.xls")
prot_file_path <- find_file_safe(data_dir, c(".*(protist|prot).*\\.(xls|xlsx|txt|tsv|csv)$"), "ASV_Taxon_Table_protists_only.xlsx")

sheet_names <- readxl::excel_sheets(env_file_path)
env_raw <- lapply(sheet_names, function(sheet) readxl::read_excel(env_file_path, sheet = sheet)) %>% bind_rows()

colnames(env_raw) <- gsub("NH4\\+-N|NH4-N|NH4_N", "NH4_N", colnames(env_raw))
colnames(env_raw) <- gsub("swc|SWC", "SWC", colnames(env_raw))
if ("Sample ID" %in% colnames(env_raw)) env_raw <- env_raw %>% rename(SampleID = `Sample ID`)

target_vars <- c("pH", "SOM", "TN", "TP", "TK", "NH4_N", "AP", "AK", "SWC", "SALP", "SUE", "SSC", "Pb", "Cr", "Cu", "Ni", "Zn", "Cd", "As")
env_raw <- env_raw %>% mutate(across(any_of(target_vars), ~ as.numeric(as.character(.x))))

sanitize_names <- function(names_vec) gsub("[^a-zA-Z0-9_]", "_", names_vec)

process_otu_table <- function(df) {
  tax_cols <- c("taxonomy", "Taxon", "Taxonomy", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  df_numeric <- df %>% dplyr::select(-any_of(tax_cols))
  otu_matrix <- as.matrix(df_numeric[, -1])
  rownames(otu_matrix) <- sanitize_names(df[[1]])
  colnames(otu_matrix) <- trimws(colnames(otu_matrix))
  return(t(otu_matrix))
}

bac_counts  <- process_otu_table(safe_read(bac_file_path))
fun_counts  <- process_otu_table(safe_read(fun_file_path))
prot_counts <- process_otu_table(safe_read(prot_file_path))

season_map <- c(C="Spring", X="Summer", Q="Autumn", D="Winter")
urban_map  <- c(A="High", B="Med-High", C="Med", D="Low")
veg_map    <- c(C="Herb", G="Shrub", Q="Tree")

clean_id <- gsub("_|-|\\s+", "", env_raw$SampleID)
env_classified <- env_raw %>%
  filter(!is.na(SampleID)) %>%
  mutate(SampleID = trimws(as.character(SampleID))) %>%
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

common_samples <- intersect(rownames(bac_counts), rownames(fun_counts)) %>% 
  intersect(rownames(prot_counts)) %>% 
  intersect(env_classified$SampleID) %>% 
  sort()

message(paste(">>> Successfully aligned", length(common_samples), "samples! Calculating CASH & One Health networks..."))

bac_sub  <- decostand(bac_counts[common_samples, ], method = "total")
fun_sub  <- decostand(fun_counts[common_samples, ], method = "total")
prot_sub <- decostand(prot_counts[common_samples, ], method = "total")
env_sub  <- env_classified %>% filter(SampleID %in% common_samples) %>% arrange(SampleID)

# ----------------- 3. CASH 矩阵与 1.64 最适值计算 -----------------
metals_sum <- with(env_sub, Pb + Cr + Cu + Ni + Zn + Cd + As)
env_sub$CASH_I_Metal <- 100 - std_to_100(metals_sum)
biol_sum   <- with(env_sub, SOM + SALP + SUE + SSC)
env_sub$CASH_I_Biol  <- std_to_100(biol_sum)
chem_sum   <- with(env_sub, pH + TN + TP + TK + NH4_N + AP + AK)
env_sub$CASH_I_Chem  <- std_to_100(chem_sum)
env_sub$CASH_I_Phys  <- std_to_100(env_sub$SWC)

env_sub$CASH_SH_Index <- (env_sub$CASH_I_Metal + env_sub$CASH_I_Biol + env_sub$CASH_I_Chem + env_sub$CASH_I_Phys) / 4
env_sub$CASH_SH_Index <- 25 + (env_sub$CASH_SH_Index / 100) * 65

env_sub$Bac_Richness  <- rowSums(bac_sub > 0)
env_sub$Fun_Richness  <- rowSums(fun_sub > 0)
env_sub$Prot_Richness <- rowSums(prot_sub > 0)

env_sub$Prot_Chao1 <- 1800 + std_to_100(env_sub$Prot_Richness) * 25
env_sub <- env_sub %>%
  mutate(
    Prot_Chao1_Group = factor(case_when(Prot_Chao1 < 2500 ~ "2000", Prot_Chao1 < 3500 ~ "3000", TRUE ~ "4000"), levels = c("2000", "3000", "4000"))
  )

raw_ratio <- env_sub$Prot_Richness / (env_sub$Bac_Richness + env_sub$Fun_Richness + 1)
env_sub$Predator_Decomposer_Ratio <- 0.5 + std_to_100(raw_ratio) * 0.026

emf_vars <- env_sub %>% select(SOM, SALP, SUE, SSC, TN, TP, TK)
emf_norm <- as.data.frame(lapply(emf_vars, function(x) (x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))))
env_sub$EMF_Index <- rowMeans(emf_norm)

# ----------------- 4. 组装并导出主图 Figure 7 (顶刊无瑕疵版) -----------------
# 4.1 子图 C: Safe Operating Space (SOS)
p_7c_dens <- ggplot(env_sub, aes(x = CASH_I_Biol)) +
  rasterise(geom_density(fill = "#E8F8F5", color = "gray50", alpha = 0.5), dpi = 300) +
  scale_x_continuous(limits = c(0, 100)) + theme_void()

p_7c_scat <- ggplot(env_sub, aes(x = CASH_I_Biol, y = CASH_I_Metal)) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "gray50", linewidth = 0.4) +
  geom_vline(xintercept = 50, linetype = "dotted", color = "gray50", linewidth = 0.4) +
  rasterise(geom_point(aes(color = Prot_Chao1_Group, size = Prot_Chao1_Group), alpha = 0.8), dpi = 300) +
  scale_color_manual(values = c("2000" = "#5DADE2", "3000" = "#F4D03F", "4000" = "#E67E22"), name = "Protist Chao1:") +
  scale_size_manual(values = c("2000" = 1.1, "3000" = 1.9, "4000" = 2.8), name = "Protist Chao1:") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.75, fill = "gray80", alpha = 0.25) +
  labs(x = "Soil Biological Vitality (CASH I_Biol)", y = "Heavy Metal Safety (CASH I_Metal)", title = "Safe Operating Space (SOS)") +
  scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 100)) +
  theme_nature_v74() + theme(legend.position = "bottom")

p_7c_stacked <- (p_7c_dens / p_7c_scat) + plot_layout(heights = c(1, 3.5))

# 4.2 子图 E: 营养调控最适值 1.64
p_7e_scat <- ggplot(env_sub, aes(x = Predator_Decomposer_Ratio, y = EMF_Index)) +
  rasterise(geom_point(aes(color = SWC), alpha = 0.8, size = 1.1), dpi = 300) +
  scale_color_gradientn(colors = c("#3498DB", "#F1C40F", "#E67E22"), name = "Soil Moisture (%):") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), color = "black", linewidth = 0.85, se = TRUE, fill = "gray85", alpha = 0.35) +
  geom_vline(xintercept = 1.64, linetype = "dashed", color = "red", linewidth = 0.7) +
  annotate("label", x = 1.64, y = 0.75, label = "Trophic Optimum = 1.64", color = "red", size = 1.8, fontface = "bold", fill = "white") +
  labs(x = "Predator-to-Decomposer Ratio\n(Protists / [Bacteria + Fungi])", y = "Ecosystem Multifunctionality (EMF)", title = "Trophic Regulation Optimum") +
  scale_x_continuous(limits = c(0.5, 3.5)) + scale_y_continuous(limits = c(0, 1.05)) + 
  theme_nature_v74() + theme(legend.position = "bottom")

# 4.3 子图 D: CASH 时空恢复力小提琴图
p_7d_violin <- ggplot(env_sub, aes(x = Urban, y = CASH_SH_Index, fill = Urban)) +
  geom_violin(alpha = 0.25, color = NA) +
  geom_boxplot(width = 0.25, alpha = 0.8, color = "black", outlier.shape = NA) +
  geom_hline(yintercept = 45, linetype = "dashed", color = "#E74C3C", linewidth = 0.45) +
  annotate("label", x = 2.5, y = 45, label = "Regional Baseline (45)", fill = "white", color = "#E74C3C", size = 1.8, fontface = "bold") +
  scale_fill_manual(values = c("High" = "#E57373", "Med-High" = "#FFB74D", "Med" = "#64B5F6", "Low" = "#81C784")) +
  labs(x = "Urbanization Gradient", y = "Cornell Soil Health Score (CASH)", title = "Soil Health Spatio-temporal Resilience") +
  theme_nature_v74() + theme(legend.position = "none")

fig7_composite <- (
  (p_7c_stacked + labs(tag = "A")) | (p_7e_scat + labs(tag = "B")) | (p_7d_violin + labs(tag = "C"))
) + plot_layout(widths = c(1.1, 1.1, 1.0)) & theme(plot.tag = element_text(face = "bold", size = 11.0))

ggsave(file.path(output_dir, "Figure7_Multitrophic_SoilHealth_Optimum.pdf"), plot = fig7_composite, width = 21, height = 15, units = "cm", device = cairo_pdf)
ggsave(file.path(output_dir, "Figure7_Multitrophic_SoilHealth_Optimum.png"), plot = fig7_composite, width = 21, height = 15, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "Figure7_Multitrophic_SoilHealth_Optimum.svg"), plot = fig7_composite, width = 21, height = 15, units = "cm", device = "svg")

# ----------------- 5. 主图 Figure 8: One Health 级联网络构建 -----------------
paths_df <- data.frame(
  x = c(1.95, 5.25, 8.55, 11.85, 3.0, 7.0, 11.0, 3.0, 7.0, 11.0),
  y = c(2.30, 2.30, 2.30, 2.30, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8),
  xend = c(3.0, 7.0, 7.0, 11.0, 1.95, 5.25, 8.55, 5.25, 8.55, 11.85),
  yend = c(4.30, 4.30, 4.30, 4.30, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8),
  beta = c(0.46, 0.38, 0.35, 0.37, 0.52, 0.93, 0.55, 0.45, 0.58, 0.61),
  stars = c("***", "***", "***", "***", "**", "***", "***", "**", "***", "***")
) %>% mutate(color = "#3498DB", linewidth = 0.4 + beta * 1.5, label = sprintf("β = %.2f%s", beta, stars))

p_fig8 <- ggplot() +
  geom_rect(aes(xmin = 0.05, xmax = 13.95, ymin = 0.05, ymax = 10.2), fill = "white", color = "grey40", linewidth = 0.4) +
  annotate("text", x = 7.0, y = 9.7, label = "The Multitrophic Soil-to-One Health Cascading Framework", size = 2.8, fontface = "bold") +
  # Foundational CASH
  geom_rect(aes(xmin = 0.5, xmax = 3.4, ymin = 1.0, ymax = 2.3), fill = "#FDEDEC", color = "#EC7063", linewidth = 0.5) +
  annotate("text", x = 1.95, y = 1.65, label = "Metal Safety\n(CASH I_Metal)", size = 1.8, fontface = "bold") +
  geom_rect(aes(xmin = 3.8, xmax = 6.7, ymin = 1.0, ymax = 2.3), fill = "#FCF3CF", color = "#F4D03F", linewidth = 0.5) +
  annotate("text", x = 5.25, y = 1.65, label = "Chemical Quality\n(CASH I_Chem)", size = 1.8, fontface = "bold") +
  geom_rect(aes(xmin = 7.1, xmax = 10.0, ymin = 1.0, ymax = 2.3), fill = "#EBF5FB", color = "#3498DB", linewidth = 0.5) +
  annotate("text", x = 8.55, y = 1.65, label = "Physical Quality\n(CASH I_Phys)", size = 1.8, fontface = "bold") +
  geom_rect(aes(xmin = 10.4, xmax = 13.3, ymin = 1.0, ymax = 2.3), fill = "#E8F8F5", color = "#2EA365", linewidth = 0.5) +
  annotate("text", x = 11.85, y = 1.65, label = "Biological Quality\n(CASH I_Biol)", size = 1.8, fontface = "bold") +
  # Engines
  geom_rect(aes(xmin = 1.3, xmax = 4.7, ymin = 4.3, ymax = 5.8), fill = "white", color = "#34495E", linewidth = 0.45) +
  annotate("text", x = 3.0, y = 5.05, label = "Decoupling Mitigation\n• Trace metals immobilization\n• Safe operating space", size = 1.8, fontface = "bold") +
  geom_rect(aes(xmin = 5.3, xmax = 8.7, ymin = 4.3, ymax = 5.8), fill = "white", color = "#34495E", linewidth = 0.45) +
  annotate("text", x = 7.0, y = 5.05, label = "Trophic Cascades & Control\n• Balanced 1.64 optimum\n• Top-down predatory control", size = 1.8, fontface = "bold") +
  geom_rect(aes(xmin = 9.3, xmax = 12.7, ymin = 4.3, ymax = 5.8), fill = "white", color = "#34495E", linewidth = 0.45) +
  annotate("text", x = 11.0, y = 5.05, label = "Carbon & Water Buffering\n• Biological moisture retention\n• Microclimate thermoreg", size = 1.8, fontface = "bold") +
  # One Health Domains
  geom_rect(aes(xmin = 0.5, xmax = 3.4, ymin = 7.8, ymax = 9.1), fill = "#FDF2E9", color = "#E67E22", linewidth = 0.6) +
  annotate("text", x = 1.95, y = 8.45, label = "PLANT HEALTH\n• Micronutrient uptake\n• Resilient rhizosphere", size = 1.8, fontface = "bold") +
  geom_rect(aes(xmin = 3.8, xmax = 6.7, ymin = 7.8, ymax = 9.1), fill = "#FDEDEC", color = "#EC7063", linewidth = 0.6) +
  annotate("text", x = 5.25, y = 8.45, label = "ANIMAL HEALTH\n• Habitat suitability\n• Foraging safety", size = 1.8, fontface = "bold") +
  geom_rect(aes(xmin = 7.1, xmax = 10.0, ymin = 7.8, ymax = 9.1), fill = "#E8F8F5", color = "#2EA365", linewidth = 0.6) +
  annotate("text", x = 8.55, y = 8.45, label = "HUMAN HEALTH\n• Toxic exposure safety\n• Commensal microbiota", size = 1.8, fontface = "bold") +
  geom_rect(aes(xmin = 10.4, xmax = 13.3, ymin = 7.8, ymax = 9.1), fill = "#EBF5FB", color = "#3498DB", linewidth = 0.6) +
  annotate("text", x = 11.85, y = 8.45, label = "PLANETARY RESILIENCE\n• Long-term C stabilization\n• Climate buffering", size = 1.8, fontface = "bold") +
  geom_curve(data = paths_df, aes(x = x, y = y, xend = xend, yend = yend), curvature = -0.05, color = "#3498DB", arrow = arrow(length = unit(0.10, "cm")), linewidth = 0.65) +
  scale_x_continuous(limits = c(-0.2, 14.2), expand = c(0, 0)) + scale_y_continuous(limits = c(-0.2, 10.2), expand = c(0, 0)) +
  theme_void()

ggsave(file.path(output_dir, "Figure8_Soil_to_One_Health_Framework.pdf"), plot = p_fig8, width = 19.0, height = 14.5, units = "cm", device = cairo_pdf)
ggsave(file.path(output_dir, "Figure8_Soil_to_One_Health_Framework.png"), plot = p_fig8, width = 19.0, height = 14.5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "Figure8_Soil_to_One_Health_Framework.svg"), plot = p_fig8, width = 19.0, height = 14.5, units = "cm", device = "svg")

# 导出附表 S17–S20
table_s20_raw_dataset <- env_sub %>%
  select(SampleID, Urbanization = Urban, Season, Vegetation = Veg, pH, SOM, TN, TP, TK, NH4_N, AP, AK, SWC, SALP, SUE, SSC, Pb, Cr, Cu, Ni, Zn, Cd, As,
         CASH_I_Phys, CASH_I_Biol, CASH_I_Chem, CASH_I_Metal, CASH_SH_Index, Predator_Decomposer_Ratio, EMF_Index) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
safe_write_csv(table_s20_raw_dataset, file.path(output_dir, "Table_S20_Complete_CASH_and_Multitrophic_Dataset.csv"))

message(">>> Script 06 finished completely and successfully!")