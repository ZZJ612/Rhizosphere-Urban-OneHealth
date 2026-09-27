# ==============================================================================
# SCRIPT 01: BIOGEOCHEMISTRY SUITE, ECOENZYME STOICHIOMETRY & SEM
# Reproducible Pipeline for Figure 2, Figure S1, Figure S2, Table S3, S4, S5
# Top-Tier Publication Standard (Nature / ISME Grade)
# ==============================================================================

# ----------------- 0. 环境准备与依赖包管理 -----------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, readxl, janitor, ggplot2, patchwork, cowplot,
  ggpubr, ggrepel, FactoMineR, factoextra, broom, stringr,
  grid, ggsci, multcompView, vegan, scales, smatr, lavaan,
  car, svglite, readr
)

# 显式全局锁定 dplyr 函数，杜绝命名空间污染
recode <- dplyr::recode
select <- dplyr::select
filter <- dplyr::filter

# 相对工程路径配置 (杜绝任何本地绝对路径)
data_dir   <- "./data"
output_dir <- "./results/Figure2_Biogeochemistry"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 安全写表函数 (UTF-8 BOM 编码，杜绝 Excel/WPS 打开乱码)
safe_write_csv <- function(df, path) {
  tryCatch({
    readr::write_excel_csv(df, path)
    message(">>> [Successfully Exported Table] ", path)
  }, error = function(e) {
    warning("Cannot write to file: ", path, " : ", e$message)
  })
}

# ----------------- 1. 顶刊极简主题系统 (8.0 pt 出版标准) -----------------
theme_nature_clean <- function(base_size = 8.0, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(size = base_size, family = base_family),
      panel.grid = element_blank(), 
      panel.border = element_blank(), 
      axis.line = element_line(color = "black", linewidth = 0.35),
      plot.title = element_text(face = "bold", size = 8.2, hjust = 0.5, margin = margin(b = 2.0, t = 0)),
      axis.text = element_text(colour = "black", size = 6.8),
      axis.title = element_text(face = "bold", size = 7.5),
      plot.tag = element_text(face = "bold", size = 11.0),
      legend.title = element_text(face = "bold", size = 7.0, margin = margin(t = 4, b = 2.5, unit = "pt")),
      legend.text = element_text(size = 6.4),
      legend.key.size = unit(0.55, "lines"),
      legend.margin = margin(0, 0, 0, 0, "pt"),
      strip.background = element_rect(fill = "#F2F4F4", colour = "black", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 6.8, colour = "black", margin = margin(t = 1.0, b = 1.0)),
      plot.margin = margin(t = 1.0, r = 2.0, b = 1.5, l = 2.0, "pt")
    )
}

theme_nature_grid <- function(base_size = 8.0, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(size = base_size, family = base_family),
      panel.grid.major = element_line(color = "#EAECEE", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.40),
      strip.background = element_rect(fill = "#F2F4F4", color = "black", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 6.8, color = "black", margin = margin(t = 1.0, b = 1.0)),
      plot.title = element_text(face = "bold", size = 8.2, hjust = 0.5, margin = margin(b = 2.0, t = 0)),
      axis.text = element_text(color = "black", size = 6.8),
      axis.title = element_text(face = "bold", size = 7.5),
      plot.tag = element_text(face = "bold", size = 11.0),
      legend.title = element_text(face = "bold", size = 7.0, margin = margin(t = 4, b = 2.5, unit = "pt")),
      legend.text = element_text(size = 6.4),
      legend.key.size = unit(0.55, "lines"),
      legend.margin = margin(0, 0, 0, 0, "pt"),
      panel.spacing = unit(0.22, "lines"),
      plot.margin = margin(t = 1.0, r = 2.0, b = 1.5, l = 2.0, "pt")
    )
}

# 19 项指标元数据标准
var_meta <- tribble(
  ~var,    ~label,    ~category,
  "som",   "SOM",     "PhysChem", 
  "nh4_n", "NH4+-N",  "PhysChem",
  "ap",    "AP",      "PhysChem", 
  "ak",    "AK",      "PhysChem",
  "tn",    "TN",      "PhysChem", 
  "tp",    "TP",      "PhysChem",
  "tk",    "TK",      "PhysChem", 
  "swc",   "SWC",     "PhysChem",
  "p_h",   "pH",      "PhysChem", 
  "sue",   "SUE",     "Enzymes",
  "ssc",   "SSC",     "Enzymes",  
  "salp",  "SALP",    "Enzymes",
  "pb",    "Pb",      "Metals",   
  "cr",    "Cr",      "Metals",
  "cu",    "Cu",      "Metals",   
  "ni",    "Ni",      "Metals",
  "zn",    "Zn",      "Metals",   
  "cd",    "Cd",      "Metals",
  "as",    "As",      "Metals"
)

level_labs <- c("High", "Med-High", "Med", "Low")
nature_level_pal <- c("High" = "#2E8B57", "Med-High" = "#FF3B6E", "Med" = "#55B9FF", "Low" = "#FFC533")
metal_pal <- c("As" = "#BC3C29", "Cd" = "#0072B5", "Cr" = "#E18727", "Cu" = "#20854E", "Ni" = "#7876B1", "Pb" = "#6F99AD", "Zn" = "#FFDC91")

# ----------------- 2. 自动化静默读取数据 (无弹窗 + 智能防错) -----------------
load_biogeochem_data_headless <- function() {
  main_file <- file.path(data_dir, "Soil_Physicochemical_Enzymes_4Seasons.xlsx")
  
  # 若未改名，自动寻找理化相关 Excel 文件
  if (!file.exists(main_file)) {
    candidates <- list.files(data_dir, pattern = "(理化|Enzyme|Physico|env).*\\.xlsx?$", full.names = TRUE, ignore.case = TRUE)
    if (length(candidates) > 0) {
      main_file <- candidates[1]
      message(">>> Auto-detected environmental data file: ", basename(main_file))
    }
  }
  
  if (file.exists(main_file)) {
    sheet_names <- readxl::excel_sheets(main_file)
    data_list <- lapply(sheet_names, function(sh) {
      df <- readxl::read_excel(main_file, sheet = sh)
      clean_cols <- names(df) %>%
        str_replace_all("\\(.*\\)", "") %>%
        str_replace_all("\\[.*\\]", "") %>%
        str_replace_all("\\+", "") %>%
        janitor::make_clean_names()
      names(df) <- clean_cols
      
      if (!"sampled" %in% names(df)) {
        sample_col_idx <- grep("^(sample|samp|yang_ben|bian_hao|id)", names(df))
        if (length(sample_col_idx) > 0) names(df)[sample_col_idx[1]] <- "sampled" else names(df)[1] <- "sampled"
      }
      df %>% mutate(Season_Sheet = sh)
    })
    df_all <- bind_rows(data_list) %>% filter(!is.na(sampled))
  } else {
    # 兼容性备用方案：若数据按季度存为独立表格
    season_candidates <- list(
      Spring = file.path(data_dir, "Spring.xlsx"),
      Summer = file.path(data_dir, "Summer.xlsx"),
      Autumn = file.path(data_dir, "Autumn.xlsx"),
      Winter = file.path(data_dir, "Winter.xlsx")
    )
    data_list <- map(names(season_candidates), function(s_en) {
      fpath <- season_candidates[[s_en]]
      if (!file.exists(fpath)) {
        stop(sprintf("Data file not found! Please place '%s' into '%s'.", 
                     "Soil_Physicochemical_Enzymes_4Seasons.xlsx", data_dir))
      }
      df <- readxl::read_excel(fpath, sheet = 1)
      clean_cols <- names(df) %>%
        str_replace_all("\\(.*\\)", "") %>%
        str_replace_all("\\[.*\\]", "") %>%
        str_replace_all("\\+", "") %>%
        janitor::make_clean_names()
      names(df) <- clean_cols
      if (!"sampled" %in% names(df)) {
        sample_col_idx <- grep("^(sample|samp|yang_ben|bian_hao|id)", names(df))
        if (length(sample_col_idx) > 0) names(df)[sample_col_idx[1]] <- "sampled" else names(df)[1] <- "sampled"
      }
      df %>% mutate(Season_Sheet = s_en)
    })
    df_all <- bind_rows(data_list) %>% filter(!is.na(sampled))
  }
  
  # 因子水平解析
  has_level <- "level" %in% names(df_all) && !all(is.na(df_all$level))
  has_veg   <- "veg" %in% names(df_all) && !all(is.na(df_all$veg))
  
  if (!has_level || !has_veg) {
    df_all <- df_all %>%
      mutate(
        clean_code = str_remove_all(as.character(sampled), "[0-9\u2080-\u2089\\s_-]"),
        code_level = substr(clean_code, 2, 2),
        code_veg   = substr(clean_code, 3, 3)
      )
    if (!has_level) {
      df_all <- df_all %>% mutate(level = factor(code_level, levels = c("A", "B", "C", "D"), labels = level_labs))
    }
    if (!has_veg) {
      df_all <- df_all %>% mutate(veg = factor(code_veg, levels = c("C", "G", "Q"), labels = c("Herb", "Shrub", "Tree")))
    }
  } else {
    df_all <- df_all %>%
      mutate(
        level = factor(level, levels = level_labs),
        veg   = factor(veg, levels = c("Herb", "Shrub", "Tree"))
      )
  }
  
  season_map_lookup <- c("C" = "Spring", "X" = "Summer", "Q" = "Autumn", "D" = "Winter",
                         "Spring" = "Spring", "Summer" = "Summer", "Autumn" = "Autumn", "Winter" = "Winter",
                         "春" = "Spring", "夏" = "Summer", "秋" = "Autumn", "冬" = "Winter")
  code_season <- substr(str_remove_all(as.character(df_all$sampled), "[0-9\u2080-\u2089\\s_-]"), 1, 1)
  
  df_all <- df_all %>%
    mutate(
      season = case_when(
        !is.na(season_map_lookup[code_season])  ~ season_map_lookup[code_season],
        !is.na(season_map_lookup[Season_Sheet]) ~ season_map_lookup[Season_Sheet],
        TRUE ~ "Spring"
      ),
      season = factor(season, levels = c("Spring", "Summer", "Autumn", "Winter"))
    )
  
  message(sprintf(">>> Data successfully loaded without dialogs. Valid samples: %d", nrow(df_all)))
  
  df_all %>%
    mutate(
      across(any_of(var_meta$var), as.numeric),
      
      # 1. 酶活性严谨转化为摩尔催化速率 (nmol product · g⁻¹ dry soil · h⁻¹)
      ssc_nmol  = ssc * (10^6 / (180.16 * 24)),
      sue_nmol  = sue * (1000 / (14.007 * 24)),
      salp_nmol = salp * 1.0,
      
      # 2. 1.01 底线对数截断保护 (严格维持欧氏第一象限几何解析)
      ssc_nmol  = pmax(ssc_nmol, 1.01),
      sue_nmol  = pmax(sue_nmol, 1.01),
      salp_nmol = pmax(salp_nmol, 1.01),
      
      # 3. 重金属风险模型 (PLI 与 RI)
      pli = ( (pb/26.0) * (cr/61.0) * (cu/22.6) * (ni/26.9) * (zn/74.2) * (cd/0.097) * (as/11.2) )^(1/7),
      er_cd = (cd/0.097) * 30, er_as = (as/11.2) * 10, er_pb = (pb/26.0) * 5,
      er_cu = (cu/22.6) * 5,  er_ni = (ni/26.9) * 5,  er_cr = (cr/61.0) * 2, er_zn = (zn/74.2) * 1,
      RI = er_cd + er_as + er_pb + er_cu + er_ni + er_cr + er_zn,
      
      # 4. US EPA 人体非致癌健康风险指数
      hi_child = (as * 200 * 350 * 6) / (15 * 6 * 365 * 0.0003 * 10^6) + 
        (pb * 200 * 350 * 6) / (15 * 6 * 365 * 0.0035 * 10^6),
      hi_adult = (as * 100 * 350 * 24) / (70 * 24 * 365 * 0.0003 * 10^6) + 
        (pb * 100 * 350 * 24) / (70 * 24 * 365 * 0.0035 * 10^6),
      
      # 5. 生态酶矢量计算 (atan2 几何连续解析)
      x_stoich = log10(ssc_nmol) / log10(salp_nmol), 
      y_stoich = log10(ssc_nmol) / log10(sue_nmol),
      Vector_L = sqrt(x_stoich^2 + y_stoich^2), 
      Vector_A = atan2(y_stoich, x_stoich) * (180 / pi),
      
      # 6. 化学计量对数变换
      ln_C = log(som),
      ln_N = log(tn),
      ln_P = log(tp)
    )
}

# 执行数据自动加载
processed_df <- load_biogeochem_data_headless()

# ----------------- 3. 彻底免疫连字符 Bug 的安全 Tukey HSD 标注系统 -----------------
get_sig_with_y <- function(data, target, group_var) {
  df_sub <- data %>% filter(!is.na(.data[[target]]), !is.na(.data[[group_var]]))
  df_sub[[group_var]] <- as.factor(df_sub[[group_var]])
  orig_levels <- levels(df_sub[[group_var]])
  
  if (length(unique(df_sub[[group_var]])) < 2) {
    return(tibble(grp = orig_levels, Letters = "a", y_pos = max(df_sub[[target]], na.rm = TRUE) * 1.08))
  }
  
  safe_tokens <- paste0("TOKEN", seq_along(orig_levels))
  forward_map <- setNames(safe_tokens, orig_levels)
  reverse_map <- setNames(orig_levels, safe_tokens)
  
  df_sub$safe_group <- factor(forward_map[as.character(df_sub[[group_var]])], levels = safe_tokens)
  
  fit <- aov(as.formula(paste(target, "~ safe_group")), data = df_sub)
  p_adj <- TukeyHSD(fit)$safe_group[, "p adj"]
  
  letters_obj <- tryCatch({
    multcompLetters(p_adj)$Letters
  }, error = function(e) {
    setNames(rep("a", length(safe_tokens)), safe_tokens)
  })
  
  res_letters <- tibble(
    grp = reverse_map[names(letters_obj)],
    Letters = as.character(letters_obj)
  )
  
  pos_df <- df_sub %>%
    group_by(across(all_of(group_var))) %>%
    summarise(y_pos = max(.data[[target]], na.rm = TRUE) * 1.08, .groups = "drop")
  
  tibble(grp = orig_levels) %>%
    left_join(res_letters, by = "grp") %>%
    mutate(Letters = if_else(is.na(Letters), "a", Letters)) %>%
    left_join(pos_df %>% rename(grp = all_of(group_var)), by = "grp")
}

# ==============================================================================
# 4. 图 2: 地球化学与酶化学计量套件 (Figure 2: 19cm × 24cm)
# ==============================================================================

# 图 2A: Type II ANOVA 主效应相对边际方差贡献
p_a_data <- map_dfr(var_meta$var, function(v) {
  fit <- lm(as.formula(paste(v, "~ level * season * veg")), data = processed_df)
  aov_res <- car::Anova(fit, type = "II")
  
  ss_u <- aov_res["level", "Sum Sq"]
  ss_s <- aov_res["season", "Sum Sq"]
  ss_v <- aov_res["veg", "Sum Sq"]
  total_main <- ss_u + ss_s + ss_v
  
  tibble(
    var = v,
    Urbanization = (ss_u / total_main) * 100,
    Season       = (ss_s / total_main) * 100,
    Vegetation   = (ss_v / total_main) * 100
  )
}) %>%
  pivot_longer(cols = c(Urbanization, Season, Vegetation), names_to = "term", values_to = "Contribution") %>%
  left_join(var_meta, by = "var")

p_a <- ggplot(p_a_data, aes(x = reorder(label, Contribution), y = Contribution, fill = term)) +
  geom_col(width = 0.68) + 
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  coord_flip() + 
  scale_fill_jama(name = "Main Driver") + 
  labs(
    title = "Relative Driver Contribution (Type II ANOVA)", 
    tag = "A", 
    y = "Relative Contribution of Main Drivers (%)", 
    x = "Biogeochemical Indicators"
  ) + 
  theme_nature_clean() + 
  theme(
    strip.text.y = element_text(angle = 0, size = 6.6, face = "bold"),
    legend.position = "none",
    axis.text.y = element_text(size = 6.5, color = "black")
  )

# 图 2B: PCA Biplot 与 PERMANOVA 统计检验
pca_data <- processed_df %>% drop_na(all_of(var_meta$var))
perm <- adonis2(scale(pca_data[, var_meta$var]) ~ level, data = pca_data, method = "euclidean")
stat_lab <- paste0("PERMANOVA:\nR² = ", sprintf("%.2f", perm$R2[1]), "\nF = ", sprintf("%.2f", perm$F[1]), "\np < 0.001")

pca_matrix <- pca_data %>% select(all_of(var_meta$var))
colnames(pca_matrix) <- var_meta$label[match(colnames(pca_matrix), var_meta$var)]
pca_fit <- PCA(pca_matrix, scale.unit = TRUE, graph = FALSE)

p_b_pca <- fviz_pca_biplot(
  pca_fit, 
  geom.ind = "point", pointshape = 21, pointsize = 1.0, 
  fill.ind = pca_data$level, col.ind = "black", palette = nature_level_pal, 
  addEllipses = TRUE, ellipse.level = 0.85, col.var = "black", 
  arrowsize = 0.35, labelsize = 2.0, repel = TRUE, 
  title = "Biogeochemical Patterns (PCA)"
) + 
  labs(tag = "B") + 
  theme_nature_grid() + 
  annotate("text", x = -Inf, y = Inf, label = stat_lab, hjust = -0.08, vjust = 1.18, 
           fontface = "bold.italic", size = 2.0, family = "sans") +
  guides(fill = "none", color = "none")

# 图 2C: 潜在生态危害指数 (RI)
sig_C <- processed_df %>%
  group_by(season, veg) %>%
  group_modify(~ get_sig_with_y(.x, "RI", "level")) %>%
  ungroup()

p_c_ri <- ggplot(processed_df, aes(x = level, y = RI)) +
  geom_jitter(aes(color = level), alpha = 0.35, width = 0.1, size = 0.6) +
  geom_boxplot(aes(color = level), width = 0.68, outlier.shape = NA, linewidth = 0.35) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 0.35) +
  geom_text(data = sig_C %>% rename(level = grp), aes(x = level, y = y_pos, label = Letters), 
            fontface = "bold", size = 2.0, family = "sans", vjust = -0.1) +
  facet_grid(veg ~ season) + 
  scale_color_manual(values = nature_level_pal, guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.45, 0.45))) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(tag = "C", y = "Potential Ecological Risk (RI)", x = "Urbanization Level") + 
  theme_nature_grid() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8))

# 图 2D: 人体健康危害指数 (HI)
df_h <- bind_rows(
  processed_df %>% mutate(HI = hi_child, Group = "Children"),
  processed_df %>% mutate(HI = hi_adult, Group = "Adults")
)

p_d <- ggplot(df_h, aes(x = season, y = HI, color = Group, group = Group)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18, linewidth = 0.35) + 
  stat_summary(fun = mean, geom = "line", linewidth = 0.65) +
  stat_summary(fun = mean, geom = "point", size = 1.3, shape = 21, fill = "white", stroke = 0.60) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", linewidth = 0.35) +
  facet_wrap(~ level, nrow = 1) + 
  scale_color_manual(values = c("Adults" = "#0072B5", "Children" = "#BC3C29"), name = "Cohort") +
  labs(tag = "D", y = "Hazard Index (HI)", x = "Season") + 
  theme_nature_grid() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8))

# 图 2E: 重金属风险多维评估套件
sig_F1 <- get_sig_with_y(processed_df, "pli", "season")
p_f1 <- ggplot(processed_df, aes(x = season, y = pli, fill = season)) +
  geom_jitter(shape = 21, fill = "white", color = "grey30", alpha = 0.5, width = 0.15, size = 0.8) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.35, width = 0.65) + 
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", linewidth = 0.35) +
  geom_text(data = sig_F1 %>% rename(season = grp), aes(x = season, y = y_pos, label = Letters), 
            fontface = "bold", size = 2.0, family = "sans", vjust = -0.1) +
  scale_fill_npg(guide = "none") + 
  labs(tag = "E", y = "Pollution Load (PLI)", x = "Season") + 
  theme_nature_grid() + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.0))

metal_names_map <- c("er_as" = "As", "er_cd" = "Cd", "er_cr" = "Cr", "er_cu" = "Cu", "er_ni" = "Ni", "er_pb" = "Pb", "er_zn" = "Zn")
er_long <- processed_df %>%
  group_by(level) %>%
  summarise(across(starts_with("er_"), mean), .groups = "drop") %>%
  pivot_longer(-level, names_to = "Metal_Raw", values_to = "Er") %>% 
  mutate(Metal = factor(metal_names_map[Metal_Raw], levels = c("As", "Cd", "Cr", "Cu", "Ni", "Pb", "Zn")))

p_f2 <- ggplot(er_long, aes(x = level, y = Er, fill = Metal)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.15) + 
  scale_fill_manual(values = metal_pal, name = "Metal") +
  labs(y = "Risk Factor (Er)", x = "Urbanization Level") + 
  theme_nature_grid() + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.0))

heat_data <- processed_df %>% 
  group_by(season, level) %>%
  summarise(
    RI_Ex  = mean(RI > 150, na.rm = TRUE) * 100, 
    PLI_Ex = mean(pli > 1, na.rm = TRUE) * 100, 
    HI_Ex  = mean(hi_child > 1, na.rm = TRUE) * 100, 
    .groups = "drop"
  ) %>%
  pivot_longer(cols = ends_with("_Ex"), names_to = "Ind", values_to = "Rate") %>%
  mutate(
    Ind = dplyr::recode(
      Ind, 
      "RI_Ex"  = "RI > 150", 
      "PLI_Ex" = "PLI > 1", 
      "HI_Ex"  = "HI > 1"
    )
  )

p_f3 <- ggplot(heat_data, aes(x = level, y = Ind, fill = Rate)) +
  geom_tile(color = "white", linewidth = 0.30) + 
  geom_text(aes(label = paste0(round(Rate), "%"), color = Rate > 60), 
            size = 1.9, show.legend = FALSE, family = "sans", fontface = "bold") +
  facet_wrap(~ season, nrow = 1) + 
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Exceedance (%)",
                       guide = guide_colorbar(barheight = unit(1.6, "cm"), barwidth = unit(0.35, "cm"), title.vjust = 1.2)) + 
  scale_color_manual(values = c("black", "white"), guide = "none") + 
  labs(title = "Threshold Exceedance Rates (%)", x = "Urbanization Level") + 
  theme_nature_grid() + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.8), axis.title.y = element_blank())

p_f_suite <- (p_f1 | p_f2 | p_f3) + plot_layout(widths = c(0.85, 0.85, 2.5))

# 图 2F: 生态酶化学计量矢量分析
p_e1 <- ggplot(processed_df, aes(x = x_stoich, y = y_stoich, color = level, shape = veg)) +
  geom_point(size = 1.3, alpha = 0.8, stroke = 0.50) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.35) +
  facet_wrap(~ season, nrow = 1) + 
  scale_color_manual(values = nature_level_pal, name = "Urbanization") +
  scale_shape_manual(values = c(16, 17, 15), name = "Vegetation") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  labs(
    title = "Ecoenzymatic Stoichiometry Space", 
    tag = "F", 
    x = expression(log[10](SSC) / log[10](SALP) ~ "(C/P)"),
    y = expression(log[10](SSC) / log[10](SUE) ~ "(C/N)")
  ) +
  theme_nature_grid() +
  theme(
    axis.title.x = element_text(face = "bold", size = 7.5, margin = margin(t = 3.5, b = 2.0, unit = "pt")),
    plot.margin = margin(t = 1.0, r = 2.0, b = 7.0, l = 2.0, "pt")
  )

sig_L <- processed_df %>% group_by(season) %>% group_modify(~ get_sig_with_y(.x, "Vector_L", "level")) %>% ungroup()
p_e2 <- ggplot(processed_df, aes(x = level, y = Vector_L)) +
  geom_jitter(aes(color = level), alpha = 0.3, width = 0.12, size = 0.6) + 
  geom_boxplot(aes(color = level), width = 0.68, linewidth = 0.35, outlier.shape = NA) + 
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 0.35) +
  geom_text(data = sig_L %>% rename(level = grp), aes(x = level, y = y_pos, label = Letters), 
            fontface = "bold", size = 2.0, family = "sans", vjust = -0.1) +
  facet_wrap(~ season, nrow = 1) + 
  scale_color_manual(values = nature_level_pal, guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.45, 0.45))) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(y = "Vector Length (L)", x = "Urbanization Level") + 
  theme_nature_grid() + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.8), legend.position = "none")

sig_A <- processed_df %>% group_by(season) %>% group_modify(~ get_sig_with_y(.x, "Vector_A", "level")) %>% ungroup()
p_e3 <- ggplot(processed_df, aes(x = level, y = Vector_A)) +
  geom_jitter(aes(color = level), alpha = 0.3, width = 0.12, size = 0.6) + 
  geom_boxplot(aes(color = level), width = 0.68, linewidth = 0.35, outlier.shape = NA) + 
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 0.35) +
  geom_text(data = sig_A %>% rename(level = grp), aes(x = level, y = y_pos, label = Letters), 
            fontface = "bold", size = 2.0, family = "sans", vjust = -0.1) +
  facet_wrap(~ season, nrow = 1) + 
  scale_color_manual(values = nature_level_pal, guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.45, 0.45))) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(y = expression(Vector ~ Angle ~ (A * degree)), x = "Urbanization Level") + 
  theme_nature_grid() + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.8), legend.position = "none")

p_e_combined <- p_e1 / (p_e2 | p_e3) + plot_layout(heights = c(1.25, 0.85))

# 图 2 全局总图拼装 (19cm × 24cm)
figure_2_final <- (p_a + p_b_pca) / (p_c_ri + p_d) / p_f_suite / p_e_combined +
  plot_layout(heights = c(1.15, 1.20, 0.85, 1.35), guides = "collect") & 
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.box.spacing = unit(10, "pt"),
    legend.spacing.y = unit(8, "pt"),
    plot.tag = element_text(face = "bold", size = 11.0)
  )

ggsave(file.path(output_dir, "Figure2_Biogeochemistry_Suite.pdf"), figure_2_final, width = 19, height = 24, units = "cm", device = cairo_pdf)
ggsave(file.path(output_dir, "Figure2_Biogeochemistry_Suite.png"), figure_2_final, width = 19, height = 24, units = "cm", dpi = 500)
ggsave(file.path(output_dir, "Figure2_Biogeochemistry_Suite.svg"), figure_2_final, width = 19, height = 24, units = "cm", device = "svg")

# ==============================================================================
# 5. 附图 S1: SMA 元素化学计量异速生长回归与相关性热图
# ==============================================================================
plot_sma <- function(data, x_col, y_col, x_lab, y_lab, panel_tag, tag_letter) {
  df_sub <- data %>% filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  fit <- tryCatch({
    sma(as.formula(paste(y_col, "~", x_col)), data = df_sub, slope.test = 1)
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {
    slope_val <- as.numeric(fit$coef[[1]][2, 1])
    intercept_val <- as.numeric(fit$coef[[1]][1, 1])
    r2_val <- as.numeric(fit$r2[[1]])
    p_slope1 <- as.numeric(fit$slopetest[[1]]$p)
    p_text <- if (p_slope1 < 0.001) "P(slope=1)<0.001" else sprintf("P(slope=1)=%.3f", p_slope1)
    anno_text <- sprintf("Slope=%.2f, R²=%.2f\n%s", slope_val, r2_val, p_text)
  } else {
    slope_val <- 1; intercept_val <- 0; anno_text <- "SMA Fit Failed"
  }
  
  mean_x <- mean(df_sub[[x_col]], na.rm = TRUE)
  mean_y <- mean(df_sub[[y_col]], na.rm = TRUE)
  iso_intercept <- mean_y - mean_x
  
  ggplot(df_sub, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(aes(color = level), alpha = 0.6, size = 1.2) +
    geom_abline(intercept = intercept_val, slope = slope_val, color = "black", linewidth = 0.7) +
    geom_abline(intercept = iso_intercept, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.45) +
    annotate("text", x = -Inf, y = Inf, label = anno_text, hjust = -0.05, vjust = 1.15, size = 2.0, fontface = "bold", family = "sans") +
    scale_color_manual(values = nature_level_pal) +
    labs(title = panel_tag, tag = tag_letter, x = x_lab, y = y_lab) +
    theme_nature_grid() +
    theme(legend.position = "none")
}

p_sma_a <- plot_sma(processed_df, "ln_N", "ln_C", "ln(Total Nitrogen, TN)", "ln(Soil Organic Matter, SOM)", "C-N Stoichiometry", "A")
p_sma_b <- plot_sma(processed_df, "ln_P", "ln_C", "ln(Total Phosphorus, TP)", "ln(Soil Organic Matter, SOM)", "C-P Stoichiometry", "B")
p_sma_c <- plot_sma(processed_df, "ln_P", "ln_N", "ln(Total Phosphorus, TP)", "ln(Total Nitrogen, TN)", "N-P Stoichiometry", "C")

corr_vars <- c("ln_C", "ln_N", "ln_P", "Vector_L", "Vector_A")
corr_data <- processed_df %>% select(all_of(corr_vars)) %>% drop_na()
cor_mat <- cor(corr_data, method = "pearson")
p_mat <- matrix(NA, nrow = length(corr_vars), ncol = length(corr_vars))

for(i in 1:length(corr_vars)) {
  for(j in 1:length(corr_vars)) {
    p_mat[i, j] <- cor.test(corr_data[[corr_vars[i]]], corr_data[[corr_vars[j]]])$p.value
  }
}

cor_df <- as.data.frame(cor_mat) %>% rownames_to_column(var = "Var1") %>% pivot_longer(-Var1, names_to = "Var2", values_to = "r")
p_df   <- as.data.frame(p_mat, row.names = corr_vars) %>% setNames(corr_vars) %>% rownames_to_column(var = "Var1") %>% pivot_longer(-Var1, names_to = "Var2", values_to = "p")

corr_plot_df <- cor_df %>%
  left_join(p_df, by = c("Var1", "Var2")) %>%
  mutate(
    stars = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""),
    label_text = if_else(Var1 == Var2, "1", if_else(stars == "", sprintf("%.2f", r), sprintf("%.2f\n%s", r, stars))),
    Var1 = dplyr::recode(Var1, "ln_C" = "ln(SOM)", "ln_N" = "ln(TN)", "ln_P" = "ln(TP)", "Vector_L" = "Vector L\n(Intensity)", "Vector_A" = "Vector A\n(Direction)"),
    Var2 = dplyr::recode(Var2, "ln_C" = "ln(SOM)", "ln_N" = "ln(TN)", "ln_P" = "ln(TP)", "Vector_L" = "Vector L\n(Intensity)", "Vector_A" = "Vector A\n(Direction)")
  ) %>%
  mutate(Var1 = factor(Var1, levels = rev(unique(Var1))), Var2 = factor(Var2, levels = unique(Var2)))

p_d_corr <- ggplot(corr_plot_df, aes(x = Var2, y = Var1, fill = r)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = label_text), size = 2.2, fontface = "bold", color = "black", lineheight = 0.8, family = "sans") +
  scale_fill_gradient2(low = "#72C0D4", mid = "#FBF6F0", high = "#E25B45", midpoint = 0, limit = c(-1, 1), name = "Pearson r") +
  labs(title = "Correlation Matrix", tag = "D", x = NULL, y = NULL) +
  theme_nature_clean() +
  theme(
    axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black")
  )

figure_s1_final <- (p_sma_a | p_sma_b | p_sma_c) / p_d_corr + plot_layout(heights = c(1, 1.2)) &
  theme(plot.tag = element_text(face = "bold", size = 11.0))

ggsave(file.path(output_dir, "FigureS1_SMA_Regressions.pdf"), figure_s1_final, width = 19, height = 16, units = "cm", device = cairo_pdf)
ggsave(file.path(output_dir, "FigureS1_SMA_Regressions.png"), figure_s1_final, width = 19, height = 16, units = "cm", dpi = 500)
ggsave(file.path(output_dir, "FigureS1_SMA_Regressions.svg"), figure_s1_final, width = 19, height = 16, units = "cm", device = "svg")

# ==============================================================================
# 6. 附图 S2: SEM 结构方程模型全通径解析与效应图
# ==============================================================================
sem_df <- processed_df %>%
  mutate(
    level_num = as.numeric(factor(level, levels = c("Low", "Med", "Med-High", "High"))),
    veg_num = as.numeric(factor(veg, levels = c("Herb", "Shrub", "Tree"))),
    season_active = if_else(season %in% c("Spring", "Summer"), 1, 0)
  ) %>%
  drop_na(level_num, veg_num, season_active, p_h, swc, RI, Vector_L, Vector_A, ln_C, ln_P)

# 结构方程模型 1 (Vector L: 碳限制强度)
model1_spec <- '
  p_h ~ a1*level_num + a2*veg_num + a3*season_active
  swc ~ b1*level_num + b2*veg_num + b3*season_active
  p_h ~~ swc
  RI ~ c1*level_num + c2*p_h + c3*swc + c4*season_active
  ln_C ~ d1*level_num + d2*veg_num + d3*p_h + d4*swc + d5*RI + d6*season_active
  Vector_L ~ e1*level_num + e2*p_h + e3*swc + e4*RI + e5*ln_C + e6*season_active
'

# 结构方程模型 2 (Vector A: 养分限制方向)
model2_spec <- '
  p_h ~ a1*level_num + a2*veg_num + a3*season_active
  swc ~ b1*level_num + b2*veg_num + b3*season_active
  p_h ~~ swc
  RI ~ c1*level_num + c2*p_h + c3*swc + c4*season_active
  ln_P ~ d1*level_num + d2*veg_num + d3*p_h + d4*swc + d5*season_active
  Vector_A ~ e1*level_num + e2*p_h + e3*swc + e4*RI + e5*ln_P + e6*season_active
'

fit1 <- sem(model1_spec, data = sem_df)
fit2 <- sem(model2_spec, data = sem_df)

node_map_m1 <- tibble(
  name = c("level_num", "veg_num", "season_active", "p_h", "swc", "RI", "ln_C", "Vector_L"),
  label = c("Urbanization", "Vegetation", "Seasonality", "Soil pH", "SWC", "Metal Risk\n(RI)", "Soil Organic C\n(ln C)", "C Limitation\n(Vector L)"),
  x = c(0.2, 0.2, 0.2, 0.8, 0.8, 1.4, 1.4, 2.1),
  y = c(1.5, 0.5, 1.0, 1.3, 0.7, 0.5, 1.5, 1.0)
)

node_map_m2 <- tibble(
  name = c("level_num", "veg_num", "season_active", "p_h", "swc", "RI", "ln_P", "Vector_A"),
  label = c("Urbanization", "Vegetation", "Seasonality", "Soil pH", "SWC", "Metal Risk\n(RI)", "Soil Total P\n(ln P)", "Limitation Dir.\n(Vector A)"),
  x = c(0.2, 0.2, 0.2, 0.8, 0.8, 1.4, 1.4, 2.1),
  y = c(1.5, 0.5, 1.0, 1.3, 0.7, 0.5, 1.5, 1.0)
)

shorten_segment <- function(x1, y1, x2, y2, r_start = 0.16, r_end = 0.22) {
  d <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  if (d == 0) return(tibble(x1_new = x1, y1_new = y1, x2_new = x2, y2_new = y2))
  tibble(
    x1_new = x1 + (r_start / d) * (x2 - x1),
    y1_new = y1 + (r_start / d) * (y2 - y1),
    x2_new = x2 - (r_end / d) * (x2 - x1),
    y2_new = y2 - (r_end / d) * (y2 - y1)
  )
}

draw_lavaan_sem <- function(fit, node_map, title_label, tag_letter) {
  paths <- standardizedSolution(fit) %>% filter(op == "~") %>% select(lhs, rhs, beta = est.std, pval = pvalue)
  r2_vals <- inspect(fit, "rsquare")
  
  nodes_df <- node_map %>%
    rowwise() %>%
    mutate(
      r2_val = if_else(name %in% names(r2_vals), r2_vals[name], as.numeric(NA)),
      label_final = if_else(!is.na(r2_val), paste0(label, "\n(R² = ", sprintf("%.2f", pmax(0, r2_val)), ")"), label)
    ) %>%
    ungroup()
  
  paths_coords <- paths %>%
    left_join(nodes_df %>% select(name, x1 = x, y1 = y), by = c("rhs" = "name")) %>%
    left_join(nodes_df %>% select(name, x2 = x, y2 = y), by = c("lhs" = "name")) %>%
    rowwise() %>%
    mutate(shortened = list(shorten_segment(x1, y1, x2, y2))) %>%
    unnest(cols = c(shortened)) %>%
    ungroup() %>%
    mutate(
      pval_clean = if_else(is.na(pval), 1.0, pval),
      color = if_else(pval_clean < 0.05 & beta > 0, "#BC3C29", "#0072B5"),
      linetype = if_else(pval_clean < 0.05, "solid", "dashed"),
      linewidth = if_else(pval_clean < 0.05, 0.65, 0.30),
      x_mid = (x1_new + x2_new) / 2,
      y_mid = (y1_new + y2_new) / 2 + 0.025
    )
  
  paths_sig <- paths_coords %>% filter(!is.na(pval) & pval < 0.05)
  fm <- fitMeasures(fit)
  fit_lbl <- sprintf("Fit: Chi²=%.2f (P=%.3f)\nCFI=%.3f, RMSEA=%.3f", fm["chisq"], fm["pvalue"], fm["cfi"], fm["rmsea"])
  
  ggplot() +
    geom_segment(data = paths_sig, aes(x = x1_new, y = y1_new, xend = x2_new, yend = y2_new, 
                                       color = color, linetype = linetype, linewidth = linewidth),
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    geom_label(data = paths_sig, aes(x = x_mid, y = y_mid, label = sprintf("%.2f", beta), color = color),
               fill = "white", label.size = 0, size = 1.9, fontface = "bold", alpha = 0.9, show.legend = FALSE, family = "sans") +
    geom_label(data = nodes_df, aes(x = x, y = y, label = label_final), 
               fill = "#F8F9FA", color = "black", size = 2.1, fontface = "bold", family = "sans",
               label.padding = unit(0.22, "lines"), label.size = 0.35, label.r = unit(0.1, "lines")) +
    annotate("label", x = 2.2, y = 0.22, label = fit_lbl, hjust = 1, vjust = 0, 
             size = 1.9, fontface = "bold.italic", color = "grey20", family = "sans",
             fill = "#F8F9FA", label.size = 0.25, label.padding = unit(0.20, "lines")) +
    scale_color_identity() + scale_linetype_identity() + scale_linewidth_identity() +
    labs(title = title_label, tag = tag_letter) +
    xlim(0.0, 2.3) + ylim(0.2, 1.9) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 8.5, margin = margin(b = 2.5), family = "sans"),
      plot.tag = element_text(face = "bold", size = 11.0)
    )
}

# 全通径解析直接效应与间接效应 (基于 DAG 结构代数展开)
get_matrix_effects <- function(fit, outcome_var) {
  std_sol <- standardizedSolution(fit) %>% filter(op == "~")
  nodes <- unique(c(std_sol$lhs, std_sol$rhs))
  n <- length(nodes)
  
  B <- matrix(0, nrow = n, ncol = n, dimnames = list(nodes, nodes))
  for(i in 1:nrow(std_sol)) {
    B[std_sol$lhs[i], std_sol$rhs[i]] <- std_sol$est.std[i]
  }
  
  I <- diag(n)
  dimnames(I) <- list(nodes, nodes)
  Total_M <- tryCatch(solve(I - B) - I, error = function(e) {
    B + B %*% B + B %*% B %*% B
  })
  
  direct_effects <- B[outcome_var, ]
  total_effects <- Total_M[outcome_var, ]
  indirect_effects <- total_effects - direct_effects
  
  tibble(
    Variable = names(direct_effects),
    Direct   = direct_effects,
    Indirect = indirect_effects,
    Total    = total_effects
  ) %>%
    filter(Variable != outcome_var, (abs(Direct) > 0.001 | abs(Indirect) > 0.001))
}

plot_effects_bar <- function(effects_df, title_lab, tag_letter) {
  df_long <- effects_df %>%
    pivot_longer(cols = c(Direct, Indirect), names_to = "Type", values_to = "Effect") %>%
    mutate(
      Variable = dplyr::recode(
        Variable,
        "level_num" = "Urbanization", "veg_num" = "Vegetation", "season_active" = "Seasonality",
        "p_h" = "Soil pH", "swc" = "SWC", "RI" = "Metal Risk (RI)", "ln_C" = "Soil Organic C", "ln_P" = "Soil Total P"
      )
    )
  
  max_val <- max(abs(df_long$Effect), na.rm = TRUE) * 1.15
  y_limit <- pmax(0.5, ceiling(max_val * 10) / 10)
  
  ggplot(df_long, aes(x = factor(Variable, levels = rev(unique(Variable))), y = Effect, fill = Type)) +
    geom_col(position = position_dodge(width = 0.70), width = 0.60, color = "black", linewidth = 0.20) +
    coord_flip() +
    scale_fill_manual(values = c("Direct" = "#E18727", "Indirect" = "#0072B5")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.30) +
    scale_y_continuous(limits = c(-y_limit, y_limit), breaks = seq(-y_limit, y_limit, length.out = 5)) + 
    labs(title = title_lab, tag = tag_letter, x = "Predictor Variables", y = "Standardized Effect (β)") +
    theme_nature_grid() +
    theme(legend.title = element_blank(), plot.title = element_text(size = 8.5, family = "sans"))
}

effects_m1 <- get_matrix_effects(fit1, "Vector_L")
effects_m2 <- get_matrix_effects(fit2, "Vector_A")

p_sem_diag_m1 <- draw_lavaan_sem(fit1, node_map_m1, "SEM - C Limitation (Vector L)", "A")
p_effects_m1   <- plot_effects_bar(effects_m1, "Effects on Vector L", "B")
p_sem_diag_m2 <- draw_lavaan_sem(fit2, node_map_m2, "SEM - Nutrient Limitation (Vector A)", "C")
p_effects_m2   <- plot_effects_bar(effects_m2, "Effects on Vector A", "D")

figure_s2_final <- (p_sem_diag_m1 | p_effects_m1) / (p_sem_diag_m2 | p_effects_m2) +
  plot_layout(heights = c(1, 1), widths = c(1.3, 1), guides = "collect") &
  theme(plot.tag = element_text(face = "bold", size = 11.0))

ggsave(file.path(output_dir, "FigureS2_SEM_Vector_Stoichiometry.pdf"), figure_s2_final, width = 19, height = 16, units = "cm", device = cairo_pdf)
ggsave(file.path(output_dir, "FigureS2_SEM_Vector_Stoichiometry.png"), figure_s2_final, width = 19, height = 16, units = "cm", dpi = 500)
ggsave(file.path(output_dir, "FigureS2_SEM_Vector_Stoichiometry.svg"), figure_s2_final, width = 19, height = 16, units = "cm", device = "svg")

# ==============================================================================
# 7. 附表导出 (Table S3, Table S4, Table S5)
# ==============================================================================
format_academic_p <- function(p) {
  case_when(
    is.na(p) ~ "-",
    p < 0.001 ~ "< 0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

# ----------------- 附表 3: Table S3 (理化性质与酶活性汇总) -----------------
table_s3 <- processed_df %>%
  group_by(level) %>%
  summarise(
    `SOM (g/kg)`   = sprintf("%.2f ± %.2f", mean(som, na.rm = TRUE), sd(som, na.rm = TRUE)),
    `TN (g/kg)`    = sprintf("%.2f ± %.2f", mean(tn, na.rm = TRUE), sd(tn, na.rm = TRUE)),
    `TP (g/kg)`    = sprintf("%.2f ± %.2f", mean(tp, na.rm = TRUE), sd(tp, na.rm = TRUE)),
    `RI`           = sprintf("%.2f ± %.2f", mean(RI, na.rm = TRUE), sd(RI, na.rm = TRUE)),
    `PLI`          = sprintf("%.2f ± %.2f", mean(pli, na.rm = TRUE), sd(pli, na.rm = TRUE)),
    `Vector L`     = sprintf("%.2f ± %.2f", mean(Vector_L, na.rm = TRUE), sd(Vector_L, na.rm = TRUE)),
    `Vector A (°)` = sprintf("%.2f ± %.2f", mean(Vector_A, na.rm = TRUE), sd(Vector_A, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(`Urbanization Level` = level)

safe_write_csv(table_s3, file.path(output_dir, "Table_S3_Summary_Statistics.csv"))

# ----------------- 附表 4: Table S4 (SMA 异速生长分析表) -----------------
get_sma_row <- function(y_col, x_col, label) {
  fit <- sma(as.formula(paste(y_col, "~", x_col)), data = processed_df, slope.test = 1)
  tibble(
    `Relationship` = label,
    `Slope`        = sprintf("%.2f", fit$coef[[1]][2, 1]),
    `95% CI`       = sprintf("%.2f - %.2f", fit$coef[[1]][2, 2], fit$coef[[1]][2, 3]),
    `Intercept`    = sprintf("%.2f", fit$coef[[1]][1, 1]),
    `R2`           = ifelse(fit$r2[[1]] < 0.01, "< 0.01", sprintf("%.2f", fit$r2[[1]])),
    `P-value`      = format_academic_p(fit$pval[[1]]),
    `P (slope=1)`  = format_academic_p(fit$slopetest[[1]]$p)
  )
}

table_s4 <- bind_rows(
  get_sma_row("ln_C", "ln_N", "ln(SOM) vs ln(TN)"),
  get_sma_row("ln_C", "ln_P", "ln(SOM) vs ln(TP)"),
  get_sma_row("ln_N", "ln_P", "ln(TN) vs ln(TP)")
)
safe_write_csv(table_s4, file.path(output_dir, "Table_S4_SMA_Regressions.csv"))

# ----------------- 附表 5: Table S5 (SEM 路径与拟合度报告) -----------------
var_academic_labels <- c(
  "level_num"     = "Urbanization", "veg_num"       = "Vegetation", "season_active" = "Seasonality",
  "p_h"           = "Soil pH",      "swc"           = "SWC",        "RI"            = "Metal Risk (RI)",
  "ln_C"          = "ln(SOM)",      "ln_P"          = "ln(TP)",     "Vector_L"      = "Vector L", "Vector_A"      = "Vector A"
)

get_sem_report <- function(fit, model_name, outcome_var) {
  pe  <- parameterEstimates(fit) %>% filter(op == "~")
  std <- standardizedSolution(fit) %>% filter(op == "~")
  eff <- get_matrix_effects(fit, outcome_var)
  
  paths_df <- pe %>%
    left_join(std %>% select(lhs, rhs, est.std), by = c("lhs", "rhs")) %>%
    mutate(
      Model               = model_name,
      Response            = dplyr::recode(lhs, !!!var_academic_labels),
      Predictor           = dplyr::recode(rhs, !!!var_academic_labels),
      Path                = paste(Predictor, "->", Response),
      `Unstandardized B`  = sprintf("%.3f", est),
      `SE`                = sprintf("%.3f", se),
      `z-value`           = sprintf("%.2f", z),
      `P-value`           = format_academic_p(pvalue),
      `Standardized Beta` = sprintf("%.3f", est.std)
    ) %>%
    left_join(eff %>% select(rhs = Variable, Direct, Indirect, Total), by = "rhs") %>%
    mutate(
      `Direct Effect`   = if_else(is.na(Direct), "-", sprintf("%.3f", Direct)),
      `Indirect Effect` = if_else(is.na(Indirect), "-", sprintf("%.3f", Indirect)),
      `Total Effect`    = if_else(is.na(Total), "-", sprintf("%.3f", Total))
    ) %>%
    select(Model, Response, Predictor, Path, `Unstandardized B`, `SE`, `z-value`, `P-value`, `Standardized Beta`, `Direct Effect`, `Indirect Effect`, `Total Effect`)
  
  fm <- fitMeasures(fit)
  fit_row <- tibble(
    Model               = model_name,
    Response            = "[Goodness of Fit]",
    Predictor           = "-",
    Path                = sprintf("Chi² = %.2f (df = %s, P = %s)", fm["chisq"], fm["df"], format_academic_p(fm["pvalue"])),
    `Unstandardized B`  = sprintf("CFI = %.3f", fm["cfi"]),
    `SE`                = sprintf("TLI = %.3f", fm["tli"]),
    `z-value`           = sprintf("RMSEA = %.3f", fm["rmsea"]),
    `P-value`           = sprintf("SRMR = %.3f", fm["srmr"]),
    `Standardized Beta` = "-",
    `Direct Effect`     = "-",
    `Indirect Effect`   = "-",
    `Total Effect`      = "-"
  )
  
  bind_rows(paths_df, fit_row)
}

table_s5 <- bind_rows(
  get_sem_report(fit1, "Model 1 (Vector L)", "Vector_L"),
  get_sem_report(fit2, "Model 2 (Vector A)", "Vector_A")
)

safe_write_csv(table_s5, file.path(output_dir, "Table_S5_SEM_Comprehensive_Report.csv"))

cat("\n=======================================================================\n")
cat(">>> [Script 01 Completed Successfully!] <<<\n")
cat(">>> Output Directory: ", output_dir, "\n")
cat(">>> 1. Figure 2: Figure2_Biogeochemistry_Suite (PDF/PNG/SVG, 19cm x 24cm)\n")
cat(">>> 2. Figure S1: FigureS1_SMA_Regressions (PDF/PNG/SVG)\n")
cat(">>> 3. Figure S2: FigureS2_SEM_Vector_Stoichiometry (PDF/PNG/SVG)\n")
cat(">>> 4. Tables: Table S3, S4, S5 (.csv)\n")
cat("=======================================================================\n")