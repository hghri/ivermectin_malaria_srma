###############################################################################
# Ivermectin–malaria SRMA — TABLE-AS-FIGURE visualisations
#
#  A. Africa study map            (from Table 1)  -> tabfig_africa_map
#  B. Heterogeneity view (I2/tau2)(from Table 2)  -> tabfig_heterogeneity
#  C. Study-characteristics panel (from Table 1)  -> tabfig_study_panel
#
# Clean journal style: short title only, dataset used as collected.
# Export: vector PDF (Cairo if available, else base pdf) + 500 dpi LZW TIFF + PNG.
#
# install.packages(c("readxl","dplyr","stringr","ggplot2","forcats","scales",
#                    "tidyr","maps","ggrepel"))
###############################################################################

library(readxl); library(dplyr); library(stringr); library(ggplot2)
library(forcats); library(scales); library(tidyr)

save_lancet <- function(plot, name, width_mm, height_mm) {
  ## Base "pdf" device: portable everywhere (Cairo is broken on some Windows R
  ## builds yet still reports as available). Swap in cairo_pdf if you have
  ## working Cairo and want embedded fonts.
  ggsave(paste0(name, ".pdf"),  plot, device = "pdf",
         width = width_mm, height = height_mm, units = "mm")
  ggsave(paste0(name, ".tiff"), plot, device = "tiff", width = width_mm,
         height = height_mm, units = "mm", dpi = 500, compression = "lzw", bg = "white")
  ggsave(paste0(name, ".png"),  plot, device = "png", width = width_mm,
         height = height_mm, units = "mm", dpi = 300, bg = "white"); invisible(NULL)
}
theme_paper <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text = element_text(colour = "black"))

## intervention palette (shared across figures)
int_levels <- c("IVM alone", "IVM + DP", "IVM + SMC", "IVM + albendazole")
int_pal <- c("IVM alone" = "#0f766e", "IVM + DP" = "#1d4ed8",
             "IVM + SMC" = "#ea7317", "IVM + albendazole" = "#7c3aed")

## Classify from the intervention ARM only. Albendazole given as a comparator
## or as a baseline dewormer (Chaccour, Foy) is NOT co-administration, so we do
## not read it from the comparator field.
classify_int <- function(pkg) {
  pkg <- coalesce(pkg, "")
  case_when(
    str_detect(pkg, regex("SMC|Seasonal", TRUE))         ~ "IVM + SMC",
    str_detect(pkg, regex("piperaquine|\\bDP\\b", TRUE)) ~ "IVM + DP",
    str_detect(pkg, regex("Albendazole", TRUE))          ~ "IVM + albendazole",
    TRUE                                                  ~ "IVM alone")
}

## ---- data (as collected) -------------------------------------------------
harm <- read_excel("ivermectin_FINAL.xlsx", sheet = "Harmonised_effects") %>%
  mutate(Coverage_pct = suppressWarnings(as.numeric(Coverage_pct)),
         intervention = factor(classify_int(Intervention_Package),
                               levels = int_levels))

###############################################################################
# A. AFRICA STUDY MAP  (one bubble per study x country; size = n estimates)
###############################################################################
library(maps); library(ggrepel)

# approximate country centroids (lon, lat) for label placement
cent <- tribble(
  ~Country,          ~lon,    ~lat,
  "Burkina Faso",   -1.75,   12.24,
  "Guinea-Bissau", -15.00,   12.00,
  "Kenya",          37.90,    0.20,
  "Mozambique",     35.50,  -18.70,
  "Senegal",       -14.50,   14.50,
  "Tanzania",       34.90,   -6.40,
  "The Gambia",    -15.40,   13.45)

# ONE message: how much evidence each country contributed. Binned for a clean,
# discrete legend (estimate counts are very skewed). Intervention is deliberately
# NOT mapped here because it varies within countries (see the study panel).
ctry <- harm %>%
  group_by(Country) %>%
  summarise(n_est = n(), n_studies = n_distinct(Author_Year), .groups = "drop") %>%
  mutate(region = recode(Country, "The Gambia" = "Gambia"),
         bin = cut(n_est, breaks = c(0, 10, 25, 50, Inf),
                   labels = c("1-10", "11-25", "26-50", "51 or more")))

africa <- map_data("world")
map_df <- left_join(africa, select(ctry, region, bin), by = "region")

lab <- ctry %>% left_join(cent, by = "Country") %>%
  mutate(txt = paste0(Country, "\n", n_studies, " studies · ", n_est, " estimates"))

fills <- c("1-10" = "#cfe8e0", "11-25" = "#84c7b6",
           "26-50" = "#2f9c8a", "51 or more" = "#0f766e")

p_map <- ggplot() +
  geom_polygon(data = map_df, aes(long, lat, group = group, fill = bin),
               colour = "grey80", linewidth = 0.15) +
  scale_fill_manual(values = fills, na.value = "grey94",
                    name = "Estimates contributed", drop = FALSE) +
  geom_point(data = lab, aes(lon, lat), size = 0.5, colour = "grey15") +
  geom_text_repel(data = lab, aes(lon, lat, label = txt),
                  size = 2.4, lineheight = 0.9, box.padding = 0.6,
                  min.segment.length = 0, segment.colour = "grey55",
                  segment.size = 0.2, seed = 1) +
  coord_quickmap(xlim = c(-20, 52), ylim = c(-37, 40)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(legend.position = c(0.16, 0.30),
        legend.key.size = unit(4.5, "mm"),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7.5))

save_lancet(p_map, "tabfig_africa_map", width_mm = 160, height_mm = 160)

###############################################################################
# B. HETEROGENEITY VIEW  (I2 lollipop per outcome; colour = tau2)  [Table 2]
###############################################################################
het <- read_excel("ivermectin_FINAL.xlsx", sheet = "Main_Table2_revised_PI") %>%
  rename(I2 = `I² (%)`, tau2 = `tau²`) %>%
  filter(Outcome != "Other outcome") %>%
  mutate(I2 = as.numeric(I2), tau2 = as.numeric(tau2),
         domain = factor(Domain,
                    levels = c("Malaria epidemiology","Entomological","Safety"))) %>%
  arrange(domain, I2) %>% mutate(Outcome = fct_inorder(Outcome))

p_het <- ggplot(het, aes(I2, Outcome)) +
  geom_vline(xintercept = c(25, 50, 75), linetype = 3, colour = "grey80", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = I2, yend = Outcome), colour = "grey75", linewidth = 0.4) +
  geom_point(aes(fill = tau2), shape = 21, size = 4, colour = "black") +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient(low = "#dbeafe", high = "#b91c1c", name = expression(tau^2)) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(x = expression("Heterogeneity, "*I^2), y = NULL) +
  theme_paper +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(6, "pt"),
        legend.position = "right", legend.key.width = unit(3, "mm"))

save_lancet(p_het, "tabfig_heterogeneity", width_mm = 168, height_mm = 105)

###############################################################################
# C. STUDY-CHARACTERISTICS PANEL  (coverage bar + intervention colour) [Table 1]
###############################################################################
panel <- harm %>%
  group_by(Author_Year, Country) %>%
  summarise(coverage = suppressWarnings(max(Coverage_pct, na.rm = TRUE)),
            n_est = n(), intervention = first(intervention), .groups = "drop") %>%
  mutate(coverage = ifelse(is.finite(coverage), coverage, NA_real_),
         label = paste0(Author_Year, " (", Country, ")")) %>%
  arrange(intervention, coverage) %>% mutate(label = fct_inorder(label))

p_panel <- ggplot(panel, aes(coverage, label, fill = intervention)) +
  geom_col(width = 0.68, na.rm = TRUE) +
  # mark "not reported" coverage with a hollow tick at the axis
  geom_point(data = subset(panel, is.na(coverage)), aes(y = label),
             x = 2, shape = 4, size = 2, colour = "grey40", inherit.aes = FALSE) +
  geom_text(aes(x = ifelse(is.na(coverage), 6, coverage + 2),
                label = ifelse(is.na(coverage),
                               paste0("NR · ", n_est, " est."),
                               paste0(round(coverage), "% · ", n_est, " est."))),
            hjust = 0, size = 2.5, colour = "black") +
  scale_fill_manual(values = int_pal, drop = FALSE, name = NULL) +
  scale_x_continuous(limits = c(0, 118), breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Median ivermectin coverage", y = NULL) +
  theme_paper +
  theme(legend.position = "top", legend.key.size = unit(3.5, "mm"),
        legend.text = element_text(size = 7.5))

save_lancet(p_panel, "tabfig_study_panel", width_mm = 168, height_mm = 95)

message("Done. Wrote tabfig_africa_map, tabfig_heterogeneity, tabfig_study_panel (pdf/tiff/png).")
