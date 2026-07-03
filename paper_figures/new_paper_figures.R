###############################################################################
# Ivermectin–malaria SRMA — NEW PUBLICATION FIGURES (additions to the paper)
#
# 1. Study-nested forest plots (per outcome)        -> newfig_forest_<outcome>
# 2. Meta-regression bubble plots                    -> newfig_metareg_<mod>
# 3. Risk-of-bias (robvis) traffic-light + summary   -> newfig_rob_*
# 4. Evidence-contribution map (studies x outcomes)  -> newfig_evidence_map
#
# Clean, journal-style: short title + legend only, no explanatory captions.
# Dataset used as collected. Export: vector PDF + 500 dpi LZW TIFF + PNG.
#
# install.packages(c("readxl","dplyr","stringr","ggplot2","forcats","scales",
#                    "metafor","robvis"))
###############################################################################

library(readxl); library(dplyr); library(stringr); library(ggplot2)
library(forcats); library(scales); library(metafor)

## ---- export helpers ------------------------------------------------------
save_lancet <- function(plot, name, width_mm, height_mm) {           # for ggplot
  ## Base "pdf" device for portability (Cairo unreliable on some Windows builds).
  ggsave(paste0(name, ".pdf"),  plot, device = "pdf",
         width = width_mm, height = height_mm, units = "mm")
  ggsave(paste0(name, ".tiff"), plot, device = "tiff", width = width_mm,
         height = height_mm, units = "mm", dpi = 500, compression = "lzw", bg = "white")
  ggsave(paste0(name, ".png"),  plot, device = "png", width = width_mm,
         height = height_mm, units = "mm", dpi = 300, bg = "white"); invisible(NULL)
}
save_base <- function(draw, name, width_mm, height_mm) {             # for base graphics (metafor)
  pdf(paste0(name, ".pdf"), width = width_mm/25.4, height = height_mm/25.4)
  draw(); dev.off()
  tiff(paste0(name, ".tiff"), width = width_mm, height = height_mm, units = "mm",
       res = 500, compression = "lzw"); draw(); dev.off()
  png(paste0(name, ".png"), width = width_mm, height = height_mm, units = "mm",
      res = 300); draw(); dev.off(); invisible(NULL)
}

theme_paper <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text = element_text(colour = "black"))

## ---- data (as collected) -------------------------------------------------
dat <- read_excel("ivermectin_FINAL.xlsx", sheet = "Poolable_post_MDA") %>%
  filter(!is.na(yi), !is.na(vi)) %>%
  mutate(yi = as.numeric(yi), vi = as.numeric(vi),
         Coverage_pct = suppressWarnings(as.numeric(Coverage_pct)))

###############################################################################
# 1. STUDY-NESTED FOREST PLOTS  (individual estimates within each outcome)
#    One file per outcome; pooled REML diamond added at the foot.
###############################################################################
forest_outcomes <- c("Mosquito mortality", "Entomological inoculation rate",
                     "Malaria incidence/infection", "Mosquito survival",
                     "Vector density", "Any adverse events")

for (oc in forest_outcomes) {
  d <- dat[dat$Outcome_measure == oc, ]
  d <- d[order(d$Author_Year, d$yi), ]
  m <- rma(yi, vi, data = d, method = "REML")
  h <- max(90, 12 + 4.2 * nrow(d))                       # scale height to n rows (mm)
  save_base(function() {
    forest(m, slab = d$Author_Year, atransf = exp,
           at = log(c(0.1, 0.25, 0.5, 1, 2, 4)),
           xlab = "Ratio (<1 favours ivermectin)", header = TRUE,
           cex = 0.62, mlab = "Pooled (REML)", col = "grey20")
  }, paste0("newfig_forest_", str_replace_all(tolower(oc), "[^a-z]+", "_")),
     width_mm = 168, height_mm = h)
}

###############################################################################
# 2. META-REGRESSION BUBBLE PLOTS
#    Bubble = one estimate (area ~ meta-analytic weight); line = fitted effect.
###############################################################################
bubble_plot <- function(df, mod, xlab, title) {
  df <- df[is.finite(df[[mod]]), ]
  df$w <- 1 / df$vi
  m  <- rma(yi, vi, mods = as.formula(paste("~", mod)), data = df, method = "REML")
  xs <- seq(min(df[[mod]]), max(df[[mod]]), length.out = 100)
  pr <- predict(m, newmods = xs)
  line <- data.frame(x = xs, ratio = exp(pr$pred), lo = exp(pr$ci.lb), hi = exp(pr$ci.ub))
  ggplot(df, aes(.data[[mod]], exp(yi))) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(data = line, aes(x, ymin = lo, ymax = hi), inherit.aes = FALSE,
                fill = "grey85", alpha = 0.6) +
    geom_line(data = line, aes(x, ratio), inherit.aes = FALSE, linewidth = 0.6) +
    geom_point(aes(size = w, colour = Author_Year), alpha = 0.75) +
    scale_y_log10() + scale_size_area(max_size = 7, guide = "none") +
    labs(x = xlab, y = "Ratio (log scale)", colour = "Study") +
    theme_paper + theme(legend.position = "bottom",
                        legend.key.size = unit(3, "mm"),
                        legend.text = element_text(size = 6))
}

# Example: entomological effect vs ivermectin coverage
ento <- dat[dat$Domain == "Entomological", ]
p_cov <- bubble_plot(ento, "Coverage_pct",
                     "Ivermectin coverage (%)",
                     "Entomological effect vs coverage")
save_lancet(p_cov, "newfig_metareg_coverage", width_mm = 130, height_mm = 95)

# By co-administered drug (categorical moderator) — boxplot-style bubbles.
# Harmonise label variants of the same co-intervention for a clean axis.
ento2 <- ento %>%
  mutate(coadmin = case_when(
    str_detect(CoAdmin_Drug, regex("piperaquine|\\bDP\\b", ignore_case = TRUE)) ~ "DP",
    str_detect(CoAdmin_Drug, regex("SMC", ignore_case = TRUE))                  ~ "SMC",
    str_detect(CoAdmin_Drug, regex("Albendazole", ignore_case = TRUE))          ~ "Albendazole",
    CoAdmin_Drug == "/"                                                          ~ "Ivermectin alone",
    TRUE ~ CoAdmin_Drug)) %>%
  filter(!is.na(coadmin))
p_coadmin <- ggplot(ento2, aes(fct_reorder(coadmin, yi), exp(yi))) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(size = 1/vi, colour = Author_Year), alpha = 0.7,
             position = position_jitter(width = 0.12, height = 0)) +
  scale_y_log10() + scale_size_area(max_size = 6, guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Ratio (log scale)", colour = "Study") +
  theme_paper + theme(legend.position = "none")
save_lancet(p_coadmin, "newfig_metareg_coadmin", width_mm = 130, height_mm = 80)

###############################################################################
# 3. RISK OF BIAS (robvis)
#    NB: requires an author-completed RoB2/ROBINS-I assessment. Template below
#    with placeholder judgements — REPLACE with your consensus scoring, then plot.
###############################################################################
library(robvis)

# ---- TEMPLATE: one row per study; fill D1..D5 with "Low"/"Some concerns"/"High"
rob <- tribble(
  ~Study,           ~D1, ~D2, ~D3, ~D4, ~D5, ~Overall, ~Weight,
  "Sylla 2010",     NA,  NA,  NA,  NA,  NA,  NA,        1,
  "Alout 2014",     NA,  NA,  NA,  NA,  NA,  NA,        1,
  "Derua 2015",     NA,  NA,  NA,  NA,  NA,  NA,        1,
  "Foy 2019",       NA,  NA,  NA,  NA,  NA,  NA,        1,
  "Smit 2018",      NA,  NA,  NA,  NA,  NA,  NA,        1,
  "Smit 2019",      NA,  NA,  NA,  NA,  NA,  NA,        1,
  "Dabira 2022",    NA,  NA,  NA,  NA,  NA,  NA,        1,
  "Soumare 2022",   NA,  NA,  NA,  NA,  NA,  NA,        1,
  "Chaccour 2025",  NA,  NA,  NA,  NA,  NA,  NA,        1,
  "Pretorius 2024", NA,  NA,  NA,  NA,  NA,  NA,        1,
  "Some 2025",      NA,  NA,  NA,  NA,  NA,  NA,        1
)
write.csv(rob, "rob2_scoring_template.csv", row.names = FALSE)   # fill this in

# Column names robvis expects for a RoB2 tool:
names(rob) <- c("Study",
  "D1: Randomisation process", "D2: Deviations from intervention",
  "D3: Missing outcome data",  "D4: Measurement of the outcome",
  "D5: Selection of reported result", "Overall", "Weight")

## Once `rob` is completed (no NAs), uncomment to render:
# tl <- rob_traffic_light(data = rob, tool = "ROB2", psize = 8)
# sm <- rob_summary(data = rob, tool = "ROB2")
# ggsave("newfig_rob_trafficlight.pdf", tl, device = cairo_pdf, width = 150, height = 200, units = "mm")
# ggsave("newfig_rob_trafficlight.tiff", tl, device = "tiff", width = 150, height = 200,
#        units = "mm", dpi = 500, compression = "lzw")
# ggsave("newfig_rob_summary.pdf", sm, device = cairo_pdf, width = 168, height = 60, units = "mm")
# ggsave("newfig_rob_summary.tiff", sm, device = "tiff", width = 168, height = 60,
#        units = "mm", dpi = 500, compression = "lzw")
message("robvis: fill rob2_scoring_template.csv, then uncomment the render block.")

###############################################################################
# 4. EVIDENCE-CONTRIBUTION MAP  (studies x outcomes; cell = number of estimates)
###############################################################################
dom_order <- c("Malaria epidemiology", "Entomological", "Safety")
cnt <- dat %>%
  count(Author_Year, Domain, Outcome_measure, name = "n") %>%
  mutate(Domain = factor(Domain, levels = dom_order))

study_order <- dat %>% count(Author_Year, name = "tot") %>% arrange(tot) %>% pull(Author_Year)

p_map <- cnt %>%
  mutate(Author_Year = factor(Author_Year, levels = study_order),
         Outcome_measure = reorder(Outcome_measure, as.integer(Domain))) %>%
  ggplot(aes(Outcome_measure, Author_Year, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n), size = 2.6,
            colour = ifelse(cnt$n > 20, "white", "black")) +
  facet_grid(. ~ Domain, scales = "free_x", space = "free_x") +
  scale_fill_gradient(low = "#dbeafe", high = "#0f766e", name = "Estimates") +
  labs(x = NULL, y = NULL) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        legend.key.size = unit(4, "mm"),
        plot.margin = margin(5, 6, 10, 5))
save_lancet(p_map, "newfig_evidence_map", width_mm = 180, height_mm = 110)

message("Done. Wrote nested forests, meta-regression bubbles, evidence map (pdf/tiff/png), robvis template.")
