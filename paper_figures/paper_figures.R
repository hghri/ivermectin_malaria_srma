###############################################################################
# Ivermectin–malaria SRMA — PUBLICATION FIGURES (reproduces paper Figures 1–3)
#
# Style: clean, journal-ready. Figures carry a SHORT title + legend only; all
#        explanatory text belongs in the manuscript figure legend, not the image.
# Export: vector PDF (Lancet-preferred) + 500 dpi LZW TIFF + PNG preview.
# Data  : ivermectin_FINAL.xlsx, used exactly as collected (nothing re-estimated).
#
# Packages:
#   install.packages(c("readxl","dplyr","stringr","ggplot2","forcats","scales",
#                      "patchwork","DiagrammeR","DiagrammeRsvg","rsvg","magick"))
###############################################################################

library(readxl); library(dplyr); library(stringr); library(ggplot2)
library(forcats); library(scales); library(patchwork)

## ---- Lancet export helper (vector PDF + high-res TIFF + PNG preview) ------
save_lancet <- function(plot, name, width_mm, height_mm) {
  ## Base "pdf" device: portable everywhere (Cairo is broken on some Windows R
  ## builds yet still reports as available). Swap in cairo_pdf for embedded fonts.
  ggsave(paste0(name, ".pdf"),  plot, device = "pdf",
         width = width_mm, height = height_mm, units = "mm")
  ggsave(paste0(name, ".tiff"), plot, device = "tiff",
         width = width_mm, height = height_mm, units = "mm",
         dpi = 500, compression = "lzw", bg = "white")
  ggsave(paste0(name, ".png"),  plot, device = "png",
         width = width_mm, height = height_mm, units = "mm",
         dpi = 300, bg = "white")
  invisible(NULL)
}

## Minimal journal theme — no subtitle, no caption.
theme_paper <- theme_classic(base_size = 9) +
  theme(plot.title   = element_text(face = "bold", size = 10),
        axis.title   = element_text(size = 9),
        axis.text    = element_text(colour = "black"),
        legend.position = "none")

###############################################################################
# FIGURE 1 — PRISMA 2020 flow diagram  (DiagrammeR / Graphviz)
###############################################################################
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(magick)

prisma_dot <- '
digraph PRISMA {
  graph [layout = dot, fontsize = 10, fontname = Helvetica, nodesep = 0.35, ranksep = 0.45]
  node  [shape = box, fontname = Helvetica, fontsize = 9, width = 2.6, style = filled,
         fillcolor = "#eef2f7", color = "#334155"]

  id  [label = "Records identified (n=289)\\lDatabases (n=206)\\lRegisters/other sources (n=83)\\l"]
  dup [label = "Records removed before screening:\\lDuplicates (n=75)\\l"]
  scr [label = "Records screened (n=214)"]
  exs [label = "Records excluded on\\ltitle/abstract (n=167)"]
  ret [label = "Reports sought for retrieval (n=47)"]
  nret[label = "Reports not retrieved (n=0)"]
  elig[label = "Reports assessed for eligibility (n=47)"]
  exf [label = "Reports excluded (n=35):\\lNo comparative data (n=9)\\lLab-only mosquito studies (n=7)\\lModelling-only (n=5)\\lPharmacokinetic-only (n=4)\\lInsufficient outcome/variance (n=6)\\lDuplicate/overlapping datasets (n=2)\\lNon-African settings (n=2)\\l"]
  inc [label = "Studies in quantitative synthesis (n=12)\\lAdditional narrative records (n=5)\\l", fillcolor = "#d9ede8"]

  id -> dup [minlen=2]; id -> scr;
  scr -> exs [minlen=2]; scr -> ret;
  ret -> nret [minlen=2]; ret -> elig;
  elig -> exf [minlen=2]; elig -> inc;

  { rank = same; id;  dup  }
  { rank = same; scr; exs  }
  { rank = same; ret; nret }
  { rank = same; elig; exf }
}'

svg <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(prisma_dot))
rsvg::rsvg_pdf(charToRaw(svg), "figure1_prisma.pdf")
rsvg::rsvg_png(charToRaw(svg), "figure1_prisma.png", width = 1800)   # ~300 dpi at 150 mm
# rsvg has no TIFF writer -> rasterise the SVG at high DPI with magick and write LZW TIFF
tmp_svg <- tempfile(fileext = ".svg"); writeLines(svg, tmp_svg)
magick::image_write(magick::image_read(tmp_svg, density = 500),
                    "figure1_prisma.tiff", format = "tiff", compression = "lzw")

###############################################################################
# Shared: read the study's published primary estimates
###############################################################################
tab <- read_excel("ivermectin_FINAL.xlsx", sheet = "Main_Table2_revised_PI") %>%
  rename(ratio_ci = `Ratio (95% CI)`, PI = `Prediction interval`) %>%
  mutate(across(c(ratio_ci, PI), ~ str_replace_all(., "[–—]", "-")),
         ratio = as.numeric(str_extract(ratio_ci, "^[0-9.]+")),
         lo    = as.numeric(str_match(ratio_ci, "\\(([0-9.]+)-")[,2]),
         hi    = as.numeric(str_match(ratio_ci, "-([0-9.]+)\\)")[,2]),
         pi_lo = as.numeric(str_match(PI, "^([0-9.]+)-")[,2]),
         pi_hi = as.numeric(str_match(PI, "-([0-9.]+)$")[,2]),
         domain = factor(Domain,
                    levels = c("Malaria epidemiology","Entomological","Safety"))) %>%
  filter(Outcome != "Other outcome")

###############################################################################
# FIGURE 2 — Primary outcome-specific REML forest plot
#   Point = pooled ratio; thick bar = 95% CI; thin bar = 95% prediction interval.
###############################################################################
f2 <- tab %>% arrange(domain, ratio) %>%
  mutate(row = fct_inorder(Outcome),
         txt = paste0(sprintf("%.2f", ratio), " (", sprintf("%.2f", lo),
                      "-", sprintf("%.2f", hi), ")"))

xb <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 4)   # axis wide enough to show every prediction interval

## LEFT panel: the forest (95% CI thick bar, prediction interval thin bar)
p_forest <- ggplot(f2, aes(y = row)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_linerange(aes(xmin = pi_lo, xmax = pi_hi),
                 linewidth = 0.4, colour = "grey60", na.rm = TRUE) +          # prediction interval
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 1.1, colour = "black") + # 95% CI
  geom_point(aes(x = ratio), shape = 15, size = 2, colour = "black") +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = xb, labels = as.character(xb)) +
  coord_cartesian(xlim = c(0.03, 4)) +
  labs(x = "Ratio (<1 favours ivermectin)", y = NULL) +
  theme_paper +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(6, "pt"))

## RIGHT panel: the estimate column, in its own panel so it is never clipped
p_txt <- ggplot(f2, aes(y = row, x = 0)) +
  geom_text(aes(label = txt), hjust = 0, size = 2.4) +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void(base_size = 9) +
  theme(strip.text = element_blank(), panel.spacing = unit(6, "pt"),
        plot.margin = margin(5, 6, 5, 2))

p2 <- p_forest + p_txt + plot_layout(widths = c(3.1, 1.2))
save_lancet(p2, "figure2_forest", width_mm = 185, height_mm = 120)

###############################################################################
# FIGURE 3 — Bayesian sensitivity: posterior probability that ratio < 1
###############################################################################
bayes <- read_excel("ivermectin_FINAL.xlsx", sheet = "Bayesian_sensitivity") %>%
  filter(str_detect(Group, ":")) %>%
  transmute(domain  = str_squish(str_extract(Group, "^[^:]+")),
            outcome = str_squish(str_remove(Group, "^[^:]+:")),
            p = as.numeric(Pr_ratio_lt_1)) %>%
  filter(outcome != "Other outcome") %>%                                   # match Figure 2
  mutate(domain = factor(domain,
           levels = c("Malaria epidemiology", "Entomological", "Safety"))) %>%
  arrange(domain, p) %>%
  mutate(outcome = fct_inorder(outcome), strong = p >= 0.95)

p3 <- ggplot(bayes, aes(p, outcome)) +
  geom_vline(xintercept = 0.5,  linetype = 3, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0.95, linetype = 2, colour = "grey30", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = p, yend = outcome), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(fill = strong), shape = 21, size = 3, colour = "black") +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y", switch = "y") +  # grouped like Fig 2
  scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, .95),
                     labels = percent_format(accuracy = 1)) +
  labs(x = "Posterior probability ratio < 1", y = NULL) +
  theme_paper +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(6, "pt"))

save_lancet(p3, "figure3_bayes", width_mm = 168, height_mm = 105)

message("Wrote figure1_prisma, figure2_forest, figure3_bayes (pdf/tiff/png).")
