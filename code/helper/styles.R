# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Colors
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.inst = c("ggthemes", "scales") %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
colors_use.10 = ggthemes::tableau_color_pal("Tableau 10")(10)
colors_use.20 = ggthemes::tableau_color_pal("Tableau 20")(20)
colors_stata =  ggthemes::stata_pal("s2color")(15)
cell.cylce.col = c("#004488", "#DDAA33", "#BB5566", "black")
names(cell.cylce.col) = c("G1", "S", "G2M", "-")

cont.col = c(
  "#412856", "#412856", "#344174", "#386293", "#5686AC",
  "#85A9C2", "#B5C0C7", "#f0f0f0", "#D8A88D", "#CC8864",
  "#B2613C", "#903A22", "#721F1E", "#5E1529", "#5E1529"
)

ct.col = c(
  "Plasma(blast)" = "#7C2529",
  "Plasmablast" = "#7C2529",
  "Plasma cell" = "#7C2529",
  "B-Cell" = "#E18A8D",
  "NK" = "#BC8400",
  "NK_CD56bright" = "#BC8400",
  "CD56 bright NK" = "#BC8400",
  "CD4 T-Cell" = "#e8d725",
  "CD8 T-Cell" = "#AFA10D",
  "gdT-Cell" = "#20581C",
  "T-Cell (cycling)" = "#CC3311",
  "Mono CD14" = "#93aeba",
  "CD14 Mono" = "#93aeba",
  "Mono CD16" = "#4B859F",
  "CD16 Mono" = "#4B859F",
  "cDC" = "#9C9BDB",
  "pDC" = "#194573",
  "other DC" = "#D1BBD7",
  "Macrophage" = "black",
  "Erythrocyte" = "#B281A6",
  "Platelet" = "#AA4499",
  "Progenitor" = "#555555",
  "Other" = "black",
  "Cycling" = "grey"
)

cd8.col = c("#0077BB", "#33BBEE", "#009988", "#EE7733", "#EE6677", "#EE3377", "#997700", "#994455")
names(cd8.col) = c("CD8.NaiveLike", "CD8.CM", "CD8.EM", "CD8.TEMRA", "CD8.TPEX", "CD8.TEX", "CD8.MAIT", "CD8.Cycling")
cd4.col = c("#b7d2e0", "#da6f6f", "#72b28a", "#e5bfaf", "#aca6e0", "#f5d39f", "#fdbfd4", "#994455")
names(cd4.col) = c("CD4.NaiveLike", "CD4.CTL_EOMES", "CD4.CTL_GNLY", "CD4.CTL_Exh", "CD4.Tfh", "CD4.Th17", "CD4.Treg", "CD4.Cycling")
ct.misc = c("#BBBBBB", "#555555", "#BBBBBB", "#555555", "#332288", "#DDDDDD", "#9e4f6c", "#4f2535", "#555555", "#CCDDAA", "#994455", "#DDDDDD")
names(ct.misc) = c("CD4.other", "CD8.other", "CD4", "CD8", "NaiveLike", "NE", "S", "G2M", "T_myeloid", "gdT", "Cycling", "Not Estimable")
til.col = c(cd8.col, cd4.col, ct.misc)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ggplot theme
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mytheme = function(base_size = 8, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      axis.ticks = element_line(colour = "black"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(size = rel(1), colour = "black"),
      strip.text.y = element_text(size = rel(1), colour = "black"),
      strip.text = element_text(size = rel(1), colour = "black"),
      axis.text = element_text(size = rel(1), colour = "black"),
      axis.title = element_text(size = rel(1), colour = "black"),
      legend.title = element_text(colour = "black", size = rel(1)),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = .3),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(hjust = 0, face = "plain", colour = "black", size = rel(1)),
      plot.subtitle = element_text(colour = "black", size = rel(.85))
    )
}

mytheme_grid = function(base_size = 8, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      axis.ticks = element_line(colour = "black"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(size = rel(1), colour = "black"),
      strip.text.y = element_text(size = rel(1), colour = "black"),
      strip.text = element_text(size = rel(1), colour = "black"),
      axis.text = element_text(size = rel(1), colour = "black"),
      axis.title = element_text(size = rel(1), colour = "black"),
      legend.title = element_text(colour = "black", size = rel(1)),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = .3),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(hjust = 0, face = "plain", colour = "black", size = rel(1)),
      plot.subtitle = element_text(colour = "black", size = rel(.85))
    )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Make grandient for ggplot (background)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
make_gradient <- function(deg = 45, n = 100, cols = blues9) {

  .cran_packages = c("grid", "ggplot2","RColorBrewer")
  .inst = .cran_packages %in% installed.packages()
  if (any(!.inst)) {
    install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
  }
  for (pack in .cran_packages) {
    suppressMessages(library(
      pack,
      quietly = TRUE,
      verbose = FALSE,
      character.only = TRUE
    ))
  }

  cols <- colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(
    data = rep(seq(0, 1, length.out = n) * cos(rad), n),
    byrow = TRUE,
    ncol = n
  ) +
    matrix(
      data = rep(seq(0, 1, length.out = n) * sin(rad), n),
      byrow = FALSE,
      ncol = n
    )
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(
    image = mat,
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = TRUE
  )
}
