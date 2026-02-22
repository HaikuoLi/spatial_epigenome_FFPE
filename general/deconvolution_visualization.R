##results from standard pipeline of SPOTlight and here we show how we visualize

library(scatterpie)

df <- read.csv("deconvolved_ref1_filtered.csv")

head(df)

ct <- colnames(df)[4:(ncol(df)-1)] 
print(ct)
ct_color <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

options(repr.plot.width = 20, repr.plot.height = 20)
p1 <- ggplot() +
  geom_scatterpie(aes(x=x_coord, y=y_coord),
                  data=df,
                  cols=ct, pie_scale = 0.5)+
  scale_fill_manual(values = ct_color) +theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  ) +
  coord_equal()
p1

keep_top_3 <- function(x) {
  ranks <- rank(-x, ties.method = "min")
  ifelse(ranks <= 3, x, 0)
}

# please change the column index, astro is column 4
df_top3 <- as.data.frame(t(apply(df[, 4:(ncol(df)-1)], 1, keep_top_3)))
df_top3$x_coord<- df$x_coord
df_top3$y_coord<- df$y_coord
head(df_top3)

p2 <- ggplot() +
  geom_scatterpie(aes(x=x_coord, y=y_coord),
                  data=df_top3,
                  cols=ct, pie_scale = 0.5)+
  scale_fill_manual(values = ct_color) +theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  ) +
  coord_equal()
p2

