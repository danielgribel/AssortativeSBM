library(ggplot2)
library(Hmisc)

sort_data <- function(data_) {
  return(data_[order(data_$likelihood), ])
}

save_plot <- function(plot, filename, w, h) {
  pdf(file = filename, width = w, height = h)
  plot(plot)
  dev.off()
}

data_ssbm <- read.table("/home/dgribel/src/ConstrainedSBM/out/P-SSBM-02.txt", head = TRUE)
data_csbm <- read.table("/home/dgribel/src/ConstrainedSBM/out/P-CSBM-02.txt", head = TRUE)

df_ssbm <- as.data.frame(data_ssbm)
df_csbm <- as.data.frame(data_csbm)

top <- 100

df <- data.frame()

u <- c(.2, .3, .4, .5, .6, .7, .8, .9, 1.)

for(i in 1:length(u)) {
	  x_ssbm <- sort_data(df_ssbm[df_ssbm$win == u[i], ])
    x_csbm <- sort_data(df_csbm[df_csbm$win == u[i], ])
    df <- rbind(df, x_ssbm[1:top, ], x_csbm[1:top, ])
}

lines <- ggplot(df, aes(x = win, y = nmi, group = model, color = model)) +
    # geom_smooth(
    #     aes(fill = model, linetype = model),
    #     level = 0.99,
    #     color = "gray6",
    #     size = 0.5,
    #     method = "auto"
    # ) + 
    geom_point(aes(shape = model), size = 2.5, alpha = 0.5, position = position_dodge(0.02)) +
    stat_summary(fun.data = "mean_cl_normal",
               geom = 'smooth', se = TRUE) +
    xlab(expression(omega["in"])) + 
    ylab("NMI") + 
    theme(legend.position = "top") +
    scale_x_continuous(breaks = round(seq(min(df$win), max(df$win), by = 0.05), 2)) +
    scale_y_continuous(breaks = round(seq(min(df$nmi), max(df$nmi), by = 0.1), 1)) +
    scale_colour_grey(start = .0, end = .4) + 
      scale_fill_grey(start = .0, end = .9) +
      theme_bw()


### Boxplot PPM

boxplot <- ggplot(df, aes(x = reorder(win, win, FUN = median), y = nmi, fill = factor(model))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), color = "grey30", size = 0.3) +
  # scale_color_manual(values = c("grey80", "grey40")) +
  # coord_flip() + 
  scale_fill_manual(values = c("grey80", "grey40"), name = "", labels = c("Standard DC-SBM", "Assortative DC-SBM")) +
  xlab(expression(omega["in"])) + 
  ylab("NMI") +
  # geom_point(position = position_jitterdodge(), alpha = 0.1) +
  # scale_color_manual(values = c("grey80", "grey40")) +
  # scale_colour_grey(start = .9, end = .5) + 
  # scale_fill_grey(start = .9, end = .5) +
  theme_light() + 
  theme(legend.position = "top")

plotName <- paste("SmoothLine", toString(top), sep = "-")
plotName <- paste(plotName, "pdf", sep = ".")

# save_plot(boxplot, plotName, 5, 4)