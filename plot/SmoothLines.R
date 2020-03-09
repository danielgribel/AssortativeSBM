library(ggplot2)
library(Hmisc)

sort_data <- function(data_) {
  return(data_[order(data_$likelihood),])
}

data_ssbm <- read.table("/home/dgribel/src/ConstrainedSBM/out/X-SSBM.txt", head = TRUE)
data_csbm <- read.table("/home/dgribel/src/ConstrainedSBM/out/X-CSBM.txt", head = TRUE)

df_ssbm <- as.data.frame(data_ssbm)
df_csbm <- as.data.frame(data_csbm)

top <- 50
u <- unique(df_ssbm$win)
frame_ssbm <- data.frame()
frame_csbm <- data.frame()

u2 <- c(0.25, 0.3, 0.35, 0.4)

for(i in 1:length(u2)) {
	x_ssbm <- sort_data(df_ssbm[df_ssbm$win == u2[i],])
    x_csbm <- sort_data(df_csbm[df_csbm$win == u2[i],])
    frame_ssbm <- rbind(frame_ssbm, x_ssbm[1:top,])
    frame_csbm <- rbind(frame_csbm, x_csbm[1:top,])
}

# df <- rbind(df_ssbm, df_csbm)
df <- rbind(frame_ssbm, frame_csbm)

plotName <- paste("SmoothLine", toString(top), sep = "-")
plotName <- paste(plotName, "pdf", sep = ".")

# pdf(file = plotName, width = 8, height = 6)

p <- ggplot(df, aes(x = win, y = nmi, group = model, color = model)) +
    # geom_smooth(
    #     aes(fill = model, linetype = model),
    #     level = 0.99,
    #     color = "gray6",
    #     size = 0.5,
    #     method = "auto"
    # ) + 
    geom_point(alpha = 0.5, position = position_dodge(0.005)) +
    stat_summary(fun.data = "mean_cl_normal",
               geom = 'smooth', se = TRUE) +
    xlab("w_in") + 
    ylab("NMI") + 
    theme(legend.position = "top") +
    scale_x_continuous(breaks = round(seq(min(df$win), max(df$win), by = 0.05), 2)) +
    scale_y_continuous(breaks = round(seq(min(df$nmi), max(df$nmi), by = 0.1), 1))

# plot(p)
# dev.off()
