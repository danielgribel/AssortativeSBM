library(ggplot2)

filter_data <- function(instance_, model_) {
  return(df[df$win == w_in & df$instance == instance_ & df$model == model_,])
}

sort_data <- function(data_) {
  return(data_[order(data_$likelihood),])
}

save_plot <- function(plot, filename, w, h) {
  pdf(file = filename, width = w, height = h)
  plot(plot)
  dev.off()
}

create_barplot <- function(df) {
  df2 <- data.frame(Dataset = c(1:length(df)), Diff = df)
  bar <- ggplot(df2, aes(Dataset, Diff, label = "")) +
      geom_text(aes(y = 0)) +
      geom_bar(stat = "identity", position = "identity", fill = "grey30") +
      scale_x_continuous(breaks = round(seq(1, nb_datasets, by = 1), 1)) +
      #  scale_fill_grey(start = .0, end = .9) +
      theme_light() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 0)) +
      ylab("Difference in NMI") +
      xlab("Datasets")
  return(bar)
}

# load input file
data_ssbm <- read.table("/home/dgribel/src/ConstrainedSBM/out/Boxplot-SSBM05-50-50.txt", head = TRUE)
data_csbm <- read.table("/home/dgribel/src/ConstrainedSBM/out/Boxplot-CSBM05-50-50.txt", head = TRUE)
data_modu <- read.table("/home/dgribel/src/ConstrainedSBM/out/Boxplot-MODU05-50-50.txt", head = TRUE)

# create dataframes
df_ssbm <- as.data.frame(data_ssbm)
df_csbm <- as.data.frame(data_csbm)
df_modu <- as.data.frame(data_modu)

# append dataframes
df <- rbind(df_ssbm, df_csbm, df_modu)

nb_datasets <- length(unique(df_csbm$instance))

w_in <- 0.5

top  <- 50

init <- 100

df2 <- data.frame()

diff_nmi_ssbm <- array(0, dim = (nb_datasets))
diff_nmi_modu <- array(0, dim = (nb_datasets))

avg_nmi_ssbm  <- array(0, dim = (nb_datasets))
avg_nmi_csbm  <- array(0, dim = (nb_datasets))
avg_nmi_modu  <- array(0, dim = (nb_datasets))

avg_ass_ssbm  <- array(0, dim = (nb_datasets))
avg_ass_csbm  <- array(0, dim = (nb_datasets))
avg_ass_modu  <- array(0, dim = (nb_datasets))

for(r in 1:nb_datasets) {
  ssbm <- data.frame()
  csbm <- data.frame()
  modu <- data.frame()

  instanceId <- init + r 

  # filter data by instancemydata
  ssbm <- filter_data(paste("A", toString(instanceId), sep = ""), "SSBM")
  csbm <- filter_data(paste("A", toString(instanceId), sep = ""), "CSBM")
  modu <- filter_data(paste("A", toString(instanceId), sep = ""), "MODU")
  
  # sort datasets
  ssbm <- sort_data(ssbm)
  csbm <- sort_data(csbm)
  modu <- sort_data(modu)

  df2 <- rbind(df2, ssbm[1:top,], csbm[1:top,], modu[1:top,])

  avg_nmi_ssbm[r] <- mean(ssbm[1:top, ]$nmi)
  avg_nmi_csbm[r] <- mean(csbm[1:top, ]$nmi)
  avg_nmi_modu[r] <- mean(modu[1:top, ]$nmi)

  avg_ass_ssbm[r] <- mean(ssbm[1:top, ]$assortativity)
  avg_ass_csbm[r] <- mean(csbm[1:top, ]$assortativity)
  avg_ass_modu[r] <- mean(modu[1:top, ]$assortativity)

  diff_nmi_ssbm[r] <- avg_nmi_csbm[r] - avg_nmi_ssbm[r]
  diff_nmi_modu[r] <- avg_nmi_csbm[r] - avg_nmi_modu[r]
}

k <- 4

### Histogram plot

histo_ssbm <- data.frame(Model = rep("SSBM", nb_datasets), Ass = k*avg_ass_ssbm)
histo_csbm <- data.frame(Model = rep("CSBM", nb_datasets), Ass = k*avg_ass_csbm)

histo <- data.frame(rbind(histo_ssbm, histo_csbm))

compare_mean <- aggregate(. ~ Model, histo, mean)

histogram <- ggplot(histo, aes(x = Ass, fill = Model, color = Model)) +
    geom_histogram(alpha = 0.2, binwidth = 0.25, position = "identity") +
    geom_vline(data = compare_mean, aes(xintercept = Ass, color = Model), linetype = "dashed", size = 1)

### Density plot

density <- ggplot(histo, aes(x = Ass, fill = Model)) +
    # geom_density(alpha = .75) +
    # scale_fill_grey(start = .7, end = .3, name = "", labels = c("Standard DC-SBM", "Assortative DC-SBM")) +
    geom_density() +
    scale_fill_manual(values = c("grey80", "grey40"), name = "", labels = c("Standard DC-SBM", "Assortative DC-SBM")) +
    theme_bw() +
    theme(legend.position = "top") +
    xlab("Number of assortative clusters") +
    ylab("Density") +
    scale_x_continuous(breaks = round(seq(0, 4, by = 0.5), 1))

### Boxplot

boxplot <- ggplot(data = df2[df2$model %in% c("CSBM", "SSBM"),],
              aes(x = reorder(instance, -nmi, FUN = median), y = nmi, fill = factor(model))) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c("grey80", "grey40"), name = "", labels = c("Standard DC-SBM", "Assortative DC-SBM")) +
    theme_light() +
    xlab("Datasets") + 
    ylab("NMI") + 
    # coord_flip() + 
    # geom_point(position = position_jitterdodge(), alpha = 0.1) +
    scale_y_continuous(breaks = round(seq(min(df$nmi), max(df$nmi), by = 0.1), 1)) +
    theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 0))

### Bar plot

diff_nmi_ssbm <- sort(diff_nmi_ssbm, decreasing = TRUE)
diff_nmi_modu <- sort(diff_nmi_modu, decreasing = TRUE)

bar_ssbm <- create_barplot(diff_nmi_ssbm)
bar_modu <- create_barplot(diff_nmi_modu)

### Lines plot

lines <- ggplot(data = df2[df2$model %in% c("CSBM", "SSBM"),],
            aes(x = reorder(instance, -nmi, FUN = median), y = nmi, group = model, color = model)) +    
    geom_point(aes(shape = model), size = 2.5, alpha = 0.9, position = position_dodge(0.55)) +
    # stat_summary(fun.data = "mean_cl_normal", geom = 'smooth', se = TRUE) +
    xlab("Datasets") + 
    ylab("NMI") + 
    theme(legend.position = "top")
    # scale_x_continuous(breaks = round(seq(min(df$win), max(df$win), by = 0.05), 2)) +
    # scale_y_continuous(breaks = round(seq(min(df$nmi), max(df$nmi), by = 0.1), 1))

lines <- lines + scale_colour_grey(start = .7, end = .5) + 
      scale_fill_grey(start = .7, end = .5) +
      theme_bw()

### Save plot to file

# save_plot(boxplot, "Boxplot-05-50.pdf", 10, 4)
# save_plot(density, "Density.pdf", 4, 4)
# save_plot(bar_nmi, "Bar-SSBM-05-5.pdf", 7, 4)
# save_plot(bar_nmi, "Bar-MODU-05-5.pdf", 7, 4)