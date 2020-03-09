library(ggplot2)

filter_data <- function(instance_, model_) {
  return(df[df$win == w_in & df$instance == instance_ & df$model == model_,])
}

sort_data <- function(data_) {
  return(data_[order(data_$likelihood),])
}

# load input file
data_ssbm <- read.table("/home/dgribel/src/ConstrainedSBM/out/Boxplot-SSBM-5050.txt", head = TRUE)
data_csbm <- read.table("/home/dgribel/src/ConstrainedSBM/out/Boxplot-CSBM-5050.txt", head = TRUE)
data_modu <- read.table("/home/dgribel/src/ConstrainedSBM/out/Boxplot-MODU-5050.txt", head = TRUE)

# create dataframes
df_ssbm <- as.data.frame(data_ssbm)
df_csbm <- as.data.frame(data_csbm)
df_modu <- as.data.frame(data_modu)

# append dataframes
df <- rbind(df_ssbm, df_csbm, df_modu)

nb_datasets <- length(unique(df_csbm$instance))

w_in <- 0.2

top  <- 5

init <- 100

df2 <- data.frame()

diff_nmi <- array(0, dim=(nb_datasets))

for(r in 1:nb_datasets) {
  csbm <- data.frame()
  modu <- data.frame()

  instanceId <- init + r 

  # filter data by instance
  ssbm <- filter_data(paste("A", toString(instanceId), sep=""), "SSBM")
  csbm <- filter_data(paste("A", toString(instanceId), sep=""), "CSBM")
  modu <- filter_data(paste("A", toString(instanceId), sep=""), "MODU")
  
  # sort datasets
  ssbm <- sort_data(ssbm)
  csbm <- sort_data(csbm)
  modu <- sort_data(modu)

  df2 <- rbind(df2, ssbm[1:top,], csbm[1:top,], modu[1:top,])

  avg_ssbm <- mean(ssbm[1:top,]$nmi)
  avg_csbm <- mean(csbm[1:top,]$nmi)
  avg_modu <- mean(modu[1:top,]$nmi)

  diff_nmi[r] <- avg_csbm - avg_ssbm
}

plotName <- paste("Boxplot", gsub('\\.', '', toString(w_in)), toString(top), sep = "-")
plotName <- paste(plotName, "pdf", sep = ".")

# pdf(file = plotName, width = 4, height = 4)

box <- ggplot(data = df2, aes(x = instance, y = nmi, fill = factor(model))) +
  geom_boxplot() +
  labs(fill = "Model: ") + 
  xlab("Dataset") + 
  ylab("NMI") + 
  geom_point(position = position_jitterdodge(), alpha = 0.2) +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = round(seq(min(df$nmi), max(df$nmi), by = 0.1), 1))

# plot(box)
# dev.off()

diff_nmi <- sort(diff_nmi, decreasing = TRUE)

dtf1 <- data.frame(Dataset = c(1:length(diff_nmi)), Diff = diff_nmi)
dtf1$colour <- ifelse(dtf1$Diff < 0, "negative","positive")

bar <- ggplot(dtf1, aes(Dataset, Diff, label = "")) +
  geom_text(aes(y = 0, colour = colour))+
  geom_bar(stat = "identity", position = "identity", aes(fill = colour))