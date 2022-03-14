library(ggplot2)
library(mclust)

# getting optimal number of clusters - getting model 
# https://mclust-org.github.io/mclust/reference/mclustModelNames.html
# https://mclust-org.github.io/mclust/reference/Mclust.html
# https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html

can_data <- read.delim("cannabis_log_ks_values.txt", header = FALSE, 
                       sep = "", dec = ".")
hop_data <- read.delim("hop_log_ks_values.txt", header = FALSE, 
                       sep = "", dec = ".")
hop_vs_can_data <- read.delim("hop_vs_can_log_ks_values.txt", 
                              header = FALSE, sep = "", dec = ".")

# create "double" data frames
can_matrix <- as.matrix(can_data)
can_ks <- as.double(can_matrix)

hop_matrix <- as.matrix(hop_data)
hop_ks <- as.double(hop_matrix)

hop_vs_can_matrix <- as.matrix(hop_vs_can_data)
hop_vs_can_ks <- as.double(hop_vs_can_matrix)


# cannabis
can_ICL <- mclustICL(can_ks)
summary(can_ICL)
plot(can_ICL)

can_ks_cluster <- Mclust(can_ks, modelName="V", G=2)
summary(can_ks_cluster)
can_ks_cluster$modelName
can_ks_cluster$G
can_ks_densityMclust <- densityMclust(can_ks, modelName="V", G=2)
plot(can_ks_densityMclust, what = "density", data = can_ks, breaks = 35)
abline(v = can_ks_cluster$parameters$mean)
plot(can_ks_densityMclust, what = "diagnostic", type = "cdf")
plot(can_ks_densityMclust, what = "diagnostic", type = "qq")

write.table(can_ks_densityMclust$density, "can_log_ks_mclust_densities.txt", 
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)
write.table(can_ks_densityMclust$data, "can_log_ks_mclust_data.txt", 
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)
write.table(can_ks_densityMclust$parameters$mean, "can_log_ks_mclust_means.txt",
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)
write.table(can_ks_densityMclust$parameters$pro, "can_log_ks_mclust_probs.txt",
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)

### hop
hop_ICL <- mclustICL(hop_ks)
summary(hop_ICL)
plot(hop_ICL)

hop_ks_cluster <- Mclust(hop_ks, modelName="V", G=3)
summary(hop_ks_cluster)
hop_ks_cluster$modelName
hop_ks_cluster$G
hop_ks_densityMclust <- densityMclust(hop_ks, modelName="V", G=3)
plot(hop_ks_densityMclust, what = "density", data = hop_ks, breaks = 35)
abline(v = hop_ks_cluster$parameters$mean)
plot(hop_ks_densityMclust, what = "diagnostic", type = "cdf")
plot(hop_ks_densityMclust, what = "diagnostic", type = "qq")


write.table(hop_ks_densityMclust$density, "hop_log_ks_mclust_densities.txt", 
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)
write.table(hop_ks_densityMclust$data, "hop_log_ks_mclust_data.txt", 
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)
write.table(hop_ks_densityMclust$parameters$mean, "hop_log_ks_mclust_means.txt",
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)
write.table(hop_ks_densityMclust$parameters$pro, "hop_log_ks_mclust_probs.txt",
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)


### hop vs can
hop_vs_can_ICL <- mclustICL(hop_vs_can_ks)
summary(hop_vs_can_ICL)
plot(hop_vs_can_ICL)

hop_vs_can_ks_cluster <- Mclust(hop_vs_can_ks, modelName="V", G=3)
summary(hop_vs_can_ks_cluster)
hop_vs_can_ks_cluster$modelName
hop_vs_can_ks_cluster$G
hop_vs_can_ks_densityMclust <- densityMclust(hop_vs_can_ks, modelName="V", G=3)
summary(hop_vs_can_ks_densityMclust)
plot(hop_vs_can_ks_densityMclust, what = "density", data = hop_vs_can_ks)
abline(v = hop_vs_can_ks_cluster$parameters$mean)

write.table(hop_vs_can_ks_densityMclust$density, "hop_vs_can_log_ks_mclust_densities.txt", 
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)
write.table(hop_vs_can_ks_densityMclust$data, "hop_vs_can_log_ks_mclust_data.txt", 
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)
write.table(hop_vs_can_ks_densityMclust$parameters$mean, "hop_vs_can_log_ks_mclust_means.txt",
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)
write.table(hop_vs_can_ks_densityMclust$parameters$pro, "hop_vs_can_log_ks_mclust_probs.txt",
            append = FALSE, sep = " ", dec = ".", row.names = FALSE, 
            col.names = FALSE)