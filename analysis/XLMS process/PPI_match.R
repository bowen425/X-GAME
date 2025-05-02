library(tidyverse)
library(dplyr)
data <- read.csv(file = "inter_470_PPI.csv")

data <- data %>% filter(.,(data$Gene1 != "" & data$Gene2 != ""))
biogrid <- read.delim('BIOGRID_human_20230723.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
biogrid <- biogrid[,c(3,4,8,9)]
colnames(biogrid) <- c("Gene1", "Gene2", "Publiction", "Throughput")
data <- data%>% as_tibble() %>%
  mutate(pair1 = paste0(data$Gene1,data$Gene2))
data <- data%>% as_tibble() %>%
  mutate(pair2 = paste0(data$Gene2,data$Gene1))
biogrid <- biogrid%>% as_tibble() %>%
  mutate(pair = paste0(biogrid$Gene1,biogrid$Gene2))
list <- cbind(biogrid[,5], biogrid[,3], biogrid[,4])
list_1 <- list[match(data$pair1,list[,"pair"]),]
list_2 <- list[match(data$pair2,list[,"pair"]),]
data <- cbind(data[,1:6], list_1[,2:3], list_2[,2:3])

data$BiogridThroughput <- "NAN"
data$BiogridPublication <- "NAN"
data[is.na(data)] <- "NAN"
for (i in 1:nrow(data)) {
  if (data[i,8] == "Low Throughput") {
    data[i,11] = "Low Throughput"
    data[i,12] = data[i,7]
  }
  else if (data[i,10] == "Low Throughput") {
    data[i,11] = "Low Throughput"
    data[i,12] = data[i,9]
  }
  else if (data[i,8] == "High Throughput"){
    data[i,11] = "High Throughput"
    data[i,12] = data[i,7]
  }
  else if (data[i,10] == "High Throughput"){
    data[i,11] = "High Throughput"
    data[i,12] = data[i,9]
  }
}
data <- data[,-(7:10)]
protein_info <- read.delim('9606_string_protein_info.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
####protein_info1 <- protein_info %>% as_tibble() %>%
####filter(., protein_info$source == "Ensembl_UniProt_AC")
protein_info2 <- protein_info %>% as_tibble() %>%
  filter(., protein_info$source == "BLAST_UniProt_GN_Name")
list <- cbind (protein_info2[,2], protein_info2[,1])
colnames(list) <- c("Gene", "STRING")
list_A=list[match(data$Gene1,list[,"Gene"]),]
list_B=list[match(data$Gene2,list[,"Gene"]),]
data2 <- cbind(data, list_A[,2], list_B[,2])
data2 <- data2%>% as_tibble() %>%
  mutate(pair1 = paste0(data2$`list_A[, 2]`,data2$`list_B[, 2]`))
data2 <- data2%>% as_tibble() %>%
  mutate(pair2 = paste0(data2$`list_B[, 2]`,data2$`list_A[, 2]`))
string <- read.delim('9606_string_protein.txt', sep = ' ', stringsAsFactors = FALSE, check.names = FALSE)
string <- string%>% as_tibble() %>%
  mutate(pair = paste0(string$protein1,string$protein2))
list <- cbind(string[,4], string[,3])
colnames(list) <- c("pair", "score")
list_1 <- list[match(data2$pair1,list[,"pair"]),]
list_2 <- list[match(data2$pair2,list[,"pair"]),]
data3 <- cbind(data2[,1:8], list_1[,2], list_2[,2])
data3 <- data3 %>% as_tibble() %>%
  mutate(StringScore = 0)
data3[is.na(data3)] = 0
for (i in 1:nrow(data3)) {
  if (data3[i,9] != 0) {
    data3[i,11] = data3[i,9]
  }
  else if (data3[i,10] != 0) {
    data3[i,11] = data3[i,10]
  }
}
data3 <- data3[,-(9:10)]
write.csv(data3, file = "inter_470_PPI_match_20250405.csv", row.names = F)


