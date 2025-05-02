library(tidyverse)
library(dplyr)
data <- read.csv(file = "PPIs.csv")
list <- read.delim('uniprot-gene-human.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
list <- list %>% separate(`Gene Names`,into="Gene Names",sep=" ")
list_A=list[match(data$Protein1,list[,"Entry"]),]
list_B=list[match(data$Protein2,list[,"Entry"]),]
data <- cbind(data[,1:2], list_A[,2], list_B[,2], data[,3])
colnames(data) <- c("Protein1", "Protein2", "Gene1", "Gene2", "Context")
data <- data %>% filter(.,(data$Gene1 != "" & data$Gene2 != ""))

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
data3 <- cbind(data2[,1:5], list_1[,2], list_2[,2])
data3 <- data3 %>% as_tibble() %>%
  mutate(StringScore = 0)
data3[is.na(data3)] = 0
for (i in 1:nrow(data3)) {
  if (data3[i,6] != 0) {
    data3[i,8] = data3[i,6]
  }
  else if (data3[i,7] != 0) {
    data3[i,8] = data3[i,7]
  }
}
data3 <- data3[,-(6:7)]
write.csv(data3, file = "PPI_match.csv", row.names = F)

