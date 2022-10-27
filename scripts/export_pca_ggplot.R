# 2020.04.16

# Function to export pca data for ggplot


library(tidyverse)
library(ggplot2)
library(cowplot)


# Functions to create eigenvalues plots
export_eig <- function(name_pca){
  eig <- data.frame(name_pca$eig)
  names(eig)[1] <- "pca_eig"
  eig$percentage = (eig[, 1]/sum(eig$pca_eig))*100
  eig$Pos <- as.numeric(rownames(eig))
  eig$percentage <- round(eig$percentage, digits = 1)
  return(eig)
}

plot_eig <- function(name_eig_df){
  ggplot(name_eig_df)+
    geom_bar(aes(Pos, percentage), stat="identity",fill = heat.colors(nrow(name_eig_df)))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    xlab("Eigenvalues")+
    ylab("Percent of variance\nexplained")+
    ggtitle("PCA eigenvalues")
}


export_pca <- function(name_pca, name_genL){
  pc1 <- name_pca$scores[,1]
  pc2 <- name_pca$scores[,2]
  pc3 <- name_pca$scores[,3]
  pc4 <- name_pca$scores[,4]
  pc5 <- name_pca$scores[,5]
  pc6 <- name_pca$scores[,6]
  pc7 <- name_pca$scores[,7]
  pc8 <- name_pca$scores[,8]
  pc9 <- name_pca$scores[,9]
  pc10 <- name_pca$scores[,10]
  
  ind_names <- name_genL@ind.names
  
  ggplot_pca <- as.data.frame(cbind(pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10)) %>%
    mutate(ind_names = ind_names) %>%
    mutate(pop = name_genL@pop)
  
  return(ggplot_pca)
}



  
plot_pca_1_2 <- function(name_exported_pca, pop_column_name){
  name_exported_pca %>% 
    rename(pop_column = {{pop_column_name}}) %>% 
    ggplot()+
    geom_hline(yintercept = 0, colour="grey70")+
    geom_vline(xintercept = 0, colour="grey70")+
    geom_point(aes(pc1, pc2, colour = factor(pop_column)), 
               size = 5) +
    ylab(paste0("PC2 (",eig$percentage[2],"%)")) +
    xlab(paste0("PC1 (",eig$percentage[1],"%)")) +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_text(face = "plain"),
      axis.title.y = element_text(face = "plain")) 
  }


plot_pca_1_3 <- function(name_exported_pca, pop_column_name){
  name_exported_pca %>% 
    rename(pop_column = {{pop_column_name}}) %>% 
    ggplot()+
    geom_hline(yintercept = 0, colour="grey70")+
    geom_vline(xintercept = 0, colour="grey70")+
    geom_point(aes(pc1, pc3, colour = factor(pop_column)), 
               size = 5) +
    ylab(paste0("PC3 (",eig$percentage[3],"%)")) +
    xlab(paste0("PC1 (",eig$percentage[1],"%)")) +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_text(face = "plain"),
      axis.title.y = element_text(face = "plain")) 
}
