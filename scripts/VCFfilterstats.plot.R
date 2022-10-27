## ------------------------------------------------------------------------- ##
# --------------------  plot VCFtools stats files ----------------------- #

# These plots were downloaded from Shannon O'Leary's github repository https://github.com/sjoleary/SNPFILT (updated 24 Jan 2018) and modified to enable rapid plotting of all panels #

library(ggplot2)
library(patchwork)

plot.vcf.ind.stats <- function(ind.file){
  
  # plot missing data per indv ----
  p1 <- ggplot(ind.file, aes(x = MISS)) +
    geom_histogram(binwidth = .01, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(MISS, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = 0.5),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "missing data per indv") 
  
  # plot Fis per indv ----
  p2 <- ggplot(ind.file, aes(x = Fis)) +
    geom_histogram(binwidth = .01, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(Fis, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = 0),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "Fis per indv")
  
  # plot read depth per indv ----
  p3 <- ggplot(ind.file, aes(x = MEAN_DEPTH)) +
    geom_histogram(binwidth = 10, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(MEAN_DEPTH, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = 20),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "mean read depth per indv") 
  
  # plot depth vs missing ----
  p4 <- ggplot(ind.file, aes(x = MEAN_DEPTH, y = MISS)) +
    geom_point() +
    geom_vline(aes(xintercept = mean(MEAN_DEPTH, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = 20),
               color = "darkblue", linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = mean(MISS, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = 0.5),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "mean depth per indv", y = "% missing data") 
  
  # plot Fis vs missing data per indv ----
  p5 <- ggplot(ind.file, aes(x = Fis, y = MISS)) +
    geom_point() +
    geom_vline(aes(xintercept = mean(Fis, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = 0),
               color = "darkblue", linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = mean(MISS, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = 0.5),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "Fis per indv", y = "% missing data") 
  
  # plot Fis vs mean depth per indv ----
  p6 <- ggplot(ind.file, aes(x = Fis, y = MEAN_DEPTH)) +
    geom_point() +
    geom_vline(aes(xintercept = mean(Fis, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = 0),
               color = "darkblue", linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = mean(MEAN_DEPTH, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = 20),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "Fis per indv", y = "mean depth per indv") 
  
plot.ind.stats <- (p1 | p2) / (p3 | p4) / (p5 | p6)  

return(plot.ind.stats)
}


plot.vcf.loc.stats <- function(loc.file){
  

  # plot distribution missing data per locus ----
  p7 <- ggplot(loc.file, aes(x = MISS)) +
    geom_histogram(binwidth = 0.01, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(MISS, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = 0.1),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "% missing data per locus") 
  
  # plot distribution mean read depth ----
  p8 <- ggplot(loc.file, aes(x = MEAN_DEPTH)) +
    geom_histogram(binwidth = 5, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(MEAN_DEPTH, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = 20),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "mean read depth per locus") 
  
  # plot read depth vs missing data ----
  p9 <- ggplot(loc.file, aes(x = MEAN_DEPTH, y = MISS)) +
    geom_point() +
    geom_vline(aes(xintercept = mean(MEAN_DEPTH, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = 20),
               color = "darkblue", linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = mean(MISS, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = 0.1),
               color = "darkblue", linetype = "dashed", size = 1) +
    labs(x = "mean depth per locus", y = "% missing data") 
  
  
  # plot no of SNPs per locus ----
  p10 <- loc.file %>%
    count(CHR) %>%
    ggplot(aes(x = n)) +
    geom_histogram(binwidth = 1, color = "black", fill = "grey95") + 
    labs(x = "number of SNPs per locus") 
  
  
  temp <- loc.file %>%
    count(CHR)
  
  # plot number of SNPs per contig vs. mean depth ----
  p11 <- left_join(temp, loc.file) %>%
    ggplot() +
    geom_point(aes(x = n, y = MEAN_DEPTH)) +
    labs(x = "number of SNPs per contig", y = "mean depth") 
  
plot.loc.stats <- (p7 | p8) / (p9 | p10) / (p11)  
return(plot.loc.stats)
  
}


