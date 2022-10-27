# KaraDicks
# This script uses the function `basic.stats` in hierfstat to generate summary stats for each population from the SNP files. It takes a hierfstat input file and calculates the following basic populations statistics: Hobs, Hs and Fis. It then summarizes these per population, providing means and confidence intervals. 


library(hierfstat)



basic_stats_summary <- function(input){ 
  
  basic_stats <- basic.stats(input)
  
  pop_Hobs <- data.frame(basic_stats$Ho) %>% 
    rownames_to_column(var = "locus") %>% 
    pivot_longer(!locus, 
                 names_to = "Pop", 
                 values_to = "Hobs")
  
  pop_Hs <- data.frame(basic_stats$Hs) %>% 
    rownames_to_column(var = "locus") %>% 
    pivot_longer(!locus, 
                 names_to = "Pop", 
                 values_to = "Hs")
  
  pop_Fis <- data.frame(basic_stats$Fis) %>% 
    rownames_to_column(var = "locus") %>% 
    pivot_longer(!locus, 
                 names_to = "Pop", 
                 values_to = "Fis")
  
  stats <- full_join(pop_Hobs, pop_Hs)
  stats <- full_join(stats, pop_Fis)
  
  
  Hobs_sum <- 
    stats %>% 
    group_by(Pop) %>% 
    summarise(mean.Hobs = mean(Hobs, na.rm = TRUE),
              sd = sd(Hobs, na.rm = TRUE),
              count.Hobs = n()) %>% 
    mutate(se = sd / sqrt(count.Hobs),
           l_ci.Hobs = (mean.Hobs + -1.96 * (sd/sqrt(count.Hobs))),
           u_ci.Hobs = (mean.Hobs + 1.96 * (sd/sqrt(count.Hobs)))) %>% 
    select(-sd, -se)
  
  
  
  Hs_sum <- 
    stats %>% 
    group_by(Pop) %>% 
    subset(!is.na(Hs)) %>% 
    summarise(mean.Hs = mean(Hs, na.rm = TRUE),
              sd = sd(Hs, na.rm = TRUE),
              count.Hs = n()) %>% 
    mutate(se = sd / sqrt(count.Hs),
           l_ci.Hs = (mean.Hs + -1.96 * (sd/sqrt(count.Hs))),
           u_ci.Hs = (mean.Hs + 1.96 * (sd/sqrt(count.Hs)))) %>% 
    select(-sd, -se)
  
  
  Fis_sum <- 
    stats %>% 
    group_by(Pop) %>% 
    subset(!is.na(Fis)) %>% 
    summarise(mean.Fis = mean(Fis, na.rm = TRUE),
              sd = sd(Fis, na.rm = TRUE),
              count.Fis = n()) %>% 
    mutate(se = sd / sqrt(count.Fis),
           l_ci.Fis = (mean.Fis + -1.96 * (sd/sqrt(count.Fis))),
           u_ci.Fis = (mean.Fis + 1.96 * (sd/sqrt(count.Fis)))) %>% 
    select(-sd, -se)
  
  
  stats_per_pop <- full_join(Hobs_sum, Hs_sum)
  stats_per_pop <- full_join(stats_per_pop, Fis_sum)
  
  return(stats_per_pop)
}