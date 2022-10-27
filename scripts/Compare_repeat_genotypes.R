# Compare the repeat samples in a plink file. 
# @KaraDicks 2020.01

# This script compares the genotypes of repeated samples contained within a plink data.

###########
#~ REQUIREMENTS:
# 
# .ped & .map file
#
#

############
#~ PERFORMS
#
# Pairwise comparison of SNP genotypes between samples with same WG_ID
# e.g. can compare NUB003a_40 with NUB003_47
# letters a & b indicating within plate repeat and _00 indicating library are dealt with. 

##########
#~ OUTPUT
# It outputs a table in the format: 
#~
# differences            n sample_pair          
# <chr>              <int> <chr>                
#~
#
#
# Where there are 4 categories of differences:
  #~ same : the genotypes are identical (Note that in the example datasets, no SNPs fail for both samples)
  #~ one_failed : one sample does not have a genotype and therfore cannot be compared at that SNP
  #~ allelic_dropout : one sample is heterozygous and the other is homozygous. It could of course be that there is a false allele, but allelic dropout is much more likely. 
  #~ alternative_allele : both samples are homozygous for alternative genotypes. A high proportion of alternative alleles would strongly suggest the two samples are from different individuals. 


# This script was created and tested using the Nubian ibex data (3.WG_Projects/WG1712_Nubian_Ibex_Mataab/WGXX_NUB_HybridCapture/1.SNP_ID)

library(tidyverse)
# library(data.table)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# convert IID into WG_ID

convert_ID <- function(input_file, ID_column){
  df <- input_file %>% 
    mutate(WG_ID = gsub("_[0-9][0-9]$", "", {{ID_column}} )) %>% 
    mutate(WG_ID = gsub("[a-z]$", "", WG_ID))
  return(df)
}


export_numbers_diffs <- function(infile){
#########################
#~ Create a genotyping file to work with 

# import the ped file
ped <- read.table(paste0(infile, ".ped"), colClasses = 'character') %>%
  select(-3, -4,-5, -6) 

# Create a list of the loci to act as column headers
n_loci <- 
  read.table(paste0(infile, ".map"), colClasses = 'character') %>% 
  nrow

loci <- 
  read.table(paste0(infile, ".map"), colClasses = 'character') %>% 
  mutate(loci = paste(V1, V2, sep=":")) %>% 
  select(loci) %>% 
  slice(rep(1:n(), each = 2)) %>% 
  mutate(rep = rep(c("a", "b"), times=as.numeric(n_loci))) %>% 
  mutate(names_cols = paste(rep, loci, sep="-")) 

names_cols <- data.frame(names_cols = loci$names_cols)
temp <- data.frame(names_cols = c("FID", "IID")) 
names_cols <- rbind(temp, names_cols)

# Add column headers to ped file
names(ped) <- names_cols$names_cols



# Convert to long format
ped_a <- 
  ped %>% 
  select(FID, IID, starts_with("a-")) %>% 
  pivot_longer(!FID:IID, 
               names_to = "locus_names",
               names_prefix = "a-",
               values_to = "allele_a")
ped_b <- 
  ped %>% 
  select(FID, IID, starts_with("b-")) %>% 
  pivot_longer(!FID:IID, 
               names_to = "locus_names",
               names_prefix = "b-",
               values_to = "allele_b")

ped2 <- full_join(ped_a, ped_b)

##################
# Start looking at the repeated individuals


# Create a list of the repeated individuals
repeat_inds <- 
  ped2 %>% 
  distinct(IID) %>% 
  convert_ID(IID) 

reps <- count(repeat_inds, WG_ID) %>% 
  subset(n > 1)

repeat_inds <- 
  subset(repeat_inds, WG_ID %in% reps$WG_ID)


#~ Subset the ped file to only those which are repeated
ped_reps <- 
  ped2 %>% 
  subset(IID %in% repeat_inds$IID) %>% 
  mutate(genotype = paste0(allele_a, allele_b)) %>% 
  convert_ID(IID)
  
#~ Create a column of genotypes per individual
ped_wide <- 
ped_reps %>% 
  select(locus_names,IID, genotype) %>% 
  pivot_wider(names_from = IID, 
              values_from = genotype)

#Create an empty dataframe for the results
comp_diffs <- data.frame(differences = as.character(),
                         n = as.numeric(),
                         sample_pair = as.character())

#################### LOOP!
for(i in unique(ped_reps$WG_ID)){
  # First find all the possible combinations of samples for individual i 
  samples <- subset(repeat_inds, WG_ID == i)
  combinations <- apply(combn(as.vector(samples$IID),2),2,paste,collapse=':')
  combinations <- data.frame(comb = combinations) %>% 
    separate(comb, 
             into=c("IID1", "IID2"),
             sep=":", remove=F)
  
  
    ##################### LOOP!
  # now run comparisons between each pair of samples
  for(i in as.vector(combinations$comb)){
  samples_i <- subset(combinations, comb == i)
  
  ped_wide <- subset(ped_reps, IID == samples_i$IID1 | IID == samples_i$IID2) %>% 
    select(locus_names,IID, genotype) %>% 
    pivot_wider(names_from = IID, 
                values_from = genotype)
  
  #~ rename the columns
  names(ped_wide)[2] <- "sample_a"
  names(ped_wide)[3] <- "sample_b"
  
  #~ compare the genotypes
  comp_geno <- 
  ped_wide %>% 
    # define if genotypes are identical, different, or a sample has failed
    mutate(comparison = 
             case_when(sample_a == sample_b ~ "same",
                       sample_a == "00" & sample_b == "00" ~ "both_failed",
                       sample_a == "00" & sample_b != "00" ~ "one_failed",
                       sample_a != "00" & sample_b == "00" ~ "one_failed",
                       TRUE ~ "other")) %>% 
    mutate(a_allele1 = substr(sample_a, 1, 1),
           a_allele2 = substr(sample_a, 2, 2),
           b_allele1 = substr(sample_b, 1, 1),
           b_allele2 = substr(sample_b, 2, 2)) %>% 
    mutate(hetz_a = ifelse(comparison != "other", NA,
                                    ifelse(a_allele1 == a_allele2, "homz",
                                           "hetz")),
           hetz_b = ifelse(comparison != "other", NA,
                                     ifelse(b_allele1 == b_allele2, "homz",
                                            "hetz"))
           ) %>%
   mutate(differences = ifelse(comparison != "other", comparison,
                                    ifelse(hetz_a == "homz" & hetz_b == "hetz", 
                                           "allelic_dropout", 
                                           ifelse(hetz_a == "hetz" & hetz_b == "homz",
                                                  "allelic_dropout", 
                                            ifelse(hetz_a == "homz" & hetz_b == "homz" & sample_a != sample_b, "alternative_allele", "other"))))) 
  
  
  #~ count the differences
  temp_out <- comp_geno %>% count(differences) %>% 
    mutate(sample_pair = i)
    
  comp_diffs <- rbind(comp_diffs, temp_out)
  }
}

return(comp_diffs)
}


