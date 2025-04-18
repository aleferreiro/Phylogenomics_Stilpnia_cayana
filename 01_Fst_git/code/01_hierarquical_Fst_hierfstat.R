# 0. Packages needed -----------------------------------------------------
# install.packages(c("vcfR", "hierfstat", "adegenet", "poppr"))
library(vcfR)
library(tidyverse)
library(hierfstat)
library(adegenet)
library(poppr)
library(ggplot2)


# 1. Read vcf file --------------------------------------------------------
  vcf <- read.vcfR("01_Fst/data/Scay_83_maf_mr80_sinMarajo.recode.vcf")

# 2. Extract genotypic data -----------------------------------------------
genind_obj <- vcfR2genind(vcf)

# 3. Add population information -------------------------------------------
Localities_df <- read.csv("00_Sampling_sites/Sampling_sites83_sinMarajo.csv", sep = ";") #|> 
  # filter(Sample != "UFG4028",
  #        Sample != "UFG5668",
  #        Sample != "UFG5753")  # !BAJA CALIDAD DE READS!!
 
pop(genind_obj) <- as.factor(Localities_df$Level_1)

#strata(genind_obj) <- data.frame(Population = as.factor(Localities_df$popNS),
#                                 Subpopulation = as.factor(Localities_df$Populations2))
genind_obj

# 4. Compute observed and expected Heterozygosity -------------------------

basic_stats_Scay <- basic.stats(genind_obj)
basic_stats_Scay

# 5. Compute Weir & Cockerman Fst -----------------------------------------
Scay_Fst_wc <- wc(genind_obj)
Scay_Fst_wc

# Estimate CI
Scay_Fst_wc_ci <- boot.ppfst(genind_obj, nboot=100, quant = c(0.025,0.975))
Scay_Fst_wc_ci

write.csv(Scay_Fst_wc_ci,
          file = "01_Fst/result/hierfstat_ppfst_ci_ll.csv")

# 6. Compute hierarchical Fst ---------------------------------------------

# Convert genind to hierfstat format
hierfstat_obj <- genind2hierfstat(genind_obj)

# Hierarchical FST calculation
fst_results <- varcomp.glob(levels = data.frame(Population = as.factor(Localities_df$Level_1)),
                            loci = hierfstat_obj[, -1])
fst_results

write.csv(fst_results$F,
          file = "01_Fst/result/hierfstat_amova.csv")

# Estimate CI of amova
fst_results_ci <- boot.vc(levels = data.frame(Population = as.factor(Localities_df$Level_1)),
                          loci = hierfstat_obj[, -1],
                          nboot = 1000)
fst_results_ci
write.csv(fst_results_ci$ci,
          file = "01_Fst/result/hierfstat_amova_CI.csv")



# Significance of the level in genetic differentiation 
fst_test_Lev1 <- test.g(hierfstat_obj[, -1], level = as.factor(Localities_df$Level_1)) 
fst_test_Lev1




# 5. Interpret the results ------------------------------------------------

# The output from varcomp.glob will provide you with the variance components
# and FST values at each level of the hierarchy.

print(fst_results)

# 6. Plotting the results


