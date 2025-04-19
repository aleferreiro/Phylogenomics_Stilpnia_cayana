# Packages needed ---------------------------------------------------------

library(tidyverse)
library(broom)
library(gt)
library(gtable)
library(gtExtras)
library(scales)

# Load bootstrap params ---------------------------------------------------
bs_params <- read.table("06_fastsimcoal2/result/Model_18_MSFSrecal/bs_params_sinPV.txt",
                        header = T)
# Convert migration rate to total migration (ind/gen)
# Uncomment lines for assymetrical migration
bs_params_recal <- bs_params |> 
  mutate(Ne_Guyana = NGUY/2,
         Ne_ParaAmapa = NPAR/2, # Get effective size by dividing by 2 as
         Ne_Marajo = NMJO/2,    # fsc works with gene number and, we are
         Ne_Caatinga = NCAA/2,  # working with a diploid organismns
         Ne_East = NEST/2,
         Ne_West = NWST/2, # 6
         Guyana_ParaAmapa = (GUYPAR*NPAR)/2, # Get total migration (ind/gen) by
         ParaAmapa_Guyana = (GUYPAR*NGUY)/2, # multipying for the effective size
         Guyana_Marajo = (GUYMJO*NMJO)/2,    # of the target population
         Marajo_Guyana = (GUYMJO*NGUY)/2,
         Guyana_Caatinga = (GUYCAA*NCAA)/2,
         Caatinga_Guyana = (GUYCAA*NGUY)/2,
         Guyana_East = (GUYEST*NEST)/2,
         East_Guyana = (GUYEST*NGUY)/2,
         Guyana_West = (GUYWST*NWST)/2,
         West_Guyana = (GUYWST*NGUY)/2, # 5 (11)
         ParaAmapa_Marajo = (PARMJO*NMJO)/2, 
         Marajo_ParaAmapa = (PARMJO*NPAR)/2,
         ParaAmapa_Caatinga = (PARCAA*NCAA)/2,
         Caatinga_ParaAmapa = (PARCAA*NPAR)/2,
         ParaAmapa_East = (PAREST*NEST)/2,
         East_ParaAmapa = (PAREST*NPAR)/2,
         ParaAmapa_West = (PARWST*NWST)/2,
         West_ParaAmapa = (PARWST*NPAR)/2, # 4 (15)
         Marajo_Caatinga = (MJOCAA*NCAA)/2,
         Caatinga_Marajo = (MJOCAA*NMJO)/2,
         Marajo_East = (MJOEST*NEST)/2,
         East_Marajo = (MJOEST*NMJO)/2,
         Marajo_West = (MJOWST*NWST)/2,
         West_Marajo = (MJOWST*NMJO)/2, # 3 (18)
         Caatinga_East = (CAAEST*NEST)/2,
         East_Caatinga = (CAAEST*NCAA)/2,
         Caatinga_West = (CAAWST*NWST)/2,
         West_Caatinga = (CAAWST*NCAA)/2, # 2 (20)
         East_West = (ESTWST*NWST)/2,
         West_East = (ESTWST*NEST)/2, # 1 (21)
         TDiv_1 = TDIV1*1.24, # Get divergence times in years based on the
         TDiv_2 = TDIV2*1.24, # generation time provided by Bird et al. 2020
         TDiv_3 = TDIV3*1.24,
         TDiv_4 = TDIV4*1.24,
         TDiv_5 = TDIV5*1.24,
  ) |> 
  select(Ne_Guyana:TDiv_5) # me quedo con las columnas con los parametros recalculados

# # Estimate confidence intervals for each parameter under normality (NOT USED!!!!)
# bs_ci <- apply(as.matrix(t(bs_params_recal)),
#                1, 
#                function(x){
#                  mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))
#                  }
#                )

# Estimate confidence intervals for each parameter without supposing normality
bs_ci <- apply(as.matrix(t(bs_params_recal)),
               1, 
               function(x){
                 alpha <- 0.05 # Define the confidence level
                 c(quantile(x, alpha / 2), quantile(x, 1 - alpha / 2))
               }
)

# Convert to data frame and name rows
bs_ci <- bs_ci |> 
  as.data.frame()
row.names(bs_ci) <- c("lower_ci","upper_ci")

# load point-estimates to the table previously obtained
params <- read.csv("06_fastsimcoal2/result/Model_18_MSFSrecal/params_recalc.csv") 
row.names(params) <- c("point_est")
# Agrego al data frame de CI para tener todo en un solo data frame
result_fsc <- rbind(params, bs_ci)
# Separate estimator according to the type of parameter
# effective size
result_fsc_Ne <- result_fsc |> 
  select(Ne_Guyana:Ne_West)
# divergence times
result_fsc_Tdiv <- result_fsc |> 
  select(TDiv_1:TDiv_5)
# migration between current populations
result_fsc_Mig <- result_fsc |> 
  select(Guyana_ParaAmapa:West_East) 

# Save results
write_csv(result_fsc_Ne,
          file = "06_fastsimcoal2/result/Model_18_MSFSrecal/bs_result_Ne_sinPV.csv")
write_csv(result_fsc_Tdiv,
          file = "06_fastsimcoal2/result/Model_18_MSFSrecal/bs_result_Tdiv_sinPV.csv")
write_csv(result_fsc_Mig,
          file = "06_fastsimcoal2/result/Model_18_MSFSrecal/bs_result_Mig_sinPV.csv")


# 4. Migration table ---------------------------------------------------------

Mig_df <- result_fsc_Mig |> 
  t() |> # transpose df
  as.data.frame() |>  # reconvert matrix to df
  mutate(across(where(is.numeric), round, 4)) |> 
  unite("CI", 2:3, sep = "-", remove = F) |> # unite CI columns
  mutate(CI = paste0(CI, ")")) |> 
  unite("estimates", 1:2, sep = " (", remove = F) # unite CI and point_est columns

Mig_df_reduced <- Mig_df |> 
  select(1)

Mig_table <- Mig_df_reduced |> 
  mutate(Mig_Source = rownames(Mig_df), .before = estimates) |> # create column with source and target mig
  separate_wider_delim(cols = Mig_Source, delim = "_",
                       names = c("From \\ To", "Mig_target")) |>    # separate in two columns
  pivot_wider(names_from = Mig_target, values_from = estimates) |>   # create a two-wat mig table
  relocate(ParaAmapa, .before = Marajo)

# Make a fancy table
Mig_gtTable <- Mig_table |> 
  gt() |> 
  tab_header(
    title = md("fastsimcoal 2.8 results"),
    subtitle = md("Total migration expressed as `ind/gen` (95% CI)")
  ) |> 
  fmt_number(columns = everything(),
             decimals = 3) |> 
  # tab_source_note(
  #   source_note = md("Model 18: Results obtained generating SFS with 6 
  #                    populations and 5 samples from each one. Using the 
  #                    full VCF with 1,424,327 sites, and RAxML topology. 
  #                    Without setting best point estimates as starting values.")
  # ) |> 
  # data_color(               # Paint cells according to degree of migration
  #   columns = c(2:6),
  #   palette = "YlOrRd"
  #   domain = c(min(Mig_table$estimates), max(Mig_table$estimates))
  # ) |>
  tab_stubhead(
    label = md("")
  ) |> 
  sub_missing(
    missing_text = "" # Remove text from NA cells
  ) |> 
  cols_align(
    align = c("center"),
    columns = everything()
  )

# View it
Mig_gtTable

# Save it as an image
gtsave(Mig_gtTable,
       "06_fastsimcoal2/plot/Model_18_SFSrecal/Mig_table_CI_sinPV.png")
gtsave(Mig_gtTable,
       "06_fastsimcoal2/plot/Model_18_SFSrecal/Mig_table_CI_sinPV.rtf")


# 5. Table with effective sizes and divergence times ----------------------

Parameters_vector <- c("Ne Guyana", "Ne Pará-Amapa","Ne Marajo",
                       "Ne Caatinga", "Ne East", "Ne West",
                       "TDiv 1", "TDiv 2", "TDiv 3", "TDiv 4", "TDiv 5")



NeTdiv_table <- result_fsc |> 
  select(!Guyana_ParaAmapa:West_East) |> # remove migration values
  t() |> # transpose df
  as.data.frame() |> # reconvert matrix to df
  mutate(Parameters = Parameters_vector, .before = point_est)


# Make a fancy table
NeTdiv_gtTable <- NeTdiv_table |> 
  gt(rowname_col = "Parameters") |> 
  tab_header(
    title = md("fastsimcoal 2.8 results"),
    #subtitle = md("Divergence times (`years`) and Effective sizes (`individuals`)")
  ) |> 
  fmt_number(columns = everything()) |> 
  # tab_source_note(
  #   source_note = md("Model 18: Results obtained generating SFS with 6 
  #                    populations and 5 samples from each one. Using the 
  #                    full VCF with 1,424,327 sites, and RAxML topology. 
  #                    Without setting best point estimates as starting values.")
  # ) |> 
  tab_stubhead(
    label = md("Parameters")
  ) |> 
  cols_width(
    everything() ~ px(120)
  ) |> 
  tab_row_group(
    label = md("Effective sizes (`individuals`)"),
    rows = 1:6
  ) |> 
  tab_row_group(
    label = md("Divergence times (`years`)") ,
    rows = 7:11
  ) |> 
  cols_label(
    point_est = "Point estimate",
    lower_ci = "lower 95% CI",
    upper_ci = "upper 95% CI"
  )

# View it
NeTdiv_gtTable

# Save it as an image
gtsave(NeTdiv_gtTable,
       "06_fastsimcoal2/plot/Model_18_SFSrecal/NeTdiv_table_CI_sinPV.png")



