# Packages needed ---------------------------------------------------------

library(tidyverse)
library(broom)
library(gt)
library(gtExtras)


# 1. Parameters rescaling -------------------------------------------------

# Load bootstrap params 
params <- read.table("06_fastsimcoal2/result/Model_18_MSFSrecal/best_params.txt",
                        header = T)
# Convert migration rate to total migration (ind/gen)
# Uncomment lines for assymetrical migration
params_recal <- params |> 
  mutate(Ne_Guyana = NGUY/2,
         Ne_ParaAmapa = NPAR/2, # Get effective size by dividing by 2 as
         Ne_Marajo = NMJO/2, # fsc works with gene number and, we are
         Ne_Caatinga = NCAA/2, # working with a diploid organismns
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

# Save results
write_csv(params_recal,
          file = "06_fastsimcoal2/result/Model_18_MSFSrecal/params_recalc.csv")


# 2. Migration table ---------------------------------------------------------

Mig_df <- params_recal |> 
  select(Guyana_ParaAmapa:West_East) |> # keep migration values
  t() |> # transpose df
  as.data.frame() |> # reconvert matrix to df
  rename(Tot_mig = V1) # rename mig values

Mig_table <- Mig_df |> 
  mutate(Mig_Source = rownames(Mig_df), .before = Tot_mig) |> # create column with source and target mig
  separate_wider_delim(cols = Mig_Source, delim = "_",
                       names = c("From | To", "Mig_target")) |>    # separate in two columns
  pivot_wider(names_from = Mig_target, values_from = Tot_mig) |>    # create a two-wat mig table
  relocate(Guyana, .before = ParaAmapa) 

# Make a fancy table
Mig_gtTable <- Mig_table |> 
  gt() |> 
  tab_header(
    title = md("fastsimcoal 2.7 results"),
    subtitle = md("Total migration `(ind/gen)`")
  ) |> 
  fmt_number(columns = everything(),
             decimals = 4) |> 
  tab_source_note(
    source_note = md("Model 18: Results obtained generating SFS with 6 populations and 5 samples from each one. 
                     Using the full VCF with 1,424,327 sites, and RAxML topology. Full vcf was obtained with a
                     pyton script provided by in the github repository of the ipyrad")
  ) |> 
  data_color(               # Paint cells according to degree of migration
    columns = c(2:7),
    palette = "YlOrRd",
    domain = c(0, max(Mig_df$Tot_mig))
  ) |> 
  tab_stubhead(
    label = md("")
  ) |> 
  sub_missing(
    missing_text = "" # Remove text from NA cells
  )

# View it
Mig_gtTable

# Save it as an image
gtsave(Mig_gtTable,
       "06_fastsimcoal2/plot/Model_18_SFSrecal/Mig_table_withoutCI.png")


# 3. Table with effective sizes  ----------------------

Parameters_vector <- c("Ne Guyana", "Ne Para-Amapa", "Ne Marajo",
                       "Ne Caatinga", "Ne East", "Ne West",
                       "TDiv 1", "TDiv 2", "TDiv 3", "TDiv 4", "TDiv 5")


NeTdiv_table <- params_recal |> 
  select(!Guyana_ParaAmapa:West_East) |> # remove migration values
  t() |> # transpose df
  as.data.frame() |> # reconvert matrix to df
  rename(Point_estimate = V1) |> # rename mig values
  mutate(Parameters = Parameters_vector)


# Make a fancy table
NeTdiv_gtTable <- NeTdiv_table |> 
  gt(rowname_col = "Parameters") |> 
  tab_header(
    title = md("fastsimcoal 2.7 results"),
    #subtitle = md("Divergence times (`years`) and Effective sizes (`individuals`)")
  ) |> 
  fmt_number(columns = everything()) |> 
  tab_source_note(
    source_note = md("Model 18: Results obtained generating SFS with 6 populations and 5 samples from each one. 
                     Using the full VCF with 1,424,327 sites, and RAxML topology. Full vcf was obtained with a
                     pyton script provided by in the github repository of the ipyrad")
  ) |> 
  tab_stubhead(
    label = md("Parameters")
  ) |> 
  cols_width(
    everything() ~ px(150)
  ) |> 
  tab_row_group(
    label = md("Effective sizes (`individuals`)"),
    rows = 1:6
  )  |> 
  tab_row_group(
    label = md("Divergence times (`years`)"),
    rows = 7:11
  ) |>
  cols_label(
    Point_estimate = "Point estimate"
  )

# View it
NeTdiv_gtTable

# Save it as an image
gtsave(NeTdiv_gtTable,
       "06_fastsimcoal2/plot/Model_16_fullVCF/NeTdiv_table_withoutCI.png")
