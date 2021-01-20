#Written by Deus Thindwa
#Pneumococcal carriage prevalence in HIV-infected adults in PCV era
#Generalized additive model
#14/12/2020 - 30/06/2021

#=======================================================================

#load required packages
pcvpa.packages <- c("tidyverse", "table1", "readstata13", "patchwork", "boot","mgcv", "devtools", "Metrics", "gratia", "MuMIn","cgam", "here")
lapply(pcvpa.packages, library, character.only = TRUE)

#load phirst datasets (household-level, individual-level, follow-up & antibiotic use)
pcvpa <- read.dta13(here("data", "PCVPA.dta"))

#select required variables
pcvpa <- rename(select(pcvpa, pid, collection_date, surv, serotype, risk_h, age_flr, sex, artdat_adj, artreg, ctx, cd4cnt, nochild5, ses_cat),
                date = collection_date, age = age_flr, artdate = artdat_adj, sescat = ses_cat)

#=======================================================================

#cleaning/recoding variables for descriptive analysis
pcvpa.des <- pcvpa

#survey number
pcvpa.des$surv <- as.integer(pcvpa.des$surv)

#serogroup
pcvpa.des$serogroup <- if_else(pcvpa.des$serotype == "NoCarriage" | pcvpa.des$serotype == "Pending", "None", 
                               if_else(pcvpa.des$serotype == "NVT", "NVT", "VT"))

pcvpa.des <- select(pcvpa.des, pid, date, surv, serotype, serogroup, age, sex, artdate, artreg, ctx, cd4cnt, nochild5, sescat)

#serotype
pcvpa.des$serotype <- if_else(pcvpa.des$serotype == "NoCarriage" | pcvpa.des$serotype == "Pending", NA_character_, pcvpa.des$serotype)

#age
pcvpa.des$age <- as.integer(pcvpa.des$age)

#sex
pcvpa.des$sex <- if_else(pcvpa.des$sex == 0, "Female", 
                         if_else(pcvpa.des$sex == 1, "Male", 
                                 if_else(pcvpa.des$sex == 2, "Female", NA_character_)))

#ART duration
pcvpa.des <- rename(pcvpa.des %>% 
                      mutate(artdate = as.integer((date-artdate)/365.25)), artdur = artdate)

#ART regimen
pcvpa.des$artreg <- if_else(pcvpa.des$artreg == "TDF/3TC/EFV" | pcvpa.des$artreg == "TDF/3TC/NVP" | pcvpa.des$artreg == "AZT/3TC/EFV", "First line",
                            if_else(pcvpa.des$artreg == "TDF/3TC/LPVr", "Second line", NA_character_))

#Cotrimoxozole
pcvpa.des$ctx <- if_else(pcvpa.des$ctx == 0, "No", 
                         if_else(pcvpa.des$ctx == 1,"Yes", NA_character_))

#CD4+ count
pcvpa.des$cd4cnt <- as.integer(if_else(pcvpa.des$cd4cnt > 1500, NA_integer_, pcvpa.des$cd4cnt))

#living with <5 years-old children
pcvpa.des$nochild5 <- as.integer(pcvpa.des$nochild5)

#social economic status
pcvpa.des$sescat <- if_else(pcvpa.des$sescat == 1, "Low", 
                            if_else(pcvpa.des$sescat == 2, "Medium", 
                                    if_else(pcvpa.des$sescat == 3, "High", NA_character_)))


#=======================================================================

#cleaning/recoding variables for modelling
pcvpa.mod <- pcvpa.des

#vaccine types
pcvpa.mod$vtcarr <- if_else(pcvpa.mod$serogroup == "VT", 1L, 
                            if_else(pcvpa.mod$serogroup == "None" | pcvpa.mod$serogroup == "NVT", 0L, NA_integer_ ))

#non-vaaccine types
pcvpa.mod$nvtcarr <- if_else(pcvpa.mod$serogroup == "NVT", 1L, 
                             if_else(pcvpa.mod$serogroup == "None" | pcvpa.mod$serogroup == "VT", 0L, NA_integer_))

#sex
pcvpa.mod$sex <- if_else(pcvpa.mod$sex == "Female", 0L, 
                         if_else(pcvpa.mod$sex == "Male", 1L, NA_integer_))

#ART duration
pcvpa.mod$artdur <- if_else(pcvpa.mod$artdur <2, 0L, 
                            if_else(pcvpa.mod$artdur >=2 & pcvpa.mod$artdur <7, 1L, 
                                    if_else(pcvpa.mod$artdur >=7 & pcvpa.mod$artdur <20, 2L, NA_integer_)))

#ART regimen
pcvpa.mod$artreg <- if_else(pcvpa.mod$artreg == "First line", 0L,
                            if_else(pcvpa.mod$artreg == "Second line", 1L, NA_integer_))

#Cotrimoxozole
pcvpa.mod$ctx <- if_else(pcvpa.mod$ctx == "No", 0L,
                         if_else(pcvpa.mod$ctx =="Yes", 1L, NA_integer_))

#cd4+ count
pcvpa.mod$cd4cnt <- if_else(pcvpa.mod$cd4cnt < 350, 0L, 
                            if_else(pcvpa.mod$cd4cnt >=350 & pcvpa.mod$cd4cnt < 3000, 1L, NA_integer_))

#adults living with children in the household
pcvpa.mod$nochild5 <- if_else(pcvpa.mod$nochild5 == 0, 0L, 
                              if_else(pcvpa.mod$nochild5 >1 & pcvpa.mod$nochild5 <5, 1L, NA_integer_))

#social economic status
pcvpa.mod$sescat <- if_else(pcvpa.mod$sescat == "Low", 0L, 
                            if_else(pcvpa.mod$sescat == "Medium", 1L, 
                                    if_else(pcvpa.mod$sescat == "High", 2L, NA_integer_)))

pcvpa.mod <- select(pcvpa.mod, pid, nvtcarr, vtcarr, surv, age, sex, artdur, artreg, ctx, cd4cnt, nochild5, sescat)

#=======================================================================

#adult characteristics for continuous variables (figure 1)
source(here("script/2_continuous_var.R"))

#adult characteristics for categorical variables (table 1)
source(here("script/3_categorical_var.R"))

#participant characteristics
source(here("script/4_model_fitting.R"))

#participant characteristics (supplementary figures 1-8)
source(here("script/5_age_surv_prev.R"))

#participant characteristics (figure 2)
source(here("script/6_age_prev_foi.R"))


