#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#=======================================================================

# load the require packages
if (!require(pacman)){
  install.packages("pacman")
}
pacman::p_load(char = c("tidyverse", "table1", "readstata13", "patchwork", "boot","mgcv", "devtools", "Metrics", 
                        "MuMIn","PropCIs", "forecast", "asbio", "scam", "here"))

#Yesoptions(stringsAsFactors = FALSE)
setwd(here::here())

#load phirst datasets (household-level, individual-level, follow-up & antibiotic use)
pcvpa <- read.dta13(here("data", "PCVPA.dta"))
micro <- read.csv(here("data", "microarray.txt"), sep="\t", header=T)

#select required variables
pcvpa <- rename(select(pcvpa, pid, labid, collection_date, surv, serotype, risk_h, age_flr, sex, artdat_adj, artreg, ctx, cd4cnt, nochild5, ses_cat),
                date = collection_date, age = age_flr, artdate = artdat_adj, sescat = ses_cat)

micro <- rename(select(micro, SampleID.1, ArraySero), labid = SampleID.1, arraysero = ArraySero)

#=======================================================================

#cleaning/recoding variables for descriptive analysis
pcvpa.des <- pcvpa

#survey number
pcvpa.des$surv <- as.integer(pcvpa.des$surv)

#serogroup
pcvpa.des$serogroup <- if_else(pcvpa.des$serotype == "NoCarriage" | pcvpa.des$serotype == "Pending", "None", 
                               if_else(pcvpa.des$serotype == "NVT", "NVT", "VT"))

pcvpa.des <- select(pcvpa.des, pid, labid, date, surv, serotype, serogroup, age, sex, artdate, artreg, ctx, cd4cnt, nochild5, sescat)

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
                            if_else(pcvpa.des$sescat == 2, "Middle", 
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

#sensitivity of ST3
pcvpa.mod$vtcarr1 <- if_else(pcvpa.mod$serotype != "NVT" & !is.na(pcvpa.mod$serotype) & pcvpa.mod$serotype != 3, 1L, 0L)

pcvpa.mod$nvtcarr1 <- if_else(pcvpa.mod$serotype == "NVT" | pcvpa.mod$serotype == "3", 1L, 0L)
pcvpa.mod$nvtcarr1[is.na(pcvpa.mod$nvtcarr1)] <- 0

#age group
pcvpa.mod$agegp <-  if_else(pcvpa.mod$age <=20, 19L,
                          if_else(pcvpa.mod$age >=21 & pcvpa.mod$age <=23, 22L,
                                if_else(pcvpa.mod$age >=24 & pcvpa.mod$age <=26, 25L,
                                      if_else(pcvpa.mod$age >=27 & pcvpa.mod$age <=29, 28L,
                                            if_else(pcvpa.mod$age >=30 & pcvpa.mod$age <=32, 31L,
                                                  if_else(pcvpa.mod$age >=33 & pcvpa.mod$age <=35, 34L,
                                                        if_else(pcvpa.mod$age >=36 & pcvpa.mod$age <=38, 37L, 40L)))))))
#sex
pcvpa.mod$sex <- if_else(pcvpa.mod$sex == "Female", 0L, 
                         if_else(pcvpa.mod$sex == "Male", 1L, NA_integer_))

#ART duration
pcvpa.mod$artdur <- if_else(pcvpa.mod$artdur <=3, 0L, 
                                    if_else(pcvpa.mod$artdur >3 & pcvpa.mod$artdur <20, 1L, NA_integer_))

#ART regimen
pcvpa.mod$artreg <- if_else(pcvpa.mod$artreg == "First line", 0L,
                            if_else(pcvpa.mod$artreg == "Second line", 1L, NA_integer_))

#Cotrimoxozole
pcvpa.mod$ctx <- if_else(pcvpa.mod$ctx == "No", 0L,
                         if_else(pcvpa.mod$ctx =="Yes", 1L, NA_integer_))

#cd4+ count
pcvpa.mod$cd4cnt <- if_else(pcvpa.mod$cd4cnt < 250, 0L, 
                            if_else(pcvpa.mod$cd4cnt >=250 & pcvpa.mod$cd4cnt < 3000, 1L, NA_integer_))

#adults living with children in the household
pcvpa.mod$nochild5 <- if_else(is.na(pcvpa.mod$nochild5), 1L, pcvpa.mod$nochild5)

pcvpa.mod$nochild5 <- if_else(pcvpa.mod$nochild5 == 0, 0L, 
                              if_else(pcvpa.mod$nochild5 >=1 & pcvpa.mod$nochild5 <5, 1L, NA_integer_))

#social economic status
pcvpa.mod$sescat <- if_else(pcvpa.mod$sescat == "Low", 0L, 
                            if_else(pcvpa.mod$sescat == "Middle", 1L, 
                                    if_else(pcvpa.mod$sescat == "High", 2L, NA_integer_)))

pcvpa.mod <- select(pcvpa.mod, pid, labid, nvtcarr, vtcarr, nvtcarr1, vtcarr1, surv, age, agegp, sex, artdur, artreg, ctx, cd4cnt, nochild5, sescat)

#=======================================================================

#descriptive characteristics of study population (figure 1)
source(here("script/b_descriptive_vars.R"))

#overall carriage prevalence & FOI (figure 2)
source(here("script/c_crude_prev_foi.R"))

#carriage prevalence & FOI by sex (figure 3)
source(here("script/d_sex_prev_foi.R"))

#carriage prevalence & FOI by ART duration (figure 4)
source(here("script/e_artdur_prev_foi.R"))

#carriage prevalence & FOI by social economic status (figure 5)
source(here("script/f_ses_prev_foi.R"))

#carriage prevalence & FOI by CD4+ count (figure 6)
source(here("script/g_cd4cnt_prev_foi.R"))

#sensitivity analysis of the spline types(figure 7)
source(here("script/h_sensitivity_spline_prev.R"))
              
#sensitivity analysis of the model residuals (figure 8)
source(here("script/i_sensitivity_model_acf.R"))
