#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#======================================================================================

micro <- micro %>% 
  select(labid, pid, st1_serotype, st2_serotype, st3_serotype) %>%
  rename("st1" = st1_serotype, "st2" = st2_serotype, "st3" = st3_serotype) %>% filter(st1 !="")

micro <- micro %>% mutate(st1 = if_else(st1 == "VT", "VT",
                                        if_else(st1 == "NT", "NVT",
                                                if_else(st1 == "NVT", "NVT", NA_character_))),
                          st2 = if_else(st2 == "VT", "VT",
                                        if_else(st2 == "NT", "NVT",
                                                if_else(st2 == "NVT", "NVT", NA_character_))),
                          st3 = if_else(st3 == "VT", "VT",
                                        if_else(st3 == "NT", "NVT",
                                                if_else(st3 == "NVT", "NVT", NA_character_)))
)

micro <- micro %>% mutate(mc = if_else(st1 == "VT" & st2 == "NVT", 1L,
                                       if_else(st1 == "NVT" & st2 == "VT", 1L,
                                               if_else(st1 == "VT" & st3 == "NVT", 1L,
                                                       if_else(st1 == "NVT" & st3 == "VT", 1L,
                                                               if_else(st2 == "VT" & st3 == "NVT", 1L,
                                                                       if_else(st2 == "NVT" & st3 == "VT", 1L, 0L))))))
)

#merge with the main dataset
micro <- micro %>% mutate(mcc = if_else(is.na(mc),0L, 1L)) %>% select(labid, pid, mcc)
pcvpa.mod <- left_join(pcvpa.mod, micro)

#recode the outcome in the main dataset baswed on multiple serotype
pcvpa.mod$nvtcarr2 <- pcvpa.mod$nvtcarr
pcvpa.mod$nvtcarr2[pcvpa.mod$nvtcarr == 0 & pcvpa.mod$mcc == 1] <- 1

pcvpa.mod$vtcarr2 <- pcvpa.mod$vtcarr
pcvpa.mod$vtcarr2[pcvpa.mod$vtcarr == 0 & pcvpa.mod$mcc == 1] <- 1

#=======================================================================================


#subset for a correct dataset
crude = pcvpa.mod %>% select(nvtcarr, nvtcarr2, vtcarr, vtcarr2, agegp, surv, sex, nochild5) %>% group_by(surv)

#fit model to obtain predictions of FOI
model_crude = scam(nvtcarr ~ s(agegp, bs="mdcv") + s(surv, bs="mdcv") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi_nvt <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

model_crude = scam(nvtcarr2 ~ s(agegp, bs="mdcv") + s(surv, bs="mdcv") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi_nvt2 <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

model_crude = scam(vtcarr ~ s(agegp, bs="mdcx") + s(surv, bs="mdcx") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi_vt <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

model_crude = scam(vtcarr2 ~ s(agegp, bs="mdcx") + s(surv, bs="mdcx") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi_vt2 <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp NVT
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi_nvt)) %>% ungroup()) %>% select(agegp, foi) %>% mutate(Category = "Original")

crude2 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr2 == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr2 != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi_nvt2)) %>% ungroup()) %>% select(agegp, foi) %>% mutate(Category = "Multiple")

#join observed and predicted datasets for survey number NVT
crude3 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi_nvt)) %>% ungroup()) %>% select(surv, foi) %>% mutate(Category = "Original")

crude4 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr2 == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr2 != 0) %>% group_by(surv) %>% summarise(foi = mean(foi_nvt2)) %>% ungroup()) %>% select(surv, foi) %>% mutate(Category = "Original")







