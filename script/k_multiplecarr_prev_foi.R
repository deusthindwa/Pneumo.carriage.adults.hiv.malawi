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

model_crude = scam(vtcarr ~ s(agegp, bs="mdcv") + s(surv, bs="mdcx") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi_vt <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

model_crude = scam(vtcarr2 ~ s(agegp, bs="mdcv") + s(surv, bs="mdcx") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi_vt2 <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi_nvt)) %>% ungroup()) %>% select(agegp, foi) %>% mutate(Detection = "Latex")

crude2 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr2 == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr2 != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi_nvt2)) %>% ungroup()) %>% select(agegp, foi) %>% mutate(Detection = "Microarray")

crude3 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi_vt)) %>% ungroup()) %>% select(agegp, foi) %>% mutate(Detection = "Latex")

crude4 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr2 == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr2 != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi_vt2)) %>% ungroup()) %>% select(agegp, foi) %>% mutate(Detection = "Microarray")

#join observed and predicted datasets for survey number
crude5 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi_nvt)) %>% ungroup()) %>% select(surv, foi) %>% mutate(Detection = "Latex")

crude6 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr2 == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr2 != 0) %>% group_by(surv) %>% summarise(foi = mean(foi_nvt2)) %>% ungroup()) %>% select(surv, foi) %>% mutate(Detection = "Microarray")

crude7 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi_vt)) %>% ungroup()) %>% select(surv, foi) %>% mutate(Detection = "Latex")

crude8 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr2 == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr2 != 0) %>% group_by(surv) %>% summarise(foi = mean(foi_vt2)) %>% ungroup()) %>% select(surv, foi) %>% mutate(Detection = "Microarray")

#plot FOI curves

A <- ggplot(data = rbind(crude1, crude2)) +
  geom_line(aes(x = agegp, y = foi, color = Detection), lty = "dashed", size = 1) +
  labs(title = "NVT(-ST3), Overall", x = "Age,y", y = "Force of infection") +
  theme_bw() +
  ylim(0, 0.02) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

B <- ggplot(data = rbind(crude5, crude6)) +
  geom_line(aes(x = surv, y = foi, color = Detection), lty = "dashed", size = 1) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  labs(title = "", x = "Survey number", y = "") +
  theme_bw() + 
  ylim(0, 0.02) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

C <- ggplot(data = rbind(crude3, crude4)) +
  geom_line(aes(x = agegp, y = foi, color = Detection), lty = "dashed", size = 1) +
  labs(title = "VT(+ST3), Overall", x = "Age,y", y = "Force of infection") +
  theme_bw() +
  ylim(0, 0.02) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

D <- ggplot(data = rbind(crude7, crude8)) +
  geom_line(aes(x = surv, y = foi, color = Detection), lty = "dashed", size = 1) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  labs(title = "", x = "Survey number", y = "") +
  theme_bw() + 
  ylim(0, 0.02) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = c(0.6, 0.7))

#=======================================================================================

ggsave(here("output", "SFig4_sens_mult_carr.tiff"),
       plot = (A | B | C | D),
       width = 14, height = 4, unit="in", dpi = 200)

