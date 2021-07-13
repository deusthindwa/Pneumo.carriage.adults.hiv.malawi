#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)
dur = 42

#======================================================================================

micro2 <- micro1 %>% 
  select(labid, pid, st1_serotype, st2_serotype, st3_serotype) %>%
  rename("st1" = st1_serotype, "st2" = st2_serotype, "st3" = st3_serotype) %>% filter(st1 !="")

micro2$st1[micro2$st1 == "NT"] <- "NVT"
micro2$st2[micro2$st2 == "NT"] <- "NVT"
micro2$st3[micro2$st3 == "NT"] <- "NVT"

micro2 <- micro2 %>% mutate(mc = if_else(st1 == "VT" & st2 == "NVT", 1L,
                                       if_else(st1 == "NVT" & st2 == "VT", 1L,
                                               if_else(st1 == "VT" & st3 == "NVT", 1L,
                                                       if_else(st1 == "NVT" & st3 == "VT", 1L,
                                                               if_else(st2 == "VT" & st3 == "NVT", 1L,
                                                                       if_else(st2 == "NVT" & st3 == "VT", 1L, 0L))))))
)

#merge with the main dataset
micro2 <- micro2 %>% select(labid, pid, mc) %>% filter(mc == 1)
pcvpa.mod <- left_join(pcvpa.mod, micro2)

 #recode the outcome in the main dataset baswed on multiple serotype
pcvpa.mod$nvtcarr2 <- pcvpa.mod$nvtcarr
pcvpa.mod$nvtcarr2[pcvpa.mod$mc == 1] <- 1

pcvpa.mod$vtcarr2 <- pcvpa.mod$vtcarr
pcvpa.mod$vtcarr2[pcvpa.mod$mc == 1] <- 1

pcvpa.mod <- pcvpa.mod %>% select(pid, labid, nvtcarr, vtcarr, nvtcarr1, vtcarr1, nvtcarr2, vtcarr2, mc, everything())

#=======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(nvtcarr, nvtcarr2, vtcarr, vtcarr2, age, year)

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24",
                                   if_else(age >24 & age <=29, "25-29",
                                           if_else(age >29 & age <=34, "30-34",
                                                   if_else(age >34 & age <=40, "35-40", NA_character_))))))

#fit model to obtain predictions of prevalence
model_crude = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$foi_vt <- model_crude$fitted.values

model_crude = gam(vtcarr2 ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$foi_vt2 <- model_crude$fitted.values

model_crude = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$foi_nvt <- model_crude$fitted.values

model_crude = gam(nvtcarr2 ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$foi_nvt2 <- model_crude$fitted.values

#plot prevalence preditions for single and multiple carriage

A <- rbind(
  rbind(crude %>% filter(vtcarr == 1) %>% group_by(year) %>% summarise(foi = mean(foi_vt)) %>% mutate(detection = "Latex", catg = "VT carriage by year") %>% rename("varx" = year),
        crude %>% filter(vtcarr2 == 1) %>% group_by(year) %>% summarise(foi = mean(foi_vt2)) %>% mutate(detection = "Microarray", catg = "VT carriage by year") %>% rename("varx" = year),
        crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% summarise(foi = mean(foi_vt)) %>% mutate(detection = "Latex", catg = "VT carriage by age group") %>% rename("varx" = agegp),
        crude %>% filter(vtcarr2 == 1) %>% group_by(agegp) %>% summarise(foi = mean(foi_vt2)) %>% mutate(detection = "Microarray", catg = "VT carriage by age group") %>% rename("varx" = agegp)),
  
  rbind(crude %>% filter(nvtcarr == 1) %>% group_by(year) %>% summarise(foi = mean(foi_nvt)) %>% mutate(detection = "Latex", catg = "NVT carriage by year") %>% rename("varx" = year),
        crude %>% filter(nvtcarr2 == 1) %>% group_by(year) %>% summarise(foi = mean(foi_nvt2)) %>% mutate(detection = "Microarray", catg = "NVT carriage by year") %>% rename("varx" = year),
        crude %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% summarise(foi = mean(foi_nvt)) %>% mutate(detection = "Latex", catg = "NVT carriage by age group") %>% rename("varx" = agegp),
        crude %>% filter(nvtcarr2 == 1) %>% group_by(agegp) %>% summarise(foi = mean(foi_nvt2)) %>% mutate(detection = "Microarray", catg = "NVT carriage by age group") %>% rename("varx" = agegp))
) %>% mutate(varx = as_factor(varx), detection = as_factor(detection), catg = as_factor(catg)) %>%
  
  ggplot() +
  geom_line(aes(x = varx, y = foi, color = detection, group = detection), lty = "dashed", size = 0.8) +
  facet_grid(.~catg, scales = "free_x") +
  theme_bw() +
  ylim(0, 0.35) +
  labs(title = "", x = "", y = "Carriage prevalence") +
  guides(color=guide_legend(title=""), fill = FALSE) +
  theme(axis.text.x = element_text(size = 7, face = "bold"), axis.text.y = element_text(size = 7, face = "bold")) +
  theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))

#=======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(nvtcarr, nvtcarr2, vtcarr, vtcarr2, age, year)

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24",
                                   if_else(age >24 & age <=29, "25-29",
                                           if_else(age >29 & age <=34, "30-34",
                                                   if_else(age >34 & age <=40, "35-40", NA_character_))))))

#fit model to obtain predictions of acquisitions
model_crude = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$foi_vt <- model_crude$fitted.values/dur

model_crude = gam(vtcarr2 ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$foi_vt2 <- model_crude$fitted.values/dur

model_crude = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$foi_nvt <- model_crude$fitted.values/dur

model_crude = gam(nvtcarr2 ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$foi_nvt2 <- model_crude$fitted.values/dur

#plot acquisition preditions for single and multiple carriage

B <- rbind(
rbind(crude %>% filter(vtcarr == 1) %>% group_by(year) %>% summarise(foi = mean(foi_vt)) %>% mutate(detection = "Latex", catg = "VT carriage by year") %>% rename("varx" = year),
      crude %>% filter(vtcarr2 == 1) %>% group_by(year) %>% summarise(foi = mean(foi_vt2)) %>% mutate(detection = "Microarray", catg = "VT carriage by year") %>% rename("varx" = year),
      crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% summarise(foi = mean(foi_vt)) %>% mutate(detection = "Latex", catg = "VT carriage by age group") %>% rename("varx" = agegp),
      crude %>% filter(vtcarr2 == 1) %>% group_by(agegp) %>% summarise(foi = mean(foi_vt2)) %>% mutate(detection = "Microarray", catg = "VT carriage by age group") %>% rename("varx" = agegp)),

rbind(crude %>% filter(nvtcarr == 1) %>% group_by(year) %>% summarise(foi = mean(foi_nvt)) %>% mutate(detection = "Latex", catg = "NVT carriage by year") %>% rename("varx" = year),
      crude %>% filter(nvtcarr2 == 1) %>% group_by(year) %>% summarise(foi = mean(foi_nvt2)) %>% mutate(detection = "Microarray", catg = "NVT carriage by year") %>% rename("varx" = year),
      crude %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% summarise(foi = mean(foi_nvt)) %>% mutate(detection = "Latex", catg = "NVT carriage by age group") %>% rename("varx" = agegp),
      crude %>% filter(nvtcarr2 == 1) %>% group_by(agegp) %>% summarise(foi = mean(foi_nvt2)) %>% mutate(detection = "Microarray", catg = "NVT carriage by age group") %>% rename("varx" = agegp))
) %>% mutate(varx = as_factor(varx), detection = as_factor(detection), catg = as_factor(catg)) %>%
  
  ggplot() +
  geom_line(aes(x = varx, y = foi, color = detection, group = detection), lty = "dashed", size = 0.6) +
  facet_grid(.~catg, scales = "free_x") +
  theme_bw() +
  ylim(0, 0.008) +
  labs(title = "", x = "", y = "Daily carriage acquisition") +
  guides(color = FALSE, fill = FALSE) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.text.x = element_text(size = 7, face = "bold"), axis.text.y = element_text(size = 7, face = "bold")) +
  theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8))

#=======================================================================================

ggsave(here("output", "FigS4_sens_mult_carr.tiff"),
       plot = (A/B),
       width = 14, height = 6, unit="in", dpi = 200)

