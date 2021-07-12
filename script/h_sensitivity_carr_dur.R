#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr, nvtcarr, age, year)

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24",
                                   if_else(age >24 & age <=29, "25-29",
                                           if_else(age >29 & age <=34, "30-34",
                                                   if_else(age >34 & age <=40, "35-40", NA_character_))))))

#fit model to individual trajectories & obtain predictions and 95%CI
model_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
model_nvt = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)

crude$prevt <- model_vt$fitted.values
crude$prenvt <- model_nvt$fitted.values
crude <- crude %>% mutate(acqvt1 = prevt/11, acqvt2 = prevt/42, acqnvt1 = prenvt/11, acqnvt2 = prenvt/42)

#join observed and predicted datasets for survey year
rbind(
rbind((crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(acq = mean(acqvt1), dur = "11 days")),
(crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(acq = mean(acqvt2), dur = "42 days"))) %>% mutate(catg = "VT carriage by year") %>% rename("varx" = year),

rbind((crude %>% filter(nvtcarr != 0) %>% group_by(year) %>% summarise(acq = mean(acqnvt1), dur = "11 days")),
(crude %>% filter(nvtcarr != 0) %>% group_by(year) %>% summarise(acq = mean(acqnvt2), dur = "42 days"))) %>% mutate(catg = "NVT carriage by year") %>% rename("varx" = year),

rbind((crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(acq = mean(acqvt1), dur = "11 days")),
(crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(acq = mean(acqvt2), dur = "42 days"))) %>% mutate(catg = "VT carriage by age") %>% rename("varx" = agegp),

rbind((crude %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(acq = mean(acqnvt1), dur = "11 days")),
(crude %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(acq = mean(acqnvt2), dur = "42 days"))) %>% mutate(catg = "NVT carriage by age") %>% rename("varx" = agegp)
) %>% mutate(varx = as_factor(varx), dur = as_factor(dur), catg = as_factor(catg), fac = if_else(varx == "2015" | varx == "2016" | varx == "2017" | varx == "2018" | varx == "2019", "year", "age")) %>%

ggplot(aes(x = varx, y = acq, color = dur, group = dur)) +
  geom_line(lty = "dashed", size = 0.7) +
  facet_grid(.~catg, scales = "free_x") +
  theme_bw() +
  ylim(0, 0.03) +
  labs(title = "", x = "", y = "Daily carriage acquisition") +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "SFig2_sens_carr_dur.tiff"),
       plot = (A | B | C | D),
       width = 13, height = 4, unit="in", dpi = 200)
