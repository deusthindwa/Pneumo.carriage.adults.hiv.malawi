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

#fit model to individual trajectories & obtain predictions and 95%CI
model_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
model_nvt = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$prevt <- model_vt$fitted.values
crude$prenvt <- model_nvt$fitted.values
crude <- crude %>% mutate(acqvt1 = prevt/11, acqvt2 = prevt/42, acqvt3 = prevt/365.25, acqnvt1 = prenvt/11, acqnvt2 = prenvt/42, acqnvt3 = prenvt/365.25)

#join observed and predicted datasets for survey year
crude1 <- crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(acqvt1 = mean(acqvt1), acqvt2 = mean(acqvt2), acqvt3 = mean(acqvt3), catg = "VT acquisition ")
crude2 <- crude %>% filter(nvtcarr != 0) %>% group_by(year) %>% summarise(acqnvt1 = mean(acqnvt1), acqnvt2 = mean(acqnvt2), acqnvt3 = mean(acqnvt3))
crude3 <- crude %>% filter(vtcarr != 0) %>% group_by(year, agegp) %>% summarise(foi1 = mean(acqvt2))
crude4 <- crude %>% filter(nvtcarr != 0) %>% group_by(year, agegp) %>% summarise(foi2 = mean(acqnvt2))
crude5 <- crude %>% filter(vtcarr != 0) %>% group_by(agegp, year) %>% summarise(foi3 = mean(acqvt2))
crude6 <- crude %>% filter(nvtcarr != 0) %>% group_by(agegp, year) %>% summarise(foi4 = mean(acqnvt2))

A <- ggplot(data = crude1) +
  geom_line(aes(x = year, y = acqvt1, color = "11d"), lty = "dashed", size = 0.7) +
  geom_line(aes(x = year, y = acqvt2, color = "42d"), lty = "dashed", size = 0.7) +
  scale_colour_manual(name = "Carriage duration", values = c("11d" = "darkblue", "42d" = "darkred")) +
  ylim(0, 0.03) +
  labs(title = "Overall VT carriage", x = "Survey year", y = "Daily carriage acquisition") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = c(0.7, 0.7))

B <- ggplot(data = crude2) +
  geom_line(aes(x = year, y = acqnvt1, color = "11d"), lty = "dashed", size = 0.7) +
  geom_line(aes(x = year, y = acqnvt2, color = "42d"), lty = "dashed", size = 0.7) +
  scale_colour_manual(name = "Carriage duration", values = c("11d" = "darkblue", "42d" = "darkred")) +
  ylim(0, 0.03) +
  labs(title = "Overall NVT carriage", x = "Survey year", y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = c(0.7, 0.7))

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "SFig2_sens_carr_dur.tiff"),
       plot = (A | B | C | D),
       width = 13, height = 4, unit="in", dpi = 200)
