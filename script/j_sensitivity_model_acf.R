#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1) 

#======================================================================================

#subset for NVT(-ST3) & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(nvtcarr, agegp, surv, sex, nochild5) %>% group_by(surv)
model_overall_nvt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

#subset for VT(+ST3) & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(vtcarr, agegp, surv, sex, nochild5) %>% group_by(surv)
model_overall_vt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

#subset for NVT(+ST3) & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(nvtcarr1, agegp, surv, sex, nochild5) %>% group_by(surv)
model_overall_nvt1 = gam(nvtcarr1 ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

#subset for VT(-ST3) & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(vtcarr1, agegp, surv, sex, nochild5) %>% group_by(surv)
model_overall_vt1 = gam(vtcarr1 ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

#subset for female & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(nvtcarr, agegp, surv, sex, nochild5) %>% filter(sex == 0) %>% group_by(surv)
model_female_nvt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

ACFdata = pcvpa.mod %>% select(vtcarr, agegp, surv, sex, nochild5) %>% filter(sex == 0) %>% group_by(surv)
model_female_vt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

#subset for male & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(nvtcarr, agegp, surv, sex, nochild5) %>% filter(sex == 1) %>% group_by(surv)
model_male_nvt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

ACFdata = pcvpa.mod %>% select(vtcarr, agegp, surv, sex, nochild5) %>% filter(sex == 1) %>% group_by(surv)
model_male_vt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

#subset for NVT ART duration & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(nvtcarr, agegp, surv, artdur, sex, nochild5) %>% filter(artdur == 0) %>% group_by(surv)
model_short_nvt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

ACFdata = pcvpa.mod %>% select(nvtcarr, agegp, surv, artdur, sex, nochild5) %>% filter(artdur == 1) %>% group_by(surv)
model_medium_nvt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

#subset for VT ART duration & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(vtcarr, agegp, surv, artdur, sex, nochild5) %>% filter(artdur == 0) %>% group_by(surv)
model_short_vt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

ACFdata = pcvpa.mod %>% select(vtcarr, agegp, surv, artdur, sex, nochild5) %>% filter(artdur == 1) %>% group_by(surv)
model_medium_vt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

#subset for NVT CD4+ & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(nvtcarr, agegp, surv, cd4cnt, sex, nochild5) %>% filter(cd4cnt == 0) %>% group_by(surv)
ACFdata = ACFdata %>% mutate(surv = if_else(surv == 8L | surv == 7L, 6L, surv))
model_cdlow_nvt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

ACFdata = pcvpa.mod %>% select(nvtcarr, agegp, surv, cd4cnt, sex, nochild5) %>% filter(cd4cnt == 1) %>% group_by(surv)
ACFdata = ACFdata %>% mutate(surv = if_else(surv == 8L | surv == 7L, 6L, surv))
model_cdhigh_nvt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

#subset for VT CD4+ & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(vtcarr, agegp, surv, cd4cnt, sex, nochild5) %>% filter(cd4cnt == 0) %>% group_by(surv)
ACFdata = ACFdata %>% mutate(surv = if_else(surv == 8L | surv == 7L, 6L, surv))
model_cdlow_vt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

ACFdata = pcvpa.mod %>% select(vtcarr, agegp, surv, cd4cnt, sex, nochild5) %>% filter(cd4cnt == 1) %>% group_by(surv)
ACFdata = ACFdata %>% mutate(surv = if_else(surv == 8L | surv == 7L, 6L, surv))
model_cdhigh_vt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

#subset for NVT SES & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(nvtcarr, agegp, surv, sescat, sex, nochild5) %>% filter(sescat == 0) %>% group_by(surv)
model_low_nvt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

ACFdata = pcvpa.mod %>% select(nvtcarr, agegp, surv, sescat, sex, nochild5) %>% filter(sescat == 1) %>% group_by(surv)
model_mid_nvt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

ACFdata = pcvpa.mod %>% select(nvtcarr, agegp, surv, sescat, sex, nochild5) %>% filter(sescat == 2) %>% group_by(surv)
model_high_nvt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps"), family = binomial(link = "cloglog"), data = ACFdata)

#subset for VT SES & fit model to obtain residuals for ACF checking
ACFdata = pcvpa.mod %>% select(vtcarr, agegp, surv, sescat, sex, nochild5) %>% filter(sescat == 0) %>% group_by(surv)
model_low_vt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

ACFdata = pcvpa.mod %>% select(vtcarr, agegp, surv, sescat, sex, nochild5) %>% filter(sescat == 1) %>% group_by(surv)
model_mid_vt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

ACFdata = pcvpa.mod %>% select(vtcarr, agegp, surv, sescat, sex, nochild5) %>% filter(sescat == 2) %>% group_by(surv)
model_high_vt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = ACFdata)

#======================================================================================

A <- ggAcf(model_overall_nvt$residuals, type = "correlation", lag.max = 8, color = "red", ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Overall NVT(-ST3)", x = "", y = "ACF") + 
  theme_bw()

B <- ggAcf(model_overall_vt$residuals, type = "correlation", lag.max = 8, color = "red", ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Overall VT(+ST3)", x = "", y = "ACF") + 
  theme_bw()   

C <- ggAcf(model_overall_nvt1$residuals, type = "correlation", lag.max = 8, ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Overall NVT(+ST3)", x = "", y = "ACF") + 
  theme_bw() 

D <- ggAcf(model_overall_vt1$residuals, type = "correlation", lag.max = 8, ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Overall VT(-ST3)", x = "", y = "ACF") + 
  theme_bw()  

E <- ggAcf(model_female_nvt$residuals, type = "correlation", lag.max = 8, color = "red", ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Female NVT", x = "", y = "ACF") + 
  theme_bw()  

F <- ggAcf(model_female_vt$residuals, type = "correlation", lag.max = 8, color = "red", ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Female VT", x = "", y = "ACF") + 
  theme_bw()  

G <- ggAcf(model_male_nvt$residuals, type = "correlation", lag.max = 8, ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Male NVT", x = "", y = "ACF") + 
  theme_bw()  

H <- ggAcf(model_male_vt$residuals, type = "correlation", lag.max = 8, ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Male VT", x = "", y = "ACF") + 
  theme_bw()  

I <- ggAcf(model_short_nvt$residuals, type = "correlation", lag.max = 7, color = "red", ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "ART ≤3y NVT", x = "", y = "ACF") + 
  theme_bw()  

J <- ggAcf(model_medium_nvt$residuals, type = "correlation", lag.max = 7, color = "red", ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "ART >3y NVT", x = "", y = "ACF") + 
  theme_bw()  

K <- ggAcf(model_short_vt$residuals, type = "correlation", lag.max = 7, ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "ART ≤3y VT", x = "", y = "ACF") + 
  theme_bw()  

L <- ggAcf(model_medium_vt$residuals, type = "correlation", lag.max = 7, ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "ART >3y VT", x = "", y = "ACF") + 
  theme_bw()  

M <- ggAcf(model_cdlow_nvt$residuals, type = "correlation", lag.max = 6, color = "red", ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Low CD4+ NVT", x = "Lag survey", y = "ACF") + 
  theme_bw()  

N <- ggAcf(model_cdhigh_nvt$residuals, type = "correlation", lag.max = 6, color = "red", ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "High CD4+ NVT", x = "Lag survey", y = "ACF") + 
  theme_bw()  

O <- ggAcf(model_cdlow_vt$residuals, type = "correlation", lag.max = 6, ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Low CD4+ VT", x = "Lag survey", y = "ACF") + 
  theme_bw()  

P <- ggAcf(model_cdhigh_vt$residuals, type = "correlation", lag.max = 6, ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "High CD4+ VT", x = "Lag survey", y = "ACF") + 
  theme_bw()  

Q <- ggAcf(model_low_nvt$residuals, type = "correlation", lag.max = 7, color = "red", ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Low-SES NVT", x = "", y = "ACF") + 
  theme_bw()  

R <- ggAcf(model_mid_nvt$residuals, type = "correlation", lag.max = 7, color = "red", ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Mid-SES NVT", x = "", y = "ACF") + 
  theme_bw()  

S <- ggAcf(model_high_nvt$residuals, type = "correlation", lag.max = 7, color = "red", ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "High-SES NVT", x = "Lag survey", y = "ACF") + 
  theme_bw()  

T <- ggAcf(model_low_vt$residuals, type = "correlation", lag.max = 7, ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Low-SES VT", x = "Lag survey", y = "ACF") + 
  theme_bw()  

U <- ggAcf(model_mid_vt$residuals, type = "correlation", lag.max = 7, ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "Mid-SES VT", x = "Lag survey", y = "ACF") + 
  theme_bw()  

V <- ggAcf(model_high_vt$residuals, type = "correlation", lag.max = 7, ylim = c(-0.15,0.15), size = 1) + 
  labs(title = "High-SES VT", x = "Lag survey", y = "ACF") + 
  theme_bw()

#=======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "SFig3_sens_model_acf.tiff"),
       plot = (A | B | C | D | E | F | G | H)/(I | J | K | L | M | N | O | P ) / (Q | R | S | T | U | V),
       width = 22, height = 10, unit="in", dpi = 200)

