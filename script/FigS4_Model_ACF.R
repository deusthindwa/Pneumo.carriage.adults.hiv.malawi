#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1) 

#======================================================================================

#subset for VT & fit model to obtain residuals for ACF checking
overall_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
female_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 1))
male_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 2))
lses_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 1))
hses_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 2))
sart_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 1))
lart_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 2))
nochild5_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 1))
yeschild5_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 2))

#subset for VT & fit model to obtain residuals for ACF checking
overall_nvt = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
female_nvt = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 1))
male_nvt = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 2))
lses_nvt = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 1))
hses_nvt = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 2))
sart_nvt = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 1))
lart_nvt = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 2))
nochild5_nvt = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 1))
yeschild5_nvt = gam(nvtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 2))

#======================================================================================

A <- ggAcf(overall_vt$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(VT) Overall", x = "", y = "ACF") + 
  theme_bw()

B <- ggAcf(female_vt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(VT) Female", x = "", y = "ACF") + 
  theme_bw() 

C <- ggAcf(male_vt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(VT) Male", x = "", y = "ACF") + 
  theme_bw()  

D <- ggAcf(lses_vt$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(VT) Low SES", x = "", y = "ACF") + 
  theme_bw()  

E <- ggAcf(hses_vt$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(VT) Middle/High SES", x = "", y = "ACF") + 
  theme_bw()  

F <- ggAcf(sart_vt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(VT) ART <3y", x = "", y = "ACF") + 
  theme_bw()  

G <- ggAcf(lart_vt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(VT) ART ≥3y", x = "", y = "ACF") + 
  theme_bw()  

H <- ggAcf(nochild5_nvt$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(VT) Without <5y child", x = "", y = "ACF") + 
  theme_bw()  

I <- ggAcf(yeschild5_nvt$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(VT) With <5y child", x = "", y = "ACF") + 
  theme_bw()  

J <- ggAcf(overall_nvt$residuals, type = "correlation", lag.max = 5,  color = "red", ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(NVT) Overall", x = "Lag year", y = "ACF") + 
  theme_bw()  

K <- ggAcf(female_nvt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(NVT) Female", x = "Lag year", y = "ACF") + 
  theme_bw()  

L <- ggAcf(male_nvt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(NVT) Male", x = "Lag year", y = "ACF") + 
  theme_bw()  

M <- ggAcf(lses_nvt$residuals, type = "correlation", lag.max = 5,  color = "red", ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(NVT) Low SES", x = "Lag year", y = "ACF") + 
  theme_bw()  

N <- ggAcf(hses_nvt$residuals, type = "correlation", lag.max = 5,  color = "red", ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(NVT) Middle/High SES", x = "Lag year", y = "ACF") + 
  theme_bw()  

O <- ggAcf(sart_nvt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(NVT) ART <3y", x = "Lag year", y = "ACF") + 
  theme_bw()  

P <- ggAcf(lart_nvt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(NVT) ART ≥3y", x = "Lag year", y = "ACF") + 
  theme_bw()  

Q <- ggAcf(nochild5_nvt$residuals, type = "correlation", lag.max = 5,  color = "red", ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(NVT) Without <5y child", x = "Lag year", y = "ACF") + 
  theme_bw()  

R <- ggAcf(yeschild5_nvt$residuals, type = "correlation", lag.max = 5,  color = "red", ylim = c(-0.10,0.10), size = 1) + 
  labs(title = "(NVT) With <5y child", x = "Lag year", y = "ACF") + 
  theme_bw()  

#=======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "FigS4_Model_ACF.tiff"),
       plot = (A | B | C | D | E | F | G | H | I)/(J | K | L | M | N | O | P | Q | R ),
       width = 23, height = 6, unit="in", dpi = 200)

