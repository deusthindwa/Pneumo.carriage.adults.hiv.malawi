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
overall_vt1 = gam(vtcarr1 ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
female_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 1))
male_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 2))
lses_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 1))
hses_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 2))
sart_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 1))
lart_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 2))
nochild5_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 1))
yeschild5_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 2))

#======================================================================================

A <- ggAcf(overall_vt$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "VT carriage", x = "", y = "") + 
  theme_bw()

B <- ggAcf(female_vt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "Female", x = "", y = "") + 
  theme_bw() 

C <- ggAcf(male_vt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "Male", x = "", y = "") + 
  theme_bw()  

D <- ggAcf(lses_vt$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "Low SES", x = "", y = "") + 
  theme_bw()  

E <- ggAcf(hses_vt$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "Middle/High SES", x = "Survey year", y = "Autocorrelation (ACF)") + 
  theme_bw()  

F <- ggAcf(sart_vt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "ART <3y", x = "Survey year", y = "") + 
  theme_bw()  

G <- ggAcf(lart_vt$residuals, type = "correlation", lag.max = 5, ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "ART ≥3y", x = "Survey year", y = "") + 
  theme_bw()  

H <- ggAcf(yeschild5_vt$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "With <5y child", x = "Survey year", y = "") + 
  theme_bw()

I <- ggAcf(nochild5_vt$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "Without <5y child", x = "Survey year", y = "") + 
  theme_bw()  

#=======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "FigS4_VT_carriage_ACF.png"),
       plot = (A | B | C | D | E)/(F | G | H | I),
       width = 11, height = 6, unit="in", dpi = 300)

