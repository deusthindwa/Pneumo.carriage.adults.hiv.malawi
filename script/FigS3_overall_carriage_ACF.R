#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1) 

#=====================================================================================

pcvpa.mod2  <- pcvpa.mod %>% mutate(ctime = seq(1:2067)) 

#subset for VT & fit model to obtain residuals for ACF checking
overall = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod2)
female = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 1))
male = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 2))
lses = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 1))
hses = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 2))
sart = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 1))
lart = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 2))
nochild5 = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 1))
yeschild5 = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 2))

#======================================================================================

A <- ggAcf(overall$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "Overall carriage", x = "", y = "Autocorrelation (ACF)") + 
  theme_bw()

B <- ggAcf(female$residuals, type = "correlation", lag.max = 5, ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "Female", x = "", y = "") + 
  theme_bw()

C <- ggAcf(male$residuals, type = "correlation", lag.max = 5, ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "Male", x = "", y = "") + 
  theme_bw()

D <- ggAcf(lses$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "Low SES", x = "", y = "") + 
  theme_bw()

E <- ggAcf(hses$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "Middle/High SES", x = "Survey year", y = "") + 
  theme_bw()

F <- ggAcf(sart$residuals, type = "correlation", lag.max = 5, ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "ART <3y", x = "Survey year", y = "Autocorrelation (ACF)") + 
  theme_bw()

G <- ggAcf(lart$residuals, type = "correlation", lag.max = 5, ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "ART ≥3y", x = "Survey year", y = "") + 
  theme_bw()

H <- ggAcf(yeschild5$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "With <5y child", x = "Survey year", y = "") + 
  theme_bw()

I <- ggAcf(nochild5$residuals, type = "correlation", lag.max = 5, color = "red", ylim = c(-0.5, 0.5), size = 1) + 
  labs(title = "Without <5y child", x = "Survey year", y = "") + 
  theme_bw()

#=======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "FigS3_overall_carriage_ACF.png"),
       plot = (A | B | C | D | E)/(F | G | H | I ),
       width = 11, height = 6, unit="in", dpi = 300)

