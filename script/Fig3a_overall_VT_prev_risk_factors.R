#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 17/12/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#===================================================================================

#subset for a dataset to store model estimates
female = filter(pcvpa.mod, sex == 1) %>% select(carr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_female = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 1), na.action = na.exclude)
female$fit = predict.gam(model_female, type = "response", se.fit = TRUE)$fit
female$fit_lci = model_female$family$linkinv(predict.gam(model_female, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_female, type = "link", se.fit = TRUE)$se.fit))
female$fit_uci = model_female$family$linkinv(predict.gam(model_female, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_female, type = "link", se.fit = TRUE)$se.fit))
female <- female %>% mutate(category = "Female", status = "Sex") 

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
male = filter(pcvpa.mod, sex == 2) %>% select(carr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_male = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 2), na.action = na.exclude)
male$fit = predict.gam(model_male, type = "response", se.fit = TRUE)$fit
male$fit_lci = model_male$family$linkinv(predict.gam(model_male, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_male, type = "link", se.fit = TRUE)$se.fit))
male$fit_uci = model_male$family$linkinv(predict.gam(model_male, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_male, type = "link", se.fit = TRUE)$se.fit))
male <- male %>% mutate(category = "Male", status = "Sex") 

#===================================================================================

#subset for a dataset to store model estimates
lses = filter(pcvpa.mod, sescat == 1) %>% select(carr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_lses = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 1), na.action = na.exclude)
lses$fit = predict.gam(model_lses, type = "response", se.fit = TRUE)$fit
lses$fit_lci = model_lses$family$linkinv(predict.gam(model_lses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lses, type = "link", se.fit = TRUE)$se.fit))
lses$fit_uci = model_lses$family$linkinv(predict.gam(model_lses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lses, type = "link", se.fit = TRUE)$se.fit))
lses <- lses %>% mutate(category = "Low SES", status = "SES")

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
hses = filter(pcvpa.mod, sescat == 2) %>% select(carr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_hses = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 2), na.action = na.exclude)
hses$fit = predict.gam(model_hses, type = "response", se.fit = TRUE)$fit
hses$fit_lci = model_hses$family$linkinv(predict.gam(model_hses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_hses, type = "link", se.fit = TRUE)$se.fit))
hses$fit_uci = model_hses$family$linkinv(predict.gam(model_hses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_hses, type = "link", se.fit = TRUE)$se.fit))
hses <- hses %>% mutate(category = "High SES", status = "SES")

#===================================================================================

#subset for a dataset to store model estimates
sart = filter(pcvpa.mod, artdur == 1) %>% select(carr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_sart = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 1), na.action = na.exclude)
sart$fit = predict.gam(model_sart, type = "response", se.fit = TRUE)$fit
sart$fit_lci = model_sart$family$linkinv(predict.gam(model_sart, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_sart, type = "link", se.fit = TRUE)$se.fit))
sart$fit_uci = model_sart$family$linkinv(predict.gam(model_sart, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_sart, type = "link", se.fit = TRUE)$se.fit))
sart <- sart %>% mutate(category = "ART <3y", status = "ART duration")

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
lart = filter(pcvpa.mod, artdur == 2) %>% select(carr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_lart = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 2), na.action = na.exclude)
lart$fit = predict.gam(model_lart, type = "response", se.fit = TRUE)$fit
lart$fit_lci = model_lart$family$linkinv(predict.gam(model_lart, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lart, type = "link", se.fit = TRUE)$se.fit))
lart$fit_uci = model_lart$family$linkinv(predict.gam(model_lart, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lart, type = "link", se.fit = TRUE)$se.fit))
lart <- lart %>% mutate(category = "ART ≥3y", status = "ART duration")

#===================================================================================

#subset for a dataset to store model estimates
nochild5 = filter(pcvpa.mod, nochild5 == 1) %>% select(carr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_nochild5 = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 1), na.action = na.exclude)
nochild5$fit = predict.gam(model_nochild5, type = "response", se.fit = TRUE)$fit
nochild5$fit_lci = model_nochild5$family$linkinv(predict.gam(model_nochild5, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_nochild5, type = "link", se.fit = TRUE)$se.fit))
nochild5$fit_uci = model_nochild5$family$linkinv(predict.gam(model_nochild5, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_nochild5, type = "link", se.fit = TRUE)$se.fit))
nochild5 <- nochild5 %>% mutate(category = "Without <5y child", status = "Cohabitation status")

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
yeschild5 = filter(pcvpa.mod, nochild5 == 2) %>% select(carr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_yeschild5 = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 2), na.action = na.exclude)
yeschild5$fit = predict.gam(model_yeschild5, type = "response", se.fit = TRUE)$fit
yeschild5$fit_lci = model_yeschild5$family$linkinv(predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$se.fit))
yeschild5$fit_uci = model_yeschild5$family$linkinv(predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$se.fit))
yeschild5 <- yeschild5 %>% mutate(category = "With <5y child", status = "Cohabitation status")

#===================================================================================

# mean values of fits

# sex
female %>% filter(carr == 1) %>% summarise(total = sum(carr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))
male %>% filter(carr == 1) %>% summarise(total = sum(carr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))

# SES
lses %>% filter(carr == 1) %>% summarise(total = sum(carr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))
hses %>% filter(carr == 1) %>% summarise(total = sum(carr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))

# ART
sart %>% filter(carr == 1) %>% summarise(total = sum(carr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))
lart %>% filter(carr == 1) %>% summarise(total = sum(carr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))

# U5 child
nochild5 %>% filter(carr == 1) %>% summarise(total = sum(carr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))
yeschild5 %>% filter(carr == 1) %>% summarise(total = sum(carr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))

#===================================================================================

# combined datasets for interaction test
interaction <- rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))

# test for interaction between age group and independent variable
summary(aov(fit ~ agegp * category, data = interaction %>% filter(status == "Sex")))
summary(aov(fit ~ agegp * category, data = interaction %>% filter(status == "SES")))
summary(aov(fit ~ agegp * category, data = interaction %>% filter(status == "ART duration")))
summary(aov(fit ~ agegp * category, data = interaction %>% filter(status == "Cohabitation status")))

# test for interaction between time and independent variable
summary(aov(fit ~ year * category, data = interaction %>% filter(status == "Sex")))
summary(aov(fit ~ year * category, data = interaction %>% filter(status == "SES")))
summary(aov(fit ~ year * category, data = interaction %>% filter(status == "ART duration")))
summary(aov(fit ~ year * category, data = interaction %>% filter(status == "Cohabitation status")))

#===================================================================================
#===================================================================================

#subset for a dataset to store model estimates
female = filter(pcvpa.mod, sex == 1) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_female = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 1), na.action = na.exclude)
female$fit = predict.gam(model_female, type = "response", se.fit = TRUE)$fit
female$fit_lci = model_female$family$linkinv(predict.gam(model_female, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_female, type = "link", se.fit = TRUE)$se.fit))
female$fit_uci = model_female$family$linkinv(predict.gam(model_female, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_female, type = "link", se.fit = TRUE)$se.fit))
female <- female %>% mutate(category = "Female", status = "Sex") 

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
male = filter(pcvpa.mod, sex == 2) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_male = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 2), na.action = na.exclude)
male$fit = predict.gam(model_male, type = "response", se.fit = TRUE)$fit
male$fit_lci = model_male$family$linkinv(predict.gam(model_male, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_male, type = "link", se.fit = TRUE)$se.fit))
male$fit_uci = model_male$family$linkinv(predict.gam(model_male, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_male, type = "link", se.fit = TRUE)$se.fit))
male <- male %>% mutate(category = "Male", status = "Sex") 

#===================================================================================

#subset for a dataset to store model estimates
lses = filter(pcvpa.mod, sescat == 1) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_lses = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 1), na.action = na.exclude)
lses$fit = predict.gam(model_lses, type = "response", se.fit = TRUE)$fit
lses$fit_lci = model_lses$family$linkinv(predict.gam(model_lses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lses, type = "link", se.fit = TRUE)$se.fit))
lses$fit_uci = model_lses$family$linkinv(predict.gam(model_lses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lses, type = "link", se.fit = TRUE)$se.fit))
lses <- lses %>% mutate(category = "Low SES", status = "SES")

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
hses = filter(pcvpa.mod, sescat == 2) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_hses = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 2), na.action = na.exclude)
hses$fit = predict.gam(model_hses, type = "response", se.fit = TRUE)$fit
hses$fit_lci = model_hses$family$linkinv(predict.gam(model_hses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_hses, type = "link", se.fit = TRUE)$se.fit))
hses$fit_uci = model_hses$family$linkinv(predict.gam(model_hses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_hses, type = "link", se.fit = TRUE)$se.fit))
hses <- hses %>% mutate(category = "High SES", status = "SES")

#===================================================================================

#subset for a dataset to store model estimates
sart = filter(pcvpa.mod, artdur == 1) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_sart = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 1), na.action = na.exclude)
sart$fit = predict.gam(model_sart, type = "response", se.fit = TRUE)$fit
sart$fit_lci = model_sart$family$linkinv(predict.gam(model_sart, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_sart, type = "link", se.fit = TRUE)$se.fit))
sart$fit_uci = model_sart$family$linkinv(predict.gam(model_sart, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_sart, type = "link", se.fit = TRUE)$se.fit))
sart <- sart %>% mutate(category = "ART <3y", status = "ART duration")

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
lart = filter(pcvpa.mod, artdur == 2) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_lart = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 2), na.action = na.exclude)
lart$fit = predict.gam(model_lart, type = "response", se.fit = TRUE)$fit
lart$fit_lci = model_lart$family$linkinv(predict.gam(model_lart, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lart, type = "link", se.fit = TRUE)$se.fit))
lart$fit_uci = model_lart$family$linkinv(predict.gam(model_lart, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lart, type = "link", se.fit = TRUE)$se.fit))
lart <- lart %>% mutate(category = "ART ≥3y", status = "ART duration")

#===================================================================================

#subset for a dataset to store model estimates
nochild5 = filter(pcvpa.mod, nochild5 == 1) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_nochild5 = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 1), na.action = na.exclude)
nochild5$fit = predict.gam(model_nochild5, type = "response", se.fit = TRUE)$fit
nochild5$fit_lci = model_nochild5$family$linkinv(predict.gam(model_nochild5, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_nochild5, type = "link", se.fit = TRUE)$se.fit))
nochild5$fit_uci = model_nochild5$family$linkinv(predict.gam(model_nochild5, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_nochild5, type = "link", se.fit = TRUE)$se.fit))
nochild5 <- nochild5 %>% mutate(category = "Without <5y child", status = "Cohabitation status")

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
yeschild5 = filter(pcvpa.mod, nochild5 == 2) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_yeschild5 = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 2), na.action = na.exclude)
yeschild5$fit = predict.gam(model_yeschild5, type = "response", se.fit = TRUE)$fit
yeschild5$fit_lci = model_yeschild5$family$linkinv(predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$se.fit))
yeschild5$fit_uci = model_yeschild5$family$linkinv(predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$se.fit))
yeschild5 <- yeschild5 %>% mutate(category = "With <5y child", status = "Cohabitation status")

#===================================================================================

# mean values of fits

# sex
female %>% filter(vtcarr == 1) %>% summarise(total = sum(vtcarr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))
male %>% filter(vtcarr == 1) %>% summarise(total = sum(vtcarr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))

# SES
lses %>% filter(vtcarr == 1) %>% summarise(total = sum(vtcarr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))
hses %>% filter(vtcarr == 1) %>% summarise(total = sum(vtcarr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))

# ART
sart %>% filter(vtcarr == 1) %>% summarise(total = sum(vtcarr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))
lart %>% filter(vtcarr == 1) %>% summarise(total = sum(vtcarr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))

# U5 child
nochild5 %>% filter(vtcarr == 1) %>% summarise(total = sum(vtcarr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))
yeschild5 %>% filter(vtcarr == 1) %>% summarise(total = sum(vtcarr), mfit = mean(fit), lfit = mean(fit_lci), ufit = mean(fit_uci))

#===================================================================================

# combined datasets for interaction test
interaction <- rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))


# test for interaction between age group and independent variable
summary(aov(fit ~ agegp * category, data = interaction %>% filter(status == "Sex")))
summary(aov(fit ~ agegp * category, data = interaction %>% filter(status == "SES")))
summary(aov(fit ~ agegp * category, data = interaction %>% filter(status == "ART duration")))
summary(aov(fit ~ agegp * category, data = interaction %>% filter(status == "Cohabitation status")))

# test for interaction between age group and independent variable
summary(aov(fit ~ year * category, data = interaction %>% filter(status == "Sex")))
summary(aov(fit ~ year * category, data = interaction %>% filter(status == "SES")))
summary(aov(fit ~ year * category, data = interaction %>% filter(status == "ART duration")))
summary(aov(fit ~ year * category, data = interaction %>% filter(status == "Cohabitation status")))

#===================================================================================


