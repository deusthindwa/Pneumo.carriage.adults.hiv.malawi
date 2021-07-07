#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#======================================================================================

#subset for a dataset to store model estimates
male = filter(pcvpa.mod, sex == 2) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_male = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 2), na.action = na.exclude)
male$fit = predict.gam(model_male, type = "response", se.fit = TRUE)$fit
male$fit_lci = model_male$family$linkinv(predict.gam(model_male, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_male, type = "link", se.fit = TRUE)$se.fit))
male$fit_uci = model_male$family$linkinv(predict.gam(model_male, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_male, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
male <- male %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
male_age <- left_join(left_join(male %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
male %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
male %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Male")

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
male_year <- left_join(left_join(male %>% group_by(year) %>% tally() %>% rename(Tot = n), 
male %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
male %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Male")

#======================================================================================

#subset for a dataset to store model estimates
female = filter(pcvpa.mod, sex == 1) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_female = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 1), na.action = na.exclude)
female$fit = predict.gam(model_female, type = "response", se.fit = TRUE)$fit
female$fit_lci = model_female$family$linkinv(predict.gam(model_female, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_female, type = "link", se.fit = TRUE)$se.fit))
female$fit_uci = model_female$family$linkinv(predict.gam(model_female, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_female, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
female <- female %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
female_age <- left_join(left_join(female %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
female %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
female %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Female")

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
female_year <- left_join(left_join(female %>% group_by(year) %>% tally() %>% rename(Tot = n), 
female %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
female %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Female")

#======================================================================================

#subset for a dataset to store model estimates
lses = filter(pcvpa.mod, sescat == 1) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_lses = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 1), na.action = na.exclude)
lses$fit = predict.gam(model_lses, type = "response", se.fit = TRUE)$fit
lses$fit_lci = model_lses$family$linkinv(predict.gam(model_lses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lses, type = "link", se.fit = TRUE)$se.fit))
lses$fit_uci = model_lses$family$linkinv(predict.gam(model_lses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lses, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
lses <- lses %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
lses_age <- left_join(left_join(lses %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
lses %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
lses %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Low SES")

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
lses_year <- left_join(left_join(lses %>% group_by(year) %>% tally() %>% rename(Tot = n), 
lses %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
lses %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Low SES")

#======================================================================================

#subset for a dataset to store model estimates
hses = filter(pcvpa.mod, sescat == 2) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_hses = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 2), na.action = na.exclude)
hses$fit = predict.gam(model_hses, type = "response", se.fit = TRUE)$fit
hses$fit_lci = model_hses$family$linkinv(predict.gam(model_hses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_hses, type = "link", se.fit = TRUE)$se.fit))
hses$fit_uci = model_hses$family$linkinv(predict.gam(model_hses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_hses, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
hses <- hses %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
hses_age <- left_join(left_join(hses %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
hses %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
hses %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Middle/High SES")

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
hses_year <- left_join(left_join(hses %>% group_by(year) %>% tally() %>% rename(Tot = n), 
hses %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
hses %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Middle/High SES")

#======================================================================================

#subset for a dataset to store model estimates
sart = filter(pcvpa.mod, artdur == 1) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_sart = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 1), na.action = na.exclude)
sart$fit = predict.gam(model_sart, type = "response", se.fit = TRUE)$fit
sart$fit_lci = model_sart$family$linkinv(predict.gam(model_sart, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_sart, type = "link", se.fit = TRUE)$se.fit))
sart$fit_uci = model_sart$family$linkinv(predict.gam(model_sart, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_sart, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
sart <- sart %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
sart_age <- left_join(left_join(sart %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
sart %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
sart %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "ART <3y")

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
sart_year <- left_join(left_join(sart %>% group_by(year) %>% tally() %>% rename(Tot = n), 
sart %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
sart %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "ART <3y")

#======================================================================================

#subset for a dataset to store model estimates
lart = filter(pcvpa.mod, artdur == 2) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_lart = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 2), na.action = na.exclude)
lart$fit = predict.gam(model_lart, type = "response", se.fit = TRUE)$fit
lart$fit_lci = model_lart$family$linkinv(predict.gam(model_lart, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lart, type = "link", se.fit = TRUE)$se.fit))
lart$fit_uci = model_lart$family$linkinv(predict.gam(model_lart, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lart, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
lart <- lart %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
lart_age <- left_join(left_join(lart %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
lart %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
lart %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "ART ≥3y")

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
lart_year <- left_join(left_join(lart %>% group_by(year) %>% tally() %>% rename(Tot = n), 
lart %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
lart %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "ART ≥3y")

#======================================================================================

#subset for a dataset to store model estimates
nochild5 = filter(pcvpa.mod, nochild5 == 1) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_nochild5 = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 1), na.action = na.exclude)
nochild5$fit = predict.gam(model_nochild5, type = "response", se.fit = TRUE)$fit
nochild5$fit_lci = model_nochild5$family$linkinv(predict.gam(model_nochild5, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_nochild5, type = "link", se.fit = TRUE)$se.fit))
nochild5$fit_uci = model_nochild5$family$linkinv(predict.gam(model_nochild5, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_nochild5, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
nochild5 <- nochild5 %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
nochild5_age <- left_join(left_join(nochild5 %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
nochild5 %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
nochild5 %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Living with <5y child")

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
nochild5_year <- left_join(left_join(nochild5 %>% group_by(year) %>% tally() %>% rename(Tot = n), 
nochild5 %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
nochild5 %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Living with <5y child")

#======================================================================================

#subset for a dataset to store model estimates
yeschild5 = filter(pcvpa.mod, nochild5 == 2) %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_yeschild5 = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, nochild5 == 2), na.action = na.exclude)
yeschild5$fit = predict.gam(model_yeschild5, type = "response", se.fit = TRUE)$fit
yeschild5$fit_lci = model_yeschild5$family$linkinv(predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$se.fit))
yeschild5$fit_uci = model_yeschild5$family$linkinv(predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_yeschild5, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
yeschild5 <- yeschild5 %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
yeschild5_age <- left_join(left_join(yeschild5 %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
yeschild5 %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
yeschild5 %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Living without <5y child")

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
yeschild5_year <- left_join(left_join(yeschild5 %>% group_by(year) %>% tally() %>% rename(Tot = n), 
yeschild5 %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
yeschild5 %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42, rf = "Living without <5y child")

#======================================================================================

pcvpa_age <- rbind(male_age, female_age, lses_age, hses_age, sart_age, lart_age, nochild5_age, yeschild5_age)
pcvpa_year <- rbind(male_year, female_year, lses_year, hses_year, sart_year, lart_year, nochild5_year, yeschild5_year)

#plot prevalence curves
A <- pcvpa_age %>% mutate(rff = factor(rf, levels=c('Male', 'Female', 'Low SES', 'Middle/High SES', 'ART <3y', 'ART ≥3y', 'Living with <5y child', 'Living without <5y child'))) %>%
  ggplot() +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1, color = rf), size = 0.7) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci, fill = rf), alpha = 0.2) +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1, color = rf), lty = "dashed", size = 0.6) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1, fill = rf), alpha = 0.2) +
  scale_y_continuous("VT carriage prevalence", sec.axis = sec_axis(~. * 0.1, name = "Daily VT carriage acquisition"), limits = c(0, 0.45)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(x = "Age group (years)") +
  theme_bw() +
  facet_grid(.~rff) +
  theme(strip.text.x = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 7, face = "bold"), axis.text.y = element_text(size = 7, face = "bold")) +
  theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8)) +
  theme(legend.position = "none")

B <- pcvpa_year %>% mutate(rff = factor(rf, levels=c('Male', 'Female', 'Low SES', 'Middle/High SES', 'ART <3y', 'ART ≥3y', 'Living with <5y child', 'Living without <5y child'))) %>%
  ggplot() +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1, color = rf), size = 0.7) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci, fill = rf), alpha = 0.2) +
  geom_line(aes(x = year, y = foi/0.1, group = 1, color = rf), lty = "dashed", size = 0.6) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1, fill = rf), alpha = 0.2) +
  scale_y_continuous("VT carriage prevalence", sec.axis = sec_axis(~. * 0.1, name = "Daily VT carriage acquisition"), limits = c(0, 0.40)) + 
  labs(x = "Survey year") +
  theme_bw() +
  facet_grid(.~rff) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.text.x = element_text(size = 7, face = "bold"), axis.text.y = element_text(size = 7, face = "bold")) +
  theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8)) +
  theme(legend.position = "none")

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig3_VT_prev_acq_risk_factors.tiff"),
       plot = (A / B ),
       width = 18, height = 6, unit="in", dpi = 200)

