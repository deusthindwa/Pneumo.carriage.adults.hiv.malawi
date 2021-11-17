#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 20/11/2021

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

#create age group and ggplot
A <- left_join(
left_join(
rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_))))) %>%
  group_by(agegp, status, category) %>% tally() %>% rename(Tot = n),


rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_))))) %>%
  filter(carr == 1) %>% group_by(agegp, status, category) %>% tally() %>% rename(Pos = n)),


rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_))))) %>%
  filter(carr != 0) %>% group_by(agegp, status, category) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% ungroup()) %>%
  group_by(category) %>%
  mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8],
         category = factor(category, levels = c('Male', 'Female', 'Low SES', 'High SES', 'ART <3y', 'ART ≥3y', 'With <5y child', 'Without <5y child')),
         status = factor(status, levels = c('Sex', 'SES', 'ART duration', 'Cohabitation status'))) %>%
  
  #plot prevalence curves
  ggplot() + 
  geom_point(aes(x = agegp, y = obs, size = Pos, color = category), shape = 1, position=position_dodge(width=0.05)) + 
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci, color = category), width = 0, size = 0.3, position=position_dodge(width=0.05)) + 
  geom_line(aes(x = agegp, y = fit, group = category, color = category), size = 1) + 
  geom_ribbon(aes(x = agegp, y = fit, group = category, fill = category, color = category, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, size = 0.1) + 
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +  
  scale_x_discrete(expand = c(0.05,0.05)) +
  labs(title = "", x = NULL, y = "Overall carriage prevalence") +
  facet_grid(.~status) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16), strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept=0.05, linetype = "dashed", color = "black", size = 0.2) +
  theme(axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  theme(legend.position = "top") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#===================================================================================

#create age group and ggplot
B <- left_join(
  left_join(
    rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
      group_by(year, status, category) %>% tally() %>% rename(Tot = n),
    
    rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
      filter(carr == 1) %>% group_by(year, status, category) %>% tally() %>% rename(Pos = n)),
  
  rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
    filter(carr != 0) %>% group_by(year, status, category) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% ungroup()) %>%
  group_by(category) %>%
  mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10],
         category = factor(category, levels = c('Male', 'Female', 'Low SES', 'High SES', 'ART <3y', 'ART ≥3y', 'With <5y child', 'Without <5y child')),
         status = factor(status, levels = c('Sex', 'SES', 'ART duration', 'Cohabitation status'))) %>%
  
  #plot prevalence curves
  ggplot() + 
  geom_point(aes(x = year, y = obs, size = Pos, color = category), shape = 1, position=position_dodge(width=0.05)) + 
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci, color = category), width = 0, size = 0.3, position=position_dodge(width=0.05)) + 
  geom_line(aes(x = year, y = fit, group = category, color = category), size = 1) + 
  geom_ribbon(aes(x = year, y = fit, group = category, fill = category, color = category, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, size = 0.1) + 
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +  
  scale_x_discrete(expand = c(0.05,0.05)) +
  labs(title = "", x = NULL, y = "") +
  facet_grid(.~status) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16), strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept=0.05, linetype = "dashed", color = "black", size = 0.2) +
  theme(axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.text.y=element_blank()) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

#======================================================================================
#======================================================================================

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

#create age group and ggplot
C <- left_join(
  left_join(
    rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
      mutate(agegp = if_else(age <=24, "18-24",
                             if_else(age >24 & age <=29, "25-29",
                                     if_else(age >29 & age <=34, "30-34",
                                             if_else(age >34 & age <=40, "35-40", NA_character_))))) %>%
      group_by(agegp, status, category) %>% tally() %>% rename(Tot = n),
    
    
    rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
      mutate(agegp = if_else(age <=24, "18-24",
                             if_else(age >24 & age <=29, "25-29",
                                     if_else(age >29 & age <=34, "30-34",
                                             if_else(age >34 & age <=40, "35-40", NA_character_))))) %>%
      filter(vtcarr == 1) %>% group_by(agegp, status, category) %>% tally() %>% rename(Pos = n)),
  
  
  rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
    mutate(agegp = if_else(age <=24, "18-24",
                           if_else(age >24 & age <=29, "25-29",
                                   if_else(age >29 & age <=34, "30-34",
                                           if_else(age >34 & age <=40, "35-40", NA_character_))))) %>%
    filter(vtcarr != 0) %>% group_by(agegp, status, category) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% ungroup()) %>%
  group_by(category) %>%
  mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8],
         category = factor(category, levels = c('Male', 'Female', 'Low SES', 'High SES', 'ART <3y', 'ART ≥3y', 'With <5y child', 'Without <5y child')),
         status = factor(status, levels = c('Sex', 'SES', 'ART duration', 'Cohabitation status'))) %>%
  
  #plot prevalence curves
  ggplot() + 
  geom_point(aes(x = agegp, y = obs, size = Pos, color = category), shape = 1, position=position_dodge(width=0.05)) + 
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci, color = category), width = 0, size = 0.3, position=position_dodge(width=0.05)) + 
  geom_line(aes(x = agegp, y = fit, group = category, color = category), size = 1) + 
  geom_ribbon(aes(x = agegp, y = fit, group = category, fill = category, color = category, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, size = 0.1) + 
  coord_cartesian(ylim = c(0, 0.45)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +  
  scale_x_discrete(expand = c(0.05,0.05)) +
  labs(title = "", x = "Age group (years)", y = "VT(+st3) carr1iage prevalence") +
  facet_grid(.~status) +
  theme_bw() + 
  theme(strip.text.x = element_blank(), strip.text.y = element_blank(), strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept=0.05, linetype = "dashed", color = "black", size = 0.2) +
  theme(axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(plot.title = element_text(size = 20), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) + 
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) 

#===================================================================================

#create age group and ggplot
D <- left_join(
  left_join(
    rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
      group_by(year, status, category) %>% tally() %>% rename(Tot = n),
    
    rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
      filter(vtcarr == 1) %>% group_by(year, status, category) %>% tally() %>% rename(Pos = n)),
  
  rbind(female, male, lses, hses, sart, lart, nochild5, yeschild5) %>% 
    filter(vtcarr != 0) %>% group_by(year, status, category) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% ungroup()) %>%
  group_by(category) %>%
  mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10],
         category = factor(category, levels = c('Male', 'Female', 'Low SES', 'High SES', 'ART <3y', 'ART ≥3y', 'With <5y child', 'Without <5y child')),
         status = factor(status, levels = c('Sex', 'SES', 'ART duration', 'Cohabitation status'))) %>%
  
  #plot prevalence curves
  ggplot() + 
  geom_point(aes(x = year, y = obs, size = Pos, color = category), shape = 1, position=position_dodge(width=0.05)) + 
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci, color = category), width = 0, size = 0.3, position=position_dodge(width=0.05)) + 
  geom_line(aes(x = year, y = fit, group = category, color = category), size = 1) + 
  geom_ribbon(aes(x = year, y = fit, group = category, fill = category, color = category, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, size = 0.1) + 
  coord_cartesian(ylim = c(0, 0.45)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +  
  #scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "", x = "Time (years)", y = "") +
  facet_grid(.~status) +
  theme_bw() + 
  theme(strip.text.x = element_blank(), strip.text.y = element_blank(), strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept=0.05, linetype = "dashed", color = "black", size = 0.2) +
  theme(axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(plot.title = element_text(size = 20), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y=element_blank()) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig3_Overall_VT_prev_risk_factors.png"),
       plot = ((A | B)/(C | D)),
       width = 24, height = 12, unit="in", dpi = 300)

