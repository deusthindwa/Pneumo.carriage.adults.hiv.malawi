#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
male = filter(pcvpa.mod, sex == 1) %>% select(vtcarr, age, year, sex, nochild5, seas)

#fit model to individual trajectories & obtain predictions and 95%CI
model_male = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + nochild5 + seas, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 1), na.action = na.exclude)
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
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
male_year <- left_join(left_join(male %>% group_by(year) %>% tally() %>% rename(Tot = n), 
male %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
male %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
A <- ggplot(male_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"], color = "gray") +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"], color = "gray") +
  scale_y_continuous("VT carriage prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.65)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "Male", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

B <- ggplot(male_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 1, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"], color = "gray") +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.65)) + 
  labs(title = "", x = "Year") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a dataset to store model estimates
female = filter(pcvpa.mod, sex == 0) %>% select(vtcarr, age, year, sex, nochild5, seas)

#fit model to individual trajectories & obtain predictions and 95%CI
model_female = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + nochild5 + seas, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 0), na.action = na.exclude)
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
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
female_year <- left_join(left_join(female %>% group_by(year) %>% tally() %>% rename(Tot = n), 
female %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
female %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
C <- ggplot(female_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"], color = "gray") +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.65)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "Female", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

D <- ggplot(female_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 1, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"], color = "gray") +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "Daily carriage acquisition"), limits = c(0, 0.65)) + 
  labs(title = "", x = "Year") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a dataset to store model estimates
lses = filter(pcvpa.mod, sescat == 0) %>% select(vtcarr, age, year, sex, nochild5, seas)

#fit model to individual trajectories & obtain predictions and 95%CI
model_lses = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + sex + nochild5 + seas, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 0), na.action = na.exclude)
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
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
lses_year <- left_join(left_join(lses %>% group_by(year) %>% tally() %>% rename(Tot = n), 
lses %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
lses %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
E <- ggplot(lses_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"], color = "gray") +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"], color = "gray") +
  scale_y_continuous("VT carriage prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.65)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "Low SES", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

F <- ggplot(lses_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 1, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"], color = "gray") +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.65)) + 
  labs(title = "", x = "Year") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a dataset to store model estimates
hses = filter(pcvpa.mod, sescat == 1) %>% select(vtcarr, age, year, sex, nochild5, seas)

#fit model to individual trajectories & obtain predictions and 95%CI
model_hses = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + sex + nochild5 + seas, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 1), na.action = na.exclude)
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
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
hses_year <- left_join(left_join(hses %>% group_by(year) %>% tally() %>% rename(Tot = n), 
hses %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
hses %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
G <- ggplot(hses_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"], color = "gray") +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.65)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "Middle/High SES", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

H <- ggplot(hses_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 1, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"], color = "gray") +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "Daily carriage acquisition"), limits = c(0, 0.65)) + 
  labs(title = "", x = "Year") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig3_VT_prev_acq_sex_ses.tiff"),
       plot = (A | B | C | D) / (E | F | G | H),
       width = 13, height = 6, unit="in", dpi = 200)
