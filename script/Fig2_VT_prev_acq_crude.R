#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
crude = pcvpa.mod %>% select(vtcarr, age, surv, sex, nochild5)

#fit model to individual trajectories & obtain predictions and 95%CI
model_crude = gam(vtcarr ~ te(age, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = pcvpa.mod, na.action = na.exclude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
crude <- crude %>% 
  mutate(agegp = if_else(age <=20, "18-20",
                         if_else(age >20 & age <=29, "21-29",
                                 if_else(age >29 & age <=35, "30-35",
                                         if_else(age >35 & age <=40, "36-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
crude_age <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#create survey group
crude <- crude %>% 
  mutate(survgp = if_else(surv == 1 | surv == 2, "1-2",
                         if_else(surv == 3 | surv == 4, "3-4",
                                 if_else(surv == 5 | surv == 6, "5-6",
                                         if_else(surv == 7 | surv == 8, "7-8", NA_character_)))))

#join observed and predicted datasets for survey number
crude_surv <- left_join(left_join(crude %>% group_by(survgp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(survgp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(survgp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
A <- ggplot(crude_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"], color = "gray") +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"], color = "gray") +
  scale_y_continuous("Carriage prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.35)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "Vaccine Serotypes (VT)", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

B <- ggplot(crude_surv) +
  geom_point(aes(x = survgp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(survgp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = survgp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = survgp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"], color = "gray") +
  geom_line(aes(x = survgp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = survgp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.35)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "", x = "Survey group") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a dataset to store model estimates
crude = pcvpa.mod %>% select(vtcarr1, age, surv, sex, nochild5)

#fit model to individual trajectories & obtain predictions and 95%CI
model_crude = gam(vtcarr1 ~ te(age, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = pcvpa.mod, na.action = na.exclude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
crude <- crude %>% 
  mutate(agegp = if_else(age <=20, "18-20",
                         if_else(age >20 & age <=29, "21-29",
                                 if_else(age >29 & age <=35, "30-35",
                                         if_else(age >35 & age <=40, "36-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
crude_age <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
                                 crude %>% filter(vtcarr1 == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
                       crude %>% filter(vtcarr1 != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
                         ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#create survey group
crude <- crude %>% 
  mutate(survgp = if_else(surv == 1 | surv == 2, "1-2",
                          if_else(surv == 3 | surv == 4, "3-4",
                                  if_else(surv == 5 | surv == 6, "5-6",
                                          if_else(surv == 7 | surv == 8, "7-8", NA_character_)))))

#join observed and predicted datasets for survey number
crude_surv <- left_join(left_join(crude %>% group_by(survgp) %>% tally() %>% rename(Tot = n), 
                                  crude %>% filter(vtcarr1 == 1) %>% group_by(survgp) %>% tally() %>% rename(Pos = n)), 
                        crude %>% filter(vtcarr1 != 0) %>% group_by(survgp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
                          ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
C <- ggplot(crude_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"], color = "gray") +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.35)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "VT without Serotype 3", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

D <- ggplot(crude_surv) +
  geom_point(aes(x = survgp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(survgp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = survgp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = survgp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"], color = "gray") +
  geom_line(aes(x = survgp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = survgp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "Daily carriage acquisition"), limits = c(0, 0.35)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "", x = "Survey group") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig2_VT_prev_acq_crude.tiff"),
       plot = (A | B | C | D),
       width = 12, height = 3, unit="in", dpi = 200)
