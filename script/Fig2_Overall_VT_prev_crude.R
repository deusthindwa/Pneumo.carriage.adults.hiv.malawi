#Written by Deus Thindwa
#Pneumococcal carriage prevalence in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 20/11/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#=================================================================================

#VT CARRIAGE WITH SEROTYPE 3

#subset for a dataset to store model estimates
crude = pcvpa.mod %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_crude = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod, na.action = na.exclude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age groups
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                                 if_else(age >34 & age <=40, "35-40", NA_character_))))))

#get predicted mean prevalence for each age group 
crude_age1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], status = "VT(+st3) carriage")

#----------------------------------------------------------------------------------

#get predicted mean prevalence for each year 
crude_year1 <- left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], status = "VT(+st3) carriage")

#=================================================================================

#VT CARRIAGE WITHOUT SEROTYPE 3

#subset for a dataset to store model estimates
crude = pcvpa.mod %>% select(vtcarr1, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_crude = gam(vtcarr1 ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod, na.action = na.exclude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24",
                                   if_else(age >24 & age <=29, "25-29",
                                           if_else(age >29 & age <=34, "30-34",
                                                   if_else(age >34 & age <=40, "35-40", NA_character_))))))

#get predicted mean prevalence for each age group
crude_age2 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr1 == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr1 != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], status = "VT(-st3) carriage")

#----------------------------------------------------------------------------------

#get predicted mean prevalence for each year
crude_year2 <- left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr1 == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr1 != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], status = "VT(-st3) carriage")

#=================================================================================

#OVERALL CARRIAGE

#subset for a dataset to store model estimates
crude = pcvpa.mod %>% select(carr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_crude = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod, na.action = na.exclude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24",
                                   if_else(age >24 & age <=29, "25-29",
                                           if_else(age >29 & age <=34, "30-34",
                                                   if_else(age >34 & age <=40, "35-40", NA_character_))))))

#get predicted mean prevalence for each age group
crude_age3 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(carr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(carr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], status = "Overall carriage")

#----------------------------------------------------------------------------------

#get predicted mean prevalence for each year
crude_year3 <- left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
crude %>% filter(carr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(carr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], status = "Overall carriage")

#=================================================================================

#plot prevalence curves
A <- ggplot(rbind(crude_age1, crude_age2, crude_age3)) + 
  geom_point(aes(x = agegp, y = obs, size = Pos, color = status), shape = 1, position=position_dodge(width=0.05)) + 
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci, color = status), width = 0, size = 0.3, position=position_dodge(width=0.05)) + 
  geom_line(aes(x = agegp, y = fit, group = status, color = status), size = 1) + 
  geom_ribbon(aes(x = agegp, y = fit, group = status, fill = status, color = status, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, size = 0.1) + 
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +  
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "A", x = "Age group (years)", y = "Carriage prevalence") +
  theme_bw() + 
  geom_hline(yintercept=0.05, linetype = "dashed", color = "black", size = 0.2) +
  theme(axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(plot.title = element_text(size = 20), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

B <- ggplot(rbind(crude_year1, crude_year2, crude_year3)) + 
  geom_point(aes(x = year, y = obs, size = Pos, color = status), shape = 1, position=position_dodge(width=0.05)) + 
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci, color = status), width = 0, size = 0.3, position=position_dodge(width=0.05)) + 
  geom_line(aes(x = year, y = fit, group = status, color = status), size = 1) + 
  geom_ribbon(aes(x = year, y = fit, group = status, fill = status, color = status, ymin = fit_lci, ymax = fit_uci), alpha = 0.3, size = 0.1) + 
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) + 
  labs(title = "B", x = "Time (years)", y = "") +
  theme_bw() + 
  geom_hline(yintercept=0.05, linetype = "dashed", color = "black", size = 0.2) +
  theme(axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(plot.title = element_text(size = 20), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y=element_blank()) + 
  theme(legend.position = "right") +
  guides(color = guide_legend(title = "Pneumococcal carriage"), size = guide_legend(title = "Sample size"), group = "none", fill = "none") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

#=================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig2_Overall_VT_prev_crude.png"),
       plot = (A | B),
       width = 10, height = 5, unit="in", dpi = 200)
