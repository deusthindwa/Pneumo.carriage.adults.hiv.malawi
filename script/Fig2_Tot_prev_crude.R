#Written by Deus Thindwa
#Pneumococcal carriage prevalence in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)
#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
crude = pcvpa.mod %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_crude = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod, na.action = na.exclude)
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

#get age group mean predicted prevalence
crude_age <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], status = "VT (with Serotype 3)")

#----------------------------------------------------------------------------------

#get year group mean predicted prevalence
crude_year <- left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], status = "VT (with Serotype 3)")

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
crude = pcvpa.mod %>% select(carr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_crude = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod, na.action = na.exclude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

summary(model_crude)
intervals(model_crude)
#----------------------------------------------------------------------------------

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24",
                                   if_else(age >24 & age <=29, "25-29",
                                           if_else(age >29 & age <=34, "30-34",
                                                   if_else(age >34 & age <=40, "35-40", NA_character_))))))

#get age group mean predicted prevalence
crude_age2 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
                                  crude %>% filter(carr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
                        crude %>% filter(carr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
                         ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], status = "Overall")

#----------------------------------------------------------------------------------

#get year group mean predicted prevalence
crude_year2 <- left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
                                  crude %>% filter(carr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
                        crude %>% filter(carr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
                          ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], status = "Overall")

#----------------------------------------------------------------------------------

#plot prevalence curves
A <- ggplot(rbind(crude_age, crude_age2)) + 
  geom_point(aes(x = agegp, y = obs, size = Pos, color = status), shape = 1) + 
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci, color = status), width = 0, size = 0.3) + 
  geom_line(aes(x = agegp, y = fit, group = status, color = status), size = 0.7) + 
  geom_ribbon(aes(x = agegp, y = fit, group = status, fill = status, ymin = fit_lci, ymax = fit_uci), alpha = 0.2) + 
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +   scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "A", x = "Age group (years)", y = "Carriage prevalence") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) + 
  theme(legend.position = "bottom") +
  guides(color=guide_legend(title="Carriage:"), size = "none", group = "none", fill = "none")

B <- ggplot(rbind(crude_year, crude_year2)) + 
  geom_point(aes(x = year, y = obs, size = Pos, color = status), shape = 1) + 
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci, color = status), width = 0, size = 0.3) + 
  geom_line(aes(x = year, y = fit, group = status, color = status), size = 0.7) + 
  geom_ribbon(aes(x = year, y = fit, group = status, fill = status, ymin = fit_lci, ymax = fit_uci), alpha = 0.2) + 
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) + 
  labs(title = "B", x = "Year", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) + 
  theme(legend.position = "bottom") +
  guides(color=guide_legend(title=""), size = "none", group = "none", fill = "none")

#======================================================================================


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

#get age group mean predicted prevalence
crude_age <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
                                 crude %>% filter(vtcarr1 == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
                       crude %>% filter(vtcarr1 != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
                         ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], status = "VT (without Serotype 3)")

#----------------------------------------------------------------------------------

#get year group mean predicted prevalence
crude_year <- left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
                                  crude %>% filter(vtcarr1 == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
                        crude %>% filter(vtcarr1 != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
                          ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], status = "VT (without Serotype 3)")

#----------------------------------------------------------------------------------

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

#get age group mean predicted prevalence
crude_age2 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
                                  crude %>% filter(carr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
                        crude %>% filter(carr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
                          ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], status = "Overall")

#----------------------------------------------------------------------------------

#get year group mean predicted prevalence
crude_year2 <- left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
                                   crude %>% filter(carr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
                         crude %>% filter(carr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
                           ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], status = "Overall")

#----------------------------------------------------------------------------------

#plot prevalence curves
C <- ggplot(rbind(crude_age, crude_age2)) + 
  geom_point(aes(x = agegp, y = obs, size = Pos, color = status), shape = 1) + 
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci, color = status), width = 0, size = 0.3) + 
  geom_line(aes(x = agegp, y = fit, group = status, color = status), size = 0.7) + 
  geom_ribbon(aes(x = agegp, y = fit, group = status, fill = status, ymin = fit_lci, ymax = fit_uci), alpha = 0.2) + 
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +   scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "C", x = "Age group (years)", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) + 
  theme(legend.position = "bottom") +
  guides(color=guide_legend(title=""), size = "none", group = "none", fill = "none")

D <- ggplot(rbind(crude_year, crude_year2)) + 
  geom_point(aes(x = year, y = obs, size = Pos, color = status), shape = 1) + 
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci, color = status), width = 0, size = 0.3) + 
  geom_line(aes(x = year, y = fit, group = status, color = status), size = 0.7) + 
  geom_ribbon(aes(x = year, y = fit, group = status, fill = status, ymin = fit_lci, ymax = fit_uci), alpha = 0.2) + 
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) + 
  labs(title = "D", x = "Year", y = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) + 
  theme(legend.position = "bottom") +
  guides(color=guide_legend(title=""), size = "none", group = "none", fill = "none")

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig2_VT_prev_acq_crude.png"),
       plot = (A | B | C | D),
       width = 17, height = 5, unit="in", dpi = 200)


#======================================================================================

## Generate a table of cell frequencies. First set the levels of the outcome
## and the exposure so the frequencies in the 2 by 2 table come out in the
## conventional format:
dat1$low <- factor(dat1$low, levels = c(1,0))
dat1$smoke <- factor(dat1$smoke, levels = c(1,0))
dat1$race <- factor(dat1$race, levels = c(1,2,3))
## Generate the 2 by 2 table. Exposure (rows) = smoke. Outcome (columns) = low.
tab1 <- table(dat1$smoke, dat1$low, dnn = c("Smoke", "Low BW"))
print(tab1)
## Compute the incidence risk ratio and other measures of association:
epi.2by2(dat = tab1, method = "cohort.count", 
         conf.level = 0.95, units = 100, outcome = "as.columns")