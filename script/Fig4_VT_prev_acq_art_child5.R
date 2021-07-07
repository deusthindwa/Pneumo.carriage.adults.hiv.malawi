#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#----------------------------------------------------------------------------------

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
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
sart_year <- left_join(left_join(sart %>% group_by(year) %>% tally() %>% rename(Tot = n), 
sart %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
sart %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
A <- ggplot(sart_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 0.7, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"]) +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.6, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"]) +
  scale_y_continuous("VT carriage prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.45)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "ART <3y", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

B <- ggplot(sart_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 0.7, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"]) +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.6, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"]) +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.45)) + 
  labs(title = "", x = "Year") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

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
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
lart_year <- left_join(left_join(lart %>% group_by(year) %>% tally() %>% rename(Tot = n), 
lart %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
lart %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
C <- ggplot(lart_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 0.7, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"]) +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.6, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"]) +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.45)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "ART ≥3y", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

D <- ggplot(lart_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 0.7, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"]) +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.6, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"]) +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "Daily carriage acquisition"), limits = c(0, 0.45)) + 
  labs(title = "", x = "Year") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

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
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
nochild5_year <- left_join(left_join(nochild5 %>% group_by(year) %>% tally() %>% rename(Tot = n), 
nochild5 %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
nochild5 %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
E <- ggplot(nochild5_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 0.7, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"]) +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.6, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"]) +
  scale_y_continuous("VT carriage prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.45)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "With <5y child", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

F <- ggplot(nochild5_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 0.7, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"]) +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.6, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"]) +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.45)) + 
  labs(title = "", x = "Year") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

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
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
yeschild5_year <- left_join(left_join(yeschild5 %>% group_by(year) %>% tally() %>% rename(Tot = n), 
yeschild5 %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
yeschild5 %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
G <- ggplot(yeschild5_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 0.7, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"]) +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.6, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"]) +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.45)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "Without <5y child", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

H <- ggplot(yeschild5_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 0.7, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"]) +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.6, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"]) +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "Daily carriage acquisition"), limits = c(0, 0.45)) + 
  labs(title = "", x = "Year") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig4_VT_prev_acq_art_child5.tiff"),
       plot = (A | B | C | D) / (E | F | G | H),
       width = 13, height = 6, unit="in", dpi = 200)
