#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#----------------------------------------------------------------------------------

#subset for a dataset to store model estimates
sart = filter(pcvpa.mod, artdur == 0) %>% select(vtcarr, age, year, sex, nochild5, seas)

#fit model to individual trajectories & obtain predictions and 95%CI
model_sart = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + nochild5 + seas, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 0), na.action = na.exclude)
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
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
A <- ggplot(sart_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"], color = "gray") +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"], color = "gray") +
  scale_y_continuous("VT carriage prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.8)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "ART <3y", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

B <- ggplot(sart_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 1, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"], color = "gray") +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.8)) + 
  labs(title = "", x = "Year") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a dataset to store model estimates
lart = filter(pcvpa.mod, artdur == 1) %>% select(vtcarr, age, year, sex, nochild5, seas)

#fit model to individual trajectories & obtain predictions and 95%CI
model_lart = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + nochild5 + seas, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 1), na.action = na.exclude)
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
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
C <- ggplot(lart_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"], color = "gray") +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.8)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "ART ≥3y", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

D <- ggplot(lart_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 1, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"], color = "gray") +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "Daily carriage acquisition"), limits = c(0, 0.8)) + 
  labs(title = "", x = "Year") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a dataset to store model estimates
lcd4 = filter(pcvpa.mod, cd4cnt == 0) %>% select(vtcarr, age, year, sex, nochild5, seas)

#fit model to individual trajectories & obtain predictions and 95%CI
model_lcd4 = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + sex + nochild5 + seas, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, cd4cnt == 0), na.action = na.exclude)
lcd4$fit = predict.gam(model_lcd4, type = "response", se.fit = TRUE)$fit
lcd4$fit_lci = model_lcd4$family$linkinv(predict.gam(model_lcd4, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lcd4, type = "link", se.fit = TRUE)$se.fit))
lcd4$fit_uci = model_lcd4$family$linkinv(predict.gam(model_lcd4, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lcd4, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
lcd4 <- lcd4 %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
lcd4_age <- left_join(left_join(lcd4 %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
lcd4 %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
lcd4 %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
lcd4_year <- left_join(left_join(lcd4 %>% group_by(year) %>% tally() %>% rename(Tot = n), 
lcd4 %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
lcd4 %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
E <- ggplot(lcd4_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"], color = "gray") +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"], color = "gray") +
  scale_y_continuous("VT carriage prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.8)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "CD4+ count <250", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

F <- ggplot(lcd4_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 1, color = brocolors("crayons")["Forest Green"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Forest Green"], color = "gray") +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Mahogany"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Mahogany"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.8)) + 
  labs(title = "", x = "Year") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a dataset to store model estimates
hcd4 = filter(pcvpa.mod, cd4cnt == 1) %>% select(vtcarr, age, year, sex, nochild5, seas)

#fit model to individual trajectories & obtain predictions and 95%CI
model_hcd4 = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + sex + nochild5 + seas, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, cd4cnt == 1), na.action = na.exclude)
hcd4$fit = predict.gam(model_hcd4, type = "response", se.fit = TRUE)$fit
hcd4$fit_lci = model_hcd4$family$linkinv(predict.gam(model_hcd4, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_hcd4, type = "link", se.fit = TRUE)$se.fit))
hcd4$fit_uci = model_hcd4$family$linkinv(predict.gam(model_hcd4, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_hcd4, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

#create age group
hcd4 <- hcd4 %>% 
  mutate(agegp = if_else(age <=24, "18-24",
                         if_else(age >24 & age <=29, "25-29",
                                 if_else(age >29 & age <=34, "30-34",
                                         if_else(age >34 & age <=40, "35-40", NA_character_)))))

#get age group mean predicted prevalence and acquisitions
hcd4_age <- left_join(left_join(hcd4 %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
hcd4 %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
hcd4 %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#join observed and predicted datasets for yearey number
hcd4_year <- left_join(left_join(hcd4 %>% group_by(year) %>% tally() %>% rename(Tot = n), 
hcd4 %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
hcd4 %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], foi = fit/42, foi_lci = fit_lci/42, foi_uci = fit_uci/42)

#----------------------------------------------------------------------------------

#plot prevalence curves
G <- ggplot(hcd4_age) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit, group = 1), size = 1, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = agegp, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"], color = "gray") +
  geom_line(aes(x = agegp, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = agegp, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.8)) + 
  scale_x_discrete(expand = c(0.04,0.04)) +
  labs(title = "CD4+ count ≥250", x = "Age group (years)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

H <- ggplot(hcd4_year) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit, group = 1), size = 1, color = brocolors("crayons")["Navy Blue"]) +
  geom_ribbon(aes(x = year, y = fit, group = 1, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = brocolors("crayons")["Navy Blue"], color = "gray") +
  geom_line(aes(x = year, y = foi/0.1, group = 1), lty = "dashed", size = 0.7, color = brocolors("crayons")["Yellow Orange"]) +
  geom_ribbon(aes(x = year, y = foi/0.1, group = 1, ymin = foi_lci/0.1, ymax = foi_uci/0.1), alpha = 0.2, fill = brocolors("crayons")["Yellow Orange"], color = "gray") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "Daily carriage acquisition"), limits = c(0, 0.8)) + 
  labs(title = "", x = "yearey group") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig4_VT_prev_acq_art_cd4.tiff"),
       plot = (A | B | C | D) / (E | F | G | H),
       width = 13, height = 6, unit="in", dpi = 200)
