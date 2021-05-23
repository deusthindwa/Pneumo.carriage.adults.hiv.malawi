#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(nvtcarr, agegp, surv, sex, nochild5) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_crude = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_crude = scam(nvtcarr ~ s(agegp, bs="mdcx") + s(surv, bs="mdcx") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp
crude3 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
crude4 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
A <- ggplot(data = cbind(crude1, crude3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "red") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("NVT prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.6)) + 
  labs(title = "NVT(-ST3), Overall", x = "Age,y") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

B <- ggplot(data = cbind(crude2, crude4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "red") +
  geom_point(aes(x = surv, y = foi/0.1), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = surv, y = foi/0.1), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "NVT force of infection"), limits = c(0, 0.6)) + 
  geom_ribbon(aes(x = surv, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  labs(title = "", x = "Survey number") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr, agegp, surv, sex, nochild5) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_crude = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_crude = scam(vtcarr ~ s(agegp, bs="mdcv") + s(surv, bs="mdcx") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp
crude3 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
crude4 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
C <- ggplot(data = cbind(crude1, crude3)) +
  geom_point(aes(x = agegp, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(x = agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "blue") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "blue") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "blue") +
  scale_y_continuous("VT prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.6)) + 
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "blue", color = "gray") +
  labs(title = "VT(+ST3), Overall", x = "Age,y") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

D <- ggplot(data = cbind(crude2, crude4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "blue") +
  geom_point(aes(x = surv, y = foi/0.1), size = 1, shape = 18, color = "blue") +
  geom_line(aes(x = surv, y = foi/0.1), lty = "dashed", size = 0.7, color = "blue") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "VT force of infection"), limits = c(0, 0.6)) + 
  geom_ribbon(aes(x = surv, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "blue", color = "gray") +
  labs(title = "", x = "Survey number") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#=======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(nvtcarr1, agegp, surv, sex, nochild5) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_crude = gam(nvtcarr1 ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr1 == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr1 != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr1 == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr1 != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_crude = scam(nvtcarr1 ~ s(agegp, bs="mdcv") + s(surv, bs="mpd") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp
crude3 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr1 == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr1 != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
crude4 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr1 == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr1 != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
E <- ggplot(data = cbind(crude1, crude3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "red") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("NVT prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.6)) + 
  labs(title = "NVT(+ST3), Overall", x = "Age,y") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

F <- ggplot(data = cbind(crude2, crude4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "red") +
  geom_point(aes(x = surv, y = foi/0.1), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = surv, y = foi/0.1), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "NVT force of infection"), limits = c(0, 0.6)) + 
  geom_ribbon(aes(x = surv, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  labs(title = "", x = "Survey number") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr1, agegp, surv, sex, nochild5) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_crude = gam(vtcarr1 ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr1 == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr1 != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr1 == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr1 != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_crude = scam(vtcarr1 ~ s(agegp, bs="mdcx") + s(surv, bs="mdcx") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp
crude3 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr1 == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr1 != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
crude4 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr1 == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr1 != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
G <- ggplot(data = cbind(crude1, crude3)) +
  geom_point(aes(x = agegp, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(x = agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "blue") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "blue") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "blue") +
  scale_y_continuous("VT prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 0.6)) + 
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "blue", color = "gray") +
  labs(title = "VT(-ST3), Overall", x = "Age,y") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

H <- ggplot(data = cbind(crude2, crude4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "blue") +
  geom_point(aes(x = surv, y = foi/0.1), size = 1, shape = 18, color = "blue") +
  geom_line(aes(x = surv, y = foi/0.1), lty = "dashed", size = 0.7, color = "blue") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "VT force of infection"), limits = c(0, 0.6)) + 
  geom_ribbon(aes(x = surv, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "blue", color = "gray") +
  labs(title = "", x = "Survey number") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig2_crude_prev_foi.tiff"),
       plot = (A | B | C | D) / (E | F | G | H),
       width = 14, height = 6, unit="in", dpi = 200)
