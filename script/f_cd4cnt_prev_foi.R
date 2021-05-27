#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1) 

#======================================================================================

#subset for a correct dataset
lowcd4cnt = pcvpa.mod %>% select(nvtcarr, agegp, surv, cd4cnt, sex, nochild5) %>% filter(cd4cnt == 0) %>% group_by(surv)
lowcd4cnt = lowcd4cnt %>% mutate(surv = if_else(surv == 8L | surv == 7L, 6L, surv))

#fit model & obtain predictions and 95%CI
model_lowcd4cnt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = lowcd4cnt)
lowcd4cnt$fit = predict.gam(model_lowcd4cnt, type = "response", se.fit = TRUE)$fit
lowcd4cnt$fit_lci = model_lowcd4cnt$family$linkinv(predict.gam(model_lowcd4cnt, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lowcd4cnt, type = "link", se.fit = TRUE)$se.fit))
lowcd4cnt$fit_uci = model_lowcd4cnt$family$linkinv(predict.gam(model_lowcd4cnt, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lowcd4cnt, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
lowcd4cnt1 <- left_join(left_join(lowcd4cnt %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
lowcd4cnt %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
lowcd4cnt %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(cd4cnt="Low CD4+ count") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
lowcd4cnt2 <- left_join(left_join(lowcd4cnt %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
lowcd4cnt %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
lowcd4cnt %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(cd4cnt="Low CD4+ count") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:6], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[7:12])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_lowcd4cnt = scam(nvtcarr ~ s(agegp, bs="micx") + s(surv, bs="micx") + sex + nochild5, family = binomial(link = "cloglog"), data = lowcd4cnt)
lowcd4cnt$foi <- ((derivative.scam(model_lowcd4cnt, deriv = 1)$d * model_lowcd4cnt$fitted.values) + (1/42*model_lowcd4cnt$fitted.values))/(1-model_lowcd4cnt$fitted.values)

#join observed and predicted datasets for agegp
lowcd4cnt3 <- left_join(left_join(lowcd4cnt %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
lowcd4cnt %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
lowcd4cnt %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
lowcd4cnt4 <- left_join(left_join(lowcd4cnt %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
lowcd4cnt %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
lowcd4cnt %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
A <- ggplot(data = cbind(lowcd4cnt1, lowcd4cnt3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "red") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("NVT prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 1)) + 
  labs(title = "CD4+ count <250", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

B <- ggplot(data = cbind(lowcd4cnt2, lowcd4cnt4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "red") +
  geom_point(aes(x = surv, y = foi/0.1), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = surv, y = foi/0.1), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 1)) + 
  geom_ribbon(aes(x = surv, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  labs(title = "", x = "Survey number") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#=======================================================================================

#subset for a correct dataset
highcd4cnt = pcvpa.mod %>% select(nvtcarr, agegp, surv, cd4cnt, sex, nochild5) %>% filter(cd4cnt == 1) %>% group_by(surv)
highcd4cnt = highcd4cnt %>% mutate(surv = if_else(surv == 8L | surv == 7L, 6L, surv))

#fit model & obtain predictions and 95%CI
model_highcd4cnt = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = highcd4cnt)
highcd4cnt$fit = predict.gam(model_highcd4cnt, type = "response", se.fit = TRUE)$fit
highcd4cnt$fit_lci = model_highcd4cnt$family$linkinv(predict.gam(model_highcd4cnt, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_highcd4cnt, type = "link", se.fit = TRUE)$se.fit))
highcd4cnt$fit_uci = model_highcd4cnt$family$linkinv(predict.gam(model_highcd4cnt, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_highcd4cnt, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
highcd4cnt1 <- left_join(left_join(highcd4cnt %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
highcd4cnt %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
highcd4cnt %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(cd4cnt="High CD4+ count") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
highcd4cnt2 <- left_join(left_join(highcd4cnt %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
highcd4cnt %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
highcd4cnt %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(cd4cnt="High CD4+ count") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:6], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[7:12])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_highcd4cnt = scam(nvtcarr ~ s(agegp, bs="mpd") + s(surv, bs="mpd") + sex + nochild5, family = binomial(link = "cloglog"), data = highcd4cnt)
highcd4cnt$foi <- ((-derivative.scam(model_highcd4cnt, deriv = 1)$d * model_highcd4cnt$fitted.values) + (1/42*model_highcd4cnt$fitted.values))/(1-model_highcd4cnt$fitted.values)

#join observed and predicted datasets for agegp
highcd4cnt3 <- left_join(left_join(highcd4cnt %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
highcd4cnt %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
highcd4cnt %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
highcd4cnt4 <- left_join(left_join(highcd4cnt %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
highcd4cnt %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
highcd4cnt %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
C <- ggplot(data = cbind(highcd4cnt1, highcd4cnt3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "blue") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "blue", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "blue") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "blue") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 1)) + 
  labs(title = "CD4+ count ≥250", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

D <- ggplot(data = cbind(highcd4cnt2, highcd4cnt4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "blue") +
  geom_point(aes(x = surv, y = foi/0.1), size = 1, shape = 18, color = "blue") +
  geom_line(aes(x = surv, y = foi/0.1), lty = "dashed", size = 0.7, color = "blue") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "NVT force of infection"), limits = c(0, 1)) + 
  geom_ribbon(aes(x = surv, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "blue", color = "gray") +
  labs(title = "", x = "Survey number") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a correct dataset
lowcd4cnt = pcvpa.mod %>% select(vtcarr, agegp, surv, cd4cnt, sex, nochild5) %>% filter(cd4cnt == 0) %>% group_by(surv)
lowcd4cnt = lowcd4cnt %>% mutate(surv = if_else(surv == 8L | surv == 7L, 6L, surv))

#fit model & obtain predictions and 95%CI
model_lowcd4cnt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = lowcd4cnt)
lowcd4cnt$fit = predict.gam(model_lowcd4cnt, type = "response", se.fit = TRUE)$fit
lowcd4cnt$fit_lci = model_lowcd4cnt$family$linkinv(predict.gam(model_lowcd4cnt, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lowcd4cnt, type = "link", se.fit = TRUE)$se.fit))
lowcd4cnt$fit_uci = model_lowcd4cnt$family$linkinv(predict.gam(model_lowcd4cnt, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lowcd4cnt, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
lowcd4cnt1 <- left_join(left_join(lowcd4cnt %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
lowcd4cnt %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
lowcd4cnt %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(cd4cnt="Low CD4+ count") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
lowcd4cnt2 <- left_join(left_join(lowcd4cnt %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
lowcd4cnt %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
lowcd4cnt %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(cd4cnt="Low CD4+ count") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:6], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[7:12])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_lowcd4cnt = scam(vtcarr ~ s(agegp, bs="mdcx") + s(surv, bs="mpd") + sex + nochild5, family = binomial(link = "cloglog"), data = lowcd4cnt)
lowcd4cnt$foi <- ((-derivative.scam(model_lowcd4cnt, deriv = 1)$d * model_lowcd4cnt$fitted.values) + (1/42*model_lowcd4cnt$fitted.values))/(1-model_lowcd4cnt$fitted.values)

#join observed and predicted datasets for agegp
lowcd4cnt3 <- left_join(left_join(lowcd4cnt %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
lowcd4cnt %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
lowcd4cnt %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
lowcd4cnt4 <- left_join(left_join(lowcd4cnt %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
lowcd4cnt %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
lowcd4cnt %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
E <- ggplot(data = cbind(lowcd4cnt1, lowcd4cnt3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "red") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("VT prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 1)) + 
  labs(title = "CD4+ count <250", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

F <- ggplot(data = cbind(lowcd4cnt2, lowcd4cnt4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "red") +
  geom_point(aes(x = surv, y = foi/0.1), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = surv, y = foi/0.1), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 1)) + 
  geom_ribbon(aes(x = surv, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  labs(title = "", x = "Survey number") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#=======================================================================================

#subset for a correct dataset
highcd4cnt = pcvpa.mod %>% select(vtcarr, agegp, surv, cd4cnt, sex, nochild5) %>% filter(cd4cnt == 1) %>% group_by(surv)
highcd4cnt = highcd4cnt %>% mutate(surv = if_else(surv == 8L | surv == 7L, 6L, surv))

#fit model & obtain predictions and 95%CI
model_highcd4cnt = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = highcd4cnt)
highcd4cnt$fit = predict.gam(model_highcd4cnt, type = "response", se.fit = TRUE)$fit
highcd4cnt$fit_lci = model_highcd4cnt$family$linkinv(predict.gam(model_highcd4cnt, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_highcd4cnt, type = "link", se.fit = TRUE)$se.fit))
highcd4cnt$fit_uci = model_highcd4cnt$family$linkinv(predict.gam(model_highcd4cnt, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_highcd4cnt, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
highcd4cnt1 <- left_join(left_join(highcd4cnt %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
highcd4cnt %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
highcd4cnt %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(cd4cnt="High CD4+ count") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
highcd4cnt2 <- left_join(left_join(highcd4cnt %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
highcd4cnt %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
highcd4cnt %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(cd4cnt="High CD4+ count") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:6], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[7:12])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_highcd4cnt = scam(vtcarr ~ s(agegp, bs="mpd") + s(surv, bs="mpd") + sex + nochild5, family = binomial(link = "cloglog"), data = highcd4cnt)
highcd4cnt$foi <- ((-derivative.scam(model_highcd4cnt, deriv = 1)$d * model_highcd4cnt$fitted.values) + (1/42*model_highcd4cnt$fitted.values))/(1-model_highcd4cnt$fitted.values)

#join observed and predicted datasets for agegp
highcd4cnt3 <- left_join(left_join(highcd4cnt %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
highcd4cnt %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
highcd4cnt %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
highcd4cnt4 <- left_join(left_join(highcd4cnt %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
highcd4cnt %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
highcd4cnt %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
G <- ggplot(data = cbind(highcd4cnt1, highcd4cnt3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "blue") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "blue", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "blue") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "blue") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 1)) + 
  labs(title = "CD4+ count ≥250", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

H <- ggplot(data = cbind(highcd4cnt2, highcd4cnt4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "blue") +
  geom_point(aes(x = surv, y = foi/0.1), size = 1, shape = 18, color = "blue") +
  geom_line(aes(x = surv, y = foi/0.1), lty = "dashed", size = 0.7, color = "blue") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = "VT force of infection"), limits = c(0, 1)) + 
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

ggsave(here("output", "Fig5_cd4+_prev_foi.tiff"),
       plot = (A | B | C | D) / (E | F | G | H),
       width = 14, height = 6, unit="in", dpi = 200)
