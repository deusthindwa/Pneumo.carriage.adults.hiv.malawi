#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#======================================================================================

#subset for a correct dataset
lowses = pcvpa.mod %>% select(nvtcarr, agegp, surv, sescat, sex, nochild5) %>% filter(sescat == 0) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_lowses = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = lowses)
lowses$fit = predict.gam(model_lowses, type = "response", se.fit = TRUE)$fit
lowses$fit_lci = model_lowses$family$linkinv(predict.gam(model_lowses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lowses, type = "link", se.fit = TRUE)$se.fit))
lowses$fit_uci = model_lowses$family$linkinv(predict.gam(model_lowses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lowses, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
lowses1 <- left_join(left_join(lowses %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
lowses %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
lowses %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
lowses2 <- left_join(left_join(lowses %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
lowses %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
lowses %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:7], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[8:14])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_lowses = scam(nvtcarr ~ s(agegp, bs="mdcv") + s(surv, bs="mpi") + sex + nochild5, family = binomial(link = "cloglog"), data = lowses)
lowses$foi <- ((-derivative.scam(model_lowses, smooth.number = 1, deriv = 1)$d * model_lowses$fitted.values) + (1/42*model_lowses$fitted.values))/(1-model_lowses$fitted.values)

#join observed and predicted datasets for agegp
lowses3 <- left_join(left_join(lowses %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
lowses %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
lowses %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
lowses4 <- left_join(left_join(lowses %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
lowses %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
lowses %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
A <- ggplot(data = cbind(lowses1, lowses3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "red") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("NVT prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 1)) + 
  labs(title = "Low SES", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

B <- ggplot(data = cbind(lowses2, lowses4)) +
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
midses = pcvpa.mod %>% select(nvtcarr, agegp, surv, sescat, sex, nochild5) %>% filter(sescat == 1) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_midses = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = midses)
midses$fit = predict.gam(model_midses, type = "response", se.fit = TRUE)$fit
midses$fit_lci = model_midses$family$linkinv(predict.gam(model_midses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_midses, type = "link", se.fit = TRUE)$se.fit))
midses$fit_uci = model_midses$family$linkinv(predict.gam(model_midses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_midses, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
midses1 <- left_join(left_join(midses %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
midses %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
midses %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16]) 

#join observed and predicted datasets for survey number
midses2 <- left_join(left_join(midses %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
midses %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
midses %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:7], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[8:14])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_midses = scam(nvtcarr ~ s(agegp, bs="mpi") + s(surv, bs="mpd") + sex + nochild5, family = binomial(link = "cloglog"), data = midses)
midses$foi <- ((-derivative.scam(model_midses, smooth.number = 1, deriv = 1)$d * model_midses$fitted.values) + (1/42*model_midses$fitted.values))/(1-model_midses$fitted.values)

#join observed and predicted datasets for agegp
midses3 <- left_join(left_join(midses %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
midses %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
midses %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
midses4 <- left_join(left_join(midses %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
midses %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
midses %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
C <- ggplot(data = cbind(midses1, midses3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "blue") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "blue", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "blue") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "blue") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 1)) + 
  labs(title = "Middle/High SES", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

D <- ggplot(data = cbind(midses2, midses4)) +
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

#=======================================================================================

#subset for a correct dataset
lowses = pcvpa.mod %>% select(vtcarr, agegp, surv, sescat, sex, nochild5) %>% filter(sescat == 0) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_lowses = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = lowses)
lowses$fit = predict.gam(model_lowses, type = "response", se.fit = TRUE)$fit
lowses$fit_lci = model_lowses$family$linkinv(predict.gam(model_lowses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lowses, type = "link", se.fit = TRUE)$se.fit))
lowses$fit_uci = model_lowses$family$linkinv(predict.gam(model_lowses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lowses, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
lowses1 <- left_join(left_join(lowses %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
lowses %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
lowses %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
lowses2 <- left_join(left_join(lowses %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
lowses %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
lowses %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:7], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[8:14])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_lowses = scam(vtcarr ~ s(agegp, bs="mdcx") + s(surv, bs="mdcv") + sex + nochild5, family = binomial(link = "cloglog"), data = lowses)
lowses$foi <- ((-derivative.scam(model_lowses, smooth.number = 1, deriv = 1)$d * model_lowses$fitted.values) + (1/42*model_lowses$fitted.values))/(1-model_lowses$fitted.values)

#join observed and predicted datasets for agegp
lowses3 <- left_join(left_join(lowses %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
lowses %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
lowses %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
lowses4 <- left_join(left_join(lowses %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
lowses %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
lowses %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
E <- ggplot(data = cbind(lowses1, lowses3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "red") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("VT prevalence", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 1)) + 
  labs(title = "Low SES", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

F <- ggplot(data = cbind(lowses2, lowses4)) +
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
midses = pcvpa.mod %>% select(vtcarr, agegp, surv, sescat, sex, nochild5) %>% filter(sescat == 1) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_midses = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = midses)
midses$fit = predict.gam(model_midses, type = "response", se.fit = TRUE)$fit
midses$fit_lci = model_midses$family$linkinv(predict.gam(model_midses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_midses, type = "link", se.fit = TRUE)$se.fit))
midses$fit_uci = model_midses$family$linkinv(predict.gam(model_midses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_midses, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
midses1 <- left_join(left_join(midses %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
midses %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
midses %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
midses2 <- left_join(left_join(midses %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
midses %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
midses %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:7], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[8:14])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_midses = scam(vtcarr ~ s(agegp, bs="mdcx") + s(surv, bs="mdcv") + sex + nochild5, family = binomial(link = "cloglog"), data = midses)
midses$foi <- ((-derivative.scam(model_midses, smooth.number = 1, deriv = 1)$d * model_midses$fitted.values) + (1/42*model_midses$fitted.values))/(1-model_midses$fitted.values)

#join observed and predicted datasets for agegp
midses3 <- left_join(left_join(midses %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
midses %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
midses %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
midses4 <- left_join(left_join(midses %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
midses %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
midses %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
G <- ggplot(data = cbind(midses1, midses3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "blue") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "blue", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.1), size = 1, shape = 18, color = "blue") +
  geom_line(aes(x = agegp, y = foi/0.1), lty = "dashed", size = 0.7, color = "blue") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.1, name = ""), limits = c(0, 1)) + 
  labs(title = "Middle/High SES", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

H <- ggplot(data = cbind(midses2, midses4)) +
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

#=======================================================================================

#turn on warnings
options(warn = defaultW)

#with middle/high SES combine
ggsave(here("output", "Fig6_ses_prev_foi.tiff"),
       plot = (A | B | C | D) / (E | F | G | H),
       width = 14, height = 6, unit="in", dpi = 200)
