#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1) 

#======================================================================================

#subset for a correct dataset
shorter = pcvpa.mod %>% select(nvtcarr, agegp, surv, artdur) %>% filter(artdur == 0) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_shorter = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps"), family = binomial(link = "cloglog"), data = shorter)
shorter$fit = predict.gam(model_shorter, type = "response", se.fit = TRUE)$fit
shorter$fit_lci = model_shorter$family$linkinv(predict.gam(model_shorter, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_shorter, type = "link", se.fit = TRUE)$se.fit))
shorter$fit_uci = model_shorter$family$linkinv(predict.gam(model_shorter, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_shorter, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
shorter1 <- left_join(left_join(shorter %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
shorter %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
shorter %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
shorter2 <- left_join(left_join(shorter %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
shorter %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
shorter %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:7], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[8:14])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_shorter = scam(nvtcarr ~ s(agegp, bs="mpi") + s(surv, bs="mpi"), family = binomial(link = "cloglog"), data = shorter)
shorter$foi <- ((-derivative.scam(model_shorter, smooth.number = 1, deriv = 1)$d * model_shorter$fitted.values) + (1/42*model_shorter$fitted.values))/(1-model_shorter$fitted.values)

#join observed and predicted datasets for agegp
shorter3 <- left_join(left_join(shorter %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
shorter %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
shorter %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#fit model to obtain predictions of FOI
model_shorter = scam(nvtcarr ~ s(agegp, bs="mpd") + s(surv, bs="mpd"), family = binomial(link = "cloglog"), data = shorter)
shorter$foi <- -derivative.scam(model_shorter, deriv = 1)$d * model_shorter$fitted.values

#join observed and predicted datasets for survey number
shorter4 <- left_join(left_join(shorter %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
shorter %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
shorter %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
A <- ggplot(data = cbind(shorter1, shorter3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "red") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.15), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = agegp, y = foi/0.15), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("NVT prevalence", sec.axis = sec_axis(~. * 0.15, name = ""), limits = c(0, 0.75)) + 
  labs(title = "ART <5y", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

B <- ggplot(data = cbind(shorter2, shorter4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "red") +
  geom_point(aes(x = surv, y = foi/0.15), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = surv, y = foi/0.15), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.15, name = ""), limits = c(0, 0.75)) + 
  geom_ribbon(aes(x = surv, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  labs(title = "", x = "Survey number") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#=======================================================================================

#subset for a correct dataset
longer = pcvpa.mod %>% select(nvtcarr, agegp, surv, artdur) %>% filter(artdur == 1) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_longer = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps"), family = binomial(link = "cloglog"), data = longer)
longer$fit = predict.gam(model_longer, type = "response", se.fit = TRUE)$fit
longer$fit_lci = model_longer$family$linkinv(predict.gam(model_longer, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_longer, type = "link", se.fit = TRUE)$se.fit))
longer$fit_uci = model_longer$family$linkinv(predict.gam(model_longer, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_longer, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
longer1 <- left_join(left_join(longer %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
longer %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
longer %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
longer2 <- left_join(left_join(longer %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
longer %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
longer %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:7], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[8:14])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_longer = scam(nvtcarr ~ s(agegp, bs="mdcv") + s(surv, bs="mdcv"), family = binomial(link = "cloglog"), data = longer)
longer$foi <- ((-derivative.scam(model_longer, smooth.number = 1, deriv = 1)$d * model_longer$fitted.values) + (1/42*model_longer$fitted.values))/(1-model_longer$fitted.values)

#join observed and predicted datasets for agegp
longer3 <- left_join(left_join(longer %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
longer %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
longer %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
longer4 <- left_join(left_join(longer %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
longer %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
longer %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
C <- ggplot(data = cbind(longer1, longer3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "darkgreen") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "green", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.15), size = 1, shape = 18, color = "darkgreen") +
  geom_line(aes(x = agegp, y = foi/0.15), lty = "dashed", size = 0.7, color = "darkgreen") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.15, name = ""), limits = c(0, 0.75)) + 
  labs(title = "ART 5+y", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

D <- ggplot(data = cbind(longer2, longer4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "darkgreen") +
  geom_point(aes(x = surv, y = foi/0.15), size = 1, shape = 18, color = "darkgreen") +
  geom_line(aes(x = surv, y = foi/0.15), lty = "dashed", size = 0.7, color = "darkgreen") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.15, name = "NVT force of infection"), limits = c(0, 0.75)) + 
  geom_ribbon(aes(x = surv, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "green", color = "gray") +
  labs(title = "", x = "Survey number") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#=======================================================================================

#subset for a correct dataset
shorter = pcvpa.mod %>% select(vtcarr, agegp, surv, artdur) %>% filter(artdur == 0) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_shorter = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps"), family = binomial(link = "cloglog"), data = shorter)
shorter$fit = predict.gam(model_shorter, type = "response", se.fit = TRUE)$fit
shorter$fit_lci = model_shorter$family$linkinv(predict.gam(model_shorter, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_shorter, type = "link", se.fit = TRUE)$se.fit))
shorter$fit_uci = model_shorter$family$linkinv(predict.gam(model_shorter, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_shorter, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
shorter1 <- left_join(left_join(shorter %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
shorter %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
shorter %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
shorter2 <- left_join(left_join(shorter %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
shorter %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
shorter %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:7], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[8:14])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_shorter = scam(vtcarr ~ s(agegp, bs="mpd") + s(surv, bs="mpd"), family = binomial(link = "cloglog"), data = shorter)
shorter$foi <- ((-derivative.scam(model_shorter, smooth.number = 1, deriv = 1)$d * model_shorter$fitted.values) + (1/42*model_shorter$fitted.values))/(1-model_shorter$fitted.values)

#join observed and predicted datasets for agegp
shorter3 <- left_join(left_join(shorter %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
shorter %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
shorter %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
shorter4 <- left_join(left_join(shorter %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
shorter %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
shorter %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
E <- ggplot(data = cbind(shorter1, shorter3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "red") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.15), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = agegp, y = foi/0.15), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("VT prevalence", sec.axis = sec_axis(~. * 0.15, name = ""), limits = c(0, 0.75)) + 
  labs(title = "ART <5y", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

F <- ggplot(data = cbind(shorter2, shorter4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "red") +
  geom_point(aes(x = surv, y = foi/0.15), size = 1, shape = 18, color = "red") +
  geom_line(aes(x = surv, y = foi/0.15), lty = "dashed", size = 0.7, color = "red") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.15, name = ""), limits = c(0, 0.75)) + 
  geom_ribbon(aes(x = surv, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "red", color = "gray") +
  labs(title = "", x = "Survey number") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#=======================================================================================

#subset for a correct dataset
longer = pcvpa.mod %>% select(vtcarr, agegp, surv, artdur) %>% filter(artdur == 1) %>% group_by(surv)

#fit model & obtain predictions and 95%CI
model_longer = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps"), family = binomial(link = "cloglog"), data = longer)
longer$fit = predict.gam(model_longer, type = "response", se.fit = TRUE)$fit
longer$fit_lci = model_longer$family$linkinv(predict.gam(model_longer, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_longer, type = "link", se.fit = TRUE)$se.fit))
longer$fit_uci = model_longer$family$linkinv(predict.gam(model_longer, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_longer, type = "link", se.fit = TRUE)$se.fit))

#join observed and predicted datasets for agegp
longer1 <- left_join(left_join(longer %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
longer %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
longer %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
longer2 <- left_join(left_join(longer %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
longer %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
longer %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>% mutate(sex="Female") %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:7], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[8:14])

#----------------------------------------------------------------------------------

#fit model to obtain predictions of FOI
model_longer = scam(vtcarr ~ s(agegp, bs="mpd") + s(surv, bs="mpd"), family = binomial(link = "cloglog"), data = longer)
longer$foi <- ((-derivative.scam(model_longer, smooth.number = 1, deriv = 1)$d * model_longer$fitted.values) + (1/42*model_longer$fitted.values))/(1-model_longer$fitted.values)

#join observed and predicted datasets for agegp
longer3 <- left_join(left_join(longer %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
longer %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
longer %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#join observed and predicted datasets for survey number
longer4 <- left_join(left_join(longer %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
longer %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
longer %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(foi = mean(foi)) %>%
ungroup()) %>% select(foi)

#----------------------------------------------------------------------------------

#plot prevalence curves
G <- ggplot(data = cbind(longer1, longer3)) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit), size = 1, color = "darkgreen") +
  geom_ribbon(aes(x = agegp, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "green", color = "gray") +
  geom_point(aes(x = agegp, y = foi/0.15), size = 1, shape = 18, color = "darkgreen") +
  geom_line(aes(x = agegp, y = foi/0.15), lty = "dashed", size = 0.7, color = "darkgreen") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.15, name = ""), limits = c(0, 0.75)) + 
  labs(title = "ART 7+y", x = "Age,y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

H <- ggplot(data = cbind(longer2, longer4)) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit), size = 1, color = "darkgreen") +
  geom_point(aes(x = surv, y = foi/0.15), size = 1, shape = 18, color = "darkgreen") +
  geom_line(aes(x = surv, y = foi/0.15), lty = "dashed", size = 0.7, color = "darkgreen") +
  scale_y_continuous("", sec.axis = sec_axis(~. * 0.15, name = "VT force of infection"), limits = c(0, 0.75)) +
  geom_ribbon(aes(x = surv, y = fit, ymin = fit_lci, ymax = fit_uci), alpha = 0.2, fill = "green", color = "gray") +
  labs(title = "", x = "Survey number") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#=======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig4_artdur_prev_foi2.tiff"),
       plot = (A | B | C | D) / (E | F | G | H ),
       width = 14, height = 6, unit="in", dpi = 200)
