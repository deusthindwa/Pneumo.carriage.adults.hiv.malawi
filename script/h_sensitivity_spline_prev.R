#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1) 

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(nvtcarr, agegp, surv) %>% group_by(surv)

#fit different models & obtain predictions and 95%CI
model_crude = gam(nvtcarr ~ te(agegp, bs="cs") + te(surv, bs="cs"), family = binomial(link = "cloglog"), data = crude)
crude$fitcs = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(nvtcarr ~ te(agegp, bs="tp") + te(surv, bs="tp"), family = binomial(link = "cloglog"), data = crude)
crude$fittp = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(nvtcarr ~ te(agegp, bs="cr") + te(surv, bs="cr"), family = binomial(link = "cloglog"), data = crude)
crude$fitcr = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps"), family = binomial(link = "cloglog"), data = crude)
crude$fitps = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(fitcs = mean(fitcs), fittp = mean(fittp), fitcr = mean(fitcr), fitps = mean(fitps)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(fitcs = mean(fitcs), fittp = mean(fittp), fitcr = mean(fitcr), fitps = mean(fitps)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#plot prevalence curves
A <- ggplot(data = crude1) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fitcs, colour = "Natural cubic spline"), size = 0.8) +
  geom_line(aes(x = agegp, y = fittp, colour = "Thin plate regression spline"), size = 0.8) +
  geom_line(aes(x = agegp, y = fitps, colour = "P-spline"), size = 0.8) +
  geom_line(aes(x = agegp, y = fitcr, colour = "Cubic regression spline"), size = 0.8) +
  scale_color_discrete(breaks=c("Natural cubic spline", "Thin plate regression spline", "P-spline", "Cubic regression spline")) +
  guides(color=guide_legend(title="Types of splines fitted"), size = FALSE) +
  labs(title = "NVT(-ST3), Overall", x = "Age,y", y = "NVT prevalence") +
  theme_bw() +
  ylim(0, 0.4) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

B <- ggplot(data = crude2) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fitcs, color = "Natural cubic spline"), size = 0.8) +
  geom_line(aes(x = surv, y = fittp, color = "Thin plate regression spline"), size = 0.8) +
  geom_line(aes(x = surv, y = fitps, color = "P-spline"), size = 0.8) +
  geom_line(aes(x = surv, y = fitcr, color = "Cubic regression spline"), size = 0.8) +
  scale_color_discrete(breaks=c("Natural cubic spline", "Thin plate regression spline", "P-spline", "Cubic regression spline")) +
  guides(color=guide_legend(title="Types of splines fitted"), size = FALSE) +
  labs(title = "", x = "Survey number", y = "") +
  theme_bw() +
  ylim(0, 0.4) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit different models & obtain predictions and 95%CI
model_crude = gam(vtcarr ~ te(agegp, bs="cs") + te(surv, bs="cs"), family = binomial(link = "cloglog"), data = crude)
crude$fitcs = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(vtcarr ~ te(agegp, bs="tp") + te(surv, bs="tp"), family = binomial(link = "cloglog"), data = crude)
crude$fittp = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(vtcarr ~ te(agegp, bs="cr") + te(surv, bs="cr"), family = binomial(link = "cloglog"), data = crude)
crude$fitcr = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps"), family = binomial(link = "cloglog"), data = crude)
crude$fitps = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fitcs = mean(fitcs), fittp = mean(fittp), fitcr = mean(fitcr), fitps = mean(fitps)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(fitcs = mean(fitcs), fittp = mean(fittp), fitcr = mean(fitcr), fitps = mean(fitps)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#plot prevalence curves
C <- ggplot(data = crude1) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fitcs, colour = "Natural cubic spline"), size = 0.8) +
  geom_line(aes(x = agegp, y = fittp, colour = "Thin plate regression spline"), size = 0.8) +
  geom_line(aes(x = agegp, y = fitps, colour = "P-spline"), size = 0.8) +
  geom_line(aes(x = agegp, y = fitcr, colour = "Cubic regression spline"), size = 0.8) +
  scale_color_discrete(breaks=c("Natural cubic spline", "Thin plate regression spline", "P-spline", "Cubic regression spline")) +
  guides(color=guide_legend(title="Types of splines fitted"), size = FALSE) +
  labs(title = "VT(+ST3), Overall", x = "Age,y", y = "VT prevalence") +
  theme_bw() +
  ylim(0, 0.4) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

D <- ggplot(data = crude2) +
  geom_point(aes(x = surv, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fitcs, color = "Natural (cubic) spline"), size = 0.8) +
  geom_line(aes(x = surv, y = fittp, color = "Thin plate regression spline"), size = 0.8) +
  geom_line(aes(x = surv, y = fitps, color = "P-spline (B-spline plus penalisation)                   "), size = 0.8) +
  geom_line(aes(x = surv, y = fitcr, color = "Cubic regression spline"), size = 0.8) +
  scale_color_discrete(breaks=c("Natural (cubic) spline", "Thin plate regression spline", "P-spline (B-spline plus penalisation)                   ", "Cubic regression spline")) +
  guides(color=guide_legend(title="Types of splines fitted"), size = FALSE) +
  labs(title = "", x = "Survey number", y = "") +
  theme_bw() +
  ylim(0, 0.4) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "right")

#=======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(nvtcarr, agegp, surv) %>% group_by(surv)

#fit different model formulations & obtain predictions and 95%CI
model1_crude = gam(nvtcarr ~ te(agegp, bs="ps"), family = binomial(link = "cloglog"), data = crude)
crude$fit1 = predict.gam(model1_crude, type = "response", se.fit = TRUE)$fit

model2_crude = gam(nvtcarr ~ te(surv, bs="ps"), family = binomial(link = "cloglog"), data = crude)
crude$fit2 = predict.gam(model2_crude, type = "response", se.fit = TRUE)$fit

model3_crude = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="cs"), family = binomial(link = "cloglog"), data = crude)
crude$fit3 = predict.gam(model3_crude, type = "response", se.fit = TRUE)$fit

model4_crude = gam(nvtcarr ~ te(agegp, bs="ps") + te(surv, bs="cs") + ti(agegp, surv, bs="ps"), family = binomial(link = "cloglog"), data = crude)
crude$fit4 = predict.gam(model4_crude, type = "response", se.fit = TRUE)$fit

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(fit1 = mean(fit1), fit2 = mean(fit2), fit3 = mean(fit3), fit4 = mean(fit4)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(fit1 = mean(fit1), fit2 = mean(fit2), fit3 = mean(fit3), fit4 = mean(fit4)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#plot prevalence curves
E <- ggplot(data = crude1) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit1, colour = "Model1"), size = 0.8) +
  geom_line(aes(x = agegp, y = fit2, colour = "Model2"), size = 0.8) +
  geom_line(aes(x = agegp, y = fit3, colour = "Model3"), size = 0.8) +
  geom_line(aes(x = agegp, y = fit4, colour = "Model4"), size = 0.8) +
  scale_color_discrete(breaks=c("Model1", "Model2", "Model3", "Model4")) +
  guides(color=guide_legend(title="Model formulation"), size = FALSE) +
  labs(title = "NVT(-ST3), Overall", x = "Age,y", y = "NVT prevalence") +
  theme_bw() +
  ylim(0, 0.4) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

F <- ggplot(data = crude2) +
  geom_point(aes(x = surv, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit1, colour = "Model1"), size = 0.8) +
  geom_line(aes(x = surv, y = fit2, colour = "Model2"), size = 0.8) +
  geom_line(aes(x = surv, y = fit3, colour = "Model3"), size = 0.8) +
  geom_line(aes(x = surv, y = fit4, colour = "Model4"), size = 0.8) +
  scale_color_discrete(breaks=c("Model1", "Model2", "Model3", "Model4")) +
  guides(color=guide_legend(title="Model formulation"), size = FALSE) +
  labs(title = "", x = "Survey number", y = "") +
  theme_bw() +
  ylim(0, 0.4) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#=======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit different model formulations & obtain predictions and 95%CI
model1_crude = gam(vtcarr ~ te(agegp, bs="ps"), family = binomial(link = "cloglog"), data = crude)
crude$fit1 = predict.gam(model1_crude, type = "response", se.fit = TRUE)$fit

model2_crude = gam(vtcarr ~ te(surv, bs="ps"), family = binomial(link = "cloglog"), data = crude)
crude$fit2 = predict.gam(model2_crude, type = "response", se.fit = TRUE)$fit

model3_crude = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="cs"), family = binomial(link = "cloglog"), data = crude)
crude$fit3 = predict.gam(model3_crude, type = "response", se.fit = TRUE)$fit

model4_crude = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="cs") + ti(agegp, surv, bs="ps"), family = binomial(link = "cloglog"), data = crude)
crude$fit4 = predict.gam(model4_crude, type = "response", se.fit = TRUE)$fit

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit1 = mean(fit1), fit2 = mean(fit2), fit3 = mean(fit3), fit4 = mean(fit4)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(fit1 = mean(fit1), fit2 = mean(fit2), fit3 = mean(fit3), fit4 = mean(fit4)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:8], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[9:16])

#plot prevalence curves
G <- ggplot(data = crude1) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit1, colour = "Model1"), size = 0.8) +
  geom_line(aes(x = agegp, y = fit2, colour = "Model2"), size = 0.8) +
  geom_line(aes(x = agegp, y = fit3, colour = "Model3"), size = 0.8) +
  geom_line(aes(x = agegp, y = fit4, colour = "Model4"), size = 0.8) +
  scale_color_discrete(breaks=c("Model1", "Model2", "Model3", "Model4")) +
  guides(color=guide_legend(title="Model formulation"), size = FALSE) +
  labs(title = "VT(+ST3), Overall", x = "Age,y", y = "VT prevalence") +
  theme_bw() +
  ylim(0, 0.4) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

H <- ggplot(data = crude2) +
  geom_point(aes(x = surv, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(surv, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = surv, y = fit1, colour = "Model1: constant + te(age)"), size = 0.8) +
  geom_line(aes(x = surv, y = fit2, colour = "Model2: constant + te(time)"), size = 0.8) +
  geom_line(aes(x = surv, y = fit3, colour = "Model3: constant + te(age) + te(time)"), size = 0.8) +
  geom_line(aes(x = surv, y = fit4, colour = "Model4: constant + te(age) + te(time) + ti(age, time)"), size = 0.8) +
  scale_color_discrete(breaks=c("Model1: constant + te(age)", "Model2: constant + te(time)", "Model3: constant + te(age) + te(time)", "Model4: constant + te(age) + te(time) + ti(age, time)")) +
  guides(color=guide_legend(title="Model formulation"), size=FALSE) +
  labs(title = "", x = "Survey number", y = "") +
  theme_bw() +
  ylim(0, 0.4) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "right")

#=======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "SFig1_sens_spline_model.tiff"),
       plot = (A | B | C | D) / (E | F | G | H),
       width = 14, height = 6, unit="in", dpi = 200)
