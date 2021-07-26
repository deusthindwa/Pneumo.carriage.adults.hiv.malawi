#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1) 

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr, age, year)

#fit different models & obtain predictions and 95%CI
model_crude = gam(vtcarr ~ te(age, bs="cs") + te(year, bs="cs") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fitcs = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(vtcarr ~ te(age, bs="tp") + te(year, bs="tp") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fittp = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(vtcarr ~ te(age, bs="cr") + te(year, bs="cr") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fitcr = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fitps = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24",
                                   if_else(age >24 & age <=29, "25-29",
                                           if_else(age >29 & age <=34, "30-34",
                                                   if_else(age >34 & age <=40, "35-40", NA_character_))))))

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fitcs = mean(fitcs), fittp = mean(fittp), fitcr = mean(fitcr), fitps = mean(fitps)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8])

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fitcs = mean(fitcs), fittp = mean(fittp), fitcr = mean(fitcr), fitps = mean(fitps)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10])

#plot prevalence curves
A <- ggplot(data = crude1) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fitcs, colour = "Natural cubic spline", group = 1), size = 0.9) +
  geom_line(aes(x = agegp, y = fittp, colour = "Thin plate regression spline", group = 1), size = 0.9) +
  geom_line(aes(x = agegp, y = fitps, colour = "P-spline", group = 1), size = 0.9) +
  geom_line(aes(x = agegp, y = fitcr, colour = "Cubic regression spline", group = 1), size = 0.9) +
  scale_color_discrete(breaks=c("Natural cubic spline", "Thin plate regression spline", "P-spline", "Cubic regression spline")) +
  guides(color=guide_legend(title="Types of fitted splines"), size = FALSE) +
  labs(title = "Overall VT carriage", x = "Age in years", y = "Prevalence") +
  theme_bw() +
  ylim(0, 0.4) +
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold")) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) +
  theme(legend.position = "none")

B <- ggplot(data = crude2) +
  geom_point(aes(x = year, y = Pos/Tot, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fitcs, color = "Natural cubic spline", group = 1), size = 0.9) +
  geom_line(aes(x = year, y = fittp, color = "Thin plate regression spline", group = 1), size = 0.9) +
  geom_line(aes(x = year, y = fitps, color = "P-spline", group = 1), size = 0.9) +
  geom_line(aes(x = year, y = fitcr, color = "Cubic regression spline", group = 1), size = 0.9) +
  scale_color_discrete(breaks=c("Natural cubic spline", "Thin plate regression spline", "P-spline", "Cubic regression spline")) +
  guides(color=guide_legend(title=""), size = FALSE) +
  labs(title = "", x = "Survey year", y = "") +
  theme_bw() +
  ylim(0, 0.4) +
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold")) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) +
  theme(legend.position = c(0.5, 0.72))

#=======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr, age, year)

#fit different model formulations & obtain predictions and 95%CI
model1_crude = gam(vtcarr ~ te(age, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fit1 = predict.gam(model1_crude, type = "response", se.fit = TRUE)$fit

model2_crude = gam(vtcarr ~ te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fit2 = predict.gam(model2_crude, type = "response", se.fit = TRUE)$fit

model3_crude = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fit3 = predict.gam(model3_crude, type = "response", se.fit = TRUE)$fit

model4_crude = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + ti(age, year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fit4 = predict.gam(model4_crude, type = "response", se.fit = TRUE)$fit

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24",
                                   if_else(age >24 & age <=29, "25-29",
                                           if_else(age >29 & age <=34, "30-34",
                                                   if_else(age >34 & age <=40, "35-40", NA_character_))))))

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit1 = mean(fit1), fit2 = mean(fit2), fit3 = mean(fit3), fit4 = mean(fit4)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8])

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit1 = mean(fit1), fit2 = mean(fit2), fit3 = mean(fit3), fit4 = mean(fit4)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10])

#plot prevalence curves
C <- ggplot(data = crude1) +
  geom_point(aes(x = agegp, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(agegp, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = agegp, y = fit1, colour = "Model1", group = 1), size = 0.9) +
  geom_line(aes(x = agegp, y = fit2, colour = "Model2", group = 1), size = 0.9) +
  geom_line(aes(x = agegp, y = fit3, colour = "Model3", group = 1), size = 0.9) +
  geom_line(aes(x = agegp, y = fit4, colour = "Model4", group = 1), size = 0.9) +
  scale_color_discrete(breaks=c("Model1", "Model2", "Model3", "Model4")) +
  guides(color=guide_legend(title="Model formulation"), size = FALSE) +
  labs(title = "Overall VT carriage", x = "Age in years", y = "Prevalence") +
  theme_bw() +
  ylim(0, 0.4) +
  
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold")) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) +
  theme(legend.position = "none")

D <- ggplot(data = crude2) +
  geom_point(aes(x = year, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(year, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = year, y = fit1, colour = "Model1: β0 + ∑(βkGi) + te(age)"), size = 0.9) +
  geom_line(aes(x = year, y = fit2, colour = "Model2: β0 + ∑(βkGi) + te(time)"), size = 0.9) +
  geom_line(aes(x = year, y = fit3, colour = "Model3: β0 + ∑(βkGi) + te(age) + te(time)"), size = 0.9) +
  geom_line(aes(x = year, y = fit4, colour = "Model4: β0 + ∑(βkGi) + te(age) + te(time) + ti(age, time)"), size = 0.9) +
  scale_color_discrete(breaks=c("Model1: β0 + ∑(βkGi) + te(age)", 
                                "Model2: β0 + ∑(βkGi) + te(time)", 
                                "Model3: β0 + ∑(βkGi) + te(age) + te(time)", 
                                "Model4: β0 + ∑(βkGi) + te(age) + te(time) + ti(age, time)")) +
  guides(color=guide_legend(title=""), size=FALSE) +
  labs(title = "", x = "Survey year", y = "") +
  theme_bw() +
  ylim(0, 0.4) +
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold")) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) +
  theme(legend.position = c(0.5, 0.72))

#=======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "SFig4_Model_Selection.tiff"),
       plot = (A | B) / (C | D),
       width = 10, height = 7, unit="in", dpi = 200)
