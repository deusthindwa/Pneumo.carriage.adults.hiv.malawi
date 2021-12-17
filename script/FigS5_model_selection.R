#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 17/12/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1) 

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(carr, age, year)

#fit different models & obtain predictions and 95%CI
model_crude = gam(carr ~ te(age, bs="cs") + te(year, bs="cs") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fitcs = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(carr ~ te(age, bs="tp") + te(year, bs="tp") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fittp = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(carr ~ te(age, bs="cr") + te(year, bs="cr") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fitcr = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

model_crude = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fitps = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24y",
                                   if_else(age >24 & age <=29, "25-29y",
                                           if_else(age >29 & age <=34, "30-34y",
                                                   if_else(age >34 & age <=40, "35-40y", NA_character_))))),
         year = as.factor(year))

A <- rbind(
#join observed and predicted datasets for agegp
left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(carr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(carr != 0) %>% group_by(agegp) %>% summarise(fitcs = mean(fitcs), fittp = mean(fittp), fitcr = mean(fitcr), fitps = mean(fitps)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], cat = "VT(+st3) carriage", status = "A, Age group (years), Spline type") %>% rename("covar" = agegp),

#join observed and predicted datasets for survey number
left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
crude %>% filter(carr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(carr != 0) %>% group_by(year) %>% summarise(fitcs = mean(fitcs), fittp = mean(fittp), fitcr = mean(fitcr), fitps = mean(fitps)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], cat = "VT(+st3) carriage", status = "B, Time (years), Spline type") %>% rename("covar" = year)
) %>%

  #plot prevalence curves
  ggplot() +
  geom_point(aes(x = covar, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(covar, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = covar, y = fitcs, colour = "Natural cubic spline", group = 1), size = 0.9) +
  geom_line(aes(x = covar, y = fittp, colour = "Thin plate regression spline", group = 1), size = 0.9) +
  geom_line(aes(x = covar, y = fitps, colour = "P-spline", group = 1), size = 0.9) +
  geom_line(aes(x = covar, y = fitcr, colour = "Cubic regression spline", group = 1), size = 0.9) +
  scale_color_discrete(breaks=c("Natural cubic spline", "Thin plate regression spline", "P-spline", "Cubic regression spline")) +
  facet_grid(.~status, scales = "free") +
  guides(color = guide_legend(title=""), size = FALSE) +
  labs(title = "", x = "", y = "Overall carriage prevalence") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14), strip.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0,0.55), labels = scales::percent_format(accuracy = 1)) +  
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold")) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) +
  theme(legend.position = "top") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#=======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(carr, age, year)

#fit different model formulations & obtain predictions and 95%CI
model1_crude = gam(carr ~ te(age, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fit1 = predict.gam(model1_crude, type = "response", se.fit = TRUE)$fit

model2_crude = gam(carr ~ te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fit2 = predict.gam(model2_crude, type = "response", se.fit = TRUE)$fit

model3_crude = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fit3 = predict.gam(model3_crude, type = "response", se.fit = TRUE)$fit

model4_crude = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + ti(age, year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$fit4 = predict.gam(model4_crude, type = "response", se.fit = TRUE)$fit

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24y",
                                   if_else(age >24 & age <=29, "25-29y",
                                           if_else(age >29 & age <=34, "30-34y",
                                                   if_else(age >34 & age <=40, "35-40y", NA_character_))))),
         year = as.factor(year))

B <- rbind(
#join observed and predicted datasets for agegp
left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(carr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(carr != 0) %>% group_by(agegp) %>% summarise(fit1 = mean(fit1), fit2 = mean(fit2), fit3 = mean(fit3), fit4 = mean(fit4)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], cat = "VT(+st3) carriage", status = "C, Age group (years), Model variant") %>% rename("covar" = agegp),

#join observed and predicted datasets for survey number
left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
crude %>% filter(carr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(carr != 0) %>% group_by(year) %>% summarise(fit1 = mean(fit1), fit2 = mean(fit2), fit3 = mean(fit3), fit4 = mean(fit4)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], cat = "VT(+st3) carriage", status = "D, Time (years), Model variant") %>% rename("covar" = year)
) %>%
  
  ggplot() +
  geom_point(aes(x = covar, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(covar, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = covar, y = fit1, colour = "β0 + ∑(βkGi) + te(age)", group = 1), size = 0.9) + 
  geom_line(aes(x = covar, y = fit2, colour = "β0 + ∑(βkGi) + te(time)", group = 1), size = 0.9) +
  geom_line(aes(x = covar, y = fit3, colour = "β0 + ∑(βkGi) + te(age) + te(time)", group = 1), size = 0.9) +
  geom_line(aes(x = covar, y = fit4, colour = "β0 + ∑(βkGi) + te(age) + te(time) + ti(age, time)", group = 1), size = 0.9) +
  scale_color_discrete(breaks=c("β0 + ∑(βkGi) + te(age)", "β0 + ∑(βkGi) + te(time)","β0 + ∑(βkGi) + te(age) + te(time)","β0 + ∑(βkGi) + te(age) + te(time) + ti(age, time)" )) +
  facet_grid(.~status, scales = "free") +
  guides(color = guide_legend(title = ""), size = FALSE) + 
  labs(title = "", x = "", y = "") +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14), strip.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 0.55), labels = scales::percent_format(accuracy = 1)) +  
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_blank()) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "top") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#=======================================================================================


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
  mutate(agegp = as.factor(if_else(age <=24, "18-24y",
                                   if_else(age >24 & age <=29, "25-29y",
                                           if_else(age >29 & age <=34, "30-34y",
                                                   if_else(age >34 & age <=40, "35-40y", NA_character_))))),
         year = as.factor(year))

C <- rbind(
#join observed and predicted datasets for agegp
left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fitcs = mean(fitcs), fittp = mean(fittp), fitcr = mean(fitcr), fitps = mean(fitps)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], cat = "VT(+st3) carriage", status = "E, Age group (years), Spline type") %>% rename("covar" = agegp),

#join observed and predicted datasets for survey number
left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fitcs = mean(fitcs), fittp = mean(fittp), fitcr = mean(fitcr), fitps = mean(fitps)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], cat = "VT(+st3) carriage", status = "F, Time (years), Spline type") %>% rename("covar" = year)
) %>%

#plot prevalence curves
ggplot() +
  geom_point(aes(x = covar, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(covar, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = covar, y = fitcs, colour = "Natural cubic spline", group = 1), size = 0.9) +
  geom_line(aes(x = covar, y = fittp, colour = "Thin plate regression spline", group = 1), size = 0.9) +
  geom_line(aes(x = covar, y = fitps, colour = "P-spline", group = 1), size = 0.9) +
  geom_line(aes(x = covar, y = fitcr, colour = "Cubic regression spline", group = 1), size = 0.9) +
  scale_color_discrete(breaks=c("Natural cubic spline", "Thin plate regression spline", "P-spline", "Cubic regression spline")) +
  facet_grid(.~status, scales = "free") +
  guides(color = guide_legend(title=""), size = FALSE) +
  labs(title = "", x = "", y = "VT carriage prevalence") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14), strip.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks = seq(0, 1, 0.05), limits = c(0,0.25), labels = scales::percent_format(accuracy = 1)) +  
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold")) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

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
  mutate(agegp = as.factor(if_else(age <=24, "18-24y",
                                   if_else(age >24 & age <=29, "25-29y",
                                           if_else(age >29 & age <=34, "30-34y",
                                                   if_else(age >34 & age <=40, "35-40y", NA_character_))))),
         year = as.factor(year))

D <- rbind(
#join observed and predicted datasets for agegp
left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit1 = mean(fit1), fit2 = mean(fit2), fit3 = mean(fit3), fit4 = mean(fit4)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], cat = "VT(+st3) carriage", status = "G, Age group (years), Model variant") %>% rename("covar" = agegp),

#join observed and predicted datasets for survey number
left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit1 = mean(fit1), fit2 = mean(fit2), fit3 = mean(fit3), fit4 = mean(fit4)) %>%
ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], cat = "VT(+st3) carriage", status = "H, Time (years), Model variant") %>% rename("covar" = year)
) %>%

ggplot() +
  geom_point(aes(x = covar, y = obs, size = Pos), shape = 1) +
  geom_errorbar(aes(covar, ymin = obs_lci, ymax = obs_uci), width = 0, size = 0.3) +
  geom_line(aes(x = covar, y = fit1, colour = "β0 + ∑(βkGi) + te(age)", group = 1), size = 0.9) + 
  geom_line(aes(x = covar, y = fit2, colour = "β0 + ∑(βkGi) + te(time)", group = 1), size = 0.9) +
  geom_line(aes(x = covar, y = fit3, colour = "β0 + ∑(βkGi) + te(age) + te(time)", group = 1), size = 0.9) +
  geom_line(aes(x = covar, y = fit4, colour = "β0 + ∑(βkGi) + te(age) + te(time) + ti(age, time)", group = 1), size = 0.9) +
  scale_color_discrete(breaks=c("β0 + ∑(βkGi) + te(age)", "β0 + ∑(βkGi) + te(time)", "β0 + ∑(βkGi) + te(age) + te(time)","β0 + ∑(βkGi) + te(age) + te(time) + ti(age, time)")) +
  facet_grid(.~status, scales = "free") +
  guides(color = guide_legend(title = ""), size = FALSE) +
  labs(title = "", x = "", y = "") +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14), strip.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks = seq(0, 1, 0.05), limits = c(0, 0.25), labels = scales::percent_format(accuracy = 1)) +  
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_blank()) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#=======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "FigS5_model_selection.png"),
       plot = ((A | B)/(C | D)),
       width = 18, height = 8, unit="in", dpi = 300)
