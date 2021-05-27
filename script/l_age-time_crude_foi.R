#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(nvtcarr, agegp, surv, sex, nochild5) %>% group_by(surv) %>%
  mutate(agegpn = if_else(agegp == 19, "18-20y",
                          if_else(agegp == 22, "21-23y",
                                  if_else(agegp == 25, "24-26y",
                                          if_else(agegp == 28, "27-29y",
                                                  if_else(agegp == 31, "30-32y",
                                                          if_else(agegp == 34, "33-35y",
                                                                  if_else(agegp == 37, "36-38y", "39-40y"))))))))

#fit model to obtain predictions of FOI
model_crude = scam(nvtcarr ~ s(agegp, bs="mdcv") + s(surv, bs="mdcv") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp across surveys
A <- crude %>% filter(nvtcarr != 0) %>% group_by(surv, agegpn) %>% summarise(foi = mean(foi)) %>%
ggplot() +
  geom_point(aes(x = surv, y = foi, color = as.factor(agegpn)), size = 1.5, shape = 18) +
  geom_line(aes(x = surv, y = foi, color = as.factor(agegpn)), size = 0.7) +
  labs(title = "NVT(-ST3), Overall", x = "Survey number", y = "Force of infection") +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  ylim(0,0.02) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr, agegp, surv, sex, nochild5) %>% group_by(surv) %>%
  mutate(agegpn = if_else(agegp == 19, "18-20y",
                          if_else(agegp == 22, "21-23y",
                                  if_else(agegp == 25, "24-26y",
                                          if_else(agegp == 28, "27-29y",
                                                  if_else(agegp == 31, "30-32y",
                                                          if_else(agegp == 34, "33-35y",
                                                                  if_else(agegp == 37, "36-38y", "39-40y"))))))))

#fit model to obtain predictions of FOI
model_crude = scam(vtcarr ~ s(agegp, bs="mdcv") + s(surv, bs="mdcx") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp across surveys
B <- crude %>% filter(vtcarr != 0) %>% group_by(surv, agegpn) %>% summarise(foi = mean(foi)) %>%
ggplot() +
  geom_point(aes(x = surv, y = foi, color = as.factor(agegpn)), size = 1.5, shape = 18) +
  geom_line(aes(x = surv, y = foi, color = as.factor(agegpn)), size = 0.7) +
  labs(title = "VT(+ST3), Overall", x = "Survey number", y = "") +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  ylim(0,0.02) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(nvtcarr1, agegp, surv, sex, nochild5) %>% group_by(surv) %>%
  mutate(agegpn = if_else(agegp == 19, "18-20y",
                          if_else(agegp == 22, "21-23y",
                                  if_else(agegp == 25, "24-26y",
                                          if_else(agegp == 28, "27-29y",
                                                  if_else(agegp == 31, "30-32y",
                                                          if_else(agegp == 34, "33-35y",
                                                                  if_else(agegp == 37, "36-38y", "39-40y"))))))))

#fit model to obtain predictions of FOI
model_crude = scam(nvtcarr1 ~ s(agegp, bs="mdcv") + s(surv, bs="mpd") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp across surveys
C <- crude %>% filter(nvtcarr1 != 0) %>% group_by(surv, agegpn) %>% summarise(foi = mean(foi)) %>%
  ggplot() +
  geom_point(aes(x = surv, y = foi, color = as.factor(agegpn)), size = 1.5, shape = 18) +
  geom_line(aes(x = surv, y = foi, color = as.factor(agegpn)), size = 0.7) +
  labs(title = "NVT(+ST3), Overall", x = "Survey number", y = "") +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  ylim(0,0.02) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr1, agegp, surv, sex, nochild5) %>% group_by(surv) %>%
  mutate(agegpn = if_else(agegp == 19, "18-20y",
                          if_else(agegp == 22, "21-23y",
                                  if_else(agegp == 25, "24-26y",
                                          if_else(agegp == 28, "27-29y",
                                                  if_else(agegp == 31, "30-32y",
                                                          if_else(agegp == 34, "33-35y",
                                                                  if_else(agegp == 37, "36-38y", "39-40y"))))))))

#fit model to obtain predictions of FOI
model_crude = scam(vtcarr1 ~ s(agegp, bs="mdcx") + s(surv, bs="mdcx") + sex + nochild5, family = binomial(link = "cloglog"), data = crude)
crude$foi <- ((-derivative.scam(model_crude, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp across surveys
D <- crude %>% filter(vtcarr1 != 0) %>% group_by(surv, agegpn) %>% summarise(foi = mean(foi)) %>%
  ggplot() +
  geom_point(aes(x = surv, y = foi, color = as.factor(agegpn)), size = 1.5, shape = 18) +
  geom_line(aes(x = surv, y = foi, color = as.factor(agegpn)), size = 0.7) +
  labs(title = "VT(-ST3), Overall", x = "Survey number", y = "") +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  ylim(0,0.02) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "right") + 
  guides(color = guide_legend(title = "Age group"))

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "SFig3_sens_age-time.tiff"),
       plot = (A | B | C | D),
       width = 14, height = 3, unit="in", dpi = 200)
