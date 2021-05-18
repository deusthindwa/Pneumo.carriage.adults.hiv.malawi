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

#fit model to obtain predictions of FOI
model_crude = scam(nvtcarr ~ s(agegp, bs="mdcv") + s(surv, bs="mdcv"), family = binomial(link = "cloglog"), data = crude)
crude$foi1 <- ((-derivative.scam(model_crude, smooth.number = 1, deriv = 1)$d * model_crude$fitted.values) + (0*model_crude$fitted.values))/(1-model_crude$fitted.values)
crude$foi2 <- ((-derivative.scam(model_crude, smooth.number = 1, deriv = 1)$d * model_crude$fitted.values) + (1/11*model_crude$fitted.values))/(1-model_crude$fitted.values)
crude$foi3 <- ((-derivative.scam(model_crude, smooth.number = 1, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(agegp) %>% summarise(afoi1 = mean(foi1), afoi2 = mean(foi2), afoi3 = mean(foi3)) %>%
ungroup()) %>% select(agegp, afoi1, afoi2, afoi3)

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(nvtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(nvtcarr != 0) %>% group_by(surv) %>% summarise(tfoi1 = mean(foi1), tfoi2 = mean(foi2), tfoi3 = mean(foi3)) %>%
ungroup()) %>% select(surv, tfoi1, tfoi2, tfoi3)

#plot prevalence curves
A <- ggplot(data = cbind(crude1, crude2)) +
  geom_line(aes(x = agegp, y = afoi1), lty = "dashed", size = 0.7, color = "darkgreen") +
  geom_line(aes(x = agegp, y = afoi2), lty = "dashed", size = 0.7, color = "darkblue") +
  geom_line(aes(x = agegp, y = afoi3), lty = "dashed", size = 0.7, color = "darkred") +
  ylim(0, 0.08) +
  labs(title = "NVT(-ST3), Overall", x = "Age,y", y = "Force of infection") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

B <- ggplot(data = cbind(crude1, crude2)) +
  geom_line(aes(x = surv, y = tfoi1), lty = "dashed", size = 0.7, color = "darkgreen") +
  geom_line(aes(x = surv, y = tfoi2), lty = "dashed", size = 0.7, color = "darkblue") +
  geom_line(aes(x = surv, y = tfoi3), lty = "dashed", size = 0.7, color = "darkred") +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  ylim(0, 0.08) +
  labs(title = "NVT(-ST3), Overall", x = "Survey number", y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit model to obtain predictions of FOI
model_crude = scam(vtcarr ~ s(agegp, bs="mdcx") + s(surv, bs="mdcx"), family = binomial(link = "cloglog"), data = crude)
crude$foi1 <- ((-derivative.scam(model_crude, smooth.number = 1, deriv = 1)$d * model_crude$fitted.values) + (0*model_crude$fitted.values))/(1-model_crude$fitted.values)
crude$foi2 <- ((-derivative.scam(model_crude, smooth.number = 1, deriv = 1)$d * model_crude$fitted.values) + (1/11*model_crude$fitted.values))/(1-model_crude$fitted.values)
crude$foi3 <- ((-derivative.scam(model_crude, smooth.number = 1, deriv = 1)$d * model_crude$fitted.values) + (1/42*model_crude$fitted.values))/(1-model_crude$fitted.values)

#join observed and predicted datasets for agegp
crude1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(afoi1 = mean(foi1), afoi2 = mean(foi2), afoi3 = mean(foi3)) %>%
ungroup()) %>% select(agegp, afoi1, afoi2, afoi3)

#join observed and predicted datasets for survey number
crude2 <- left_join(left_join(crude %>% group_by(surv) %>% tally() %>% rename(Tot = n), 
crude %>% filter(vtcarr == 1) %>% group_by(surv) %>% tally() %>% rename(Pos = n)), 
crude %>% filter(vtcarr != 0) %>% group_by(surv) %>% summarise(tfoi1 = mean(foi1), tfoi2 = mean(foi2), tfoi3 = mean(foi3)) %>%
ungroup()) %>% select(surv, tfoi1, tfoi2, tfoi3)

#plot prevalence curves
C <- ggplot(data = cbind(crude1, crude2)) +
  geom_line(aes(x = agegp, y = afoi1), lty = "dashed", size = 0.7, color = "darkgreen") +
  geom_line(aes(x = agegp, y = afoi2), lty = "dashed", size = 0.7, color = "darkblue") +
  geom_line(aes(x = agegp, y = afoi3), lty = "dashed", size = 0.7, color = "darkred") +
  ylim(0, 0.08) +
  labs(title = "VT(+ST3), Overall", x = "Age,y", y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none")

D <- ggplot(data = cbind(crude1, crude2)) +
  geom_line(aes(x = surv, y = tfoi1, color = "Immunising infection"), lty = "dashed", size = 0.7) +
  geom_line(aes(x = surv, y = tfoi2, color = "11 days (Thindwa et al.)"), lty = "dashed", size = 0.7) +
  geom_line(aes(x = surv, y = tfoi3, color = "42 days (Kalata et al.)"), lty = "dashed", size = 0.7) +
  scale_colour_manual(name = "Carriage duration", values = c("Immunising infection" = "darkgreen", "11 days (Thindwa et al.)" = "darkblue", "42 days (Kalata et al.)"="darkred")) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  ylim(0, 0.08) +
  labs(title = "VT(+ST3), Overall", x = "Survey number", y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = c(0.6, 0.7))

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "SFig3_sens_carr_dur.tiff"),
       plot = (A | B | C | D),
       width = 13, height = 4, unit="in", dpi = 200)
