#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit model & obtain prevalence predictions and 95%CI
model_crude = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = pcvpa.mod, na.action = na.exclude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#divide by carriage duration to obtaine acquisitions and 95%CI
crude$foi <- crude$fit/42; crude$foi_lci <- crude$fit_lci/42; crude$foi_uci <- crude$fit_uci/42
crude1 <- crude %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(agegp) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "VT(+ST3)") %>% ungroup()
crude2 <- crude %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(surv) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "VT(+ST3)") %>% ungroup()

#----------------------------------------------------------------------------------

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr1, agegp, surv) %>% group_by(surv)

#fit model & obtain prevalence predictions and 95%CI
model_crude = gam(vtcarr1 ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = pcvpa.mod, na.action = na.exclude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#divide by carriage duration to obtaine acquisitions and 95%CI
crude$foi <- crude$fit/42; crude$foi_lci <- crude$fit_lci/42; crude$foi_uci <- crude$fit_uci/42
crude3 <- crude %>% filter(vtcarr1 == 1 & !is.na(fit)) %>% group_by(agegp) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "VT(-ST3)") %>% ungroup()
crude4 <- crude %>% filter(vtcarr1 == 1 & !is.na(fit)) %>% group_by(surv) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "VT(-ST3)") %>% ungroup()

#----------------------------------------------------------------------------------

#plot prevalence curves
A <- ggplot(data = rbind(crude1, crude3)) +
  geom_line(aes(x = agegp, y = foi*365.25, color = strat, lty = strat), size = 1) +
  geom_ribbon(aes(x = agegp, y = foi*365.25, ymin = foi_lci*365.25, ymax = foi_uci*365.25, color = strat, fill = strat, lty = strat), alpha = 0.3) +
  labs(title = "Overall carriage", x = "Age,y", y = "Annual carriage acquisition") +
  theme_bw() +
  ylim(0, 6) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none") +
  guides(color = guide_legend(title = ""), lty = guide_legend(title = ""), fill = guide_legend(title = ""))

B <- ggplot(data = rbind(crude2, crude4)) +
  geom_line(aes(x = surv, y = foi*365.25, color = strat, lty = strat), size = 1) +
  geom_ribbon(aes(x = surv, y = foi*365.25, ymin = foi_lci*365.25, ymax = foi_uci*365.25, color = strat, fill = strat, lty = strat), alpha = 0.3) +
  labs(title = "", x = "Survey number", y = "") +
  theme_bw() +
  ylim(0, 6) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = c(0.6, 0.6)) +
  guides(color = guide_legend(title = ""), lty = guide_legend(title = ""), fill = guide_legend(title = ""))

#======================================================================================

#subset for a correct dataset
male = filter(pcvpa.mod, sex == 1) %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit model & obtain prevalence predictions and 95%CI
model_male = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 1), na.action = na.exclude)
male$fit = predict.gam(model_male, type = "response", se.fit = TRUE)$fit
male$fit_lci = model_male$family$linkinv(predict.gam(model_male, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_male, type = "link", se.fit = TRUE)$se.fit))
male$fit_uci = model_male$family$linkinv(predict.gam(model_male, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_male, type = "link", se.fit = TRUE)$se.fit))

#divide by carriage duration to obtaine acquisitions and 95%CI
male$foi <- male$fit/42; male$foi_lci <- male$fit_lci/42; male$foi_uci <- male$fit_uci/42
crude1 <- male %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(agegp) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "Male") %>% ungroup()
crude2 <- male %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(surv) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "Male") %>% ungroup()

#----------------------------------------------------------------------------------

#subset for a correct dataset
female = filter(pcvpa.mod, sex == 0) %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit model & obtain prevalence predictions and 95%CI
model_female = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sex == 0), na.action = na.exclude)
female$fit = predict.gam(model_female, type = "response", se.fit = TRUE)$fit
female$fit_lci = model_female$family$linkinv(predict.gam(model_female, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_female, type = "link", se.fit = TRUE)$se.fit))
female$fit_uci = model_female$family$linkinv(predict.gam(model_female, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_female, type = "link", se.fit = TRUE)$se.fit))

#divide by carriage duration to obtaine acquisitions and 95%CI
female$foi <- female$fit/42; female$foi_lci <- female$fit_lci/42; female$foi_uci <- female$fit_uci/42
crude3 <- female %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(agegp) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "Female") %>% ungroup()
crude4 <- female %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(surv) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "Female") %>% ungroup()

#----------------------------------------------------------------------------------

#plot prevalence curves
C <- ggplot(data = rbind(crude1, crude3)) +
  geom_line(aes(x = agegp, y = foi*365.25, color = strat, lty = strat), size = 1) +
  geom_ribbon(aes(x = agegp, y = foi*365.25, ymin = foi_lci*365.25, ymax = foi_uci*365.25, color = strat, fill = strat, lty = strat), alpha = 0.3) +
  labs(title = "Sex", x = "Age,y", y = "Annual carriage acquisition") +
  theme_bw() +
  ylim(0, 6) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none") +
  guides(color = guide_legend(title = ""), lty = guide_legend(title = ""), fill = guide_legend(title = ""))

D <- ggplot(data = rbind(crude2, crude4)) +
  geom_line(aes(x = surv, y = foi*365.25, color = strat, lty = strat), size = 1) +
  geom_ribbon(aes(x = surv, y = foi*365.25, ymin = foi_lci*365.25, ymax = foi_uci*365.25, color = strat, fill = strat, lty = strat), alpha = 0.3) +
  labs(title = "", x = "Survey number", y = "") +
  theme_bw() +
  ylim(0, 6) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = c(0.6, 0.6)) +
  guides(color = guide_legend(title = ""), lty = guide_legend(title = ""), fill = guide_legend(title = ""))

#======================================================================================

#subset for a correct dataset
lart = filter(pcvpa.mod, artdur == 1) %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit model & obtain prevalence predictions and 95%CI
model_lart = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 1), na.action = na.exclude)
lart$fit = predict.gam(model_lart, type = "response", se.fit = TRUE)$fit
lart$fit_lci = model_lart$family$linkinv(predict.gam(model_lart, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lart, type = "link", se.fit = TRUE)$se.fit))
lart$fit_uci = model_lart$family$linkinv(predict.gam(model_lart, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lart, type = "link", se.fit = TRUE)$se.fit))

#divide by carriage duration to obtaine acquisitions and 95%CI
lart$foi <- lart$fit/42; lart$foi_lci <- lart$fit_lci/42; lart$foi_uci <- lart$fit_uci/42
crude1 <- lart %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(agegp) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "ART >3y") %>% ungroup()
crude2 <- lart %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(surv) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "ART >3y") %>% ungroup()

#----------------------------------------------------------------------------------

#subset for a correct dataset
sart = filter(pcvpa.mod, artdur == 0) %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit model & obtain prevalence predictions and 95%CI
model_sart = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, artdur == 0), na.action = na.exclude)
sart$fit = predict.gam(model_sart, type = "response", se.fit = TRUE)$fit
sart$fit_lci = model_sart$family$linkinv(predict.gam(model_sart, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_sart, type = "link", se.fit = TRUE)$se.fit))
sart$fit_uci = model_sart$family$linkinv(predict.gam(model_sart, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_sart, type = "link", se.fit = TRUE)$se.fit))

#divide by carriage duration to obtaine acquisitions and 95%CI
sart$foi <- sart$fit/42; sart$foi_lci <- sart$fit_lci/42; sart$foi_uci <- sart$fit_uci/42
crude3 <- sart %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(agegp) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "ART ≤3y") %>% ungroup()
crude4 <- sart %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(surv) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "ART ≤3y") %>% ungroup()

#----------------------------------------------------------------------------------

#plot prevalence curves
E <- ggplot(data = rbind(crude1, crude3)) +
  geom_line(aes(x = agegp, y = foi*365.25, color = strat, lty = strat), size = 1) +
  geom_ribbon(aes(x = agegp, y = foi*365.25, ymin = foi_lci*365.25, ymax = foi_uci*365.25, color = strat, fill = strat, lty = strat), alpha = 0.3) +
  labs(title = "ART duration", x = "Age,y", y = "") +
  theme_bw() +
  ylim(0, 6) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none") +
  guides(color = guide_legend(title = ""), lty = guide_legend(title = ""), fill = guide_legend(title = ""))

F <- ggplot(data = rbind(crude2, crude4)) +
  geom_line(aes(x = surv, y = foi*365.25, color = strat, lty = strat), size = 1) +
  geom_ribbon(aes(x = surv, y = foi*365.25, ymin = foi_lci*365.25, ymax = foi_uci*365.25, color = strat, fill = strat, lty = strat), alpha = 0.3) +
  labs(title = "", x = "Survey number", y = "") +
  theme_bw() +
  ylim(0, 6) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = c(0.6, 0.6)) +
  guides(color = guide_legend(title = ""), lty = guide_legend(title = ""), fill = guide_legend(title = ""))

#======================================================================================

#subset for a correct dataset
hcd4 = filter(pcvpa.mod, cd4cnt == 1) %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit model & obtain prevalence predictions and 95%CI
model_hcd4 = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, cd4cnt == 1), na.action = na.exclude)
hcd4$fit = predict.gam(model_hcd4, type = "response", se.fit = TRUE)$fit
hcd4$fit_lci = model_hcd4$family$linkinv(predict.gam(model_hcd4, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_hcd4, type = "link", se.fit = TRUE)$se.fit))
hcd4$fit_uci = model_hcd4$family$linkinv(predict.gam(model_hcd4, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_hcd4, type = "link", se.fit = TRUE)$se.fit))

#divide by carriage duration to obtaine acquisitions and 95%CI
hcd4$foi <- hcd4$fit/42; hcd4$foi_lci <- hcd4$fit_lci/42; hcd4$foi_uci <- hcd4$fit_uci/42
crude1 <- hcd4 %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(agegp) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "≥250 cells/μL") %>% ungroup()
crude2 <- hcd4 %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(surv) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "≥250 cells/μL") %>% ungroup()

#----------------------------------------------------------------------------------

#subset for a correct dataset
lcd4 = filter(pcvpa.mod, cd4cnt == 0) %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit model & obtain prevalence predictions and 95%CI
model_lcd4 = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, cd4cnt == 0), na.action = na.exclude)
lcd4$fit = predict.gam(model_lcd4, type = "response", se.fit = TRUE)$fit
lcd4$fit_lci = model_lcd4$family$linkinv(predict.gam(model_lcd4, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lcd4, type = "link", se.fit = TRUE)$se.fit))
lcd4$fit_uci = model_lcd4$family$linkinv(predict.gam(model_lcd4, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lcd4, type = "link", se.fit = TRUE)$se.fit))

#divide by carriage duration to obtaine acquisitions and 95%CI
lcd4$foi <- lcd4$fit/42; lcd4$foi_lci <- lcd4$fit_lci/42; lcd4$foi_uci <- lcd4$fit_uci/42
crude3 <- lcd4 %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(agegp) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "<250 cells/μL") %>% ungroup()
crude4 <- lcd4 %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(surv) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "<250 cells/μL") %>% ungroup()

#----------------------------------------------------------------------------------

#plot prevalence curves
G <- ggplot(data = rbind(crude1, crude3)) +
  geom_line(aes(x = agegp, y = foi*365.25, color = strat, lty = strat), size = 1) +
  geom_ribbon(aes(x = agegp, y = foi*365.25, ymin = foi_lci*365.25, ymax = foi_uci*365.25, color = strat, fill = strat, lty = strat), alpha = 0.3) +
  labs(title = "CD4 count", x = "Age,y", y = "Annual carriage acquisition") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 6)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none") +
  guides(color = guide_legend(title = ""), lty = guide_legend(title = ""), fill = guide_legend(title = ""))

H <- ggplot(data = rbind(crude2, crude4)) +
  geom_line(aes(x = surv, y = foi*365.25, color = strat, lty = strat), size = 1) +
  geom_ribbon(aes(x = surv, y = foi*365.25, ymin = foi_lci*365.25, ymax = foi_uci*365.25, color = strat, fill = strat, lty = strat), alpha = 0.3) +
  labs(title = "", x = "Survey number", y = "") +
  theme_bw() +
  ylim(0, 6) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = c(0.6, 0.6)) +
  guides(color = guide_legend(title = ""), lty = guide_legend(title = ""), fill = guide_legend(title = ""))

#======================================================================================

#subset for a correct dataset
hses = filter(pcvpa.mod, sescat == 1) %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit model & obtain prevalence predictions and 95%CI
model_hses = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 1), na.action = na.exclude)
hses$fit = predict.gam(model_hses, type = "response", se.fit = TRUE)$fit
hses$fit_lci = model_hses$family$linkinv(predict.gam(model_hses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_hses, type = "link", se.fit = TRUE)$se.fit))
hses$fit_uci = model_hses$family$linkinv(predict.gam(model_hses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_hses, type = "link", se.fit = TRUE)$se.fit))

#divide by carriage duration to obtaine acquisitions and 95%CI
hses$foi <- hses$fit/42; hses$foi_lci <- hses$fit_lci/42; hses$foi_uci <- hses$fit_uci/42
crude1 <- hses %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(agegp) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "High/Middle SES") %>% ungroup()
crude2 <- hses %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(surv) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "High/Middle SES") %>% ungroup()

#----------------------------------------------------------------------------------

#subset for a correct dataset
lses = filter(pcvpa.mod, sescat == 0) %>% select(vtcarr, agegp, surv) %>% group_by(surv)

#fit model & obtain prevalence predictions and 95%CI
model_lses = gam(vtcarr ~ te(agegp, bs="ps") + te(surv, bs="ps") + sex + nochild5, family = binomial(link = "cloglog"), data = filter(pcvpa.mod, sescat == 0), na.action = na.exclude)
lses$fit = predict.gam(model_lses, type = "response", se.fit = TRUE)$fit
lses$fit_lci = model_lses$family$linkinv(predict.gam(model_lses, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_lses, type = "link", se.fit = TRUE)$se.fit))
lses$fit_uci = model_lses$family$linkinv(predict.gam(model_lses, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_lses, type = "link", se.fit = TRUE)$se.fit))

#divide by carriage duration to obtaine acquisitions and 95%CI
lses$foi <- lses$fit/42; lses$foi_lci <- lses$fit_lci/42; lses$foi_uci <- lses$fit_uci/42
crude3 <- lses %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(agegp) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "Low SES") %>% ungroup()
crude4 <- lses %>% filter(vtcarr == 1 & !is.na(fit)) %>% group_by(surv) %>% summarise(foi = mean(foi), foi_lci = mean(foi_lci), foi_uci = mean(foi_uci)) %>% mutate(strat = "Low SES") %>% ungroup()

#----------------------------------------------------------------------------------

#plot prevalence curves
I <- ggplot(data = rbind(crude1, crude3)) +
  geom_line(aes(x = agegp, y = foi*365.25, color = strat, lty = strat), size = 1) +
  geom_ribbon(aes(x = agegp, y = foi*365.25, ymin = foi_lci*365.25, ymax = foi_uci*365.25, color = strat, fill = strat, lty = strat), alpha = 0.3) +
  labs(title = "CD4 count", x = "Age,y", y = "") +
  theme_bw() +
  ylim(0, 8) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "none") +
  guides(color = guide_legend(title = ""), lty = guide_legend(title = ""), fill = guide_legend(title = ""))

J <- ggplot(data = rbind(crude2, crude4)) +
  geom_line(aes(x = surv, y = foi*365.25, color = strat, lty = strat), size = 1) +
  geom_ribbon(aes(x = surv, y = foi*365.25, ymin = foi_lci*365.25, ymax = foi_uci*365.25, color = strat, fill = strat, lty = strat), alpha = 0.3) +
  labs(title = "", x = "Survey number", y = "") +
  theme_bw() +
  ylim(0, 8) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = c(0.6, 0.6)) +
  guides(color = guide_legend(title = ""), lty = guide_legend(title = ""), fill = guide_legend(title = ""))

#======================================================================================

(A | B | plot_layout(ncol = 2, width = c(1,1,1,1))) / (C | D | E | F) / (G | H | I | J)

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig2_CrudePrevFOI.tiff"),
       plot = (A | B | C | D) / (E | F | G | H),
       width = 14, height = 6, unit="in", dpi = 200)
