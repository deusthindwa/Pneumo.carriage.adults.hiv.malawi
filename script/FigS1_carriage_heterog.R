#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr, nvtcarr, age, year)

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24y",
                                   if_else(age >24 & age <=29, "25-29y",
                                           if_else(age >29 & age <=34, "30-34y",
                                                   if_else(age >34 & age <=40, "35-40y", NA_character_))))))

#fit model to individual trajectories & obtain predictions and 95%CI
model_vt = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)

crude$prevt <- model_vt$fitted.values
crude$prevtl = model_vt$family$linkinv(predict.gam(model_vt, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_vt, type = "link", se.fit = TRUE)$se.fit))
crude$prevtu = model_vt$family$linkinv(predict.gam(model_vt, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_vt, type = "link", se.fit = TRUE)$se.fit))

crude <- crude %>% mutate(acqvt1 = prevt/11, 
                          acqvt1l = prevtl/11,
                          acqvt1u = prevtu/11,
                          acqvt2 = prevt/42, 
                          acqvt2l = prevtl/42,
                          acqvt2u = prevtu/42)

#join observed and predicted datasets for survey year
A <- rbind(
rbind((crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(acq = mean(acqvt1), acql = mean(acqvt1l), acqu = mean(acqvt1u), dur = "11 days")),
(crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(acq = mean(acqvt2), acql = mean(acqvt2l), acqu = mean(acqvt2u), dur = "42 days"))) %>% mutate(catg = "VT carriage by year") %>% rename("varx" = year),

rbind((crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(acq = mean(acqvt1), acql = mean(acqvt1l), acqu = mean(acqvt1u), dur = "11 days")),
(crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(acq = mean(acqvt2), acql = mean(acqvt2l), acqu = mean(acqvt2u), dur = "42 days"))) %>% mutate(catg = "VT carriage by age") %>% rename("varx" = agegp)

) %>% mutate(varx = as_factor(varx), dur = as_factor(dur), catg = as_factor(catg)) %>%

  ggplot() +
  geom_line(aes(x = varx, y = acq, color = dur, group = dur), lty = "dashed", size = 0.9) +
  geom_ribbon(aes(x = varx, y = acq, group = dur, ymin = acql, ymax = acqu, fill = dur), alpha = 0.2) +
  facet_grid(.~catg, scales = "free_x") +
  theme_bw() +
  ylim(0, 0.04) +
  labs(title = "", x = "", y = "Daily carriage acquisition") +
  guides(color=guide_legend(title="Assumed carriage duration"), fill = FALSE) +
  theme(strip.text.x = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold")) +
  theme(axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) +
  theme(legend.position = c(0.8, 0.75))

#======================================================================================

#join observed and predicted datasets for agegp across surveys
B <- crude %>% filter(vtcarr != 0) %>% group_by(year, agegp) %>% summarise(foi = mean(acqvt2)) %>%
  ggplot() +
  geom_point(aes(x = year, y = foi, color = as.factor(agegp)), size = 3, shape = 18) +
  geom_line(aes(x = year, y = foi, color = as.factor(agegp)), size = 0.9) +
  labs(title = "", x = "Survey number", y = "Daily carriage acquisition") +
  theme_bw() +
  ylim(0, 0.007) +
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold")) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) +
  theme(legend.position = c(0.8, 0.75), legend.title=element_blank())

C <- crude %>% filter(vtcarr != 0) %>% group_by(agegp, year) %>% summarise(foi = mean(acqvt2)) %>%
  ggplot() +
  geom_point(aes(x = agegp, y = foi, color = as.factor(year), group = year), size = 3, shape = 18) +
  geom_line(aes(x = agegp, y = foi, color = as.factor(year), group = year), size = 0.9) +
  labs(title = "", x = "Age group (years)", y = "") +
  theme_bw() +
  ylim(0, 0.007) +
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_blank()) +
  theme(plot.title = element_text(size = 14), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11)) +
  theme(legend.position = c(0.8, 0.75), legend.title=element_blank())

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "FigS1_carriage_heterog.tiff"),
       plot = (A) / (B | C),
       width = 11, height = 8, unit="in", dpi = 200)
