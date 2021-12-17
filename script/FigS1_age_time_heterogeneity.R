#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 17/12/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(carr, vtcarr, age, year)

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24y",
                                   if_else(age >24 & age <=29, "25-29y",
                                           if_else(age >29 & age <=34, "30-34y",
                                                   if_else(age >34 & age <=40, "35-40y", NA_character_))))))

#fit model to individual trajectories & obtain predictions and 95%CI for overall carriage
model0 = gam(carr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$prev0 <- model0$fitted.values

#fit model to individual trajectories & obtain predictions and 95%CI for VT carriage
model1 = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$prev1 <- model1$fitted.values

#======================================================================================

A <- rbind(
  crude %>% filter(carr != 0) %>% group_by(year, agegp) %>% summarise(prev = mean(prev0)) %>% mutate(cat = "A, Overall carriage across age"),
  crude %>% filter(vtcarr != 0) %>% group_by(year, agegp) %>% summarise(prev = mean(prev1)) %>% mutate(cat = "B, VT carriage across age")
) %>%
  
  ggplot() +
  geom_point(aes(x = agegp, y = prev, color = as.factor(year), group = year), size = 3, shape = 18) +
  geom_line(aes(x = agegp, y = prev, color = as.factor(year), group = year), size = 0.9) +
  labs(title = "", x = "Age group (years)", y = "Carriage prevalence") +
  theme_bw() +
  facet_grid(.~cat) +
  theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16), strip.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +  
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold")) +
  theme(plot.title = element_text(size = 20), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "top", legend.title=element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))


B <- rbind(
crude %>% filter(carr != 0) %>% group_by(agegp, year) %>% summarise(prev = mean(prev0))  %>% mutate(cat = "C, Overall carriage across time"),
crude %>% filter(vtcarr != 0) %>% group_by(agegp, year) %>% summarise(prev = mean(prev1)) %>% mutate(cat = "D, VT carriage across time")
) %>%

  ggplot() +
  geom_point(aes(x = year, y = prev, color = as.factor(agegp)), size = 3, shape = 18) +
  geom_line(aes(x = year, y = prev, color = as.factor(agegp)), size = 0.9) +
  labs(title = "", x = "Time (in years)", y = "") +
  theme_bw() +
  facet_grid(.~cat) +
  theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16), strip.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +  
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_blank()) +
  theme(plot.title = element_text(size = 20), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
  theme(legend.position = "top", legend.title=element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "FigS1_age_time_heterogeneity.png"),
       plot = (A | B),
       width = 15, height = 5, unit="in", dpi = 300)
