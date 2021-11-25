#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)
dur = 42

#======================================================================================

micro2 <- micro1 %>% 
  select(labid, pid, st1_serotype, st2_serotype, st3_serotype) %>%
  rename("st1" = st1_serotype, "st2" = st2_serotype, "st3" = st3_serotype) %>% filter(st1 !="")

micro2$st1[micro2$st1 == "NT"] <- "NVT"
micro2$st2[micro2$st2 == "NT"] <- "NVT"
micro2$st3[micro2$st3 == "NT"] <- "NVT"

micro2 <- micro2 %>% mutate(mc = if_else(st1 == "VT" & st2 == "NVT", 1L,
                                       if_else(st1 == "NVT" & st2 == "VT", 1L,
                                               if_else(st1 == "VT" & st3 == "NVT", 1L,
                                                       if_else(st1 == "NVT" & st3 == "VT", 1L,
                                                               if_else(st2 == "VT" & st3 == "NVT", 1L,
                                                                       if_else(st2 == "NVT" & st3 == "VT", 1L, 0L))))))
)

#merge with the main dataset
micro2 <- micro2 %>% select(labid, pid, mc) %>% filter(mc == 1)
pcvpa.mod <- left_join(pcvpa.mod, micro2)

#recode the outcome in the main dataset baswed on multiple serotype
pcvpa.mod$nvtcarr2 <- pcvpa.mod$nvtcarr
pcvpa.mod$nvtcarr2[pcvpa.mod$mc == 1] <- 1

pcvpa.mod$vtcarr2 <- pcvpa.mod$vtcarr
pcvpa.mod$vtcarr2[pcvpa.mod$mc == 1] <- 1

pcvpa.mod <- pcvpa.mod %>% select(pid, labid, nvtcarr, vtcarr, nvtcarr1, vtcarr1, nvtcarr2, vtcarr2, mc, everything())

#=======================================================================================

#subset for a correct dataset
crude = pcvpa.mod %>% select(vtcarr, vtcarr2, age, year)

#create age group
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24y",
                                   if_else(age >24 & age <=29, "25-29y",
                                           if_else(age >29 & age <=34, "30-34y",
                                                   if_else(age >34 & age <=40, "35-40y", NA_character_))))))

#fit model to obtain predictions of prevalence
model_crude = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$foi_vt <- model_crude$fitted.values

model_crude = gam(vtcarr2 ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod)
crude$foi_vt2 <- model_crude$fitted.values

#=======================================================================================

#plot prevalence predictions for single and multiple VT carriage detection
ggsave(here("output", "FigS2_multiple_carriage.png"),
       
plot = rbind(
  crude %>% filter(vtcarr == 1) %>% group_by(year) %>% summarise(foi = mean(foi_vt)) %>% mutate(detection = "Latex", catg = "B, VT(+st3) carriage across time (years)") %>% rename("varx" = year),
  crude %>% filter(vtcarr2 == 1) %>% group_by(year) %>% summarise(foi = mean(foi_vt2)) %>% mutate(detection = "Microarray", catg = "B, VT(+st3) carriage across time (years)") %>% rename("varx" = year),
  crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% summarise(foi = mean(foi_vt)) %>% mutate(detection = "Latex", catg = "A, VT(+st3) carriage across age group (years)") %>% rename("varx" = agegp),
  crude %>% filter(vtcarr2 == 1) %>% group_by(agegp) %>% summarise(foi = mean(foi_vt2)) %>% mutate(detection = "Microarray", catg = "A, VT(+st3) carriage across age group (years)") %>% rename("varx" = agegp)
  ) %>% 

 ggplot() +
  geom_point(aes(x = varx, y = foi, color = detection, group = detection), size = 3, shape = 11) +
  geom_line(aes(x = varx, y = foi, color = detection, group = detection), size = 1) +
  facet_grid(.~catg, scales = "free") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16), strip.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks = seq(0, 1, 0.05), labels = scales::percent_format(accuracy = 1)) +  
  labs(title = "", x = "", y = "Carriage prevalence") +
  guides(color=guide_legend(title=""), fill = FALSE) +
  theme(axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 10, face = "bold")) +
  theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) + 
  theme(legend.position = c(0.40, 0.85)) +
  guides(color = guide_legend(title = "Detection method")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)),

width = 10, height = 4, unit="in", dpi = 200)

#=======================================================================================

