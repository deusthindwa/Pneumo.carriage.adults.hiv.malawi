#Written by Deus Thindwa
#Pneumococcal carriage prevalence in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#Age and time-specific prevalence and force of infection estimations

#==========================================================================

#select appropriate dataset to fit
pcvpa.sex0 <- filter(pcvpa.mod, surv == 1 & sex == 0)
pcvpa.sex1 <- filter(pcvpa.mod, surv == 1 & sex == 1)

model.sex0 = gam(nvtcarr ~ te(age), family = binomial(link = "cloglog"), data = pcvpa.sex0)
model.sex1 = gam(nvtcarr ~ te(age), family = binomial(link = "cloglog"), data = pcvpa.sex1)

#estimate prevalence and 95%CI from optimal model
newDF <- with(pcvpa.sex, data.frame(age,sex))

W <- predict(model.sex, 
             newDF, 
             type = "response", 
             se.fit = TRUE)

pcvpa.sex$fit <- W$fit

pcvpa.sex$fit_der <- iSpline(pcvpa.sex$age, degree = 2, intercept = TRUE, derivs = 1)[,3]


#plot observed & fitted prevalences, and FoI
pcvpa.sex0 <- pcvpa.sex0 %>% group_by(age, nvtcarr) %>% tally() %>% mutate(prev = n/sum(n)) %>% filter(nvtcarr !=0)
pcvpa.sex1 <- pcvpa.sex1 %>% group_by(age, nvtcarr) %>% tally() %>% mutate(prev = n/sum(n)) %>% filter(nvtcarr !=0)
pcvpa.sex0$sex <- "female"; pcvpa.sex1$sex <- "male"
pcvpa.sex <- rbind(pcvpa.sex0, pcvpa.sex1)

pcvpa.mod %>% group_by(age, sex, nvtcarr) %>% tally() %>% mutate(prev = n/sum(n)) %>% filter(nvtcarr !=0) %>%  
ggplot() +
  geom_point(aes(x = age, y = prev, size = n, color=as.factor(sex)), shape = 1) +
  geom_smooth(aes(x = age, y = prev, color=as.factor(sex)), se=FALSE, method = gam, size = 0.5) +
  labs(title = "A", x = "Age in years", y = "Prevalence") +
  theme_bw() +
  coord_cartesian(ylim=c(0,0.6)) +
  scale_y_continuous(breaks = seq(0, 0.6, 0.2)) +
  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "right")

#===============================================================================




