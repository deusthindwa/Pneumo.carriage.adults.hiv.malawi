ggplot(data = pcvpa.des) +
  #geom_density(aes(x = age)) + 
  geom_density(aes(x = artdur))

#duration on ART distribution
A <- pcvpa.des %>%  
  ggplot() + 
  geom_histogram(aes(x = age), bins = 23, color = "black", fill =  brocolors("crayons")["Goldenrod"]) +
  geom_density(aes(x = age, y = ..density../0.001), alpha = 0.3, size = 1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(18, 40, 2)) +
  scale_y_continuous("Number of participants", sec.axis = sec_axis(~. * 0.001, name = "Probability density"), limits = c(0, 170)) + 
  labs(title = "A", x = "Participant age (years)") +
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 0, vjust = 0.5, hjust = 0.3), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "none")

#duration on ART distribution
B <- pcvpa.des %>%  
  ggplot() + 
  geom_histogram(aes(x = artdur), bins = 23, color = "black", fill =  brocolors("crayons")["Goldenrod"]) +
  geom_density(aes(x = artdur, y = ..density../0.001), alpha = 0.3, size = 1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 18, 2)) +
  scale_y_continuous("Number of participants", sec.axis = sec_axis(~. * 0.001, name = "Probability density"), limits = c(0, 400)) + 
  labs(title = "B", x = "Duration on ART (years)") +
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 0, vjust = 0.5, hjust = 0.3), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "none")

#duration on ART by age groups 
C <- pcvpa.des %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24y",
                                   if_else(age >24 & age <=29, "25-29y",
                                           if_else(age >29 & age <=34, "30-34y",
                                                   if_else(age >34 & age <=40, "35-40y", NA_character_)))))) %>%
  ggplot(aes(x = agegp, y = artdur, fill = agegp)) + 
  geom_boxplot(notch = TRUE, color = "black", size = 1) +
  theme_bw() +
  labs(title = "C", x = "Age group (years)", y = "Duration on ART (years)") +
  #scale_y_continuous(breaks = seq(0, 1400, 200)) +
  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "none")


#combined plots
ggsave(here("output", "FigS6_age_ART_duration.png"),
       plot = (A | B | C),
       width = 13, height = 5, unit="in", dpi = 300)
