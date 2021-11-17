#Written by Deus Thindwa
#Pneumococcal carriage prevalence in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#=======================================================================

#descriptive analysis plots

#serotype frequency
A <- ggplot(data = filter(pcvpa.des, serotype != "NVT" & !is.na(serotype)), mapping = aes(x = factor(serotype,levels(factor(serotype))[c(1, 7, 8, 9, 10, 11, 12, 13, 2, 3, 4, 5, 6)]), fill = serotype)) + 
  geom_bar(color = "black") + 
  theme_bw() + 
  labs(title = "A", x = "Pneumococcal congugate vaccine serotype (VT)", y = "Frequency") +
  scale_fill_grey(start = 0.9, end = 0.2) +
  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold",size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y=element_text(face = "bold", size = 12)) +
  theme(legend.position = "none")

#serogroup frequency
X <- ggplot(data = pcvpa.des, mapping = aes(x = serogroup, fill = serogroup)) + 
  geom_bar(color = "black") + 
  theme_bw() + 
  labs(title = "", x = "Serogroup", y = "Frequency") +
  scale_y_continuous(breaks = seq(0, 1200, 300)) +
  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "none")

#surveys summary
B <- pcvpa.des %>% group_by(age, surv) %>% tally() %>%
  ggplot(mapping = aes(x = age, y = surv, fill = n)) +
  geom_point(aes(size = n, color = n), shape = 21) +
  theme_bw() +
  labs(title = "B", x = "Age (years)", y = "Survey number") +
  scale_fill_continuous(type = "viridis") +
  scale_color_continuous(type = "viridis") +
  scale_x_continuous(breaks = seq(18, 40, 2)) +
  scale_y_continuous(breaks = seq(1, 8, 1)) +
  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  guides(fill = guide_legend(title = ""), size = guide_legend(title = "")) + 
  theme(legend.position = "bottom")

#pneumococcal serogroup and number of children within household
C <-  pcvpa.des %>% group_by(nochild5, serogroup) %>% tally() %>%
  ggplot(mapping = aes(x = nochild5, y = serogroup, fill = n)) +
  geom_point(aes(size = n, color = n), shape = 21) +
  theme_bw() +
  labs(title = "C", x = "Number of children in the house", y = "Serogroup in HIV-infected adults") +
  scale_fill_continuous(type = "viridis") +
  scale_color_continuous(type = "viridis") +  
  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size=12), axis.title.y = element_text(face = "bold", size = 12)) +
  guides(fill = guide_legend(title = ""), size = guide_legend(title = "")) + 
  theme(legend.position = "bottom")

#pneumococcal serogroup and age
D <- ggplot(data = pcvpa.des, mapping = aes(x = serogroup, y = age, fill=serogroup)) + 
  geom_boxplot(notch = TRUE, color = "black", size = 1) + 
  theme_bw() +
  labs(title = "D", x = "Serogroup", "Age (years)") +
  scale_y_continuous(breaks = seq(0, 40, 4)) +
  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "none")

#pneumococcal serogroup and ART duration
E <- ggplot(data = pcvpa.des, mapping = aes(x = serogroup, y = artdur, fill = serogroup)) + 
  geom_boxplot(notch = TRUE, color = "black", size = 1) + 
  theme_bw() +
  labs(title = "E", x = "Serogroup", y = "ART duration (years)") +
  scale_y_continuous(breaks = seq(0, 18, 3)) +
  theme(axis.text.x = element_text(face = "bold", size = 12),axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "none")

#pneumococcal serogroup and CD4+ cell count
F <- ggplot(data = pcvpa.des, mapping = aes(x = serogroup, y = cd4cnt, fill = serogroup)) + 
  geom_boxplot(notch = TRUE, color = "black", size = 1) +
  theme_bw() +
  labs(title = "F", x = "Serogroup", y = "CD4+ count (cells/mm3)") +
  scale_y_continuous(breaks = seq(0, 1400, 200)) +
  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "none")

#sex
filter(pcvpa.des, !is.na(sex) & serogroup !="None") %>% group_by(sex, serogroup) %>% tally() %>% mutate(sex_p = n/sum(n))
G <- filter(pcvpa.des, !is.na(sex)) %>% group_by(serogroup, sex) %>% tally() %>% mutate(sex_p = n/sum(n)) %>% 
  ggplot(mapping = aes(x = serogroup, y = sex_p, fill = sex)) + 
  geom_bar(stat = "identity", color = "black", size = 1) +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme_bw() +
  labs(title = "G", x = "", y = "Proportion in sex groups") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(face = "bold", size = 12), legend.box.spacing = unit(0, 'cm'))

#social economic status
filter(pcvpa.des, !is.na(sescat) & serogroup !="None") %>% group_by(sescat, serogroup) %>% tally() %>% mutate(ses_p = n/sum(n))
H <- filter(pcvpa.des, !is.na(sescat)) %>% group_by(serogroup, sescat) %>% tally() %>% mutate(ses_p = n/sum(n)) %>%
  ggplot(mapping = aes(x = serogroup, y = ses_p, fill = factor(sescat,levels(factor(sescat))[c(2, 3, 1)]))) + 
  geom_bar(stat = "identity", color = "black", size = 1) +  
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme_bw() +
  labs(title = "H", x = "", y = "Proportion in social economic status") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(face = "bold", size = 12), legend.box.spacing = unit(0, 'cm'))

#ART regimen
filter(pcvpa.des, !is.na(artreg) & serogroup !="None") %>% group_by(artreg, serogroup) %>% tally() %>% mutate(artreg_p = n/sum(n))
I <- filter(pcvpa.des, !is.na(artreg)) %>% group_by(serogroup, artreg) %>% tally() %>% mutate(artreg_p = n/sum(n)) %>%
  ggplot(mapping = aes(x = serogroup, y = artreg_p, fill = artreg)) + 
  geom_bar(stat = "identity", color = "black", size = 1) +  
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme_bw() +
  labs(title = "I", x = "", y = "Proportion in ART regimen") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(face = "bold", size = 12), legend.box.spacing = unit(0, 'cm'))

#cotrimoxazole
filter(pcvpa.des, !is.na(ctx) & serogroup !="None") %>% group_by(ctx, serogroup) %>% tally() %>% mutate(ctx_p = n/sum(n))
J <- filter(pcvpa.des, !is.na(ctx)) %>% group_by(serogroup, ctx) %>% tally() %>% mutate(ctx_p = n/sum(n)) %>%
  ggplot(mapping = aes(x = serogroup, y = ctx_p, fill = ctx)) + 
  geom_bar(stat = "identity", color = "black", size = 1) +  
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme_bw() +
  labs(title = "J", x = "", y = "Proportion in cotrimoxazole status") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +  theme(axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(face = "bold", size = 12), legend.box.spacing = unit(0, 'cm'))

#combined plots
ggsave(here("output", "Fig1_descriptive.png"),
       plot = (A + inset_element(X, right = 0.7, bottom = 0.4, left = 0.3, top = 0.9) + B + C) / (D | E | F | G | H | I | J),
       width = 22, height = 13, unit="in", dpi = 250)

#=======================================================================

#continuous variables notched box plot statistics

#age
filter(pcvpa.des, !is.na(age)) %>% group_by(serogroup) %>% summarise(median(age), quantile(age, 1/4), quantile(age, 3/4))
filter(pcvpa.des, serogroup == "NVT" | serogroup == "VT" & !is.na(age)) %>% summarise(median(age), quantile(age, 1/4), quantile(age, 3/4))
filter(pcvpa.des, !is.na(age)) %>% summarise(median(age), quantile(age, 1/4), quantile(age, 3/4))
median_se(subset(pcvpa.des, !is.na(age) & serogroup == "VT")$age)
median_se(subset(pcvpa.des, !is.na(age) & serogroup == "NVT")$age)

#ART duration
filter(pcvpa.des, !is.na(artdur)) %>% group_by(serogroup) %>% summarise(median(artdur), quantile(artdur, 1/4), quantile(artdur, 3/4))
filter(pcvpa.des, (serogroup == "NVT" | serogroup == "VT") & !is.na(artdur)) %>% summarise(median(artdur), quantile(artdur, 1/4), quantile(artdur, 3/4))
filter(pcvpa.des, !is.na(artdur)) %>% summarise(median(artdur), quantile(artdur, 1/4), quantile(artdur, 3/4))
median_se(subset(pcvpa.des, !is.na(artdur) & serogroup == "VT")$artdur)
median_se(subset(pcvpa.des, !is.na(artdur) & serogroup == "NVT")$artdur)

#CD4+ count
filter(pcvpa.des, !is.na(cd4cnt)) %>% group_by(serogroup) %>% summarise(median(cd4cnt), quantile(cd4cnt, 1/4), quantile(cd4cnt, 3/4))
filter(pcvpa.des, (serogroup == "NVT" | serogroup == "VT") & !is.na(cd4cnt)) %>% summarise(median(cd4cnt), quantile(cd4cnt, 1/4), quantile(cd4cnt, 3/4))
filter(pcvpa.des, !is.na(cd4cnt)) %>% summarise(median(cd4cnt), quantile(cd4cnt, 1/4), quantile(cd4cnt, 3/4))
median_se(subset(pcvpa.des, !is.na(cd4cnt) & serogroup == "VT")$cd4cnt)
median_se(subset(pcvpa.des, !is.na(cd4cnt) & serogroup == "NVT")$cd4cnt)

