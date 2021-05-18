#Written by Deus Thindwa
#Pneumococcal carriage prevalence & FOI in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#======================================================================================

micro <- micro %>% mutate(arrayplus = str_detect(arraysero, '[+]'),
         arraystar = str_detect(arraysero, '[*]'),
         multicarr = if_else(arrayplus == TRUE & arraystar == FALSE, 1, 0)) %>%
         select(labid, multicarr) 

pcvpa.mod <- left_join(pcvpa.mod, micro)

#=======================================================================================

pcvpa.mod %>% filter(!is.na(multicarr)) %>% group_by(nvtcarr, age, multicarr) %>% tally() %>%
ggplot(aes(x = age, y = n/sum(n), color = as.factor(multicarr))) +
  geom_line()

pcvpa.mod %>% filter(!is.na(multicarr)) %>% group_by(nvtcarr, surv, multicarr) %>% tally() %>%
  ggplot(aes(x = surv, y = n/sum(n), color = as.factor(multicarr))) +
  geom_line()

#=======================================================================================

#turn on warnings
options(warn = defaultW)

ggsave(here("output", "Fig4_artdur_prev_foi.tiff"),
       plot = (A | B | C | D | E | F) / (G | H | I | J | K | L),
       width = 16, height = 6, unit="in", dpi = 200)
