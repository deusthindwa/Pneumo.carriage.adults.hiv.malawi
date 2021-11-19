#Written by Deus Thindwa
#Pneumococcal carriage prevalence in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 20/11/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#=================================================================================

#VT CARRIAGE WITH SEROTYPE 3

#subset for a dataset to store model estimates
crude = pcvpa.mod %>% select(vtcarr, age, year)

#fit model to individual trajectories & obtain predictions and 95%CI
model_crude = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod, na.action = na.exclude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))

#----------------------------------------------------------------------------------

str(model_crude$smooth)




#create age groups
crude <- crude %>% 
  mutate(agegp = as.factor(if_else(age <=24, "18-24",
                                   if_else(age >24 & age <=29, "25-29",
                                           if_else(age >29 & age <=34, "30-34",
                                                   if_else(age >34 & age <=40, "35-40", NA_character_))))))

#get predicted mean prevalence for each age group 
crude_age1 <- left_join(left_join(crude %>% group_by(agegp) %>% tally() %>% rename(Tot = n), 
                                  crude %>% filter(vtcarr == 1) %>% group_by(agegp) %>% tally() %>% rename(Pos = n)), 
                        crude %>% filter(vtcarr != 0) %>% group_by(agegp) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
                          ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:4], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[5:8], status = "VT(+st3) carriage")

#----------------------------------------------------------------------------------

#get predicted mean prevalence for each year 
crude_year1 <- left_join(left_join(crude %>% group_by(year) %>% tally() %>% rename(Tot = n), 
                                   crude %>% filter(vtcarr == 1) %>% group_by(year) %>% tally() %>% rename(Pos = n)), 
                         crude %>% filter(vtcarr != 0) %>% group_by(year) %>% summarise(fit = mean(fit), fit_lci = mean(fit_lci), fit_uci = mean(fit_uci)) %>%
                           ungroup()) %>% mutate(obs = Pos/Tot, obs_lci = exactci(Pos, Tot, 0.95)$conf.int[1:5], obs_uci = exactci(Pos, Tot, 0.95)$conf.int[6:10], status = "VT(+st3) carriage")

# Table 2.3, page 20
MI <- matrix(c(crude_age1$fit, 104, 10845, 10933), nrow = 2)
dimnames(MI) <- list("Group" = c("Placebo","Aspirin"), "MI" = c("Yes","No"))
MI
##          MI
## Group     Yes    No
##   Placebo 189 10845
##   Aspirin 104 10933



