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
model_crude = gam(vtcarr ~ s(age, bs="ps") + s(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = pcvpa.mod, na.action = na.exclude)
crude$fit = predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
crude$fit_lci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit - (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))
crude$fit_uci = model_crude$family$linkinv(predict.gam(model_crude, type = "link", se.fit = TRUE)$fit + (2 * predict.gam(model_crude, type = "link", se.fit = TRUE)$se.fit))





