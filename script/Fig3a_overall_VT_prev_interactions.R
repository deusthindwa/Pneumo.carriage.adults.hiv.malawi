#Written by Deus Thindwa
#Pneumococcal carriage prevalence & acquisition in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 17/12/2021

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1)

#===================================================================================

#add age group to dataset
interaction <-
  pcvpa.mod %>%
  mutate(agegp = if_else(age <=24, 1L,
                       if_else(age >24 & age <=29, 2L,
                               if_else(age >29 & age <=34, 3L,
                                       if_else(age >34 & age <=40, 4L, NA_integer_)))))

#=============

model_a = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(agegp, bs="ps", by = sex) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(agegp, bs="ps", by = sescat) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(agegp, bs="ps", by = artdur) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(agegp, bs="ps", by = nochild5) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

#=============

model_a = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(year, bs="ps", by = sex) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(year, bs="ps", by = sescat) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(year, bs="ps", by = artdur) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(carr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(year, bs="ps", by = nochild5) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

#=============

model_a = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(agegp, bs="ps", by = sex) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(agegp, bs="ps", by = sescat) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(agegp, bs="ps", by = artdur) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(agegp, bs="ps", by = nochild5) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

#=============

model_a = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(year, bs="ps", by = sex) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(year, bs="ps", by = sescat) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(year, bs="ps", by = artdur) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

model_a = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
model_b = gam(vtcarr ~ te(agegp, bs="ps") + te(year, bs="ps") + ti(year, bs="ps", by = nochild5) + seas + sex + artdur + nochild5 + sescat, family = binomial(link = "cloglog"), data = interaction, na.action = na.exclude)
AICc(model_a, model_b)

#=============




