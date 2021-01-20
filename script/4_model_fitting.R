#Written by Deus Thindwa
#Pneumococcal carriage prevalence in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#model fitting

#==========================================================================

#non-VT carriage: fitting 6 distinct models based on parametric, semi-parametric, non-parametric and non-proportion hazard assumptions

model1 = gam(nvtcarr ~ age + surv + sex + artdur + artreg + ctx + cd4cnt + nochild5 + sescat, 
             family = binomial(link = "cloglog"), 
             na.action=na.exclude,
             data = pcvpa.mod,
             control = gam.control(trace = TRUE))

model2 = gam(nvtcarr ~ te(age) + surv + sex + artdur + artreg + ctx + cd4cnt + nochild5 + sescat, 
             family = binomial(link = "cloglog"), 
             na.action=na.exclude,
             data = pcvpa.mod,
             control = gam.control(trace = TRUE))

model3 = gam(nvtcarr ~ te(age) + te(surv) + sex + artdur + artreg + ctx + cd4cnt + nochild5 + sescat, 
             family = binomial(link = "cloglog"), 
             na.action=na.exclude,
             data = pcvpa.mod,
             control = gam.control(trace = TRUE))

model4 = gam(nvtcarr ~ ti(age, surv) + sex + artdur + artreg + ctx + cd4cnt + nochild5 + sescat, 
             family = binomial(link = "cloglog"), 
             na.action=na.exclude,
             data = pcvpa.mod,
             control = gam.control(trace = TRUE))

#model comparisons using AICc
AICc(model1, model2, model3, model4)
BIC(model1, model2, model3, model4)

#==========================================================================

#VT carriage: fitting 6 distinct models based on parametric, semi-parametric, non-parametric and non-proportion hazard assumptions

model5 = gam(vtcarr ~ age + surv + sex + artdur + artreg + ctx + cd4cnt + nochild5 + sescat, 
             family = binomial(link = "cloglog"), 
             na.action=na.exclude,
             data = pcvpa.mod,
             control = gam.control(trace = TRUE))

model6 = gam(vtcarr ~ te(age) + surv + sex + artdur + artreg + ctx + cd4cnt + nochild5 + sescat, 
             family = binomial(link = "cloglog"), 
             na.action=na.exclude,
             data = pcvpa.mod,
             control = gam.control(trace = TRUE))

model7 = gam(vtcarr ~ te(age) + te(surv) + sex + artdur + artreg + ctx + cd4cnt + nochild5 + sescat, 
             family = binomial(link = "cloglog"), 
             na.action=na.exclude,
             data = pcvpa.mod,
             control = gam.control(trace = TRUE))

model8 = gam(vtcarr ~ ti(age, surv) + sex + artdur + artreg + ctx + cd4cnt + nochild5 + sescat, 
              family = binomial(link = "cloglog"), 
              na.action=na.exclude,
              data = pcvpa.mod,
              control = gam.control(trace = TRUE))

#model comparisons using AICc
AICc(model5, model6, model7, model8)
BIC(model5, model6, model7, model8)
