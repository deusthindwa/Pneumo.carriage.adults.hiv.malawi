slope = function(data, indices){
  data = data[indices,]
  vtcarr = data[,1]
  age = data[,2]
  year = data[,3]
  sex = data[,4]
  nochild5 = data[,5]
  seas = data[,6]
  model_crude = gam(vtcarr ~ te(age, bs="ps") + te(year, bs="ps") + sex + nochild5 + seas, family = binomial(link = "cloglog"))
  predict.gam(model_crude, type = "response", se.fit = TRUE)$fit
}

data = crude
set.seed(1988)
slopeb = as_tibble(t(boot(data, slope, R = 1000, parallel = "multicore", ncpus = 3)$t))
slopec <- slopeb %>% 
  rowwise() %>% 
  mutate(fit = mean(c_across(c(V1:V1000)), na.rm = TRUE), 
         se = (sd(c_across(c(V1:V1000)), na.rm = TRUE))/length(c_across(c(V1:V1000))),
         fit_lci = fit - qt(1-(0.05/2), length(c_across(c(V1:V1000)))-1)*se,
         fit_uci = fit + qt(1-(0.05/2), length(c_across(c(V1:V1000)))-1)*se) %>%
  select(fit, fit_lci, fit_uci)

crude <- cbind(crude, slopec)