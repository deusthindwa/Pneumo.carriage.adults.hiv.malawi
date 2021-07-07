#imputing 3 variables: ART duration, CD4 count, SES
summary(pcvpa.mod)
set.seed(1988)
pcvpa.mod1 <- pcvpa.mod %>% select(vtcarr, vtcarr1, nvtcarr, nvtcarr1, year, age, seas, sex, artdur, cd4cnt, nochild5, sescat) %>% mutate(artdur = as_factor(artdur), cd4cnt = as_factor(cd4cnt), sescat = as_factor(sescat))
pcvpa.mod2 <- missForest(pcvpa.mod1, ntree = 100)
pcvpa.mod2$OOBerror
pcvpa.err <- mixError(pcvpa.mod2$ximp, pcvpa.mod1, pcvpa.mod)

pcvpa.mod2 <- pcvpa.mod2$ximp




