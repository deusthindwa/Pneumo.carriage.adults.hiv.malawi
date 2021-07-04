#imputing 3 variables: ART duration, CD4 count, SES
summary(pcvpa.mod)

pcvpa.modi <- aregImpute(~ age + surv + seas + sex + artdur + cd4cnt + nochild5 + sescat, data = pcvpa.mod, n.impute = 100)

pcvpa.modi$imputed$artdur
pcvpa.modi$imputed$cd4cnt
pcvpa.modi$imputed$sescat
