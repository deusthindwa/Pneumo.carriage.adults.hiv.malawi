#Written by Deus Thindwa
#Pneumococcal carriage prevalence in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#age-survey prevalence estimations

#==========================================================================

#plot age and time-dependent overall, NVT, and VT carriage prevalences
png(filename="output/Sfig1_overall.tiff", width = 2000, height = 900, res = 200)
par(mfrow=c(1,2))
par(mar=c(5,5,2,2))
vis.gam(model4, view = c("age", "surv"), n.grid=100, contour.col = "black", type = "response", plot.type = "contour", color = 'heat', main = "NVT prevalence", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, contour.col = "black", type = "response", plot.type = "contour", color = "heat", main = "VT prevalence", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

#plot age and time-dependent NVT, and VT carriage prevalences by sex
png(filename="output/Sfig2_sex.tiff", width = 4500, height = 1500, res = 400)
par(mfrow=c(1,4))
par(mar=c(5,5,2,2))
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(sex=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, females", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(sex=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, males", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(sex=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, females", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(sex=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, males", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

#plot age and time-dependent NVT, and VT carriage prevalences by ART duration
png(filename="output/Sfig3_artdur.tiff", width = 3000, height = 2000, res = 300)
par(mfrow=c(2,3))
par(mar=c(5,5,2,2))
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(artdur=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, ART <2y", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(artdur=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, ART 2-6y", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(artdur=2), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, ART 7+y", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(artdur=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, ART <2y", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(artdur=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, ART 2-6y", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(artdur=2), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, ART 7+y", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

#plot age and time-dependent NVT, and VT carriage prevalences by ART regimen
png(filename="output/Sfig4_artreg.tiff", width = 4500, height = 1500, res = 400)
par(mfrow=c(1,4))
par(mar=c(5,5,2,2))
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(artreg=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, first line ART", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(artreg=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, second line ART", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(artreg=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, first line ART", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(artreg=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, second line ART", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

#plot age and time-dependent NVT and VT carriage prevalences by cotrimoxazole
png(filename="output/Sfig5_cotri.tiff", width = 4500, height = 1500, res = 400)
par(mfrow=c(1,4))
par(mar=c(5,5,2,2))
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(ctx=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, cotrimoxazole", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(ctx=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, no cotrimoxazole", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(ctx=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, cotrimoxazole", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(ctx=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, no cotrimoxazole", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

#plot age and time-dependent NVT and VT carriage prevalences by CD4+ count
png(filename="output/Sfig6_cd4.tiff", width = 4500, height = 1500, res = 400)
par(mfrow=c(1,4))
par(mar=c(5,5,2,2))
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(cd4cnt=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, low CD4+ count", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(cd4cnt=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, high CD4+ count", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(cd4cnt=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, low CD4+ count", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(cd4cnt=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, high CD4+ count", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

#plot age and time-dependent NVT and VT carriage prevalences by cohabiting children
png(filename="output/Sfig7_child.tiff", width = 4500, height = 1500, res = 400)
par(mfrow=c(1,4))
par(mar=c(5,5,2,2))
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(nochild5=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, households without <5y children", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(nochild5=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, households with <5y children", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(nochild5=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, households without <5y children", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(nochild5=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, households with <5y children", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

#plot age and time-dependent NVT and VT carriage prevalences by social economic status
png(filename="output/Sfig8_ses.tiff", width = 3000, height = 2000, res = 300)
par(mfrow=c(2,3))
par(mar=c(5,5,2,2))
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(sescat=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, low SES", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(sescat=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, medium SES", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model4, view = c("age", "surv"), n.grid=100, cond = list(sescat=2), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "NVT, high SES", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(sescat=0), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, low SES", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(sescat=1), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, medium SES", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
vis.gam(model8, view = c("age", "surv"), n.grid=100, cond = list(sescat=2), contour.col = "black", type = "response", plot.type = "contour", color = "heat", too.far = 0, main = "VT, high SES", xlab = "Age,y", ylab = "Survey number", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

