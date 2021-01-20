#Written by Deus Thindwa
#Pneumococcal carriage prevalence in HIV-infected adults in PCV era
#Generalized additive model.
#14/12/2020 - 30/06/2021

#participants characteristics (Table 1)


#categorical variables
pcvpa.render <- function(x, name,...){
  if(!is.numeric(x)) 
    return(render.categorical.default(na.omit(x)))
}

table1(~sex + artreg + ctx + sescat | serogroup, data = pcvpa.des, topclass = "Rtable1-grid Rtable1-shade Rtable1-times", render = pcvpa.render)
