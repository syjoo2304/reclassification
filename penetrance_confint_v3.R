library(binom)

##from https://www.cureffi.org/2016/10/19/estimation-of-penetrance/
penetrance = function(af_case, af_control, baseline_risk) {
  calculated_penetrance = af_case * baseline_risk / af_control
  estimated_penetrance = pmin(1,pmax(0,calculated_penetrance)) # trim to [0,1] support
  return (estimated_penetrance)
}
##from https://www.cureffi.org/2016/10/19/estimation-of-penetrance/
penetrance_confint = function (ac_case, n_case, ac_control, n_control, baseline_risk) {
  # for a genotypic model, use 1*n_case; for allelic, use 2*n_case
  case_confint = binom.confint(x=ac_case,n=2*n_case,method='wilson')
  control_confint = binom.confint(x=ac_control,n=2*n_control,method='wilson')
  lower_bound = penetrance(case_confint$lower,control_confint$upper,baseline_risk)
  best_estimate = penetrance(case_confint$mean,control_confint$mean,baseline_risk)
  upper_bound = penetrance(case_confint$upper,control_confint$lower,baseline_risk)
  return ( c(lower_bound, best_estimate, upper_bound) )
}

summarize_penetrance_variants <- function(penetrance_list_named) {
  # 입력: list("variant1" = list(x, y, z), "variant2" = list(...), ...)
  summary_list <- lapply(names(penetrance_list_named), function(variant_name) {
    variant_data <- penetrance_list_named[[variant_name]]
    
    lowers <- sapply(variant_data, function(x) x[1])
    estimates <- sapply(variant_data, function(x) x[2])
    uppers <- sapply(variant_data, function(x) x[3])
    
    data.frame(
      variant = variant_name,
      lower = mean(lowers),
      estimate = mean(estimates),
      upper = mean(uppers)
    )
  })
  
  summary_df <- do.call(rbind, summary_list)
  return(summary_df)
}

#MYO7A R675C
x=penetrance_confint(4,1302,8,5144,0.033) 
y=penetrance_confint(4,1302,5,8548,0.033) 
z=penetrance_confint(4,1302,32,108604,0.033) 

#TJP2 A122T
x1=penetrance_confint(3,1302,73,10562,0.033) 
y1=penetrance_confint(3,1302,54,8548,0.033) 
z1=penetrance_confint(3,1302,1306,108604,0.033) 

#TJP2 p.T1188A
x2=penetrance_confint(4,1302,32,10478,0.033)
y2=penetrance_confint(4,1302,29,8548,0.033) 
z2=penetrance_confint(4,1302,706,108604,0.033) 

#P2RX2 p.D273Y
x3=penetrance_confint(4,1302,20,10468,0.033)
y3=penetrance_confint(4,1302,27,8548,0.033) 
z3=penetrance_confint(4,1302,132,108604,0.033) 

#MYH14 p.R1726W
x4=penetrance_confint(2,1302,19,10006,0.033)
y4=penetrance_confint(2,1302,19,8548,0.033) 
z4=penetrance_confint(2,1302,62,108604,0.033) 

#MYH9 p.R1730C
x5=penetrance_confint(3,1302,16,9944,0.033)
y5=penetrance_confint(3,1302,20,8548,0.033) 
z5=penetrance_confint(3,1302,62,108602,0.033) 

#KCNQ4 p.F182L
x6=penetrance_confint(2,1302,14,9722,0.033)
y6=penetrance_confint(2,1302,19,8548,0.033) 
z6=penetrance_confint(2,1302,423,108604,0.033) 

#MYH9 p.R802W
x7=penetrance_confint(1,1302,10,10496,0.033)
y7=penetrance_confint(1,1302,8,8548,0.033) 
z7=penetrance_confint(1,1302,43,108604,0.033) 

#DIAPH1 p.I700N
x8=penetrance_confint(3,1302,20,6528,0.033)
y8=penetrance_confint(3,1302,11,8548,0.033) 
z8=penetrance_confint(3,1302,451,106882,0.033) 

#MYH14 p.G697S
x9=penetrance_confint(1,1302,18,10388,0.033)
y9=penetrance_confint(1,1302,18,8548,0.033) 
z9=penetrance_confint(1,1302,368,108604,0.033) 

#GSDME p.R261*
x10=penetrance_confint(2,1302,13,10550,0.033)
y10=penetrance_confint(2,1302,11,8548,0.033) 
z10=penetrance_confint(2,1302,136,108604,0.033) 

#GRHL2 p.Q445R
x11=penetrance_confint(3,1302,13,10586,0.033)
y11=penetrance_confint(3,1302,3,8548,0.033) 
z11=penetrance_confint(3,1302,285,108602,0.033) 

#MYO6 p.R1166L
x12=penetrance_confint(2,1302,8,10546,0.033)
y12=penetrance_confint(2,1302,4,8548,0.033) 
z12=penetrance_confint(2,1302,66,108604,0.033) 

#MYH14 p.K898Q
x13=penetrance_confint(3,1302,4,9788,0.033)
y13=penetrance_confint(3,1302,5,8548,0.033) 
z13=penetrance_confint(3,1302,82,108604,0.033) 

#
df <- summarize_penetrance_variants(list("MYO7A p.R675C"=list(x,y,z),"TJP2 p.A122T"=list(x1,y1,z1),"TJP2 p.T1188A"=list(x2,y2,z2),
"P2RX2 p.D273Y"=list(x3,y3,z3),"MYH14 p.R1726W"=list(x4,y4,z4),"MYH9 p.R1730C"=list(x5,y5,z5),"KCNQ4 p.F182L"=list(x6,y6,z6),
"MYH9 p.R802W"=list(x7,y7,z7),"DIAPH1 p.I700N"=list(x8,y8,z8),"MYH14 p.G697S"=list(x9,y9,z9),"GSDME p.R261*"=list(x10,y10,z10),
"GRHL2 p.Q445R"=list(x11,y11,z11),"MYO6 p.R1166L"=list(x12,y12,z12),"MYH14 p.K898Q"=list(x13,y13,z13)))

write.table(df,"reclassification_combine.txt",quote=F,sep='\t')
