
func_PhenotypicAge = function(dataset){
  attach(dataset)
  xb= -19.9067-0.0336*(Albumin*10)+0.0095*(Creatinine*88.4)+0.1953*(Glucose*0.0555) + # all reqd: Albumin - g/L, Creat - umol/L, Gluc - mmol/L
      0.0954*log(hsCRP*0.1) -0.012*Lymp_pct +
      0.0268*MCV +0.3306*RDW  + #units RDW % , MCV fL (1e-15*L) reqd.
      0.00188*AlkPhos+0.0554*WBC+0.0804*Age #AlkPhos - U/L, WBC - 1000cells/uL
  
  detach(dataset)
  
  M = 1-exp(-1.51714*exp(xb)/0.0076927)
  pheno_age=141.50225+log(-0.00553*log(1-M))/0.09165 # or 0.090165?
  
  # Note: If not working, might have to determine gamma (Gompertz coefs) based on dataset instead
  # of taking it as 0.0076927... perhaps other coefs pertaining to xb covariates too..
  
  return(pheno_age)
}
