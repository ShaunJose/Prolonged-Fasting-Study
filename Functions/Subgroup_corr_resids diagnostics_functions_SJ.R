# IMPORTANT FUNCTIONS

## Assigns 0 to correlations under threshold for easy visualisation of strong correlations
func_sig.corr_matrices = function(corr_mat, threshold = 0.7){
  corr_mat[abs(corr_mat)<threshold] = 0 # or 'NA' depending on visualisation used.
  return(corr_mat)
}

## Subgrouping function for delta fractional and delta basic datasets
func_delta_frac_subgroups = function(delta_dataset, grp_var_name,
                                     same_dir=FALSE, basic_diff=FALSE,
                                     thresh =0.5 # 50% magnitude change
                                     ) {
  attach(delta_dataset)
  if (basic_diff){
    incr_grp=delta_dataset[grp_var_name >= 0,]
    decr_grp=delta_dataset[grp_var_name < 0,]
    num_incr=dim(incr_grp)[1]
    num_decr=dim(decr_grp)[1]
    output_list=list("increased_grp" = incr_grp, "decreased_grp" = decr_grp,
                     "num_increased_grp" = num_incr,"num_decreased_grp" = num_decr)
  } else {
    if (same_dir) {
      # when biomarkers go in same direction for everyone
      high_grp=delta_dataset[abs(grp_var_name-1) >= frac_thresh,]
      low_grp=delta_dataset[abs(grp_var_name-1) < frac_thresh,] # here the group that doesn't go as high
      num_highr=dim(high_grp)[1]
      num_lowr=dim(low_grp)[1]
      output_list=list("higher_grp" = high_grp, "lower_grp" = low_grp,
                       "num_higher_grp" = num_highr,"num_lower_grp" = num_lowr)
      
    } else {
      incr_grp=delta_dataset[grp_var_name >= 1,]
      decr_grp=delta_dataset[grp_var_name < 1,]
      num_incr=dim(incr_grp)[1]
      num_decr=dim(decr_grp)[1]
      output_list=list("increased_grp" = incr_grp, "decreased_grp" = decr_grp,
                       "num_increased_grp" = num_incr,"num_decreased_grp" = num_decr)
    }
  }
  detach(delta_dataset)
  return(output_list)
}

## Subgrouping function for delta log datasets
func_delta_log_subgroups = function(delta_log_dataset, grp_var_name,
                                     same_dir=FALSE, 
                                     thresh =0.5) {
  attach(delta_log_dataset)
  if (same_dir) {
    # when biomarkers go in same direction for everyone
    high_grp=delta_log_dataset[abs(exp(grp_var_name)-1) >= thresh,]
    low_grp=delta_log_dataset[abs(exp(grp_var_name)-1) < thresh,] # here the group that doesn't go as high
    num_highr=dim(high_grp)[1]
    num_lowr=dim(low_grp)[1]
    output_list=list("higher_grp" = high_grp, "lower_grp" = low_grp,
                     "num_higher_grp" = num_highr,"num_lower_grp" = num_lowr)
    
  } else {
    incr_grp=delta_log_dataset[grp_var_name >= 0,]
    decr_grp=delta_log_dataset[grp_var_name < 0,]
    num_incr=dim(incr_grp)[1]
    num_decr=dim(decr_grp)[1]
    output_list=list("increased_grp" = incr_grp, "decreased_grp" = decr_grp,
                     "num_increased_grp" = num_incr,"num_decreased_grp" = num_decr)
  }
detach(delta_log_dataset)
return(output_list)
}



## Subgrouping function for delta fractional datasets
func_DeltaFracSubgroups = function(delta_dataset, grp_var_name, same_dir=FALSE, thresh =0.5){
  attach(delta_dataset)
  if (same_dir) {
    # when biomarkers go in same direction for everyone
    high_grp=delta_dataset[abs(grp_var_name-1) >= thresh,]
    low_grp=delta_dataset[abs(grp_var_name-1) < thresh,] # here the group that doesn't go as high
    num_highr=dim(high_grp)[1]
    num_lowr=dim(low_grp)[1]
    # Get group names from grp_var_name, to paste it to get into plot legend
    high_name <- paste("High", deparse(substitute(grp_var_name)), sep = " ") ##SJ
    low_name <- paste("Low", deparse(substitute(grp_var_name)), sep = " ") ##SJ
    high_name <- gsub("_.*", "\\1", high_name) ##AM removes _ + all char. after
    low_name <- gsub("_.*", "\\1", low_name) ##AM removes _ + all char. after
    # Code new group var
    df_highr_grp_ID <- tibble::rownames_to_column(high_grp, "ID")
    df_highr_grp_ID <- add_column(df_highr_grp_ID, 
                                  Group = high_name, .after = 1) ##AM high_name does not work ##SJ will work now
    df_highr_grp_ID$Group <- as.character(high_name)
    df_lowr_grp_ID <- tibble::rownames_to_column(low_grp, "ID")
    df_lowr_grp_ID <- add_column(df_lowr_grp_ID, 
                                 Group = low_name, .after = 1)
    df_lowr_grp_ID$Group <- as.character(low_name)
    df_both_grp_ID <- rbind(df_highr_grp_ID, df_lowr_grp_ID)
    df_both_grp_ID[1:2] <- lapply(df_both_grp_ID[1:2], factor)
    output_list=list("higher_grp" = high_grp, "lower_grp" = low_grp,
                     "both_grp" = df_both_grp_ID, 
                     "num_higher_grp" = num_highr,
                     "num_lower_grp" = num_lowr)
    
  } else {
    incr_grp=delta_dataset[grp_var_name >= 1,]
    decr_grp=delta_dataset[grp_var_name < 1,]
    num_incr=dim(incr_grp)[1]
    num_decr=dim(decr_grp)[1]
    # Get group names ### passing column name not working atm
    inc_name <- paste("Incr", deparse(substitute(grp_var_name)), sep = " ")
    dec_name <- paste("Decr", deparse(substitute(grp_var_name)), sep = " ")
    inc_name <- gsub("_.*", "\\1", inc_name) ##AM removes _ + all char. after
    dec_name <- gsub("_.*", "\\1", dec_name) ##AM removes _ + all char. after
    # Code new group var
    df_incr_grp_ID <- tibble::rownames_to_column(incr_grp, "ID")
    df_incr_grp_ID <- add_column(df_incr_grp_ID, 
                                 Group = inc_name, .after = 1)
    # df_highr_grp_ID$Group <- as.character(high_name)
    df_decr_grp_ID <- tibble::rownames_to_column(decr_grp, "ID")
    df_decr_grp_ID <- add_column(df_decr_grp_ID, 
                                 Group = dec_name, .after = 1)
    # df_lowr_grp_ID$Group <- as.character(low_name)
    df_both_grp_ID <- rbind(df_incr_grp_ID, df_decr_grp_ID)
    df_both_grp_ID[1:2] <- lapply(df_both_grp_ID[1:2], factor)
    output_list=list("increased_grp" = incr_grp, "decreased_grp" = decr_grp,
                     "both_grp" = df_both_grp_ID,
                     "num_increased_grp" = num_incr,"num_decreased_grp" = num_decr)
  }
  detach(delta_dataset)
  return(output_list)
}



## Subgrouping function for delta fractional datasets (AM)
func_DeltaFracSubgroups_AM = function(delta_dataset, grp_var_name, same_dir=FALSE, thresh =0.5){
  attach(delta_dataset)
  if (same_dir) {
    # when biomarkers go in same direction for everyone
    high_grp=delta_dataset[abs(grp_var_name-1) >= thresh,]
    low_grp=delta_dataset[abs(grp_var_name-1) < thresh,] # here the group that doesn't go as high
    num_highr=dim(high_grp)[1] 
    num_lowr=dim(low_grp)[1] 
    # Get group names from grp_var_name, to paste it to get into plot legend
    high_name <- paste("High", deparse(substitute(grp_var_name)), sep = ".") ##SJ
    low_name <- paste("Low", deparse(substitute(grp_var_name)), sep = ".") ##SJ
    # Code new group var
    df_highr_grp_ID <- tibble::rownames_to_column(high_grp, "ID")
    df_highr_grp_ID <- add_column(df_highr_grp_ID, 
                                  Group = high_name, .after = 1)
    # df_highr_grp_ID$Group <- as.character(high_name)
    df_lowr_grp_ID <- tibble::rownames_to_column(low_grp, "ID")
    df_lowr_grp_ID <- add_column(df_lowr_grp_ID, 
                                 Group = low_name, .after = 1)
    # df_lowr_grp_ID$Group <- as.character(low_name)
    df_both_grp_ID <- rbind(df_highr_grp_ID, df_lowr_grp_ID)
    df_both_grp_ID[1:2] <- lapply(df_both_grp_ID[1:2], factor)
    output_list=list("higher_grp" = high_grp, "lower_grp" = low_grp,
                     "both_grp" = df_both_grp_ID, 
                     "num_higher_grp" = num_highr,
                     "num_lower_grp" = num_lowr)
    
  } else {
    incr_grp=delta_dataset[grp_var_name >= 1,]
    decr_grp=delta_dataset[grp_var_name < 1,]
    num_incr=dim(incr_grp)[1]
    num_decr=dim(decr_grp)[1]
    # Get group names ### passing column name not working atm
    inc_name <- paste("Incr", deparse(substitute(grp_var_name)), sep = ".")
    dec_name <- paste("Decr", deparse(substitute(grp_var_name)), sep = ".")
    # Code new group var
    df_incr_grp_ID <- tibble::rownames_to_column(incr_grp, "ID")
    df_incr_grp_ID <- add_column(df_incr_grp_ID, 
                                 Group = inc_name, .after = 1)
    # df_highr_grp_ID$Group <- as.character(high_name)
    df_decr_grp_ID <- tibble::rownames_to_column(decr_grp, "ID")
    df_decr_grp_ID <- add_column(df_decr_grp_ID, 
                                 Group = dec_name, .after = 1)
    # df_lowr_grp_ID$Group <- as.character(low_name)
    df_both_grp_ID <- rbind(df_incr_grp_ID, df_decr_grp_ID)
    df_both_grp_ID[1:2] <- lapply(df_both_grp_ID[1:2], factor)
    output_list=list("increased_grp" = incr_grp, "decreased_grp" = decr_grp,
                     "both_grp" = df_both_grp_ID,
                     "num_increased_grp" = num_incr,"num_decreased_grp" = num_decr)
  }
  detach(delta_dataset)
  return(output_list)
}

## ggScatter plot of selected biomarkers
func_cor_scatter_plot= function(dataset,selected_cols=1:ncol(dataset), title=""){
  ggpairs(dataset, columns = selected_cols,
          lower = list(continuous ="smooth"), progress = F)+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
}

# Single scatter plot
# sp <- ggscatter(df, x = "wt", y = "mpg",
#    add = "reg.line",  # Add regressin line
#    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#    conf.int = TRUE # Add confidence interval
#    )
# # Add correlation coefficient
# sp + stat_cor(method = "pearson", p.accuracy=0.001, r.accuracy=0.001)

## More strict than sig.corr_matrics by only looking at significant p-values.
func_cor_pval_plot = function(dataset,method="spearman", thresh = 0.7, sig_level=0.05,title=""){
  corr_matrix = func_sig.corr_matrices(cor(dataset,method=method), threshold = thresh)
  corr_pmat = cor_pmat(dataset, method=method)
  ggcorrplot(corr_matrix, p.mat = corr_pmat, type="lower",
             insig = "blank", sig.level = sig_level, title=title) %>% ggplotly(width=850, height=750)
}


## For descriptive statistics (see below)
func_Mode = function(x) {
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

func_neat_descriptive_stats = function(datatable,title="", sig.figs=3){
  
  stats_table=data.frame(round(cbind(apply(datatable,2,mean), apply(datatable,2,sd),
                                     apply(datatable,2,Mode),
                                     t(apply(datatable,2,quantile))), sig.figs))
  stats_table=pander(stats_table, keep.line.breaks = TRUE,style="grid",caption=title,
                     col.names=c("Mean", "Std. Dev", "Mode",
                                 "Min", "1st Q", "Median\n (2nd Q)",
                                 "3rd Q","Max"),
                     row.names=c("Body weight", "LDL Cholestrol\n (mg/dL)", 
                                 "HDL Cholestrol\n (mg/dL)",
                                 "Triglyceride, Tg\n (mg/dL)",
                                 "Glucose-Fasting\n (mg/dL)", 
                                 "Insulin-Fasting\n (uIU/mL)",
                                 "C-Reactive Protein\n (ug/mL)",
                                 "TNF-alpha \n (pg/mL)",
                                 "Triiodothyronine, T3\n (ng/dL)", 
                                 "Cortisol\n (ug/dL)", "IGF-1\n (ng/mL)",
                                 "IGFBP-1\n (pg/mL)", "IGFBP-3\n (ng/mL)"),
                     split.table=Inf)
  return (stats_table)
}

func_elbow_test = function(criterion){
  attach(criterion)
  plot(1:dim(criterion[1]),AdjR2*100,pch=20,type='b',main="Elbow test", xlab="Number of explanatory variables", ylab="AdjR2")
  plot(1:dim(exhaustive_models[1]),Cp, pch=20, type="b", col="red",
       main="Elbow test", xlab="Number of explanatory variables", ylab="Mallow's Cp")
  plot(1:dim(exhaustive_models[1]),BIC,  pch=20, type="b", col="blue",
       main="Elbow test", xlab="Number of explanatory variables", ylab="BIC")
  detach(criterion)
  
}

func_resids_diagnostics = function(fit,explanatory_vars,title="", ident=FALSE,
                                   subtitle="", cooks_title=""){
  par(mfcol=c(2,2),oma = c(2,0,1,0))
  n = dim(explanatory_vars)[1]
  
  plot(c(1:n),hat(explanatory_vars),ylim=c(0,max(hat(w))), 
       main="(a) Leverage for datapoints", 
       xlab="Observation Number",
       ylab="Leverage",pch=20)
  # segments(c(1:n),0,c(1:n),hat(w))
  abline(h=(2*3)/n,col="red") #2(p+1)/n
  
  # par(pty="s")
  p2=qqnorm(residuals(fit),
            main="(b) Quantile-Quantile Plot",
            xlab="Gaussian Quantiles",
            ylab="Residuals",pch=20)
  qqline(residuals(fit), col="red")
  abline(v=0,col="Blue",lty=2)
  abline(h=0,col="Blue",lty=2)
  if (ident){
    identify(p2, plot = TRUE)
  }
  
  # par(pty="m")
  
  plot(fitted(fit),residuals(fit),
       main="(c) Residual Plot",
       xlab="Fitted Values",
       ylab="Residuals",pch=20)
  abline(h=0,col="blue",lty=2)
  abline(h=c(1.96,-1.96),lty=2,col="red")
  lines(lowess(fitted(fit),residuals(fit)), col="red")
  # identify(fitted(fit),residuals(fit))
  
  plot(fitted(fit),abs(residuals(fit)),
       main="(d) Absolute Residual Plot",
       xlab="Fitted Values",
       ylab="Absolute Residuals",pch=20)
  lines(lowess(fitted(fit),abs(residuals(fit))), col="red")
  # identify(fitted(fit),abs(residuals(fit)))
  mtext(title, col = "red",
        side=3,line=-1.5,outer=T,cex=1.2,font=2)
  mtext(subtitle,
        side=1,line=1,outer=T,cex=1.1,col='blue')
  
  par(mfrow=c(1,1), pty="m")
  plot(cooks.distance(fit),pch=20, 
       main=cooks_title, 
       col.main="red", ylab="Cook's Distance", xlab="Observation")
  abline(h=0.01,col='blue', lty=2)
  mtext(subtitle,
        side=1,line=0.5,outer=T,cex=1.1,col='blue')
  if (ident){
    identify(cooks.distance(fit), plot = TRUE, col='red')
  }
}