
### Correlation Functions ###

# Corr Matrix function (SJ)
func_SigCorrMmatrices <- function(corr_mat, threshold = 0.4){
  corr_mat[abs(corr_mat)<threshold] = 0 # or 'NA' depending on visualisation used.
  return(corr_mat)
}

# Get only lower triangle of the correlation matrix
# This removes diagonal ones and duplicates
func_GetLowerTri<-function(cormat){
  cormat[upper.tri(cormat, diag = TRUE)] <- NA
  return(cormat)
}

# Function to output formatted table to html via kable
func_Kable <- function(df) {
  knitr::kable(
    df, booktabs = TRUE, digits = 3, 
    align = "c", caption = table_title) %>% 
    kable_styling(bootstrap_options =c(
      "striped", "scale_down","hover", "condensed", "responsive"))
}

##### 
# SJ edits

# Corr function
func_Cor_SJ <- function(df, method = "pearson") {
  M <- Hmisc::rcorr(as.matrix(df), type = method)
  Md <- map(M, ~func_GetLowerTri(.x)) #rm diagonal and duplicates
  Mdf <- map(Md, ~data.frame(.x))
  return(Mdf)
}

# Cor output to df
func_FormattedCor_df_SJ <- function(df, p_val = 0.05, corr_method = "pearson"){
  func_Cor_SJ(df, method = corr_method) %>%
    map(~rownames_to_column(.x, var="Marker1")) %>%
    map(~pivot_longer(.x, -Marker1, "Marker2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    rename(N = n) %>%
    rename(p = P) %>%
    mutate(sig_p = ifelse(p < p_val, T, F), # set significance
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA),
           Method = ifelse(corr_method=="pearson", "Pearson's", "Spearman's"))
}
# Cor final table output
# Usage: func_FormattedCor_SJ(df_BL_centr, p_val, r_threshold, methods)
func_FormattedCor_SJ <- function(df, p_val=0.05, r_threshold, 
                                 methods = c("pearson", "spearman"), timepts.option=FALSE,
                                 timepts.exclude=c("_blef","_bsl")){
  df_cor_table = NULL
  for (i in 1:length(methods)){
    ## Baseline only cor
    df_corr <- func_FormattedCor_df_SJ(df, p_val, corr_method = methods[i])
    # Drop rows containing NA's (self-cor, non-unique, non-significant p-values)
    df_corr <- df_corr %>% drop_na()
    if (timepts.option==T){
      # Only BLEF (default) in Marker1 column
      df_corr <- df_corr[!grepl(timepts.exclude[1],df_corr$Marker1),]
      # Only BSL (default) in Marker2 column
      df_corr <- df_corr[!grepl(timepts.exclude[2],df_corr$Marker2),]
    }
    # Sort by descending absolute R value
    df_corr <- df_corr %>% arrange(desc(abs(r)))
    # Remove rows with R lower than R_threshold
    df_corr <- df_corr %>% subset(abs(r) >= r_threshold)
    # Remove extra columns
    df_corr <- df_corr %>% select(-c(6:8))
    # Join dfs
    df_cor_table <- base::rbind(df_cor_table, df_corr)
  }
  #  and remove duplicated pairs
  df_cor_table <- distinct(df_cor_table, Marker1, Marker2, .keep_all = TRUE)
  # Switch Marker1 Marker2 names
  df_cor_table <- df_cor_table %>% rename(Marker1 = Marker2, Marker2 = Marker1)
  # Add row numbers, re-order
  df_cor_table <- df_cor_table %>%
    mutate(Nr = row_number()) %>% 
    select(Nr, Marker1, Marker2, r, p, Method, N) # reorder
  return(df_cor_table)
}
#####