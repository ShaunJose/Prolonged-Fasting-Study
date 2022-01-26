
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

# Pearsons corr function
func_CorPears <- function(df) {
  M <- Hmisc::rcorr(as.matrix(df))
  Md <- map(M, ~func_GetLowerTri(.x)) #rm diagonal and duplicates
  Mdf <- map(Md, ~data.frame(.x))
  return(Mdf)
}

# Spearmans corr function
func_CorSpearm <- function(df) {
  M <- Hmisc::rcorr(as.matrix(df), type = "spearman")
  Md <- map(M, ~func_GetLowerTri(.x))
  Mdf <- map(Md, ~data.frame(.x))
  return(Mdf)
}

# Cor output to df - Pearsons
func_FormattedCorPears <- function(df, p_val){
  func_CorPears(df) %>%
    map(~rownames_to_column(.x, var="Marker1")) %>%
    map(~pivot_longer(.x, -Marker1, "Marker2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    rename(N = n) %>%
    rename(p = P) %>%
    mutate(sig_p = ifelse(p < p_val, T, F), # set significance
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA),
           pearsons = ifelse(p < 1, T, F))
}

# Cor output to df - Spearmans
func_FormattedCorSpearm <- function(df, p_val){
  func_CorSpearm(df) %>%
    map(~rownames_to_column(.x, var="Marker1")) %>%
    map(~pivot_longer(.x, -Marker1, "Marker2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    rename(N = n) %>%
    rename(p = P) %>%
    mutate(sig_p = ifelse(p < p_val, T, F), # set significance
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA),
           spearmans = ifelse(p < 1, T, F)) 
}

# Cor final table output
# Usage: func_FormattedCor(df_BL_centr, p_val, r_threshold)
func_FormattedCor <- function(df, p_val, r_threshold){
  ## Baseline only cor
  df_pears <- func_FormattedCorPears(df, p_val)
  df_spearm <- func_FormattedCorSpearm(df, p_val)
  # Drop rows containing NA's (self-cor, non-unique, non-significant p-values)
  df_pears <- df_pears %>% drop_na()    
  df_spearm <- df_spearm %>% drop_na() 
  # Sort by descending absolute R value
  df_pears <- df_pears %>% arrange(desc(abs(r)))
  df_spearm <- df_spearm %>% arrange(desc(abs(r)))
  # Remove rows with R lower than R_threshold
  df_pears <- df_pears %>% subset(abs(r) >= r_threshold)
  df_spearm <- df_spearm %>% subset(abs(r) >= r_threshold)
  # Add method label, remove extra columns
  df_pears <- df_pears %>% mutate(
    Method =  case_when(pearsons = TRUE ~ "Pearson's")) %>% 
    select(-c(6:9))
  df_spearm <- df_spearm %>% mutate(
    Method =  case_when(spearmans = TRUE ~ "Spearman's")) %>% 
    select(-c(6:9))
  # Join dfs and remove non-unique spearman
  df_cor_table <- base::rbind(df_pears, df_spearm)
  df_cor_table <- distinct(
    df_cor_table, Marker1, Marker2, .keep_all = TRUE)
  df_cor_table <- df_cor_table %>%
    mutate(Nr = row_number()) %>% 
    select(Nr, Marker1, Marker2, r, p, Method, N) # reorder
  return(df_cor_table)
}

# Cor final table output 2 - for BL vs BLEF fold change
# Usage: func_FormattedCor(df_BL_centr, p_val, r_threshold)
func_FormattedCor2 <- function(df, p_val, r_threshold){
  ## Baseline only cor
  df_pears <- func_FormattedCorPears(df, p_val)
  df_spearm <- func_FormattedCorSpearm(df, p_val)
  # Drop rows containing NA's (self-cor, non-unique, non-significant p-values)
  df_pears <- df_pears %>% drop_na()    
  df_spearm <- df_spearm %>% drop_na()
  # Drop Baseline from Marker1 column
  df_pears <- 
    df_pears[!grepl("_bsl",df_pears$Marker1),]
  df_spearm <- 
    df_spearm[!grepl("_bsl",df_spearm$Marker1),]
  # Drop BLEF from Marker2 column
  df_pears <- 
    df_pears[!grepl("_blef_fold",df_pears$Marker2),]
  df_spearm <- 
    df_spearm[!grepl("_blef_fold",df_spearm$Marker2),]
  # Sort by descending absolute R value
  df_pears <- df_pears %>% arrange(desc(abs(r)))
  df_spearm <- df_spearm %>% arrange(desc(abs(r)))
  # Remove rows with R lower than R_threshold
  df_pears <- df_pears %>% subset(abs(r) >= r_threshold)
  df_spearm <- df_spearm %>% subset(abs(r) >= r_threshold)
  # Add method label, remove extra columns
  df_pears <- df_pears %>% mutate(
    Method =  case_when(pearsons = TRUE ~ "Pearson's")) %>% 
    select(-c(6:9))
  df_spearm <- df_spearm %>% mutate(
    Method =  case_when(spearmans = TRUE ~ "Spearman's")) %>% 
    select(-c(6:9))
  # Join dfs and remove non-unique spearman
  df_cor_table <- base::rbind(df_pears, df_spearm)
  df_cor_table <- distinct(
    df_cor_table, Marker1, Marker2, .keep_all = TRUE)
  # Switch Marker1 Marker2 names
  df_cor_table <- df_cor_table %>% rename(Marker1 = Marker2, Marker2 = Marker1 )
  # Add row numbers, re-order
  df_cor_table <- df_cor_table %>%
    mutate(Nr = row_number()) %>% 
    select(Nr, Marker1, Marker2, r, p, Method, N) # reorder
  return(df_cor_table)
}

# Cor final table output 3 - for BL vs BLEF absolute change
# Usage: func_FormattedCor(df_BL_centr, p_val, r_threshold)
func_FormattedCor3 <- function(df, p_val, r_threshold){
  ## Baseline only cor
  df_pears <- func_FormattedCorPears(df, p_val)
  df_spearm <- func_FormattedCorSpearm(df, p_val)
  # Drop rows containing NA's (self-cor, non-unique, non-significant p-values)
  df_pears <- df_pears %>% drop_na()    
  df_spearm <- df_spearm %>% drop_na()
  # Drop Baseline from Marker1 column
  df_pears <- 
    df_pears[!grepl("_bsl",df_pears$Marker1),]
  df_spearm <- 
    df_spearm[!grepl("_bsl",df_spearm$Marker1),]
  # Drop BLEF from Marker2 column
  df_pears <- 
    df_pears[!grepl("_blef_abs",df_pears$Marker2),]
  df_spearm <- 
    df_spearm[!grepl("_blef_abs",df_spearm$Marker2),]
  # Sort by descending absolute R value
  df_pears <- df_pears %>% arrange(desc(abs(r)))
  df_spearm <- df_spearm %>% arrange(desc(abs(r)))
  # Remove rows with R lower than R_threshold
  df_pears <- df_pears %>% subset(abs(r) >= r_threshold)
  df_spearm <- df_spearm %>% subset(abs(r) >= r_threshold)
  # Add method label, remove extra columns
  df_pears <- df_pears %>% mutate(
    Method =  case_when(pearsons = TRUE ~ "Pearson's")) %>% 
    select(-c(6:9))
  df_spearm <- df_spearm %>% mutate(
    Method =  case_when(spearmans = TRUE ~ "Spearman's")) %>% 
    select(-c(6:9))
  # Join dfs and remove non-unique spearman
  df_cor_table <- base::rbind(df_pears, df_spearm)
  df_cor_table <- distinct(
    df_cor_table, Marker1, Marker2, .keep_all = TRUE)
  # Switch Marker1 Marker2 names
  df_cor_table <- df_cor_table %>% rename(Marker1 = Marker2, Marker2 = Marker1 )
  # Add row numbers, re-order
  df_cor_table <- df_cor_table %>%
    mutate(Nr = row_number()) %>% 
    select(Nr, Marker1, Marker2, r, p, Method, N) # reorder
  return(df_cor_table)
}

# Function to output formatted table to html via kable
func_Kable <- function(df) {
  knitr::kable(
    df, booktabs = TRUE, digits = 3, 
    align = "c", caption = table_title) %>% 
    kable_styling(bootstrap_options =c(
      "striped", "scale_down","hover", "condensed", "responsive"))
}