library(EnhancedVolcano)
library(dplyr)
library(forecast)
library(UniprotR)
# library(pathfindR)

BoxCoxTransf = function(x){
  lambda = BoxCox.lambda(x)
  if(lambda == 1){y = log(x)} else {
    y = x^(lambda-1)/lambda}
  if(cor(x,y, use = "pairwise.complete.obs") >= 0){return(y)} else {return(-y)}}

MWAS = function(x,y=group,set1,set2){
  x_new = BoxCoxTransf(x)
  x_new_set1 = x_new[which(y == set1)]
  x_new_set2 = x_new[which(y == set2)]
  p.ttest = t.test(x_new_set1,x_new_set2)$p.val
  fc = mean(x[which(y == set1)]) / mean(x[which(y == set2)])
  # fc = mean(x[which(y == set1)]/x[which(y == set2)])
  log2fc = log(fc,2)
  output = c(log2fc,p.ttest)
  return(output)}

create_results_table = function(df,group1,group2){
  res = t(apply(df,1,MWAS,set1=group1,set2=group2))
  res = data.frame(res)
  colnames(res) = c("log2FC","pval")
  res$FDR = p.adjust(res$pval, method = "fdr") #Default is Holm's
  # Eg. (a<-c(rep(3,4),rep(5,5), rep(1, 11))/100); p.adjust(a); p.adjust(a, "fdr")
  res = res[order(res$pval),]
  return(res)}



func_merge_metabolomics_biomarker <-function(metabolomics_data, biomarker_data, biomarker){
  # Transpose to match on rownames (i.e. Timepoint + Patient ID)
  df1<-data.frame(t(metabolomics_data)); df2 <- data.frame(t(biomarker_data))
  attach(df2)
  named_biomarker=data.frame(biomarker, row.names = row.names(df2))
  detach(df2)
  
  # Clean merged output
  df_merged <- merge(named_biomarker, df1, by=0) #biomarker is first column
  row.names(df_merged) <- df_merged$Row.names
  df_merged <- df_merged[,-1]
  
  df_merged <- data.frame(t(df_merged))
  return (df_merged)
}

# To get pearsons correlations and p-val for biomarker with metabolomics
# NOTE: Default: 1st row is biomarker, followed by other metabolites. Change row_num argument if it's not the case.
func_corr_results_table <- function(merged_dataset,biomarker_row_num = 1){
  # Empty space to fill in
  results_table <- NULL
  
  # Segregating biomarker and metabolites
  biomarker <- as.numeric(merged_dataset[biomarker_row_num,])
  metabolites <- merged_dataset[-biomarker_row_num,]
  
  # Getting correlations and p-values for each biomarker-metabolite pair
  for (i in 1:dim(metabolites)[1]){
    test <- cor.test(biomarker,as.numeric(metabolites[i,]), method="pearson")
    cor_pearsons <- test$estimate
    pval <- test$p.value
    
    # Updating table to include all pairs
    results_table = rbind(results_table, 
                          cbind(row.names(metabolites)[i],cor_pearsons,pval))
  }
  # Cleaning results table
  rownames(results_table) <- results_table[,1]
  results_table <- data.frame(results_table[,-1])
  results_table %>% mutate_if(is.character,as.numeric) -> results_table
  results_table$FDR = p.adjust(results_table$pval, method = "fdr") #Default is Holm's
  
  return (results_table[order(results_table$FDR),])
}

func_volcano_plot <- function(correlation_table, biomarker, main_title = "Biomarker Volcano Plot"){
  EnhancedVolcano(correlation_table, pCutoff = 0.01, FCcutoff = 0.5,
                  lab = rownames(correlation_table), labSize = 2, 
                  xlab = bquote(~ rho["Pearsons"] ~ .(biomarker)), xlim=c(-1,1),
                  x = 'cor_pearsons', y = 'FDR', title=main_title,
                  subtitle="Using FDR (adjusted) p-values",
                  selectLab = rownames(correlation_table)[which(abs(correlation_table$cor_pearsons)>=0.5 & correlation_table$pval<=0.05)],
                  legendPosition = 'right',legendLabSize = 9,legendIconSize = 3,
                  boxedLabels = TRUE, labFace='bold',
                  drawConnectors = TRUE,
                  widthConnectors = 0.75)
}



# To get accession numbers from UniProt website for pathway enrichment analysis
# Inputs protein short forms (with or w/o dots) to return Accs_no.
# If uniprot_target_name = T --> using Target column
# Defaults to Target_dots column
# eg. func_Accs_no_proteins("IL.10.Ra", UniprotDict (see under section Datasets)) --> "Q13651"
# eg. func_Accs_no_proteins("IL-10 Ra", UniprotDict (see under section Datasets), uniprot_target_name = T) --> "Q13651"
func_Accs_no_proteins <- function(protein_list, df_UniProt_IDs, uniprot_target_name=F){
  prot_list<-func_prots_X_numeric(protein_list)
  if (uniprot_target_name){
    accs=na.omit(df_UniProt_IDs[match(prot_list,df_UniProt_IDs$Target),]$UniProt)
  } else {
    accs<-na.omit(df_UniProt_IDs[match(prot_list,df_UniProt_IDs$Target_dots),]$UniProt)
  }
  return(accs)
}

# Protein Accs No (UniprotID) to Protein Name (or DottedName)
# Defaults to Actual Protein Name
# eg. func_AccsToProtNames(Q13651, Uniprot Dict) --> "IL.10.Ra"
# eg. func_AccsToProtNames(Q13651, Uniprot Dict, uniprot_target_name = T) --> "IL-10 Ra"
func_AccsToProtNames <- function(accs_list, df_UniProt_IDs, uniprot_target_name=F){
  if (uniprot_target_name){
    prots=na.omit(df_UniProt_IDs[match(accs_list,df_UniProt_IDs$UniProt),]$Target)
  } else {
    prots<-na.omit(df_UniProt_IDs[match(prot_list,df_UniProt_IDs$UniProt),]$Target_dots)
  }
  return(prots)
}

# Prot names (with special characters subbed in by ".") to Actual Protein Name 
# (eg. LEAP.1 to LEAP-1, IL.10.Ra to IL-10 Ra)
# Inputs: Vector of dottedNamed Proteins; Uniprot Dict (see below under section Datasets)
func_NamesDotsToTargetNames <- function(protDots_list, df_UniProt_IDs){
  prot_list<-func_prots_X_numeric(protDots_list)
  prots=na.omit(df_UniProt_IDs[match(prot_list,df_UniProt_IDs$Target_dots),]$Target)
  return(prots)
}


# Removes X as prefix if 2nd character is numeric.
func_prots_X_numeric <- function(protein_list){
  for (i in 1:length(protein_list)){
    list_chars<-unlist(strsplit(protein_list[i],""))
    # Take proteins that start with X and followed by a number, removes X in front.
    if (list_chars[1]=="X" & is.numeric(as.numeric(list_chars[2])) & 
        !is.na(as.numeric(list_chars[2]))){
      protein_list[i] <- paste(list_chars[-1],collapse = "")
    }
  }
  return (protein_list)
}
