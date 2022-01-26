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


# To select pathway analysis required (options: Pathway.Enr, BP, CC, MF)
func_path_enr <- function(protein_Accs_list, type = "Path.Enr"){
  if (type == "Path.Enr"){
    res <- Pathway.Enr(protein_Accs_list)
  } else if (type == "BP"){
    res<-Enrichment.BP(protein_Accs_list)
  } else if (type == "CC"){
    res<-Enrichment.CC(protein_Accs_list)
  } else if (type == "MF"){
    res<-Enrichment.MF(protein_Accs_list)
  }
  return (res)
}

# Returns pvalue, pathway, common proteins belonging to pathway for Enrichment.obj
func_extract_relevant_paths<- function(enr.obj){
  data_enr.obj<-enr.obj$data
  res<-data.frame(p_val = data_enr.obj$p_value,
                  pathway = data_enr.obj$term_name,
                  common_proteins = data_enr.obj$intersection)
  return (res)
}

# To get count of common elements (ordered highest to lowest)
# Although works well for pathways, can be used for other cases.
func_count_common <- function(pathways){
  count_table<-as.data.frame(table(pathways))
  ordered_table<-count_table[order(count_table$Freq, decreasing=T),]
  return (ordered_table)
}



# Like func_count_common but specific to finding number of paths that each protein is involved in
func_count_proteins_paths<- function(pathways_obj, uniprotDict_obj,uniprotTargetName = T){
  comb_accs_list<-paste(pathways_obj$common_proteins, collapse=",") # combine all accs_nos across all pathways into one long character separated by ",".
  vect_accs_list<-unlist(strsplit(comb_accs_list,",")) # separate to 1 AccsNo per vector entry.
  freq_table <- func_count_common(vect_accs_list) # Count how many times a single protein appears
  freq_table[,1]<-func_AccsToProtNames(freq_table$pathways,uniprotDict_obj, uniprot_target_name = uniprotTargetName) # to Actual Protein names (in Uniprot)
  return (freq_table)
}


# To get into format ready for network plot: pathfindR::term_gene_graph()
func_TableGraphNetwork_Omics <- function(named_cor_table_obj,pathways_obj, uniprotDict_obj, uniprotTargetName = T){
  
  
  # Uncomment either the following chunk or the one below
  # cor_positive_prots<-NULL
  # cor_negative_prots<-NULL
  # rownames(named_cor_table_obj)<-func_NamesDotsToTargetNames(rownames(named_cor_table_obj),uniprotDict_obj)
  # for (i in 1:dim(named_cor_table_obj)[1]){
  #   if (named_cor_table_obj$cor_pearsons[i] >=0){
  #     cor_positive_prots<-cbind(cor_positive_prots,
  #                               rownames(named_cor_table_obj)[i])
  #     # element goes in "positive" list
  #   } else {
  #     cor_negative_prots<-cbind(cor_negative_prots,
  #                               rownames(named_cor_table_obj)[i])
  #     # element goes in "negative" list
  #   }
  # }
  
  # Uncomment either the following chunk or the one above
  
  # Separating strongly correlated proteins into 2 lists - positive and negative
  cor_positive_prots<-NULL
  cor_negative_prots<-NULL
  # Transforming names to Actual Protein names.
  rownames(named_cor_table_obj)<-func_NamesDotsToTargetNames(rownames(named_cor_table_obj),
                                                             uniprotDict_obj)
  # Segregating into positive and negative lists
  for (element in rownames(named_cor_table_obj)){
    index<-which(rownames(named_cor_table_obj) == element) # identify index corresponding to protein
    if (named_cor_table_obj[index,] >= 0){
      cor_positive_prots<-cbind(cor_positive_prots,element)
      # element goes in "positive" list
    } else {
      cor_negative_prots<-cbind(cor_negative_prots,element)}
    # element goes in "negative" list
  }
  
  # Now we do not need the pearsons correlations anymore...
  
  # Create shallow copy which will be output
  df_network_plot<-with(pathways_obj,
                        data.frame(lowest_p = p_val,
                                   Term_Description = pathway,
                                   Up_regulated = character(length(p_val)), 
                                   Down_regulated = character(length(p_val))))
  
  # Each row is a string of proteins belonging to a specific pathway. (hence no duplicates possible)
  accs_pathways<-pathways_obj$common_proteins
  for (i in 1:length(accs_pathways)){
    expanded_accs_path_list<-unlist(strsplit(accs_pathways[i], ","))
    prots_pathway<-func_AccsToProtNames(expanded_accs_path_list,
                                        uniprotDict_obj,
                                        uniprot_target_name = uniprotTargetName)
    
    # Separating positive proteins from negative ones for every pathway
    pos_prots_pathway<-intersect(prots_pathway, cor_positive_prots)
    df_network_plot$Up_regulated[i]<-paste(pos_prots_pathway, collapse = ", ")
    
    neg_prots_pathway<-intersect(prots_pathway, cor_negative_prots)
    df_network_plot$Down_regulated[i]<-paste(neg_prots_pathway, collapse = ", ")
    
    if ((length(neg_prots_pathway) + length(pos_prots_pathway))!=length(prots_pathway)){
      stop("Some proteins in pathway not accounted for.")
    }
  }
  return (df_network_plot)
}

# Gene_Graph_Plot for Protein Network and Identifying Hubs
# Filename: .jpg
func_NetworkPlot_GeneGraph_Omics <- function(TableGraphNetwork_obj, filename,
                                             Res = 200, Width = 4000, Height = 2000,
                                             numTerms = nrow(TableGraphNetwork_obj)){
  jpeg(filename,width=Width,height=Height,res=Res)
  pathfindR::term_gene_graph(TableGraphNetwork_obj,
                             num_terms = numTerms,
                             layout = "stress",
                             use_description = TRUE,
                             node_size = "num_genes")
}