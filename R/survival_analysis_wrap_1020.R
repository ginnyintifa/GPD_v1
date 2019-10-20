

###  Association analysis between somatic units and OS


#' This function performs survival analysis for mutation counts mapped to PIU, LU and NCU.
#' 
#' @param piu_filename Filename for the PIU mapping results.
#' @param lu_filename Filename for the LU mapping results. 
#' @param ncu_filename Filename for the NCU mapping results. 
#' @param clinical_df Clinical information data frame. 
#' @param gender_as_covariate  Boolean variable indicating gender should be taken as a covariable. 
#' @param race_group_min The minimum number of patients with the same race enough for a race group.
#' @param min_surv_days The minimum number of survival time (in days) a patient has to be included in analysis.
#' @param min_surv_people The minimum number of patients survived (or censored) in a cohort that a survival analysis should be conducted.
#' @param patient_sum_min The minimum number of patients having mutation mapping to the unit for the unit to be analysed.
#' @param mutation_type A tag indicating the type of mutation; somatic or germline.
#' @param output_dir The directory you would like to have your output files in.
#' @import dplyr magrittr data.table qvalue survival
#' @export
#' @details 
#' @examples 
#' univariate_cox_model(piu_filename = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/example/piu_mapping_count.tsv",
#'                      lu_filename = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/example/lu_summarising_count.tsv",
#'                      ncu_filename = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/example/ncu_summarising_count.tsv",
#'                      barcode_stage_filename = "a file with stage information"
#'                      clinical_df = sel_example_cdr,
#'                      race_group_min = 6,
#'                      min_surv_days = 90,
#'                      min_surv_people = 5,
#'                      patient_sum_min = 3,
#'                      mutation_type = "somatic",
#'                      output_dir = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/example/")
#'                       


univariate_cox_model = function(piu_filename,
                                lu_filename,
                                ncu_filename,
                                barcode_stage_filename,
                                clinical_df,
                             #   gender_as_covariate = T,
                                race_group_min = 6,
                                min_surv_days = 90,
                                min_surv_people = 5,
                                patient_sum_min = 3,
                                mutation_type = "somatic",
                                output_dir)
  
{
  
  get_barcode_stage = fread(barcode_stage_filename, stringsAsFactors = F)
  
  #### piu
  
  if(file.exists(piu_filename))
  {
    piu_unite = piu_counts_cdr_clinical_unite(
      piu_count_filename = piu_filename,
      cdr_clinical = cdr_clinical,
      patient_sum_min = patient_sum_min,
      output_dir = output_dir,
      output_name = paste0(mutation_type,"_piu_cdr_clinical_unite.tsv"))
    
    if(file.exists(paste0(output_dir,paste0(mutation_type,"_piu_cdr_clinical_unite.tsv"))))
    {
      surv_piu_info = cdr_tidy_up_for_model(
        interest_variable_info = piu_unite[[3]],
        unite_data = piu_unite[[1]],
        race_group_min = race_group_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_piu_survival_info.tsv"))
      
    }
  }
    
    ####bpiu
  
    if(file.exists(bpiu_filename))
    {
      
      ### an additional step for selecting valid genes 
      
      
      bpiu_unite = gene_counts_cdr_clinical_unite(
        gene_count_filename = bpiu_filename,
        cdr_clinical = cdr_clinical,
        patient_sum_min = patient_sum_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_bpiu_cdr_clinical_unite.tsv"))
      
      
      
      if(file.exists(paste0(output_dir,paste0(mutation_type,"_bpiu_cdr_clinical_unite.tsv"))))
      {
        surv_bpiu_info = cdr_tidy_up_for_model(
          interest_variable_info = bpiu_unite[[2]],
          unite_data = bpiu_unite[[1]],
          race_group_min = race_group_min,
          output_dir = output_dir,
          output_name = paste0(mutation_type,"_bpiu_survival_info.tsv"))
        
      }
    }
  
  ####npc 
  
  
  
  if(file.exists(npc_filename))
  {
    

    
    npc_unite = gene_counts_cdr_clinical_unite(
      gene_count_filename = npc_filename,
      cdr_clinical = cdr_clinical,
      patient_sum_min = patient_sum_min,
      output_dir = output_dir,
      output_name = paste0(mutation_type,"_npc_cdr_clinical_unite.tsv"))
    
    
    
    if(file.exists(paste0(output_dir,paste0(mutation_type,"_npc_cdr_clinical_unite.tsv"))))
    {
      surv_npc_info = cdr_tidy_up_for_model(
        interest_variable_info = npc_unite[[2]],
        unite_data = npc_unite[[1]],
        race_group_min = race_group_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_npc_survival_info.tsv"))
      
      
    }
  }
  
  
  
  surv_names = grep("ENSG", colnames(surv_piu_info), value = T, invert = T)
  
  piu_names = grep("ENSG", colnames(surv_piu_info), value = T)
  
  bpiu_names = grep("ENSG", colnames(surv_bpiu_info), value = T)
  
  npc_names = grep("ENSG", colnames(surv_npc_info), value = T)
  
  if(length(piu_names) >0)
  {
    colnames(surv_piu_info) = c(surv_names,paste0("_piu_", piu_names))
    
  }
  
  
  if(length(bpiu_names)>0)
  {
    colnames(surv_bpiu_info) = c(surv_names,paste0("_bpiu_", bpiu_names))
    
  }
  
  if(length(npc_names)>0)
  {
    
    colnames(surv_npc_info) = c(surv_names,paste0("_npc_", npc_names))
    
  }
  
  
  
  surv_info = surv_piu_info%>%
    dplyr::left_join(surv_bpiu_info, by = surv_names)%>%
    dplyr::left_join(surv_npc_info, by = surv_names)%>%
    dplyr::filter(barcode%in%cdr_no_hyper$bcr_patient_barcode)
  
  
  
  filter_surv_os_data = filter_survival_data(surv_info_data = surv_info,
                                             surv_status = quo(OS),
                                             surv_time = quo(OS.time),
                                             surv_race = quo(os_race),
                                             min_surv_days = 90)
  
  
  cancer_mc3 = fread(paste0("/data/ginny/tcga_pancan/TCGA_all/",
                            this_cancer,"_somatic/",
                            this_cancer, "_somatic_mc3.tsv"),
                     stringsAsFactors = F)
  
  
  ns_mc3 = cancer_mc3%>%
    dplyr::filter(! Variant_Classification == "Silent")
  
  
  all_barcodes = unlist(lapply(1:nrow(ns_mc3), function(x) {
    
    
    this_row = ns_mc3$agg_sample_id[x]
    
    sep_id = unlist(strsplit(this_row, split = "_", fixed = T))
    
    
    
  }))
  
  
  
  barcode_df = as.data.frame(table(all_barcodes))%>%
    dplyr::arrange(desc(Freq))
  
  barcode_totalMut = data.frame(barcode = as.character(barcode_df$all_barcodes), 
                                total_mutation = barcode_df$Freq,
                                stringsAsFactors = F)
  
  
  
  surv_with_totalMut = barcode_totalMut %>%
    dplyr::left_join(filter_surv_os_data, by = "barcode")%>%
    dplyr::left_join(get_barcode_stage, by = "barcode")%>%
    na.omit()
  
  test_units = univariate_survival_significance_adjust_totalMut_stage(filter_surv_data = surv_with_totalMut,
                                                                output_dir = output_dir,
                                                                surv_name = "os_significance",
                                                                data_fold = "whole")
  
  
  
  write.table(test_units, paste0(output_dir,"all_units_results.tsv"),
              quote = F, row.names = F, sep = "\t")
  
  
  
  
}
  
  


  