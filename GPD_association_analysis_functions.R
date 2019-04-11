### functions needed for association implementation 



univariate_cox_model_for_piu = function(piu_filename,
                                                cdr_clinical,
                                                gender_as_covariate = T,
                                                race_group_min = 6,
                                                min_surv_days = 90,
                                                min_surv_people = 5,
                                                patient_sum_min = 3,
                                                mutation_type = "somatic",
                                                output_dir)
{
   
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
      piu_info = cdr_tidy_up_for_model(
        interest_variable_info = piu_unite[[3]],
        unite_data = piu_unite[[1]],
        race_group_min = race_group_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_piu_survival_info.tsv"))
      
      if(gender_as_covariate == T)
      {
        fit_survival_model_v2_test(
          surv_info_data = piu_info,
          interest_variable_info = piu_unite[[3]],
          min_surv_time = min_surv_days,
          min_surv_people = min_surv_people,
          output_dir = output_dir,
          output_name = paste0(mutation_type,"_piu_cdr_univariate_test.tsv"))
        
      }else{
        fit_survival_model_no_gender_v2_test(
          surv_info_data = piu_info,
          interest_variable_info = piu_unite[[3]],
          min_surv_time = min_surv_days,
          min_surv_people = min_surv_people,
          output_dir = output_dir,
          output_name = paste0(mutation_type,"_piu_cdr_univariate_test.tsv"))
      }
    }else{
      cat("...this type of PIU level data not available.", "\n")
      
    }
    
    
  }
  
}



univariate_cox_model_for_bpiu = function(bpiu_filename,
                                                 cdr_clinical,
                                                 gender_as_covariate = T,
                                                 race_group_min = 6,
                                                 min_surv_days = 90,
                                                 min_surv_people = 5,
                                                 patient_sum_min = 3,
                                                 mutation_type = "somatic",
                                                 output_dir)
{
  
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
      bpiu_info = cdr_tidy_up_for_model(
        interest_variable_info = bpiu_unite[[2]],
        unite_data = bpiu_unite[[1]],
        race_group_min = race_group_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_bpiu_survival_info.tsv"))
      
      cat("cancer patient gene dimension: ", dim(bpiu_info), "\n")  
      
      
      if(gender_as_covariate == T)
      {
        fit_survival_model_v2_test(
          surv_info_data = bpiu_info,
          interest_variable_info = bpiu_unite[[2]],
          min_surv_time = min_surv_days,
          min_surv_people = min_surv_people,
          output_dir = output_dir,
          output_name = paste0(mutation_type,"_bpiu_cdr_univariate_test.tsv"))
        
      }else{
        fit_survival_model_no_gender_v2_test(
          surv_info_data = bpiu_info,
          interest_variable_info = bpiu_unite[[2]],
          min_surv_time = min_surv_days,
          min_surv_people = min_surv_people,
          output_dir = output_dir,
          output_name = paste0(mutation_type,"_bpiu_cdr_univariate_test.tsv"))
      }
    }else{
      cat("...this bPIU data not available.", "\n")
      
    }
    
    #################
    
    cat("1/1...univariate Cox regression model on bPIU counts fitted!", "\n")
    
  }else{
    cat("1/1...bPIU data not available.", "\n")
    
  }
  
  cat("Univariate Cox regression models fitted!", "\n")
  
}





univariate_cox_model_for_npc = function(npc_filename,
                                                cdr_clinical,
                                                gender_as_covariate = T,
                                                race_group_min = 6,
                                                min_surv_days = 90,
                                                min_surv_people = 5,
                                                patient_sum_min = 3,
                                                mutation_type = "somatic",
                                                output_dir)
{
  
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
      npc_info = cdr_tidy_up_for_model(
        interest_variable_info = npc_unite[[2]],
        unite_data = npc_unite[[1]],
        race_group_min = race_group_min,
        output_dir = output_dir,
        output_name = paste0(mutation_type,"_npc_survival_info.tsv"))
      
      cat("cancer patient gene dimension: ", dim(npc_info), "\n")  
      
      if(length(npc_unite[[2]])>0)
      {
        
        if(gender_as_covariate == T)
        {
          fit_survival_model_v2_test(
            surv_info_data = npc_info,
            interest_variable_info = npc_unite[[2]],
            min_surv_time = min_surv_days,
            min_surv_people = min_surv_people,
            output_dir = output_dir,
            output_name = paste0(mutation_type,"_npc_cdr_univariate_test.tsv"))
          
        }else{
          fit_survival_model_no_gender_v2_test(
            surv_info_data = npc_info,
            interest_variable_info = npc_unite[[2]],
            min_surv_time = min_surv_days,
            min_surv_people = min_surv_people,
            output_dir = output_dir,
            output_name = paste0(mutation_type,"_npc_cdr_univariate_test.tsv"))
        }
        
      }
      
      
    }else{
      cat("...this nPC data not available.", "\n")
      
    }
    

    cat("1/1...univariate Cox regression model on nPC counts fitted!", "\n")
    
  }else{
    cat("1/1...nPC data not available.", "\n")
    
  }
  
  cat("Univariate Cox regression models fitted!", "\n")
  
}


####


nonzero = function(x) sum(x != 0)

cdr_tidy_up_for_model = function(interest_variable_info, unite_data, race_group_min, 
                                 output_dir,
                                 output_name)
{
  
 
    patient_count = unite_data %>%
    dplyr::select(barcode, age_at_initial_pathologic_diagnosis,gender,one_of(interest_variable_info)) %>%
    dplyr::rename(age = age_at_initial_pathologic_diagnosis)
  
  unite_data$race[grep("\\[", unite_data$race)] = "OTHER"
  
  os_data = unite_data %>%
    dplyr::select(barcode, race, OS, OS.time) %>%
    na.omit()
  os_race_freq = as.data.frame(table(os_data$race))
  os_race_minor = os_race_freq %>%
    dplyr::filter(Freq<race_group_min)
  os_data$race[which(os_data$race %in% os_race_minor$Var1)] = "OTHER"
  colnames(os_data)[which(colnames(os_data)=="race")] = "os_race"
  
  
  dss_data = unite_data %>%
    dplyr::select(barcode, race, DSS, DSS.time) %>%
    na.omit()
  dss_race_freq = as.data.frame(table(dss_data$race))
  dss_race_minor = dss_race_freq %>%
    dplyr::filter(Freq<race_group_min)
  dss_data$race[which(dss_data$race %in% dss_race_minor$Var1)] = "OTHER"
  colnames(dss_data)[which(colnames(dss_data)=="race")] = "dss_race"
  
  
  dfi_data = unite_data %>%
    dplyr::select(barcode, race, DFI, DFI.time) %>%
    na.omit()
  
  dfi_race_freq = as.data.frame(table(dfi_data$race))
  dfi_race_minor = dfi_race_freq %>%
    dplyr::filter(Freq<race_group_min)
  dfi_data$race[which(dfi_data$race %in% dfi_race_minor$Var1)] = "OTHER"
  colnames(dfi_data)[which(colnames(dfi_data)=="race")] = "dfi_race"
  
  
  pfi_data = unite_data %>%
    dplyr::select(barcode, race, PFI, PFI.time) %>%
    na.omit()
  pfi_race_freq = as.data.frame(table(pfi_data$race))
  pfi_race_minor = pfi_race_freq %>%
    dplyr::filter(Freq<race_group_min)
  pfi_data$race[which(pfi_data$race %in% pfi_race_minor$Var1)] = "OTHER"
  colnames(pfi_data)[which(colnames(pfi_data)=="race")] = "pfi_race"
  
  
  
  
  patient_count_survival = patient_count %>%
    dplyr::left_join(os_data, by = "barcode") %>%
    dplyr::left_join(dss_data, by = "barcode") %>%
    dplyr::left_join(dfi_data, by =  "barcode") %>%
    dplyr::left_join(pfi_data, by = "barcode") %>%
    dplyr::select(barcode, age, gender,
                  os_race, OS, OS.time,
                  dss_race, DSS, DSS.time,
                  dfi_race, DFI, DFI.time,
                  pfi_race, PFI, PFI.time,
                  one_of(interest_variable_info))
  

  
  write.table(patient_count_survival, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  return(patient_count_survival)
  
}



piu_counts_cdr_clinical_unite = function(piu_count_filename,
                                         cdr_clinical,
                                         patient_sum_min,
                                         output_dir,
                                         output_name)
  
{
  
  piu_count_df = fread(piu_count_filename,
                       stringsAsFactors = F)
  
  which_count = grep("TCGA", colnames(piu_count_df))
  non_zero = apply(piu_count_df[,..which_count],1,nonzero)
  piu_count_fsel = piu_count_df[which(non_zero>=patient_sum_min),]
  
  
  piu_count_sel = piu_count_fsel%>%
    dplyr::mutate(piu_info = paste(uniprot_accession, start_position,end_position,unit_label, unit_name, gene_name, gene_id, sep = "_"))%>%
    dplyr::mutate(gene_info = paste(gene_name, gene_id,sep = "_"))%>%
    dplyr::select(uniprot_accession, start_position, end_position, unit_name, gene_name, gene_id, unit_label,gene_info, piu_info, row_sum,
                  everything())
  

  if(nrow(piu_count_sel)>0)
  {
    which_barcode = grep("TCGA", colnames(piu_count_sel))
    
    piu_matrix = t(as.matrix(piu_count_sel[,which_barcode]))
    
    piu_count = data.frame(barcode = colnames(piu_count_sel)[which_barcode],
                           piu_matrix,
                           stringsAsFactors = F)
    
    colnames(piu_count) = c("barcode", piu_count_sel$piu_info)
    rownames(piu_count) = NULL
    
    piu_gene_df = data.frame(piu_info = piu_count_sel$piu_info,
                             gene_info = piu_count_sel$gene_info,
                             stringsAsFactors = F)
    
    
    piu_clinical_unite_data = piu_count %>%
      dplyr::left_join(cdr_clinical, by = c("barcode" = "bcr_patient_barcode")) %>%
      dplyr::select(barcode,colnames(cdr_clinical)[3:34],everything())
    
    
    return_list = vector(mode = "list", length = 4)
    
    return_list[[1]] = piu_clinical_unite_data
    return_list[[2]] = piu_gene_df
    return_list[[3]] = piu_count_sel$piu_info
    return_list[[4]] = piu_count_sel$gene_id
    
    
    write.table(return_list[[1]], paste0(output_dir, output_name),
                quote = F, row.names = F, sep = "\t")
    
    return(return_list)
    
  }else{
    cat("No item on this type of PIU.","\n")
  }
  
  
}



gene_counts_cdr_clinical_unite = function(gene_count_filename,
                                          cdr_clinical,
                                          patient_sum_min,
                                          output_dir,
                                          output_name)
  
{
  
  gene_count_df = fread(gene_count_filename,
                        stringsAsFactors = F)
  
  
  if("gene_info" %in% colnames(gene_count_df))
  {
    gene_count_df_info = as_tibble(gene_count_df)
  }else{
    gene_count_df_info = gene_count_df %>%
      dplyr::mutate(gene_info = paste(gene_name, gene_id, sep = "_"))
    
  }
  
  which_count = grep("TCGA", colnames(gene_count_df_info))
  non_zero = apply(gene_count_df_info[,which_count],1,nonzero)
  gene_count_sel = gene_count_df_info[which(non_zero>=patient_sum_min),]
  
  #piu_info = piu_count_sel$piu_info
  
  
  which_barcode = grep("TCGA", colnames(gene_count_sel))
  
  gene_matrix = t(as.matrix(gene_count_sel[,which_barcode]))
  
  gene_count = data.frame(barcode = colnames(gene_count_sel)[which_barcode],
                          gene_matrix,
                          stringsAsFactors = F)
  
  colnames(gene_count) = c("barcode", gene_count_sel$gene_info)
  rownames(gene_count) = NULL
  
  
  gene_clinical_unite_data = gene_count %>%
    dplyr::left_join(cdr_clinical, by = c("barcode" = "bcr_patient_barcode")) %>%
    dplyr::select(barcode,colnames(cdr_clinical)[3:34],everything())
  
  
  return_list = vector(mode = "list", length = 2)
  
  return_list[[1]] = gene_clinical_unite_data
  return_list[[2]] = gene_count_sel$gene_info
  
  
  write.table(return_list[[1]], paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
  
  return(return_list)
  
}






fit_survival_model_v2_test= function(surv_info_data,
                                     interest_variable_info,
                                     min_surv_time,
                                     min_surv_people,
                                     output_dir,
                                     output_name)
{
  paste0(mutation_type,"_npc_cdr_univariate_test.tsv")
  #   
  endpoint_flag = c(T,T)
  v_status = c("OS","PFI")
  v_time = c("OS.time","PFI.time")
  v_race = c("os_race","pfi_race")
  
  
  q_status = quos(OS,PFI)
  q_time = quos(OS.time,PFI.time)
  q_race = quos(os_race,pfi_race)
  
  
  list_surv_models = vector(mode = "list", length = 2)
  
  for(i in 1:2)
  {
    
    
    this_race = q_race[[i]]
    this_status = q_status[[i]]
    this_time = q_time[[i]]
    
    
    this_race_v = v_race[i]
    this_status_v = v_status[i]
    this_time_v = v_time[i]
    
    
    
    this_surv_data = surv_info_data %>%
      dplyr::select(barcode, age,gender, !!this_race, !!this_status, !!this_time) %>%
      na.omit()%>%
      dplyr::arrange(!!this_status, !!this_time) %>%
      dplyr::filter(!!this_time >= min_surv_time)
    
    if(nrow(this_surv_data) == 0)
    {
      endpoint_flag[i] = F
    }else{
      
      this_surv_table = as.data.frame(table(this_surv_data[[this_status_v]]))
      
      if(min(this_surv_table$Freq)< min_surv_people)
      {
        endpoint_flag[i] = F
      }else{
        this_age_data = as.numeric(this_surv_data$age)
        this_gender_data = relevel(as.factor(this_surv_data$gender), ref = unique(this_surv_data$gender)[1])
        this_race_data = relevel(as.factor(this_surv_data[[this_race_v]]), ref = unique(this_surv_data[[this_race_v]])[1])
        this_surv_object = Surv(time = this_surv_data[[this_time_v]], event = this_surv_data[[this_status_v]])
        
      }
    }
    
    
    this_surv_result_df = rbindlist(lapply(1:length(interest_variable_info), function(x)
    {
      
      this_count_df = surv_info_data %>%
        dplyr::select(barcode, one_of(interest_variable_info[x]))
      this_count = this_surv_data %>%
        dplyr::left_join(this_count_df, by = "barcode")
      this_num_patients = sum(this_count[,7]!= 0)
      this_total_patients = nrow(this_count)
      this_count_coeff = NA
      this_count_exp_coeff = NA
      this_count_pval = NA
      this_assum_test = NA
      
      if(endpoint_flag[i] == T)
      {
        if(length(unique(this_gender_data))>1)
        {
          if(length(unique(this_race_data))>1)
          {
            this_model = coxph(this_surv_object ~  this_age_data  + this_gender_data + this_race_data + 
                                 this_count[,7])   
          }else{
            this_model = coxph(this_surv_object ~  this_age_data  + this_gender_data +
                                 this_count[,7])   
          }
          
        }else{
          if(length(unique(this_race_data))>1)
          {
            this_model = coxph(this_surv_object ~  this_age_data  + this_race_data + 
                                 this_count[,7])   
          }else{
            this_model = coxph(this_surv_object ~  this_age_data +
                                 this_count[,7])   
          }
          
        }
        
        
        this_test = cox.zph(this_model)
        this_table = this_test$table
        
        
        
        
        this = summary(this_model)
        this_coef = this$coefficients
        r_this =  length(unique(this_gender_data)) +length(unique(this_race_data))
        this_count_coeff = this_coef[r_this,1]
        this_count_exp_coeff = this_coef[r_this,2]
        this_count_pval = this_coef[r_this,5]
        this_assum_test = this_table[r_this,3]
        
        
        
      }
      
      one_surv_result = data.frame(count_info = interest_variable_info[x],
                                   this_num_patients, this_total_patients,
                                   this_count_coeff, this_count_exp_coeff, this_count_pval,
                                   this_assum_test,
                                   this_count_qval = NA,
                                   stringsAsFactors = F)
      
      if(x%%100 == 0)
        cat(x, i, "\n")
      
      return(one_surv_result)
      
    }))
    
    
    
    if(length(interest_variable_info)>5)
    {
      if(endpoint_flag[i] == T)
      {
        if(length(this_surv_result_df$this_count_pval)>300 & min(this_surv_result_df$this_count_pval,na.rm = T)<0.05 & max(this_surv_result_df$this_count_pval,na.rm = T)>0.95)
        {
          this_qval = qvalue(this_surv_result_df$this_count_pval)
        }else{
          this_qval = qvalue(this_surv_result_df$this_count_pval, pi0 = 1)
          
        }
        this_surv_result_df$this_count_qval = this_qval$qvalues
        
      }
    }else{
      cat("Small unit size, no q-val estimation.", "\n")
      
    }
    
    
    old_colnames = colnames(this_surv_result_df)
    new_colnames = gsub("this", v_status[i], old_colnames)
    colnames(this_surv_result_df) = new_colnames
    
    
    list_surv_models[[i]] = this_surv_result_df
    
    
  }
  
  
  
  final_result = list_surv_models %>%
    Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="count_info"), .) %>%
    dplyr::arrange(desc(OS_num_patients)) %>%
    replace(is.na(.),"")
  
  
  write.table(final_result, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
}









fit_survival_model_no_gender_v2_test= function(surv_info_data,
                                               interest_variable_info,
                                               min_surv_time,
                                               min_surv_people,
                                               output_dir,
                                               output_name)
{

  
  endpoint_flag = c(T,T)
  v_status = c("OS","PFI")
  v_time = c("OS.time","PFI.time")
  v_race = c("os_race","pfi_race")
  
  
  q_status = quos(OS,PFI)
  q_time = quos(OS.time,PFI.time)
  q_race = quos(os_race,pfi_race)
  
  
  list_surv_models = vector(mode = "list", length = 2)
  
  for(i in 1:2)
  {
    # i = 1
    
    this_race = q_race[[i]]
    this_status = q_status[[i]]
    this_time = q_time[[i]]
    
    
    this_race_v = v_race[i]
    this_status_v = v_status[i]
    this_time_v = v_time[i]
    
    
    
    this_surv_data = surv_info_data %>%
      dplyr::select(barcode, age,gender, !!this_race, !!this_status, !!this_time) %>%
      na.omit()%>%
      dplyr::arrange(!!this_status, !!this_time) %>%
      dplyr::filter(!!this_time >= min_surv_time)
    
    if(nrow(this_surv_data) == 0)
    {
      endpoint_flag[i] = F
    }else{
      
      this_surv_table = as.data.frame(table(this_surv_data[[this_status_v]]))
      
      if(min(this_surv_table$Freq)< min_surv_people)
      {
        endpoint_flag[i] = F
      }else{
        this_age_data = as.numeric(this_surv_data$age)
        this_gender_data = relevel(as.factor(this_surv_data$gender), ref = unique(this_surv_data$gender)[1])
        this_race_data = relevel(as.factor(this_surv_data[[this_race_v]]), ref = unique(this_surv_data[[this_race_v]])[1])
        this_surv_object = Surv(time = this_surv_data[[this_time_v]], event = this_surv_data[[this_status_v]])
        
      }
    }
    
    
    this_surv_result_df = rbindlist(lapply(1:length(interest_variable_info), function(x)
    {
      # x = 1
      
      this_count_df = surv_info_data %>%
        dplyr::select(barcode, one_of(interest_variable_info[x]))
      this_count = this_surv_data %>%
        dplyr::left_join(this_count_df, by = "barcode")
      this_num_patients = sum(this_count[,7]!= 0)
      this_total_patients = nrow(this_count)
      this_count_coeff = NA
      this_count_exp_coeff = NA
      this_count_pval = NA
      this_assum_test = NA
      
      
      if(endpoint_flag[i] == T)
      {
        
        if(length(unique(this_race_data))>1)
        {
          this_model = coxph(this_surv_object ~  this_age_data + this_race_data + 
                               this_count[,7])   
          
          
        }else{
          this_model = coxph(this_surv_object ~  this_age_data +
                               this_count[,7])   
        }
        
        
        this_test = cox.zph(this_model)
        this_table = this_test$table
        
        this = summary(this_model)
        this_coef = this$coefficients
        r_this =  1 +length(unique(this_race_data))
        this_count_coeff = this_coef[r_this,1]
        this_count_exp_coeff = this_coef[r_this,2]
        this_count_pval = this_coef[r_this,5]
        this_assum_test = this_table[r_this,3]
        
        
      }
      
      one_surv_result = data.frame(count_info = interest_variable_info[x],
                                   this_num_patients, this_total_patients,
                                   this_count_coeff, this_count_exp_coeff, this_count_pval,
                                   this_assum_test,
                                   this_count_qval = NA,
                                   stringsAsFactors = F)
      
      if(x%%100 == 0)
        cat(x, i, "\n")
      
      return(one_surv_result)
      
    }))
    
    
    
    if(length(interest_variable_info)>5)
    {
      if(endpoint_flag[i] == T)
      {
        if(length(this_surv_result_df$this_count_pval)>300 & min(this_surv_result_df$this_count_pval,na.rm = T)<0.05 & max(this_surv_result_df$this_count_pval,na.rm = T)>0.95)
        {
          this_qval = qvalue(this_surv_result_df$this_count_pval)
        }else{
          this_qval = qvalue(this_surv_result_df$this_count_pval, pi0 = 1)
          
        }
        this_surv_result_df$this_count_qval = this_qval$qvalues
        
      }
    }else{
      
      cat("Small unit size, no q-val estimation.", "\n")
      
    }
    
    
    ### change column names of this df 
    old_colnames = colnames(this_surv_result_df)
    new_colnames = gsub("this", v_status[i], old_colnames)
    colnames(this_surv_result_df) = new_colnames
    
    
    list_surv_models[[i]] = this_surv_result_df
    
    
  }
  
  
  
  final_result = list_surv_models %>%
    Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="count_info"), .) %>%
    dplyr::arrange(desc(OS_num_patients)) %>%
    replace(is.na(.),"")
  
  
  write.table(final_result, paste0(output_dir, output_name),
              quote = F, row.names = F, sep = "\t")
  
  
}



interaction_cox_model_for_somatic_units= function(univariate_somatic_piu_results,
                                                  univariate_somatic_bpiu_results,
                                                  univariate_somatic_npc_results,
                                                  somatic_piu_filename,
                                                  somatic_bpiu_filename,
                                                  somatic_npc_filename,
                                                  germline_piu_filename,
                                                  germline_bpiu_filename,
                                                  germline_npc_filename,
                                                  cdr_clinical,
                                                  gender_as_covariate,
                                                  race_group_min = 6,
                                                  min_surv_days = 90,
                                                  min_surv_people = 5,
                                                  patient_sum_min = 3,
                                                  mutation_type = "interaction",
                                                  sig_level = 0.05,
                                                  q_pvalue,
                                                  q_status,
                                                  q_time,
                                                  q_race,
                                                  v_status,
                                                  v_time,
                                                  v_race,
                                                  output_dir = interaction_output_dir)

{
  
  somatic_piu_result = data.frame(unit_info = character(0),
                                  gene_info = character(0),
                                  unit_type = character(0),
                                  stringsAsFactors = F)
  if(file.exists(univariate_somatic_piu_results))
  {
    univariate_somatic_piu = fread(univariate_somatic_piu_results, stringsAsFactors = F)
    
    sig_piu = univariate_somatic_piu%>%
      dplyr::filter(!!q_pvalue <= sig_level)
    
    if(nrow(sig_piu)>0)
    {
      
      get_gene_info = unlist(lapply(1:nrow(sig_piu), function(x){
        
        
        this_break = unlist(strsplit(sig_piu$count_info[x], split = "_"))
        this_gene_info = paste(this_break[(length(this_break)-1):length(this_break)],
                               collapse = "_")
        
        return(this_gene_info)
        
      }))
      
      somatic_piu_result = data.frame(unit_info = sig_piu$count_info,
                                      gene_info = get_gene_info,
                                      unit_type = "PIU",
                                      stringsAsFactors = F)
      
      
      
      
    }
    
    
  }
  
  
  
  
  somatic_bpiu_result = data.frame(unit_info = character(0),
                                   gene_info = character(0),
                                   unit_type = character(0),
                                   stringsAsFactors = F)
  if(file.exists(univariate_somatic_bpiu_results))
  {
    univariate_somatic_bpiu = fread(univariate_somatic_bpiu_results, stringsAsFactors = F)
    
    sig_bpiu = univariate_somatic_bpiu%>%
      dplyr::filter(!!q_pvalue <= sig_level)
    
    if(nrow(sig_bpiu)>0)
    {
      
      get_gene_info = sig_bpiu$count_info
      
      
      somatic_bpiu_result = data.frame(unit_info = sig_bpiu$count_info,
                                       gene_info = get_gene_info,
                                       unit_type = "bPIU",
                                       stringsAsFactors = F)
      
      
    }
    
    
  }
  
  
  
  
  
  somatic_npc_result = data.frame(unit_info = character(0),
                                  gene_info = character(0),
                                  unit_type = character(0),
                                  stringsAsFactors = F)
  if(file.exists(univariate_somatic_npc_results))
  {
    univariate_somatic_npc = fread(univariate_somatic_npc_results, stringsAsFactors = F)
    
    sig_npc = univariate_somatic_npc%>%
      dplyr::filter(!!q_pvalue <= sig_level)
    
    if(nrow(sig_npc)>0)
    {
      
      get_gene_info = sig_npc$count_info
      
      
      somatic_npc_result = data.frame(unit_info = sig_npc$count_info,
                                      gene_info = get_gene_info,
                                      unit_type = "nPC",
                                      stringsAsFactors = F)
      
      
    }
    
    
  }
  
  
  somatic_sig_unit = rbind(somatic_piu_result, somatic_bpiu_result, somatic_npc_result)
  

  if(nrow(somatic_sig_unit)>0)
  {
    
    
    somatic_piu_count = fread(somatic_piu_filename, stringsAsFactors = F)
    
    piu_patient_name = grep("TCGA", colnames(somatic_piu_count), value = T)
    
    
    if(length(which(somatic_sig_unit$unit_type=="PIU"))>0)
    {
      grab_piu = somatic_piu_count%>%
        dplyr::mutate(piu_info = paste(uniprot_accession, start_position, end_position,
                                       unit_label,unit_name, gene_name, gene_id, sep = "_"))%>%
        dplyr::mutate(gene_info = paste(gene_name, gene_id, sep = "_"))%>%
        dplyr::filter(piu_info %in% somatic_piu_result$unit_info)%>%
        dplyr::mutate(unit_type = "PIU")%>%
        dplyr::select(piu_info, unit_type, gene_info, one_of(piu_patient_name), row_sum)%>%
        dplyr::rename(unit_info = piu_info)
      
    }else{
      grab_piu = data.frame(matrix(ncol = 4+length(piu_patient_name), nrow = 0), stringsAsFactors = F)
      colnames(grab_piu) = c("unit_info","unit_type", "gene_info", piu_patient_name, "row_sum")
    }
    
    somatic_bpiu_count = fread(somatic_bpiu_filename, stringsAsFactors = F)
    
    bpiu_patient_name = grep("TCGA", colnames(somatic_bpiu_count), value = T)
    
    if(length(which(somatic_sig_unit$unit_type=="bPIU"))>0)
    {
      grab_bpiu = somatic_bpiu_count%>%
        dplyr::mutate(gene_info  = paste(gene_name, gene_id, sep = "_"))%>%
        dplyr::filter(gene_info%in%somatic_bpiu_result$unit_info)%>%
        dplyr::mutate(unit_type = "bPIU")%>%
        dplyr::mutate(unit_info = gene_info)%>%
        dplyr::select(unit_info, unit_type, gene_info, one_of(bpiu_patient_name), row_sum)
    }else{
      grab_bpiu = data.frame(matrix(ncol = 4+length(bpiu_patient_name), nrow = 0), stringsAsFactors = F)
      colnames(grab_bpiu) = c("unit_info","unit_type", "gene_info", bpiu_patient_name, "row_sum")
    }
    
    
    
    somatic_npc_count = fread(somatic_npc_filename, stringsAsFactors = F)
    npc_patient_name = grep("TCGA", colnames(somatic_npc_count), value = T)
    
    if(length(which(somatic_sig_unit$unit_type=="nPC"))>0)
    {
      grab_npc = somatic_npc_count%>%
        dplyr::mutate(gene_info  = paste(gene_name, gene_id, sep = "_"))%>%
        dplyr::filter(gene_info%in%somatic_npc_result$unit_info)%>%
        dplyr::mutate(unit_type = "nPC")%>%
        dplyr::mutate(unit_info = gene_info)%>%
        dplyr::select(unit_info, unit_type, gene_info, one_of(npc_patient_name), row_sum)
    }else{
      grab_npc = data.frame(matrix(ncol = 4+length(npc_patient_name), nrow = 0), stringsAsFactors = F)
      colnames(grab_npc) = c("unit_info","unit_type", "gene_info", npc_patient_name, "row_sum")
    }
    
    grab_sig_somatic = bind_rows(grab_piu, grab_bpiu, grab_npc)
    
    write.table(grab_sig_somatic, paste0(output_dir, v_status,"_sig_somatic_counts.tsv"),
                sep = "\t", row.names = F, quote = F)
    
    
    gene_to_grab = unique(somatic_sig_unit$gene_info)
    
    
    germline_piu_count = fread(germline_piu_filename, stringsAsFactors = F)
    
    
    germline_piu_patient_name = grep("TCGA", colnames(germline_piu_count), value = T)
    
    
    if(nrow(grab_sig_somatic)>0)
    {
      grab_germline_piu = germline_piu_count%>%
        dplyr::mutate(piu_info = paste(uniprot_accession, start_position, end_position,
                                       unit_label,unit_name, gene_name, gene_id, sep = "_"))%>%
        dplyr::mutate(gene_info = paste(gene_name, gene_id, sep = "_"))%>%
        dplyr::filter(gene_info %in% gene_to_grab)%>%
        dplyr::mutate(unit_type = "PIU")%>%
        dplyr::select(piu_info, unit_type, gene_info, one_of(piu_patient_name), row_sum)%>%
        dplyr::rename(unit_info = piu_info)
      
    }else{
      grab_germline_piu = data.frame(matrix(ncol = 4+length(germline_piu_patient_name), nrow = 0), stringsAsFactors = F)
      colnames(grab_germline_piu) = c("unit_info","unit_type", "gene_info", germline_piu_patient_name, "row_sum")
    }
    
    
    ###
    
    germline_bpiu_count = fread(germline_bpiu_filename, stringsAsFactors = F)
    
    germline_bpiu_patient_name = grep("TCGA", colnames(germline_bpiu_count), value = T)
    
    if(nrow(grab_sig_somatic)>0)
    {
      grab_germline_bpiu = germline_bpiu_count%>%
        dplyr::mutate(gene_info  = paste(gene_name, gene_id, sep = "_"))%>%
        dplyr::filter(gene_info%in%gene_to_grab)%>%
        dplyr::mutate(unit_type = "bPIU")%>%
        dplyr::mutate(unit_info = gene_info)%>%
        dplyr::select(unit_info, unit_type, gene_info, one_of(bpiu_patient_name), row_sum)
    }else{
      grab_germline_bpiu = data.frame(matrix(ncol = 4+length(germline_bpiu_patient_name), nrow = 0), stringsAsFactors = F)
      colnames(grab_germline_bpiu) = c("unit_info","unit_type", "gene_info", germline_bpiu_patient_name, "row_sum")
    }
    
    ###
    
    
    germline_npc_count = fread(germline_npc_filename, stringsAsFactors = F)
    
    germline_npc_patient_name = grep("TCGA", colnames(germline_npc_count), value = T)
    
    if(nrow(grab_sig_somatic)>0)
    {
      grab_germline_npc = germline_npc_count%>%
        dplyr::mutate(gene_info  = paste(gene_name, gene_id, sep = "_"))%>%
        dplyr::filter(gene_info%in%gene_to_grab)%>%
        dplyr::mutate(unit_type = "nPC")%>%
        dplyr::mutate(unit_info = gene_info)%>%
        dplyr::select(unit_info, unit_type, gene_info, one_of(npc_patient_name), row_sum)
    }else{
      grab_germline_npc = data.frame(matrix(ncol = 4+length(germline_npc_patient_name), nrow = 0), stringsAsFactors = F)
      colnames(grab_germline_npc) = c("unit_info","unit_type", "gene_info", germline_npc_patient_name, "row_sum")
    }
    
    
    grab_need_germline = bind_rows(grab_germline_piu, grab_germline_bpiu, grab_germline_npc)
    
    
    
    write.table(grab_need_germline, paste0(output_dir, v_status,"_need_germline_counts.tsv"),
                sep = "\t", row.names = F, quote = F)
   
     
    
    
    
  }
}







