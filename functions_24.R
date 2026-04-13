# Predefine names of columns for different datasets

colnames_metrics=c("seed","nb_ind","nb_loci","nb_dim","mut_rate","mut_rate_MC","mig_rate","rec_rate","rec_rate_LB","rec_rate_CL","var_mut","fit_alt","ep_FL","var_pref","var_cost","nb_dmi","ep_dmi","initial_pref","initial_cost","origin","life_cycle","opt_one_DI","opt_one_DII","opt_one_DIII","opt_two_DI","opt_two_DII","opt_two_DIII","scenario","rep","generations","sex_ratio_1","fit_eco_f1","fit_dmi_f1","pref_f1","cost_f1","fit_f1","pheno_f1_DI","pheno_f1_DII","pheno_f1_DIII","fit_eco_m1","fit_dmi_m1","pref_m1","cost_m1","fit_m1","pheno_m1_DI","pheno_m1_DII","pheno_m1_DIII","sex_ratio_2","fit_eco_f2","fit_dmi_f2","pref_f2","cost_f2","fit_f2","pheno_f2_DI","pheno_f2_DII","pheno_f2_DIII","fit_eco_m2","fit_dmi_m2","pref_m2","cost_m2","fit_m2","pheno_m2_DI","pheno_m2_DII","pheno_m2_DIII","F1_m1f2","F1_m2f1","dmi_scheme")

colnames_metrics_bm=c("seed","nb_ind","nb_loci","nb_dim","mut_rate","mut_rate_MC","mu_back","mig_rate","rec_rate","rec_rate_LB","rec_rate_CL","var_mut","fit_alt","ep_FL","var_pref","var_cost","nb_dmi","ep_dmi","initial_pref","initial_cost","origin","life_cycle","opt_one_DI","opt_one_DII","opt_one_DIII","opt_two_DI","opt_two_DII","opt_two_DIII","scenario","rep","generations","sex_ratio_1","fit_eco_f1","fit_dmi_f1","pref_f1","cost_f1","fit_f1","pheno_f1_DI","pheno_f1_DII","pheno_f1_DIII","fit_eco_m1","fit_dmi_m1","pref_m1","cost_m1","fit_m1","pheno_m1_DI","pheno_m1_DII","pheno_m1_DIII","sex_ratio_2","fit_eco_f2","fit_dmi_f2","pref_f2","cost_f2","fit_f2","pheno_f2_DI","pheno_f2_DII","pheno_f2_DIII","fit_eco_m2","fit_dmi_m2","pref_m2","cost_m2","fit_m2","pheno_m2_DI","pheno_m2_DII","pheno_m2_DIII","F1_m1f2","F1_m2f1","dmi_scheme")

colnames_assay_raw=c("nb_ind","nb_loci","nb_dim","mut_rate","mut_rate_MC","mig_rate","rec_rate","rec_rate_LB","rec_rate_CL","var_mut","fit_alt","ep_FL","var_pref","var_cost","nb_dmi","ep_dmi","initial_pref","initial_cost","life_cycle","origin","opt_one_DI","opt_one_DII","opt_one_DIII","opt_two_DI","opt_two_DII","opt_two_DIII","scenario","rep","fit_alt_eval","ep_dmi_eval","same_conditions","barriers","mig_pop1","pheno_m1","eco_fit_m1","dmi_fit_m1","pref_m1","nb_offspring1","nb_mating1","eco_fi1","hyb_fit1","rel_fit1","mig_pop2","pheno_m2","eco_fit_m2","dmi_fit_m2","pref_m2","nb_offspring2","nb_mating2","eco_fi2","hyb_fit2","rel_fit2","lmk_1","t_lmk_1","fmk_1","t_fmk_1","lmk_2","t_lmk_2","fmk_2","t_fmk_2")

colnames_assay=c("nb_ind","nb_loci","nb_dim","mut_rate","mut_rate_MC","mig_rate","rec_rate","rec_rate_LB","rec_rate_CL","var_mut","fit_alt","ep_FL","var_pref","var_cost","nb_dmi","ep_dmi","initial_pref","initial_cost","origin","life_cycle","opt_one_DI","opt_one_DII","opt_one_DIII","opt_two_DI","opt_two_DII","opt_two_DIII","scenario","rep","fit_alt_eval","ep_dmi_eval","same_conditions","barriers","pheno_par_m1","eco_fit_par_m1","dmi_fit_par_m1","pref_par_m1","nb_offspring_m1_mean","nb_mating_m1_mean","eco_fit_m1_mean","hyb_fit_m1_mean","rel_fit_m1_mean","pheno_par_m2","eco_fit_par_m2","dmi_fit_par_m2","pref_par_m2","nb_offspring_m2_mean","nb_mating_m2_mean","eco_fit_m2_mean","hyb_fit_m2_mean","rel_fit_m2_mean","pheno_par_f1","eco_fit_par_f1","dmi_fit_par_f1","pref_par_f1","nb_offspring_f1_mean","nb_mating_f1_mean","eco_fit_f1_mean","hyb_fit_f1_mean","rel_fit_f1_mean","pheno_par_f2","eco_fit_par_f2","dmi_fit_par_f2","pref_par_f2","nb_offspring_f2_mean","nb_mating_f2_mean","eco_fit_f2_mean","hyb_fit_f2_mean","rel_fit_f2_mean","nb_seg_lmk_m1","nb_seg_lmk_m2","nb_seg_fmk_m1","nb_seg_fmk_m2","nb_seg_lmk_f1","nb_seg_lmk_f2","nb_seg_fmk_f1","nb_seg_fmk_f2","t_loss_lmk_m1","t_loss_lmk_m2","t_loss_fmk_m1","t_loss_fmk_m2","t_loss_lmk_f1","t_loss_lmk_f2","t_loss_fmk_f1","t_loss_fmk_f2","cor_ext_m1","cor_ext_m2","cor_ext_f1","cor_ext_f2","t_med_loss_lmk_m1","t_med_loss_lmk_m2","t_med_loss_fmk_m1","t_med_loss_fmk_m2","t_med_loss_lmk_f1","t_med_loss_lmk_f2","t_med_loss_fmk_f1","t_med_loss_fmk_f2","t_q95_loss_lmk_m1","t_q95_loss_lmk_m2","t_q95_loss_fmk_m1","t_q95_loss_fmk_m2","t_q95_loss_lmk_f1","t_q95_loss_lmk_f2","t_q95_loss_fmk_f1","t_q95_loss_fmk_f2","nb_trials","nb_paring_for_F1_metrics","mean_fit_eco_f1","mean_dmi_f1","mean_fit_f1","mean_pref_f1","mean_cost_f1","mean_pheno_f1_t0","mean_pheno_f1_t1","mean_pheno_f1_t2","mean_fit_eco_f2","mean_dmi_f2","mean_fit_f2","mean_pref_f2","mean_cost_f2","mean_pheno_f2_t0","mean_pheno_f2_t1","mean_pheno_f2_t2","mean_fit_eco_m1","mean_dmi_m1","mean_fit_m1","mean_pref_m1","mean_cost_m1","mean_pheno_m1_t0","mean_pheno_m1_t1","mean_pheno_m1_t2","mean_fit_eco_m2","mean_dmi_m2","mean_fit_m2","mean_pref_m2","mean_cost_m2","mean_pheno_m2_t0","mean_pheno_m2_t1","mean_pheno_m2_t2")

colnames_assay_geom=c("seed","nb_ind","nb_loci","nb_dim","mut_rate","mut_rate_MC","mig_rate","rec_rate","rec_rate_LB","rec_rate_CL","var_mut","fit_alt","ep_FL","var_pref","var_cost","nb_dmi","ep_dmi","initial_pref","initial_cost","origin","life_cycle","opt_one_DI","opt_one_DII","opt_one_DIII","opt_two_DI","opt_two_DII","opt_two_DIII","scenario","rep","fit_alt_eval","ep_dmi_eval","same_conditions","barriers","pheno_par_m1","eco_fit_par_m1","dmi_fit_par_m1","pref_par_m1","nb_offspring_m1_mean","nb_mating_m1_mean","eco_fit_m1_mean","hyb_fit_m1_mean","rel_fit_m1_mean","pheno_par_m2","eco_fit_par_m2","dmi_fit_par_m2","pref_par_m2","nb_offspring_m2_mean","nb_mating_m2_mean","eco_fit_m2_mean","hyb_fit_m2_mean","rel_fit_m2_mean","pheno_par_f1","eco_fit_par_f1","dmi_fit_par_f1","pref_par_f1","nb_offspring_f1_mean","nb_mating_f1_mean","eco_fit_f1_mean","hyb_fit_f1_mean","rel_fit_f1_mean","pheno_par_f2","eco_fit_par_f2","dmi_fit_par_f2","pref_par_f2","nb_offspring_f2_mean","nb_mating_f2_mean","eco_fit_f2_mean","hyb_fit_f2_mean","rel_fit_f2_mean","nb_seg_lmk_m1","nb_seg_lmk_m2","nb_seg_fmk_m1","nb_seg_fmk_m2","nb_seg_lmk_f1","nb_seg_lmk_f2","nb_seg_fmk_f1","nb_seg_fmk_f2","t_loss_lmk_m1","t_loss_lmk_m2","t_loss_fmk_m1","t_loss_fmk_m2","t_loss_lmk_f1","t_loss_lmk_f2","t_loss_fmk_f1","t_loss_fmk_f2","cor_ext_m1","cor_ext_m2","cor_ext_f1","cor_ext_f2","t_med_loss_lmk_m1","t_med_loss_lmk_m2","t_med_loss_fmk_m1","t_med_loss_fmk_m2","t_med_loss_lmk_f1","t_med_loss_lmk_f2","t_med_loss_fmk_f1","t_med_loss_fmk_f2","t_q95_loss_lmk_m1","t_q95_loss_lmk_m2","t_q95_loss_fmk_m1","t_q95_loss_fmk_m2","t_q95_loss_lmk_f1","t_q95_loss_lmk_f2","t_q95_loss_fmk_f1","t_q95_loss_fmk_f2","th_loss_lmk_m1","th_loss_lmk_m2","th_loss_fmk_m1","th_loss_fmk_m2","th_loss_lmk_f1","th_loss_lmk_f2","th_loss_fmk_f1","th_loss_fmk_f2","tg_loss_lmk_m1","tg_loss_lmk_m2","tg_loss_fmk_m1","tg_loss_fmk_m2","tg_loss_lmk_f1","tg_loss_lmk_f2","tg_loss_fmk_f1","tg_loss_fmk_f2","ie_lmk_m1","ie_lmk_m2","ie_lmk_f1","ie_lmk_f2","ie_fmk_m1","ie_fmk_m2","ie_fmk_f1","ie_fmk_f2","nb_trials","nb_paring_for_F1_metrics","mean_fit_eco_f1","mean_dmi_f1","mean_fit_f1","mean_pref_f1","mean_cost_f1","mean_pheno_f1_t0","mean_pheno_f1_t1","mean_pheno_f1_t2","mean_fit_eco_f2","mean_dmi_f2","mean_fit_f2","mean_pref_f2","mean_cost_f2","mean_pheno_f2_t0","mean_pheno_f2_t1","mean_pheno_f2_t2","mean_fit_eco_m1","mean_dmi_m1","mean_fit_m1","mean_pref_m1","mean_cost_m1","mean_pheno_m1_t0","mean_pheno_m1_t1","mean_pheno_m1_t2","mean_fit_eco_m2","mean_dmi_m2","mean_fit_m2","mean_pref_m2","mean_cost_m2","mean_pheno_m2_t0","mean_pheno_m2_t1","mean_pheno_m2_t2")

# create vectors for the list of the different scenarios, origins or life cycle, becasue they are coded as numbers in the simulations
list_scenario=c("all","LA+DMI","LA+MC","LA","MC+DMI","DMI","MC","neu")
list_origin=c("Ancestral","Diverged","Neutral SC","Diverged SC","Extension")
list_life_cycle=c("sel-mig","mig-sel")
list_dmi_scenario=c("default","neutral","network 100L","network 050L")

# Define list of origins
list_origin=c("Ancestral","Diverged","Neutral SC","Diverged SC","Extension","ExtensionSC","Specific_barrier_strength")

# Function that reads the metric output file of the simulation program and extract the parameters as well as the simulations outcome and add them to a dataframe.It also accounts for the type of DMIs scheme, and correct the generation times for extension of previous runs.
load_metrics=function(file_name,dataframe){
  data_test=read.table(file_name,nrows=1)
  rep=as.numeric(strsplit(strsplit(strsplit(file_name,"rep")[[1]][2],".txt")[[1]],"_")[[1]][1])
  
  if(regexpr("neuDMI",file_name)!=-1){
    dmi_sc=2
  } else {
    dmi_sc=1
  }
  if(regexpr("ntwk100L",file_name)!=-1){
    dmi_sc=3
  } else if(regexpr("ntwk50L",file_name)!=-1){
    dmi_sc=4
  }
  is_extension=FALSE
  if(regexpr("_ext.txt",file_name)!=-1){
    is_extension=TRUE  
  }
  is_extensionSC=FALSE
  if(regexpr("_extSC.txt",file_name)!=-1){
    is_extensionSC=TRUE
  }
  
  scenario=0
  if (data_test$V24<1&data_test$V34>0&data_test$V38==0){
    scenario=1
  } else if (data_test$V24<1&data_test$V34>0&data_test$V38==1){
    scenario=2
  }else if (data_test$V24<1&data_test$V34==0&data_test$V38==0){
    scenario=3
  }else if (data_test$V24<1&data_test$V34==0&data_test$V38==1){
    scenario=4
  }else if (data_test$V24==1&data_test$V34>0&data_test$V38==0){
    scenario=5
  }else if (data_test$V24==1&data_test$V34>0&data_test$V38==1){
    scenario=6
  }else if (data_test$V24==1&data_test$V34==0&data_test$V38==0){
    scenario=7
  }else if (data_test$V24==1&data_test$V34==0&data_test$V38==1){
    scenario=8
  }
  if (data_test$V8==1){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48)],NA,NA,data_test[50],NA,NA,scenario,rep)
  } else if (data_test$V8==2){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48,50)],NA,data_test[c(51,52)],scenario,rep)
  } else if (data_test$V8==3){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48,49,50)],data_test[c(52,53,54)],scenario,rep)
  }
  
  #temp_data=dataframe[,1:27]
  #temp_data=rbind(temp_data,unlist(param))
  #print(tail(duplicated(temp_data),1))
  # while (tail(duplicated(temp_data),1)){
  #  temp_data[dim(temp_data)[1],22]=as.numeric(temp_data[dim(temp_data)[1],22])+1
  #  param[[22]]=param[[22]]+1
  #}
  
  data_all=read.table(file_name,skip=1)
  data_all=data_all[data_all$V1 %% (param$V4/10)<2,]
  #data_all=cbind(rep,data_all)
  data_all=cbind(param,data_all)
  
  if (data_test$V8==1){
    colnames(data_all)=colnames_metrics[-c(37,38,45,46,54,55,62,63,66)-1]
    data_all$pheno_f1_DII=NA
    data_all$pheno_f1_DIII=NA
    data_all$pheno_m1_DII=NA
    data_all$pheno_m1_DIII=NA
    data_all$pheno_f2_DII=NA
    data_all$pheno_f2_DIII=NA
    data_all$pheno_m2_DII=NA
    data_all$pheno_m2_DIII=NA
  } else if (data_test$V8==2){
    colnames(data_all)=colnames_metrics[-c(38,46,55,63,66)-1]
    data_all$pheno_f1_DIII=NA
    data_all$pheno_m1_DIII=NA
    data_all$pheno_f2_DIII=NA
    data_all$pheno_m2_DIII=NA
  } else if (data_test$V8==3){
    colnames(data_all)=colnames_metrics[-66-1]
  }
  data_all$dmi_scheme=dmi_sc
  if(is_extension){
    data_all$generations=500000+data_all$generations
  }
  if(is_extensionSC){
    data_all$generations=500000+data_all$generations
    data_all$origin=5
  }

  rbind(dataframe,data_all)
}

# Function that reads the metric output file of the simulation program and extract the parameters as well as the simulations outcome and add them to a dataframe.It also accounts for the type of DMIs scheme, and correct the generation times for extension of previous runs. This functions handles the case where a different mutation rate for back mutation was possible.
load_metrics_bm=function(file_name,dataframe){
  data_test=read.table(file_name,nrows=1)
  rep=as.numeric(strsplit(strsplit(strsplit(file_name,"rep")[[1]][2],".txt")[[1]],"_")[[1]][1])
  
  if(regexpr("neuDMI",file_name)!=-1){
    dmi_sc=2
  } else {
    dmi_sc=1
  }
  if(regexpr("ntwk100L",file_name)!=-1){
    dmi_sc=3
  } else if(regexpr("ntwk50L",file_name)!=-1){
    dmi_sc=4
  }
  is_extension=FALSE
  if(regexpr("_ext.txt",file_name)!=-1){
    is_extension=TRUE  
  }
  
  scenario=0
  if (data_test$V26<1&data_test$V36>0&data_test$V40==0){
    scenario=1
  } else if (data_test$V26<1&data_test$V36>0&data_test$V40==1){
    scenario=2
  }else if (data_test$V26<1&data_test$V36==0&data_test$V40==0){
    scenario=3
  }else if (data_test$V26<1&data_test$V36==0&data_test$V40==1){
    scenario=4
  }else if (data_test$V26==1&data_test$V36>0&data_test$V40==0){
    scenario=5
  }else if (data_test$V26==1&data_test$V36>0&data_test$V40==1){
    scenario=6
  }else if (data_test$V26==1&data_test$V36==0&data_test$V40==0){
    scenario=7
  }else if (data_test$V26==1&data_test$V36==0&data_test$V40==1){
    scenario=8
  }
  if (data_test$V8==1){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,42,44,46,48,50)],NA,NA,data_test[52],NA,NA,scenario,rep)
  } else if (data_test$V8==2){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,42,44,46,48,50,51)],NA,data_test[c(53,54)],scenario,rep)
  } else if (data_test$V8==3){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,42,44,46,48,50,51,52)],data_test[c(54,55,56)],scenario,rep)
  }
  
  #temp_data=dataframe[,1:27]
  #temp_data=rbind(temp_data,unlist(param))
  #print(tail(duplicated(temp_data),1))
  # while (tail(duplicated(temp_data),1)){
  #  temp_data[dim(temp_data)[1],22]=as.numeric(temp_data[dim(temp_data)[1],22])+1
  #  param[[22]]=param[[22]]+1
  #}
  
  data_all=read.table(file_name,skip=1)
  data_all=data_all[data_all$V1 %% (param$V4/10)<2,]
  #data_all=cbind(rep,data_all)
  data_all=cbind(param,data_all)
  
  if (data_test$V8==1){
    colnames(data_all)=colnames_metrics_bm[-c(38,39,46,47,55,56,63,64,67)-1]
    data_all$pheno_f1_DII=NA
    data_all$pheno_f1_DIII=NA
    data_all$pheno_m1_DII=NA
    data_all$pheno_m1_DIII=NA
    data_all$pheno_f2_DII=NA
    data_all$pheno_f2_DIII=NA
    data_all$pheno_m2_DII=NA
    data_all$pheno_m2_DIII=NA
  } else if (data_test$V8==2){
    colnames(data_all)=colnames_metrics_bm[-c(39,47,56,64,67)-1]
    data_all$pheno_f1_DIII=NA
    data_all$pheno_m1_DIII=NA
    data_all$pheno_f2_DIII=NA
    data_all$pheno_m2_DIII=NA
  } else if (data_test$V8==3){
    colnames(data_all)=colnames_metrics_bm[-67-1]
  }
  data_all$dmi_scheme=dmi_sc
  if(is_extension){
    data_all$generations=500000+data_all$generations
  }
  rbind(dataframe,data_all)
}

#Function that reads the output file of the RI measure program and extract the parameters as well as the simulations outcome and add them to a dataframe.
load_RI_measure_raw=function(file_name,dataframe){
  if(FALSE){
    file_name="./../multiloci_24/output_Rwkm/pop_LC_0_LA0.1_MC1_DMI0.1_rep4_assays_ED0.1_PZ_0_MC_0.txt"
    
  }
  data_test=read.table(file_name,nrows=1)
  rep=as.numeric(strsplit(strsplit(strsplit(file_name,"rep")[[1]][2],".txt")[[1]],"_")[[1]][1])
  scenario=0
  if (data_test$V24<1&data_test$V34>0&data_test$V38==0){
    scenario=1
  } else if (data_test$V24<1&data_test$V34>0&data_test$V38==1){
    scenario=2
  }else if (data_test$V24<1&data_test$V34==0&data_test$V38==0){
    scenario=3
  }else if (data_test$V24<1&data_test$V34==0&data_test$V38==1){
    scenario=4
  }else if (data_test$V24==1&data_test$V34>0&data_test$V38==0){
    scenario=5
  }else if (data_test$V24==1&data_test$V34>0&data_test$V38==1){
    scenario=6
  }else if (data_test$V24==1&data_test$V34==0&data_test$V38==0){
    scenario=7
  }else if (data_test$V24==1&data_test$V34==0&data_test$V38==1){
    scenario=8
  }
  if (data_test$V8==1){
    param=c(data_test[c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48)],NA,NA,data_test[50],NA,NA,scenario,rep)
  } else if (data_test$V8==2){
    param=c(data_test[c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48,49)],NA,data_test[c(51,52)],scenario,rep)
  } else if (data_test$V8==3){
    param=c(data_test[c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48,49,50)],data_test[c(52,53,54)],scenario,rep)
  }
  
  
  data_assay=read.table(file_name,nrows=1,skip=1)
  
  if (data_test$V8==1){
    if(all(data_assay[c(3:18,20:23,25:33,37:42)]==data_test[c(3:18,20:23,25:33,45:50)])){
      same_conditions=TRUE
    } else {
      same_conditions=FALSE
    }
  } else if(data_test$V8==2){
    if(all(data_assay[c(3:18,20:23,25:33,37:44)]==data_test[c(3:18,20:23,25:33,45:52)])){
      same_conditions=TRUE
    } else {
      same_conditions=FALSE
    }
  } else if(data_test$V8==3){
    if(all(data_assay[c(3:18,20:23,25:33,37:46)]==data_test[c(3:18,20:23,25:33,45:54)])){
      same_conditions=TRUE
    } else {
      same_conditions=FALSE
    }
  }
  param=c(param,data_assay$V24,data_assay$V34)
  
  data_test=data_assay
  scenario=0
  if (data_test$V24<1&data_test$V34>0&data_test$V36==0){
    scenario=1
  } else if (data_test$V24<1&data_test$V34>0&data_test$V36==1){
    scenario=2
  }else if (data_test$V24<1&data_test$V34==0&data_test$V36==0){
    scenario=3
  }else if (data_test$V24<1&data_test$V34==0&data_test$V36==1){
    scenario=4
  }else if (data_test$V24==1&data_test$V34>0&data_test$V36==0){
    scenario=5
  }else if (data_test$V24==1&data_test$V34>0&data_test$V36==1){
    scenario=6
  }else if (data_test$V24==1&data_test$V34==0&data_test$V36==0){
    scenario=7
  }else if (data_test$V24==1&data_test$V34==0&data_test$V36==1){
    scenario=8
  }
  
  
  #print(paste(c(list_scenario[scenario],data_test$V22,data_test$V32,data_test$V34,file_name)),collapse="  ", sep="  ")
  param=unlist(c(param,same_conditions,scenario)) ## colnames_assay[2:30]
  #param=cbind(param,list_scenario[scenario])
  
  data_test=read.table(file_name,skip=2,nrows=1)
  test22=read.table(file_name,skip=3,flush=TRUE,fill=NA)
  test22=as.data.frame(test22)
  #test22=cbind(rep(i,length(test22[,1])),test22)
  if (data_assay$V8==1){
    colnames(test22)=c("mig_pop1","pheno_m1","eco_fit_m1","dmi_fit_m1","pref_m1","nb_offspring1","nb_mating1","eco_fi1","hyb_fit1","rel_fit1","mig_pop2","pheno_m2","eco_fit_m2","dmi_fit_m2","pref_m2","nb_offspring2","nb_mating2","eco_fi2","hyb_fit2","rel_fit2","lmk_1","t_lmk_1","fmk_1","t_fmk_1","lmk_2","t_lmk_2","fmk_2","t_fmk_2")
  }
  #test_all=rbind(test_all,test22)
  #hist(test22$nb_offspring1/test22$nb_mating1,breaks=seq(0,1,.005))
  temp_d=cbind(matrix(param,nrow=1,ncol=length(param),byrow = TRUE),test22)
  colnames(temp_d)=colnames_assay_raw
  dataframe=rbind(dataframe,temp_d)
  return(dataframe)
}

# Function that reads the metric output file of the simulation program and extract the parameters as well as the simulations outcome and add them to a dataframe.It also accounts for the type of DMIs scheme, and correct the generation times for extension of previous runs. This functions handles the case where the phenotype for mate choice was different than for fitness.
load_metrics_mc=function(file_name,dataframe){
  data_test=read.table(file_name,nrows=1)
  rep=as.numeric(strsplit(strsplit(strsplit(file_name,"rep")[[1]][2],".txt")[[1]],"_")[[1]][1])
  
  if(regexpr("neuDMI",file_name)!=-1){
    dmi_sc=2
  } else {
    dmi_sc=1
  }
  if(regexpr("ntwk100L",file_name)!=-1){
    dmi_sc=3
  } else if(regexpr("ntwk50L",file_name)!=-1){
    dmi_sc=4
  }
  is_extension=FALSE
  if(regexpr("_ext.txt",file_name)!=-1){
    is_extension=TRUE  
  }
  is_extensionSC=FALSE
  if(regexpr("_extSC.txt",file_name)!=-1){
    is_extensionSC=TRUE
  }
  
  scenario=0
  if (data_test$V24<1&data_test$V34>0&data_test$V38==0){
    scenario=1
  } else if (data_test$V24<1&data_test$V34>0&data_test$V38==1){
    scenario=2
  }else if (data_test$V24<1&data_test$V34==0&data_test$V38==0){
    scenario=3
  }else if (data_test$V24<1&data_test$V34==0&data_test$V38==1){
    scenario=4
  }else if (data_test$V24==1&data_test$V34>0&data_test$V38==0){
    scenario=5
  }else if (data_test$V24==1&data_test$V34>0&data_test$V38==1){
    scenario=6
  }else if (data_test$V24==1&data_test$V34==0&data_test$V38==0){
    scenario=7
  }else if (data_test$V24==1&data_test$V34==0&data_test$V38==1){
    scenario=8
  }
  if (data_test$V8==2){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48)],NA,NA,data_test[50],NA,NA,scenario,rep)
  } else if (data_test$V8==3){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48,50)],NA,data_test[c(51,52)],scenario,rep)
  } else if (data_test$V8==4){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48,49,50)],data_test[c(52,53,54)],scenario,rep)
  }
  
  #temp_data=dataframe[,1:27]
  #temp_data=rbind(temp_data,unlist(param))
  #print(tail(duplicated(temp_data),1))
  # while (tail(duplicated(temp_data),1)){
  #  temp_data[dim(temp_data)[1],22]=as.numeric(temp_data[dim(temp_data)[1],22])+1
  #  param[[22]]=param[[22]]+1
  #}
  
  data_all=read.table(file_name,skip=1)
  data_all=data_all[data_all$V1 %% (param$V4/10)<2,]
  #data_all=cbind(rep,data_all)
  data_all=cbind(param,data_all)
  
  if (data_test$V8==1){
    colnames(data_all)=colnames_metrics[-c(37,38,45,46,54,55,62,63,66)-1]
    data_all$pheno_f1_DII=NA
    data_all$pheno_f1_DIII=NA
    data_all$pheno_m1_DII=NA
    data_all$pheno_m1_DIII=NA
    data_all$pheno_f2_DII=NA
    data_all$pheno_f2_DIII=NA
    data_all$pheno_m2_DII=NA
    data_all$pheno_m2_DIII=NA
  } else if (data_test$V8==2){
    colnames(data_all)=colnames_metrics[-c(38,46,55,63,66)-1]
    data_all$pheno_f1_DIII=NA
    data_all$pheno_m1_DIII=NA
    data_all$pheno_f2_DIII=NA
    data_all$pheno_m2_DIII=NA
  } else if (data_test$V8==3){
    colnames(data_all)=colnames_metrics[-66-1]
  }
  data_all$dmi_scheme=dmi_sc
  if(is_extension){
    data_all$generations=500000+data_all$generations
  }
  if(is_extensionSC){
    data_all$generations=500000+data_all$generations
    data_all$origin=5
  }
  rbind(dataframe,data_all)
}

# Function that reads the output file of the RI measure program and extract the parameters, and compute the mean of the loss time and the invasion probobability and add them to a dataframe.
load_RI_measure=function(file_name,dataframe){
  if(FALSE){
    file_name="./../multiloci_24/output_default/pop_LC_0_LA0.1_MC1_DMI0.1_rep9_assays_ED0.1_PZ_0_MC_0.txt"
    #file_name="output_LC1/pop_LC_1_LA0.1_MC1_DMI0.1_rep9_assays_ED0.1_PZ_0_MC_1.txt"
  }
  data_test=read.table(file_name,nrows=1)
  rep=as.numeric(strsplit(strsplit(strsplit(file_name,"rep")[[1]][2],".txt")[[1]],"_")[[1]][1])
  scenario=0
  if (data_test$V24<1&data_test$V34>0&data_test$V38==0){
    scenario=1
  } else if (data_test$V24<1&data_test$V34>0&data_test$V38==1){
    scenario=2
  }else if (data_test$V24<1&data_test$V34==0&data_test$V38==0){
    scenario=3
  }else if (data_test$V24<1&data_test$V34==0&data_test$V38==1){
    scenario=4
  }else if (data_test$V24==1&data_test$V34>0&data_test$V38==0){
    scenario=5
  }else if (data_test$V24==1&data_test$V34>0&data_test$V38==1){
    scenario=6
  }else if (data_test$V24==1&data_test$V34==0&data_test$V38==0){
    scenario=7
  }else if (data_test$V24==1&data_test$V34==0&data_test$V38==1){
    scenario=8
  }
  if (data_test$V8==1){
    param=c(data_test[c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48)],NA,NA,data_test[50],NA,NA,scenario,rep)
  } else if (data_test$V8==2){
    param=c(data_test[c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48,49)],NA,data_test[c(51,52)],scenario,rep)
  } else if (data_test$V8==3){
    param=c(data_test[c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48,49,50)],data_test[c(52,53,54)],scenario,rep)
  }
  
  
  data_assay=read.table(file_name,nrows=1,skip=1)
  
  if (data_test$V8==1){
    if(all(data_assay[c(3:18,20:23,25:33,37:42)]==data_test[c(3:18,20:23,25:33,45:50)])){
      same_conditions=TRUE
    } else {
      same_conditions=FALSE
    }
  } else if(data_test$V8==2){
    if(all(data_assay[c(3:18,20:23,25:33,37:44)]==data_test[c(3:18,20:23,25:33,45:52)])){
      same_conditions=TRUE
    } else {
      same_conditions=FALSE
    }
  } else if(data_test$V8==3){
    if(all(data_assay[c(3:18,20:23,25:33,37:46)]==data_test[c(3:18,20:23,25:33,45:54)])){
      same_conditions=TRUE
    } else {
      same_conditions=FALSE
    }
  }
  
  param=c(param,data_assay$V24,data_assay$V34)
  
  
  data_test=data_assay
  scenario=0
  if (data_test$V24<1&data_test$V34>0&data_test$V36==0){
    scenario=1
  } else if (data_test$V24<1&data_test$V34>0&data_test$V36==1){
    scenario=2
  }else if (data_test$V24<1&data_test$V34==0&data_test$V36==0){
    scenario=3
  }else if (data_test$V24<1&data_test$V34==0&data_test$V36==1){
    scenario=4
  }else if (data_test$V24==1&data_test$V34>0&data_test$V36==0){
    scenario=5
  }else if (data_test$V24==1&data_test$V34>0&data_test$V36==1){
    scenario=6
  }else if (data_test$V24==1&data_test$V34==0&data_test$V36==0){
    scenario=7
  }else if (data_test$V24==1&data_test$V34==0&data_test$V36==1){
    scenario=8
  }
  
  
  #print(paste(c(list_scenario[scenario],data_test$V22,data_test$V32,data_test$V34,file_name)),collapse="  ", sep="  ")
  param=unlist(c(param,same_conditions,scenario)) ## colnames_assay[2:30]
  #param=cbind(param,list_scenario[scenario])
  
  data_test=read.table(file_name,skip=2,nrows=1)
  test22=read.table(file_name,skip=3,flush=TRUE,fill=NA)
  test22=as.data.frame(test22)
  #test22=cbind(rep(i,length(test22[,1])),test22)
  if (data_assay$V8==1){
    colnames(test22)=c("mig_pop1","pheno_m1","eco_fit_m1","dmi_fit_m1","pref_m1","nb_offspring1","nb_mating1","eco_fi1","hyb_fit1","rel_fit1","mig_pop2","pheno_m2","eco_fit_m2","dmi_fit_m2","pref_m2","nb_offspring2","nb_mating2","eco_fi2","hyb_fit2","rel_fit2","lmk_1","t_lmk_1","fmk_1","t_fmk_1","lmk_2","t_lmk_2","fmk_2","t_fmk_2")
  }
  #test_all=rbind(test_all,test22)
  #hist(test22$nb_offspring1/test22$nb_mating1,breaks=seq(0,1,.005))
  test22$t_lmk_1[test22$t_lmk_1==0]=param[1]
  test22$t_lmk_2[test22$t_lmk_2==0]=param[1]
  test22$t_fmk_1[test22$t_fmk_1==0]=param[1]
  test22$t_fmk_2[test22$t_fmk_2==0]=param[1]
  
  filter=test22$mig_pop1=="m1m"
  if (all(sapply(c(2:10,12:28),function(x) is.numeric(test22[,x])))){
    vec=c(colMeans(test22[filter,c(2:10,12:20)]),colMeans(test22[!filter,c(2:10,12:20)]),
          mean(test22$lmk_1[filter]>0),mean(test22$lmk_2[filter]>0),mean(test22$fmk_1[filter]>0),mean(test22$fmk_2[filter]>0),mean(test22$lmk_1[!filter]>0),mean(test22$lmk_2[!filter]>0),mean(test22$fmk_1[!filter]>0),mean(test22$fmk_2[!filter]>0),
          mean(test22$t_lmk_1[test22$t_lmk_1>0&filter]),mean(test22$t_lmk_2[test22$t_lmk_2>0&filter]),mean(test22$t_fmk_1[test22$t_fmk_1>0&filter]),mean(test22$t_fmk_2[test22$t_fmk_2>0&filter]),mean(test22$t_lmk_1[test22$t_lmk_1>0&!filter]),mean(test22$t_lmk_2[test22$t_lmk_2>0&!filter]),mean(test22$t_fmk_1[test22$t_fmk_1>0&!filter]),mean(test22$t_fmk_2[test22$t_fmk_2>0&!filter]),
          cor(test22$t_lmk_1[filter&test22$t_lmk_1>0&test22$t_fmk_1>0],test22$t_fmk_1[filter&test22$t_lmk_1>0&test22$t_fmk_1>0]),cor(test22$t_lmk_2[filter&test22$t_lmk_2>0&test22$t_fmk_2>0],test22$t_fmk_2[filter&test22$t_lmk_2>0&test22$t_fmk_2>0]),cor(test22$t_lmk_1[!filter&test22$t_lmk_1>0&test22$t_fmk_1>0],test22$t_fmk_1[!filter&test22$t_lmk_1>0&test22$t_fmk_1>0]),cor(test22$t_lmk_2[!filter&test22$t_lmk_2>0&test22$t_fmk_2>0],test22$t_fmk_2[!filter&test22$t_lmk_2>0&test22$t_fmk_2>0]))
    
    
    
    vec=c(vec,median(test22$t_lmk_1[filter]),median(test22$t_lmk_2[filter]),median(test22$t_fmk_1[filter]),median(test22$t_fmk_2[filter]),median(test22$t_lmk_1[!filter]),median(test22$t_lmk_2[!filter]),median(test22$t_fmk_1[!filter]),median(test22$t_fmk_2[!filter]))
    vec=c(vec,quantile(test22$t_lmk_1[filter],0.95),quantile(test22$t_lmk_2[filter],0.95),quantile(test22$t_fmk_1[filter],0.95),quantile(test22$t_fmk_2[filter],0.95),quantile(test22$t_lmk_1[!filter],0.95),quantile(test22$t_lmk_2[!filter],0.95),quantile(test22$t_fmk_1[!filter],0.95),quantile(test22$t_fmk_2[!filter],0.95))
  } else{
    print(file_name)
    #file.remove(file_name)
    vec=rep(NA,80)
  }
  param=c(param,unlist(vec))
  
  
  
  if (param[3]==1){
    param=c(param,unlist(c(data_test[seq(2,16,2)],NA,NA,data_test[seq(18,28,2)],NA,NA,data_test[seq(30,40,2)],NA,NA,data_test[seq(42,52,2)],NA,NA)))
  } else if (param[3]==2){
    param=c(param,unlist(c(data_test[seq(2,18,2)],NA,data_test[seq(20,32,2)],NA,data_test[seq(34,46,2)],NA,data_test[seq(48,60,2)],NA)))
  } else if (param[3]==3){
    param=c(param,unlist(data_test[c(F,T)]))
  }
  
  
  
  dataframe=rbind(dataframe,param)
  return(dataframe)
}


# Function that reads the output file of the RI measure program and extract the parameters as well as the simulations outcome and add them to a dataframe.It however uses the harmonic (despite the name) means instead of the arithmetic ones to compute the mean accross all replicates.
load_RI_measure_geom=function(file_name,dataframe){
  if(FALSE){
    file_name="./../multiloci_24/output_default/pop_LC_0_LA0.1_MC1_DMI0.1_rep9_assays_ED0.1_PZ_0_MC_0.txt"
    file_name="./../multiloci_24/output_default/pop_LC_0_LA0.1_MC1_DMI0.1_rep9_assays_ED0.1_PZ_0_MC_0.txt"
    file_name="./../multiloci_24/output_iso_50percent/pop_LC_0_LA0.5_MC1_DMI0_rep6_assays_ED1_PZ_0_MC_1.txt"
    file_name="./../multiloci_24/output_iso_10percent/pop_LC_0_LA0.1_MC1_DMI0_rep1_assays_ED0.1_PZ_0.045_MC_0.txt"
    file_name="./../multiloci_24/output_250DMIopt/pop_LC_0_LA0.1_MC1_DMI0.1_rep1_assays_ED0.1_PZ_0.1_MC_0.txt"
    #file_name="output_LC1/pop_LC_1_LA0.1_MC1_DMI0.1_rep9_assays_ED0.1_PZ_0_MC_1.txt"
  }
  if (length(unlist(strsplit(i,"partial")))>1){
    print(c("skipped",file_name))
    return(dataframe)
  }
  data_test=read.table(file_name,nrows=1)
  rep=as.numeric(strsplit(strsplit(strsplit(file_name,"rep")[[1]][2],".txt")[[1]],"_")[[1]][1])
  scenario=0
  if (data_test$V24<1&data_test$V34>0&data_test$V38==0){
    scenario=1
  } else if (data_test$V24<1&data_test$V34>0&data_test$V38==1){
    scenario=2
  }else if (data_test$V24<1&data_test$V34==0&data_test$V38==0){
    scenario=3
  }else if (data_test$V24<1&data_test$V34==0&data_test$V38==1){
    scenario=4
  }else if (data_test$V24==1&data_test$V34>0&data_test$V38==0){
    scenario=5
  }else if (data_test$V24==1&data_test$V34>0&data_test$V38==1){
    scenario=6
  }else if (data_test$V24==1&data_test$V34==0&data_test$V38==0){
    scenario=7
  }else if (data_test$V24==1&data_test$V34==0&data_test$V38==1){
    scenario=8
  }
  if (data_test$V8==1){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48)],NA,NA,data_test[50],NA,NA,scenario,rep)
  } else if (data_test$V8==2){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48,49)],NA,data_test[c(51,52)],scenario,rep)
  } else if (data_test$V8==3){
    param=c(data_test[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,40,42,44,46,48,49,50)],data_test[c(52,53,54)],scenario,rep)
  }
  
  
  data_assay=read.table(file_name,nrows=1,skip=1)
  
  if (data_test$V8==1){
    if(all(data_assay[c(3:12,15:18,20:23,25:33,37:42)]==data_test[c(3:12,15:18,20:23,25:33,45:50)])){
      same_conditions=TRUE
    } else {
      same_conditions=FALSE
    }
  } else if(data_test$V8==2){
    if(all(data_assay[c(3:12,15:18,20:23,25:33,37:44)]==data_test[c(3:12,15:18,20:23,25:33,45:52)])){
      same_conditions=TRUE
    } else {
      same_conditions=FALSE
    }
  } else if(data_test$V8==3){
    if(all(data_assay[c(3:12,15:18,20:23,25:33,37:46)]==data_test[c(3:12,15:18,20:23,25:33,45:54)])){
      same_conditions=TRUE
    } else {
      same_conditions=FALSE
    }
  }
  
  param=c(param,data_assay$V24,data_assay$V34)
  
  
  data_test=data_assay
  scenario=0
  if (data_test$V24<1&data_test$V34>0&data_test$V36==0){
    scenario=1
  } else if (data_test$V24<1&data_test$V34>0&data_test$V36==1){
    scenario=2
  }else if (data_test$V24<1&data_test$V34==0&data_test$V36==0){
    scenario=3
  }else if (data_test$V24<1&data_test$V34==0&data_test$V36==1){
    scenario=4
  }else if (data_test$V24==1&data_test$V34>0&data_test$V36==0){
    scenario=5
  }else if (data_test$V24==1&data_test$V34>0&data_test$V36==1){
    scenario=6
  }else if (data_test$V24==1&data_test$V34==0&data_test$V36==0){
    scenario=7
  }else if (data_test$V24==1&data_test$V34==0&data_test$V36==1){
    scenario=8
  }
  
  
  #print(paste(c(list_scenario[scenario],data_test$V22,data_test$V32,data_test$V34,file_name)),collapse="  ", sep="  ")
  param=unlist(c(param,same_conditions,scenario)) ## colnames_assay[2:30]
  #param=cbind(param,list_scenario[scenario])
  

  
  data_test=read.table(file_name,skip=2,nrows=1)
  test22=read.table(file_name,skip=3,flush=TRUE,fill=NA)
  test22=as.data.frame(test22)
  if (dim(test22)[1]!=2000){
    print(c(file_name,dim(test22)))
  }
  #test22=cbind(rep(i,length(test22[,1])),test22)
  if (data_assay$V8==1){
    colnames(test22)=c("mig_pop1","pheno_m1","eco_fit_m1","dmi_fit_m1","pref_m1","nb_offspring1","nb_mating1","eco_fi1","hyb_fit1","rel_fit1","mig_pop2","pheno_m2","eco_fit_m2","dmi_fit_m2","pref_m2","nb_offspring2","nb_mating2","eco_fi2","hyb_fit2","rel_fit2","lmk_1","t_lmk_1","fmk_1","t_fmk_1","lmk_2","t_lmk_2","fmk_2","t_fmk_2")
  }
  #test_all=rbind(test_all,test22)
  #hist(test22$nb_offspring1/test22$nb_mating1,breaks=seq(0,1,.005))
  test22$t_lmk_1[test22$t_lmk_1==0]=param[2]
  test22$t_lmk_2[test22$t_lmk_2==0]=param[2]
  test22$t_fmk_1[test22$t_fmk_1==0]=param[2]
  test22$t_fmk_2[test22$t_fmk_2==0]=param[2]
  
  filter=test22$mig_pop1=="m1m"
  if (all(sapply(c(2:10,12:28),function(x) is.numeric(test22[,x])))){
    vec=c(colMeans(test22[filter,c(2:10,12:20)]),colMeans(test22[!filter,c(2:10,12:20)]),
          mean(test22$lmk_1[filter]>0),mean(test22$lmk_2[filter]>0),mean(test22$fmk_1[filter]>0),mean(test22$fmk_2[filter]>0),mean(test22$lmk_1[!filter]>0),mean(test22$lmk_2[!filter]>0),mean(test22$fmk_1[!filter]>0),mean(test22$fmk_2[!filter]>0),
          mean(test22$t_lmk_1[test22$t_lmk_1>0&filter]),
          mean(test22$t_lmk_2[test22$t_lmk_2>0&filter]),
          mean(test22$t_fmk_1[test22$t_fmk_1>0&filter]),
          mean(test22$t_fmk_2[test22$t_fmk_2>0&filter]),
          mean(test22$t_lmk_1[test22$t_lmk_1>0&!filter]),
          mean(test22$t_lmk_2[test22$t_lmk_2>0&!filter]),
          mean(test22$t_fmk_1[test22$t_fmk_1>0&!filter]),
          mean(test22$t_fmk_2[test22$t_fmk_2>0&!filter]),
          cor(test22$t_lmk_1[filter&test22$t_lmk_1>0&test22$t_fmk_1>0],test22$t_fmk_1[filter&test22$t_lmk_1>0&test22$t_fmk_1>0]),cor(test22$t_lmk_2[filter&test22$t_lmk_2>0&test22$t_fmk_2>0],test22$t_fmk_2[filter&test22$t_lmk_2>0&test22$t_fmk_2>0]),cor(test22$t_lmk_1[!filter&test22$t_lmk_1>0&test22$t_fmk_1>0],test22$t_fmk_1[!filter&test22$t_lmk_1>0&test22$t_fmk_1>0]),cor(test22$t_lmk_2[!filter&test22$t_lmk_2>0&test22$t_fmk_2>0],test22$t_fmk_2[!filter&test22$t_lmk_2>0&test22$t_fmk_2>0]))
    
    
    
    vec=c(vec,median(test22$t_lmk_1[filter]),median(test22$t_lmk_2[filter]),median(test22$t_fmk_1[filter]),median(test22$t_fmk_2[filter]),median(test22$t_lmk_1[!filter]),median(test22$t_lmk_2[!filter]),median(test22$t_fmk_1[!filter]),median(test22$t_fmk_2[!filter]))
    vec=c(vec,quantile(test22$t_lmk_1[filter],0.95),quantile(test22$t_lmk_2[filter],0.95),quantile(test22$t_fmk_1[filter],0.95),quantile(test22$t_fmk_2[filter],0.95),quantile(test22$t_lmk_1[!filter],0.95),quantile(test22$t_lmk_2[!filter],0.95),quantile(test22$t_fmk_1[!filter],0.95),quantile(test22$t_fmk_2[!filter],0.95))
    vec=c(vec,1/mean(1/test22$t_lmk_1[test22$t_lmk_1>0&filter]),
          1/mean(1/test22$t_lmk_2[test22$t_lmk_2>0&filter]),
          1/mean(1/test22$t_fmk_1[test22$t_fmk_1>0&filter]),
          1/mean(1/test22$t_fmk_2[test22$t_fmk_2>0&filter]),
          1/mean(1/test22$t_lmk_1[test22$t_lmk_1>0&!filter]),
          1/mean(1/test22$t_lmk_2[test22$t_lmk_2>0&!filter]),
          1/mean(1/test22$t_fmk_1[test22$t_fmk_1>0&!filter]),
          1/mean(1/test22$t_fmk_2[test22$t_fmk_2>0&!filter]))
    vec=c(vec,exp(mean(log(test22$t_lmk_1[test22$t_lmk_1>0&filter]))),
    exp(mean(log(test22$t_lmk_2[test22$t_lmk_2>0&filter]))),
    exp(mean(log(test22$t_fmk_1[test22$t_fmk_1>0&filter]))),
    exp(mean(log(test22$t_fmk_2[test22$t_fmk_2>0&filter]))),
    exp(mean(log(test22$t_lmk_1[test22$t_lmk_1>0&!filter]))),
    exp(mean(log(test22$t_lmk_2[test22$t_lmk_2>0&!filter]))),
    exp(mean(log(test22$t_fmk_1[test22$t_fmk_1>0&!filter]))),
    exp(mean(log(test22$t_fmk_2[test22$t_fmk_2>0&!filter]))))
    vec=c(vec,mean(test22$t_lmk_1[filter]==1),mean(test22$t_lmk_2[filter]==1),mean(test22$t_lmk_1[!filter]==1),mean(test22$t_lmk_2[!filter]==1),mean(test22$t_fmk_1[filter]==1),mean(test22$t_fmk_2[filter]==1),mean(test22$t_fmk_1[!filter]==1),mean(test22$t_fmk_2[!filter]==1))
  } else{
    print(c("sthg wrong",file_name))
    #file.remove(file_name)
    vec=rep(NA,97)
  }
  param=c(param,unlist(vec))
  
  
  
  if (param[4]==1){
    param=c(param,unlist(c(data_test[seq(2,16,2)],NA,NA,data_test[seq(18,28,2)],NA,NA,data_test[seq(30,40,2)],NA,NA,data_test[seq(42,52,2)],NA,NA)))
  } else if (param[4]==2){
    param=c(param,unlist(c(data_test[seq(2,18,2)],NA,data_test[seq(20,32,2)],NA,data_test[seq(34,46,2)],NA,data_test[seq(48,60,2)],NA)))
  } else if (param[4]==3){
    param=c(param,unlist(data_test[c(F,T)]))
  }
  
  
 # if(!same_conditions){
  #  print(c("change condition",file_name))
  #}
  

  dataframe=rbind(dataframe,param)
  if(length(unlist(strsplit(file_name,"percent")))>1){
    dataframe$origin[length(dataframe$origin)]=6
    dataframe$life_cycle[length(dataframe$life_cycle)]=as.numeric(unlist(strsplit(unlist(strsplit(file_name,"LC_"))[2],"_"))[1])
  }
  return(dataframe)
}


# Function that reads the population output file of the simulation program and extract the distribution of phenotype of the 2 populations
check_LB_pheno_dist=function(name,dataset,factor,nb_dim=1,dim=1){
  for (i in 1:10){
    load_pop=read.table(paste(name,i,".txt",sep=""),fill=NA,skip=1,colClasses = "factor")
    load_pop1=load_pop[1:which(load_pop$V1=="Female_pop2"),]
    n1f=as.numeric(as.character(load_pop1[1,2]))
    n1m=as.numeric(as.character(load_pop1[2+2*n1f,2]))
    
    load_pop1=load_pop1[load_pop1[,6]!="" ,c(2,3+dim,3+1*(nb_dim+1)+dim,3+2*(nb_dim+1)+dim,3+3*(nb_dim+1)+dim,3+4*(nb_dim+1)+dim)]
    load_pop1=as.data.frame(lapply(lapply(load_pop1,as.character),as.numeric))
    load_pop1$pheno=rowSums(load_pop1[,-1])
    
    
    load_pop2=read.table(paste(name,i,".txt",sep=""),fill=NA,skip=which(load_pop$V1=="Female_pop2"),colClasses = "factor")
    n2f=as.numeric(as.character(load_pop2[1,2]))
    n2m=as.numeric(as.character(load_pop2[2+2*n2f,2]))
    load_pop2=load_pop2[load_pop2[,6]!="" ,c(2,3+dim,3+1*(nb_dim+1)+dim,3+2*(nb_dim+1)+dim,3+3*(nb_dim+1)+dim,3+4*(nb_dim+1)+dim)]
    load_pop2=as.data.frame(lapply(lapply(load_pop2,as.character),as.numeric))
    load_pop2$pheno=rowSums(load_pop2[,-1])
    
    
    break_ph=seq(-3,3,0.1)
    h1=hist(c(load_pop1[,2],load_pop2[,2]),breaks=break_ph,plot=FALSE)
    h2=hist(c(load_pop1[,2]),breaks=break_ph,plot=FALSE)
    
    for( j in seq(5)){
      test=kmeans(c(load_pop1[,1+j],load_pop2[,1+j]),centers = 2)
      dataset=rbind(dataset,c(i,j,factor,sort(test$centers),test$betweenss/test$totss))
    }
    
  }
  
  if(length(strsplit(file_name,"percent"))>1){
    param[19]=6
  }
  dataset
}

# Function to compute the relative sojourn time of mutations (for all 3 pithagorian means)
compute_relative_sojourn_time=function(assay_temp_v,merge_assay,barriers,background,metric="arith"){
  
  ref=c("neu","LA","MC","DMI","LA+MC","LA+DMI","MC+DMI","all")
  ref_mat=matrix(c(NA,"LA","MC","DMI","LA+MC","LA+DMI","MC+DMI","all",
                   NA,NA,NA,NA,"MC","DMI",NA,"MC+DMI",
                   NA,NA,NA,NA,"LA",NA,"DMI","LA+DMI",
                   NA,NA,NA,NA,NA,"LA","MC","LA+MC",
                   NA,NA,NA,NA,NA,NA,NA,"DMI",
                   NA,NA,NA,NA,NA,NA,NA,"MC",
                   NA,NA,NA,NA,NA,NA,NA,"LA",
                   NA,NA,NA,NA,NA,NA,NA,NA),nrow=8,byrow = TRUE)
  
  colnames(assay_temp_v)=c("rep","mig_rate","barriers","nb_trials","eco_fit","hyb_fit","nb_offspring","nb_mating","t_loss_linked","t_loss_free")
  assay_temp3=assay_temp_v[assay_temp_v$barriers==barriers,]
  assay_temp_bck=assay_temp_v[assay_temp_v$barriers==background,]
  if ( background=="neu"){
    if(metric=="arith"){
      assay_temp3$t_loss_linked=assay_temp3$t_loss_linked/mean(assay_temp_bck$t_loss_linked)
      assay_temp3$t_loss_free=assay_temp3$t_loss_free/mean(assay_temp_bck$t_loss_free)
    }
    if(metric=="harm"){
      assay_temp3$t_loss_linked=assay_temp3$t_loss_linked/(sum(assay_temp_bck$nb_trials)/sum(assay_temp_bck$nb_trials/assay_temp_bck$t_loss_linked))
      
      assay_temp3$t_loss_free=assay_temp3$t_loss_free/(sum(assay_temp_bck$nb_trials)/sum(assay_temp_bck$nb_trials/assay_temp_bck$t_loss_free))
    }
    if(metric=="geom"){
      assay_temp3$t_loss_linked=assay_temp3$t_loss_linked/(exp(sum(assay_temp_bck$nb_trials*log(assay_temp_bck$t_loss_linked))/sum(assay_temp_bck$nb_trials)))
      
      assay_temp3$t_loss_free=assay_temp3$t_loss_free/(exp(sum(assay_temp_bck$nb_trials*log(assay_temp_bck$t_loss_free))/sum(assay_temp_bck$nb_trials)))
    }
  }  else if( all(dim(assay_temp3)==dim(assay_temp_bck))){
    assay_temp_bck=assay_temp_bck[order(assay_temp_bck$rep),]
    assay_temp3=assay_temp3[order(assay_temp3$rep),]
    
    assay_temp3$t_loss_linked=assay_temp3$t_loss_linked/assay_temp_bck$t_loss_linked
    assay_temp3$t_loss_free=assay_temp3$t_loss_free/assay_temp_bck$t_loss_free
    
  } else {
    #print(c(assay_temp_v[,2][1],assay_temp_v$barriers[1],dim(assay_temp3),dim(assay_temp_bck)))
    #print("wrong number of replicate")
    return(merge_assay)
  }
    
    

    assay_temp3$fullbarriers=assay_temp3$barriers
    assay_temp3$barriers=ref_mat[which(ref==background),which(ref==barriers)]
    assay_temp3$background=background
    assay_temp3$cond_prob_mating=assay_temp3$nb_offspring/assay_temp3$nb_mating
    assay_temp3$uncond_prob_mating=assay_temp3$nb_offspring/((1000-assay_temp3$nb_offspring)*100+assay_temp3$nb_mating)
    
    #View(assay_temp3)
    
    if(!((barriers %in% c("LA","LA+DMI","LA+MC","all"))|(background %in% c("LA","LA+DMI","LA+MC")))){
      assay_temp3$eco_fit=NA
    }
    if(!((barriers %in% c("MC","MC+DMI","LA+MC","all"))|(background %in% c("MC","MC+DMI","LA+MC")))){
      assay_temp3$cond_prob_mating=NA
      assay_temp3$uncond_prob_mating=NA
    }
    if(!((barriers %in% c("DMI","LA+DMI","MC+DMI","all"))|(background %in% c("DMI","MC+DMI","LA+DMI")))){
      assay_temp3$hyb_fit =NA
    }
    
    merge_assay=rbind(merge_assay,assay_temp3)

  merge_assay
}

# Displayed function that illustrate the mean sojourn time for the different RI barriers, and in different background.
do_plot_decomposition=function(assay_all_temp,mean_type="arith"){
  if(mean_type=="arith"){
    name_y="Relative (arithmetic) mean sojourn time"
  } else if (mean_type=="harm"){
    name_y="Relative (harmonic) mean sojourn time"
  } else if (mean_type=="geom"){
    name_y="Relative (geoemtric) mean sojourn time"
  } 
  min_y=min(assay_all_temp$t_loss_free[assay_all_temp$mean==mean_type],assay_all_temp$t_loss_linked[assay_all_temp$mean==mean_type])
  max_y=max(assay_all_temp$t_loss_free[assay_all_temp$mean==mean_type],assay_all_temp$t_loss_linked[assay_all_temp$mean==mean_type])

  show(ggplot(assay_all_temp[assay_all_temp$background=="neu"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("Neutral background; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","red","blue","black","black","black"))+scale_color_manual(values=c("black","black","black","purple","red","blue"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1)))
  show(ggplot(assay_all_temp[assay_all_temp$background=="LA"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("LA background; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","red","blue","black","black","black"))+scale_color_manual(values=c("black","black","black","purple","red","blue"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1)))
  show(ggplot(assay_all_temp[assay_all_temp$background=="MC"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("MC background; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","red","blue","black","black","black"))+scale_color_manual(values=c("black","black","black","purple","red","blue"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1)))
  show(ggplot(assay_all_temp[assay_all_temp$background=="DMI"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("DMI background; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","red","blue","black","black","black"))+scale_color_manual(values=c("black","black","black","purple","red","blue"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1)))
  
  show(ggplot(assay_all_temp[assay_all_temp$background=="neu"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("Neutral background; unlinked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","red","blue","black","black","black"))+scale_color_manual(values=c("black","black","black","purple","red","blue"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1)))
  show(ggplot(assay_all_temp[assay_all_temp$background=="LA"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("LA background; unlinked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","red","blue","black","black","black"))+scale_color_manual(values=c("black","black","black","purple","red","blue"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1)))
  show(ggplot(assay_all_temp[assay_all_temp$background=="MC"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("MC background; unlinked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","red","blue","black","black","black"))+scale_color_manual(values=c("black","black","black","purple","red","blue"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1)))
  show(ggplot(assay_all_temp[assay_all_temp$background=="DMI"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("DMI background; unlinked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","red","blue","black","black","black"))+scale_color_manual(values=c("black","black","black","purple","red","blue"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1)))
}

# 
compute_relative_mean_sojourn_time=function(assay_temp,var="mig_rate",merge_tag=FALSE){
  assay_male_temp=c()
  assay_female_temp=c()
  assay_male_temp_harm=c()
  assay_female_temp_harm=c()
  assay_male_temp_geom=c()
  assay_female_temp_geom=c()
  assay_temp1=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_m1_mean","hyb_fit_m1_mean","nb_offspring_m1_mean","nb_mating_m1_mean","t_loss_lmk_m1","t_loss_fmk_m1")]
  assay_temp2=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_m2_mean","hyb_fit_m2_mean","nb_offspring_m2_mean","nb_mating_m2_mean","t_loss_lmk_m2","t_loss_fmk_m2")]
  assay_temp3=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_f1_mean","hyb_fit_f1_mean","nb_offspring_f1_mean","nb_mating_f1_mean","t_loss_lmk_f1","t_loss_fmk_f1")]
  assay_temp4=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_f2_mean","hyb_fit_f2_mean","nb_offspring_f2_mean","nb_mating_f2_mean","t_loss_lmk_f2","t_loss_fmk_f2")]
  
  assay_temp5=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_m1_mean","hyb_fit_m1_mean","nb_offspring_m1_mean","nb_mating_m1_mean","th_loss_lmk_m1","th_loss_fmk_m1")]
  assay_temp6=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_m2_mean","hyb_fit_m2_mean","nb_offspring_m2_mean","nb_mating_m2_mean","th_loss_lmk_m2","th_loss_fmk_m2")]
  assay_temp7=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_f1_mean","hyb_fit_f1_mean","nb_offspring_f1_mean","nb_mating_f1_mean","th_loss_lmk_f1","th_loss_fmk_f1")]
  assay_temp8=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_f2_mean","hyb_fit_f2_mean","nb_offspring_f2_mean","nb_mating_f2_mean","th_loss_lmk_f2","th_loss_fmk_f2")]
  
  assay_temp9=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_m1_mean","hyb_fit_m1_mean","nb_offspring_m1_mean","nb_mating_m1_mean","tg_loss_lmk_m1","tg_loss_fmk_m1")]
  assay_temp10=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_m2_mean","hyb_fit_m2_mean","nb_offspring_m2_mean","nb_mating_m2_mean","tg_loss_lmk_m2","tg_loss_fmk_m2")]
  assay_temp11=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_f1_mean","hyb_fit_f1_mean","nb_offspring_f1_mean","nb_mating_f1_mean","tg_loss_lmk_f1","tg_loss_fmk_f1")]
  assay_temp12=assay_temp[,c("rep",var,"barriers","nb_trials","eco_fit_f2_mean","hyb_fit_f2_mean","nb_offspring_f2_mean","nb_mating_f2_mean","tg_loss_lmk_f2","tg_loss_fmk_f2")]
  
  vec_barrier=c("DMI","LA","MC","LA+DMI","LA+MC","MC+DMI","all","LA+DMI","LA+DMI","MC+DMI","MC+DMI","LA+MC","LA+MC","all","all","all","all","all","all")
  vec_background=c("neu","neu","neu","neu","neu","neu","neu","LA","DMI","MC","DMI","LA","MC","LA+MC","LA+DMI","MC+DMI","LA","MC","DMI")
  
  assay_all_temp=c()
  print(unique(assay_temp2[,2]))
  for (i in unique(assay_temp2[,2])){
    assay_male_temp=c()
    assay_female_temp=c()
    assay_male_temp_harm=c()
    assay_female_temp_harm=c()
    assay_male_temp_geom=c()
    assay_female_temp_geom=c()
    for(j in 1:length(vec_barrier)){
      assay_male_temp=compute_relative_sojourn_time(assay_temp1[assay_temp1[,2]==i,],assay_male_temp,vec_barrier[j],vec_background[j],metric = "arith")
      assay_male_temp=compute_relative_sojourn_time(assay_temp2[assay_temp2[,2]==i,],assay_male_temp,vec_barrier[j],vec_background[j],metric = "arith")
      
      assay_female_temp=compute_relative_sojourn_time(assay_temp3[assay_temp3[,2]==i,],assay_female_temp,vec_barrier[j],vec_background[j],metric = "arith")
      assay_female_temp=compute_relative_sojourn_time(assay_temp4[assay_temp4[,2]==i,],assay_female_temp,vec_barrier[j],vec_background[j],metric = "arith")
      
      assay_male_temp_harm=compute_relative_sojourn_time(assay_temp5[assay_temp5[,2]==i,],assay_male_temp_harm,vec_barrier[j],vec_background[j],metric = "harm")
      assay_male_temp_harm=compute_relative_sojourn_time(assay_temp6[assay_temp6[,2]==i,],assay_male_temp_harm,vec_barrier[j],vec_background[j],metric = "harm")
      
      assay_female_temp_harm=compute_relative_sojourn_time(assay_temp7[assay_temp7[,2]==i,],assay_female_temp_harm,vec_barrier[j],vec_background[j],metric = "harm")
      assay_female_temp_harm=compute_relative_sojourn_time(assay_temp8[assay_temp8[,2]==i,],assay_female_temp_harm,vec_barrier[j],vec_background[j],metric = "harm")
      
      assay_male_temp_geom=compute_relative_sojourn_time(assay_temp9[assay_temp9[,2]==i,],assay_male_temp_geom,vec_barrier[j],vec_background[j],metric = "geom")
      assay_male_temp_geom=compute_relative_sojourn_time(assay_temp10[assay_temp10[,2]==i,],assay_male_temp_geom,vec_barrier[j],vec_background[j],metric = "geom")
      
      assay_female_temp_geom=compute_relative_sojourn_time(assay_temp11[assay_temp11[,2]==i,],assay_female_temp_geom,vec_barrier[j],vec_background[j],metric = "geom")
      assay_female_temp_geom=compute_relative_sojourn_time(assay_temp12[assay_temp12[,2]==i,],assay_female_temp_geom,vec_barrier[j],vec_background[j],metric = "geom")
      
    }
    if(merge_tag){
      assay_male_temp$sex=paste("male lvl=",i)
      assay_female_temp$sex=paste("female lvl=",i)
      assay_male_temp_harm$sex=paste("male lvl=",i)
      assay_female_temp_harm$sex=paste("female lvl=",i)
      assay_male_temp_geom$sex=paste("male lvl=",i)
      assay_female_temp_geom$sex=paste("female lvl=",i)
    } else{
      
      assay_male_temp$sex="male"
      assay_male_temp_harm$sex="male"
      assay_male_temp_geom$sex="male"
      assay_female_temp$sex="female"
      assay_female_temp_harm$sex="female"
      assay_female_temp_geom$sex="female"
    }
    assay_male_temp$mean="arith"
    assay_female_temp$mean="arith"
    assay_male_temp_harm$mean="harm"
    assay_female_temp_harm$mean="harm"
    assay_male_temp_geom$mean="geom"
    assay_female_temp_geom$mean="geom"
    
    if(is.null(assay_all_temp)){
      assay_all_temp=rbind(assay_male_temp,assay_female_temp)
    } else {
      assay_all_temp=rbind(assay_all_temp,assay_male_temp)
      assay_all_temp=rbind(assay_all_temp,assay_female_temp)
    }
    assay_all_temp=rbind(assay_all_temp,assay_male_temp_harm)
    assay_all_temp=rbind(assay_all_temp,assay_female_temp_harm)
    assay_all_temp=rbind(assay_all_temp,assay_female_temp_geom)
    assay_all_temp=rbind(assay_all_temp,assay_male_temp_geom)

    }

  assay_all_temp$mean_sex=as.factor(paste(assay_all_temp$sex,assay_all_temp$mean))
  
  assay_all_temp$fullbarriers=factor(as.character(assay_all_temp$fullbarriers),levels=c("neu","LA","MC","DMI","LA+MC","LA+DMI","MC+DMI","all"))
  assay_all_temp
}

do_plot_decomposition_migpanel_onlyNeu=function(assay_all_temp,mean_type="arith",ncol=3){
  if(mean_type=="arith"){
    name_y="Relative (arithmetic) mean sojourn time"
  } else if (mean_type=="harm"){
    name_y="Relative (harmonic) mean sojourn time"
  } else if (mean_type=="geom"){
    name_y="Relative (geoemtric) mean sojourn time"
  } 
  min_y=min(assay_all_temp$t_loss_free[assay_all_temp$mean==mean_type],assay_all_temp$t_loss_linked[assay_all_temp$mean==mean_type])
  max_y=max(assay_all_temp$t_loss_free[assay_all_temp$mean==mean_type],assay_all_temp$t_loss_linked[assay_all_temp$mean==mean_type])
  
  show(ggplot(assay_all_temp[assay_all_temp$life_cycle=="sel-mig"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("A/ Adult dispersal; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  show(ggplot(assay_all_temp[assay_all_temp$life_cycle=="sel-mig"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("A/ Adult dispersal")+theme()+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol)+theme_bw()+ theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  show(ggplot(assay_all_temp[assay_all_temp$life_cycle=="mig-sel"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("B/ Juvenile dispersal; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  show(ggplot(assay_all_temp[assay_all_temp$life_cycle=="mig-sel"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("B/ Juvenile dispersal")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
}

do_plot_decomposition_migpanel_onlyNeu_V2=function(assay_all_temp,mean_type="arith",ncol=3){
  if(mean_type=="arith"){
    name_y="Relative (arithmetic) mean sojourn time"
  } else if (mean_type=="harm"){
    name_y="Relative (harmonic) mean sojourn time"
  } else if (mean_type=="geom"){
    name_y="Relative (geoemtric) mean sojourn time"
  } 
  min_y=min(assay_all_temp$t_loss_free[assay_all_temp$mean==mean_type],assay_all_temp$t_loss_linked[assay_all_temp$mean==mean_type])
  max_y=max(assay_all_temp$t_loss_free[assay_all_temp$mean==mean_type],assay_all_temp$t_loss_linked[assay_all_temp$mean==mean_type])
  
  show(ggplot(assay_all_temp[assay_all_temp$life_cycle=="sel-mig"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.25)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("A/ Adult dispersal; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("purple4","darkorange4"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  show(ggplot(assay_all_temp[assay_all_temp$life_cycle=="sel-mig"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.25)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("A/ Adult dispersal")+theme()+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("purple4","darkorange4"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol)+theme_bw()+ theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  show(ggplot(assay_all_temp[assay_all_temp$life_cycle=="mig-sel"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.25)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("B/ Juvenile dispersal; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("purple4","darkorange4"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
  show(ggplot(assay_all_temp[assay_all_temp$life_cycle=="mig-sel"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.25)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("B/ Juvenile dispersal")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("purple4","darkorange4"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
}

do_plot_decomposition_migpanel=function(assay_all_temp,mean_type="arith",ncol=3){
  if(mean_type=="arith"){
    name_y="Relative (arithmetic) mean sojourn time"
  } else if (mean_type=="harm"){
    name_y="Relative (harmonic) mean sojourn time"
  } else if (mean_type=="geom"){
    name_y="Relative (geoemtric) mean sojourn time"
  } 
  min_y=min(assay_all_temp$t_loss_free[assay_all_temp$mean==mean_type],assay_all_temp$t_loss_linked[assay_all_temp$mean==mean_type])
  max_y=max(assay_all_temp$t_loss_free[assay_all_temp$mean==mean_type],assay_all_temp$t_loss_linked[assay_all_temp$mean==mean_type])
  
  show(ggplot(assay_all_temp[assay_all_temp$background=="neu"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("Neutral background; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol))
  show(ggplot(assay_all_temp[assay_all_temp$background=="LA"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("LA background; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol))
  show(ggplot(assay_all_temp[assay_all_temp$background=="MC"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("MC background; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol))
  show(ggplot(assay_all_temp[assay_all_temp$background=="DMI"&assay_all_temp$mean==mean_type,],aes(y=t_loss_linked,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("DMI background; linked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol))
  
  show(ggplot(assay_all_temp[assay_all_temp$background=="neu"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("Neutral background; unlinked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol))
  show(ggplot(assay_all_temp[assay_all_temp$background=="LA"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("LA background; unlinked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol))
  show(ggplot(assay_all_temp[assay_all_temp$background=="MC"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("MC background; unlinked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol))
  show(ggplot(assay_all_temp[assay_all_temp$background=="DMI"&assay_all_temp$mean==mean_type,],aes(y=t_loss_free,x=fullbarriers,fill=sex,col=sex))+geom_boxplot(outlier.shape = NA)+geom_point(position=position_jitterdodge(0.15),alpha=0.15)+coord_trans(y="log10")+ylab(name_y)+xlab("Genetic barriers")+ggtitle("DMI background; unlinked marker")+theme_bw()+theme(legend.position = "none")+scale_fill_manual(values=c("purple","orange"))+scale_color_manual(values=c("black","black"))+geom_abline(slope=0,intercept = 1,col="black",linetype=2)+scale_y_continuous(limits=c(min_y,max_y),breaks=c(0.01,0.05,0.2,0.5,1))+facet_wrap(~mig_rate,ncol = ncol))
}


read_genome=function(pop,i){
  vec=rep(0,500)
  vec[100:1]=as.numeric(unlist(strsplit(as.character(pop[i,1]),"")))
  vec[200:101]=as.numeric(unlist(strsplit(as.character(pop[i,2]),"")))
  vec[300:201]=as.numeric(unlist(strsplit(as.character(pop[i,3]),"")))
  vec[400:301]=as.numeric(unlist(strsplit(as.character(pop[i,4]),"")))
  vec[500:401]=as.numeric(unlist(strsplit(as.character(pop[i,5]),"")))
  vec
}

random_reflective_walk=function(var,n,start){
  step=rnorm(n,0,var)
  path=rep(start,n)
  for (i in 2:n){
      path[i]=path[i-1]+step[i]
      if (path[i]<0){path[i]=0}
      else if(path[i]>0.5){path[i]=0.5}
  }
  path
}

mutational_bias=function(mutation_file,nb_LB,opt=5){
  mutations=read.table(mutation_file)
  mutations=as.numeric(unlist(mutations))
  nb_loci=mutations[3]
  nb_dmi=mutations[5]
  sel=mutations[6:(6+nb_loci-1)]
  
  
  sel=cbind(rep(seq(1:nb_LB),each=nb_loci/nb_LB),sel)
  sel=as.data.frame(sel)
 
  
  sel_coef=sel$sel
  if (opt>0){
    one_option=which(sel_coef>0)[which(cumsum(2*sel_coef[sel_coef>0])<(opt))]
    for(i in seq(length(sel_coef))[-one_option]){
      if(sel_coef[i]>0&(sum(2*sel_coef[one_option])+2*sel_coef[i])<(opt)){
        one_option=c(one_option,i)
      }
    }
  }
  mut_av=sel_coef[-one_option]
  
  
  chosen1=sample(500,10000,replace=TRUE)
  chosen2=unlist(sapply(seq(10000), function(x) sample(seq(500)[-chosen1[x]],1)))
  
  c1=sum(2*sel_coef[one_option])-as.numeric(chosen1 %in% one_option)*sel_coef[chosen1]+as.numeric(!(chosen1 %in% one_option))*sel_coef[chosen1]
  c2=sum(2*sel_coef[one_option])-as.numeric(chosen1 %in% one_option)*sel_coef[chosen1]+as.numeric(!(chosen1 %in% one_option))*sel_coef[chosen1]-as.numeric(chosen2 %in% one_option)*sel_coef[chosen2]+as.numeric(!(chosen2 %in% one_option))*sel_coef[chosen2]                                     
  print(paste(c("Mean effect of 1 step:",round((mean(c1)-sum(2*sel_coef[one_option]))*1000)/1000,", of 2 steps:",round((mean(c2)-sum(2*sel_coef[one_option]))*1000)/1000,"\nProportion of mutations increasing z in 1 step:",round(mean(c1>sum(2*sel_coef[one_option]))*1000)/1000,", in 2 steps:",round(mean(c2>sum(2*sel_coef[one_option]))*1000)/1000),sep="",collapse=""))
  
  
  h1=hist( c1,breaks=seq(floor(min(c1,c2)*100)/100,ceiling(max(c1,c2)*100)/100,0.005),plot = FALSE)
  h2=hist( c2,breaks=seq(floor(min(c1,c2)*100)/100,ceiling(max(c1,c2)*100)/100,0.005),plot=FALSE)
  plot(h1,col="red",border="red",main="",xlab="Phenotype reachable through a single (or a pair of) mutation")
  plot(h2,col=alpha("blue",0.5),border=alpha("blue",0.5),add=TRUE)
  abline(v=opt,col="black",lwd=1.5)

}

display_dmi_map_v2=function(mutation_file,nb_LB,opt=5){
  mutations=read.table(mutation_file)
  mutations=as.numeric(unlist(mutations))
  nb_loci=mutations[3]
  nb_dmi=mutations[5]
  sel=mutations[6:(6+nb_loci-1)]
  dmi=mutations[(6+nb_loci):(6+nb_loci+2*nb_dmi-1)]
  dmi=unlist(dmi)
  
  sel=cbind(rep(seq(1:nb_LB),each=nb_loci/nb_LB),sel)
  sel=as.data.frame(sel)
  print(ggplot(sel,aes(x=sel,fill=as.factor(V1)))+geom_histogram()+facet_wrap(~V1)+scale_fill_manual(values=c("orange","red","magenta","blue","cyan"))+labs(fill = "Linkage block")+ggtitle("Distribution of mutational effects"))
  
  mean_sel=matrix(NA, nrow=nb_LB, ncol=4)
  for (i in 1:nb_LB){
    mean_sel[i,]=c(i,sum(sel[sel[,1]==i,2]),sum(sel[sel[,1]==i&sel[,2]>0,2]),sum(sel[sel[,1]==i&sel[,2]<0,2]))
  }
  colnames(mean_sel)=c("block","mean","min","max")
  mean_sel=as.data.frame(mean_sel)
  
  
  sel_nd=sel[!seq(0,(nb_loci-1))%in%dmi,]
  mean_sel_nd=matrix(NA, nrow=nb_LB, ncol=4)
  for (i in 1:nb_LB){
    mean_sel_nd[i,]=c(i,sum(sel_nd[sel_nd[,1]==i,2]),sum(sel_nd[sel_nd[,1]==i&sel_nd[,2]>0,2]),sum(sel_nd[sel_nd[,1]==i&sel_nd[,2]<0,2]))
  }
  colnames(mean_sel_nd)=c("block","mean","min","max")
  mean_sel_nd=as.data.frame(mean_sel_nd)
  print(ggplot(mean_sel,aes(x=block,y=mean))+geom_point(col="blue")+geom_point(aes(x=block,y=min),col="red")+geom_point(aes(x=block,y=max),col="orange") +geom_point(data=mean_sel_nd,aes(x=block+.1,y=mean),pch=15,col="blue")+geom_point(data=mean_sel_nd,aes(x=block+.1,y=min),col="red",pch=15)+geom_point(data=mean_sel_nd,aes(x=block+.1,y=max),col="orange",pch=15))
  
  sel_ndmi=sel$sel[seq(nb_loci)[-dmi]]
  sel_ndmi=sel_ndmi[sel_ndmi!=0]
  if (opt>0){
    one_option=which(sel_ndmi>0)[which(cumsum(sel_ndmi[sel_ndmi>0])<opt)]
    for(i in seq(length(sel_ndmi))[-one_option]){
      if(sel_ndmi[i]>0&(sum(sel_ndmi[one_option])+sel_ndmi[i])<opt){
        one_option=c(one_option,i)
      }
    }
  } else{
    one_option=which(sel_ndmi<0)[which(cumsum(sel_ndmi[sel_ndmi<0])>opt)]
    for(i in seq(length(sel_ndmi))[-one_option]){
      if(sel_ndmi[i]<0&(sum(sel_ndmi[one_option])+sel_ndmi[i])>opt){
        one_option=c(one_option,i)
      }
    }
  }
  mut_av=sel_ndmi[-one_option]
  
  
  
  c1=sum(sel_ndmi[one_option])+mut_av[sample(seq(length(mut_av)),10000,replace = TRUE)]
  c2=sum(sel_ndmi[one_option])+mut_av[sample(seq(length(mut_av)),10000,replace = TRUE)]+mut_av[sample(seq(length(mut_av)),10000,replace = TRUE)]                                      
  h1=hist( c1,breaks=seq(floor(min(c1,c2)*100)/100,ceiling(max(c1,c2)*100)/100,0.005),plot = FALSE)
  h2=hist( c2,breaks=seq(floor(min(c1,c2)*100)/100,ceiling(max(c1,c2)*100)/100,0.005),plot=FALSE)
  plot(h1,col="red",main=paste(c("Mean effect of 1 step:",round((mean(c1)-sum(sel_ndmi[one_option]))*1000)/1000,", of 2 steps:",round((mean(c2)-sum(sel_ndmi[one_option]))*1000)/1000,"\nProportion of mutations increasing z in 1 step:",round(mean(c1>sum(sel_ndmi[one_option]))*1000)/1000,", in 2 steps:",round(mean(c2>sum(sel_ndmi[one_option]))*1000)/1000),sep="",collapse=""))
  plot(h2,col=alpha("blue",0.5),add=TRUE)
  abline(v=sum(sel_ndmi[one_option]),col="black")
  abline(v=opt,col="gray")
  
  
  sel_coef=sel$sel
  if (opt>0){
    one_option=which(sel_coef>0)[which(cumsum(sel_coef[sel_coef>0])<opt)]
    for(i in seq(length(sel_coef))[-one_option]){
      if(sel_coef[i]>0&(sum(sel_coef[one_option])+sel_coef[i])<opt){
        one_option=c(one_option,i)
      }
    }
  }
  mut_av=sel_coef[-one_option]
  
  
  chosen1=sample(500,10000,replace=TRUE)
  chosen2=unlist(sapply(seq(10000), function(x) sample(seq(500)[-chosen1[x]],1)))
  
  c1=sum(sel_coef[one_option])-as.numeric(chosen1 %in% one_option)*sel_coef[chosen1]+as.numeric(!(chosen1 %in% one_option))*sel_coef[chosen1]
  c2=sum(sel_coef[one_option])-as.numeric(chosen1 %in% one_option)*sel_coef[chosen1]+as.numeric(!(chosen1 %in% one_option))*sel_coef[chosen1]-as.numeric(chosen2 %in% one_option)*sel_coef[chosen2]+as.numeric(!(chosen2 %in% one_option))*sel_coef[chosen2]                                     
  print(paste(c("Mean effect of 1 step:",round((mean(c1)-sum(sel_coef[one_option]))*1000)/1000,", of 2 steps:",round((mean(c2)-sum(sel_coef[one_option]))*1000)/1000,"\nProportion of mutations increasing z in 1 step:",round(mean(c1>sum(sel_coef[one_option]))*1000)/1000,", in 2 steps:",round(mean(c2>sum(sel_coef[one_option]))*1000)/1000),sep="",collapse=""))
  
  
  h1=hist( c1,breaks=seq(floor(min(c1,c2)*100)/100,ceiling(max(c1,c2)*100)/100,0.005),plot = FALSE)
  h2=hist( c2,breaks=seq(floor(min(c1,c2)*100)/100,ceiling(max(c1,c2)*100)/100,0.005),plot=FALSE)
  plot(h1,col="red",border="red",main="",xlab="Phenotype reachable through a single (or a pair of) mutation")
  plot(h2,col=alpha("blue",0.5),border=alpha("blue",0.5),add=TRUE)
  abline(v=sum(sel_ndmi[one_option]),col="black")
  abline(v=opt,col="black")
  
  #pts=seq(1,mutations[3]+nb_LB*10,1)/(mutations[3]+nb_LB*10)*(2*pi)
  nb_loci_block=mutations[3]/nb_LB
  pts=(10+seq(0,mutations[3]-1)+trunc(seq(0,mutations[3]-1)/nb_loci_block)*10+nb_LB)/(mutations[3]+nb_LB*10+20)*(2*pi)
  
  
  dmi_modified=10+dmi+trunc(dmi/nb_loci_block)*10+nb_LB
  
  dmi_apir=cbind(dmi_modified/(nb_loci+nb_LB*10+20)*2*pi, rep(1:mutations[5],each=2))
  colnames(dmi_apir)=c("pos","grp")
  dmi_apir=as.data.frame(dmi_apir)
  print(ggplot(data = data.frame(x = sin(pts), y = cos(pts)), aes(x, y))+geom_point()+geom_line(data=dmi_apir,aes(x=sin(dmi_apir$pos),y=cos(dmi_apir$pos),group=as.factor(dmi_apir$grp)),col="red"))
  
  pair_dmi=cbind(trunc(dmi/100)[c(T,F)],trunc(dmi/100)[c(F,T)])
  pair_dmi=pair_dmi[,1]+10*pair_dmi[,2]
  pair_dmi[pair_dmi==43]=34
  pair_dmi[pair_dmi==42]=24
  pair_dmi[pair_dmi==41]=14
  pair_dmi[pair_dmi==40]=4
  pair_dmi[pair_dmi==32]=23
  pair_dmi[pair_dmi==31]=13
  pair_dmi[pair_dmi==30]=3
  pair_dmi[pair_dmi==21]=12
  pair_dmi[pair_dmi==20]=2
  pair_dmi[pair_dmi==42]=24
  pair_dmi[pair_dmi==41]=14
  pair_dmi[pair_dmi==10]=1
  
  #ggplot(data = data.frame(x = sin(pts), y = cos(pts)), aes(x, y))+geom_point()+geom_line(aes(data=dmi_apir,x=sin(dmi_apir$pos),y=cos(dmi_apir$pos),group=as.factor(dmi_apir$grp)),col="red")
  #hist((1-(1-rho)^abs(dmi[c(T,F)]-dmi[c(F,T)])), "Distribution of the probability of a pair of alleles, involved in a DMI, to be herited together (no CO)")
  pair=t(apply(ceiling(cbind(dmi[c(T,F)],dmi[c(F,T)])/nb_loci_block),1,sort))
  pair_count=matrix(NA,nrow=nb_LB+1,ncol=nb_LB)
  for(i in 1:nb_LB){
    pair_count[i,1:nb_LB]=hist(pair[pair[,1]==i,2],breaks=0:nb_LB,plot=FALSE)$counts
  }
  pair_count[nb_LB+1,]=hist(dmi,breaks=seq(0,nb_loci,nb_loci_block),plot=FALSE)$counts
  rownames(pair_count)=c("LB_1","LB_2","LB_3","LB_4","LB_5","Total")
  colnames(pair_count)=c("LB_1","LB_2","LB_3","LB_4","LB_5")
  print(kable(pair_count))
  return(mean_sel)
}

display_constant=function(data_temp){
  list_col=c("nb_ind","nb_loci","nb_dim","mut_rate","mig_rate","rec_rate","rec_rate_LB","rec_rate_CL","var_mut","fit_alt","ep_FL","var_pref","var_cost","nb_dmi","ep_dmi","initial_pref","initial_cost","origin","scenario","life_cycle","opt_one_DI","opt_two_DI","dmi_scheme")
  is_constant=rep(FALSE,length(list_col))
  for (i in 1:length(list_col)){
    if(length(unique(data_temp[,list_col[i]]))==1){
      is_constant[i]=TRUE
    } else {
      print(paste(c(list_col[i],sort(unique(data_temp[,list_col[i]]))),sep=" ",collapse=" "))
    }
  }
  knitr::kable(data_temp[1,list_col[is_constant]])
}

display_constant_assay=function(data_temp){
  list_col=c("nb_ind","nb_loci","nb_dim", "mut_rate","mig_rate","rec_rate","rec_rate_LB", "rec_rate_CL", "var_mut","fit_alt","ep_FL",  "var_pref","var_cost","nb_dmi","ep_dmi", "initial_pref","initial_cost","life_cycle","origin", "opt_one_DI","opt_two_DI","scenario",   "same_conditions"  ,"barriers")
  is_constant=rep(FALSE,length(list_col))
  for (i in 1:length(list_col)){
    if(length(unique(data_temp[,list_col[i]]))==1){
      is_constant[i]=TRUE
    } else {
      print(paste(c(list_col[i],unique(data_temp[,list_col[i]])),collapse=""))
    }
  }
  knitr::kable(data_temp[1,list_col[is_constant]])
}

simulate_dyn_Ab=function(nb_pairs=50,pop_size=100,nb_gen=50,isSecContact=FALSE){
  if(isSecContact){
    ind_locus=rep(c(-1,1),nb_pairs/2)
  } else{
    ind_locus=rep(0,nb_pairs)
  }
  for(i in 1:(nb_gen*pop_size)){
    mut=runif(nb_pairs)
    ind=seq(nb_pairs)[mut<1/pop_size]
    if (length(ind)!=0){
      ind_locus[ind]=ind_locus[ind]+sample(c(-1,1),length(ind),replace = TRUE)
    }
    ind_locus[ind_locus< -1]=-1
    ind_locus[ind_locus>1]=1
  }
  #prob_AB=matrix(0,nrow=50000,ncol=50)
  #for(i in 1:50){
  #  prob_AB[,i]=freq_t[,4*(i-1)+1]*freq_t[,4*(i-1)+4]+freq_t[,4*(i-1)+2]*freq_t[,4*(i-1)+3]
  ind_locus
}
