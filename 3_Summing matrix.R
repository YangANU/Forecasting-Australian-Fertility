####################################################################
# Summing matrix 1: fixed COB, disaggregate by geographical factors
####################################################################

Smat_fun_fix_cob <- function(kj, age, no_area = 47, cob = 1)
{
  age_ind = match(age, age_all)
  
  ###################
  # Level_1 national
  ###################
  
  level_1 =  matrix(0,no_area,(11-kj))
  for(iw in 1:no_area)
  {
    level_1[iw,] = get(pop_ratio_Area_to_AUS[iw])[kj,age_ind,1:(11-kj), cob]
  }
  
  #################################  
  # Level 2 disaggregate by region
  #################################
  
  level_2_r_1 = level_2_r_2 = level_2_r_3 = level_2_r_4 = level_2_r_5 = 
    level_2_r_6 = level_2_r_7 = level_2_r_8 = level_2_r_9 = level_2_r_10 = 
    level_2_r_11 = matrix(0,no_area,(11-kj))
  
  # R1
  level_2_r_1[region_list_ind[[1]],] = get(pop_ratio_Area_to_R1[1])[kj,age_ind,1:(11-kj), cob]
  
  # R2
  for(ij in 1:length(pop_ratio_Area_to_R2))
  {
    row_ind = region_list_ind[[2]][ij]
    level_2_r_2[row_ind,] = get(pop_ratio_Area_to_R2[ij])[kj,age_ind,1:(11-kj), cob]
  }
  
  # R3
  level_2_r_3[region_list_ind[[3]],] = get(pop_ratio_Area_to_R3[1])[kj,age_ind,1:(11-kj), cob]
  
  # R4
  for(ij in 1:length(pop_ratio_Area_to_R4))
  {
    row_ind = region_list_ind[[4]][ij]
    level_2_r_4[row_ind,] = get(pop_ratio_Area_to_R4[ij])[kj,age_ind,1:(11-kj), cob]
  }
  
  # R5
  level_2_r_5[region_list_ind[[5]],] = get(pop_ratio_Area_to_R5[1])[kj,age_ind,1:(11-kj), cob]
  
  # R6
  level_2_r_6[region_list_ind[[6]],] = get(pop_ratio_Area_to_R6[1])[kj,age_ind,1:(11-kj), cob]
  
  # R7
  level_2_r_7[region_list_ind[[7]],] = get(pop_ratio_Area_to_R7[1])[kj,age_ind,1:(11-kj), cob]
  
  # R8
  level_2_r_8[region_list_ind[[8]],] = get(pop_ratio_Area_to_R8[1])[kj,age_ind,1:(11-kj), cob]
  
  # R9
  level_2_r_9[region_list_ind[[9]],] = get(pop_ratio_Area_to_R9[1])[kj,age_ind,1:(11-kj), cob]
  
  # R10
  for(ij in 1:length(pop_ratio_Area_to_R10))
  {
    row_ind = region_list_ind[[10]][ij]
    level_2_r_10[row_ind,] = get(pop_ratio_Area_to_R10[ij])[kj,age_ind,1:(11-kj), cob]
  }
  
  # R11
  for(ij in 1:length(pop_ratio_Area_to_R11))
  {
    row_ind = region_list_ind[[11]][ij]
    level_2_r_11[row_ind,] = get(pop_ratio_Area_to_R11[ij])[kj,age_ind,1:(11-kj), cob]
  }
  
  ###############################  
  # Level 3 disaggregate by area
  ###############################
  
  level_3 = diag(no_area)
  
  S_mat = array(0, dim = c(59, no_area, (11-kj)))
  for(ik in 1:(11-kj))
  {
    S_mat[,,ik] = rbind(level_1[,ik], level_2_r_1[,ik], level_2_r_2[,ik], level_2_r_3[,ik], 
                        level_2_r_4[,ik], level_2_r_5[,ik], level_2_r_6[,ik], level_2_r_7[,ik],
                        level_2_r_8[,ik], level_2_r_9[,ik], level_2_r_10[,ik], level_2_r_11[,ik],
                        level_3)
  }
  
  return(S_mat)
}

test = Smat_fun(kj = 1, age = 15, no_area = 47, cob = 2)




############################################
# W_h matrix for MinT reconciliation method
############################################

wh_fun_fix_cob <- function(kj, age, cob_ind)
{
  age_ind = match(age, age_all)
  
  lowerD = hts:::lowerD
  shrink.estim = hts:::shrink.estim
  
  eh_mat = matrix(NA, nrow = 59, ncol = (20+kj))
  
  # level 1
  eh_mat[1,] = get(national_train_residual_fix_COB[1])[[cob_ind]][[kj]][age_ind,]
  
  # level 2 (disaggregate by region)
  for (ikw in 1:11)
  {
    eh_mat[(1+ikw),] = get(region_train_residual_fix_COB[ikw])[[cob_ind]][[kj]][age_ind,]
  }
  
  # level 4 (disaggregate by area)
  for (ikw in 1:47)
  {
    if(ikw == 23)
    {
      eh_mat[(12+ikw),] = rep(1e-6, (20+kj))
    } else {
      eh_mat[(12+ikw),] = get(area_train_residual_fix_COB[(ikw)])[[cob_ind]][[kj]][age_ind,]
    }
  }

  target = lowerD(t(eh_mat))
  shrink = shrink.estim(t(eh_mat), target)
  wh_mat = shrink[[1]]
  
  return(wh_mat = wh_mat)
}


ind_wh_fun_fix_cob <- function(kj, age, cob_ind)
{
  age_ind = match(age, age_all)
  
  lowerD = hts:::lowerD
  shrink.estim = hts:::shrink.estim
  
  eh_mat = matrix(NA, nrow = 59, ncol = (20+kj))
  
  # level 1
  eh_mat[1,] = get(ind_national_train_residual_fix_COB[1])[[cob_ind]][[kj]][age_ind,]
  
  # level 2 (disaggregate by region)
  for (ikw in 1:11)
  {
    eh_mat[(1+ikw),] = get(ind_region_train_residual_fix_COB[ikw])[[cob_ind]][[kj]][age_ind,]
  }
  
  # level 4 (disaggregate by area)
  for (ikw in 1:47)
  {
    if(ikw == 23)
    {
      eh_mat[(12+ikw),] = rep(1e-6, (20+kj))
    } else {
      eh_mat[(12+ikw),] = get(ind_area_train_residual_fix_COB[(ikw)])[[cob_ind]][[kj]][age_ind,]
    }
  }
  
  target = lowerD(t(eh_mat))
  shrink = shrink.estim(t(eh_mat), target)
  wh_mat = shrink[[1]]
  
  return(wh_mat = wh_mat)
}


###########################################################
# Summing matrix 2: fixed GEO, disaggregate by birthplaces
###########################################################

Smat_fun_fix_geo <- function(kj, age, no_cob = 19, geo = 1)
{
  age_ind = match(age, age_all)
  
  #################
  # Level_1 global
  #################
  
  level_1 =  matrix(0,no_cob,(11-kj))
  for(iw in 1:no_cob)
  {
    level_1[iw,] = get(pop_ratio_cob_to_national[iw])[kj,age_ind,1:(11-kj), geo]
  }
  
  ##################################  
  # Level 2 disaggregate by C_group
  ##################################
  
  level_2_r_1 = level_2_r_2 = level_2_r_3 = level_2_r_4 = level_2_r_5 = 
    level_2_r_6 = level_2_r_7 = level_2_r_8 = level_2_r_9 = level_2_r_10 = matrix(0,no_cob,(11-kj))
  
  # C_group 1
  for(ij in 1:length(pop_ratio_cob_to_C1))
  {
    row_ind = COB_list[[1]][ij]
    level_2_r_2[row_ind,] = get(pop_ratio_cob_to_C1[ij])[kj,age_ind,1:(11-kj), geo]
  }

  # C_group 2
  for(ij in 1:length(pop_ratio_cob_to_C2))
  {
    row_ind = COB_list[[2]][ij]
    level_2_r_2[row_ind,] = get(pop_ratio_cob_to_C2[ij])[kj,age_ind,1:(11-kj), geo]
  }
  
  # C_group 3
  level_2_r_3[COB_list[[3]],] = get(pop_ratio_cob_to_C3[1])[kj,age_ind,1:(11-kj), geo]
  
  # C_group 4
  level_2_r_4[COB_list[[4]],] = get(pop_ratio_cob_to_C4[1])[kj,age_ind,1:(11-kj), geo]
  
  # C_group 5
  for(ij in 1:length(pop_ratio_cob_to_C5))
  {
    row_ind = COB_list[[5]][ij]
    level_2_r_5[row_ind,] = get(pop_ratio_cob_to_C5[ij])[kj,age_ind,1:(11-kj), geo]
  }

  # C_group 6
  for(ij in 1:length(pop_ratio_cob_to_C6))
  {
    row_ind = COB_list[[6]][ij]
    level_2_r_6[row_ind,] = get(pop_ratio_cob_to_C6[ij])[kj,age_ind,1:(11-kj), geo]
  }
  
  # C_group 7
  for(ij in 1:length(pop_ratio_cob_to_C7))
  {
    row_ind = COB_list[[7]][ij]
    level_2_r_7[row_ind,] = get(pop_ratio_cob_to_C7[ij])[kj,age_ind,1:(11-kj), geo]
  }
  
  # C_group 8
  level_2_r_8[COB_list[[8]],] = get(pop_ratio_cob_to_C8[1])[kj,age_ind,1:(11-kj), cob]
  
  # C_group 9
  level_2_r_9[COB_list[[9]],] = get(pop_ratio_cob_to_C9[1])[kj,age_ind,1:(11-kj), cob]
  
  # C_group 10
  level_2_r_10[COB_list[[10]],] = get(pop_ratio_cob_to_C10[1])[kj,age_ind,1:(11-kj), cob]
  
  ###############################  
  # Level 3 disaggregate by area
  ###############################
  
  level_3 = diag(no_cob)
  
  S_mat = array(0, dim = c(30, no_cob, (11-kj)))
  for(ik in 1:(11-kj))
  {
    S_mat[,,ik] = rbind(level_1[,ik], level_2_r_1[,ik], level_2_r_2[,ik], level_2_r_3[,ik], 
                        level_2_r_4[,ik], level_2_r_5[,ik], level_2_r_6[,ik], level_2_r_7[,ik],
                        level_2_r_8[,ik], level_2_r_9[,ik], level_2_r_10[,ik], level_3)
  }
  
  return(S_mat)
}

test = Smat_fun_fix_geo(kj = 1, age = 15, no_cob = 19, geo = 2)


############################################
# W_h matrix for MinT reconciliation method
############################################

wh_fun_fix_geo <- function(kj, age, geo_ind)
{
  age_ind = match(age, age_all)
  
  lowerD = hts:::lowerD
  shrink.estim = hts:::shrink.estim
  
  eh_mat = matrix(NA, nrow = 30, ncol = (20+kj))
  
  # level 1
  eh_mat[1,] = get(global_train_residual_fix_GEO[1])[[geo_ind]][[kj]][age_ind,]
  
  # level 2 (disaggregate by region)
  for (ikw in 1:10)
  {
    if(ikw == 6 && geo_ind == 23)
    {
      eh_mat[(1+ikw),] = rep(1e-6, (20+kj))
    } else {
      eh_mat[(1+ikw),] = get(C_group_train_residual_fix_GEO[ikw])[[geo_ind]][[kj]][age_ind,]
    }
  }
  
  # level 4 (disaggregate by area)
  for (ikw in 1:19)
  {
    if(ikw == 14 && geo_ind == 23)
    {
      eh_mat[(11+ikw),] = rep(1e-6, (20+kj))
    } else {
      eh_mat[(11+ikw),] = get(cob_train_residual_fix_GEO[(ikw)])[[geo_ind]][[kj]][age_ind,]
    }
  }
  
  target = lowerD(t(eh_mat))
  shrink = shrink.estim(t(eh_mat), target)
  wh_mat = shrink[[1]]
  
  return(wh_mat = wh_mat)
}

test = wh_fun_fix_geo(kj = 1, age = 25, geo_ind = 23)

ind_wh_fun_fix_geo <- function(kj, age, geo_ind)
{
  age_ind = match(age, age_all)
  
  lowerD = hts:::lowerD
  shrink.estim = hts:::shrink.estim
  
  eh_mat = matrix(NA, nrow = 30, ncol = (20+kj))
  
  # level 1
  eh_mat[1,] = get(ind_global_train_residual_fix_GEO[1])[[geo_ind]][[kj]][age_ind,]
  
  # level 2 (disaggregate by region)
  for (ikw in 1:10)
  {
    if(ikw == 6 && geo_ind == 23)
    {
      eh_mat[(1+ikw),] = rep(1e-6, (20+kj))
    } else {
      eh_mat[(1+ikw),] = get(ind_c_group_train_residual_fix_GEO[ikw])[[geo_ind]][[kj]][age_ind,]
    }
  }
  
  # level 4 (disaggregate by area)
  for (ikw in 1:19)
  {
    if(ikw == 14 && geo_ind == 23)
    {
      eh_mat[(11+ikw),] = rep(1e-6, (20+kj))
    } else {
      eh_mat[(11+ikw),] = get(ind_cob_train_residual_fix_GEO[(ikw)])[[geo_ind]][[kj]][age_ind,]
    }
  }
  
  target = lowerD(t(eh_mat))
  shrink = shrink.estim(t(eh_mat), target)
  wh_mat = shrink[[1]]
  
  return(wh_mat = wh_mat)
}

test_2 = ind_wh_fun_fix_geo(kj = 1, age = 25, geo_ind = 23)









