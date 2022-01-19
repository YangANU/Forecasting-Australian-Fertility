##########################################################
# Construct summing matrix S by forecast exposure to risk
##########################################################

Smat_fun <- function(kj, age, no_area = 47)
{
  level_0 = level_1_a = level_1_b = matrix(NA, 2*no_area, (6-kj))
  for(iw in 1:no_area)
  {
    # Level 0 
    
    level_0[(2*iw-1),] = get(pop_ratio_A_F_to_T[iw])[age, kj:5]
    level_0[(2*iw),]   = get(pop_ratio_A_M_to_T[iw])[age, kj:5]
    
    # Level 1 (disaggregate by sex)
    
    level_1_a[(2*iw-1),] = get(pop_ratio_A_F_to_F[iw])[age, kj:5]
    level_1_a[(2*iw),]   = rep(0,(6-kj))
    
    level_1_b[(2*iw-1),] = rep(0,(6-kj))
    level_1_b[(2*iw),]   = get(pop_ratio_A_M_to_M[iw])[age, kj:5]
  }
  

  ###################################  
  # Level 2 (disaggregate by region)
  ###################################
  
  level_2_a_1 = level_2_a_2 = level_2_a_3 = level_2_a_4 = level_2_a_5 = level_2_a_6 = 
    level_2_a_7 = level_2_a_8 = level_2_a_9 = level_2_a_10 = level_2_a_11 = matrix(NA, 2*no_area, (6-kj))
  
  # R01
  
  level_2_a_1[1,] = get(pop_ratio_A_F_to_R1_T[1])[age, kj:5]  
  level_2_a_1[2,] = get(pop_ratio_A_M_to_R1_T[1])[age, kj:5]
  level_2_a_1[3:(2*no_area),] = matrix(rep(0,(6-kj)*92), ncol=(6-kj))
  
  # R02
  
  level_2_a_2[1:2, ] = matrix(rep(0,(6-kj)*2), ncol=(6-kj))
  for(iw in 1:3)
  {
    level_2_a_2[(2*iw+1),] = get(pop_ratio_A_F_to_R2_T[iw])[age, kj:5]
    level_2_a_2[(2*iw+2),] = get(pop_ratio_A_M_to_R2_T[iw])[age, kj:5]
  }  
  level_2_a_2[9:(2*no_area),] = matrix(rep(0,(6-kj)*86), ncol=(6-kj))
  
  # R03
  
  level_2_a_3[1:8,] = matrix(rep(0,(6-kj)*8), ncol=(6-kj))  
  
  level_2_a_3[9,] = get(pop_ratio_A_F_to_R3_T[1])[age, kj:5]
  level_2_a_3[10,] = get(pop_ratio_A_M_to_R3_T[1])[age, kj:5]
  
  level_2_a_3[11:(2*no_area),] = matrix(rep(0,(6-kj)*84), ncol=(6-kj))
  
  # R04
  
  level_2_a_4[1:10,] = matrix(rep(0,(6-kj)*10), ncol=(6-kj))
  for(iw in 1:5)
  {
    level_2_a_4[(2*(iw+4)+1),] = get(pop_ratio_A_F_to_R4_T[iw])[age, kj:5]
    level_2_a_4[(2*(iw+4)+2),] = get(pop_ratio_A_M_to_R4_T[iw])[age, kj:5]
  }
  level_2_a_4[21:(2*no_area),] = matrix(rep(0,(6-kj)*74), ncol=(6-kj))
  
  # R05
  
  level_2_a_5[1:20,] = matrix(rep(0,(6-kj)*20), ncol=(6-kj))

  level_2_a_5[21,] = get(pop_ratio_A_F_to_R5_T[1])[age, kj:5]
  level_2_a_5[22,] = get(pop_ratio_A_M_to_R5_T[1])[age, kj:5]
  
  level_2_a_5[23:(2*no_area),] = matrix(rep(0,(6-kj)*72), ncol=(6-kj))
  
  # R06
  
  level_2_a_6[1:22,] = matrix(rep(0,(6-kj)*22), ncol=(6-kj))
  
  level_2_a_6[23,] = get(pop_ratio_A_F_to_R6_T[1])[age, kj:5]
  level_2_a_6[24,] = get(pop_ratio_A_M_to_R6_T[1])[age, kj:5]
  
  level_2_a_6[25:(2*no_area),] = matrix(rep(0,(6-kj)*70), ncol=(6-kj))
  
  # R07
  
  level_2_a_7[1:24,] = matrix(rep(0,(6-kj)*24), ncol=(6-kj))
  
  level_2_a_7[25,] = get(pop_ratio_A_F_to_R7_T[1])[age, kj:5]
  level_2_a_7[26,] = get(pop_ratio_A_M_to_R7_T[1])[age, kj:5]
  
  level_2_a_7[27:(2*no_area),] = matrix(rep(0,(6-kj)*68), ncol=(6-kj))
  
  # R08
  
  level_2_a_8[1:26,] = matrix(rep(0,(6-kj)*26), ncol=(6-kj))
  
  level_2_a_8[27,] = get(pop_ratio_A_F_to_R8_T[1])[age, kj:5]
  level_2_a_8[28,] = get(pop_ratio_A_M_to_R8_T[1])[age, kj:5]
  
  level_2_a_8[29:(2*no_area),] = matrix(rep(0,(6-kj)*66), ncol=(6-kj))
  
  # R09
  
  level_2_a_9[1:28,] = matrix(rep(0,(6-kj)*28), ncol=(6-kj))
  
  level_2_a_9[29,] = get(pop_ratio_A_F_to_R9_T[1])[age, kj:5]
  level_2_a_9[30,] = get(pop_ratio_A_M_to_R9_T[1])[age, kj:5]
  
  level_2_a_9[31:(2*no_area),] = matrix(rep(0,(6-kj)*64), ncol=(6-kj))
  
  # R10
  
  level_2_a_10[1:30,] = matrix(rep(0,(6-kj)*30), ncol=(6-kj))
  for(iw in 1:19)
  {
    level_2_a_10[(2*(iw+14)+1),] = get(pop_ratio_A_F_to_R10_T[iw])[age, kj:5]
    level_2_a_10[(2*(iw+14)+2),] = get(pop_ratio_A_M_to_R10_T[iw])[age, kj:5]
  }
  level_2_a_10[69:(2*no_area),] = matrix(rep(0,(6-kj)*26), ncol=(6-kj))
  
  # R11
  
  level_2_a_11[1:68,] = matrix(rep(0,(6-kj)*68), ncol=(6-kj))
  for(iw in 1:13)
  {
    level_2_a_11[(2*(iw+33)+1),] = get(pop_ratio_A_F_to_R11_T[iw])[age, kj:5]
    level_2_a_11[(2*(iw+33)+2),] = get(pop_ratio_A_M_to_R11_T[iw])[age, kj:5]
  }
  
  
  ##################################################
  # Level 2 (disaggregate by region + sex (female))
  ##################################################
  
  level_2_b_1 = level_2_b_2 = level_2_b_3 = level_2_b_4 = level_2_b_5 = level_2_b_6 = 
    level_2_b_7 = level_2_b_8 = level_2_b_9 = level_2_b_10 = level_2_b_11 = matrix(NA, (2*no_area),(6-kj))
  
  # R01
  
  level_2_b_1[1,] = get(pop_ratio_A_F_to_R1_F[1])[age, kj:5]  
  level_2_b_1[2:(2*no_area),] = matrix(rep(0,(6-kj)*93), ncol=(6-kj))
  
  # R02
  
  level_2_b_2[1:2, ] = matrix(rep(0,(6-kj)*2), ncol=(6-kj))
  for(iw in 1:3)
  {
    level_2_b_2[(2*iw+1),] = get(pop_ratio_A_F_to_R2_F[iw])[age, kj:5]
    level_2_b_2[(2*iw+2),] = rep(0, (6-kj))
  }  
  level_2_b_2[9:(2*no_area),] = matrix(rep(0,(6-kj)*86), ncol=(6-kj))
  
  # R03
  
  level_2_b_3[1:8,] = matrix(rep(0,(6-kj)*8), ncol=(6-kj))  

  level_2_b_3[9,] = get(pop_ratio_A_F_to_R3_F[1])[age, kj:5]
  level_2_b_3[10,] = rep(0, (6-kj))
  
  level_2_b_3[11:(2*no_area),] = matrix(rep(0,(6-kj)*84), ncol=(6-kj))
  
  # R04
  
  level_2_b_4[1:10,] = matrix(rep(0,(6-kj)*10), ncol=(6-kj))
  for(iw in 1:5)
  {
    level_2_b_4[(2*(iw+4)+1),] = get(pop_ratio_A_F_to_R4_F[iw])[age, kj:5]
    level_2_b_4[(2*(iw+4)+2),] = rep(0, (6-kj))
  }
  level_2_b_4[21:(2*no_area),] = matrix(rep(0,(6-kj)*74), ncol=(6-kj))
  
  # R05
  
  level_2_b_5[1:20,] = matrix(rep(0,(6-kj)*20), ncol=(6-kj))
  
  level_2_b_5[21,] = get(pop_ratio_A_F_to_R5_F[1])[age, kj:5]
  level_2_b_5[22,] = rep(0, (6-kj))
  
  level_2_b_5[23:(2*no_area),] = matrix(rep(0,(6-kj)*72), ncol=(6-kj))
  
  # R06
  
  level_2_b_6[1:22,] = matrix(rep(0,(6-kj)*22), ncol=(6-kj))
  
  level_2_b_6[23,] = get(pop_ratio_A_F_to_R6_F[1])[age, kj:5]
  level_2_b_6[24,] = rep(0, (6-kj))
  
  level_2_b_6[25:(2*no_area),] = matrix(rep(0,(6-kj)*70), ncol=(6-kj))
  
  # R07
  
  level_2_b_7[1:24,] = matrix(rep(0,(6-kj)*24), ncol=(6-kj))
  
  level_2_b_7[25,] = get(pop_ratio_A_F_to_R7_F[1])[age, kj:5]
  level_2_b_7[26,] = rep(0, (6-kj))
  
  level_2_b_7[27:(2*no_area),] = matrix(rep(0,(6-kj)*68), ncol=(6-kj))
  
  # R08
  
  level_2_b_8[1:26,] = matrix(rep(0,(6-kj)*26), ncol=(6-kj))
  
  level_2_b_8[27,] = get(pop_ratio_A_F_to_R8_F[1])[age, kj:5]
  level_2_b_8[28,] = rep(0, (6-kj))
  
  level_2_b_8[29:(2*no_area),] = matrix(rep(0,(6-kj)*66), ncol=(6-kj))
  
  # R09
  
  level_2_b_9[1:28,] = matrix(rep(0,(6-kj)*28), ncol=(6-kj))
  
  level_2_b_9[29,] = get(pop_ratio_A_F_to_R9_F[1])[age, kj:5]
  level_2_b_9[30,] = rep(0, (6-kj))
  
  level_2_b_9[31:(2*no_area),] = matrix(rep(0,(6-kj)*64), ncol=(6-kj))
  
  # R10
  
  level_2_b_10[1:30,] = matrix(rep(0,(6-kj)*30), ncol=(6-kj))
  for(iw in 1:19)
  {
    level_2_b_10[(2*(iw+14)+1),] = get(pop_ratio_A_F_to_R10_F[iw])[age, kj:5]
    level_2_b_10[(2*(iw+14)+2),] = rep(0, (6-kj))
  }
  level_2_b_10[69:(2*no_area),] = matrix(rep(0,(6-kj)*26), ncol=(6-kj))
  
  # R11
  
  level_2_b_11[1:68,] = matrix(rep(0,(6-kj)*68), ncol=(6-kj))
  for(iw in 1:13)
  {
    level_2_b_11[(2*(iw+33)+1),] = get(pop_ratio_A_F_to_R11_F[iw])[age, kj:5]
    level_2_b_11[(2*(iw+33)+2),] = rep(0, (6-kj))
  }
  
  
  ################################################
  # Level 2 (disaggregate by region + sex (male))
  ################################################
  
  level_2_c_1 = level_2_c_2 = level_2_c_3 = level_2_c_4 = level_2_c_5 = level_2_c_6 = 
    level_2_c_7 = level_2_c_8 = level_2_c_9 = level_2_c_10 = level_2_c_11 = matrix(NA, (2*no_area), (6-kj))
  
  # R01
  
  level_2_c_1[1,] = rep(0,(6-kj))
  level_2_c_1[2,] = get(pop_ratio_A_M_to_R1_M[1])[age, kj:5]
  level_2_c_1[3:(2*no_area),] = matrix(rep(0,(6-kj)*92), ncol=(6-kj))
  
  # R02
  
  level_2_c_2[1:2, ] = matrix(rep(0,(6-kj)*2), ncol=(6-kj))
  for(iw in 1:3)
  {
    level_2_c_2[(2*iw+1),] = rep(0,(6-kj))
    level_2_c_2[(2*iw+2),] = get(pop_ratio_A_M_to_R2_M[iw])[age, kj:5]
  }  
  level_2_c_2[9:(2*no_area),] = matrix(rep(0,(6-kj)*86), ncol=(6-kj))
  
  # R03
  
  level_2_c_3[1:8,] = matrix(rep(0,(6-kj)*8), ncol=(6-kj))  
  for(iw in 1:7)
  
  level_2_c_3[9,] = rep(0,(6-kj))
  level_2_c_3[10,] = get(pop_ratio_A_M_to_R3_M[1])[age, kj:5]
  
  level_2_c_3[11:(2*no_area),] = matrix(rep(0,(6-kj)*84), ncol=(6-kj))
  
  # R04
  
  level_2_c_4[1:10,] = matrix(rep(0,(6-kj)*10), ncol=(6-kj))
  for(iw in 1:5)
  {
    level_2_c_4[(2*(iw+4)+1),] = rep(0,(6-kj))
    level_2_c_4[(2*(iw+4)+2),] = get(pop_ratio_A_M_to_R4_M[iw])[age, kj:5]
  }
  level_2_c_4[21:(2*no_area),] = matrix(rep(0,(6-kj)*74), ncol=(6-kj))
  
  # R05
  
  level_2_c_5[1:20,] = matrix(rep(0,(6-kj)*20), ncol=(6-kj))
  
  level_2_c_5[21,] = rep(0,(6-kj))
  level_2_c_5[22,] = get(pop_ratio_A_M_to_R5_M[1])[age, kj:5]
  
  level_2_c_5[23:(2*no_area),] = matrix(rep(0,(6-kj)*72), ncol=(6-kj))
  
  # R06
  
  level_2_c_6[1:22,] = matrix(rep(0,(6-kj)*22), ncol=(6-kj))
  
  level_2_c_6[23,] = rep(0,(6-kj))
  level_2_c_6[24,] = get(pop_ratio_A_M_to_R6_M[1])[age, kj:5]
  
  level_2_c_6[25:(2*no_area),] = matrix(rep(0,(6-kj)*70), ncol=(6-kj))
  
  # R07
  
  level_2_c_7[1:24,] = matrix(rep(0,(6-kj)*24), ncol=(6-kj))

  level_2_c_7[25,] = rep(0,(6-kj))
  level_2_c_7[26,] = get(pop_ratio_A_M_to_R7_M[1])[age, kj:5]
  
  level_2_c_7[27:(2*no_area),] = matrix(rep(0,(6-kj)*68), ncol=(6-kj))
  
  # R08
  
  level_2_c_8[1:26,] = matrix(rep(0,(6-kj)*26), ncol=(6-kj))
  
  level_2_c_8[27,] = rep(0,(6-kj))
  level_2_c_8[28,] = get(pop_ratio_A_M_to_R8_M[1])[age, kj:5] 
  
  level_2_c_8[29:(2*no_area),] = matrix(rep(0,(6-kj)*66), ncol=(6-kj))
  
  # R09
  
  level_2_c_9[1:28,] = matrix(rep(0,(6-kj)*28), ncol=(6-kj))
  
  level_2_c_9[29,] = rep(0, (6-kj))
  level_2_c_9[30,] = get(pop_ratio_A_M_to_R9_M[1])[age, kj:5]
  
  level_2_c_9[31:(2*no_area),] = matrix(rep(0,(6-kj)*64), ncol=(6-kj))
  
  # R10
  
  level_2_c_10[1:30,] = matrix(rep(0,(6-kj)*30), ncol=(6-kj))
  for(iw in 1:19)
  {
    level_2_c_10[(2*(iw+14)+1),] = rep(0, (6-kj))
    level_2_c_10[(2*(iw+14)+2),] = get(pop_ratio_A_M_to_R10_M[iw])[age, kj:5]
  }
  level_2_c_10[69:(2*no_area),] = matrix(rep(0,(6-kj)*26), ncol=(6-kj))
  
  # R11
  
  level_2_c_11[1:68,] = matrix(rep(0,(6-kj)*68), ncol=(6-kj))
  for(iw in 1:13)
  {
    level_2_c_11[(2*(iw+33)+1),] = rep(0, (6-kj))
    level_2_c_11[(2*(iw+33)+2),] = get(pop_ratio_A_M_to_R11_M[iw])[age, kj:5]
  }
  
  
  #############################################################
  # Level 3 (disaggregate by prefecture + sex (female & male))
  #############################################################
  
  level_3 = array(0, dim = c((2*no_area), (6-kj), 47))
  for(iw in 1:47)
  {
    level_3[2*iw-1,,iw] = get(pop_ratio_A_F_to_A_T[iw])[age, kj:5]
    level_3[2*iw,,iw]   = get(pop_ratio_A_M_to_A_T[iw])[age, kj:5]
  }
  
  S_mat = array(NA, dim = c(177, (2*no_area), (6-kj)))
  for(ik in 1:(6-kj))
  {
    S_mat[,,ik] = rbind(level_0[,ik],     level_1_a[,ik],   level_1_b[,ik],   level_2_a_1[,ik],
                        level_2_a_2[,ik], level_2_a_3[,ik], level_2_a_4[,ik], level_2_a_5[,ik],
                        level_2_a_6[,ik], level_2_a_7[,ik], level_2_a_8[,ik], level_2_a_9[,ik],
                        level_2_a_10[,ik], level_2_a_11[,ik], 
                        level_2_b_1[,ik], level_2_b_2[,ik], level_2_b_3[,ik], level_2_b_4[,ik], level_2_b_5[,ik],
                        level_2_b_6[,ik], level_2_b_7[,ik], level_2_b_8[,ik], level_2_b_9[,ik], level_2_b_10[,ik],
                        level_2_b_11[,ik], 
                        level_2_c_1[,ik], level_2_c_2[,ik], level_2_c_3[,ik], level_2_c_4[,ik], level_2_c_5[,ik],
                        level_2_c_6[,ik], level_2_c_7[,ik], level_2_c_8[,ik], level_2_c_9[,ik], level_2_c_10[,ik],
                        level_2_c_11[,ik],
                        t(level_3[,ik,]),
                        diag(94))
  }
  return(S_mat)
}



############################################
# W_h matrix for MinT reconciliation method
############################################

# Univariate MinT reconciliation function

wh_fun <- function(kj, age)
{
  lowerD = hts:::lowerD
  shrink.estim = hts:::shrink.estim
  
  eh_mat = matrix(NA, nrow = 177, ncol = (18+kj))
  
  # level 0
  eh_mat[1,] = get(ind_Area_train_residual_total[1])[[kj]][age,]
  
  # level 1 (disaggregate by sex)
  eh_mat[2,] = get(ind_Area_train_residual_female[1])[[kj]][age,]
  eh_mat[3,] = get(ind_Area_train_residual_male[1])[[kj]][age,]
  
  # level 2 & 3 (disaggregate by Region & sex)
  for (ikw in 1:11)
  {
    eh_mat[(3+ikw),] = get(ind_Region_train_residual_total[ikw])[[kj]][age,]
    eh_mat[(14+ikw),] = get(ind_Region_train_residual_female[ikw])[[kj]][age,]
    eh_mat[(25+ikw),] = get(ind_Region_train_residual_male[ikw])[[kj]][age,]
  }
  
  # level 4 (disaggregate by state)
  for (ikw in 1:47)
  {
    eh_mat[(36+ikw),] = get(ind_Area_train_residual_total[(ikw+1)])[[kj]][age,]
  }
  
  # level 5 (disaggregate by state & sex)
  for (ikw in 2:48)
  {
    eh_mat[80+(2*ikw),] = get(ind_Area_train_residual_female[(ikw)])[[kj]][age,]
    eh_mat[81+(2*ikw),] = get(ind_Area_train_residual_male[(ikw)])[[kj]][age,]
  }
  
  target = lowerD(t(eh_mat))
  shrink = shrink.estim(t(eh_mat), target)
  wh_mat = shrink[[1]]
  
  return(wh_mat = wh_mat)
}

# Multivariate MinT reconciliation function

wh_fun_mfts <- function(kj, age)
{
  lowerD = hts:::lowerD
  shrink.estim = hts:::shrink.estim
  
  eh_mat = matrix(NA, nrow = 177, ncol = (18+kj))
  
  # level 0
  eh_mat[1,] = get(mfts_Area_train_residual_total[1])[[kj]][age,]
  
  # level 1 (disaggregate by sex)
  eh_mat[2,] = get(mfts_Area_train_residual_female[1])[[kj]][age,]
  eh_mat[3,] = get(mfts_Area_train_residual_male[1])[[kj]][age,]
  
  # level 2 & 3 (disaggregate by Region & sex)
  for (ikw in 1:11)
  {
    eh_mat[(3+ikw),] = get(mfts_Region_train_residual_total[ikw])[[kj]][age,]
    eh_mat[(14+ikw),] = get(mfts_Region_train_residual_female[ikw])[[kj]][age,]
    eh_mat[(25+ikw),] = get(mfts_Region_train_residual_male[ikw])[[kj]][age,]
  }
  
  # level 4 (disaggregate by Area)
  for (ikw in 1:47)
  {
    eh_mat[(36+ikw),] = get(mfts_Area_train_residual_total[(ikw+1)])[[kj]][age,]
  }
  
  # level 5 (disaggregate by Area & sex)
  for (ikw in 2:48)
  {
    eh_mat[80+(2*ikw),] = get(mfts_Area_train_residual_female[(ikw)])[[kj]][age,]
    eh_mat[81+(2*ikw),] = get(mfts_Area_train_residual_male[(ikw)])[[kj]][age,]
  }
  
  target = lowerD(t(eh_mat))
  shrink = shrink.estim(t(eh_mat), target)
  wh_mat = shrink[[1]]
  
  return(wh_mat = wh_mat)
}


wh_fun_gender <- function(kj, age)
{
  lowerD = hts:::lowerD
  shrink.estim = hts:::shrink.estim
  
  eh_mat = matrix(NA, nrow = 177, ncol = (18+kj))
  
  # level 0
  eh_mat[1,] = get(mfts_Area_train_residual_total[1])[[kj]][age,]
  
  # level 1 (disaggregate by sex)
  eh_mat[2,] = get(mfts_Area_train_residual_female[1])[[kj]][age,]
  eh_mat[3,] = get(mfts_Area_train_residual_male[1])[[kj]][age,]
  
  # level 2 & 3 (disaggregate by Region & sex)
  for (ikw in 1:8)
  {
    eh_mat[(3+ikw),] = get(mfts_Region_train_residual_total[ikw])[[kj]][age,]
    eh_mat[(14+ikw),] = get(mfts_gender_Region_train_residual_female[ikw])[[kj]][age,]
    eh_mat[(25+ikw),] = get(mfts_gender_Region_train_residual_male[ikw])[[kj]][age,]
  }
  
  # level 4 (disaggregate by Area)
  for (ikw in 1:47)
  {
    eh_mat[(36+ikw),] = get(mfts_Area_train_residual_total[(ikw+1)])[[kj]][age,]
  }
  
  # level 5 (disaggregate by Area & sex)
  for (ikw in 1:47)
  {
    eh_mat[(80+(2*ikw-1)),] = get(mfts_gender_Area_train_residual_female[(ikw+1)])[[kj]][age,]
    eh_mat[(80+(2*ikw)),] = get(mfts_gender_Area_train_residual_male[(ikw+1)])[[kj]][age,]
  }
  
  target = lowerD(t(eh_mat))
  shrink = shrink.estim(t(eh_mat), target)
  wh_mat = shrink[[1]]
  
  return(wh_mat = wh_mat)
}






















