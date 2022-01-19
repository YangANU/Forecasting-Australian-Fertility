#################################################
# Reconciliation of multivariate point forecasts
#################################################

########################################
# Point forecast errors; fmethod = "BU" 
########################################

# create storing objects:

Level_0_err     = Level_F_err       = Level_M_err       = matrix(NA,5,3)
Level_T_R_err   = Level_F_R_err     = Level_M_R_err     = array(NA, dim = c(5,3,11))
Level_State_err = Level_State_F_err = Level_State_M_err = array(NA, dim = c(47,5,3))

for(ikw in 1:5)
{
  dum = BU_optim_err(ik = ikw, hier_method = "BU", 
                     total_data_raw_rate_region = mfts_Region_forc_total,
                     female_data_raw_rate_region = mfts_Region_forc_female,
                     male_data_raw_rate_region = mfts_Region_forc_male,
                     female_data_raw_rate = mfts_Area_forc_female, 
                     male_data_raw_rate = mfts_Area_forc_male, 
                     total_data_raw_rate = mfts_Area_forc_total)
  
  # Total + Sex
  
  Level_0_err[ikw,1] = dum$me_total_err
  Level_F_err[ikw,1] = dum$me_female_err
  Level_M_err[ikw,1] = dum$me_male_err
  
  Level_0_err[ikw,2] = dum$mae_total_err
  Level_F_err[ikw,2] = dum$mae_female_err
  Level_M_err[ikw,2] = dum$mae_male_err
  
  Level_0_err[ikw,3] = dum$rmse_total_err
  Level_F_err[ikw,3] = dum$rmse_female_err
  Level_M_err[ikw,3] = dum$rmse_male_err
  
  # Region + Sex
  
  Level_T_R_err[ikw,1,] = dum$me_total_R_err
  Level_F_R_err[ikw,1,] = dum$me_female_R_err
  Level_M_R_err[ikw,1,] = dum$me_male_R_err
  
  Level_T_R_err[ikw,2,] = dum$mae_total_R_err
  Level_F_R_err[ikw,2,] = dum$mae_female_R_err
  Level_M_R_err[ikw,2,] = dum$mae_male_R_err
  
  Level_T_R_err[ikw,3,] = dum$rmse_total_R_err
  Level_F_R_err[ikw,3,] = dum$rmse_female_R_err
  Level_M_R_err[ikw,3,] = dum$rmse_male_R_err
  
  # Prefecture + Sex
  
  Level_State_err[,ikw,1]   = dum$me_Area_err
  Level_State_F_err[,ikw,1] = dum$me_bottom_female_err
  Level_State_M_err[,ikw,1] = dum$me_bottom_male_err
  
  Level_State_err[,ikw,2]   = dum$mae_Area_err
  Level_State_F_err[,ikw,2] = dum$mae_bottom_female_err
  Level_State_M_err[,ikw,2] = dum$mae_bottom_male_err
  
  Level_State_err[,ikw,3]   = dum$rmse_Area_err
  Level_State_F_err[,ikw,3] = dum$rmse_bottom_female_err
  Level_State_M_err[,ikw,3] = dum$rmse_bottom_male_err
}

# summary of results

BU_all_level_err_me_mfts = cbind(Level_0_err[,1], rowMeans(cbind(Level_F_err[,1], Level_M_err[,1])),
                            rowMeans(Level_T_R_err[,1,]), rowMeans(cbind(rowMeans(Level_F_R_err[,1,]), rowMeans(Level_M_R_err[,1,]))),
                            colMeans(Level_State_err[,,1]), colMeans(rbind(Level_State_F_err[,,1], Level_State_M_err[,,1])))

BU_all_level_err_mae_mfts = cbind(Level_0_err[,2], rowMeans(cbind(Level_F_err[,2], Level_M_err[,2])),
                             rowMeans(Level_T_R_err[,2,]), rowMeans(cbind(rowMeans(Level_F_R_err[,2,]), rowMeans(Level_M_R_err[,2,]))),
                             colMeans(Level_State_err[,,2]), colMeans(rbind(Level_State_F_err[,,2], Level_State_M_err[,,2])))

BU_all_level_err_rmse_mfts = cbind(Level_0_err[,3], rowMeans(cbind(Level_F_err[,3], Level_M_err[,3])),
                              rowMeans(Level_T_R_err[,3,]), rowMeans(cbind(rowMeans(Level_F_R_err[,3,]), rowMeans(Level_M_R_err[,3,]))),
                              colMeans(Level_State_err[,,3]), colMeans(rbind(Level_State_F_err[,,3], Level_State_M_err[,,3])))
colnames(BU_all_level_err_me_mfts) = colnames(BU_all_level_err_mae_mfts) = colnames(BU_all_level_err_rmse_mfts) = c("Total", "Sex", "Region", "Region + Sex", "Area", "Area + Sex")

##############################################
# Point forecast errors; fmethod = "comb_OLS" 
##############################################

# create storing objects:

optim_Level_0_err     = optim_Level_F_err       = optim_Level_M_err       = matrix(NA,5,3)
optim_Level_T_R_err   = optim_Level_F_R_err     = optim_Level_M_R_err     = array(NA, dim = c(5,3,11))
optim_Level_State_err = optim_Level_State_F_err = optim_Level_State_M_err = array(NA, dim = c(47,5,3))

for(ikw in 1:5)
{
  dum = BU_optim_err(ik = ikw, hier_method = "comb_OLS", 
                     total_data_raw_rate_region = mfts_Region_forc_total,
                     female_data_raw_rate_region = mfts_Region_forc_female,
                     male_data_raw_rate_region = mfts_Region_forc_male,
                     female_data_raw_rate = mfts_Area_forc_female, 
                     male_data_raw_rate = mfts_Area_forc_male, 
                     total_data_raw_rate = mfts_Area_forc_total)
  
  # Total + Sex
  
  optim_Level_0_err[ikw,1] = dum$me_total_err
  optim_Level_F_err[ikw,1] = dum$me_female_err
  optim_Level_M_err[ikw,1] = dum$me_male_err
  
  optim_Level_0_err[ikw,2] = dum$mae_total_err
  optim_Level_F_err[ikw,2] = dum$mae_female_err
  optim_Level_M_err[ikw,2] = dum$mae_male_err
  
  optim_Level_0_err[ikw,3] = dum$rmse_total_err
  optim_Level_F_err[ikw,3] = dum$rmse_female_err
  optim_Level_M_err[ikw,3] = dum$rmse_male_err
  
  # Region + Sex
  
  optim_Level_T_R_err[ikw,1,] = dum$me_total_R_err
  optim_Level_F_R_err[ikw,1,] = dum$me_female_R_err
  optim_Level_M_R_err[ikw,1,] = dum$me_male_R_err
  
  optim_Level_T_R_err[ikw,2,] = dum$mae_total_R_err
  optim_Level_F_R_err[ikw,2,] = dum$mae_female_R_err
  optim_Level_M_R_err[ikw,2,] = dum$mae_male_R_err
  
  optim_Level_T_R_err[ikw,3,] = dum$rmse_total_R_err
  optim_Level_F_R_err[ikw,3,] = dum$rmse_female_R_err
  optim_Level_M_R_err[ikw,3,] = dum$rmse_male_R_err
  
  # Prefecture + Sex
  
  optim_Level_State_err[,ikw,1]   = dum$me_Area_err
  optim_Level_State_F_err[,ikw,1] = dum$me_bottom_female_err
  optim_Level_State_M_err[,ikw,1] = dum$me_bottom_male_err
  
  optim_Level_State_err[,ikw,2]   = dum$mae_Area_err
  optim_Level_State_F_err[,ikw,2] = dum$mae_bottom_female_err
  optim_Level_State_M_err[,ikw,2] = dum$mae_bottom_male_err
  
  optim_Level_State_err[,ikw,3]   = dum$rmse_Area_err
  optim_Level_State_F_err[,ikw,3] = dum$rmse_bottom_female_err
  optim_Level_State_M_err[,ikw,3] = dum$rmse_bottom_male_err
}

# summary of results

optim_all_level_err_me_mfts = cbind(optim_Level_0_err[,1], rowMeans(cbind(optim_Level_F_err[,1], optim_Level_M_err[,1])),
                               rowMeans(optim_Level_T_R_err[,1,]), rowMeans(cbind(rowMeans(optim_Level_F_R_err[,1,]), rowMeans(optim_Level_M_R_err[,1,]))),
                               colMeans(optim_Level_State_err[,,1]), colMeans(rbind(optim_Level_State_F_err[,,1], optim_Level_State_M_err[,,1])))

optim_all_level_err_mae_mfts = cbind(optim_Level_0_err[,2], rowMeans(cbind(optim_Level_F_err[,2], optim_Level_M_err[,2])),
                                rowMeans(optim_Level_T_R_err[,2,]), rowMeans(cbind(rowMeans(optim_Level_F_R_err[,2,]), rowMeans(optim_Level_M_R_err[,2,]))),
                                colMeans(optim_Level_State_err[,,2]), colMeans(rbind(optim_Level_State_F_err[,,2], optim_Level_State_M_err[,,2])))

optim_all_level_err_rmse_mfts = cbind(optim_Level_0_err[,3], rowMeans(cbind(optim_Level_F_err[,3], optim_Level_M_err[,3])),
                                 rowMeans(optim_Level_T_R_err[,3,]), rowMeans(cbind(rowMeans(optim_Level_F_R_err[,3,]), rowMeans(optim_Level_M_R_err[,3,]))),
                                 colMeans(optim_Level_State_err[,,3]), colMeans(rbind(optim_Level_State_F_err[,,3], optim_Level_State_M_err[,,3])))
colnames(optim_all_level_err_me_mfts) = colnames(optim_all_level_err_mae_mfts) = colnames(optim_all_level_err_rmse_mfts) = c("Total", "Sex", "Region", "Region + Sex", "Area", "Area + Sex")

###############################################
# Point forecast errors; fmethod = "mint_mfts" 
###############################################

# create storing objects:

optim_mint_Level_0_err     = optim_mint_Level_F_err       = optim_mint_Level_M_err       = matrix(NA,5,3)
optim_mint_Level_T_R_err   = optim_mint_Level_F_R_err     = optim_mint_Level_M_R_err     = array(NA, dim = c(5,3,11))
optim_mint_Level_State_err = optim_mint_Level_State_F_err = optim_mint_Level_State_M_err = array(NA, dim = c(47,5,3))

for(ikw in 1:5)
{
  dum = BU_optim_err(ik = ikw, hier_method = "mint_mfts", 
                     total_data_raw_rate_region = mfts_Region_forc_total,
                     female_data_raw_rate_region = mfts_Region_forc_female,
                     male_data_raw_rate_region = mfts_Region_forc_male,
                     female_data_raw_rate = mfts_Area_forc_female, 
                     male_data_raw_rate = mfts_Area_forc_male, 
                     total_data_raw_rate = mfts_Area_forc_total)
  
  # Total + Sex
  
  optim_mint_Level_0_err[ikw,1] = dum$me_total_err
  optim_mint_Level_F_err[ikw,1] = dum$me_female_err
  optim_mint_Level_M_err[ikw,1] = dum$me_male_err
  
  optim_mint_Level_0_err[ikw,2] = dum$mae_total_err
  optim_mint_Level_F_err[ikw,2] = dum$mae_female_err
  optim_mint_Level_M_err[ikw,2] = dum$mae_male_err
  
  optim_mint_Level_0_err[ikw,3] = dum$rmse_total_err
  optim_mint_Level_F_err[ikw,3] = dum$rmse_female_err
  optim_mint_Level_M_err[ikw,3] = dum$rmse_male_err
  
  # Region + Sex
  
  optim_mint_Level_T_R_err[ikw,1,] = dum$me_total_R_err
  optim_mint_Level_F_R_err[ikw,1,] = dum$me_female_R_err
  optim_mint_Level_M_R_err[ikw,1,] = dum$me_male_R_err
  
  optim_mint_Level_T_R_err[ikw,2,] = dum$mae_total_R_err
  optim_mint_Level_F_R_err[ikw,2,] = dum$mae_female_R_err
  optim_mint_Level_M_R_err[ikw,2,] = dum$mae_male_R_err
  
  optim_mint_Level_T_R_err[ikw,3,] = dum$rmse_total_R_err
  optim_mint_Level_F_R_err[ikw,3,] = dum$rmse_female_R_err
  optim_mint_Level_M_R_err[ikw,3,] = dum$rmse_male_R_err
  
  # Prefecture + Sex
  
  optim_mint_Level_State_err[,ikw,1]   = dum$me_Area_err
  optim_mint_Level_State_F_err[,ikw,1] = dum$me_bottom_female_err
  optim_mint_Level_State_M_err[,ikw,1] = dum$me_bottom_male_err
  
  optim_mint_Level_State_err[,ikw,2]   = dum$mae_Area_err
  optim_mint_Level_State_F_err[,ikw,2] = dum$mae_bottom_female_err
  optim_mint_Level_State_M_err[,ikw,2] = dum$mae_bottom_male_err
  
  optim_mint_Level_State_err[,ikw,3]   = dum$rmse_Area_err
  optim_mint_Level_State_F_err[,ikw,3] = dum$rmse_bottom_female_err
  optim_mint_Level_State_M_err[,ikw,3] = dum$rmse_bottom_male_err
}

# summary of results

optim_mint_all_level_err_me_mfts = cbind(optim_mint_Level_0_err[,1], rowMeans(cbind(optim_mint_Level_F_err[,1], optim_mint_Level_M_err[,1])),
                                    rowMeans(optim_mint_Level_T_R_err[,1,]), rowMeans(cbind(rowMeans(optim_mint_Level_F_R_err[,1,]), rowMeans(optim_mint_Level_M_R_err[,1,]))),
                                    colMeans(optim_mint_Level_State_err[,,1]), colMeans(rbind(optim_mint_Level_State_F_err[,,1], optim_mint_Level_State_M_err[,,1])))

optim_mint_all_level_err_mae_mfts = cbind(optim_mint_Level_0_err[,2], rowMeans(cbind(optim_mint_Level_F_err[,2], optim_mint_Level_M_err[,2])),
                                     rowMeans(optim_mint_Level_T_R_err[,2,]), rowMeans(cbind(rowMeans(optim_mint_Level_F_R_err[,2,]), rowMeans(optim_mint_Level_M_R_err[,2,]))),
                                     colMeans(optim_mint_Level_State_err[,,2]), colMeans(rbind(optim_mint_Level_State_F_err[,,2], optim_mint_Level_State_M_err[,,2])))

optim_mint_all_level_err_rmse_mfts = cbind(optim_mint_Level_0_err[,3], rowMeans(cbind(optim_mint_Level_F_err[,3], optim_mint_Level_M_err[,3])),
                                      rowMeans(optim_mint_Level_T_R_err[,3,]), rowMeans(cbind(rowMeans(optim_mint_Level_F_R_err[,3,]), rowMeans(optim_mint_Level_M_R_err[,3,]))),
                                      colMeans(optim_mint_Level_State_err[,,3]), colMeans(rbind(optim_mint_Level_State_F_err[,,3], optim_mint_Level_State_M_err[,,3])))
colnames(optim_mint_all_level_err_me_mfts) = colnames(optim_mint_all_level_err_mae_mfts) = colnames(optim_mint_all_level_err_rmse_mfts) = c("Total", "Sex", "Region", "Region + Sex", "Area", "Area + Sex")

##################################################
# Simple averaging of reconcilied point forecasts
##################################################

Comb_simple_average_mfts = function(total_data_raw_rate_region, female_data_raw_rate_region,
                               male_data_raw_rate_region, female_data_raw_rate, male_data_raw_rate,
                               total_data_raw_rate, ik, comb_method = c("Av", "AvInt"))
{
  comb_method = match.arg(comb_method)
  
  Comb_simple_average_hier = array(NA, dim = c(18,177,(6-ik)))
  if(comb_method == "Av")
  {
    for(ik_age in 1:18)
    {
      Comb_simple_average_hier[ik_age,,] = (BU_optim_hier_mfts(total_data_raw_rate_region, female_data_raw_rate_region,
                                                               male_data_raw_rate_region,female_data_raw_rate, 
                                                               male_data_raw_rate, total_data_raw_rate,      
                                                               age = ik_age, kj = ik, hier_method = "BU") + 
                                              BU_optim_hier_mfts(total_data_raw_rate_region, female_data_raw_rate_region,
                                                                 male_data_raw_rate_region,female_data_raw_rate, 
                                                                 male_data_raw_rate, total_data_raw_rate,      
                                                                 age = ik_age, kj = ik, hier_method = "comb_OLS") +
                                              BU_optim_hier_mfts(total_data_raw_rate_region, female_data_raw_rate_region,
                                                                 male_data_raw_rate_region,female_data_raw_rate, 
                                                                 male_data_raw_rate, total_data_raw_rate,      
                                                                 age = ik_age, kj = ik, hier_method = "mint_mfts"))/3
    }
  } 
  
  if(comb_method == "AvInt")
  {
    stop("Method needs update")
  }
  
  #####################################
  # Errors, including ME, MAE and RMSE
  #####################################
  me = ftsa:::me; mae = ftsa:::mae; rmse = ftsa:::rmse
  
  # Level 0 (Total)
  me_total_err    = me(Comb_simple_average_hier[,1,],    extract.years(get(Area[1]), years = (2011+ik):2016)$rate$total)
  mae_total_err   = mae(Comb_simple_average_hier[,1,],   extract.years(get(Area[1]), years = (2011+ik):2016)$rate$total)
  rmse_total_err  = rmse(Comb_simple_average_hier[,1,],  extract.years(get(Area[1]), years = (2011+ik):2016)$rate$total)
  
  # Level 1 (Sex)
  me_female_err  = me(Comb_simple_average_hier[,2,],  extract.years(get(Area[1]), years = (2011+ik):2016)$rate$female)
  me_male_err    = me(Comb_simple_average_hier[,3,],  extract.years(get(Area[1]), years = (2011+ik):2016)$rate$male)
  
  mae_female_err = mae(Comb_simple_average_hier[,2,],  extract.years(get(Area[1]), years = (2011+ik):2016)$rate$female)
  mae_male_err   = mae(Comb_simple_average_hier[,3,],  extract.years(get(Area[1]), years = (2011+ik):2016)$rate$male)
  
  rmse_female_err = rmse(Comb_simple_average_hier[,2,],  extract.years(get(Area[1]), years = (2011+ik):2016)$rate$female)
  rmse_male_err   = rmse(Comb_simple_average_hier[,3,],  extract.years(get(Area[1]), years = (2011+ik):2016)$rate$male)
  
  # Level 2 (Region + Total, Female and Male)
  me_total_R_err = me_female_R_err = me_male_R_err = mae_total_R_err = mae_female_R_err = mae_male_R_err = rmse_total_R_err = rmse_female_R_err = rmse_male_R_err = vector("numeric",1)
  for(iw in 1:11)
  {
    data_tesing_total = extract.years(get(Region[iw]), years = (2011+ik):2016)$rate$total
    testing_ind_total = data_tesing_total != 0
    
    data_tesing_female = extract.years(get(Region[iw]), years = (2011+ik):2016)$rate$female
    testing_ind_female = data_tesing_female != 0
    
    data_tesing_male = extract.years(get(Region[iw]), years = (2011+ik):2016)$rate$male
    testing_ind_male = data_tesing_male != 0
    
    me_total_R_err[iw]  = me(Comb_simple_average_hier[,(iw+3),][testing_ind_total],   data_tesing_total[testing_ind_total])
    me_female_R_err[iw] = me(Comb_simple_average_hier[,(iw+14),][testing_ind_female], data_tesing_female[testing_ind_female])
    me_male_R_err[iw]   = me(Comb_simple_average_hier[,(iw+25),][testing_ind_male],   data_tesing_male[testing_ind_male])
    
    mae_total_R_err[iw]  = mae(Comb_simple_average_hier[,(iw+3),][testing_ind_total],   data_tesing_total[testing_ind_total])
    mae_female_R_err[iw] = mae(Comb_simple_average_hier[,(iw+14),][testing_ind_female], data_tesing_female[testing_ind_female])
    mae_male_R_err[iw]   = mae(Comb_simple_average_hier[,(iw+25),][testing_ind_male],   data_tesing_male[testing_ind_male])
    
    rmse_total_R_err[iw]  = rmse(Comb_simple_average_hier[,(iw+3),][testing_ind_total],   data_tesing_total[testing_ind_total])
    rmse_female_R_err[iw] = rmse(Comb_simple_average_hier[,(iw+14),][testing_ind_female], data_tesing_female[testing_ind_female])
    rmse_male_R_err[iw]   = rmse(Comb_simple_average_hier[,(iw+25),][testing_ind_male],   data_tesing_male[testing_ind_male])
  }
  
  # Level 3 (Prefecture + Total)
  
  me_Area_err = bottom_female_me =  bottom_male_me = mae_Area_err = bottom_female_mae =  bottom_male_mae = rmse_Area_err = bottom_female_rmse = bottom_male_rmse = vector("numeric",47)
  
  for(iwk in 2:48)
  {
    data_testing_total = extract.years(get(Area[iwk]), years = (2011+ik):2016)$rate$total
    testing_ind_total  = data_testing_total != 0
    
    me_Area_err[iwk-1]    = me(Comb_simple_average_hier[,(iwk+35),][testing_ind_total],    data_testing_total[testing_ind_total])
    mae_Area_err[iwk-1]   = mae(Comb_simple_average_hier[,(iwk+35),][testing_ind_total],   data_testing_total[testing_ind_total])
    rmse_Area_err[iwk-1]  = rmse(Comb_simple_average_hier[,(iwk+35),][testing_ind_total],  data_testing_total[testing_ind_total])
  }
  
  for(iwk in 2:48)
  {
    data_testing_female = extract.years(get(Area[iwk]), years = (2011+ik):2016)$rate$female
    testing_ind_female  = data_testing_female != 0
    
    data_testing_male = extract.years(get(Area[iwk]), years = (2011+ik):2016)$rate$male
    testing_ind_male  = data_testing_male != 0
    
    bottom_female_me[iwk-1]   = me(Comb_simple_average_hier[,(80+2*iwk),][testing_ind_female],    data_testing_female[testing_ind_female])
    bottom_female_mae[iwk-1]  = mae(Comb_simple_average_hier[,(80+2*iwk),][testing_ind_female],   data_testing_female[testing_ind_female])
    bottom_female_rmse[iwk-1] = rmse(Comb_simple_average_hier[,(80+2*iwk),][testing_ind_female],  data_testing_female[testing_ind_female])
    
    bottom_male_me[iwk-1]     = me(Comb_simple_average_hier[,(81+2*iwk),][testing_ind_male],    data_testing_male[testing_ind_male])
    bottom_male_mae[iwk-1]    = mae(Comb_simple_average_hier[,(81+2*iwk),][testing_ind_male],   data_testing_male[testing_ind_male])
    bottom_male_rmse[iwk-1]   = rmse(Comb_simple_average_hier[,(81+2*iwk),][testing_ind_male],  data_testing_male[testing_ind_male])
  }      
  
  return(list(me_total_err = me_total_err, me_female_err = me_female_err, 
              me_male_err  = me_male_err,  me_total_R_err = me_total_R_err,
              me_female_R_err = me_female_R_err, me_male_R_err = me_male_R_err,
              me_Area_err  = me_Area_err, me_bottom_female_err = bottom_female_me, 
              me_bottom_male_err = bottom_male_me, 
              
              mae_total_err = mae_total_err, mae_female_err = mae_female_err, 
              mae_male_err  = mae_male_err,  mae_total_R_err = mae_total_R_err,
              mae_female_R_err = mae_female_R_err, mae_male_R_err = mae_male_R_err,
              mae_Area_err  = mae_Area_err, mae_bottom_female_err = bottom_female_mae, 
              mae_bottom_male_err = bottom_male_mae, 
              
              rmse_total_err = rmse_total_err, rmse_female_err = rmse_female_err, 
              rmse_male_err  = rmse_male_err,  
              rmse_total_R_err = rmse_total_R_err,
              rmse_female_R_err = rmse_female_R_err,
              rmse_male_R_err = rmse_male_R_err,
              rmse_Area_err  = rmse_Area_err,
              rmse_bottom_female_err = bottom_female_rmse, rmse_bottom_male_err = bottom_male_rmse))
  
}

###########################################
# Point forecast errors; comb_method= "Av" 
###########################################

# create storing objects:

Comb_Av_Level_0_err     = Comb_Av_Level_F_err       = Comb_Av_Level_M_err       = matrix(NA,5,3)
Comb_Av_Level_T_R_err   = Comb_Av_Level_F_R_err     = Comb_Av_Level_M_R_err     = array(NA, dim = c(5,3,11))
Comb_Av_Level_State_err = Comb_Av_Level_State_F_err = Comb_Av_Level_State_M_err = array(NA, dim = c(47,5,3))

for(ikw in 1:5)
{
  dum = Comb_simple_average_mfts(ik = ikw, comb_method = "Av", 
                            total_data_raw_rate_region = mfts_Region_forc_total,
                            female_data_raw_rate_region = mfts_Region_forc_female,
                            male_data_raw_rate_region = mfts_Region_forc_male,
                            female_data_raw_rate = mfts_Area_forc_female, 
                            male_data_raw_rate = mfts_Area_forc_male, 
                            total_data_raw_rate = mfts_Area_forc_total)
  
  # Total + Sex
  
  Comb_Av_Level_0_err[ikw,1] = dum$me_total_err
  Comb_Av_Level_F_err[ikw,1] = dum$me_female_err
  Comb_Av_Level_M_err[ikw,1] = dum$me_male_err
  
  Comb_Av_Level_0_err[ikw,2] = dum$mae_total_err
  Comb_Av_Level_F_err[ikw,2] = dum$mae_female_err
  Comb_Av_Level_M_err[ikw,2] = dum$mae_male_err
  
  Comb_Av_Level_0_err[ikw,3] = dum$rmse_total_err
  Comb_Av_Level_F_err[ikw,3] = dum$rmse_female_err
  Comb_Av_Level_M_err[ikw,3] = dum$rmse_male_err
  
  # Region + Sex
  
  Comb_Av_Level_T_R_err[ikw,1,] = dum$me_total_R_err
  Comb_Av_Level_F_R_err[ikw,1,] = dum$me_female_R_err
  Comb_Av_Level_M_R_err[ikw,1,] = dum$me_male_R_err
  
  Comb_Av_Level_T_R_err[ikw,2,] = dum$mae_total_R_err
  Comb_Av_Level_F_R_err[ikw,2,] = dum$mae_female_R_err
  Comb_Av_Level_M_R_err[ikw,2,] = dum$mae_male_R_err
  
  Comb_Av_Level_T_R_err[ikw,3,] = dum$rmse_total_R_err
  Comb_Av_Level_F_R_err[ikw,3,] = dum$rmse_female_R_err
  Comb_Av_Level_M_R_err[ikw,3,] = dum$rmse_male_R_err
  
  # Prefecture + Sex
  
  Comb_Av_Level_State_err[,ikw,1]   = dum$me_Area_err
  Comb_Av_Level_State_F_err[,ikw,1] = dum$me_bottom_female_err
  Comb_Av_Level_State_M_err[,ikw,1] = dum$me_bottom_male_err
  
  Comb_Av_Level_State_err[,ikw,2]   = dum$mae_Area_err
  Comb_Av_Level_State_F_err[,ikw,2] = dum$mae_bottom_female_err
  Comb_Av_Level_State_M_err[,ikw,2] = dum$mae_bottom_male_err
  
  Comb_Av_Level_State_err[,ikw,3]   = dum$rmse_Area_err
  Comb_Av_Level_State_F_err[,ikw,3] = dum$rmse_bottom_female_err
  Comb_Av_Level_State_M_err[,ikw,3] = dum$rmse_bottom_male_err
  
  rm(dum)
}

# summary of results

Comb_Av_all_level_err_me_mfts = cbind(Comb_Av_Level_0_err[,1], rowMeans(cbind(Comb_Av_Level_F_err[,1], Comb_Av_Level_M_err[,1])),
                                     rowMeans(Comb_Av_Level_T_R_err[,1,]), rowMeans(cbind(rowMeans(Comb_Av_Level_F_R_err[,1,]), rowMeans(Comb_Av_Level_M_R_err[,1,]))),
                                     colMeans(Comb_Av_Level_State_err[,,1]), colMeans(rbind(Comb_Av_Level_State_F_err[,,1], Comb_Av_Level_State_M_err[,,1])))

Comb_Av_all_level_err_mae_mfts = cbind(Comb_Av_Level_0_err[,2], rowMeans(cbind(Comb_Av_Level_F_err[,2], Comb_Av_Level_M_err[,2])),
                                      rowMeans(Comb_Av_Level_T_R_err[,2,]), rowMeans(cbind(rowMeans(Comb_Av_Level_F_R_err[,2,]), rowMeans(Comb_Av_Level_M_R_err[,2,]))),
                                      colMeans(Comb_Av_Level_State_err[,,2]), colMeans(rbind(Comb_Av_Level_State_F_err[,,2], Comb_Av_Level_State_M_err[,,2])))

Comb_Av_all_level_err_rmse_mfts = cbind(Comb_Av_Level_0_err[,3], rowMeans(cbind(Comb_Av_Level_F_err[,3], Comb_Av_Level_M_err[,3])),
                                       rowMeans(Comb_Av_Level_T_R_err[,3,]), rowMeans(cbind(rowMeans(Comb_Av_Level_F_R_err[,3,]), rowMeans(Comb_Av_Level_M_R_err[,3,]))),
                                       colMeans(Comb_Av_Level_State_err[,,3]), colMeans(rbind(Comb_Av_Level_State_F_err[,,3], Comb_Av_Level_State_M_err[,,3])))
colnames(Comb_Av_all_level_err_me_mfts) = colnames(Comb_Av_all_level_err_mae_mfts) = colnames(Comb_Av_all_level_err_rmse_mfts) = c("Total", "Sex", "Region", "Region + Sex", "Area", "Area + Sex")



########################################
# Summary of reconciled point forecasts
########################################

point_accuracy_me_mfts   = cbind(mfts_all_level_err_me,   BU_all_level_err_me_mfts,   optim_all_level_err_me_mfts,   optim_mint_all_level_err_me_mfts,   Comb_Av_all_level_err_me_mfts)
point_accuracy_mae_mfts  = cbind(mfts_all_level_err_mae,  BU_all_level_err_mae_mfts,  optim_all_level_err_mae_mfts,  optim_mint_all_level_err_mae_mfts,  Comb_Av_all_level_err_mae_mfts)
point_accuracy_rmse_mfts = cbind(mfts_all_level_err_rmse, BU_all_level_err_rmse_mfts, optim_all_level_err_rmse_mfts, optim_mint_all_level_err_rmse_mfts, Comb_Av_all_level_err_rmse_mfts)

point_forecast_accuracy_me_mfts   = rbind(point_accuracy_me_mfts,   apply(point_accuracy_me_mfts,   2, mean), apply(point_accuracy_me_mfts,   2, median))
point_forecast_accuracy_mae_mfts = rbind(point_accuracy_mae_mfts,  apply(point_accuracy_mae_mfts,  2, mean), apply(point_accuracy_mae_mfts,  2, median))
point_forecast_accuracy_rmse_mfts = rbind(point_accuracy_rmse_mfts, apply(point_accuracy_rmse_mfts, 2, mean), apply(point_accuracy_rmse_mfts, 2, median))

rownames(point_forecast_accuracy_me_mfts) = rownames(point_forecast_accuracy_mae_mfts) = rownames(point_forecast_accuracy_rmse_mfts) = c(1:5, "Mean", "Median")









































