#####################################################################
# List of R codes for forecasting subnational mortality in Australia
#####################################################################

###########################################################################
# Import Australian subnational mortality data and define a summing matrix
###########################################################################

# 1_Data_conversion.R: Read in raw Australian subnational mortality and convert the data into demogdata-friendly format.
# 2_Exposure_to_risk.R: forecast exposure to risk at each level of hierarchy of Australian subnational mortality rates.
# 3_Summing_matrix.R: construct summing matrix using forecast exposure to risk.


###############################################
# Univariate point forecast relevant functions
###############################################

# 4_Univariate_point_forecast.R: obtain point forecasts by the grouped univariate forecating method without reconciliation.
# 5_Reconciled_univariate_point_forecast.R: reconcile univariate point forecasts using both bottom-up and optimal combination methods and calculate errors.

#################################################
# Multivariate point forecast relevant functions
#################################################

# 6_Multivariate_point_forecast.R: obtain point forecasts by the grouped mulivariate forecating method without reconciliation.
# 7_Reconciled_multivariate_point_forecast: reconcile multivariate point forecasts using both bottom-up and optimal combination methods and calculate errors.

# 8_Summary_plots.R: plot of results.


