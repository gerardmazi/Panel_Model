########################################################################################################
#
#                                             MODEL SEGMENTATION ANALYSIS
#
#########################################################################################################

wd = '..........'
code = paste0(wd, '/CODE')
input = paste0(wd, '/INPUT')
output = paste0(wd, '/OUTPUT')

source(paste0(code,'/Functions.R'))

#=========================================================================================================
# KS and MW Tests
#=========================================================================================================

# Read DRD dataset
ur = 
  fread(
    paste0(input,'/output_ur.csv'), 
    stringsAsFactors = F, 
    integer64 = 'double'
    )
ur_all_types = 
  fread(
    paste0(input,'/output_ur_all_types.csv'), 
    stringsAsFactors = F, 
    integer64 = 'double'
    )

#==================== By Industry ====================#

group_var = 'FINAL_INDUSTRY'
industry_matrices = get_ks_mw_matrix(group_var, ur)
industry_ks.p_matrix = industry_matrices$ks.p
industry_mw.p_matrix = industry_matrices$mw.p

#==================== By All Loan Types ====================#
# Note that the data used here is "Moody's data with all types" (No TYPE filter applied)

group_var = 'TYPE'
all_type_matrices = get_ks_mw_matrix(group_var, ur_all_types)
all_type_ks.p_matrix = all_type_matrices$ks.p
all_type_mw.p_matrix = all_type_matrices$mw.p

#==================== By Loan Types ====================#
# Note that the data used here is "Moody's data with all types" (No TYPE filter applied)

group_var = 'TYPE'
type_matrices = get_ks_mw_matrix(group_var, ur)
type_ks.p_matrix = type_matrices$ks.p
type_mw.p_matrix = type_matrices$mw.p

#==================== By Seniority ====================#

group_var = 'SENIORITY'
seniority_matrices = get_ks_mw_matrix(group_var, ur)
seniority_ks.p_matrix = seniority_matrices$ks.p
seniority_mw.p_matrix = seniority_matrices$mw.p

#==================== By Collateral ====================#

group_var = 'COLLATERAL'
collateral_matrices = get_ks_mw_matrix(group_var, ur)
collateral_ks.p_matrix = collateral_matrices$ks.p
collateral_mw.p_matrix = collateral_matrices$mw.p

#==================== By Seniority_Collateral ====================#

group_var = 'SENIORITY_COLLATERAL'
SENIORITY_COLLATERAL_matrices = get_ks_mw_matrix(group_var, ur)
SENIORITY_COLLATERAL_ks.p_matrix = SENIORITY_COLLATERAL_matrices$ks.p
SENIORITY_COLLATERAL_mw.p_matrix = SENIORITY_COLLATERAL_matrices$mw.p

#==================== By Broad_Seniority_Collateral ====================#

group_var = 'BROAD_SENIORITY_COLLATERAL'
BROAD_SENIORITY_COLLATERAL_matrices = get_ks_mw_matrix(group_var, ur)
BROAD_SENIORITY_COLLATERAL_ks.p_matrix = BROAD_SENIORITY_COLLATERAL_matrices$ks.p
BROAD_SENIORITY_COLLATERAL_mw.p_matrix = BROAD_SENIORITY_COLLATERAL_matrices$mw.p

#==================== By Services ====================#

group_var = 'Services'
services_matrices = get_ks_mw_matrix(group_var, ur)
services_ks.p_matrix = services_matrices$ks.p
services_mw.p_matrix = services_matrices$mw.p


#==================================================================================================
# RUN THE FOLLOWING TWO KS/MW TESTS ON ONLY TERM LOANS AND REVOLVERS
#==================================================================================================

#==================== By Collat_Type ====================#

group_var = 'COLLAT_TYPE'
Collat_Type_matrices = get_ks_mw_matrix(group_var, ur)
Collat_Type_ks.p_matrix = Collat_Type_matrices$ks.p
Collat_Type_mw.p_matrix = Collat_Type_matrices$mw.p

#==================== By Seniority_Collat_Type ====================#

group_var = 'SENIORITY_COLLAT_TYPE'
SENIORITY_COLLAT_TYPE_matrices = get_ks_mw_matrix(group_var, ur)
SENIORITY_COLLAT_TYPE_ks.p_matrix = SENIORITY_COLLAT_TYPE_matrices$ks.p
SENIORITY_COLLAT_TYPE_mw.p_matrix = SENIORITY_COLLAT_TYPE_matrices$mw.p


#=================================================================================================
#                                           CORRELATION ANALYSIS
#=================================================================================================

# Import complete set of quarters for constructing a time series
quarters = 
  read.csv(
    paste0(input,'/quarters.csv'),
    stringsAsFactors = F, 
    col.names='Quarter'
    )
source(paste0(code,'/Functions.R'))

#==================== ALL ASSETS ====================#

group_var = 'ALL_ASSETS'
seniority_corr_matrices_p = 
  get_corr_matrix(
    group_var, 
    ur, 
    quarters, 
    'pearson'
    )
seniority_corr_matrices_k = 
  get_corr_matrix(
    group_var, 
    ur,
    quarters, 
    'kendall'
    )
seniority_corr_matrices_s = 
  get_corr_matrix(
    group_var, 
    ur, 
    quarters,
    'spearman'
    )

#==================== TYPE ====================#

group_var = 'TYPE'
seniority_corr_matrices_p = 
  get_corr_matrix(
    group_var, 
    ur, 
    quarters,
    'pearson'
    )
seniority_corr_matrices_k = 
  get_corr_matrix(
    group_var, 
    ur, 
    quarters, 
    'kendall'
    )
seniority_corr_matrices_s = 
  get_corr_matrix(
    group_var, 
    ur, 
    quarters, 
    'spearman'
    )

#==================== SENIORITY ====================#

group_var = 'SENIORITY'
seniority_corr_matrices_p = 
  get_corr_matrix(
    group_var, 
    ur, 
    quarters, 
    'pearson'
    )
seniority_corr_matrices_k = 
  get_corr_matrix(
    group_var, 
    ur, 
    quarters,
    'kendall'
    )
seniority_corr_matrices_s = 
  get_corr_matrix(
    group_var,
    ur,
    quarters,
    'spearman'
    )

#==================== COLLATERAL ====================#

group_var = 'COLLATERAL'
collateral_corr_matrices_p = 
  get_corr_matrix(
    group_var, 
    ur, 
    quarters, 
    'pearson'
    )
collateral_corr_matrices_k =
  get_corr_matrix(
    group_var, 
    ur, 
    quarters,
    'kendall'
    )
collateral_corr_matrices_s =
  get_corr_matrix(
    group_var, 
    ur, 
    quarters,
    'spearman'
    )

#==================== SENIORITY COLLATERAL ====================#

group_var = 'SENIORITY_COLLATERAL'
seniority_collateral_corr_matrices_p = 
  get_corr_matrix(
    group_var, 
    ur,
    quarters, 
    'pearson'
    )
seniority_collateral_corr_matrices_k = 
  get_corr_matrix(
    group_var, 
    ur, 
    quarters,
    'kendall'
    )
seniority_collateral_corr_matrices_s = 
  get_corr_matrix(
    group_var, 
    ur,
    quarters,
    'spearman'
    )

#==================== BROAD SENIORITY COLLATERAL ====================#

group_var = 'BROAD_SENIORITY_COLLATERAL'
BROAD_SENIORITY_COLLATERAL_corr_matrices_p = 
  get_corr_matrix(
    group_var,
    ur,
    quarters,
    'pearson'
    )
BROAD_SENIORITY_COLLATERAL_corr_matrices_k = 
  get_corr_matrix(
    group_var, 
    ur, 
    quarters, 
    'kendall'
    )
BROAD_SENIORITY_COLLATERAL_corr_matrices_s = 
  get_corr_matrix(
    group_var, 
    ur,
    quarters, 
    'spearman'
    )


