#################################################################################################
#
#                         MODEL DEVELOPMENT USING INTERNAL DATA
#
#################################################################################################

wd = '.....'
code = paste0(wd, '/CODE')
input = paste0(wd, '/INPUT')
output = paste0(wd, '/OUTPUT')

source(paste0(code,'/LGD Functions.R'))

#=================================================================================================
# PREPARE REGRESSION DATASET
#=================================================================================================

#================ Read input data ================#

# Load macro data 
internal = 
  fread(
    paste0(input, '/internal_out.csv'), 
    stringsAsFactors = F, 
    integer64 = 'double', 
    select=c('Quarter','lgd','BROAD_SENIORITY_COLLATERAL','ObligorID')
    )
macro = 
  fread(
    paste0(input, '/macro_out.csv'), 
    stringsAsFactors = F, 
    integer64 = 'double'
    )
macro_adf_res = 
  fread(
    paste0(output, '/macro_adf_out.csv'), 
        stringsAsFactors = F, 
        integer64 = 'double'
    )    # ADF result

setcolorder(internal, c('Quarter','lgd','ObligorID','BROAD_SENIORITY_COLLATERAL'))
setorder(internal, Quarter)    # Sort ascending by Quarter

#================ Stationarity Test for LGD ================#

# Dickey-Fuller test with single mean
lgd.adf.p = 
  adfTest(
    x = internal$lgd, 
    lags = 2, 
    type = 'ct'
    )@test$p.value    # 0.01; stationary

#================ Construct regression data ================#

# Set up time horizon
start = '2007:Q4'
end = '2015:Q2'

# Prepare regression data
internal = 
  internal[
    head(which(Quarter==start),1):tail(which(Quarter==end),1)
    ]
internal = 
  merge(
    internal, 
    macro, 
    by='Quarter', 
    all.x=T, 
    ally=F
    )    # Merge lgd and macro data
internal = internal[!is.na(lgd)]    # Eliminate rows with missing lgd

# Segment internal data into ss and sub
ss = internal[BROAD_SENIORITY_COLLATERAL=='Senior Secured']
sub = internal[BROAD_SENIORITY_COLLATERAL=='Junior/Unsecured']
equip = internal[BROAD_SENIORITY_COLLATERAL=='Equipment']
ss = ss[, BROAD_SENIORITY_COLLATERAL:=NULL]
sub = sub[, BROAD_SENIORITY_COLLATERAL:=NULL]
equip = equip[, BROAD_SENIORITY_COLLATERAL:=NULL]

# Label "lgd" variable to segment-specific lgd
setnames(ss, 'lgd', 'ss_lgd')
setnames(sub, 'lgd', 'sub_lgd')
setnames(equip, 'lgd', 'equip_lgd')

# # Save regression datasets
# write.csv(ss, file=paste0(output, '/ss.csv'), row.names=F, na='')
# write.csv(sub, file=paste0(output, '/sub.csv'), row.names=F, na='')
# write.csv(equip, file=paste0(output, '/equip.csv'), row.names=F, na='')

#=============================================================================================
# UNIVARIATE REGRESSION ANALYSIS
#=============================================================================================

#================ Run univariate regressions ================#

ynames_equip = 'equip_lgd'
ynames_ss = 'ss_lgd'
ynames_sub = 'sub_lgd'
xnames = setdiff(names(macro), 'Quarter')
univ_reg_res_equip = univ_reg(equip, ynames_equip, xnames)
univ_reg_res_ss = univ_reg(ss, ynames_ss, xnames)
univ_reg_res_sub = univ_reg(sub, ynames_sub, xnames)

#================ Correlation Analysis ================#

lgd_vars = c('equip_lgd_4Q_Ma','ss_lgd_4Q_Ma','sub_lgd_4Q_Ma')

# Read (4Q Count weighted MA LGD) created from macro transformation code
internal_ss_sub_macro = 
  fread(
    paste0(input, '/internal_ss_sub_macro_out.csv'),
    stringsAsFactors = F, 
    integer64 = 'double', 
    select=c('Quarter', lgd_vars)
    )

# Merge count weighted 4Q MA LGD into macro data
corr_dt = 
  as.data.table(
    merge(
      macro, 
      internal_ss_sub_macro,
      by='Quarter', 
      all.x=T, 
      all.y=F
      )
    )
corr_dt = 
  corr_dt[
    head(which(Quarter==start),1):tail(which(Quarter==end),1)
    ]
setcolorder(
  corr_dt, 
  c(
    'Quarter', 
    lgd_vars, 
    setdiff(
      names(corr_dt), 
      c('Quarter', lgd_vars)
      )
    )
  )

# Get correlation
corr = 
  reshape2::melt(
    cor(
      corr_dt[, -c('Quarter'), with=F], 
      use = "pairwise.complete.obs", 
      method = "pearson"
      )
    )
corr = as.data.table(corr)
corr = 
  corr[
    which(corr$Var1 %in% lgd_vars & !(corr$Var2 %in% lgd_vars)),
    .(Var1, Var2, value, Sign=sign(value))
    ]
setnames(
  corr, 
  c('LGD_Var','Macro_Var','Correlation','Sign')
  )
setorder(corr, 'LGD_Var')

# Rename variables to make it distinguish from similar vars
setnames(corr, c('LGD_Var','Sign'), c('LGD_4Q_MA','Sign.Corr'))
corr[LGD_4Q_MA=='equip_lgd_4Q_Ma', LGD_4Q_MA:='equip_lgd']
corr[LGD_4Q_MA=='ss_lgd_4Q_Ma', LGD_4Q_MA:='ss_lgd']
corr[LGD_4Q_MA=='sub_lgd_4Q_Ma', LGD_4Q_MA:='sub_lgd']

#================ Merge expected sign info into univariate regression results  ================#

# Merge regression results with expected sign
expected_sign = 
  fread(
    paste0(input, '/con_Var_Tracker.csv'), 
    stringsAsFactors = F, 
    integer64 = 'double',
    select=c('Base','Expected.Sign')
    )
univ_reg_res_equip = 
  merge(
    univ_reg_res_equip, 
    expected_sign, 
    by.x='X.Var', 
    by.y='Base', 
    all.x=T, 
    all.y=F
    )
univ_reg_res_ss = 
  merge(
    univ_reg_res_ss, 
    expected_sign, 
    by.x='X.Var', 
    by.y='Base', 
    all.x=T, 
    all.y=F
    )
univ_reg_res_sub = 
  merge(
    univ_reg_res_sub, 
    expected_sign, 
    by.x='X.Var', 
    by.y='Base', 
    all.x=T, 
    all.y=F
    )

# Make expected sign for transformations same as and consistent with base
load(paste0(input, '/macro_base.Rdata'))
for (macro_name in macro_base) {
  
  sign_temp = 
    expected_sign[
      Base==macro_name, 
      Expected.Sign
      ]
  
  univ_reg_res_equip[
    grep(pattern = macro_name, x = X.Var), 
    Expected.Sign:=sign_temp
    ]
  
  univ_reg_res_ss[
    grep(pattern = macro_name, x = X.Var), 
    Expected.Sign:=sign_temp
    ]
  
  univ_reg_res_sub[
    grep(pattern = macro_name, x = X.Var), 
    Expected.Sign:=sign_temp
    ]
}

#================ Merge correlation results into univariate regression results ================#

# Merge regression results with correlation result
univ_reg_res_equip = 
  merge(
    univ_reg_res_equip, 
    corr, 
    by.x = c('Y.Var','X.Var'), 
    by.y=c('LGD_4Q_MA', 'Macro_Var'),
    all.x=T, 
    all.y=F
    )
univ_reg_res_ss = 
  merge(
    univ_reg_res_ss, 
    corr, 
    by.x = c('Y.Var','X.Var'), 
    by.y=c('LGD_4Q_MA', 'Macro_Var'), 
    all.x=T, 
    all.y=F
    )
univ_reg_res_sub = 
  merge(
    univ_reg_res_sub, 
    corr, 
    by.x = c('Y.Var','X.Var'), 
    by.y=c('LGD_4Q_MA', 'Macro_Var'),
    all.x=T, 
    all.y=F
    )

# Add "Sign.Consistent" flag if Sign.Corr is same as expected
univ_reg_res_equip[, Sign.Consistent := Sign.X==Expected.Sign]
univ_reg_res_ss[, Sign.Consistent := Sign.X==Expected.Sign]
univ_reg_res_sub[, Sign.Consistent := Sign.X==Expected.Sign]

#================ Merge ADF results into univariate regression results ================#

# Merge regression results with adf test results
univ_reg_res_equip = 
  merge(
    univ_reg_res_equip, 
    macro_adf_res[, -'ADF.p', with=F], 
    by.x='X.Var', 
    by.y='Var.Name', 
    all.x=T, 
    all.y=F
    )
univ_reg_res_ss = 
  merge(
    univ_reg_res_ss, 
    macro_adf_res[, -'ADF.p', with=F], 
    by.x='X.Var', 
    by.y='Var.Name', 
    all.x=T, 
    all.y=F
    )
univ_reg_res_sub = 
  merge(
    univ_reg_res_sub, 
    macro_adf_res[, -'ADF.p', with=F], 
    by.x='X.Var',
    by.y='Var.Name', 
    all.x=T, 
    all.y=F
    )

# Remove key from merging step and reset column order to start with Y
setkey(univ_reg_res_equip, NULL)
setkey(univ_reg_res_ss, NULL)
setkey(univ_reg_res_sub, NULL)
setcolorder(
  univ_reg_res_equip, 
  c('Y.Var', setdiff(names(univ_reg_res_equip), 'Y.Var'))
  )
setcolorder(
  univ_reg_res_ss, 
  c('Y.Var', setdiff(names(univ_reg_res_ss), 'Y.Var'))
  )
setcolorder(
  univ_reg_res_sub, 
  c('Y.Var', setdiff(names(univ_reg_res_sub), 'Y.Var'))
  )

#================ Select candidate factors based on univariate results ================#

# Call "select_univ_factors" function to do factor selection
selected_equip_factors = 
  select_univ_factors(
    univ_reg_res_equip, 
    cor.thresh = 0.2, 
    p.thresh = 0.1, 
    adf.p.thresh=0.05
    )
selected_ss_factors = 
  select_univ_factors(
    univ_reg_res_ss, 
    cor.thresh = 0.2, 
    p.thresh = 0.1, 
    adf.p.thresh=0.05
    )
selected_sub_factors = 
  select_univ_factors(
    univ_reg_res_sub, 
    cor.thresh = 0.2, 
    p.thresh = 0.1, 
    adf.p.thresh=0.05
    )

#==============================================================================================
# MULTIVARIATE REGRESSION ANALYSIS
#==============================================================================================

source(paste0(code,'/Functions.R'))

#================ Input parameters ================#

# Declare inputs
yname_equip = 'equip_lgd'
yname_ss = 'ss_lgd'
yname_sub = 'sub_lgd'
xnames_equip = unique(selected_equip_factors$X.Var)
xnames_ss = unique(selected_ss_factors$X.Var)
xnames_sub = unique(selected_sub_factors$X.Var)
factor_subset_list = c(2,3,4)
max_factor_subset = 4
load(paste0(input, '/macro_base.Rdata'))
macro_sign_file = 
  fread(
    paste0(input, '/con_Var_Tracker.csv'), 
    stringsAsFactors = F, 
    integer64 = 'double', 
    select=c('Base','Expected.Sign')
    )
rolling.start.qtr='2007:Q4'
rolling.end.qtr=c('2012:Q4','2013:Q4','2014:Q4')
order = 2
index = c('Quarter','ObligorID')

# Prepare count_weight_file
internal_ss_sub_macro_out = 
  fread(
    paste0(input,'/internal_ss_sub_macro_out.csv'), 
    stringsAsFactors = F, 
    integer64 = 'double'
    )
count_weight_file = 
  internal_ss_sub_macro_out[
    ,
    .(Quarter, count_non_equip=sum(count_ss, count_sub, na.rm=T)), 
    by=1:nrow(internal_ss_sub_macro_out)
    ]
count_weight_file[, nrow:=NULL]
setnames(
  count_weight_file, 
  'count_non_equip', 
  'Count'
  )    # IMPORTANT; colnames must match what is used in function

# Downturn shift factor inputs
downturn_periods=c('2008:Q4','2009:Q1','2009:Q2','2009:Q3')
cycle_start_qtr ='2008:Q1'
cycle_end_qtr = '2015:Q2'
cycle_end_qtr_sub = '2015:Q1'

#================ Get candidate models ================#

# Call "do_multiv_regression" function to get all 2-factor models
internal_best_equip2factors = 
  do_multiv_regression(
    equip, 
    macro_file=macro, 
    macro_base, 
    macro_sign_file,
    yname=yname_equip, 
    xnames=xnames_equip, 
    factor_subset=2, 
    max_factor_subset=4, 
    vif.thresh=5, 
    p.thresh=0.05, 
    hetero.p.thresh=0.05,
    rolling.start.qtr, 
    rolling.end.qtr, 
    K_fold=5,
    count_weight_file, 
    downturn_periods, 
    cycle_start_qtr, 
    cycle_end_qtr, 
    order, 
    index
    )

# Call "do_multiv_regreequipion" function to get all 3-factor models
internal_best_equip3factors = 
  do_multiv_regreequipion(
    equip, 
    macro_file=macro, 
    macro_base, 
    macro_sign_file, 
    yname=yname_equip, 
    xnames=xnames_equip, 
    factor_subset=3, 
    max_factor_subset=4, 
    vif.thresh=5, 
    p.thresh=0.05, 
    hetero.p.thresh=0.05,
    rolling.start.qtr, 
    rolling.end.qtr, 
    K_fold=5,
    count_weight_file, 
    downturn_periods, 
    cycle_start_qtr, 
    cycle_end_qtr,
    order, 
    index
    )

# # Call "do_multiv_regreequipion" function to get all 4-factor models
internal_best_equip4factors = 
  do_multiv_regreequipion(
    equip, 
    macro_file=macro, 
    macro_base, 
    macro_sign_file,
    yname=yname_equip,
    xnames=xnames_equip,
    factor_subset=4, 
    max_factor_subset=4,
    vif.thresh=5,
    p.thresh=0.05,
    hetero.p.thresh=0.05,
    rolling.start.qtr, 
    rolling.end.qtr, 
    K_fold=5,
    count_weight_file, 
    downturn_periods, 
    cycle_start_qtr, 
    cycle_end_qtr,
    order, 
    index
    )

