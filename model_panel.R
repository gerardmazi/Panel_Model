###################################################################################################
#
#                                                 MODEL DEVELOPMENT CODE
#
####################################################################################################

wd = '.......'
code = paste0(wd, '/CODE')
input = paste0(wd, '/INPUT')
output = paste0(wd, '/OUTPUT')

source(paste0(code,'/Functions.R'))


#=====================================================================================================
# PREPARE REGRESSION DATASET
#=====================================================================================================

#================ Read input data ================#

ur = 
  fread(
    paste0(input, '/output_ur.csv'), 
    stringsAsFactors = F, 
    integer64 = 'double'
    )    # Ultimate recovery
macro = 
  fread(
    paste0(input,'/macro_out.csv'), 
    stringsAsFactors = F, 
    integer64 = 'double'
    )    # Transformed macro

#================ Stationarity Test  ================#

# Dickey-Fuller test with single mean
lgd.adf.p = 
  adfTest(
    x = ur$lgd,
    lags = 2, 
    type = 'ct'
    )@test$p.value    # 0.01; stationary

#================ Stationarity Test for Macroeconomic Variables ===========#

# Test for stationarity using Augmented Dickey Fuller test
adf.p = list()
k = 2    # Lag order
for (i in 2:length(names(macro))) { # Start from 2 b/c the first column is "Quarter"
  nm = names(macro)[i]
  temp_series = macro[[nm]]
  temp_series = temp_series[!is.na(temp_series)]
  # Note below "nc" means "no constant" or "zero mean"; 
  # "c" means "constant" or "single mean"; 
  #"ct" means "constant
  # with time trend"
  adf.p[[i-1]] = 
    data.table(
      as.numeric(
        adfTest(
          temp_series,
          lags=k, 
          type='nc'
          )@test$statistic),
      as.numeric(
        adfTest(
          temp_series,
          lags=k, 
          type='nc'
          )@test$p.value
        ),
      as.numeric(
        adfTest(
          temp_series,
          lags=k, 
          type='c'
          )@test$statistic
        ),
      as.numeric(
        adfTest(
          temp_series, 
          lags=k,
          type='c'
          )@test$p.value
        ),
      as.numeric(
        adfTest(
          temp_series, 
          lags=k, 
          type='ct'
          )@test$statistic
        ),
      as.numeric(
        adfTest(
          temp_series,
          lags=k, 
          type='ct'
          )@test$p.value
        )
      )
  rm(temp_series, nm)
}
macro_adf_res = 
  data.table(
    setdiff(names(macro),'Quarter'),
    rbindlist(adf.p)
    )
setnames(
  macro_adf_res, 
  c('Var.Name','ADF.nc.DF.stat','ADF.nc.p','ADF.c.DF.stat','ADF.c.p','ADF.ct.DF.stat','ADF.ct.p')
  )

# Assume type of ADF test for all variables is "Zero Mean"
macro_adf_res[, ADF.Type := 'Zero Mean']
macro_adf_res[
  , 
  ADF.p:=ifelse(
    ADF.Type=='Zero Mean', 
    ADF.nc.p, 
    ifelse(
      ADF.Type=='Single Mean', 
      ADF.c.p, 
      ADF.ct.p
      )
    ),
  by=1:nrow(macro_adf_res)
  ]

# Save adf test results
write.csv(
  macro_adf_res,
  file=paste0(output, '/macro_adf_out.csv'), 
  row.names=F,
  na=''
  )

#================ Construct regression data ================#

# Set up time horizon
start = '1996:Q1'
end = '2016:Q2'

# Prepare regression data
ur = ur[, .(Quarter, INSTRUMENT_ID, lgd, BROAD_SENIORITY_COLLATERAL)]    # Keep only Quarter and lgd
setorder(ur, Quarter)    # Sort ascending by Quarter
ur = ur[head(which(Quarter==start),1):tail(which(Quarter==end),1)]    # Keep only declared  horizon
ur = merge(ur, macro, by='Quarter', all.x=T, ally=F)    # Merge data
ur = ur[!is.na(lgd)]    # Eliminate rows with missing lgd

# Segment ur data into ss and sub
ss = ur[BROAD_SENIORITY_COLLATERAL=='Senior Secured']
sub = ur[BROAD_SENIORITY_COLLATERAL=='Junior/Unsecured']
ss = ss[, BROAD_SENIORITY_COLLATERAL:=NULL]
sub = sub[, BROAD_SENIORITY_COLLATERAL:=NULL]

# Label "lgd" variable to segment-specific lgd
setnames(ss, 'lgd', 'ss_lgd')
setnames(sub, 'lgd', 'sub_lgd')

#====================================================================================================
# UNIVARIATE REGRESSION ANALYSIS
#====================================================================================================

#================ Run univariate regressions ================#

ynames_ss = 'ss_lgd'
ynames_sub = 'sub_lgd'
xnames = setdiff(names(macro), 'Quarter')
univ_reg_res_ss = univ_reg(ss, ynames_ss, xnames)
univ_reg_res_sub = univ_reg(sub, ynames_sub, xnames)

#================ Correlation Analysis ================#

# Read (4Q Count weighted MA LGD) created from macro transformation code
ss_sub_macro = 
  fread(
    paste0(input, '/ss_sub_macro_out.csv'), 
    stringsAsFactors = F,
    integer64 = 'double'
    )

# Merge count weighted 4Q MA LGD into macro data
corr_dt = 
  merge(
    macro, 
    ss_sub_macro[, .(Quarter, ss_lgd_4Q_Ma, sub_lgd_4Q_Ma)], 
    by='Quarter', 
    all.x=T,
    all.y=F
    )
corr_dt = as.data.table(corr_dt)

# Subset corr_dt to the desired time horizon
corr_dt = 
  corr_dt[
    head(which(Quarter==start),1):tail(which(Quarter==end),1)
    ]

# Reset column order to start with Quarter and then lgd variables
lgd_vars = c('ss_lgd_4Q_Ma','sub_lgd_4Q_Ma')
setcolorder(
  corr_dt, 
  c(
    'Quarter',
    lgd_vars, 
    setdiff(names(corr_dt), c('Quarter', lgd_vars))
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
setnames(corr, c('LGD_Var','Macro_Var','Correlation','Sign'))
setorder(corr, 'LGD_Var')

# Rename variables to make it distinguish from similar vars
setnames(corr, c('LGD_Var','Sign'), c('LGD_4Q_MA','Sign.Corr'))
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
univ_reg_res_ss[, Sign.Consistent := Sign.X==Expected.Sign]
univ_reg_res_sub[, Sign.Consistent := Sign.X==Expected.Sign]

#================ Merge ADF results into univariate regression results ================#

# Merge regression results with adf test results
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
setkey(univ_reg_res_ss, NULL)
setkey(univ_reg_res_sub, NULL)
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

#=====================================================================================================
# MULTIVARIATE REGRESSION ANALYSIS
#=====================================================================================================

#================ Input parameters ================#

# Declare inputs
yname_ss = 'ss_lgd'
yname_sub = 'sub_lgd'
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
rolling.start.qtr='1996:Q1'
rolling.end.qtr=c('2010:Q4','2011:Q4','2013:Q1','2014:Q1','2015:Q1')
order = 2
index = c('Quarter', 'INSTRUMENT_ID')

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

# Declare downturn shift factor inputs
downturn_periods=c('2008:Q4','2009:Q1','2009:Q2','2009:Q3')
cycle_start_qtr='2008:Q1'
cycle_end_qtr_ss='2015:Q2'
cycle_end_qtr_sub='2015:Q2'

#================ Get candidate models ================#

source(paste0(code,'/Functions.R'))

# Call "do_multiv_regression" function to get all 2-factor candidate models for ss
best_ss2factors = 
  do_multiv_regression(
    ss, 
    macro_file=macro,
    macro_base,
    macro_sign_file, 
    yname=yname_ss, 
    xnames=xnames_ss,
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
    cycle_end_qtr_ss,
    order, 
    index
    )

# Call "do_multiv_regression" function to get all 3-factor models
best_ss3factors = 
  do_multiv_regression(
    ss,
    macro_file=macro, 
    macro_base,
    macro_sign_file, 
    yname=yname_ss,
    xnames=xnames_ss,
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
    cycle_end_qtr_ss,
    order,
    index
    )

# # Call "do_multiv_regression" function to get all 4-factor models
best_ss4factors = 
  do_multiv_regression(
    ss, 
    macro_file=macro, 
    macro_base,
    macro_sign_file,
    yname=yname_ss,
    xnames=xnames_ss,
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
    cycle_end_qtr_ss, 
    order, 
    index
    )

# Call "do_multiv_regression" function to get all 1-factor models for sub
best_sub1factors = 
  do_multiv_regression(
    sub, 
    macro_file=macro,
    macro_base, 
    macro_sign_file, 
    yname=yname_sub, 
    xnames=xnames_sub, 
    factor_subset=1, 
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
    cycle_end_qtr_sub,
    order, 
    index
    )

#================ Combine model results and select final candidate models ================#

# Combine ss model for 2 & 3 factor subsets
best_ss = rbind(best_ss2factors, best_ss3factors)
setorder(best_ss, -Adj.R2)
# setorder(best_ss2factors, -Adj.R2)
setorder(best_sub1factors, -Adj.R2)

# Subset ss models to exclude CBC and SP500
best_ss_final = best_ss[!grep('CBC', Formula)][!grep('SP500', Formula)]
# best_ss2factors_final = best_ss2factors[!grep('CBC', Formula)][!grep('SP500', Formula)]
best_sub_final = copy(best_sub1factors)