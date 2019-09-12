##################################################################################################
#
#                                                       FUNCTIONS
#
##################################################################################################

library(lubridate)
library(data.table)
library(stats)        # ks.test and wilson.test
library(zoo)
library(Hmisc)
library(fUnitRoots)   # adfTest function to test for stationarity
library(R.utils)      # For progress bar
library(car)          # load vif() function
library(boot)         # load cv.glm() function
library(lmtest)       # BP test and White test
library(sandwich)     # White sd correction
library(nortest)      # Normality test
library(plm)          # Autocorrelation test for panel regression


#=================================================================================================
# Generic functions
#=================================================================================================

# Calculate 4-Qtr moving average
calc_4qtr_ma = function(var_data){
  new = c()
  new[1] = var_data[1]
  new[2] = mean(var_data[1:2], na.rm = T)
  new[3] = mean(var_data[1:3], na.rm = T)
  for(i in 4:length(var_data)){
    new[i] = mean(var_data[(i-3):i], na.rm=T)
  }
  return(new)
}

# Calculate 4-Qtr moving sum
calc_4qtr_ms = function(var_data){
  new = c()
  new[1] = var_data[1]
  new[2] = sum(var_data[1:2], na.rm = T)
  new[3] = sum(var_data[1:3], na.rm = T)
  for(i in 4:length(var_data)){
    new[i] = sum(var_data[(i-3):i], na.rm=T)
  }
  return(new)
}


#====================================================================================================
# KS and MW tests
#====================================================================================================

# This function takes 2 data series as input, does KS test, and returns [stat, p]
ks_stat_p = function(dt1, dt2) {
  library(stats)
  ks = ks.test(dt1, dt2, exact=T)
  output = 
    c(
      round(
        as.numeric(ks$statistic),
        4
        ), 
      round(
        as.numeric(ks$p.value),
        4)
      )
  return(output)
}

# This function takes 2 data series as input and does MW tes, returns [stat, p]
mw_stat_p = function(dt1, dt2) {
  mw = wilcox.test(dt1, dt2)
  output = 
    c(
      round(
        as.numeric(mw$statistic),
        4
        ), 
      round(
        as.numeric(mw$p.value),
        4
        )
      )
  return(output)
}

#' Get KS and Mann-Whiteney test results for all unique groups under input grouping variable
#'
#' @author Gerard Mazi
#' @param group_var A character string, specifies the grouping variable/column
#' @param ur_data A dataset of at least two columns: 1) pred, and 2) the group_var
#' 
#' @return Return a list of four datasets: 
#' 1) ks_p_matrix, 
#' 2) ks_stat_matrix, 
#' 3) mw_p_matrix, 
#' 4) mw_stat_matrix
#' @export
#'
#' @examples
#' TBD
#'
get_ks_mw_matrix = function(group_var, ur_data) {
  
  options(warn = -1)    # Turn off warnings temporarily b/c of ties when computing ks
  
  data = as.data.table(copy(ur_data))
  groups = unique(data[[group_var]])    # Get unique "group" values
  l = length(groups)
  combo = expand.grid(groups, groups)    # Get pairwise combinations for all unique groups
  setnames(data, group_var, 'GROUP_VAR')
  
  # Do KS and MW tests for each pair of groups
  ks_results = NULL
  mw_results = NULL
  for(i in 1:nrow(combo)) {
    ks_results = 
      rbind(
        ks_results, 
        ks_stat_p(
          data[GROUP_VAR==combo[i,1], lgd], 
          data[GROUP_VAR==combo[i,2], lgd]
          )
        )
    mw_results = 
      rbind(
        mw_results,
        mw_stat_p(
          data[GROUP_VAR==combo[i,1], lgd], 
          data[GROUP_VAR==combo[i,2], lgd]
          )
        )
  }
  
  combo_other = as.data.table(cbind(groups, 'All_Others'))
  setnames(combo_other, names(combo))
  combo = rbind(combo, combo_other)
  
  # Do KS and MW tests for each of one group vs. all other groups
  ks_results_others = NULL
  mw_results_others = NULL
  for(j in 1:l) {
    ks_results_others = 
      rbind(
        ks_results_others,
        ks_stat_p(
          data[GROUP_VAR==groups[j], lgd],
          data[GROUP_VAR!=groups[j], lgd]
          )
        )
    mw_results_others = 
      rbind(
        mw_results_others,
        mw_stat_p(
          data[GROUP_VAR==groups[j], lgd], 
          data[GROUP_VAR!=groups[j], lgd]
          )
        )
  }
  
  # Aggregate all test results, including test-stats and p-values
  combo = 
    cbind(
      combo, 
      rbind(
        ks_results, 
        ks_results_others),
      rbind(
        mw_results, 
        mw_results_others
        )
      )
  setnames(
    combo, 
    c(
      paste0(group_var,'.x'), 
      paste0(group_var,'.y'),
      'KS.Stat',
      'KS.p.value',
      'MW.Stat',
      'MW.p.value'
      )
    )
  
  # Transform results into matrix format
  ks_stat_matrix = 
    as.data.table(
      matrix(
        combo$KS.Stat, 
        nrow = l+1, 
        ncol = l, 
        byrow = T
        )
      )
  ks_p_matrix = 
    as.data.table(
      matrix(
        combo$KS.p.value, 
        nrow = l+1, 
        ncol = l, 
        byrow = T
        )
      )
  mw_stat_matrix = 
    as.data.table(
      matrix(
        combo$MW.Stat, 
        nrow = l+1, 
        ncol = l, 
        byrow = T
        )
      )
  mw_p_matrix =
    as.data.table(
      matrix(
        combo$MW.p.value,
        nrow = l+1, 
        ncol = l, 
        byrow = T
        )
      )
  
  # Delete diagonal repetitive values
  for(k in 1:l-1) {
    ks_p_matrix[k, (k+1):l] = NA
    ks_stat_matrix[k, (k+1):l] = NA
    mw_p_matrix[k, (k+1):l] = NA
    mw_stat_matrix[k, (k+1):l] = NA
  }
  
  # Set names for the output matrices
  setnames(ks_stat_matrix, groups)
  setnames(ks_p_matrix, groups)
  setnames(mw_stat_matrix, groups)
  setnames(mw_p_matrix, groups)
  row.names(ks_stat_matrix) = c(groups, 'All Others')
  row.names(ks_p_matrix) = c(groups, 'All Others')
  row.names(mw_stat_matrix) = c(groups, 'All Others')
  row.names(mw_p_matrix) = c(groups, 'All Others')
  
  output = list(ks_p_matrix, ks_stat_matrix, mw_p_matrix, mw_stat_matrix)
  names(output) = c('ks.p','ks.stat','mw.p','mw.stat')
  
  options(warn = 0)
  return(output)
}

#=======================================================================================================
# Correlation Matrix
#=======================================================================================================

#' Get correlation matrix for all unique groups under the input grouping variable
#'
#' @author Gerard Mazi
#' @param group_var A character string, specifies the grouping variable/column
#' @param ur A dataset consisting at least two columns: 
#' 1) Quarter, 
#' 2) pred, and 
#' 3) the group_var
#' @param quarters A dataset with only the Quarter column
#' @param method One of "pearson", "kendall", or "spearman"; a character string
#' 
#' @return Return a list of 2 matrices: 
#' 1) correlation and 
#' 2) p-value
#' @export
#'
#' @examples
#' TBD
#'
get_corr_matrix = 
  function(
    group_var, 
    ur,
    quarters,
    method='pearson'
    ) {
  
  # Warning msg will show up when using 'kendall' and 'spearman' due to ties
  options(warn = -1)
  
  data = 
    as.data.table(
      copy(
        ur[,c('Quarter', 'lgd', group_var), 
           with=F]
        )
      )
  groups = unique(data[[group_var]])    # Get unique "group" values t
  l = length(groups)
  setnames(data, group_var, 'GROUP_VAR')
  
  #---------------- Get count-weighted 4Q MA for each segment in the group_var ----------------#
  
  cor.dt = as.data.table(copy(quarters))
  for (g in 1:l) {
    # Get moving average for current segment
    temp_dt = data[
      GROUP_VAR==groups[g],
      .(avg.lgd=mean(lgd, na.rm=T), .N),
      by=Quarter
      ]
    temp_dt = 
      merge(
        temp_dt, 
        quarters, 
        by='Quarter', 
        all=T
        )
    temp_dt[is.na(N), N:=0]
    temp_dt[, CW.4QMA := calc_4qtr_ms(avg.lgd*N)/calc_4qtr_ms(N)]
    cor.dt[, groups[g]] = temp_dt$CW.4QMA
    
    # Get moving average for all other segments combined
    temp_dt_others = 
      data[
        GROUP_VAR!=groups[g], 
        .(avg.lgd=mean(lgd, na.rm=T), .N),
        by=Quarter
        ]
    temp_dt_others = 
      merge(
        temp_dt_others,
        quarters, 
        by='Quarter',
        all=T
        )
    temp_dt_others[is.na(N), N:=0]
    temp_dt_others[, CW.4QMA := calc_4qtr_ms(avg.lgd*N)/calc_4qtr_ms(N)]
    cor.dt[, paste0(groups[g],'.Others')] = temp_dt_others$CW.4QMA
    
    rm(temp_dt, temp_dt_others)
  }
  
  #---------------- Calculate correlation ----------------#
  
  # Get pairwise combinations for all unique groups
  combo = expand.grid(groups, groups)
  setnames(
    combo, 
    c('Var1','Var2'), c('V1','V2')
    )    # Renaming is required here for the rbind step below
  combo = 
    rbind(
      combo, 
      data.table(
        c(groups), 
        paste0(groups, ".Others")
        )
      )
  combo = as.matrix(combo)
  
  # Calculate correlation for each pair of groups
  corr_results_r = c()
  corr_results_p = c()
  for(i in 1:nrow(combo)) {
    temp.x = as.numeric(unlist(cor.dt[,combo[i,1],with=F]))
    temp.y = as.numeric(unlist(cor.dt[,combo[i,2],with=F]))
    temp.corr.res = 
      cor.test(
        temp.x,
        temp.y, 
        ethod = method
        )
    corr_results_r = 
      c(
        corr_results_r, 
        as.numeric(temp.corr.res$estimate)
        )
    corr_results_p = 
      c(
        corr_results_p,
        as.numeric(temp.corr.res$p.value)
        )
  }
  
  #---------------- Format results into matrices ----------------#
  
  # Transform results into matrix format
  corr_r_matrix = 
    as.data.table(
      matrix(
        corr_results_r,
        nrow = l+1, 
        ncol = l, 
        byrow = T
        )
      )
  corr_p_matrix =
    as.data.table(
      matrix(
        corr_results_p,
        nrow = l+1,
        ncol = l,
        byrow = T
        )
      )
  
  # Delete diagonal repetitive values
  for(k in 1:l-1) {
    corr_r_matrix[k, (k+1):l] = NA
    corr_p_matrix[k, (k+1):l] = NA
  }
  
  # Set names for the output matrices
  setnames(corr_r_matrix, groups)
  setnames(corr_p_matrix, groups)
  row.names(corr_r_matrix) = c(groups, 'All Others')
  row.names(corr_p_matrix) = c(groups, 'All Others')
  
  output = list(corr_r_matrix, corr_p_matrix)
  names(output) = c('correlation','p.value')
  
  options(warn = 0)
  return(output)
}

#=======================================================================================================
# Macroeconomic Factor transformation
#=======================================================================================================

#' Do transformations for a single factor
#'
#' @author Gerard Mazi
#' @param factor A numeric factor which is to be transformed
#' @param is.percentage A logical flag for percentage formatting
#' 
#' @return Return a data table with all transformed columns; 
#' @export
#'
#' @examples
#' TBD
#'
univ_transform = function(factor, is.percentage) {
  
  library(Hmisc)
  
  #Apply QoQ level difference
  factor_QQ_D = factor - Lag(factor)
  
  #Apply QoQ percentage growth
  factor_QQ_G = (factor/Lag(factor)) - 1
  
  #Apply YoY level difference
  factor_YY_D = factor - Lag(factor, 4)
  
  #Apply YoY percentage growth
  factor_YY_G = (factor/Lag(factor, 4)) - 1
  
  #Apply 1Q lag to all of the transformed factors
  factor_L1 = Lag(factor)
  factor_QQ_D_L1 = Lag(factor_QQ_D)
  factor_QQ_G_L1 = Lag(factor_QQ_G)
  factor_YY_D_L1 = Lag(factor_YY_D)
  factor_YY_G_L1 = Lag(factor_YY_G)
  
  #Apply 2Q lag to all of the transformed factors
  factor_L2 = Lag(factor,2)
  factor_QQ_D_L2 = Lag(factor_QQ_D,2)
  factor_QQ_G_L2 = Lag(factor_QQ_G,2)
  factor_YY_D_L2 = Lag(factor_YY_D,2)
  factor_YY_G_L2 = Lag(factor_YY_G,2)
  
  #Apply 1Q lead to all of the transformed factors
  factor_LL1 = Lag(factor,-1)
  factor_QQ_D_LL1 = Lag(factor_QQ_D,-1)
  factor_QQ_G_LL1 = Lag(factor_QQ_G,-1)
  factor_YY_D_LL1 = Lag(factor_YY_D,-1)
  factor_YY_G_LL1 = Lag(factor_YY_G,-1)
  
  #Apply 2Q lead to all of the transformed factors
  factor_LL2 = Lag(factor,-2)
  factor_QQ_D_LL2 = Lag(factor_QQ_D,-2)
  factor_QQ_G_LL2 = Lag(factor_QQ_G,-2)
  factor_YY_D_LL2 = Lag(factor_YY_D,-2)
  factor_YY_G_LL2 = Lag(factor_YY_G,-2)
  
  # Combine all transformations
  output = 
    data.table(
      factor_QQ_D, factor_QQ_G, factor_YY_D, factor_YY_G, factor_L1, factor_QQ_D_L1,
      factor_QQ_G_L1, factor_YY_D_L1, factor_YY_G_L1,factor_L2, factor_QQ_D_L2, 
      factor_QQ_G_L2, factor_YY_D_L2, factor_YY_G_L2,factor_LL1, factor_QQ_D_LL1, 
      factor_QQ_G_LL1, factor_YY_D_LL1, factor_YY_G_LL1,factor_LL2, factor_QQ_D_LL2, 
      factor_QQ_G_LL2, factor_YY_D_LL2, factor_YY_G_LL2
      )
  
  if (is.percentage==1) { # Exclude growth if the factor is already in percentage
    output = output[, !grep('_G', names(output)), with=F]
  }
  if (is.percentage==0) { # Exclude level difference if the factor is not in percentage
    output = output[, !grep('_D', names(output)), with=F]
  }
  
  return(output)
}

#' Call univ_transform function and transform all macros
#'
#' @author Gerard Mazi
#' @param macro_dt A data frame/table consisting macros to be transformed
#' @param pct_flag_matrix A mapping table with two columns: 
#' 1) macro factor name,  
#' 2) corresponding "is.percentage" flag
#' 
#' @return Return a data frame that combines macro_dt and all transformations
#' @export
#'
#' @examples
#' TBD
#'
do_transformation = function (macro_dt, pct_flag_matrix) {
  
  # Make a copy of input data to avoid changing it after calling function
  transformed_macro = copy(as.data.frame(macro_dt))
  
  # Run a loop to apply function "univ_transform" to all factors of interest
  for (i in 1:nrow(pct_flag_matrix)) {
    factor_name = as.data.frame(pct_flag_matrix)[i, 1]
    percentage_flag = as.data.frame(pct_flag_matrix)[i, 2]
    new_transform = 
      univ_transform(
        macro[,factor_name], 
        percentage_flag
        )
    setnames(
      new_transform, 
      gsub(
        pattern = 'factor', 
        replacement = factor_name, 
        x = names(new_transform)
        )
      )
    
    # Append transformations of new base factor to the output transformed_macro dataset
    transformed_macro = 
      data.frame(
        transformed_macro, 
        new_transform
        )
  }
  
  return(transformed_macro)
}


#=========================================================================================================
# Univariate regression analysis
#=========================================================================================================

#' Display regression summary table
#'
#' @author Gerard Mazi
#' @param model A model object, e.g., output from lm or glm
#' @param dt Regression dataset which includes all x and y variables and Quarter 
#'
#' @return Return the summary table of a single univariate model; a data table
#' @export
#'
#' @examples
#' TBD
#'
return_univ_reg_summary = function(model) {
  
  options(warn = -1)
  
  # Get each piece of information needed from the regression model
  Y.Var = names(model$model)[1]
  X.Var = names(model$model)[-1]
  formula = 
    paste0(
      names(model$model)[1],
      ' ~ ', 
      paste0(
        names(model$model)[-1], 
        collapse='+'
        )
      )
  coef.intercept = as.numeric(model$coefficients[1])
  coef.x = as.numeric(model$coefficients[2])
  sign.x = sign(coef.x)
  t.stats = summary(model)$coefficients[2,3]   
  p.value = summary(model)$coefficients[2,4]   
  adj.r2 = summary(model)$adj.r.square
  rmse = sqrt(mean(residuals(model)^2))
  N.Total = nrow(model$model)
  
  # Combine all information
  output = 
    data.table(
      Y.Var, 
      X.Var, 
      coef.intercept, 
      coef.x, 
      sign.x, 
      adj.r2,
      t.stats,
      p.value, 
      rmse, 
      N.Total
      )
  setnames(
    output, 
    c(
      'Y.Var',
      'X.Var',
      'Coef.Intercept',
      'Coef.X',
      'Sign.X',
      'Adj.R2',
      't.stats',
      'p.value',
      'RMSE', 
      'N.Total'
      )
    )
  
  options(warn = 0)
  return(output)
}

#' Univariate regression analysis
#'
#' @author Gerard Mazi
#' @param dt Regression dataset which includes all x and y variables, Quarter
#' @param ynames Names of all depedent variables; a character vector
#' @param xnames Names of all independent variables that will be used to do regression
#' @param type Type of regression; must be one of "ols" and "logistic"
#'
#' @return Return summary table of all univariate regressions; a data table
#' @export
#'
#' @examples
#' TBD
#'
univ_reg = function(dt, ynames, xnames) {
  
  # Constuct univariate regression formulas
  yx_combo = 
    expand.grid(
      ynames, 
      xnames
      )
  formulas = 
    paste(
      yx_combo$Var1, 
      yx_combo$Var2, 
      sep = '~'
      )
  
  # Get a list of all regression models
  reg_res = 
    lapply(
      formulas, 
      lm, 
      data=dt
      )
  
  # Call "return_reg_summary" function above and outputs regression summary for all models
  reg_res_summary =
    rbindlist(
      lapply(
        reg_res,
        return_univ_reg_summary
        )
      )
  
  return(reg_res_summary)
}


#========================================================================================================
# Univariate factor selection
#========================================================================================================

#' Univariate factor selection
#'
#' @author Gerard Mazi
#' @param reg_res A summary table which is merged from: 
#' 1) univariate regression results (the output from "univ_reg"),
#' 2) correlation results (4 cols: Y, X, Cor, Sign.Cor), 
#' 3) expected correlation sign table (2 cols: X, Expected.Sign),
#' 4) ADF stationarity test
#' @param cor.thresh Threshold of correlation; double number
#' @param p.thresh Threshold of p-value; double number
#' @param adf.p.thresh Threshold of p-value in ADF stationarity test; double number
#'
#' @return Return summary table of selected univariate regressions; a data table
#' @export
#'
#' @examples
#' TBD
#'
# This functions select univariate factors based on a series of criterion
select_univ_factors = 
  function (
    reg_res, 
    cor.thresh=0.2, 
    p.thresh=0.1, 
    adf.p.thresh=0.05
    ) {
  
  # Make a copy of input dataset to avoid changing the original input data if it is data.table
  output = copy(as.data.table(reg_res))
  
  # 1. Filter out non-stationary factors
  output = 
    output[
      , 
      ADF.p := 
        ifelse(
          ADF.Type=='Zero.Mean', 
          ADF.nc.p, 
          ifelse(
            ADF.Type=='Single.Mean', 
            ADF.c.p, 
            ADF.ct.p
            )
          ), 
      by=1:nrow(output)
      ]
  
  output = output[ADF.p < adf.p.thresh]
  cat(
    paste0(
      nrow(output), 
      ' factors are left after excluding non-stationary factors\n\n')
    )
  output[, ADF.p := NULL]
  if (nrow(output)==0) return(output)
  
  # 2. Filter out factors that do not have same coefficient sign as expected
  output = output[Sign.Consistent==TRUE]
  cat(
    paste0(
      nrow(output), 
      ' factors are left after checking coefficient sign\n\n'
      )
    )
  if (nrow(output)==0) return(output)
  
  # 3. Filter out factors that have less than 20% correlation with Y
  output = output[abs(Correlation) > cor.thresh]
  cat(
    paste0(
      nrow(output), 
      ' factors are left after excluding correlation less than ', cor.thresh, '\n\n'
      )
    )
  if (nrow(output)==0) return(output)
  
  # 4. Filter out factors that have less than 0.1 p-value
  output = output[p.value < p.thresh]
  cat(
    paste0(
      nrow(output), 
      ' factors are left after excluding p-value greater than ', 
      p.thresh, 
      '\n\n'
      )
    )
  
  return(output)
}


#=========================================================================================================
# Multivariate regression analysis
#=========================================================================================================

#' Create all regression formulas that exclude factor transformations from the same macro base variable;
#' a much faster version compared to "create_formulas"
#'
#' @author Gerard Mazi
#' @param factor_subset The number of factors to be included in the model
#' @param max_factor_subset The max number of factors that is allowed
#' @param yname Name of single dependent variable
#' @param xnames Names of all independent variables; a character vector
#'
#' @return Return all constructed formulas; a vector of character
#' @export
#'
#' @examples
#' TBD
#'
create_formulas_v2 = 
  function(
    factor_subset, 
    max_factor_subset, 
    yname, xnames
    ) {
  
  # Get macro base for selected xnames
  xbases = 
    unique(
      gsub("^(.*?)_.*", "\\1", xnames)
      )
  
  # Put xnames into distinct groups based on their corresponding macro_base
  xname_list = list(length(xbases))
  xname_list_length = c()
  for (i in 1:length(xbases)) {
    xname_list[[i]] = 
      grep(
        xbases[i], 
        xnames, 
        value=T
        )
    xname_list_length[i] = length(xname_list[[i]])
  }
  
  # Generate combinations of groups which we will choose xnames from; 
  # Note the numeric values here indicates list element index in xname_list
  if (length(xbases) < factor_subset) {
    message(
      paste0(
        'Number of distinct xbases is less than specified factor subset; ',
        'factor subset was reset to the number of distinct xbases'
        )
      )
    factor_subset = length(xbases)
  }
  group_combo = combn(1:length(xbases), factor_subset)

  # This function gets formulas within each xnames groups, which is defined above
  get_ingroup_formulas = 
    function(
      group, 
      max_factor_subset, 
      factor_subset, 
      xname_list
      ) {
    formulas = list()
    formulas[[1]] = 
      paste0(
        yname, 
        ' ~ ', 
        xname_list[[group[1]]]
        )
    j=2
    while (j <= max_factor_subset) {
      if (factor_subset >= j) {
        temp_formulas = list()
        for (m in 1:length(xname_list[[group[j]]])) {
          temp_formulas[[m]] = 
            paste0(
              formulas[[j-1]], 
              ' + ', 
              xname_list[[group[j]]][m]
              )
        }
        formulas[[j]] = unlist(temp_formulas)
      }
      j=j+1
    }
    # Delete univariate formulas
    output = formulas[[factor_subset]]
    
    # Return formula
    return(output)
  }
  
  # Call above function to get all formulas under the specified factor subset
  formulas = 
    unique(
      unlist(
        apply(
          group_combo, 
          2, 
          get_ingroup_formulas,
          max_factor_subset, 
          factor_subset, 
          xname_list
          )
        )
      )
  cat(
    paste0(
      '\n', 
      length(formulas), 
      ' formulas were constructed for factor subset ', 
      factor_subset, '\n'
      )
    )
  return(formulas)
}


#' Create all regression formulas that exclude factor transformations from the same macro base variable;
#' a slower version compared to "create_formulas_v2"
#'
#' @author Gerard Mazi
#' @param factor_subset The number of factors to be included in the model
#' @param yname Name of single dependent variable
#' @param xnames Names of all independent variables; a character vector
#'
#' @return Return all constructed formulas; a vector of character
#' @export
#'
#' @examples
#' TBD
#'
# Create all regression formulas
create_formulas = 
  function(
    factor_subset, 
    yname, 
    xnames
    ) {
  
  temp_combo = 
    as.data.table(
      combn(
        xnames, 
        factor_subset
        )
      )    # All possible combination of x_var
  
  # Assign indexes to x_combos that contain 2 or more transformations from same macro base
  dup_macro_index = list()
  pb = 
    txtProgressBar(
      min = 0, 
      max = length(macro_base), 
      style = 3
      )
  for (i in 1:length(macro_base)) {
    dup_macro_index[[i]] = 
      c(
        dup_macro_index, 
        which(
          as.numeric(
            temp_combo[
              , 
              lapply(
                .SD, 
                function(x){
                  length(
                    grep(macro_base[i], x)
                    )
                  }
                )
              ]
            )>1)
        )
    setTxtProgressBar(pb, i)
  }
  close(pb)
  # Exclude above indexes and get x_combos with transformations from distinct macro_base
  x_combo = 
    as.data.table(t(temp_combo))[!unique(unlist(dup_macro_index))]
  yx_combo = 
    expand.grid(
      yname, 
      apply(
        x_combo[, .SD],
        1, 
        paste, 
        collapse=" + "
        )
      )
  formulas = 
    paste(
      yx_combo$Var1,
      yx_combo$Var2,
      sep=' ~ '
      )
  
  return(formulas)
}

# This function calculates the R-Squared in lm and glm model
R2gauss<- function(y,model){
  moy<-mean(y)
  N<- length(y)
  p<-length(model$coefficients)-1
  SSres<- sum((y-predict(model))^2)
  SStot<-sum((y-moy)^2)
  R2<-1-(SSres/SStot)
  Rajust<-1-(((1-R2)*(N-1))/(N-p-1))
  return(data.frame(R2,Rajust,SSres,SStot))
}

# This function conducts White test for heteroskedasticity
do_WhiteTest = function (model, data) {
  library(lmtest)
  
  yname = names(model$model)[1]    # Y variable name
  xnames = setdiff(names(model$model), yname)
  square_terms = paste0('I(', xnames, '^2)')
  
  if (length(xnames)>1) {
    xnames_combo = combn(xnames, 2)
    cross_product_terms = 
      apply(
        xnames_combo, 
        2,
        paste, 
        collapse=' * '
        )
    auxiliary_formula = 
      formula(
        paste0(
          yname, 
          ' ~ ', 
          paste(
            c(xnames, square_terms, cross_product_terms), 
            collapse = ' + '
            )
          )
        )
  } else {
    auxiliary_formula = 
      formula(
        paste0(
          yname, 
          ' ~ ', 
          paste(
            c(xnames, square_terms), 
            collapse = ' + '
            )
          )
        )
  }
  bp_res = 
    bptest(
      model,
      auxiliary_formula, 
      data = data
      )
  # white.p = bp_res$p.value
  return(bp_res)
}

#' Display multivariate regression result for a single model
#'
#' @author Gerard Mazi
#' @param model A model object obtained from lm, glm, etc
#' @param data Regression dataset
#' @param max_factor_subset The max number of factors that is allowed
#' @param macro_sign_file A dataset with two columns: 
#' 1) macro base var names 
#' 2) corresponding intuitive signs
#' @param hetero.p.thresh P-value threshold for heteroskedasticity; default is 5%
#'
#' @return Return regression result of a single model; a dataset with only one row
#'
#' @examples
#' 
return_multiv_reg_summary = 
  function(
    model,
    data, 
    max_factor_subset=4, 
    macro_sign_file, 
    hetero.p.thresh=0.05
    ) {
  
  library(car)    # load vif() function
  library(data.table)
  library(lmtest)    # BP test
  library(sandwich)    # White sd correction
  
  #---------------- Outputs that can be obtained by one step ----------------#
  
  n.xvar = ncol(model$model)-1    # Number of X variables
  Y.Var = names(model$model)[1]    # Y variable name
  xnames = setdiff(names(model$model), Y.Var)
  X.Var = c(
    xnames, 
    rep(NA, max_factor_subset-ncol(model$model)+1)
    )
  formula = 
    paste0(
      names(model$model)[1], 
      ' ~ ', 
      paste0(
        names(model$model)[-1], 
        collapse=' + '
        )
      )    # Regression formula
  coef.intcpt = as.numeric(model$coefficients[1])    # Intercept coefficient
  t.intcpt = summary(model)$coefficients[1,3]    # Intercept t-stats
  p.intcpt = summary(model)$coefficients[1,4]    # Intercept p-value
  adj.r2 = as.numeric(R2gauss(model$model[,Y.Var], model)[2])
  rmse = sqrt(mean(residuals(model)^2))
  max.vif = ifelse(length(xnames)>1, max(vif(model)), 1)    # VIF
  model.aic = AIC(model)    # AIC
  model.bic = BIC(model)    # BIC
  N.Total = nrow(model$model)    # Total number of observations used to run regression
  rmse.cv = as.numeric(NA)
  
  # White test for heteroskedasticity
  white.res = do_WhiteTest(model, data)
  white.test.p = white.res$p.value
  white.corrected.stats = 
    coeftest(
      model, 
      vcov=vcovHC(model, type='HC3')
      )
  # white.coef.intcpt = white.corrected.stats[1,1]
  white.p.intcpt = white.corrected.stats[1,4]
  hetero.flag = as.numeric(white.test.p<hetero.p.thresh)
  
  #---------------- Declare outputs that need multiple steps to get  ----------------#
  
  coef.x = as.numeric(rep(NA, max_factor_subset))    # Coefficients of X variables
  white.coef.x = as.numeric(rep(NA, max_factor_subset))
  sign.as.expected = as.numeric(rep(NA, max_factor_subset))    # Sign-as-expected flag
  t.stats = as.numeric(rep(NA, max_factor_subset))    # t-statistics
  p.value = as.numeric(rep(NA, max_factor_subset))    # p-value
  white.p = as.numeric(rep(NA, max_factor_subset))    # p-value 

  #---------------- Run loop to get the above outputs  ----------------#
  # Note placeholders are created for up to max_factor_subset x vars; if less x vars, output NA
  
  for (i in 1:n.xvar) { # Loop starts from 2 b/c 1 is Y
    coef.x[i] = as.numeric(model$coefficients[i+1])
    sign.x = sign(coef.x[i])    # a temp var to store sign of coef.x in current loop
    x.base = gsub("^(.*?)_.*", "\\1", X.Var[i])    # a temp var to store base of x in current loop
    sign.expected = macro_sign_file[Base==x.base, Expected.Sign]    # a temp var to store expected sign of x
    sign.as.expected[i] = as.numeric(sign.x==sign.expected)
    t.stats[i] = summary(model)$coefficients[i+1,3]   # row# is the row of x var, and col# is where t-stats is
    p.value[i] = summary(model)$coefficients[i+1,4]   # row# is the row of x var, and col# is where p-value is
    white.p[i] = white.corrected.stats[i+1,4]
  }
  
  #---------------- Construct output dataset  ----------------#
  
  output = 
    data.frame(
      formula, 
      Y.Var, 
      t(X.Var), 
      coef.intcpt, 
      t(coef.x), 
      t(sign.as.expected),
      t.intcpt, 
      t(t.stats), 
      p.intcpt, 
      t(p.value),
      adj.r2, 
      rmse, 
      rmse.cv, 
      max.vif,
      model.aic, 
      model.bic, 
      white.test.p,
      hetero.flag, 
      white.p.intcpt, 
      t(white.p), 
      N.Total, 
      stringsAsFactors = F
      )
  setnames(
    output, 
    c(
      'Formula', 
      'Y', 
      paste0('X',1:max_factor_subset),
      'Coef.Intcpt', 
      paste0('Coef.X',1:max_factor_subset),
      paste0('Sign.Match.X',1:max_factor_subset), 
      't.Intcpt', 
      paste0('t.X',1:max_factor_subset),
      'p.Intcpt', 
      paste0('p.X',1:max_factor_subset),
      'Adj.R2', 
      'RMSE', 
      'RMSE.CV', 
      'VIF', 
      'AIC', 
      'BIC', 
      'Hetero.Test.p', 
      'Hetero.Flag', 
      'White.p.Intcpt', 
      paste0('White.p.X',1:max_factor_subset), 
      'N.Total'
      )
    )
  return(output)
}

#' Display multivariate regression results for all input formulas
#'
#' @author Gerard Mazi
#' @param factor_subset
#' @param max_factor_subset The max number of factors that is allowed
#' @param data Regression dataset
#' @param formulas All regression formulas;
#' @param macro_base A character vector of macro base variable names
#' @param macro_sign_file A dataset with two columns: 
#' 1) macro base var names 
#' 2) corresponding intuitive signs
#' @param hetero.p.thresh P-value threshold for heteroskedasticity; default is 5%
#'
#' @return Return a dataset with all regression results; each row corresponds to a model/formula
#'
#' @examples
#' 
multiv_reg = 
  function(
    factor_subset,
    max_factor_subset=4, 
    data, 
    formulas,
    macro_base, 
    macro_sign_file, 
    hetero.p.thresh=0.05
    ) {
  
  library(R.utils)    # For progress bar
  library(data.table)
  library(car)    # load vif() function
  options(warn = -1)
  
  reg_res = list()
  cat(
    paste0(
      '\nThere are ',
      length(formulas),
      ' formulas in factor subset ',
      factor_subset,' run\n'
      )
    )
  pb = 
    txtProgressBar(
      min = 0, 
      max = length(formulas), 
      initial=NA,
      style = 3
      )    # Progress Bar
  for (i in 1:length(formulas)) {
    fitted.model = 
      glm(
        formulas[i], 
        data, 
        family = 'gaussian'
        )
    reg_res[[i]] = 
      return_multiv_reg_summary(
        fitted.model, 
        data, 
        max_factor_subset, 
        macro_sign_file, 
        hetero.p.thresh
        )
    setTxtProgressBar(pb, i)
  }
  close(pb)
  reg_res_summary = rbindlist(reg_res)
  
  options(warn=0)
  return(reg_res_summary)
}


#=========================================================================================================
# Multivaraite factor selection
#=========================================================================================================

#' Multivariate factor selection based on 
#' 1) intuitive signs, 
#' 2) heteroskedasticity and (White corrected) p-value, 
#' 3) VIF
#'
#' @author Gerard Mazi
#' @param reg_res A dataset with regression results of all models; output from multiv_reg()
#' @param vif.thresh VIF threshold for model selection; default is 5
#' @param p.thresh P-value threshold for model selection; default is 5%
#' @param hetero.p.thresh P-value threshold to detect heteroskedasticity; default is 5%
#' 
#' @return Return a dataset of survived models with their summary stats
#' @export
#'
#' @examples
#' TBD
#'
select_multiv_factors = 
  function (
    reg_res, 
    vif.thresh=5, 
    p.thresh=0.05,
    hetero.p.thresh=0.05
    ) {
  
  # Make a copy of input dataset to avoid changing the original input data if it is data.table
  output = copy(as.data.table(reg_res))
  cat(
    paste0(
      '\nThere are ',
      nrow(output),
      ' models before best subset starts\n\n'
      )
    )
  
  #---------------- 1. Check intuition of coefficient sign ----------------#
  
  sign.check = 
    output[
      ,c(grep('Sign.Match.X', names(output))), 
      with=F
      ]   # Subset data to Sign-Match columns
  keep_index = c()
  cat('Checking coefficient signs...\n')
  pb = 
    txtProgressBar(
      min = 0, 
      max = nrow(sign.check), 
      initial=NA, 
      style = 3
      )    # Progress Bar
  for (i in 1:nrow(sign.check)) {
    if (all(sign.check[i]==1, na.rm=T)) keep_index[i] = i
    else keep_index[i] = NA
    setTxtProgressBar(pb,i)
  }
  close(pb)
  keep_index = keep_index[!is.na(keep_index)]
  output = output[keep_index]
  rm(sign.check)
  rm(keep_index)
  cat(paste0(nrow(output), ' models are left after excluding unintuitive signs\n\n'))
  if (nrow(output)==0) return(output)
  
  #---------------- 2. Check heteroskedasticity and p-value ----------------#
  
  p.check = 
    output[
      , 
      c('p.Intcpt', grep('^p.X', names(output), value=T)),
      with=F
      ]
  White.p.check = 
    output[
      , 
      c('White.p.Intcpt', grep('White.p.X', names(output), value=T)),
      with=F
      ]
  keep_index = c()
  cat('Checking heteroskedasticity and p-value...\n')
  pb = 
    txtProgressBar(
      min = 0, 
      max = nrow(p.check), 
      initial=NA, 
      style = 3
      )    # Progress Bar
  for (i in 1:nrow(p.check)) {
    hetero.p = output[i, Hetero.Test.p]
    if (hetero.p>hetero.p.thresh & all(p.check[i]<p.thresh, na.rm=T)) {
      keep_index[i] = i
    } else if (hetero.p<=hetero.p.thresh & all(White.p.check[i]<p.thresh, na.rm=T)) {
      keep_index[i] = i
    } else keep_index[i] = NA
    setTxtProgressBar(pb,i)
  }
  close(pb)
  keep_index = keep_index[!is.na(keep_index)]
  output = output[keep_index]
  rm(p.check)
  rm(White.p.check)
  rm(keep_index)
  cat(
    paste0(
      nrow(output), 
      ' models are left after excluding (White corrected) p-value greater than ', 
      p.thresh, 
      '\n\n'
      )
    )
  if (nrow(output)==0) return(output)
  
  #---------------- 3. Check VIF ----------------#
  
  cat('Checking VIF...\n')
  output = output[VIF<=vif.thresh]
  cat(
    paste0(
      nrow(output), 
      ' models are left after excluding VIF higher than ',
      vif.thresh, 
      '\n\n'
      )
    )

  return(output)
}

#' Add a new column named "Pass.Rolling.Test" to the survived regression result dataset, indicating
#' whether the model passed the forward chain rolling back test; exclude models that fail the rolling test
#'
#' @author Gerard Mazi
#' @param selected_models A dataset with all survived regression model results
#' @param data Same regression dataset as what was used for modeling
#' @param max_factor_subset The max number of factors that is allowed/
#' @param hetero.p.thresh p-value threshold to detect heteroskedasticity; default is 5%
#' @param p.thresh p-value threshold for model selection; default is 5%
#' @param start_qtr Starting quarter of the rolling test period; default is "1996:Q1"
#' @param end_qtr A vector of ending quarters of the rolling test period; default is 
#' c('2010:Q4','2011:Q4','2012:Q4','2014:Q1','2015:Q1')
#' 
#' @return Return exactly the same dataset as the input "selected_models", 
#' @export
#'
#' @examples
#' TBD
#'
forward_chain_rolling = 
  function (
    selected_models, 
    data, 
    max_factor_subset=4, 
    hetero.p.thresh=0.05, 
    p.thresh=0.05,
    start_qtr='1996:Q1', 
    end_qtr=c('2010:Q4','2011:Q4','2012:Q4','2014:Q1','2015:Q1')
    ) {
  # Check if input calculation periods are reasonable
  if (!start_qtr %in% data$Quarter) stop('Start Quarter is not an element in the Quarter list of data!')
  if (!all(end_qtr %in% data$Quarter)) stop('Not all End Quarters are in the Quarter list of data!')
  models = copy(as.data.table(selected_models))
  setorder(data, Quarter)    # Order data by quarter
  
  cat('Doing forward-chain-rolling-back test...\n')
  pb = 
    txtProgressBar(
      min=0, 
      max=nrow(models), 
      initial=NA, 
      style=3
      )
  for (i in 1:nrow(models)) {
    # Get coefficient signs from original models
    models.coef = as.data.table(models)[i, c(grep('^Coef.X', names(models))), with=F]
    models.coef.signs = models.coef[, lapply(.SD, sign)]
    models.coef.signs = as.numeric(models.coef.signs)[!is.na((models.coef.signs))]
    
    # Get regression formula from regression result
    formula = formula(models$Formula[i])
    
    # Get "Pass.Rolling.Test" flag and add it to regression result
    for (qtr in end_qtr) {
      train_dt = as.data.table(data)[(which(Quarter==start_qtr)[1]):(last(which(Quarter==qtr)))]
      loop_model = lm(formula, train_dt)
      # print(summary(loop_model))
      loop_signs = as.numeric(sign(loop_model$coefficients[-1]))
      loop_p_values = as.numeric(summary(loop_model)$coefficients[,4])
      if (!all(loop_signs==models.coef.signs) | !all(loop_p_values<p.thresh)) {
        pass.rolling.test = 0
        break
      } else pass.rolling.test = 1
    }
    models[i, Passed.Rolling.Test := pass.rolling.test]
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  # Return output
  output = copy(models)
  output = output[Passed.Rolling.Test==1]    # Keep only models that passed the rolling back test
  cat(paste0(nrow(output), ' models passed the forward-chain-rolling-back test\n\n'))
  rm(models)
  return(output)
}

#' Conduct K-fold validation to survived regression models and calculate the cross-validation RMSE
#'
#' @author Gerard Mazi
#' @param selected_models A dataset with all survived regression model results
#' @param data Same regression dataset as what was used for modeling
#' @param K-fold Number of folds to be done in cross validation; default is 5
#' 
#' @return Add "RMSE.CV" added
#' @export
#'
#' @examples
#' TBD
#'
calculate_rmse.cv = function (selected_models, data, K_fold=5) {
  
  library(boot)    # load cv.glm() function
  
  models = copy(as.data.table(selected_models))
  
  # Do cross validation
  cat('Doing Cross-Validation...\n')
  pb = 
    txtProgressBar(
      min=0, 
      max=nrow(models), 
      initial=NA, 
      style=3
      )
  for (i in 1:nrow(models)) {
    formula = models$Formula[i]
    model.glm = 
      glm(
        formula, 
        data=data, 
        family=gaussian()
        )
    cv.res = 
      cv.glm(
        data=data, 
        glmfit=model.glm, 
        K=K_fold
        )
    rmse.cv = cv.res$delta[1]
    models[i, RMSE.CV:=rmse.cv]
    setTxtProgressBar(pb,i)
  }
  close(pb)
  cat('RMSE.CV values are populated in the regression output\n\n')
  
  output = copy(models)
  rm(models)
  return (output)
}


#' Add two new columns to the input regression results: 
#' 1) Unweighted.Shift.Factor 
#' 2) CW.Shift.Factor;
#' Note that col 
#' 1) is obtained by directly calculating predicted  using the macro data in each quarter, 
#' whereas col 
#' 2) is a count-weighted version
#'
#' @author Gerard Mazi
#' @param selected_models A dataset with all survived regression model results
#' @param data Same regression dataset as what was used for modeling
#' @param macro_file A dataset with macro data for each quarter, which is used to predict 
#' @param count_weight_file A dataset with two columns: 1) Quarter and 2) Count
#' @param downturn_periods A character vector specifying selected downturn periods
#' @param denominator_start_qtr Starting quarter of entire calculation period; default is '2008:Q1'
#' @param denominator_end_qtr Ending quarter of entire calculation period; default is '2016:Q2'
#' 
#' @return Return "Downturn.Shift.Factor" and "Downturn.Shift.Factor.CW" 
#' @export
#'
#' @examples
#' TBD
#'
calculate_downturn_shift_factor = 
  function (
    selected_models, 
    data, macro_file, 
    count_weight_file,
    downturn_periods=c('2008:Q4','2009:Q1','2009:Q2','2009:Q3'),
    denominator_start_qtr='2008:Q1', denominator_end_qtr='2015:Q2') {
  # Check if input calculation periods are reasonable
  if (!all(downturn_periods %in% data$Quarter)) stop('Not all downturn periods are in the Quarter list of data!')
  if (!denominator_start_qtr %in% data$Quarter) stop('Start Quarter in denominator is not in the Quarter list of data!')
  if (!denominator_end_qtr %in% data$Quarter) stop('End Quarter in denominator is not in the Quarter list of data!')
  
  models = copy(as.data.table(selected_models))
  overall_periods = 
    macro[
      which(Quarter==denominator_start_qtr):which(Quarter==denominator_end_qtr), Quarter
      ]
  
  cat('Calculating Downturn Shift Factor...\n')
  pb = 
    txtProgressBar(
      min=0, 
      max=nrow(models), 
      initial=NA, 
      style=3
      )
  for (i in 1:nrow(models)) {
    formula = models$Formula[i]
    model.glm = 
      glm(
        formula,
        data=data, 
        family=gaussian()
        )
    
    # Get predicted 
    fitted.y = 
      data.table(
        Quarter=
          macro_file[Quarter %in% overall_periods, Quarter],
        Fitted.y = 
          predict.glm(
            model.glm,
            newdata=macro_file[Quarter %in% overall_periods]
            )
        )
    fitted.y[Fitted.y<0] = 0
    fitted.y[Fitted.y>1] = 1
    downturn_fitted.y = fitted.y[Quarter %in% downturn_periods, Fitted.y]
    overall_fitted.y = fitted.y[Quarter %in% overall_periods, Fitted.y]
    
    # Get unweighted downturn shift factor by using directly predicted quarterly pred
    downturn_lgd = mean(downturn_fitted.y, na.rm=T)
    overall_lgd = mean(overall_fitted.y, na.rm=T)
    unweighted_shift_factor = downturn_lgd/overall_lgd - 1
    models[i, Unweighted.Shift.Factor:=unweighted_shift_factor]
    
    # Get count-weighted downturn shift factor
    downturn_counts = count_weight_file[Quarter %in% downturn_periods, Count]
    overall_counts = count_weight_file[Quarter %in% overall_periods, Count]
    downturn_cw_lgd = sum(downturn_fitted.y*downturn_counts, na.rm=T)/sum(downturn_counts, na.rm=T)
    overall_cw_lgd = sum(overall_fitted.y*overall_counts, na.rm=T)/sum(overall_counts, na.rm=T)
    count_weighted_shift_factor = downturn_cw_lgd/overall_cw_lgd - 1
    models[i, CW.Shift.Factor:=count_weighted_shift_factor]
    
    setTxtProgressBar(pb,i)
  }
  close(pb)
  cat('\n')
  
  output = copy(models)
  rm(models)
  return (output)
}

#' Conduct normality tests to the input regression results
#'
#' @author Gerard Mazi
#' @param selected_models A dataset with all survived regression model results
#' @param data Same regression dataset as what was used for modeling
#' 
#' @return Return 4 new columns
#' added: 
#' 1) SW.Normality.p, 
#' 2) KS.Normality.p, 
#' 3) CVM.Normality.p, 
#' 4) AD.Normality.p
#' @export
#'
#' @examples
#' TBD
#'
# This function does normality test to the input data series
do_normality_test = function (selected_models, data) {
  
  library(nortest)
  options(warn = -1)    # Suppress warnings b/c there can easily be ties in KS test
  models = copy(as.data.table(selected_models))
  
  cat('Conducting normality test...\n')
  pb = 
    txtProgressBar(
      min=0, 
      max=nrow(models),
      initial=NA, 
      style=3
      )
  for (i in 1:nrow(models)) {
    formula = models$Formula[i]
    model.glm = glm(
      formula, 
      data=data,
      family=gaussian()
      )
    epsilon = as.numeric(model.glm$residuals)
    
    # Do normality test
    sw.p = shapiro.test(epsilon)$p.value    # Shapiro-Wilk test
    ks.p = ks.test(epsilon, 'pnorm')$p.value    # Kolmogorov-Smirnov test
    cvm.p = cvm.test(epsilon)$p.value    # Cramer-von Mises test
    ad.p = ad.test(epsilon)$p.value    # Anderson-Darling test
    
    models[i, c('SW.Normality.p', 'KS.Normality.p', 'CVM.Normality.p', 'AD.Normality.p')
           := .(sw.p, ks.p, cvm.p, ad.p)]
    setTxtProgressBar(pb,i)
  }
  close(pb)
  cat('\n')
  
  output = copy(models)
  rm(models)
  options(warn = 0)
  return (output)
}

#' Conduct autocorrelation tests to model residuals
#'
#' @author Gerard Mazi
#' @param selected_models A dataset with all survived regression model results
#' @param data Same regression dataset as what was used for modeling
#' @param order Maximal order of serial correlation to be tested; an integer
#' @param max_factor_subset The max number of factors that is allowed/
#' @param index Input "time" and "group" variables
#' 
#' @return Return 6 new columns added:
#'         1) Autocorrelation.test.p (small p indicates there exists autocorrelation), 
#'         2) NW.p.intcpt (Newey-West corrected p-value for intercept; similar concept as White p), 
#'         3) NW.p.1 (Newey-West corrected p-value for the 1st x variable), 
#'         4) NW.p.2, 
#'         5) NW.p.3, and 
#'         6) NW.p.4
#' @export
#'
#' @examples
#' TBD
#'
do_autocorrelation_test = 
  function (
    selected_models, 
    data, order=2,
    max_factor_subset=4, 
    index=c('Quarter','INSTRUMENT_ID')) {
  
  library(plm)
  options(warn = -1)    # Suppress warnings b/c there can easily be ties in KS test
  models = copy(as.data.table(selected_models))
  
  cat('Conducting autocorrelation test...\n')
  pb = 
    txtProgressBar(
      min=0,
      max=nrow(models), 
      initial=NA, 
      style=3
      )
  
  NW.p.list = list()
  for (i in 1:nrow(models)) {
    # Get autocorrelation test p-value
    formula = formula(models$Formula[i])
    model = plm(formula, data=data, model='pooling', index=index)
    bg.res = pbgtest(model, order=order)
    bg.p = bg.res$p.value
    models[i, Autocorrelation.test.p := bg.p]
    
    # Use Newey-West standard error to correct p-value for autocorrelation (simliar idea as White corrected p)
    NW.corrected.stats = coeftest(model, vcov = vcovNW(model, type='HC3'))
    models[i, NW.p.intcpt:=NW.corrected.stats[1,4]]
    
    NW.p = rep(NA, max_factor_subset)
    for (j in 1:(nrow(NW.corrected.stats)-1)) {
      NW.p[j] = NW.corrected.stats[j+1,4]
    }
    NW.p = data.table(t(NW.p))
    setnames(NW.p, paste0('NW.p.', 1:max_factor_subset))
    NW.p.list[[i]] = NW.p
    
    setTxtProgressBar(pb,i)
  }
  close(pb)
  cat('\n')
  
  NW.p = rbindlist(NW.p.list)
  output = cbind(models, NW.p)
  rm(models)
  
  options(warn = 0)
  return (output)
}

#' Run regressions and do model selection
#'
#' @author Gerard Mazi
#' 
#' @return Return a dataset with models and corresponding key statistics
#' @export
#'
#' @examples
#' TBD
#'
do_multiv_regression = 
  function (
    data, 
    macro_file, 
    macro_base, 
    macro_sign_file, 
    yname, 
    xnames, 
    factor_subset,
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
    order=2, 
    index=c('Quarter','INSTRUMENT_ID')
    ) {
  
  #-------------------------------- Run regressions and get results --------------------------------#
  
  message('\nRunning regressions...')
  formulas = 
    create_formulas_v2(
      factor_subset, 
      max_factor_subset,
      yname, 
      xnames
      )
  multiv_reg_res = 
    multiv_reg(
      factor_subset,
      max_factor_subset, 
      data, 
      formulas, 
      macro_base,
      macro_sign_file,
      hetero.p.thresh
      )
  
  #-------------------------------- Model selection --------------------------------#
  
  message('\n\nDoing model selection...')
  
  # Step 1: Check heteroskedasticity, p-value and VIF
  best_subset_v1 = 
    select_multiv_factors(
      multiv_reg_res, 
      vif.thresh, 
      p.thresh, 
      hetero.p.thresh
      )
  if(nrow(best_subset_v1)==0) {
    message('\nNo model survives in the best subset model selection process\n\n')
    return(best_subset_v1)
  }
  
  # Step 2: Do forward-chain-rolling-back test
  best_subset_v2 = 
    forward_chain_rolling(
      best_subset_v1, 
      data, 
      max_factor_subset,
      hetero.p.thresh, 
      p.thresh,
      rolling.start.qtr, 
      rolling.end.qtr
      )
  if(nrow(best_subset_v2)==0) {
    message('\nNo model survives in the best subset model selection process\n\n')
    return(best_subset_v2)
  }
  
  # Step 3: Do cross-validation
  best_subset_v3 = 
    calculate_rmse.cv(
      best_subset_v2, 
      data, 
      K_fold
      )
  if(nrow(best_subset_v3)==0) {
    message('\nNo model survives in the best subset model selection process\n\n')
    return(best_subset_v3)
  }
  
  # Step 4: Calculate downturn shift factors
  best_subset_v4 = 
    calculate_downturn_shift_factor(
      best_subset_v3, 
      data,
      macro_file,
      count_weight_file,
      downturn_periods, 
      cycle_start_qtr, 
      cycle_end_qtr
      )
  if(nrow(best_subset_v4)==0) {
    message('\nNo model survives in the best subset model selection process\n\n')
    return(best_subset_v4)
  }
  
  # Step 5: Do normality test
  best_subset_v5 = 
    do_normality_test(
      best_subset_v4,
      data
      )
  if(nrow(best_subset_v5)==0) {
    message('\nNo model survives in the best subset model selection process\n\n')
    return(best_subset_v5)
  }
  
  # Step 6: Do autocorrelation test
  best_subset_v6 = 
    do_autocorrelation_test(
      best_subset_v5, 
      data,
      order, 
      max_factor_subset, 
      index=index
      )
  if(nrow(best_subset_v6)==0) {
    message('\nNo model survives in the best subset model selection process\n\n')
    return(best_subset_v6)
  }
  
  #-------------------------------- Rank candidate models --------------------------------#
  
  setorder(best_subset_v6, -Adj.R2)    # Rank by adj.r2
  message(
    paste0(
      '\nThere are a total of ',
      nrow(best_subset_v6),
      ' candidate models left after model selection\n\n'
      )
    )
  
  return(best_subset_v6)
}