#############################################################################################################
#
#                                                 MACRO DATA PROCESSING
#
##############################################################################################################

wd = '.....'
code = paste0(wd, '/CODE')
input = paste0(wd, '/INPUT')
output = paste0(wd, '/OUTPUT')

source(paste0(code,'/Functions.R'))

#==============================================================================================================
# MACRO FACTOR TRANSFORMATION
#==============================================================================================================

#==================== Declare parameters ====================#

# Factors of interest
macro_base = c("BKB","CBC","GDPR","IP","LBR","SP500","BAA7YSPREAD","GDPN")
save(macro_base, file=paste0(output,'/macro_base.Rdata'))   

#==================== Import and clean macro data ====================#

macro = 
  read.csv(
    file = paste0(input,'/ST_ECON_DATA_QTR.csv'), 
    stringsAsFactors = F
    )
names(macro) = gsub('_US','',names(macro))    # Simplify names)

# Assign quarter from As_of_Date
macro$Quarter = as.yearqtr(macro$As_of_Date, "%m/%d/%Y")
macro$Quarter = 
  paste(
    substr(macro$Quarter, 1, 4), ':', substr(macro$Quarter, 6, 7), 
    sep=""
    )

# Keep only actual data and remove forecasted
macro = macro[macro$Scenario == "Actual",]

# Keep only factors of interest as declared above, plus Quarter
macro = macro[, c('Quarter', macro_base)]

# Format as numeric for further data processing;
# Note that warning message will pop up due to missing values; it's OK to ignore
macro[, macro_base] = lapply(macro[,macro_base], as.numeric)

#==================== Do transformations ====================#

# Get the percentage flag matrix, which is an input for function "do_transformation"
var_tracker = 
  read.csv(
    paste0(input,'/con_Var_Tracker.csv'), 
    stringsAsFactors = F
    )
pct_flag_matrix = var_tracker[var_tracker$Base %in% macro_base, c('Base','Percent')]
rm(var_tracker)

# Call function "do_transformation"
transformed_macro = do_transformation(macro, pct_flag_matrix)

# Save transformed macro data
write.csv(
  transformed_macro, 
  paste0(output,'/macro_out.csv'), 
  row.names = F, 
  na=''
  )


#=======================================================================================================
# TIME SERIES ANALYSIS
#=======================================================================================================

ur = read.csv(
  paste0(input,'/output_ur.csv'), 
  stringsAsFactors = F
  )

start = '1996:Q1'
end = '2016:Q2'

#==================== Segment data and create time series respectively ====================#

# Create Segments
ss = ur[ur$BROAD_SENIORITY_COLLATERAL=="Senior Secured",]
sub = ur[ur$BROAD_SENIORITY_COLLATERAL=="Junior/Unsecured",]

# Sort by Quarter
setorder(ss, Quarter)
setorder(sub, Quarter)

# Create time series of LGD
ss_series = aggregate(x = ss$lgd, by=list(ss$Quarter), mean)
colnames(ss_series) = c("Quarter", "lgd_ss")
sub_series = aggregate(x = sub$lgd, by=list(sub$Quarter), mean)
colnames(sub_series) = c("Quarter", "lgd_sub")

# Add LGD observation count for ss, sub
ss_count = aggregate(x = ss$lgd, by=list(ss$Quarter), length)
colnames(ss_count) = c("Quarter", "count_ss")
sub_count = aggregate(x = sub$lgd, by=list(sub$Quarter), length)
colnames(sub_count) = c("Quarter", "count_sub")

# Add missing quarters in order to have a complete time series
quarters = 
  read.csv(
    file = paste0(input,'/quarters.csv'), 
    stringsAsFactors = F, 
    col.names = 'Quarter'
    )
ss_series = 
  merge(
    x = quarters, 
    y = ss_series, 
    by = "Quarter", 
    all.x = T, 
    all.y = F
    )
ss_series = 
  merge(
    x = ss_series, 
    y = ss_count, 
    by = "Quarter", 
    all.x = T, 
    all.y = F
    )
sub_series = 
  merge(
    x = quarters,
    y = sub_series, 
    by = "Quarter", 
    all.x = T, 
    all.y = F
    )
sub_series =
  merge(
    x = sub_series, 
    y = sub_count, 
    by = "Quarter",
    all.x = T, 
    all.y = F
    )

#================ Merge time series together and calculate count-weighted MA ====================#

# Merge ss and sub in one dataset for futher analysis
ss_sub = 
  merge(
    ss_series, 
    sub_series, 
    by= "Quarter"
    )

# Merge with macro data to plot. Plots done in Excel are based on this data
ss_sub_macro = 
  merge(
    ss_sub, 
    macro, 
    by="Quarter", 
    all.x = T, 
    all.y = F
    )

# Calculate 4 quarter count weighted moving average LGD
ss_sub_macro$ss_lgd_4Q_Ma = 
  calc_4qtr_ms(ss_sub_macro[, "lgd_ss"] * ss_sub_macro[,"count_ss"])/
  calc_4qtr_ms(ss_sub_macro$count_ss)
ss_sub_macro$sub_lgd_4Q_Ma = 
  calc_4qtr_ms(ss_sub_macro[, "lgd_sub"] * ss_sub_macro[,"count_sub"])/
  calc_4qtr_ms(ss_sub_macro$count_sub)

# Keep only the non-transformed factors to perform time series plots
ss_sub_macro = 
  ss_sub_macro[,
               c("Quarter", "ss_lgd_4Q_Ma", "count_ss", "sub_lgd_4Q_Ma", "count_sub", macro_base)
               ]

# Apply time horizon
ss_sub_macro = 
  ss_sub_macro[
    match(start, ss_sub_macro$Quarter):match(end, ss_sub_macro$Quarter),
    ]

# Save output and delete intermediate datasets
write.csv(
  ss_sub_macro, 
  paste0(output,"/ss_sub_macro_out.csv"), 
  row.names = F, 
  na = ''
  )
rm(quarters, ss, ss_count, ss_series, ss_sub, sub, sub_count, sub_series, ur)

