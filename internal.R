###################################################################################################################
#
#                                            INTERNAL LGD ANALYSIS
#
###################################################################################################################

wd = '...'
code = paste0(wd, '/CODE')
input = paste0(wd, '/INPUT')
output = paste0(wd, '/OUTPUT')

source(paste0(code,'/LGD Functions.R'))


#==================================================================================================================
# INTERNAL LGDs
#==================================================================================================================

# Load internal data
internal = 
  read.csv(
    paste0(input, '/Full_V1_Dataset - handover.csv'), 
    header = T, 
    stringsAsFactors = F
    )

# Cap LGD to 100%
internal$lgd = pmin(internal$Accounting.LGD, 1)
internal$lgd = pmax(internal$lgd, 0)

# Create three segments: Equipment, Secured, Unsecured
internal$BROAD_SENIORITY_COLLATERAL = as.character(internal$Final_collateral_type_grouped)
internal[
  internal$BROAD_SENIORITY_COLLATERAL == "Blank", 
  "BROAD_SENIORITY_COLLATERAL"
  ] = "Junior/Unsecured"
internal[
  internal$BROAD_SENIORITY_COLLATERAL == "Unsecured", 
  "BROAD_SENIORITY_COLLATERAL"
  ] = "Junior/Unsecured"
internal[
  internal$BROAD_SENIORITY_COLLATERAL == "Machinery and Equipment", 
  "BROAD_SENIORITY_COLLATERAL"
  ] = "Equipment"
internal[
  !internal$BROAD_SENIORITY_COLLATERAL %in% c("Equipment","Junior/Unsecured"),
  "BROAD_SENIORITY_COLLATERAL"
  ] = "Senior Secured"

# Add a flag for quarter th
internal$Quarter = as.numeric(substr(internal$Default_Year_Month_ID, 5 ,6))
internal[internal$Quarter %in% c(1,2,3),"Quarter"] = "Q1"
internal[internal$Quarter %in% c(4,5,6),"Quarter"] = "Q2"
internal[internal$Quarter %in% c(7,8,9),"Quarter"] = "Q3"
internal[internal$Quarter %in% c(10,11,12),"Quarter"] = "Q4"
internal$Quarter = 
  paste(
    substr(internal$Default_Year_Month_ID, 1, 4), ':', internal$Quarter, 
    sep = ""
    )
internal = internal[order(internal$Quarter),]

write.csv(
  internal, 
  file=paste0(output,'/internal_out.csv'), 
  row.names = F, 
  na = ""
  )


#============================================================================================================
# TIME SERIES ANALYSIS
#============================================================================================================
# Create a 4 quarter count-weighted time series LGD for each segment for correlation analysis later on

# Create Segments
ss = internal[internal$BROAD_SENIORITY_COLLATERAL=="Senior Secured",]
sub = internal[internal$BROAD_SENIORITY_COLLATERAL=="Junior/Unsecured",]
equip = internal[internal$BROAD_SENIORITY_COLLATERAL=="Equipment",]
non_equip = internal[internal$BROAD_SENIORITY_COLLATERAL!="Equipment",]

# Sort by Quarter
ss = ss[order(ss$Quarter),]
sub = sub[order(sub$Quarter),]
equip = equip[order(equip$Quarter),]
non_equip = non_equip[order(non_equip$Quarter),]

# Create time series of LGD
ss_series = 
  aggregate(
    x = ss$lgd, 
    by=list(ss$Quarter), 
    mean
    )
sub_series = 
  aggregate(
    x = sub$lgd, 
    by=list(sub$Quarter), 
    mean
    )
equip_series = 
  aggregate(
    x = equip$lgd, 
    by=list(equip$Quarter), 
    mean
    )
non_equip_series = 
  aggregate(
    x = non_equip$lgd, 
    by=list(non_equip$Quarter),
    mean
    )
colnames(ss_series) = c("Quarter", "lgd_ss")
colnames(sub_series) = c("Quarter", "lgd_sub")
colnames(equip_series) = c("Quarter", "lgd_equip")
colnames(non_equip_series) = c("Quarter", "lgd_non_equip")

# Add LGD observation count for segments
ss_count = 
  aggregate(
    x = ss$lgd, 
    by=list(ss$Quarter),
    length
    )
sub_count = 
  aggregate(
    x = sub$lgd, 
    y=list(sub$Quarter), 
    length
    )
equip_count = 
  aggregate(
    x = equip$lgd, 
    by=list(equip$Quarter),
    length
    )
non_equip_count = 
  aggregate(
    x = non_equip$lgd,
    by=list(non_equip$Quarter),
    length
    )
colnames(ss_count) = c("Quarter", "count_ss")
colnames(sub_count) = c("Quarter", "count_sub")
colnames(equip_count) = c("Quarter", "count_equip")
colnames(non_equip_count) = c("Quarter", "count_non_equip")

# Add missing quarters in order to have a complete time series
quarters =
  read.csv(
    file=paste0(input,'/quarters.csv'),
    stringsAsFactors = F,
    col.names = 'Quarter'
    )

# Merge mean lgd and count series for each segment
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
equip_series = 
  merge(
    x = quarters, 
    y = equip_series,
    by = "Quarter", 
    all.x = T, 
    all.y = F
    )
equip_series = 
  merge(
    x = equip_series, 
    y = equip_count, 
    by = "Quarter", 
    all.x = T, 
    all.y = F
    )
non_equip_series = 
  merge(
    x = quarters, 
    y = non_equip_series, 
    by = "Quarter", 
    all.x = T, 
    all.y = F
    )
non_equip_series = 
  merge(
    x = non_equip_series, 
    y = non_equip_count, 
    by = "Quarter", 
    all.x = T, 
    all.y = F
    )

# Merge into one dataset for futher analysis
internal_ss_sub = 
  merge(
    ss_series, 
    sub_series, 
    by= "Quarter"
    )
internal_ss_sub = 
  merge(
    internal_ss_sub, 
    equip_series, 
    by= "Quarter"
    )
internal_ss_sub = 
  merge(
    internal_ss_sub, 
    non_equip_series, 
    by= "Quarter"
    )

# Calculate 4 quarter count weighted moving average LGD
internal_ss_sub$ss_lgd_4Q_Ma = 
  calc_4qtr_ms(internal_ss_sub[, "lgd_ss"] * internal_ss_sub[,"count_ss"])/
  calc_4qtr_ms(internal_ss_sub$count_ss)
internal_ss_sub$sub_lgd_4Q_Ma = 
  calc_4qtr_ms(internal_ss_sub[, "lgd_sub"] * internal_ss_sub[,"count_sub"])/
  calc_4qtr_ms(internal_ss_sub$count_sub)
internal_ss_sub$equip_lgd_4Q_Ma = 
  calc_4qtr_ms(internal_ss_sub[, "lgd_equip"] * internal_ss_sub[,"count_equip"])/
  calc_4qtr_ms(internal_ss_sub$count_equip)
internal_ss_sub$non_equip_lgd_4Q_Ma = 
  calc_4qtr_ms(internal_ss_sub[, "lgd_non_equip"] * internal_ss_sub[,"count_non_equip"])/
  calc_4qtr_ms(internal_ss_sub$count_non_equip)

# Keep only target variables
internal_ss_sub = 
  internal_ss_sub[,c("Quarter","ss_lgd_4Q_Ma","count_ss","sub_lgd_4Q_Ma","count_sub","equip_lgd_4Q_Ma",
                                     "count_equip","non_equip_lgd_4Q_Ma","count_non_equip")]

# Set up time horizon
start = '2007:Q4'
end = '2015:Q2'
internal_ss_sub = 
  internal_ss_sub[
    match(start, internal_ss_sub$Quarter):match(end, internal_ss_sub$Quarter),
    ]

# Save output
write.csv(
  internal_ss_sub, 
  paste0(output, "/internal_ss_sub_macro_out.csv"), 
  row.names = F, 
  na = ''
  )

rm(quarters, ss, ss_count, ss_series, sub, sub_count, sub_series, equip, equip_count, equip_series, 
   non_equip, non_equip_count, non_equip_series)


