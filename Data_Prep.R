##################################################################################################################
#
#                                              BASE DATA PREPARATION
#
##################################################################################################################
wd = '...'
input = paste0(wd, '/INPUT')
output = paste0(wd, '/OUTPUT')

# Load packages
library(data.table)
library(lubridate)    # To get year and month from date variables

#=================================================================================================================
#=================================================================================================================

#==================== Read and merge raw ultimate reovery data ====================#

# Variable names to be selected
recovery_class_names = 
  c(
    'CLASS_ID','CLASS_NAM','CLASS_TAG','EVENT_ID','PREFERRED_METHOD'
    )
recovery_instrument_names = 
  c(
    'CLASS_ID','DEBT_ABOVE','DEBT_BELOW','DEBT_CURRENT','DISCOUNT_RECOMMENDED','EFFECTIVE_INTEREST_RATE',
    'INSTRUMENT_ID','LAST_CASH_PAID_DATETIME','MAST_ISSU_NUM','TYPE','DESC','DID_NOT_DEFAULT_F',
    'NOMINAL_SETTLEMENT_TOTAL','DISCOUNT_SETTLEMENT_TOTAL','NOMINAL_TRADING_PRICE','DISCOUNT_TRADING_PRICE',
    'NOMINAL_LIQUIDITY_TOTAL', 'DISCOUNT_LIQUIDITY_TOTAL','COLLATERAL_RANK','COLLATERAL_TYPE','PERCENT_ABOVE',
    'PERCENT_BELOW','PERCENT_CURRENT','PRINCIPAL_AMOUNT_AT_DEFAULT','TRADING_PRICE_30_DAY',
    'TRADING_PRICE_30_DAY_DATETIME','TRADING_PRICE_EMERGENCE','TRADING_PRICE_EMERGENCE_DT'
    )
recovery_oblgior_names = 
  c(
    'MAST_ISSR_NUM','NAM','OBLIGOR_ID','CUSIP','TICKER','SIC_CODE','SIC_DESC'
    )
recovery_event_names = 
  c(
    'DEF_NUM','EVENT_ID','FAMILY_RECOVERY','EMERGENCE_DATETIME','ISSUER_DEFAULT_DATETIME','OBLIGOR_ID',
    'OUTCOME','TOTAL_DEBT','DEFAULT_EVENT_DESC'
    )

# Read ultimate recovery (ur) data and read only selected variable names
recovery_class = 
  fread(
    paste0(input,'/RECOVERY_CLASS.txt'), 
    stringsAsFactors = F, 
    integer64 = 'double', 
    select = recovery_class_names
    )
recovery_instrument = 
  fread(
    paste0(input,'/RECOVERY_INSTRUMENT.txt'),
    stringsAsFactors = F, 
    integer64 = 'double', 
    select = recovery_instrument_names
    )
recovery_obligor = 
  fread(
    paste0(input,'/RECOVERY_OBLIGOR.txt'), 
    stringsAsFactors = F, 
    integer64 = 'double', 
    select = recovery_oblgior_names
    )
recovery_event = 
  fread(
    paste0(input,'/RECOVERY_EVENT.txt'), 
    stringsAsFactors = F, 
    integer64 = 'double', 
    select = recovery_event_names
    )

# Merge data read above into a single recovery dataset called "ur_info"
ur_info = 
  merge(
    recovery_instrument, 
    recovery_class, 
    by='CLASS_ID', 
    all.x=T, 
    all.y=F
    )
ur_info = 
  merge(
    ur_info, 
    recovery_event, 
    by='EVENT_ID', 
    all.x=T, 
    all.y=F
    )
ur_info = 
  merge(
    ur_info, 
    recovery_obligor, 
    by='OBLIGOR_ID', 
    all.x=T, 
    all.y=F
    )

# Change dates to "Date" format
ur_info[, c('LAST_CASH_PAID_DATETIME', 'TRADING_PRICE_30_DAY_DATETIME', 'TRADING_PRICE_EMERGENCE_DT', 
            'ISSUER_DEFAULT_DATETIME', 'EMERGENCE_DATETIME') 
        := lapply(.(LAST_CASH_PAID_DATETIME, TRADING_PRICE_30_DAY_DATETIME, TRADING_PRICE_EMERGENCE_DT, 
                    ISSUER_DEFAULT_DATETIME, EMERGENCE_DATETIME), as.Date, format='%m/%d/%Y')]

#==================== Add new variables which are needed in later analysis ====================#

# Add default year and quarter info
ur_info[, ISSUER_DEFAULT_YEAR:=year(ISSUER_DEFAULT_DATETIME)]
ur_info[, ISSUER_DEFAULT_MONTH:=month(ISSUER_DEFAULT_DATETIME)]
ur_info[!is.na(ISSUER_DEFAULT_YEAR), Quarter := 
          paste0(
            ISSUER_DEFAULT_YEAR, ':', ifelse(
              ISSUER_DEFAULT_MONTH %in% c(1,2,3), 
              'Q1',
              ifelse(
                ISSUER_DEFAULT_MONTH %in% c(4,5,6), 
                'Q2', 
                ifelse(
                  ISSUER_DEFAULT_MONTH %in% c(7,8,9), 
                  'Q3', 
                  'Q4'
                  )
                )
              )
            )
        ]

# Add recommended recovery info based on perferred method
ur_info[, RECOMMENDED_NOMINAL_RECOVERY:=
          ifelse(
            PREFERRED_METHOD=='Trading Price', 
            NOMINAL_TRADING_PRICE,
            ifelse(
              PREFERRED_METHOD=='Liquidity',
              NOMINAL_LIQUIDITY_TOTAL,
              ifelse(
                PREFERRED_METHOD=='Settlement',
                NOMINAL_SETTLEMENT_TOTAL, 
                NA
                )
              )
            ), 
        by=1:nrow(ur_info)]
ur_info[, RECOMMENDED_DISCOUNT_RECOVERY:=
          ifelse(
            PREFERRED_METHOD=='Trading Price',
            DISCOUNT_TRADING_PRICE,
            ifelse(
              PREFERRED_METHOD=='Liquidity',
              DISCOUNT_LIQUIDITY_TOTAL,
              ifelse(
                PREFERRED_METHOD=='Settlement',
                DISCOUNT_SETTLEMENT_TOTAL, 
                NA
                )
              )
            ), 
        by=1:nrow(ur_info)]

# Bound LGD between [0, 1]; cap recoveries to 100 because anything > 100 is workout cost reimbursement
ur_info$RECOMMENDED_NOMINAL_RECOVERY = 
  pmin(
    ur_info$RECOMMENDED_NOMINAL_RECOVERY, 
    100
    )
ur_info$RECOMMENDED_NOMINAL_RECOVERY = 
  pmax(
    ur_info$RECOMMENDED_NOMINAL_RECOVERY, 
    0
    )

# Calculate LGD
ur_info$lgd = 1 - (ur_info$RECOMMENDED_NOMINAL_RECOVERY/100)

# Calculate the number of days between the last cash paid by the obligor and the date of defaulted
ur_info[, days_passed_btw_lastcash_and_dflt:=
          as.numeric(ISSUER_DEFAULT_DATETIME-LAST_CASH_PAID_DATETIME)]

# Calculate the number of days it takes the obligor to go from default to emergence
ur_info[, days_passed_btw_dflt_and_emergence:=
          as.numeric(EMERGENCE_DATETIME-ISSUER_DEFAULT_DATETIME)]

# Sort column in alphabetical order so it is easier to locate variable
setcolorder(ur_info, sort(names(ur_info)))
setkey(ur_info, NULL)


#===================================================================================================================
# APPLY FILTERS
#===================================================================================================================

#==================== Exclude foreign obligors ====================#

# Get domain(country) info from Moody's DRD data
mast_issr = 
  fread(
    paste0(input,'/MAST_ISSR.txt'), 
    stringsAsFactors = F, 
    integer64 = 'double', 
    select = c('MAST_ISSR_NUM','MDY_DOMN_NUM')
    )
govt_domain = 
  fread(
    paste0(input,'/GOVT_DOMAIN.txt'),
    stringsAsFactors = F, 
    integer64 = 'double', 
    select = c('DOMN_NAM','MDY_DOMN_NUM')
    )
domain_info = 
  merge(
    mast_issr, 
    govt_domain, 
    by='MDY_DOMN_NUM', 
    all.x=T, 
    all.y=F
    )

# Merge Domicile into ur_info
ur_info = 
  merge(
    ur_info, 
    unique(domain_info[,.(MAST_ISSR_NUM, DOMN_NAM)]), 
    by='MAST_ISSR_NUM', 
    all.x=T, 
    all.y=F
    )

# Exclude foreign obligors, but keep N/A because N/A's have been verified to be US
ur_info = ur_info[DOMN_NAM=='UNITED STATES'|is.na(DOMN_NAM)]

#==================== Apply standard cleanup filters ====================#

ur_info = ur_info[DID_NOT_DEFAULT_F!=1]
ur_info = ur_info[PRINCIPAL_AMOUNT_AT_DEFAULT>10000 | is.na(PRINCIPAL_AMOUNT_AT_DEFAULT)]
# Exclude less than 90 days delinquencies for consistency with Basel default, unless they are bankruptcy, 
# destressed exchange or other restructuring
ur_info = 
  ur_info[
    DEFAULT_EVENT_DESC %in% c("Bankruptcy","Distressed Exchange", "Other Restructuring") |
      days_passed_btw_lastcash_and_dflt >= 90
    ]
# Assuming payment schedule is monthly, exclude default and cure facilities that are less than 120 (30 + 90) days
# delinquent; otherwise if we assume quarterly payments, then use 180 (90 + 90) days
ur_info = 
  ur_info[
    !(DEFAULT_EVENT_DESC=='Default and Cure' & 
        days_passed_btw_lastcash_and_dflt < 120)
    ]
ur_info = ur_info[COLLATERAL_TYPE!='Real Estate']

# Save ur_info before further filtering
write.csv(
  ur_info, 
  file=paste0(output, '/output_ur_all_types.csv'), 
  row.names=F, 
  na=''
  )
ur_info[, .(.N, Avg.LGD=mean(lgd, na.rm=T)), by='TYPE']    # Quick summary of counts and average LGD

# Keep only Term Loans or Revolvers
ur_info = ur_info[TYPE %in% c('Term Loan', 'Revolver')]


#===================================================================================================================
# ADDITION OF FIELDS FOR SEGMENTATION ANALYSIS
#===================================================================================================================

#==================== Add Industry information ====================#

# Map SIC code to Moodys NDY and industry descriptions
sic_mapping = 
  fread(
    paste0(input,'/SIC_Mapping.csv'), 
    stringsAsFactors=F, 
    integer64 = 'double'
    )
sic_mapping[, SIC:=as.character(SIC)]
ur_info[, SIC_CODE:=as.character(SIC_CODE)]
ur_info = 
  merge(
    ur_info, 
    sic_mapping, 
    by.x='SIC_CODE', 
    by.y='SIC', 
    all.x=T, 
    all.y=F
    )

# Add FINAL_INDUSTRY; Carve out Energy from the standard SIC Broad categories to craete an Energy segment.
ur_info[, FINAL_INDUSTRY:=Sector4.0]
ur_info[SIC.2 %in% c('12','13','29','46'), FINAL_INDUSTRY:='Energy']

#==================== Add Seniority information ====================#

# Assign a dummy variable to all senior and subordinate debts
ur_info[, SENIORITY:=ifelse(DEBT_ABOVE>0, 'Junior', 'Senior')]

# Assign a dummy variable to collateral types for Secured and Unsecured
ur_info[, COLLATERAL:=COLLATERAL_TYPE]
ur_info[
  COLLATERAL_TYPE %in% c('Accounts Receivable',
                         'All Assets',
                         'All Current Assets',
                         'All Non-current Assets',
                         'Cash',
                         'Equipment',
                         'Inventory',
                         'Inventory and Accounts Receivable',
                         'Most Assets',
                         'Oil and Gas Properties',
                         'PP&E',
                         'Real Estate'),
  COLLATERAL:='Secured'
  ]

ur_info[, COLLATERAL:=
          ifelse(
            COLLATERAL!='Secured', 
            'Unsecured', 
            'Secured'
            )
        ]

#==================== Add segments for All Assets and Other Secured ====================#

ur_info[, ALL_ASSETS :=COLLATERAL_TYPE]
ur_info[
  COLLATERAL_TYPE %in% c('All Assets','All Current Assets','All Non-current Assets','Most Assets'),
  ALL_ASSETS:='All_Assets'
  ]
ur_info[
  COLLATERAL_TYPE %in% c('Accounts Receivable','Cash','Equipment','Inventory',
                         'Inventory and Accounts Receivable','Oil and Gas Properties','PP&E','Real Estate'),
  ALL_ASSETS:='Other_Assets'
  ]

ur_info[
  COLLATERAL_TYPE %in% c('Unsecured','Capital Stock','Second Lien','Third Lien','Guarantees',
                         'Intellectual Property'),  
  ALL_ASSETS:='Unsecured'
  ]

#==================== Add segments for Services and Non-Services ====================#

service_desc_spec = 
  c('Railroad Transportation',
    'Local And Suburban Transit And Interurban Highway Passenger Transportation',
    'Water Transportation',
    'Transportation By Air',
    'Transportation Services',
    'Communications',
    'Eating And Drinking Places',
    'Depository Institutions',
    'Non-depository Credit Institutions',
    'Security And Commodity Brokers, Dealers, Exchanges, And Services',
    'Real Estate',
    'Holding And Other Investment Offices',
    'Hotels, Rooming Houses, Camps, And Other Lodging Places',
    'Personal Services',
    'Business Services',
    'Automotive Repair, Services, And Parking',
    'Miscellaneous Repair Services',
    'Motion Pictures',
    'Amusement And Recreation Services',
    'Health Services',
    'Legal Services',
    'Social Services',
    'Engineering, Accounting, Research, Management, And Related Services',
    'Miscellaneous Services',
    'Nonclassifiable Establishments'
    )
ur_info[, Services := 'Non-Services']
ur_info[SIC.2_DESC_SPEC %in% service_desc_spec, Services := 'Services']

#==================== Create new variables based on different variable combination ====================#

# Create combinations of Seniority and COllateral (i.e., Senior Secured, etc.)
ur_info[
  SENIORITY=='Senior' & COLLATERAL=='Secured', 
  SENIORITY_COLLATERAL:="Senior Secured"
  ]
ur_info[
  SENIORITY=='Senior' & COLLATERAL=='Unsecured', 
  SENIORITY_COLLATERAL:="Senior Unsecured"
  ]
ur_info[
  SENIORITY=='Junior' & COLLATERAL=='Secured', 
  SENIORITY_COLLATERAL:="Junior Secured"
  ]
ur_info[
  SENIORITY=='Junior' & COLLATERAL=='Unsecured', 
  SENIORITY_COLLATERAL:="Junior Unsecured"
  ]

# Create combination of 1) Senior Secured and 2) all Juniors or all Unsecured
ur_info[
  SENIORITY_COLLATERAL=='Senior Secured', 
  BROAD_SENIORITY_COLLATERAL:="Senior Secured"
  ]
ur_info[
  is.na(BROAD_SENIORITY_COLLATERAL), 
  BROAD_SENIORITY_COLLATERAL:="Junior/Unsecured"
  ]

# Create combination of Type and Collateral (i.e., Secured Revolver, Unsecured Term Loan, etc.)
ur_info[
  COLLATERAL=='Secured' & TYPE=='Revolver', 
  COLLAT_TYPE:="Secured Revolver"
  ]
ur_info[
  COLLATERAL=='Secured' & TYPE=='Term Loan', 
  COLLAT_TYPE:="Secured Term Loan"
  ]
ur_info[
  COLLATERAL=='Unsecured' & TYPE=='Revolver', 
  COLLAT_TYPE:="Unsecured Revolver"
  ]
ur_info[
  COLLATERAL=='Unsecured' & TYPE=='Term Loan', 
  COLLAT_TYPE:="Unsecured Term Loan"
  ]
ur_info[
  is.na(COLLAT_TYPE), COLLAT_TYPE:="Other"]

# Create combinations of all three dummy variables (i.e., Senior Secured Revolver, etc.)
ur_info[
  SENIORITY_COLLATERAL=='Senior Secured' & TYPE=='Revolver', 
  SENIORITY_COLLAT_TYPE:="Senior Secured Revolver"
  ]
ur_info[
  SENIORITY_COLLATERAL=='Senior Secured' & TYPE=='Term Loan', 
  SENIORITY_COLLAT_TYPE:="Senior Secured Term Loan"
  ]
ur_info[
  SENIORITY_COLLATERAL=='Senior Unsecured' & TYPE=='Revolver', 
  SENIORITY_COLLAT_TYPE:="Senior Unsecured Revolver"
  ]
ur_info[
  SENIORITY_COLLATERAL=='Senior Unsecured' & TYPE=='Term Loan', 
  SENIORITY_COLLAT_TYPE:="Senior Unsecured Term Loan"
  ]
ur_info[
  SENIORITY_COLLATERAL=='Junior Secured' & TYPE=='Revolver', 
  SENIORITY_COLLAT_TYPE:="Junior Secured Revolver"
  ]
ur_info[
  SENIORITY_COLLATERAL=='Junior Secured' & TYPE=='Term Loan', 
  SENIORITY_COLLAT_TYPE:="Junior Secured Term Loan"
  ]
ur_info[
  SENIORITY_COLLATERAL=='Junior Unsecured' & TYPE=='Revolver', 
  SENIORITY_COLLAT_TYPE:="Junior Unsecured Revolver"
  ]
ur_info[
  SENIORITY_COLLATERAL=='Junior Unsecured' & TYPE=='Term Loan', 
  SENIORITY_COLLAT_TYPE:="Junior Unsecured Term Loan"
  ]
ur_info[
  is.na(SENIORITY_COLLAT_TYPE), 
  SENIORITY_COLLAT_TYPE:="Other"
  ]


#==================================================================================================================
# SAVE FINAL ULTIMATE RECOVERY DATA FOR LATER USE
#==================================================================================================================

#==================== Save all variables as in ur_info ====================#

setorderv(ur_info, c('Quarter','ISSUER_DEFAULT_DATETIME'))
write.csv(
  ur_info, 
  paste0(output,'/output_ur_all.csv'), 
  row.names = F, 
  na = ''
  )

#==================== Subset ur_info and keep only selected variables for further analysis ====================#

keep_list = 
  c(
    'Quarter','ISSUER_DEFAULT_DATETIME','lgd','INSTRUMENT_ID','FINAL_INDUSTRY','TYPE','SENIORITY',
    'COLLATERAL', 'ALL_ASSETS', 'Services','SENIORITY_COLLATERAL','BROAD_SENIORITY_COLLATERAL',
    'COLLAT_TYPE','SENIORITY_COLLAT_TYPE'
    )
ur = ur_info[, keep_list, with=F]
setorderv(ur, c('Quarter','ISSUER_DEFAULT_DATETIME'))
write.csv(
  ur, 
  paste0(output,'/output_ur.csv'), 
  row.names = F, 
  na = ''
  )

