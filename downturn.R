##########################################################################################
#
#                                       DOWNTURN
#
##########################################################################################

wd = '...'
code = paste0(wd, '/CODE')
input = paste0(wd, '/INPUT')
output = paste0(wd, '/OUTPUT')

source(paste0(code,'/LGD Functions.R'))

#==========================================================================================
# CORP DEFAULT RATE
# APPROACH 1: ALL SPECULATIVE US ISSUERS
#==========================================================================================

def = 
  read.csv(
    paste0(input,'corp_default_data.csv'), 
    header = F, 
    skip = 10, 
    nrows = 7
    )
def = def[,-1]
def = t(def)

colnames(def) = def[1,]
def = as.data.frame(def[-1,])

def$date = 
  as.Date(
    def$`Up/Down`, 
    format = "%m/%d/%Y"
    )

# Assign quarter
def$Quarter = as.yearqtr(def$date, "%m/%d/%Y")
def$Quarter = 
  paste(
    substr(def$Quarter, 1, 4), ':', substr(def$Quarter, 6, 7), 
    sep=""
    )

# Format as numeric for further data processing
index = sapply(def, is.factor)
def[index] = 
  lapply(
    def[index], 
    function(x) as.numeric(as.character(x))
    )

# Re-arrange the column order with quarter first
def = cbind(date=def$date, Quarter=def$Quarter, def[,2:7])

# Calculate the historical default rate (Defaults/Totals for each quarter)
def$def_rate = def$Defaulted/(def$Totals - def$`Non-Rated`)

#================ Downturn period calculation ================#

# Standard Deviation approach (def rate > 1 standard deviation above mean)
sd = 
  mean(
    def$def_rate, 
    na.rm = TRUE
    ) + 
  sd(
    def$def_rate, 
    na.rm = TRUE
    )
def$downturn.sd = as.numeric(def$def_rate > sd)

# 80th Percentile approach (def rate > 80th percentile def rate)
pct = 
  as.double(
    quantile(
      def$def_rate,
      probs =  0.8, 
      na.rm = TRUE
      )
    )
def$downturn.pct = as.numeric(def$def_rate > pct)

# Median approach (def rate > 2 x Median)
med = median(def$def_rate, na.rm  = TRUE) * 2
def$downturn.med = as.numeric(def$def_rate > med)


#===================================================================================
# S&P CORP DEFAULT RATE
# APPROACH 2: RATING WEIGHTED DEFAULT RATES
#===================================================================================

#' Conduct autocorrelation tests to model residuals
#'
#' @author Gerard Mazi
#' @param data Raw dataset for a specific rating
#' 
#' @return Return a cleaned and formotted dataset with default rate column added
#' @export
#'
#' @examples
#' TBD
#'
get_default_rate = function(data, rating) {
  options(warn = -1)
  data = data[,-1]
  data = t(data)
  
  colnames(data) = data[1,]
  data = as.data.frame(data[-1,])
  
  data$date = as.Date(data$`Up/Down`, format = "%m/%d/%Y")
  
  # Assign quarter
  data$Quarter = as.yearqtr(data$date, "%m/%d/%Y")
  data$Quarter = 
    paste(
      substr(data$Quarter, 1, 4), ':', substr(data$Quarter, 6, 7), 
      sep=""
      )
  
  # Format as numeric for further data processing
  index = sapply(data, is.factor)
  data[index] = 
    lapply(
      data[index], 
      function(x) as.numeric(as.character(x))
      )
  
  # Re-arrange the column order with quarter first
  data = 
    cbind(
      date=data$date, 
      Quarter=data$Quarter, 
      data[,2:7]
      )
  
  # Calculate the historical default rate (Defaults/Totals for each quarter)
  data$data_def_rate = data$Defaulted/(data$Totals - data$`Non-Rated`)
  setnames(data, 'data_def_rate', paste0(rating, '_def_rate'))
  
  options(warn = 0)
  return(data)
}

#================ Get default rate for each rating  ================#

rating_list = c('AAA','AA','A','BBB','BB','B','CCC','CC','C')

dt = list()
for (rating in rating_list) {
  dt[[rating]] = 
    read.csv(
      paste0(input,'/', rating, '.csv'), 
      header = F, 
      skip = 10, 
      nrows = 7, 
      stringsAsFactors = F
      )
  dt[[rating]] = get_default_rate(dt[[rating]], rating)
}

AAA = dt$AAA
AA = dt$AA
A = dt$A
BBB = dt$BBB
BB = dt$BB
B = dt$B
CCC = dt$CCC
CC = dt$CC
C = dt$C

Defaulted = CCC$Defaulted + CC$Defaulted + C$Defaulted
Non_Rated = CCC$`Non-Rated` + CC$`Non-Rated` + C$`Non-Rated`
Totals = CCC$Totals + CC$Totals + C$Totals

CCC_C_default_rate = Defaulted / (Totals - Non_Rated)

CCC_C = data.frame(Quarter = CCC$Quarter, CCC_C_default_rate)

# Bring all the rating-specific default rates in one data frame
rating_def = 
  data.frame(
    Quarter = AAA$Quarter,
    AAA = AAA$AAA_def_rate, 
    AA = AA$AA_def_rate, 
    A = A$A_def_rate,
    BBB = BBB$BBB_def_rate, 
    BB = BB$BB_def_rate, 
    B = B$B_def_rate, 
    CCC_C = CCC_C$CCC_C_default_rate, 
    stringsAsFactors = F
    )

#================ Get rating weighted default rates  ================#

# COF Portfolio Weights
weights = 
  read.table(
    paste0(input, '/weight_201612.txt'), 
    stringsAsFactors = F
    )
colnames(weights) = c("COF_rating", "count")

# Aggregate COF portfolio by the following mapping:
# AAA-1, AA-2, A-3, BBB-4~6, BB-7~11, B-12~15, CCC/CC/C-16~18
rating_def$COF_AAA_weight = 0
rating_def$COF_AA_weight = 0
rating_def$COF_A_weight = sum(weights[weights$COF_rating==3,2])/sum(weights[,2])
rating_def$COF_BBB_weight = sum(weights[weights$COF_rating %in% c(4,5,6),2])/sum(weights[,2])
rating_def$COF_BB_weight = sum(weights[weights$COF_rating %in% c(7,8,9,10,11),2])/sum(weights[,2])
rating_def$COF_B_weight = sum(weights[weights$COF_rating %in% c(12,13,14,15),2])/sum(weights[,2])
rating_def$COF_CCC_C_weight = sum(weights[weights$COF_rating %in% c(16,17,18),2])/sum(weights[,2])

# Calculate weighted average default rate based on COF portfolio weights
rating_def$average_def_rate = 
  rating_def$AAA * rating_def$COF_AAA_weight +
  rating_def$AA * rating_def$COF_AA_weight +
  rating_def$A * rating_def$COF_A_weight +
  rating_def$BBB * rating_def$COF_BBB_weight +
  rating_def$BB * rating_def$COF_BB_weight +
  rating_def$B * rating_def$COF_B_weight +
  rating_def$CCC_C * rating_def$COF_CCC_C_weight

# Cleanup  
rating_def = rating_def[,c("Quarter","average_def_rate")]
rm(AAA,AA,A,BBB,BB,B,CCC,CC,C,CCC_C,weights)

#================ Calculate downturn period  ================#

def$average_def_rate = rating_def$average_def_rate

# Standard Deviation approach (def rate > 1 standard deviation above mean)
sd_w = 
  mean(
    rating_def$average_def_rate, 
    na.rm = TRUE
    ) + 
  sd(
    rating_def$average_def_rate, 
    na.rm = TRUE
    )
def$downturn.sd_w = as.numeric(rating_def$average_def_rate > sd_w)

# 80th Percentile approach (def rate > 80th percentile def rate)
pct_w = 
  as.double(
    quantile(
      rating_def$average_def_rate,
      probs =  0.8, 
      na.rm = TRUE
      )
    )
def$downturn.pct_w = as.numeric(rating_def$average_def_rate > pct_w)

# Median approach (def rate > 2 x Median)
med_w = 
  median(
    rating_def$average_def_rate, 
    na.rm  = TRUE
    ) * 2
def$downturn.med_w = as.numeric(rating_def$average_def_rate > med_w)

write.csv(
  def, 
  file=paste0(output, '/def_out.csv'), 
  row.names = F
  )

