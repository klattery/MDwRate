stacker <- function(kdata, head_cols){
  # General stacking function
  kdata <- as.matrix(kdata) # force kdata as matrix
  result <- do.call(rbind,
                    lapply(1:nrow(kdata), function(i){
                      head <- kdata[i,head_cols] 
                      values <- kdata[i,-head_cols] 
                      counter <- 1:length(values)
                      head_s <- do.call(rbind, lapply(counter, function(x) head)) # repeat heading
                      result <- cbind(head_s, counter, values)
                    })) 
  
}



merge_function <- function(data_master, des_version, module_name) {
#######################################################
#  Setup MaxDiff with Design File + Choices
#######################################################

####################################################################
##       1.  SPECIFY DIRECTORIES, DATA, & COLUMNS                 ##
####################################################################

# directory <- "C:/Users/m.kang/SKIM/AT&T - Master Database and Dashboard/Analysis/0. Data/"

####################################################################
##       2.  SPECIFY DESIGN FILE AND SURVEY DATA FILE              ##
####################################################################
# des_version <- as.matrix(read.csv(paste0(directory, file = "F0000_M_Design.csv"), as.is = TRUE)) # Version
# data_master <- read.csv(paste0(directory, file = "F0000_data_for_analysis.csv"), as.is = TRUE)
des_version <- as.data.frame(des_version)
data_master <- as.data.frame(data_master)

#####################################################################
##       3.  GET ID, VERSION, BEST PICKS, WORST PICKS, Anchor Data ##
#####################################################################
id <- data_master[,1]                           # id variable
moduleversiongrab = paste0("sys_MXDVersion_",module_name)
ver <- data_master[,moduleversiongrab] # version variable

bests <- paste0(module_name,".*_b")
worsts <- paste0(module_name,".*_w")

MD_best <- data_master[,grep(bests, colnames(data_master))] # Best choices in: M_*_b
MD_worst <- data_master[,grep(worsts, colnames(data_master))] # worst choices in: M_*_w
# You should be sure that the 2 files above have the right variables

#############################################################################
##    4.  RUN REST OF CODE                                        
##     AFTER LINE 80 CHECK THAT YOUR DESIGN AND DATA FILE MATCH
##     FILE bad SHOULD HAVE 0 ROWS OF DATA
##     OTHERWISE YOUR DESIGN FILE DOES NOT MATCH YOUR DATA FILE ON bad 
#############################################################################

#######################################################################
#          PREPARE MAXDIFF DATA
######################################################################
des_version_long <- stacker(des_version, 1:2)
colnames(des_version_long) <- c("version", "task", "concept", "items")
# Stack tasks as long data

ver_split <- split(as.data.frame(des_version_long), des_version_long[,1]) # Each version in list
id_ver <- cbind(id, ver)
resp_design <- do.call(rbind, lapply(1:nrow(id_ver), function(i){
  id <- id_ver[i,1]
  ver <- id_ver[i,2]
  shown <- ver_split[[ver]]
  id <- rep(id, nrow(shown))
  result <- cbind(id, ver_split[[ver]])  
})) # What each respondent saw

best_long <- stacker(cbind(id, MD_best), 1)
worst_long <- stacker(cbind(id, MD_worst), 1)
resp_picks <- as.data.frame(cbind(best_long, worst_long[,3]))
colnames(resp_picks) <- c("id", "task", "best", "worst")
# resp_picks has stacked MD choices

data_merge <- merge(resp_design, resp_picks, by.x = c("id", "task"), by.y = c("id", "task"), all.x = TRUE)
korder <- order(data_merge$id, data_merge$task, data_merge$concept)
data_merge <- data_merge[korder,]
# SQL Code
#data_merge <- sqldf("select resp_design.*, resp_picks.best, resp_picks.worst from
#              resp_design left join resp_picks on
#              (resp_design.id = resp_picks.id) and
#              (resp_design.task = resp_picks.task)") # Merge with design
dep_bw <- cbind(1 * (data_merge$items == data_merge$best), 1 * (data_merge$items == data_merge$worst))
colnames(dep_bw) <- c("dep_best", "dep_worst")
data_merge <- cbind(data_merge, dep_bw)
check <- aggregate(dep_bw, list(id = data_merge$id, task = data_merge$task), sum)
bad <- check[rowSums(check[,3:4]) < 2,] #Will show all rows that do not match
print(bad)
if (nrow(bad) == 0){
  MD_stacked <- data_merge[,!(colnames(data_merge) %in% c("best", "worst"))]
  cat("********************************")
  cat("\n")
  cat("*  Data Files Match Perfectly  *")
  cat("\n")
  cat("********************************")
  cat("\n")
} else message ("Data Files Do Not Match! See file named bad for all mismatches.")
return(MD_stacked)
}


ratestack_function <- function(data_master, anchor_name) {
#######################################################################
#          PREPARE PI, ANCHOR, OR OTHER RATING DATA (IF YOU HAVE)
######################################################################
data_master <- as.data.frame(data_master)
anchors <- paste0(anchor_name, "_r*")
id <- data_master[,1]
Rate_data <- data_master[,grep(anchors, colnames(data_master))] # PI Ratings 
Rate_stacked <- stacker(cbind(id, Rate_data), 1)
Rate_stacked <- Rate_stacked[!is.na(Rate_stacked[,3]),] # Remove NA
colnames(Rate_stacked) <- c("id", "items", "rating")
return(Rate_stacked)
}

