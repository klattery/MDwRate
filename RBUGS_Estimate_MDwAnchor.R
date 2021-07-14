rbugs_function <- function(MD_stacked, Anchor_stacked, anchor_bounds) {
########################################################
#  1. Specify Directories and Data                     #
#    
#  MD_stacked has these first 7 variables in this order:
#   id, task, version, concept, items, dep_best, dep_worst
#   version and concept are not used and can be set to anything
#
#  Rate_stacked has these first 3 variables in this order:
#   id, items, rating
#   Rating can also be a binary choice 0 = No, 1 = Yes
########################################################
# dir_data <- "C:/Users/k.lattery/OneDrive - SKIM/Development/Conjoint/MaxDiff/Data_Both3/"
# dir_out <- "C:/Users/k.lattery/OneDrive - SKIM/Development/Conjoint/MaxDiff/Data_Both3/Ratings/"

MD_stacked <- as.data.frame(MD_stacked)
Anchor_stacked <- as.data.frame(Anchor_stacked)
table(Anchor_stacked[,3])

########################################################
#  2. Specify Type of Anchor with rate_bounds
#   Examples of 2 options:
#     anchor_bounds <- c(.5, 5.5) # 5 Pt scale with extra .5 on ends
#     anchor_bounds <- "Binary"
########################################################

# Bounds for scaling of MD utilities: MD * k
kbounds <- c(.1, 3) # Recommend: (.1, 3) for ratings in [0,10].  

# Resources >> Tools & Templates >>  Tools >> 03 Analysis >> HB_RBUGS
# source("RBUGS_Functions_v1.1.R")
# source("RBUGS_Model_MDwAnchor_v1.R")

# Prep Data
ind_all <- as.matrix(c(MD_stacked[,5], Anchor_stacked[,2]))
ind_type <- rep(2, nrow(ind_all))
ind_type[1:length(MD_stacked[,5])] <- 1
colnames(ind_all) <-"items"

indcode_list <- list()
indcode_list[[1]] <- catcode(ind_all, 1, codetype = 1) # Column 5 has MD items

make_codefiles(indcode_list) # Create: code_master, indcode, indprior 
# write.table(cbind(rownames(code_master), code_master), file = paste0(dir_out, "code_master.csv"), sep = ",", na = ".", row.names = FALSE)

data_all <- prep_file_MD(idtask = MD_stacked[,1:2], # id, task
                         best = MD_stacked[,6], 
                         worst = MD_stacked[,7],
                         ind = indcode[ind_type == 1,])

# Code ratings file same way (use code_master)
data_all$rate$ind <- indcode[ind_type == 2,]
data_all$rate$dep <- as.matrix(Anchor_stacked[,3])
data_all$rate$resp_id <- as.vector(unique(Anchor_stacked[,1]))
data_all$rate$match_betasrow <- as.vector(match(data_all$rate$resp_id, data_all$resp_id)) # subset of betas that match
data_all$rate$match_betascol <- TRUE
data_all$rate$match_id <- match(Anchor_stacked[,1], data_all$rate$resp_id)

rm(ind_all, ind_type, indcode_list)
if (is.numeric(anchor_bounds)){
  data_all$rate$type <- "Rate"
  data_all$rate$lb <- anchor_bounds[1]
  data_all$rate$ub <- anchor_bounds[2]
  if ((anchor_bounds[2] < max(data_all$rate$dep)) | (anchor_bounds[1] > min(data_all$rate$dep))){
    message("*********************************************************")
    message("\n")
    message("* WARNING!!!  Some values in your data are outside anchor_bounds")     
  }
} else {
  data_all$rate$type <- "Binary"
  data_all$rate$lb <- 0
  data_all$rate$ub <- 1
}


# Useful constants
npar <- ncol(data_all$ind) #number of parameters
nbetas <- length(data_all$resp_id) # number of betas
betas_0 <- matrix(rnorm(nbetas * npar)/100, nrow = nbetas, ncol = npar)
alpha_0 <- colMeans(betas_0) + rnorm(npar)/1000
priorcov <- indprior
ntasks2 <- nrow(Anchor_stacked)/length(unique(Anchor_stacked[,1])) # Tasks per resp in rating

pred_func <- function(data_list, pred_env){
  # MaxDiff + Ratings
  # parameter environment has betas, k, threshold
  # data_list has standard plus $rate with ind, dep, matching stuff, lb, ub
  
  betas <- pred_env$betas
  
  # Best
  V <- rowSums(data_list$ind * betas[data_list$match_id,])
  U <- exp(V) #Utility   
  tasksum <- rowsum(U, data_list$idtask_r) # Sum U each task
  esum <- (tasksum[data_list$idtask_r,]) # Map Sum to rows
  predprob1 <- U/esum
  
  # Worst
  U <- 1/U # Equivalent to exp(-V)
  tasksum <- rowsum(U, data_list$idtask_r) # Sum U each task
  esum <- (tasksum[data_list$idtask_r,]) # Map Sum to rows
  predprob2 <- U/esum  
  
  # Ratings
  betas_rate <- sweep(betas, 1, pred_env$k, "*")[data_list$rate$match_betasrow, data_list$rate$match_betascol]
  V <- rowSums(data_list$rate$ind * betas_rate[data_list$rate$match_id,])
  predrate <- data_list$rate$lb + ((data_list$rate$ub - data_list$rate$lb) *
                                     1/(1 + exp(0 - V))) 
  return(list(predprob1 = predprob1, predprob2 = predprob2, predrate = predrate))
}


LL_func_rate <- function(data_list, pred_env) { # Standard likelihood
  predlist <- pred_env$pred_func(data_list, pred_env)
  LL1 <- rowsum(log(predlist[[1]]) * data_list$dep1 * data_list$wts, data_list$match_id)
  LL2 <- rowsum(log(predlist[[2]]) * data_list$dep2 * data_list$wts, data_list$match_id)
  LL3 <- rowsum(OLS_LL(predlist[[3]], data_list$rate$dep), data_list$rate$match_id)
  LL3_match <- as.matrix(rep(0, nrow(LL1)))
  LL3_match[data_list$rate$match_betasrow,1] <- LL3
  return(as.vector(LL1 + LL2 + model_env$LLwt * LL3_match)) # Multiply rate LL because there are only a few tasks
}

LL_func_binary <- function(data_list, pred_env) { # Standard likelihood
  predlist <- pred_env$pred_func(data_list, pred_env)
  LL1 <- rowsum(log(predlist[[1]]) * data_list$dep1 * data_list$wts, data_list$match_id)
  LL2 <- rowsum(log(predlist[[2]]) * data_list$dep2 * data_list$wts, data_list$match_id)
  pred3 <- pmin(pmax(predlist[[3]], .0000001),.999999)
  LL3 <- rowsum(log(pred3) * data_list$rate$dep + 
                  log(1 - pred3) * (1- data_list$rate$dep),
                data_list$rate$match_id)
  LL3_match <- as.matrix(rep(0, nrow(LL1)))
  LL3_match[data_list$rate$match_betasrow,1] <- LL3
  return(as.vector(LL1 + LL2 + model_env$LLwt * LL3_match)) # Multiply rate LL because there are only a few tasks
}

model_pars <- list(
  # Set of parameters and prediction function
  # par_func are parameters used by pred_func
  par_func = list(
    betas = betas_0,
    k = matrix(1 + rnorm(nbetas)/100, nbetas), # Must start k inside bounds
    thresh = matrix(rnorm(nbetas)/100, nbetas)
  ),
  par_others = list(
    b_alpha =  alpha_0,
    b_cov =  priorcov,
    b_cov_IW = priorcov,  # Used to update cov as IW 
    # b ~ MVN(b_alpha, b_cov)
    
    k_alpha = 1,
    k_cov = matrix((kbounds[2] - kbounds[1])/3)^2,
    # k ~ N(k_alpha, SD = (UB - LB)/3)  
    
    #thresh_alpha = 0,
    #thresh_cov = matrix(2), # Vague Prior cov 
    # thresh ~ N(k_alpha, SD = 5)
    
    b_prop = priorcov,
    norm_1 = matrix(1),
    norm_2 = matrix(2),
    # Extra Parameters Used  
    
    LLwt = max(1, 9/ntasks2) # Weight of the LL values for ratings
  )
)

cond_betas <- list(
  method = "MH",
  par = "betas",
  constraints = NULL,
  prior_func = HB_R$Prior_MVN,
  prior_data = list(
    alpha = "b_alpha",
    cov = "b_cov"
  ),
  prop_func = HB_R$Proposal_RW,
  prop_data = list(
    cov = "b_prop", # b_cov normally
    rho = 2.4 / sqrt(npar),
    update_rho = TRUE
  ),
  accept = list(
    target = .35, # Target 
    current = .35 # Initial Value updated
  )
)

cond_b_alpha <- list(
  method = "Norm_DataCov",
  par = "b_alpha",
  prior = list(
    data = "betas",
    cov = "b_cov")
)

cond_b_cov <- list(
  method = "IW_Cov",
  par = "b_cov",
  prior = list(
    data = "betas",
    alpha = "b_alpha",
    cov = "b_cov_IW",
    df = 5)
)

cond_k <- list(
  method = "MH",
  par = "k",
  constraints = NULL,
  prior_func = HB_R$Prior_MVN,
  prior_data = list(
    alpha = "k_alpha",
    cov = "k_cov"
  ),
  prop_func = Proposal_RW_NT,
  prop_data = list(
    bound = kbounds,
    rho = 2.4,
    update_rho = TRUE
  ),
  accept = list(
    target = .45, # Target 
    current = .45 # Initial Value updated
  )
)
cond_k_alpha <- list(
  method = "NormTrunc_DataCov",
  par = "k_alpha",
  prior = list(
    data = "k",
    cov = "k_cov"),
  bound = kbounds # Truncated draws need bounds
)

cond_thresh <- list(
  method = "MH",
  par = "thresh",
  constraints = NULL,
  prior_func = HB_R$Prior_MVN,
  prior_data = list(
    alpha = "thresh_alpha",
    cov = "thresh_cov"
  ),
  prop_func = HB_R$Proposal_RW,
  prop_data = list(
    cov = "norm_1",
    rho = 2.4, 
    update_rho = TRUE
  ),
  accept = list(
    target = .45, # Target
    current = .45 # Initial 
  )
)

cond_thresh_alpha <- list(
  method = "Norm_DataCov",
  par = "thresh_alpha",
  prior = list(
    data = "thresh",
    cov = "thresh_cov")
)

cond_thresh_cov <- list(
  method = "IW_Cov",
  par = "thresh_cov",
  prior = list(
    data = "thresh",
    alpha = "thresh_alpha",
    cov = "norm_1",
    df = 5)
)

model_methods <- list(cond_betas = cond_betas,
                      cond_b_alpha = cond_b_alpha,
                      cond_b_cov = cond_b_cov,
                      cond_k = cond_k,
                      cond_k_alpha = cond_k_alpha
) #Each element describes updates of parameters

Report_MDwRate <- function(data_list){
  if (!is.null(dev.list())) dev.off()
  #dev.new()
  m <- rbind(c(0, 0.5, 0.55, 1), c(0.5, 1, 0.55, 1),
             c(0, 0.5, 0, 0.55), c(0.5, 1, 0, 0.55), c(.5,.55,.57,.66))
  split.screen(m)
  screen(1)
  par(mar = c(4, 4, 1, 1))
  screen(2)
  par(mar = c(4, 4, 1, 1))
  screen(3)
  par(mar = c(4, 4, 2, 1))
  screen(4)
  par(mar = c(4, 4, 2, 1))
  Pred_PtEst <- pred_env$pred_func(data_list, pred_env)    
  check <- cbind(data_list$dep1, Pred_PtEst[[1]])
  check2 <- cbind(data_list$dep2, Pred_PtEst[[2]])
  check3 <- cbind(data_list$rate$dep, Pred_PtEst[[3]])
  RLH1 <- exp(sum(check[,1]* log(check[,2]))/sum(check[,1]))
  RLH2 <- exp(sum(check2[,1]* log(check2[,2]))/sum(check2[,1]))
  MuP1 <- sum(check[,1]* check[,2])/sum(check[,1])
  MuP2 <- sum(check2[,1]* check2[,2])/sum(check2[,1])
  MSE <- mean((check3[,1] - check3[,2])^2)
  corr <- cor(check3[,1],check3[,2])
  
  screen(1)
  hist(check[check[,1] == 1,2], breaks = 100, col = "cadetblue", main = "Predicted Probs for Best Choices", xlab = "")
  mtext(paste0("RLH/mu ", "\n", format(RLH1, digits = 3, nsmall = 3)), side = 3, adj = .05, padj = 1.5)
  mtext(format(MuP1, digits = 3), side = 3, adj = .05, padj = 5)
  screen(2)
  hist(check2[check2[,1] == 1,2], col = "indianred1", breaks = 100, main = "Predicted Probs for Worst Choices", xlab = "", ylab = "")
  mtext(paste0("RLH/mu ", "\n", format(RLH2, digits = 3, nsmall = 3)), side = 3, adj = .05, padj = 1.5)
  mtext(format(MuP2, digits = 3), side = 3, adj = .05, padj = 5)
  screen(3)
  rate_err <- check3[,2] - check3[,1]
  hist(rate_err,  col = "orange", breaks = 100, main = "Ratings: Predicted - Observed", xlab = "")
  mtext(paste0("MSE/R ", "\n", format(MSE, digits = 3)), side = 3, adj = .05, padj = 1.5)
  mtext(format(corr, digits = 3), side = 3, adj = .05, padj = 5)
  screen(4)
  tpi <- (Sys.time() - model_env$time_beg) /model_env$iter
  tleft <- tpi * (model_env$tot_iter - model_env$iter)
  units(tleft) <- "mins"
  LLPlot <- model_env$LL_Hist
  LLPlot[,1] <- model_env$LL_Hist[,1]/1000
  LLPlot[,2] <- model_env$LL_Hist[,2]/model_env$scale_y
  ymin <- round(min(LLPlot[,2]) * 1.05, 1)
  ymax <- round(max(LLPlot[,2]) * .95, 1)
  interval <- max(round((ymax - ymin)/4,1),.5)
  ymax <- ymin + (4 * interval)
  plot(x = 0, y = 0, type = "n", main = paste0("MCMC Fit -- Time Left: ", format(tleft, digits = 2)),
       xlim = c(0,model_env$tot_iter/1000), xlab = "Iterations x 1000", 
       ylim = c(ymin, ymax), ylab = paste0("Log Likelihood / ",model_env$scale_y), axes = FALSE, cex = 0.5)
  #mtext(paste0("Time to completion: ", format(tleft, digits = 3)), side = 3)
  segments(hb_control$iter_burn/1000, ymin, hb_control$iter_burn/1000, ymax, col = "red", 
           lty = 2, lwd = 2)
  segments(0, 0, model_env$tot_iter, 0, col = "gray", lty = 1, lwd = 1)
  axis(1, at = seq(from = 0, to = model_env$tot_iter/1000, by = (model_env$tot_iter)/10000))
  axis(2, at = seq(from = ymin, to = ymax, by = interval))
  lines(LLPlot, type = 'l', col = 'blue', cex = .5)
}

Report_MDwBinary <- function(data_list){
  if (!is.null(dev.list())) dev.off()
  #dev.new()
  m <- rbind(c(0, 0.5, 0.55, 1), c(0.5, 1, 0.55, 1),
             c(0, 0.5, 0, 0.55), c(0.5, 1, 0, 0.55), c(.5,.55,.57,.66))
  split.screen(m)
  screen(1)
  par(mar = c(4, 4, 1, 1))
  screen(2)
  par(mar = c(4, 4, 1, 1))
  screen(3)
  par(mar = c(4, 4, 2, 1))
  screen(4)
  par(mar = c(4, 4, 2, 1))
  Pred_PtEst <- pred_env$pred_func(data_list, pred_env)    
  check <- cbind(data_list$dep1, Pred_PtEst[[1]])
  check2 <- cbind(data_list$dep2, Pred_PtEst[[2]])
  check3 <- cbind(data_list$rate$dep, Pred_PtEst[[3]])
  check3 <- rbind(check3, 1 - check3)
  RLH1 <- exp(sum(check[,1]* log(check[,2]))/sum(check[,1]))
  RLH2 <- exp(sum(check2[,1]* log(check2[,2]))/sum(check2[,1]))
  RLH3 <- exp(sum(check3[,1]* log(check3[,2]))/sum(check3[,1]))
  
  MuP1 <- sum(check[,1]* check[,2])/sum(check[,1])
  MuP2 <- sum(check2[,1]* check2[,2])/sum(check2[,1])
  MuP3 <- sum(check3[,1]* check3[,2])/sum(check3[,1])
  
  screen(1)
  hist(check[check[,1] == 1,2], breaks = 100, col = "cadetblue", main = "Predicted Probs for Best Choices", xlab = "")
  mtext(paste0("RLH/mu ", "\n", format(RLH1, digits = 3, nsmall = 3)), side = 3, adj = .05, padj = 1.5)
  mtext(format(MuP1, digits = 3), side = 3, adj = .05, padj = 5)
  screen(2)
  hist(check2[check2[,1] == 1,2], col = "indianred1", breaks = 100, main = "Predicted Probs for Worst Choices", xlab = "", ylab = "")
  mtext(paste0("RLH/mu ", "\n", format(RLH2, digits = 3, nsmall = 3)), side = 3, adj = .05, padj = 1.5)
  mtext(format(MuP2, digits = 3), side = 3, adj = .05, padj = 5)
  screen(3)
  hist(check3[check3[,1] == 1,2], col = "orange", breaks = 100, main = "Predicted Probs for Binary Anchor", xlab = "", ylab = "")
  mtext(paste0("RLH/mu ", "\n", format(RLH3, digits = 3, nsmall = 3)), side = 3, adj = .05, padj = 1.5)
  mtext(format(MuP3, digits = 3), side = 3, adj = .05, padj = 5)
  screen(4)
  tpi <- (Sys.time() - model_env$time_beg) /model_env$iter
  tleft <- tpi * (model_env$tot_iter - model_env$iter)
  units(tleft) <- "mins"
  LLPlot <- model_env$LL_Hist
  LLPlot[,1] <- model_env$LL_Hist[,1]/1000
  LLPlot[,2] <- model_env$LL_Hist[,2]/model_env$scale_y
  ymin <- round(min(LLPlot[,2]) * 1.05, 1)
  ymax <- round(max(LLPlot[,2]) * .95, 1)
  interval <- max(round((ymax - ymin)/4,1),.5)
  ymax <- ymin + (4 * interval)
  plot(x = 0, y = 0, type = "n", main = paste0("MCMC Fit -- Time Left: ", format(tleft, digits = 2)),
       xlim = c(0,model_env$tot_iter/1000), xlab = "Iterations x 1000", 
       ylim = c(ymin, ymax), ylab = paste0("Log Likelihood / ",model_env$scale_y), axes = FALSE, cex = 0.5)
  #mtext(paste0("Time to completion: ", format(tleft, digits = 3)), side = 3)
  segments(hb_control$iter_burn/1000, ymin, hb_control$iter_burn/1000, ymax, col = "red", 
           lty = 2, lwd = 2)
  segments(0, 0, model_env$tot_iter, 0, col = "gray", lty = 1, lwd = 1)
  axis(1, at = seq(from = 0, to = model_env$tot_iter/1000, by = (model_env$tot_iter)/10000))
  axis(2, at = seq(from = ymin, to = ymax, by = interval))
  lines(LLPlot, type = 'l', col = 'blue', cex = .5)
}

if (data_all$rate$type == "Binary"){
  data_all$LL_func <- LL_func_binary
  data_all$report_function <- Report_MDwBinary
} else {
  data_all$LL_func <- LL_func_rate
  data_all$report_function <- Report_MDwRate
}

message("*********************************************************")
message("\n")
message("* Data Fusion Model for MaxDiff and Anchor Loaded")
message("\n")
message(paste0("* Anchor Type is ", data_all$rate$type, " with range [", data_all$rate$lb, ", ", data_all$rate$ub, "]"))
message("\n")
message(paste0("* Created list data_all with total of ", nbetas, " respondents in MaxDiff."))
message("\n")
message(paste0("* We found a mean of ", ntasks2, " anchors per respondent."))
message("\n")
message(paste0("* Anchors weighted at ", model_pars$par_others$LLwt, " in model_pars$par_others$LLwt"))
message("\n")
message("*********************************************************")





# Creates model_pars, model_methods, pred_func, LL_func

# Set parameters for run
hb_control <- list(
  iter_burn = 2000,
  iter_sample = 1000,
  thin = 1,
  report_interval = 50,
  report_function = data_all$report_function,
  track = c(
    "betas",
    "k",
    "b_alpha",
    "k_alpha"
  )
)

#model_pars$par_others$LLwt <- 1 to change weight of ratings
HB_EST(data_list = data_all,
       pred_func = pred_func,
       LL_func = data_all$LL_func,
       model_pars = model_pars, # LL function and subset of parameters used in that function 
       model_methods = model_methods,
       hb_control = hb_control,
       ContinueOld = FALSE)

#############################################

# Check Convergence
close.screen(all = TRUE) # Close split screen
draws_alpha <- do.call(rbind,model_env$track$b_alpha)
GR1 <- GR(list(draws_alpha, draws_alpha)) # Convergence of MD

draws_ratep <- as.matrix(unlist(model_env$track$k_alpha))
GR2 <- GR(list(draws_ratep, draws_ratep)) # Convergence of Ratings
plot(draws_ratep[,1], type = "l", ylab = "Scaling Alpha", xlab ="Iteration")
mean(draws_ratep)

# Compute Point estimates
betas <- Reduce("+",model_env$track$betas)/length(model_env$track$betas) #coded
betas_r <- betas %*% t(code_master) #back code
k <- (Reduce("+",model_env$track$k)/length(model_env$track$k))
hist(k, breaks = 100, main = "Respondent Utility: Scaling Factor K")

# Compute Predicted Ratings for convenience
betas_rate <- sweep(betas_r, 1, k, "*")
rating <- data_all$rate$lb + (data_all$rate$ub - data_all$rate$lb)/(1 + exp(-betas_rate))

# Export
export_utils <- data.frame(cbind(id = data_all$resp_id, betas_r))
export_utils$k <- k
export_ratings <- data.frame(cbind(id = data_all$resp_id, rating))
return(export_utils)
# write.table(export_ratings, file = paste0(dir_out, "export_ratings.csv"), sep = ",", na = ".", row.names = FALSE)
}
