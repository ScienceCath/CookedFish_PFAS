#' @title pred_interval_esmeans
#' @description Function to get prediction intervals (crediblity intervals) from esmeans objects (metafor)
#' @param model rma.mv object
#' @param esmeans result from emmeans::emmeans object
#' @param ... other arguments passed to function
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export
pred_interval_esmeans <- function(model, mm, ...){
  
  tmp <- summary(mm)
  test.stat <- qt(0.975, tmp$df)
  
  if(length(model$tau2) <= 1){
    sigmas <- sum(model$sigma2)
    PI <- test.stat * sqrt(tmp$SE^2 + sigmas)
  } else {
    sigmas <- sum(model$sigma2)
    taus   <- model$tau2
    w <- model$g.levels.k
    
    if(pred == "1"){
      tau <- weighted_var(taus, weights = w)
      PI <- test.stat * sqrt(tmp$SE^2 + sigmas + tau)
      
    } else {
      PI <- test.stat * sqrt(tmp$SE^2 + sigmas + taus)
    }
  }
  
  tmp$lower.PI <- tmp$emmean - PI
  tmp$upper.PI <- tmp$emmean + PI
  
  return(tmp)
}

#' @title marginalised_means
#' @description Function to to get marginalised means from met-regression models with single or multiple moderator variables that are both continuous or categorical.
#' @param model rma.mv object
#' @param data data frame used to fit rma.mv model
#' @param pred predictor variable of interest that one wants marginalised means for.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @example \dontrun{
#'warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat,                               control=list(optimizer="optim", optmethod="Nelder-Mead"))
#'   overall <- marginalised_means(model, data = warm_dat)
#' across_trait <- marginalised_means(model, data = warm_dat, pred = "trait.type")
#' across_trait_by_degree_diff <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif")
#' across_trait_by_degree_diff_at_treat_end_days10 <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = 10), by = "deg_dif")
#' across_trait_by_degree_diff_at_treat_end_days10And50 <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "deg_dif")
#' across_trait_by_treat_end_days10And50 <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days")
#' across_trait_by_treat_end_days10And50_ordinaryMM <- marginalised_means(model, data = warm_dat, pred = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", weights = "prop")
#' }
#' @export
#'
#'
# We will ned to make sure people use "1" pr "moderator_names"
marginalised_means <- function(model, data, pred = "1", by = NULL, at = NULL, ...){
  model$data <- data
  
  grid <- emmeans::qdrg(object = model, at = at)
  mm <- emmeans::emmeans(grid, specs = pred, df = model$df, by = by, ...)
  mm_pi <- pred_interval_esmeans(model, mm, pred = pred)
  
  
  if(is.null(by)){
    mod_table <- tibble::tibble(name = mm_pi[,1], estimate = mm_pi[,"emmean"], lowerCL = mm_pi[,"lower.CL"], upperCL = mm_pi[,"upper.CL"], lowerPI = mm_pi[,"lower.PI"], upperPI = mm_pi[,"upper.PI"])
    
  } else{
    mod_table <- tibble::tibble(name = mm_pi[,1], mod = mm_pi[,2], estimate = mm_pi[,"emmean"], lowerCL = mm_pi[,"lower.CL"], upperCL = mm_pi[,"upper.CL"], lowerPI = mm_pi[,"lower.PI"], upperPI = mm_pi[,"upper.PI"])
    
  }
  
  output <- list(mod_table = mod_table,
                 data = data)
  
  class(output) <- "orchard"
  
  return(output)
}

