LogLikelihood <- function(init, num_nmodel, vec_modelname, envs,
                          lst_option, lst_psdparams, # model information
                          psd_preexpr, psd_timedepexpr, # parsed expression
                          vec_parameter, vec_variable, # variable name
                          vec_innername){ # inner variable name (nullvalue, weight, time1, time2, status)
  ret <- 0
  for(i in 1:num_nmodel){
    if(vec_modelname[i] == "cox"){
      
      ret <- ret + LogCoxLikelihood(init=init, envs=envs[[i]], lst_option=lst_option[[i]],
                             lst_psdparams=lst_psdparams[[i]], psd_preexpr=psd_preexpr[[i]],
                             psd_timedepexpr=psd_timedepexpr[[i]],
                             vec_parameter=vec_parameter, vec_variable=vec_variable,
                             vec_innername=vec_innername)
      
    } else {

      if(is.environment(envs[[i]][[1]])){ # no random effect

        ret <- ret +
          InnerLogLikelihood(init=init, modelname=vec_modelname[i], envs=envs[[i]],
                             lst_option=lst_option[[i]], lst_psdparams=lst_psdparams[[i]],
                             psd_preexpr=psd_preexpr[[i]],
                             vec_parameter=vec_parameter, vec_variable=vec_variable,
                             vec_innername=vec_innername)
        
      } else { # with random effect
        stop("reached unreachable code area in LogLikelihood function")
      }
    }
  }
  
  return(ret)
}

LogCoxLikelihood <- function(init, envs, lst_option, lst_psdparams, psd_preexpr,
                             psd_timedepexpr, vec_parameter, vec_variable, vec_innername){

  ret <- 0
  ties <- ""
  
  if(is.null(lst_option)){
    ties <- "efron"
  } else {
    ties <- lst_option$"ties"
    if(is.null(ties))
      ties <- "efron"
  }

  if(!ties %in% c("efron","breslow","average","discrete")){
    stop("Invalid tie specification in cox regression")
  }
  
  for(strata in 1:length(envs)){

    ## Assign parameter value
    for(i in 1:length(vec_parameter))
      assign(vec_parameter[i], init[i], envir=envs[[strata]])

    time1 <- get(vec_innername[2], envir=envs[[strata]])
    time2 <- get(vec_innername[3], envir=envs[[strata]])
    status <- get(vec_innername[4], envir=envs[[strata]])

    nsubject <- length(status)
    result <- 0
    riskset <- numeric(nsubject)
    pretime <- time2[1]
    duringtie <- FALSE
    tiebegin <- 0
    
    if(!exists(vec_innername[1], mode="numeric", envir=envs[[strata]])){
      weight <- rep(1, length(status))
    } else {
      if(ties=="discrete" || ties=="average")
        stop("weight is not supported for ties=\"", ties,"\" specification")
      weight <- get(vec_innername[1], envir=envs[[strata]])
    }
    
    if(length(psd_timedepexpr) > 0){
      
      ## setting inner variable "T" to the first time point
      ## calculate riskset first for all datapoint
      ## assign(vec_innername[6], rep(time2[1], length(time2)), envir=envs[[strata]])
      assign(vec_innername[5], time2[1], envir=envs[[strata]])
      eval(psd_timedepexpr, envir=envs[[strata]])
      hazard <- eval(lst_psdparams[[1]], envir=envs[[strata]])
      whazard <- hazard*weight
      phazard <- numeric(nsubject)
      phazard[1] <- hazard[1]^weight[1]
      riskset[1] <- sum(as.numeric((time1 < time2[1]) & (time2[1] <= time2))*whazard)
      
    } else {
      hazard <- eval(lst_psdparams[[1]], envir=envs[[strata]])
      whazard <- hazard*weight
      phazard <- hazard^weight
      riskset[1] <- sum(as.numeric((time1 < time2[1]) & (time2[1] <= time2))*whazard)
    }
    
    ## Main calculation loop begins here
    for(i in 2:nsubject){

      if(length(psd_timedepexpr) > 0){
        assign(vec_innername[5], time2[i], envir=envs[[strata]])    
        eval(psd_timedepexpr, envir=envs[[strata]])
        hazard <- eval(lst_psdparams[[1]], envir=envs[[strata]])
        whazard <- hazard*weight
        phazard[i] <- hazard[i]^weight[i]
        riskset[i] <- sum(as.numeric((time1 < time2[i]) & (time2[i] <= time2))*whazard)
      } else {
        riskset[i] <- sum(as.numeric((time1 < time2[i]) & (time2[i] <= time2))*whazard)
      }
      
      ## in case of tie
      if(time2[i] == pretime){
        
        pretime <-  time2[i]
        
        if(!duringtie){
          duringtie <- TRUE
          tiebegin <- i-1
        }
        
        if(i==nsubject){
          ## Assume events occure before censoring
          phazard[tiebegin:nsubject] <- phazard[tiebegin:nsubject][order(status[tiebegin:nsubject], decreasing=TRUE)]
          status[tiebegin:nsubject] <- status[tiebegin:nsubject][order(status[tiebegin:nsubject], decreasing=TRUE)]
          
          tieevent <- sum(status[tiebegin:nsubject]) # number of event in tie
          
          ## Efron approximation
          if(ties=="efron"){
            
            if(tieevent > 0){
              tiehazard <- sum(whazard[tiebegin:(tiebegin+tieevent-1)])
              
              ## Calculate riskset for efron: tie is between (tiebegin ... nsubject)
              for(j in (tieevent-1):0){
                riskset[tiebegin+j] <- (riskset[tiebegin] - j/tieevent*tiehazard)^(sum(weight[tiebegin:(tiebegin+tieevent-1)])/tieevent)
              }
            }
            
          } else if(ties=="breslow"){ # do nothing
            riskset[tiebegin:(tiebegin+tieevent-1)] <- riskset[tiebegin:(tiebegin+tieevent-1)]^weight[tiebegin:(tiebegin+tieevent-1)]
            
            ## Average partial likelihood 
          } else if(ties=="average"){
            
            eventhazard <- 0
            
            if(tieevent > 0){
              eventhazard <- sum(whazard[tiebegin:(tiebegin+tieevent-1)])
              
              lik <- function(t){
                res <- 1
                for(k in tiebegin:(tiebegin+tieevent-1)){
                  res <- res*(1 - exp(-(hazard[k]*t)/(riskset[k]-eventhazard)))^weight[k]
                }
                return(res*exp(-t))
              }
              
              ## put all tie combination log partial likelihood into the last tie position
              phazard[tiebegin] <- integrate(lik, 0, Inf)$value
              status[tiebegin] <- 1 # regard as event
              riskset[tiebegin] <- 1
              status[(tiebegin+1):nsubject] <- 0
            }
            
            ## "discrete" in SAS or "exact" in R
          } else if(ties=="discrete"){
            
            if(tieevent > 0){
              phazard[tiebegin] <- prod(phazard[tiebegin:(tiebegin+tieevent-1)])/.Call(Rf_select, tieevent, nsubject-tiebegin+1, phazard[tiebegin:nsubject])
              status[tiebegin] <- 1 # regard as event
              riskset[tiebegin] <- 1
              status[(tiebegin+1):(i-1)] <- 0
            }
          }
          
        } else { # i=nsubject ends
          ## in the tie
          next # skip until end of tie
        }
      } else { # not in case of tie
        pretime <- time2[i]
        
        ## when tie ends (tie is between tiebegin and (i-1))
        ##    tiebegin ... eventpos       (i-1)   i
        ## no    tie   ...   tie    ...    tie    no  (N of tie is (tiebegin-i))
        ## ...   evt   ...   evt    cen    cen    ...
        if(duringtie){
          
          duringtie <- FALSE
          
          ## Assume event occurs before censoring
          phazard[tiebegin:(i-1)] <- phazard[tiebegin:(i-1)][order(status[tiebegin:(i-1)], decreasing=TRUE)]
          status[tiebegin:(i-1)] <- status[tiebegin:(i-1)][order(status[tiebegin:(i-1)], decreasing=TRUE)]
          
          tieevent <- sum(status[tiebegin:(i-1)]) # number of event in tie
          
          ## Efron approximation
          if(ties=="efron"){
            
            if(tieevent > 0){
              tiehazard <- sum(whazard[tiebegin:(tiebegin+tieevent-1)])
              
              for(j in (tieevent-1):0){
                riskset[tiebegin+j] <- (riskset[tiebegin] - j/tieevent*tiehazard)^(sum(weight[tiebegin:(tiebegin+tieevent-1)])/tieevent)
              }
            }
            
          } else if(ties=="breslow"){ # do nothing
            riskset[tiebegin:(tiebegin+tieevent-1)] <- riskset[tiebegin:(tiebegin+tieevent-1)]^weight[tiebegin:(tiebegin+tieevent-1)]
          } else if(ties=="average"){
            
            eventhazard <- 0
            
            if(tieevent > 0){
              eventhazard <- sum(whazard[tiebegin:(tiebegin+tieevent-1)])
              
              lik <- function(t){
                res <- 1
                for(k in tiebegin:(tiebegin+tieevent-1)){
                  res <- res*(1 - exp(-hazard[k]*t/(riskset[k]-eventhazard)))^weight[k]
                }
                return(res*exp(-t))
              }
              
              phazard[tiebegin] <- integrate(lik, 0, Inf)$value
              status[tiebegin] <- 1 # regard as event
              riskset[tiebegin] <- 1
              status[(tiebegin+1):(i-1)] <- 0
            }
            
          } else if(ties=="discrete"){
            
            ## when tie ends (tie is between tiebegin and (i-1))
            ##    tiebegin ... eventpos       (i-1)   i
            ## no    tie   ...   tie    ...    tie    no  (N of tie is (tiebegin-i))
            ## ...   evt   ...   evt    cen    cen    ...
            
            if(tieevent > 0){
              phazard[tiebegin] <- prod(phazard[tiebegin:(tiebegin+tieevent-1)])/.Call("Rf_select", tieevent, nsubject-tiebegin+1, phazard[tiebegin:nsubject])
              status[tiebegin] <- 1 # regard as event
              riskset[tiebegin] <- 1
              status[(tiebegin+1):(i-1)] <- 0
            }
          }
        } # end of during tie
        ## no tie...nothing to do
      }
    }

    ret <- ret - sum(log(phazard/riskset)*(as.numeric(status==1)))
  }
  return(ret)
}

InnerLogLikelihood <- function(init, modelname, envs, lst_option, lst_psdparams, psd_preexpr,
                               vec_parameter, vec_variable, vec_innername){
  env <- envs[[1]]
  
  for(i in 1:length(vec_parameter))
    assign(vec_parameter[i], init[i], envir=env)

  if(length(psd_preexpr) > 0)
    eval(psd_preexpr, envir=env)

  if(modelname=="general"){
    if(exists(vec_innername[1], envir=env)){
      return(-sum(eval(lst_psdparams[[1]], envir=env)*get(vec_innername[1], envir=env)))
    } else {
      return(-sum(eval(lst_psdparams[[1]], envir=env)))
    }
  } else {
    func <- get(paste("d", modelname, sep=""), mode="function",
                envir=loadNamespace("stats"))
    if(exists(vec_innername[1], envir=env)){

      # results of eval are stored as list with additional option
      return(-sum(do.call(func, c(lapply(as.expression(lst_psdparams), eval, envir=env),
                                  log=TRUE))*get(vec_innername[1], envir=envs[[1]])))
    } else {

      return(-sum(do.call(func, c(lapply(as.expression(lst_psdparams), eval, envir=env),
                                  log=TRUE))))
    }
  }
}


## obtain the first environment from "envs" list
GetFirstEnvironment <- function(envs){
  while(!is.environment(envs))
    envs <- "[["(envs, 1)
  return(envs)
}

AssignNumericVector <- function(vec_value, env){
  vec_varname <- names(vec_value)
  if(is.null(vec_varname)){
    warning("Names are not set in AssignNumericVector")
  } else {
    for(i in 1:length(vec_varname))
      assign(vec_varname[i], vec_value[i], envir=env)
  }
}

CopyVariableBetweenEnvironment <- function(vec_varname, dstenv, srcenv, index, check_length){
  if(check_length){
    n <- length(index)
    for(i in 1:length(vec_varname)){
      tmp <- get(vec_varname[i], envir=srcenv)
      if(length(tmp) != n)
        warning("Length of variable \"", vec_varname[i], "\" are different from others")
      dim(tmp) <- NULL
      assign(vec_varname[i], tmp[index], dstenv)
    }
  } else {
    for(i in 1:length(vec_varname)){
      assign(vec_varname[i], get(vec_varname[i], envir=srcenv)[index], dstenv)
    }
  }
}

# Make equations which solve dependency
SolveDependence <- function(vec_depvar, lst_eqnassigned, lst_eqndepend, vec_eqns){
  unsolved <- vec_depvar
  resfml <- character(0) # result formula
  neqns <- length(vec_eqns)
  flag <- TRUE
  while(length(unsolved) > 0){
    target <- unsolved[1]
    flag <- TRUE
    for(i in neqns:1){
      if(target %in% lst_eqnassigned[[i]]){ # find
        flag <- FALSE
        resfml <- c(vec_eqns[i], resfml)
        unsolved <- unsolved[unsolved != target] # remove target
        target <- unsolved[1]
        if(length(lst_eqndepend[[i]])==1 && nchar(lst_eqndepend[[i]]) > 0){ # solve further dependency
          unsolved <- unsolved[unsolved != lst_eqndepend[[i]]]
          resfml <- c(SolveDependence(lst_eqndepend[[i]], lst_eqnassigned, lst_eqndepend, vec_eqns), resfml)
        } else if(length(lst_eqndepend[[i]]) > 1){
          unsolved <- RemoveVariableName(unsolved, lst_eqndepend[[i]])
          resfml <- c(SolveDependence(lst_eqndepend[[i]], lst_eqnassigned, lst_eqndepend, vec_eqns), resfml)
        }
        break
      } # find end
    } # eqns for loop end
    if(flag){
      stop(paste("Dependency cannot be solved for", target))
    }
  }
  return(resfml)
}

## Replace formula in language object
## eg., mu <- a + b*x; y = log(mu) --> y = log(a + b*x)
InsertFormula <- function(psd_target, vec_depfml){

  for(i in length(vec_depfml):1){
    if(length(grep("<-", vec_depfml[i])) > 0){
      chr_var <- trim(strsplit(vec_depfml[i], "<-")[[1]][1])
      chr_fml <- trim(strsplit(vec_depfml[i], "<-")[[1]][2])
    } else if(length(grep("=", vec_depfml[i])) > 0){
      chr_var <- trim(strsplit(vec_depfml[i], "=")[[1]][1])
      chr_fml <- trim(strsplit(vec_depfml[i], "=")[[1]][2])
    } else if(length(grep("->", vec_depfml[i])) > 0){
      chr_var <- trim(strsplit(vec_depfml[i], "=")[[1]][2])
      chr_fml <- trim(strsplit(vec_depfml[i], "=")[[1]][1])
    } else {
      stop("vec_depfml must be assigning sentence in InsertFormula function")
    }
    
    psd_fml <- ParseLine(chr_fml)
  
    psd_target[[1]] <- InnerInsertFormula(psd_target[[1]], chr_var, psd_fml)
  }
  return(psd_target)
}

InnerInsertFormula <- function(psd_target, chr_var, psd_fml){
  if(is.symbol(psd_target)){
    if(as.character(psd_target) == chr_var){
      return(psd_fml)
    } else {
      return(psd_target)
    }
  }
  for(i in 2:length(psd_target)){
    if(is.symbol(psd_target[[i]])){
      if(as.character(psd_target[[i]]) == chr_var){
        psd_target[[i]] <- psd_fml
      }
    } else {
      psd_target[[i]] <- InnerInsertFormula(psd_target[[i]], chr_var, psd_fml)
    }
  }
  return(psd_target)
}

## remove some variable names from variable list
RemoveVariableName <- function(varlist, remove){
  if(length(remove) == 0)
    return(varlist);
  flag <- rep(TRUE, length(varlist))
  for(i in 1:length(remove)){
    for(j in 1:length(varlist)){
      if(varlist[j]==remove[i]){
        flag[j]=FALSE
      }
    }
  }
  return(varlist[flag])
}

GetOptionString <- function(separator, src){

  ntotal <- nchar(src)
  nsep <- nchar(separator)
  range <- ntotal - nsep

  if(range > 0){
    for(i in 1:range){
      if(identical(separator, substr(src, i, i + nsep -1)))
        return(c(substr(src, 1, i - 1), substr(src, i + nsep, ntotal)))
    }
  }
  return(c("", ""))
}


## Obtain option list
GetOptions <- function(modelstr){
  ret <- list(NULL)
  name <- character(0)
  options <- strsplit(modelstr, "/")[[1]][2]
  options <- strsplit(options, ",")[[1]]
  if(is.na(options[[1]])){
    return(ret)
  } else {
    for(i in 1:length(options)){
      tmp <- GetOptionString("=", options[i])
      ret[[i]] <- tmp[2]
      name <- c(name, tmp[1])
    }
    names(ret) <- name
    return(ret)
  }    
}

## Print Information depending on verbatim
CatVerbatim <- function(level, verbatim, ...){
  content <- list(...)
  if(verbatim >= level){
    do.call(cat, content)
  }
}

## Print Information depending on verbatim
PrintVerbatim <- function(level, verbatim, variable, ...){
  content <- list(...)
  if(verbatim >= level){
    do.call(cat, content)
    print(variable)
  }
}

## return variable list (variable, pre-parameter)
## called from ClassifyParameter
ClassifyVariable <- function(vec_varlist, envir){
  res <- list(character(0), character(0))
  for(i in 1:length(vec_varlist)){
    if(exists(vec_varlist[i], mode="numeric", envir=envir))
      res[[1]] <- c(res[[1]], vec_varlist[i])
    else
      res[[2]] <- c(res[[2]], vec_varlist[i])
  }
  return(res)
}

## return variable list (parameter, variable)
## call ClassifyVariable
ClassifyParameter <- function(vec_varlist, envir, vec_remove=""){
  tmp <- ClassifyVariable(vec_varlist, envir)
  res <- list(character(0), character(0))
  res[[2]] <- tmp[[1]]

  index <- rep(TRUE, length(tmp[[2]]))

  if(length(tmp[[2]]) == 0){
    return(res)
  }
  
  ## constants
  for(i in 1:length(index)){
    ## Built-in constants
    if(tmp[[2]][i] %in% c("LETTERS", "letters", "month.abb",
                           "month.name", "pi", vec_remove)){
      index[i] <- FALSE
      next
    }
       
       tryCatch({
         as.numeric(tmp[[2]][i])
         index[i] <- FALSE
       },
                warning=function(e){},
                error=function(e){})
  }
  res[[1]] <- tmp[[2]][index]
  return(res)
}

## remove all other variables not included in varpool
LimitVarlist <- function(vec_varlist, vec_varpool){
  if(length(vec_varpool)==0){
    vec_varpool=""
  }
  res <- character(0)
  for(i in 1:length(vec_varpool)){
    if(vec_varpool[i] %in% vec_varlist){
      res <- c(res, vec_varpool[i])
    }
  }
  return(res)
}

## Obtain parameter position as integer subset list
GetParamPosition <- function(param, paramlist){
  sapply(param, function(x){
    for(i in 1:length(paramlist)){
      if(x == paramlist[i])
        return(i)
    }
  }, USE.NAMES=FALSE)
}

## Make epifit result object from optim function
MakeResultFromOptim <- function(result, ans, hessian, nulllik){
  result$coefficients <- ans$par
  result$loglik <- c(-nulllik, -ans$value)
  result$var <- ginv(hessian)
  result$iter <- ans$counts
  if(ans$convergence==0) result$convergence <- 0
  else if(ans$convergence==1) result$convergence <- 4
  else if(ans$convergence==10) result$convergence <- 6
  else if(ans$convergence==51) result$convergence <- 7
  else if(ans$convergence==52) result$convergence <- 8
  result$wald.test <- t(ans$par)%*%(ans$hessian)%*%(ans$par)
  return(result)
}

## Make epifit result object from nlm function
MakeResultFromNlm <- function(result, ans, hessian, nulllik){
  result$coefficients <- ans$estimate
  result$loglik <- c(-nulllik, -ans$minimum)
  result$var <- ginv(hessian)
  result$iter <- ans$iterations
  result$convergence <- ans$code
  result$wald.test <- t(ans$estimate)%*%(ans$hessian)%*%(ans$estimate)
  return(result)
}

## Call InnerListVariable for multiple expression
ListVariable <- function(expression){
  res <- list(NULL)
  if(is.call(expression) || is.name(expression)){
    res <- InnerListVariable(expression)
  } else if(is.list(expression)){
    if(length(expression) == 1){
      res[[1]] <- ListVariable(expression[[1]])
    } else {
      res <- lapply(expression, ListVariable)
    }
  } else if(is.expression(expression)){
    if(length(expression) == 1){
      res[[1]] <- ListVariable(expression[[1]])
    } else {
      res <- lapply(expression, ListVariable)
    }
  } else if(is.character(expression)){
    if(length(expression) == 1){
      res <- ListVariable(parse(text=expression))
    } else {
      res <- lapply(expression, ListVariable)
    }
  }
  return(res)
}

## Supported operator
getSupportedOperator <- function(){
  return(c("-","+","*","/","^","<",">","==",">=","<=", "&","|","(",")",
           "abs","acos","acosh", "as.integer","as.numeric","asin","asinh","atan",
           "atanh","cos","cosh", "digamma","exp","expm1","factorial","floor",
           "gamma","ifelse","lgamma", "lfactorial","log","log10","log1p","log2",
           "logb","pmax","pmax.int", "pmin","pmin.int",
           "sin","sinh", "sum", "tan","tanh","trigamma",
           "[", "length", "print", "cat", "rep")) ## mainly for displaying inner result
}

## Checking supported operator is included
## List assigned and used variables in 
InnerListVariable <- function(tree){
  res <- list(character(0), character(0))
  if(length(tree) == 1){
    res[[2]] <- c(res[[2]], as.character(tree))
    return(res)
  }
  op <- as.character(tree[[1]])
  if(op == "<-" || op == "="){
    res[[1]] <- c(res[[1]], as.character(tree[[2]]))
    tmp <- InnerListVariable(tree[[3]])
    res[[1]] <- c(res[[1]], tmp[[1]])
    res[[2]] <- c(res[[2]], tmp[[2]])
  } else if(op == "->"){
    res[[1]] <- c(res[[1]], as.character(tree[[3]]))
    tmp <- InnerListVariable(tree[[2]])
    res[[1]] <- c(res[[1]], tmp[[1]])
    res[[2]] <- c(res[[2]], tmp[[2]])
  } else if(!op %in% getSupportedOperator()){
    stop(paste("unsupported operator is found:", as.character(tree[[1]])), sep=" ")
  } else {
    for(j in 2:length(tree)){
      tmp <- InnerListVariable(tree[[j]])
      res[[1]] <- c(res[[1]], tmp[[1]])
      res[[2]] <- c(res[[2]], tmp[[2]])
    }
  }
  return(res)
}

## Parse one line text
ParseLine <- function(expression){
  return(parse(text=expression)[[1]])
}

## Remove spaces
trim <- function(char){
  ret <- sub("^[[:blank:]]+", "", char)
  return(sub("[[:blank:]]+$", "", ret))
}

## Divide R expressions into each line
DivideExpression <- function(string){
  fml_ret <- character(0)
  for(i in 1:length(string)){
    tmpeq <- strsplit(string[i], "\n")[[1]]
    for(j in 1:length(tmpeq)){
      fml_ret <- c(fml_ret, trim(strsplit(tmpeq, ";")[[j]]))
    }
  }
  return(fml_ret[rep("",length(fml_ret))!=fml_ret])
}

mumul <- function(vec){
  return(.Call(Rf_mumul, vec))
  #if(length(vec) < 2)
  #  return(vec[1])
  #ret <- vec[1]
  #for(i in 2:length(vec))
  #  ret <- ret*vec[i]
  #return(ret)
}

# make combination of initial parameter candidate values
makeInitVector <- function(...){
  param <- list(...) # param values
  nparam <- length(param) # number of parameters
  name <- names(param) # names of parameters
  ninit <- numeric(nparam) # number of init param candidate
  ncomb <- 1 # number of combinations
  
  for(i in 1:nparam){
    if(!is.vector(param[[i]]))
      stop("Initial parameter specification must be numeric vector")
    ninit[i] <- length(param[[i]])
    ncomb <- ncomb * length(param[[i]])
  }

  mat_ret <- matrix(0, nparam, ncomb)
  
  if(nparam == 1)
    return(param[[1]])

  dim(param[[1]]) <- NULL
  mat_ret[1,] <- kronecker(param[[1]], rep(1, mumul(ninit[-1])))
  
  if(nparam > 2){
    for(i in 2:(nparam-1)){
      dim(param[[i]]) <- NULL
      mat_ret[i,] <- rep(kronecker(param[[i]], rep(1, mumul(ninit[(i+1):nparam]))), mumul(ninit[1:(i-1)]))
    }
  }
  
  mat_ret[nparam,] <- rep(param[[nparam]], mumul(ninit[1:(nparam-1)]))
  lst_ret <- list()

  for(i in 1:ncomb){
    lst_ret[[i]] <- mat_ret[,i]
    names(lst_ret[[i]]) <- name
  }
  return(lst_ret)
}

mkSeqVariable <- function(varname="", start=0, end=0, step=1, dstenvir=NULL, srcenvir=NULL, func=mean, toLowerCategory=TRUE, defaultval=0){
  if(!is.environment(srcenvir))
    stop("environment must be specified in srcenvir argument in mkSeqVariable")
  if(!is.environment(dstenvir))
    stop("environment must be specified in dstenvir argument in mkSeqVariable")
  
  allvars <- ls(envir=srcenvir)
  idx <- grep(varname, allvars)
  
  if(length(idx) == 0){
    stop("No variables found")
  }
  
  varlist <- allvars[idx]

  if(start >= end){
    stop("start must be specified with smaller value than end")
  }
  
  resvar <- .Call(Rf_scanVarname, varlist)

  if(is.na(resvar[[2]][1]))
    stop(resvar[[1]])

  point <- resvar[[2]][order(resvar[[2]])]

  pos <- 1 # current position in pos
  resvarname <- character(length(point))
  nsubject <- length(get(paste(resvar[[1]], point[pos], sep=""), envir=srcenvir))
  
  while(point[pos] < start && pos < length(point)){
    pos <- pos + 1
  }
  
  for(i in seq(start, end, by=step)){
    varname <- paste(resvar[[1]], i, sep="")    
    resvarname[i] <- varname
    
    if(i <= point[pos]){
      npos <- 0 # number of position included in interval [i, i+step)
      while(point[pos + npos] < i + step && pos < length(point))
        npos <-  npos + 1
      if(npos > 0){
        pointdata <- matrix(0, nsubject, npos)
        value <- numeric(nsubject)
        for(j in 1:npos) # obtain all data between interval [i, i+step)
          pointdata[,j] <- get(paste(resvar[[1]], point[pos + j - 1], sep=""), envir=srcenvir)

        for(j in 1:nsubject)
          value[j] <- func(pointdata[j,]) # take mean or some functions for each point data
        assign(varname, value, envir=dstenvir)
      }

    } else {
      assign(paste(resvar[[1]], i, sep=""), rep(defaultval, nsubject), envir=dstenvir)
    }
  }
  return(resvarname)
}

# time...time point, timename...name of time variable, param...time parameter value
# paramname...name of parameter, psd_func...parsed function of weight
getVariableWeight <- function(time, timename, param, paramname, psd_func, envir){
  for(i in 1:length(param))
    assign(paramname[i], param[i], envir=envir)
  assign(timename, time)
  return(eval(psd_func, envir=envir))
}

calcWeightedVariable <- function(param, paramname, varlist=NULL, time, timename, psd_func, timepoint, envir=NULL){
  pointweight <- getVariableWeight(timepoint, timename, param, paramname, psd_func, envir)
  totalweight <- sum(pointweight)
  ntimepoint <- length(timepoint)
  nsubject <- length(get(varlist[1], envir=envir))
  pos <- findInterval(time, timepoint)
  
  ret <- numeric(nsubject)

  if(pos < 1)
    return(ret)
  for(i in 1:pos)
    ret <- ret + pointweight[i]*get(varlist[i], envir=envir)*ntimepoint/totalweight
  return(ret)
}

cumtdc <- function(time, value, varname, data1, delay, na.rm, default, epstime, options){

  if(is.na(time) || is.na(value))
    return(data1)

  varlist <- names(data1)
  startpos <- grep("tstart", varlist)
  stoppos <- grep("tstop", varlist)
  varpos <- grep(varname, varlist)
  time <- time + delay
  
  if(length(startpos) == 0 || length(stoppos) == 0)
    stop("neither a tstop argument nor an initial event argument was found")

  data1 <- data1[order(data1[, startpos]),]
  
  tstart <- data1[,startpos]
  tstop <- data1[,stoppos]
  pos <- (tstart <= time) & (time < tstop)

  if(sum(pos) > 1)
    stop("Dupicate data for one id")

  posnum <- sum((1:nrow(data1))*as.numeric(pos))

  if(sum(pos) == 1){ # interval data exist

    if(data1[posnum, startpos] == time){ # no need to divide
      
      if(length(varpos) == 0){ # no pre-value

        if(posnum > 1)
          data1[1:(posnum-1), varname] <- default
        
        data1[posnum:nrow(data1), varname] <- value
        
      } else { # pre-value exist

        data1[posnum:nrow(data1), varname] <- data1[posnum:nrow(data1), varname] + value
      }

    } else { # need to divide interval
      if(length(varpos) == 0){ # no pre-value
        #                         tstart tstop cumcov
        # pre-range(posnum - 1)      *     *
        # target-range(posnum)       *   time  value-added
        # NEW DATA inserted        time    *
        # post-range(posnum + 1)     *     *   value-added

        newdata <- data1[posnum,] # upper divided range created as new data
        data1[1:posnum, varname] <- default # lower range included
        if(nrow(data1) > posnum) # when later interval exists
          data1[(posnum+1):nrow(data1), varname] <- value
        newdata[, varname] <- value
        
        newdata[, startpos] <- time 
        data1[posnum, stoppos] <- time
        
        data1 <- rbind(data1, newdata)
        
      } else { # pre-value exist
        
        #                         tstart tstop cumcov
        # pre-range(posnum - 1)      *     *     *
        # NEW DATA inserted          *    time   *
        # target-range(posnum)     time    *   value-added
        # post-range(posnum + 1)     *     *   value-added
        
        newdata <- data1[posnum,] # lower divided range created as new data
        data1[posnum:nrow(data1), varname] <- data1[posnum:nrow(data1), varname] + value

        newdata[, stoppos] <- time
        data1[posnum, startpos] <- time
        
        data1 <- rbind(data1, newdata)
      }
    }
    
  } else if(time < tstart[1]){ # before the first interval

    if(length(varpos) == 0)
      data1[, varname] <- value
    else
      data1[, varname] <- data1[, varname] + value
    
  } else if(tstop[length(pos)] == time){ # covariate data added at the end of last interval

      #                         tstart    tstop   cumcov
      # pre-range(posnum - 1)      *       *        *
      # target-range(last)         *    time-eps value-added
      # NEW DATA(inserted)     time-eps   time   value-added
      
    if(length(varpos) == 0){ # no pre-data
      
      newdata <- data1[length(pos),]
      data1[, varname] <- default

      while(epstime > data1[length(pos), stoppos] - data1[length(pos), startpos]){
        epstime <- epstime/2
      }
      
      data1[length(pos), stoppos] <- data1[length(pos), stoppos] - epstime
      newdata[, varname] <- value
      newdata[, startpos] <- data1[length(pos), stoppos]
      if(data1[length(pos), startpos] >= data1[length(pos), stoppos])
        stop("epstime is too large (this error should not occure, maybe due to bug)")
      
      data1 <- rbind(data1, newdata)
      
    } else { # pre-data exist
    
      newdata <- data1[length(pos),]
      newdata[, varname] <- data1[length(pos), varname] + value

      while(epstime > data1[length(pos), stoppos] - data1[length(pos), startpos]){
        epstime <- epstime/2
      }
      
      newdata[, startpos] <-  data1[length(pos), stoppos] - epstime
      data1[length(pos), stoppos] <- newdata[, startpos]
      
      if(data1[length(pos), startpos] >= data1[length(pos), stoppos])
        stop("epstime is too large (this error should not occure, maybe due to bug)")

      data1 <- rbind(data1, newdata)
    }
  } else { # upper outside
    
    if(length(varpos) == 0) # case of no pre-value
      data1[, varname] <- default
  }
  
  return(data1)
}

event <- function(time, value, varname, data1, delay, na.rm, default, epstime, options){

  if(is.na(time) || is.na(value))
    return(data1)
  
  varlist <- names(data1)
  startpos <- grep("tstart", varlist)
  stoppos <- grep("tstop", varlist)
  time <- time + delay
  
  # interval exists
  if(length(startpos) > 0){
    
    tstart <- data1[,startpos]
    tstop <- data1[,stoppos]
    pos <- (tstart < time) & (time <= tstop)
    posnum <- sum((1:nrow(data1))*as.numeric(pos))
    
    if(sum(pos) > 1)
      stop("Dupicate range for one id")

    # inside of interval
    if(sum(pos) == 1){
    
      if(data1[posnum, stoppos] == time){ # event at the end of interval

        data1[      , varname] <- default
        data1[posnum, varname] <- value
        
      } else { # event at the middle of interval
        
        newdata <- data1[posnum,]
        data1[posnum, stoppos] <- time
        data1[posnum, varname] <- value
        newdata[, startpos] <- time
        newdata[, varname] <- default
        data1 <- rbind(data1, newdata)
        
      }
      
    } else { # interval exists, but time is outside of them
      
      data1 <- data1[order(data1[, startpos]),]
      
      if(data1[1, startpos] == time){ # event at the first time point

        newdata <- data1[1,]
        newdata[, startpos] <- data1[1, startpos] - epstime
        newdata[, stoppos] <- data1[1, startpos]
        newdata[, varname] <- value
        data1[, varname] <- default
        data1 <- rbind(newdata, data1)
        
      }
    }
    # do nothing in case of outside
    
  } else { # no interval exists
    data1[, "tstart"] <- 0
    data1[, "tstop"] <- time
    data1[, varname] <- value
  }
  
  return(data1)
}

tdispatch <- function(id, id1, id2, data1, data2, eqns, varnames, delay, na.rm, default, epstime, options){
  data1 <- data1[id1==id, , drop=FALSE]
  data2 <- data2[id2==id, , drop=FALSE]
  if(nrow(data2) > 1)
    stop("more elements supplied than there are to replace")

  for(i in 1:length(eqns)){
    if(length(eqns[[1]]) < 3){
      if(length(eqns[[1]]) == 2 && eqns[[1]][[1]] == "event")
        eqns[[i]][[length(eqns[[i]]) + 1]] <- 1
      else
        stop("incorrect number of argument")
    }
    eqns[[i]][[length(eqns[[i]]) + 1]] <- varnames[i]
    eqns[[i]][[length(eqns[[i]]) + 1]] <- as.data.frame(data1)
    eqns[[i]][[length(eqns[[i]]) + 1]] <- delay
    eqns[[i]][[length(eqns[[i]]) + 1]] <- na.rm
    eqns[[i]][[length(eqns[[i]]) + 1]] <- default
    eqns[[i]][[length(eqns[[i]]) + 1]] <- epstime
    eqns[[i]][[length(eqns[[i]]) + 1]] <- options
    ret <- eval(eqns[[i]], data2)
  }
  return(ret)
}

tdmerge <- function(data1, data2, id, ..., tstart, tstop, options){

  Call <- match.call()

  if(missing(data1) || missing(data2) || missing(id))
    stop("data1, data2, and id arguments are required")
  if(!is.data.frame(data1) || !is.data.frame(data2))
    stop("data1 and data2 must be a data.frame")
    
  id2 <- eval(Call[["id"]], data2, enclos = emptyenv())
  id <- unique(id2[order(id2)])
  id1 <- eval(Call[["id"]], data1, enclos = emptyenv())

  if(!missing(tstart)){
    
    if(missing(tstop))
      stop("either tstart or tstop cannot be specified")

    data1[, "tstart"] <-  eval(Call[["tstart"]], data2)

    data1[, "tstop"] <- eval(Call[["tstop"]], data2)

  } else {
    if(!missing(tstop))
      stop("either tstart or tstop cannot be specified")

    tstart <- NULL
    tstop <- NULL
  }

  eqns <- substitute(as.list(...))[-1]
  varnames <- names(eqns)
  na.rm <- TRUE
  delay <- 0
  default <- 0
  epstime <- 0.1
  
  if(!missing(options)){
    
    options <- as.list(options)
    if(!is.null(options[["cmd"]]) || is.character(options[["cmd"]])){
      eqns <- list(eqns, parse(text=options[["cmd"]]))
    }
    if(!is.null(options[["na.rm"]]))
      na.rm <- options[["na.rm"]]
    if(!is.null(options[["delay"]]))
      delay <- options[["delay"]]
    if(!is.null(options[["default"]]))
      default <- options[["default"]]
    if(!is.null(options[["epstime"]]))
      epstime <- options[["epstime"]]
    
  } else {
    options <- list()
  }
  
  res <- lapply(id, tdispatch, id1=id1, id2=id2, data1=data1, data2=data2,
                eqns=eqns, varnames=varnames, delay=delay, na.rm=na.rm,
                default=default, epstime=epstime, options=options)

  rowlengths <- sapply(res, nrow)
  totlength <- sum(rowlengths)
  listlength <- length(res)
  reslist <- names(res[[1]])
  
  ret <- list()
  for(var in reslist){
    ret[[var]] <- res[[1]][,var]
    if(listlength > 1){
      length(ret[[var]]) <- totlength
      curlength <- rowlengths[1]
      for(i in 2:listlength){
        ret[[var]][(curlength+1):(curlength + rowlengths[i])] <- res[[i]][,var]
        curlength <- curlength + rowlengths[i]
      }
    }
  }
  return(as.data.frame(ret))
}
