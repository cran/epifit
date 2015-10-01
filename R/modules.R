## Likelihood function
LogLikelihood <- function(init, num_nmodel, envs, weight, itereq,
                          lst_psdresponse, lst_psdrespvar, lst_psdparam, lst_option, # model
                          vec_ranvar, vec_rparameter, vec_randens, lst_psdrparam, # random
                          vec_parameter, lst_parameter # variable list
                          ){
  ret <- 0
  
  for(i in 1:num_nmodel){

    if(lst_psdresponse[[i]][[1]] == "cox"){
      
      if(is.null(lst_option[[i]]))
        ties <- "efron"
      else
        ties <- lst_option[[i]]$"ties"

      ret <- ret + LogCoxLikelihood(init[GetParamPosition(lst_parameter[[i]], vec_parameter)],
                                    lst_parameter[[i]], lst_psdrespvar[[i]], itereq[i],
                                    envs[[i]], "time1_inner_", "time2_inner_", "status_inner_",
                                    ties=ties, "weight_inner_")

    } else {

      if(is.environment(envs[[i]][[1]])){ # no random effect

        ret <- ret +
          InnerLogLikelihood(init[GetParamPosition(lst_parameter[[i]],vec_parameter)],
                          lst_parameter[[i]], lst_psdresponse[[i]][[1]], lst_psdparam[[i]],
                          lst_psdrespvar[[i]], envs[[i]][[1]])
        
      } else { # with random effect
        stop("reached unreachable code area in LogLikelihood function")
      }
    }
  }
  
  return(ret)
}

InnerLogLikelihood <- function(init, parameters, distname, lst_psdparam, psdrespvar, env){
  for(i in 1:length(parameters))
    assign(parameters[i], init[i], envir=env)

  if(distname=="pois"){
    if(exists("weight_inner_", envir=env)){
      return(-sum(x=dpois(eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env), log=TRUE)*
                  get("weight_inner_", envir=env)))
    } else {
      return(-sum(x=dpois(eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env), log=TRUE)))
    }
  } else if(distname=="norm"){

    if(exists("weight_inner_", envir=env)){
      return(-sum(dnorm(x=eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env),
                         eval(lst_psdparam[[2]], envir=env), log=TRUE)*get("weight_inner_", envir=env)))
    } else {
      return(-sum(dnorm(x=eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env),
                         eval(lst_psdparam[[2]], envir=env), log=TRUE)))
    }
    
  } else if(distname=="binom"){
    
    if(exists("weight_inner_", envir=env)){
      return(-sum(dbinom(x=eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env),
                         eval(lst_psdparam[[2]], envir=env), log=TRUE)*get("weight_inner_", envir=env)))
    } else {
      return(-sum(dbinom(x=eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env),
                         eval(lst_psdparam[[2]], envir=env), log=TRUE)))
    }
    
  } else if(distname=="nbinom"){
    if(exists("weight_inner_", envir=env)){
      return(-sum(dnbinom(x=eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env),
                         eval(lst_psdparam[[2]], envir=env), log=TRUE)*get("weight_inner_", envir=env)))
    } else {
      return(-sum(dnbinom(x=eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env),
                         eval(lst_psdparam[[2]], envir=env), log=TRUE)))
    }
    
  } else if(distname=="gamma"){
    
    if(exists("weight_inner_", envir=env)){
      return(-sum(dgamma(x=eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env),
                 eval(lst_psdparam[[2]], envir=env), log=TRUE)*get("weight_inner_", envir=env)))
    } else {
      return(-sum(dgamma(x=eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env),
                         eval(lst_psdparam[[2]], envir=env), log=TRUE)))
    }
    
  } else if(distname=="weibull"){

    if(exists("weight_inner_", envir=env)){
      return(-sum(dweibull(x=eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env),
                 eval(lst_psdparam[[2]], envir=env), log=TRUE)*get("weight_inner_", envir=env)))
    } else {
      return(-sum(dweibull(x=eval(psdrespvar, envir=env), eval(lst_psdparam[[1]], envir=env),
                         eval(lst_psdparam[[2]], envir=env), log=TRUE)))
    }
    
  } else if(distname=="general"){
    if(exists("weight_inner_", envir=env)){
      return(-sum(eval(psdrespvar, envir=env)*get("weight_inner_", envir=env)))
    } else {
      return(-sum(eval(psdrespvar, envir=env)))
    }
  }
}

## obtain the first environment from "envs" list
GetFirstEnvironment <- function(envs){
  while(!is.environment(envs))
    envs <- "[["(envs, 1)
  return(envs)
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
          unsolved <- RemoveVariable(unsolved, lst_eqndepend[[i]])
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
RemoveVariable <- function(varlist, remove){
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
      tmp <- strsplit(options[i], "\\(")[[1]]
      ret[[i]] <- strsplit(tmp[2], "\\)")[[1]][1]
      name <- c(name, tmp[1])
    }
    names(ret) <- name
    return(ret)
  }
}

## Inner function for "modeleq" argument
GetResponse <- function(modelstr){
  tmp <- strsplit(modelstr, "~")[[1]][1]
  return(strsplit(tmp, "/")[[1]][1])
}

## Inner function for "modeleq" argument
GetRespvar <- function(modelstr){
  return(strsplit(modelstr, "~")[[1]][2])
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
  })
}

## Make epifit result object from optim function
MakeResultFromOptim <- function(result, ans, nulllik){
  result$coefficients <- ans$par
  result$loglik <- c(-nulllik, -ans$value)
  result$var <- ginv(ans$hessian)
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
MakeResultFromNlm <- function(result, ans, nulllik){
  result$coefficients <- ans$estimate
  result$loglik <- c(-nulllik, -ans$minimum)
  result$var <- ginv(ans$hessian)
  result$iter <- ans$iterations
  result$convergence <- ans$code
  result$wald.test <- t(ans$estimate)%*%(ans$hessian)%*%(ans$estimate)
  return(result)
}

## Checking supported operator is included
InnerListVariable <- function(tree){
  supported <- c("-","+","*","/","^","<",">","==",">=","<=", "&","|","(",")",
                 "abs","acos","acosh", "as.integer","as.numeric","asin","asinh","atan",
                 "atanh","cos","cosh", "digamma","exp","expm1","factorial","floor",
                 "gamma","ifelse","lgamma", "lfactorial","log","log10","log1p","log2",
                 "logb","pmax","pmax.int", "pmin","pmin.int","print",
                 "sin","sinh","tan","tanh","trigamma")  
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
  } else if(!op %in% supported){
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

## obtain variable list summarized over sentences ({assigned}, {arguments})
ListVariable <- function(expression){
  expression <- DivideExpression(expression)
  res <- list(character(0), character(0))
  for(i in 1:length(expression)){
    tryCatch({tree <- ParseLine(expression[i])},
             warning=function(e){print(e)},
             error=function(e){stop("Invalid expression: check equation")},
             finally={}
             )
    tmp <- InnerListVariable(tree)
    res[[1]] <- c(res[[1]], tmp[[1]])
    res[[2]] <- c(res[[2]], tmp[[2]])
  }
  return(res)
}

## Parse one line text
ParseLine <- function(expression){
  return(parse(text=expression)[[1]])
}

## Obtain result per sentence
ListVariablePerSentence <- function(expression){
  expression <- DivideExpression(expression)
  res <- list(NULL)
  for(i in 1:length(expression)){
    res[[i]] <- InnerListVariable(expression[i])
  }
  return(res)
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

## Partial likelihood function for Cox regressions
LogCoxLikelihood <- function(init, parameters, equations, itereq, envs, time1name, time2name, statusname, ties=c("efron","breslow","average","discrete"), weightname){

  ret <- 0
  ties <- match.arg(ties)

  for(strata in 1:length(envs)){

    ## Assign parameter value
    for(i in 1:length(parameters))
      assign(parameters[i], init[i], envir=envs[[strata]])

    time2 <- get(time2name, envir=envs[[strata]])
    status <- get(statusname, envir=envs[[strata]])    
    if(!exists(time1name, envir=envs[[strata]])){
      time1 <- rep(0, length(time2))
    } else {
      time1 <- get(time1name, envir=envs[[strata]])
    }
    
    if(nchar(itereq) > 0){
      assign("T", rep(time2[1], length(time2)), envir=envs[[strata]])    
      eval(itereq, envir=envs[[strata]])
    }
    
    hazard <- eval(equations, envir=envs[[strata]])
    whazard <- NULL
    phazard <- NULL

    if(!exists(weightname, mode="numeric", envir=envs[[strata]])){
      weight <- rep(1, length(hazard))
      whazard <- hazard
      phazard <- hazard
    } else {
      if(ties=="discrete" || ties=="average")
        stop("weight is not supported for ties=\"", ties,"\" specification")
      weight <- get(weight, envir=envs[[strata]])
      whazard <- hazard*weight
      phazard <- hazard^weight
    }
  
    nsubject <- length(hazard)
    result <- 0
    riskset <- numeric(nsubject)
    pretime <- time2[1]
    duringtie <- FALSE
    tiebegin <- 0
    
    riskset[1] <- sum(as.numeric((time1 < time2[1]) & (time2[1] <= time2))*whazard)
    
    ## Main calculation loop begins here
    for(i in 2:nsubject){
      
      if(nchar(itereq) > 0){
        assign("T", rep(time2[i], nsubject), envir=envs[[strata]])    
        eval(itereq, envir=envs[[strata]])
        hazard <- eval(equations, envir=envs[[strata]])
        phazard <- hazard^weight
        whazard <- hazard*weight
      }
      
      riskset[i] <- sum(as.numeric((time1 < time2[i]) & (time2[i] <= time2))*whazard)
      
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
              phazard[tiebegin] <- prod(phazard[tiebegin:(tiebegin+tieevent-1)])/.Call("Rf_select", tieevent, nsubject-tiebegin+1, phazard[tiebegin:nsubject])
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

## obsolete
## Check operator in equations (currently not yet used)
## Included in InnerListVariable
CheckOperator <- function(expression){
  supported <- c("-","+","*","/","^","<",">","==",">=","<=","<-", "=", "&","|",
                 "abs","acos","acosh", "as.integer","as.numeric","asin","asinh","atan",
                 "atanh","cos","cosh", "digamma","exp","expm1","factorial","floor",
                 "gamma","ifelse","lgamma", "lfactorial","log","log10","log1p","log2",
                 "logb","pmax","pmax.int", "pmin","pmin.int","print",
                 "sin","sinh","tan","tanh","trigamma","(",")")
  tryCatch({tree <- ParseLine(expression)},
           warning=function(e){print(e)},
           error=function(e){stop("Invalid expression: check equation")},
           finally={}
           )
  if(length(tree) == 1){
    return()
  }
  if(!(as.character(tree[[1]]) %in% supported)){
    stop(paste("unsupported operator is found:", as.character(tree[[1]])), sep=" ")
  }
  for(i in 2:length(tree)){
    CheckOperator(deparse(tree[[i]]))
  }
}

## return assigned and varlist
InnerListVariableOld <- function(expression){
  tryCatch({tree <- parse(text=expression)[[1]]},
           warning=function(e){print(e)},
           error=function(e){stop("Invalid expression: check equation")},
           finally={}
           )
  res <- list(character(0), character(0))
  if(length(tree) == 1){
    res[[2]] <- c(res[[2]], as.character(tree))
    return(res)
  }
  if(as.character(tree[[1]]) == "<-" || as.character(tree[[1]]) == "="){
    res[[1]] <- c(res[[1]], as.character(tree[[2]]))
    tmp <- InnerListVariableOld(deparse(tree[[3]]))
    res[[1]] <- c(res[[1]], tmp[[1]])
    res[[2]] <- c(res[[2]], tmp[[2]])
  } else if(as.character(tree[[1]]) == "->"){
    res[[1]] <- c(res[[1]], as.character(tree[[3]]))
    tmp <- InnerListVariableOld(deparse(tree[[2]]))
    res[[1]] <- c(res[[1]], tmp[[1]])
    res[[2]] <- c(res[[2]], tmp[[2]])
  } else {
    for(j in 2:length(tree)){
      tmp <- InnerListVariableOld(deparse(tree[[j]]))
      res[[1]] <- c(res[[1]], tmp[[1]])
      res[[2]] <- c(res[[2]], tmp[[2]])
    }
  }
  return(res)
}

## obtain variable list summarized over sentences ({assigned}, {arguments})
ListVariableOld <- function(expression){
  expression <- DivideExpression(expression)
  res <- list(character(0), character(0))
  for(i in 1:length(expression)){
    tmp <- InnerListVariable(expression[i])
    res[[1]] <- c(res[[1]], tmp[[1]])
    res[[2]] <- c(res[[2]], tmp[[2]])
  }
  return(res)
}
