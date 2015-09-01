#' This function maximizes an arbitrary likelihood including generalized linear models and Cox partial likelihood.
#'
#' This function provides flexible model fitting. The main model specification is written in \code{modeleq}. \code{modeleq} consisits of two parts separated by \sQuote{~}. The distribution is specified in the first part, and the second part includes variable name which follows the specified distribution in the first part. Available distributional specifications are \sQuote{cox}(Cox partial likelihood), \sQuote{pois}(Poisson distribution), \sQuote{norm}(normal distribution), \sQuote{binom}(binomial distribution), \sQuote{nbinom}(negative binomial distribution), \sQuote{gamma}(gamma distribution) and \sQuote{weibull}(Weibull distribution). Options can be specified for some distribution specification after \sQuote{/} in the distribution specification part. Multiple options separated by \sQuote{,} can also be specified.
#'
#' For Cox regressions, time and status variable must be specified in parentheses like \code{cox(time, status)}. Some options are also available for Cox regressions, and \sQuote{efron}, \sQuote{breslow}, \sQuote{discrete} and \sQuote{average} is available for tie handling method. \code{ties(discrete)} specification corresponds to \sQuote{exact} ties specification in \code{coxph} function, and \code{ties(average)} corresonds to \sQuote{exact} specification in SAS PHREG procecure. See references for further details. Strata option which specifies a variable indicating strata is also available in Cox regressions. Subset option which has same functinality as subset argument below is also available for Cox regressions and other distribution specifications. For other distribution specifications, parameters must be specified in parentheses. For poisson distribution, mean parameter must be specified as \code{pois(mean)}. Note that each parameter specificaiton can be a variable or R equation. For other distributions, \code{norm(mean, variance)}, \code{binom(size, probability)}, \code{nbinom(size, probability)}, \code{gamma(shape, scale)}, \code{weibull(shape, scale)}.
#'
#' When distributions are specified (not Cox regressions), additional R expressions can be specified in \code{equations} argument. R expressions are parsed to make variable list. Variables which exist in data.frame or the global environment must be vector, and the rest of variables are regarded as parameters. If you define variable \sQuote{mu} in \code{equations}, you can use \sQuote{mu} in \code{modeleq} argument. Refer Poisson regression examples below.
#' 
#' @title Model fitting function for epifit package
#' @param modeleq a character string specifying the model. See \sQuote{Details}.
#' @param equations a character string specifying R expressions executed before \sQuote{modeleq}. Multiple expressions separated by \sQuote{;} is allowed. See \sQuote{Details}.
#' @param initeq a character string specifying R expressions executed once at the first time. Typical use is specifying initial values for parameters using R expressions. Multiple expressions separated by \sQuote{;} is allowed.
#' @param itereq a character string specifying R expressions executed during Cox regression at each time of event occur. Typical use is incorporating time-dependent variables. Not yet implemented completely.
#' @param endeq a character string specifying R expressions executed at the end of calculations.
#' @param data a data.frame in which to interpret the variables named in the formula, or in the subset and the weights argument.
#' @param subset expression indicating which subset of the rows of data should be used in the fit. All observations are included by default.
#' @param weight vector of case weights. If weights is a vector of integers, then the estimated coefficients are equivalent to estimating the model from data with the individual cases replicated as many times as indicated by weights. Only supported \sQuote{breslow} and \sQuote{efron} ties specification in Cox regression models and other likelihood specifications.
# @param random a character string specifying random effects. See \sQuote{Details}.
#' @param na.action a missing-data filter function. This is applied when data.frame is supplied as \sQuote{data} parameter. Default is \code{options()$na.action}.
#' @param opt a character string specifying the method for optimization. When sQuote{newrap} is specified, \code{nlm} function that uses Newton-type algorithm used for obtaining maximum likelihood estimate. For the rest of specifications (\sQuote{BFGS} for quasi-Newton method, \sQuote{CG} for conjugate gradients method, \sQuote{SANN} for simulated annealing, \sQuote{Nelder-Mead} for Nelder and Mead simplex method), \code{optim} is used. The default is \sQuote{newrap}.
#' @param tol1 a numeric value specifying \code{gtol} in nlm, \code{abstol} in optim.
#' @param tol2 a numeric value specifying \code{stol} in nlm, \code{reltol} in optim.
#' @param maxiter a integer value specifying the maximum number of iterations. Defaults is 200.
#' @param init Initial values for the parameters specified in the form of vector.
#' @param verbatim a integer value from 0 (minimum) to 2 (maximum) controlling the amount of information printed during calculation.
#' @param ... For the arguments used in the inner functions (currently not used).
#' @return A list containing the result of model fitting including parameter estimates, variance of parameter estimates, log likelihood and so on.
#' @references DeLong, D. M., Guirguis, G.H., and So, Y.C. (1994). Efficient computation of subset selection probabilities with application to Cox regression. \emph{Biometrika} \strong{81}, 607-611.
#' @references Gail, M. H., Lubin, J. H., and Rubinstein, L. V. (1981). Likelihood calculations for matched case-control studies and survival studies with tied death times. \emph{Biometrika} \strong{68}, 703-707.
#' @examples
#' library(survival)
#' 
#' # Make sample data
#' set.seed(123)
#' nsub <- 20
#' follow <- 5
#' x <- rnorm(nsub)
#' rate <- exp(-0.5 + x)
#' etime <- rweibull(nsub, 1, 1/rate)
#' status <- as.integer(etime < follow)
#' time <- pmin(follow, etime)
#' dat <- data.frame(status, x, time, ltime=log(time))
#'
#' coxph(Surv(time, status)~x, data=dat)
#' modeleq <- "cox(time,status)~exp(beta*x)"
#' epifit(modeleq=modeleq, data=dat)
#'
#' glm(status ~ x + offset(ltime), family=poisson(), data=dat)
#' equations <- "mu <- time*exp(beta0 + beta1*x)"
#' modeleq <- "pois(mu) ~ status"
#' epifit(modeleq=modeleq, equations=equations, data=dat)
#'
#' # The simplest test data set from coxph function
#' test1 <- list(time=c(4,3,1,1,2,2,3),
#'               status=c(1,1,1,0,1,1,0),
#'               x=c(0,2,1,1,1,0,0),
#'               sex=c(0,0,0,0,1,1,1))
#'
#' # Cox regressions with strata
#' coxph(Surv(time, status) ~ x + strata(sex), data=test1)
#' modeleq <- "cox(time,status)/strata(sex)~exp(beta*x)"
#' epifit(modeleq=modeleq, data=test1)
#'
#' # Tie specification example in Cox regressions
#' coxph(Surv(time, status) ~ x + strata(sex), data=test1, ties="breslow")
#' modeleq <- "cox(time,status)/strata(sex),ties(breslow)~exp(beta*x)"
#' epifit(modeleq=modeleq, data=test1)
#' 
#' # Average partial likelihood
#' modeleq <- "cox(time,status)/strata(sex),ties(average)~exp(beta*x)"
#' epifit(modeleq=modeleq, data=test1)
#'
#' # Conditional logistic regression for matched case-control studies
#' # hypothetical data
#' conlog <- data.frame(strata=c(1,1,2,2,3,3,4,4,5,5), outcome=c(1,0,1,0,1,0,1,0,1,0),
#'                      cov=c(1,3,2,1,5,2,4,2,2,2))
#' # Make dummy survival time so that all the cases in a matched set have the same survival
#' # time value, and the corresponding controls are censored at later times
#' conlog <- cbind(conlog, dummy=(2 - conlog$outcome))
#' coxph(Surv(dummy, outcome)~cov + strata(strata), ties="exact", data=conlog)
#' modeleq <- "cox(dummy,outcome)/ties(discrete),strata(strata)~exp(beta*cov)"
#' epifit(modeleq=modeleq, data=conlog)
#' @export
epifit <- function(modeleq, equations="", initeq="", itereq="", endeq="",
                   data, subset, weight, na.action,
                   opt=c("newrap", "BFGS", "CG", "Nelder-Mead", "SANN"),
                   tol1=1e-8, tol2=1e-8, maxiter=200, init, verbatim=0, ...){

  ## argument processing
  args <- list(...)
  Call <- match.call(expand.dots=FALSE)
  opt <- match.arg(opt)
  
  if(missing(modeleq))
    stop("modeleq cannot be omitted")

  if(missing(na.action)){
    na.action <- options()$na.action
  }

  if(na.action=="na.fail"){
    
    if(missing(data)){
      data <- .GlobalEnv
    } else {
      data <- list2env(na.fail(data), parent=.GlobalEnv)
    }

  } else if(na.action=="na.omit"){

    if(missing(data)){
      data <- .GlobalEnv
    } else {
      data <- list2env(na.omit(data), parent=.GlobalEnv)
    }
    
  } else if(na.action=="na.exclude"){

    if(missing(data)){
      data <- .GlobalEnv
    } else {
      data <- list2env(na.exclude(data), parent=.GlobalEnv)
    }
    
  } else if(na.action=="na.pass"){

    if(missing(data)){
      data <- .GlobalEnv
    } else {
      data <- list2env(na.pass(data), parent=.GlobalEnv)
    }

  } else {
      stop("invalid na.action specification")
  }

  ## number of models
  num_nmodel <- length(modeleq)
  
  ## response variable
  lst_option <- list(NULL) # option list
  lst_psdparam <- list(NULL) # parsed 1st, 2nd 3rd element of modeleq *******
  lst_psdresponse <- list(NULL) # parsed response including model specification
  lst_psdrespvar <- list(NULL) # parsed respvar
  
  ## temporary variables
  vec_response <- character(num_nmodel) # left part of modeleq per model (temporary)
  vec_respvar <- character(num_nmodel) # right part of modeleq per model (temporary)
  vec_fmlparam <- character(0) # parameter formula (temporary)
  vec_respdepend <- character(0) # respvar depend variable list (temporary)
  eqn_equation <- character(0) # divided equations (temporary)
  eqn_initeq <- character(0) # divided initieq (temporary)
  eqn_itereq <- character(0) # divided itereq (temporary)
  eqn_endeq <- character(0) # divided endeq (temporary)

  ## random effects related variables
  psd_random <- NULL # parsed random effect (temporary)
  lst_psdrparam <- list(NULL) # random effect param formula (1st..mu, 2nd..sigma for normal)
  lst_roption <- list(NULL) # option list for random effect
  vec_subject <- numeric(0) # variable specifying subject in random effect
  vec_ranvar <- character(0) # random effect variable
  vec_randens <- character(0) # random effect density function

  ## initialize variable lists (varlist=(variables+parameters), assigned)
  ## for current model
  vec_varlist <- character(0) # all varlist including parameter (temporary working)
  vec_depend <- character(0) # dependent variable per parameter (temporary working)
  
  vec_assigned <- character(0) # all assigned variable
  vec_variable <- character(0) # all variables
  vec_parameter <- character(0) # all parameters
  ## per sentence
  lst_varlist <- list(NULL) # reserve variable list per sentence
  ## per model
  lst_parameter <- list(NULL) # parameter list per model (do not exist in the data)
  lst_variable <- list(NULL) # variable list per model (exist in the data)
  ## random effect
  vec_rparameter <- character(0)
  vec_rvariable <- character(0)

  ## dependent variables in "equations"
  lst_eqnassigned <- list(NULL)
  lst_eqndepend <- list(NULL)

  if(num_nmodel > 1){
    if(length(equations) != num_nmodel){
      if(length(equations) != 1){
        stop("number of equations is not equal to that of models")
      }
      equations <- rep(equations, num_nmodel)
    }
    if(length(initeq) != num_nmodel){
      if(length(initeq) != 1){
        stop("number of initeq is not equal to that of models")
      }
      initeq <- rep(initeq, num_nmodel)      
    }
    if(length(itereq) != num_nmodel){
      if(length(itereq) != 1){
        stop("number of itereq is not equal to that of models")
      }
      itereq <- rep(itereq, num_nmodel)      
    }
    if(length(endeq) != num_nmodel){
      if(length(endeq) != 1){
        stop("number of endeq is not equal to that of models")
      }
      endeq <- rep(endeq, num_nmodel)      
    }
  }

  ## random effect (enabled when multivariate improper integral is implemented)
  ## if(nchar(random) > 0){
  ##   lst_roption <- GetOptions(random)
  ##   vec_ranvar <- trim(GetRespvar(random))
  ##   psd_random <- parse(text=GetResponse(random))
    
  ##   if(psd_random[[1]][[1]] == "normal"){ # random effect distribution
  ##     vec_randens <- "dnorm"
  ##   } else {
  ##     stop("not supported distribution in random effect specification")
  ##   }

  ##   if(is.null(lst_roption$"subject"))
  ##     stop("subject option must be specified in random effect")
  ##   vec_subject <- eval(parse(text=lst_roption$"subject"), envir=data)

  ##   lst_psdrparam[[1]] <- list(NULL)
  ##   for(i in 2:length(psd_random[[1]])){
  ##     lst_psdrparam[[1]][[i-1]] <- as.character(psd_random[[1]][[i]])
  ##     vec_varlist <- c(vec_varlist, InnerListVariable(psd_random[[1]][[i]])[[2]])
  ##   }

  ##   vec_varlist <- unique(vec_varlist)

  ##   tmp <- ClassifyParameter(vec_varlist, data)
  ##   vec_rparameter <- c(vec_rparameter, as.character(tmp[[1]]))
  ##   vec_rvariable <- c(vec_rvariable, as.character(tmp[[2]]))

  ##   vec_varlist <- character(0) # initialize for later use
  ## } # random effect end

  for(i in 1:num_nmodel){

    ## left term is response variable specification
    tmpstr <- strsplit(modeleq[i], "~")[[1]]
    
    if(length(tmpstr) != 2)
      stop("\"modeleq\" must include response variable as left term separated by \"~\"")

    ## left side of "modeleq" with option specifications
    vec_response[i] <- trim(tmpstr[1])
    
    ## right term of "modeleq" should be evaluated
    vec_respvar[i] <- trim(tmpstr[2])
    
    lst_option[[i]] <- GetOptions(vec_response[i])
    vec_response[i] <- strsplit(vec_response[i], "/")[[1]][1] # remove option specification
    
    ## Divide argument "equations"
    if(nchar(equations[i]) > 0){
      eqn_equation <- DivideExpression(equations[i])
    } else {
      eqn_equation <- character(0)
    }
    
    ## Divide argument "itereq"
    if(nchar(itereq[i]) > 0){
      eqn_itereq <- DivideExpression(itereq)
    } else {
      eqn_itereq <- character(0)
    }
    
    ## Divide respvar
    if(length(grep(";", vec_respvar[i])) > 0)
      stop("right side of modeleq must be one sentence")
    eqn_respvar <- vec_respvar[i]
    
    ## eval vec_response
    lst_psdresponse[[i]] <- parse(text=vec_response[i])[[1]]
    
    ## get parameter formula
    vec_fmlparam <- character(length(lst_psdresponse[[i]])-1)
    
    ## check number of parameters
    if(length(vec_fmlparam) == 0){
      stop("parameter should be specified in modeleq")
    } else if(length(vec_fmlparam) == 1){
      if(!as.character(lst_psdresponse[[i]][[1]]) %in% c("pois")){
        stop("number of parameters specified in \"modeleq\" is incorrect")
      }
    } else if(length(vec_fmlparam) == 2){
      if(!as.character(lst_psdresponse[[i]][[1]]) %in% c("cox", "norm", "binom",
                                                      "nbinom", "gamma", "weibull")){
        stop("number of parameters specified in \"modeleq\" is incorrect")
      }
    } else if(length(vec_fmlparam) == 3){
      if(!as.character(lst_psdresponse[[i]][[1]]) %in% c("cox")){
        stop("number of parameters specified in \"modeleq\" is incorrect")
      }
    } else {
        stop("number of parameters specified in \"modeleq\" is incorrect")
    }
    
    for(j in 2:length(lst_psdresponse[[i]])){
      tmp <- deparse(lst_psdresponse[[i]][[j]])
      if(length(grep(";", tmp) > 0))
        stop("each parameter must be conposed of one sentence")
      vec_fmlparam[j-1] <- tmp
    }

    ## conbine all equations
    eqns <- c(eqn_equation, eqn_itereq, eqn_respvar, vec_fmlparam)
    tmp <- ListVariable(eqns)
    vec_varlist <- tmp[[2]]

    tmp <- ClassifyParameter(vec_varlist, data, tmp[[1]]) # remove assigned variable(tmp[[1]])
    lst_parameter[[i]] <- tmp[[1]]
    lst_variable[[i]] <- tmp[[2]]
    
    ## make dependent variable list for equations
    ## (lst_eqndepend, lst_eqnassigned)
    if(nchar(equations[i]) > 0){
      for(j in 1:length(eqn_equation)){
        tmp <- InnerListVariable(eqn_equation[j])
        if(length(tmp[[1]]) > 0)
          lst_eqnassigned[[j]] <- tmp[[1]]
        if(length(tmp[[2]]) > 0)
          lst_varlist[[j]] <- tmp[[2]]
      }
      for(j in 1:length(eqn_equation)){
        if(length(LimitVarlist(unlist(lst_varlist[[j]]), unlist(lst_eqnassigned[[j]]))) > 0)
          lst_eqndepend[[j]] <- LimitVarlist(unlist(lst_varlist[[j]]), unlist(lst_eqnassigned[[j]]))
        else
          lst_eqndepend[[j]] <- c("")
      }

      ## warn when the data is changed
      vec_assigned <- unlist(lst_eqnassigned)
      if(length(vec_assigned) != 0 && sum(vec_assigned %in% lst_variable[[i]]))
        warning("Values of some objects will be changed during calculations")
      
      ## Make dependence list
      vec_varlist <- InnerListVariable(eqn_respvar)[[2]]
      vec_depend <- LimitVarlist(vec_varlist, vec_assigned)
    } else {
      vec_depend <- character(0)
    }

    if(length(vec_depend) > 0){
      vec_depfml <- SolveDependence(vec_depend, lst_eqnassigned, lst_eqndepend, eqn_equation)
      lst_psdrespvar[[i]] <- InsertFormula(parse(text=eqn_respvar), vec_depfml)
    } else {
      lst_psdrespvar[[i]] <- parse(text=eqn_respvar)
    }
    
    ## random effect in respvar
    ##if(length(vec_ranvar) > 0){
      ## not yet implemented
    ##}

    lst_psdparam[[i]] <- list(NULL)
    for(j in 1:length(vec_fmlparam)){
      vec_varlist <- InnerListVariable(vec_fmlparam[j])[[2]]
      vec_depend <- LimitVarlist(vec_varlist, vec_assigned)
      if(length(vec_depend) > 0){
        vec_depfml <- SolveDependence(vec_depend, lst_eqnassigned, lst_eqndepend, eqn_equation)
        lst_psdparam[[i]][[j]] <- InsertFormula(parse(text=vec_fmlparam[j]), vec_depfml)
      } else {
        lst_psdparam[[i]][[j]] <- parse(text=vec_fmlparam[j])
      }
      ## random effect in parameter
      ##if(length(vec_ranvar) > 0){
        ## not yet implemented
      ##}
    }
  }

  ## summarize parameters and variables 
  vec_parameter <- unique(unlist(lst_parameter))
  vec_variable <- unique(unlist(lst_variable), vec_rvariable)

  if(length(vec_parameter) < 1){
    stop("there are no parameters to estimate")
  }
  
  ## obtain number of subjects from the first variable
  nsubject <- length(get(vec_variable[1], envir=data))

  if(!missing(weight)){
    if(is.character(weight)){
      weight <- get(weight, envir=data)
    } else {
      weight <- eval(weight, envir=data)
    }
  } else {
    weight <- NULL
  }
  
  ## Get intial value from init argument
  if(!missing(init)){
    init <- rep(init, length(vec_parameter))
  } else {
    init <- rep(0, length(vec_parameter))
  }

  envs <- list(NULL)
  envstrata <- list(NULL)
  envrand <- list(NULL)
  
  ## make environments envs[[model, strata, random]]
  for(i in 1:num_nmodel){

    ## subset processing (commented out for debug)
    if(!missing(subset)){
      index <- rep(TRUE, nsubject) & as.logical(eval(parse(text=subset), envir=data))
    } else {
      index <- rep(TRUE, nsubject)
    }
    
    ## Cox regression
    if(as.character(lst_psdresponse[[i]][[1]])=="cox"){
      
      if(!is.null(lst_option[[i]]$"subset")){
        idx <- index & as.logical(eval(parse(text=lst_option[[i]]$"subset"), envir=data))
      } else {
        idx <- index
      }

      ## Counting process type input
      ## TODO: check dimension of variables
      if(length(lst_psdresponse[[i]])==4){
        time1 <- eval(parse(text=lst_psdresponse[[i]][[2]]), envir=data)
        time2 <- eval(parse(text=lst_psdresponse[[i]][[3]]), envir=data)
        status <- eval(parse(text=lst_psdresponse[[i]][[4]]), envir=data)
      } else if(length(lst_psdresponse[[i]])==3){
        time1 <- rep(0, nsubject)
        time2 <- eval(parse(text=lst_psdresponse[[i]][[2]]), envir=data)
        status <- eval(parse(text=lst_psdresponse[[i]][[3]]), envir=data)
      } else { # length(resp) is not between 3 and 4
        stop("unsupported model specification for Cox regression")
      }
      
      if(min(status) < 0 || max(status) > 1){
        stop("range of status variable is not between 0 and 1")
      }      

      ## Make environments (strata/random variable)
      if(!is.null(lst_option[[i]]$"strata")){ # with strata
        strata <- eval(parse(text=lst_option[[i]]$"strata"), envir=data)[idx]
        stratalist <- unique(strata[order(strata)])
        
        for(j in 1:length(stratalist)){
          
          idx2 <- idx & (strata==stratalist[j]) # strata restriction within model
          
          if(length(vec_ranvar)==0){ # no random effect
            
            envstrata[[j]] <- new.env()
            orderedtime <- order(time2[idx2])

            for(k in 1:length(vec_variable)){
              assign(vec_variable[k], get(vec_variable[k], envir=data)[idx2][orderedtime], envir=envstrata[[j]])
            }

            ## assign null model parameter
            ## null parameter value "nullpara" can be changed
            ## by assigning other values to nullpar in "initeq"
            assign("nullpara", rep(0, length(vec_parameter)), envir=envstrata[[j]])

            assign("time1_inner_", time1[idx2][orderedtime], envir=envstrata[[j]])
            assign("time2_inner_", time2[idx2][orderedtime], envir=envstrata[[j]])
            assign("status_inner_", status[idx2][orderedtime], envir=envstrata[[j]])
            
            if(!is.null(weight))
              assign("weight_inner_", weight[idx2][orderedtime], envir=envstrata[[j]])

            if(nchar(initeq[i] > 0)){
              eval(parse(text=initeq[i]), envir=envstrata[[j]])
            }
                        
          } else { # with random effect with strata
            stop("reached invalid code region")
            #stop("random effect in Cox regression models is not supported")
          }
        }
        
      } else { # without strata

        if(length(vec_ranvar)==0){

          envstrata[[1]] <- new.env()
          orderedtime <- order(time2[idx])
          for(k in 1:length(vec_variable)){
            assign(vec_variable[k], get(vec_variable[k], envir=data)[idx][orderedtime], envir=envstrata[[1]])
          }

          assign("time1_inner_", time1[idx][orderedtime], envir=envstrata[[1]])
          assign("time2_inner_", time2[idx][orderedtime], envir=envstrata[[1]])
          assign("status_inner_", status[idx][orderedtime], envir=envstrata[[1]])
          if(!is.null(weight))
            assign("weight_inner_", weight[idx][orderedtime], envir=envstrata[[1]])
          
          ## assign null model parameter
          ## null parameter value "nullpara" can be changed
          ## by assigning other values to nullpar in "initeq"
          assign("nullpara", rep(0, length(vec_parameter)), envir=envstrata[[1]])
          if(nchar(initeq[i]) > 0){
            eval(parse(text=initeq[i]), envir=envstrata[[1]])
          }
          
        } else { # with random effect without strata
          stop("random effect in Cox regression models is not supported")
        }
      }
      
      envs[[i]] <- envstrata
      
    ## other likelihood models      
    } else {
      
      if(!is.null(lst_option[[i]]$"strata")){
        warning("strata option is ignored because supported only in Cox regression model")
      }
      
      if(!is.null(lst_option[[i]]$"subset")){
        idx <- index & eval(parse(text=lst_option[[i]]$"subset"), envir=data)
        
      } else {
        idx <- index
      }

      if(length(vec_ranvar)==0){ # no random effect
        
        envs[[i]] <- list(new.env())

        for(j in 1:length(vec_variable)){
          assign(vec_variable[j], get(vec_variable[j], envir=data)[idx], envs[[i]][[1]])
        }
        
        ## assign null model parameter
        ## null parameter value "nullpara" can be changed
        ## by assigning other values to nullpar in "initeq"
        assign("nullpara", rep(0, length(vec_parameter)), envir=envs[[i]][[1]])
        if(!is.null(weight))
          assign("weight_inner_", weight[idx], envir=envs[[i]][[1]])
        
        if(nchar(initeq[i] > 0)){
          eval(parse(text=initeq[i]), envir=envs[[i]][[1]])
        }
        
      } else { # with random effect
        ## will be implemented later
      }
    }

    ## some work for every model
    if(nchar(itereq[i]) > 0)
      itereq[i] <- parse(text=itereq[i]) # better to prepare another parsed variable
  }

  result <- NULL

  ## Overwrite initial value by "initeq" specification from the first env
  for(i in 1:length(vec_parameter)){
    if(exists(vec_parameter[i], mode="numeric", envir=GetFirstEnvironment(envs)))
      init[i] <- get(vec_parameter[i], envir=GetFirstEnvironment(envs))
  }

  nulllik <- LogLikelihood(init=get("nullpara", envir=GetFirstEnvironment(envs)),
                           num_nmodel, envs, weight, itereq,
                           lst_psdresponse, lst_psdrespvar, lst_psdparam, lst_option,
                           vec_ranvar, vec_rparameter, vec_randens, lst_psdrparam,
                           vec_parameter, lst_parameter)

  result <- list(n=nsubject, call=Call, modelformula=modeleq,
                 parameters=vec_parameter, variables=vec_variable)
      
  class(result) <- "epifit"
  attr(result, "df") <- length(vec_parameter)
    
  if(opt=="newrap"){
    
    ans <- nlm(f=LogLikelihood, p=init, iterlim=maxiter, hessian=TRUE, gradtol=tol1,
               steptol=tol2, print.level=verbatim,
               num_nmodel=num_nmodel, envs=envs, weight=weight, itereq=itereq,
               lst_psdresponse=lst_psdresponse, lst_psdrespvar=lst_psdrespvar,
               lst_psdparam=lst_psdparam, lst_option=lst_option,
               vec_ranvar=vec_ranvar, vec_rparameter=vec_rparameter,
               vec_randens=vec_randens, lst_psdrparam=lst_psdrparam,
               vec_parameter=vec_parameter, lst_parameter=lst_parameter)
    
    result <- MakeResultFromNlm(result, ans, nulllik)
    
  } else {
    
    ans <- optim(par=init, fn=LogLikelihood, method=opt, hessian=TRUE,
                 control=list(reltol=tol1, maxit=maxiter, trace=as.integer(verbatim==2),
                   abstol=tol2),
                 num_nmodel=num_nmodel, envs=envs, weight=weight, itereq=itereq,
                 lst_psdresponse=lst_psdresponse, lst_psdrespvar=lst_psdrespvar,
                 lst_psdparam=lst_psdparam, lst_option=lst_option,
                 vec_ranvar=vec_ranvar, vec_rparameter=vec_rparameter,
                 vec_randens=vec_randens, lst_psdrparam=lst_psdrparam,
                 vec_parameter=vec_parameter, lst_parameter=lst_parameter)
    
    result <- MakeResultFromOptim(result, ans, nulllik)
    
  }

  return(result)
}
