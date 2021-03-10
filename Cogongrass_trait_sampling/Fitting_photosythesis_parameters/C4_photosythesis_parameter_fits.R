##'
##' This is an adaptation of the existing photosynthesis module to c4 photosynthesis using 
##' Xiaohui (aka Sunny's) photosynthesis fits. I belive the citation is
##' 
##' Feng, X. & Dietze, M.C. (2013). Scale dependence in the effects of leaf ecophysiological traits on photosynthesis: Bayesian parameterization of photosynthesis models. New Phytol., 200, 1132–1144.
##'   
##' @name fitA
##' @title fitA
##' @author Mike Dietze
##' @author Xiaohui Feng
##' @author Tempest McCabe
##' @export
##' 
##' @param flux.data  data.frame of Licor data, concatenated by rows, and with a leading column 'fname' that is used to count the number of curves and match to covariates
##' @param cov.data   data.frame of covariate data. Column names used in formulas
##' @param model      list including at least 6 components: the fixed effects model for alpha (a.fixed) and Vcmax (V.fixed), the random effects for these (a.random, V.random), the variable used to match the gas-exchange and covariate data (match), and the number of MCMC interations (n.iter). Additional optional arguments: TPU = TRUE turns on TPU limitation; Temp == 'Bernacchi01' turns on the Bernacchi et al 2001 temperature correction. If this is turned on all parameters are estimated for 25C, otherwise no temperature correction is applied. Setting Temp = 'June2004' will turn on the June et al 2004 Funct Plant Biol temperature correction to Jmax. Note: these two corrections are not mutually exclusive, you can set Temp = c('June2004','Bernacchi2001')
##' 
##' Right now the fixed effects are specified as a string using the standard R lm formula syntax, but without the LHS variable (e.g. '~ SLA + chl + SLA:chl'). The tilde is optional. For random effects, the two options right now are just 'leaf' for leaf-level random effects and NULL. 'model' has a default that sets all effects to NULL (fit one curve to all data) and n.iter=1000.
##' 
fitA_c4 <- function(flux.data, cov.data = NULL, model = NULL) {
  
  
  library(rjags)
  
  if (is.null(model)) {
    model <- list(a.fixed = NULL, a.random = NULL, V.fixed = NULL, V.random = NULL, 
                  n.iter = 25000, match = "fname")
  }
  out.variables <- c("r","vmax","alpha", "k","pmean", "pA")
  
  a.fixed <- model$a.fixed
  a.random <- model$a.random
  V.fixed <- model$V.fixed
  V.random <- model$V.random
  if (is.null(model$match)) {
    model$match <- "fname"
  }
  
  dat <- as.data.frame(flux.data)
  
  id         <- dat[, model$match]
  n.curves   <- length(unique(id))
  curve.id   <- as.numeric(as.factor(id))
  curve.code <- tapply(as.character(id), curve.id, unique)
  
  ## Vcmax design matrix
  if (is.null(V.fixed)) {
    XV <- NULL
  } else {
    break
    if (is.null(cov.data)) {
      print("Vcmax formula provided but covariate data is absent:", V.fixed)
    }
    if (length(grep("~", V.fixed)) == 0) {
      V.fixed <- paste("~", V.fixed)
    }
    XV      <- with(cov.data, model.matrix(formula(V.fixed)))
    XV.cols <- colnames(XV)
    XV.cols <- XV.cols[XV.cols != "(Intercept)"]
    XV      <- as.matrix(XV[, XV.cols])
    colnames(XV) <- XV.cols
    Vcenter <- apply(XV, 2, mean, na.rm = TRUE)
    XV      <- t(t(XV) - Vcenter)
  }
  
  ## alpha design matrix
  if (is.null(a.fixed)) {
    Xa <- NULL
  } else {
    break
    if (is.null(cov.data)) {
      print("alpha formula provided but covariate data is absent:", a.fixed)
    }
    a.fixed <- ifelse(length(grep("~", a.fixed)) == 0, paste("~", a.fixed), a.fixed)
    Xa      <- with(cov.data, model.matrix(formula(a.fixed)))
    Xa      <- as.matrix(Xa[, -which(colnames(Xa) == "(Intercept)")])
    acenter <- apply(Xa, 2, mean, na.rm = TRUE)
    Xa      <- t(t(Xa) - acenter)
  }
  
  ## C4 model
  
  my.model  <-
   "model{

    ## Priors
    alpha ~ dlnorm(-3.21,3.7) 	    	## initial slope of photosynthesis light response
    vmax~dlnorm(3,3)                  ## maximum rubisco capacity
    r ~ dlnorm(-0.2,2.8)              ## leaf respiration
    k ~ dlnorm(11.5, 2.8)             ## initial slope of photosynthetic CO2 response
    tau ~ dgamma(0.1,0.1)
    
    for(i in 1:n  ){                 ## process model

      al[i]<- alpha*q[i]                                  ## light limited without covariates
      ac[i]<-k*pi[i]/100000                               ## CO2 limited without covariates
      ae[i]<-vmax                                         ## rubisco limited without covariates
      pmean[i] <- min(min(al[i],ac[i]),ae[i]) - r
      an[i] ~ dnorm(pmean[i],tau)                         ## likelihood
      pA[i] ~ dnorm(pmean[i],tau)                         ## prediction
    }
  }
  "

## prep data
sel <- seq_len(nrow(dat))  #which(dat$spp == s)
if (!any(names(dat) == "Tleaf")) {
  dat$Tleaf <- rep(25 + 273.15, nrow(dat))  ## if leaf temperature is absent, assume 25C
}
mydat <- list(an = dat$Photo[sel], 
              pi = dat$Ci[sel], 
              q = dat$PARi[sel],
              #T = dat$Tleaf, # I 
              n = length(sel), 
              Kc = 46, 
              Ko = 22000, 
              po = 21000, 
              rep = curve.id, 
              nrep = n.curves)
#  Kc<-46                          ## Michaelis constant CO2 (Pa)
#  Ko<-33000                       ## Michaelis constant O2  (Pa)
#  po<-21000                       ## partial pressure of O2  (Pa)

## TPU Limitation
if ("TPU" %in% names(model)) {
  if (model$TPU == TRUE) {
    my.model <- gsub(pattern = "#TPU", " ", my.model)
    out.variables <- c(out.variables, "tpu")
  }
}

## Temperature scaling
Vformula <- NULL
if ("Temp" %in% names(model)) {
  if ("Bernacchi01" %in% model$Temp) {
    my.model <- gsub(pattern = "##B01", " ", my.model)
  }
  if ("June2004" %in% model$Temp) {
    my.model <- gsub(pattern = "##J04", " ", my.model)
  }
}

## VCmax Formulas
Vformula <- NULL
if ("leaf" %in% V.random) {
  Vformula <- " + Vleaf[rep[i]]"
  my.model <- gsub(pattern = "#RLEAF.V", " ", my.model)
  out.variables <- c(out.variables, "tau.Vleaf")
}

if (!is.null(XV)) {
  Vnames <- gsub(" ", "_", colnames(XV))
  Vformula <- paste(Vformula,
                    paste0("+ betaV", Vnames, "*XV[rep[i],", seq_along(XV), "]", collapse = " "))
  Vpriors <- paste0("     betaV", Vnames, "~dnorm(0,0.001)", collapse = "\n")
  my.model <- sub(pattern = "## Vcmax BETAS", Vpriors, my.model)
  mydat[["XV"]] <- XV
  out.variables <- c(out.variables, paste0("betaV", Vnames))
}
if (!is.null(Vformula)) {
  my.model <- sub(pattern = "#VFORMULA", Vformula, my.model)
} 

## alpha Formulas
Aformula <- NULL
if ("leaf" %in% a.random) {
  Aformula <- " + Aleaf[rep[i]]"
  my.model <- gsub(pattern = "#RLEAF.A", "", my.model)
  out.variables <- c(out.variables, "tau.Aleaf")
}

if (!is.null(Xa)) {
  Anames <- gsub(" ", "_", colnames(Xa))
  Aformula <- paste(Aformula, paste0("+ betaA", Anames, "*Xa[rep[i],", 1:ncol(Xa), 
                                     "]", collapse = " "))
  apriors <- paste0("betaA", Anames, "~dnorm(0,0.001)", collapse = "\n")
  my.model <- sub(pattern = "## alpha BETAs", apriors, my.model)
  mydat[["Xa"]] <- Xa
  out.variables <- c(out.variables, paste0("betaA", Anames))
}
if (!is.null(Aformula)) {
  my.model <- sub(pattern = "#AFORMULA", Aformula, my.model)
}

## Define initial conditions
init <- list()
init[[1]]<-list(r=0.8, vmax=30,alpha=0.03, tau=10, k=0.7*100000)   ## ,tau.Vleaf=300,tau.Kleaf=1e-10, beta1=4, beta2=1,beta5=3,tau.Vmon=10
init[[2]]<-list(r=1, vmax=20, alpha=0.07, tau=20, k=0.8*100000)    ## ,tau.Vleaf=200,tau.Kleaf=2e-9,beta1=1,beta2=1,beta5=-1,tau.Vmon=20
init[[3]]<-list(r=2, vmax=15,alpha=0.06, tau=20, k=0.2*1000000)    ## ,tau.Vleaf=100,tau.Kleaf=3e-8,beta1=1,beta2=2,beta5=2,tau.Vmon=3

mc3 <- jags.model(file = textConnection(my.model), data = mydat, inits = init, n.chains = 3)

mc3.out <- coda.samples(model = mc3, variable.names = out.variables, n.iter = model$n.iter)

## split output
out         <- list(params = NULL, predict = NULL, model = my.model)
mfit        <- as.matrix(mc3.out, chains = TRUE)
pred.cols   <- union(grep("pA", colnames(mfit)), grep("pmean", colnames(mfit)))
chain.col   <- which(colnames(mfit) == "CHAIN")
out$predict <- mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
out$params  <- mat2mcmc.list(mfit[, -pred.cols])
return(out)
} # fitA


##' @name read_Licor
##' @title read_Licor
##' 
##' @author Mike Dietze
##' @export
##' 
##' @param filename  name of the file to read
##' @param sep       file delimiter. defaults to tab
##' @param ...       optional arguements forwarded to read.table
read_Licor <- function(filename, sep = "\t", ...) {
  fbase <- sub(".txt", "", tail(unlist(strsplit(filename, "/")), n = 1))
  print(fbase)
  full <- readLines(filename)
  ## remove meta=data
  start <- grep(pattern = "OPEN", full)
  skip <- grep(pattern = "STARTOFDATA", full)
  for (i in rev(seq_along(start))) {
    full <- full[-(start[i]:(skip[i] + 1 * (i > 1)))]  # +1 is to deal with second header
  }
  full <- full[grep("\t", full)]  ## skip timestamp lines
  dat <- read.table(textConnection(full), header = TRUE, blank.lines.skip = TRUE,
                    sep = sep, ...)
  fname <- rep(fbase, nrow(dat))
  dat <- as.data.frame(cbind(fname, dat))
  return(dat)
} # read_Licor


mat2mcmc.list <- function(w) {
  temp <- list()
  chain.col <- which(colnames(w) == "CHAIN")
  for (i in unique(w[, "CHAIN"])) {
    temp[[i]] <- as.mcmc(w[w[, "CHAIN"] == i, -chain.col])
  }
  return(as.mcmc.list(temp))
} # mat2mcmc.list

