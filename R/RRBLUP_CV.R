#' @title Ridge Regression, Leave-one-out and Leave-D-out Cross-validation
#' 
#' 
#' @description This function solves univariate linear mixed models by likelihood methods. The optimization methods is based on Efficient Mixed Model Association (Kang et al. 2008).
#' 
#' 
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame containing the variables in the model.
#' @param Z SNP markers matrix with animals in rows and markers in columns.
#' @param CV (logical) if TRUE the Leave-one-out Cross-validation is performed, default TRUE.
#' @param folds is a vector with validation sets used to perform k-folds cross-validation. Should be NULL or a numeric vector. If it is NULL, LOOCV will be performed.
#' @param weights is an optional vector to allow for heterogeneous error variance: \eqn{Var[\varepsilon_i] = R_i \sigma^2_e}. Should be NULL or a numeric vector. The length of weights must be equal to the number of individuals. 
#' 
#' 
#' @return A list with genetic and residual variances; a vector with BLUPs for random effects; a vector for BLUEs of fixed effects; the log-likelihood. AIC and BIC statistics. If CV is TRUE, a vector with BLUPs for random effects using LOOCV method.
#' 
#' @references Hyun Min Kang, Noah A. Zaitlen, Claire M. Wade, Andrew Kirby, David Heckerman, Mark J. Daly and Eleazar Eskin, 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723. doi:10.1534/genetics.107.080101.
#' 
#' @references Gianola D, Schon C-C, 2016. Cross-Validation Without Doing Cross-Validation in Genome-Enabled Prediction. G3: Genes|Genomes|Genetics. 6(10):3107-3128. doi:10.1534/g3.116.033381.
#' 
#' 
#' @examples
#' ## Not to run ##
#' 
#' ## RRBLUP_CV(Phen ~ Effect, data=Data, K=M)
#'
#' ## End(Not run)
#' 
#' @export RRBLUP_CV
#' 
#' @import stats
#' @useDynLib ezCV
#' @importFrom Rcpp evalCpp
RRBLUP_CV <- function(formula, data, Z, CV=TRUE, folds=NULL, weights=NULL){
  cat("\n")
  centerText <- function() {
    width <- getOption("width")
    A <- ("                                                       ._____.    \n")
    B <- ("    _/////_            Fernando Brito Lopes           _|_____|_   \n")
    C <- ("   (' o o ')     Animal Scientist (Zootechnician)     (' o o ')   \n")
    D <- ("__ooO_(_)_Ooo_____ Animal Breeding and Genetics _____ooO_(_)_Ooo__\n")
    E <- ("                    e-mail: <camult@gmail.com>                    \n")
    ws <- rep(" ", floor((width - nchar(A))/2))
    cat(ws,A,ws,B,ws,C,ws,D,ws,E,ws,sep = "")
  }
  centerText()
  cat("\n")
  form <- as.formula(formula)
  phenName <- as.character(form[[2]])
  Y <- as.numeric(data[, phenName])
  not.miss <- which(!is.na(Y))
  resid <- rep(NA, length(Y))
  if (length(not.miss) < length(Y)) {
    data <- data[not.miss, ]
    Y <- Y[not.miss]
    if (!is.null(weights)) {
      weights <- weights[not.miss]
    }
  }
  n <- length(Y)
  X <- model.matrix(form, data = data)
  not.miss.gid <- as.character(unique(data[, 1]))
  gid <- rownames(Z)
  mid <- colnames(Z)
  ix.pheno <- match(not.miss.gid, gid)
  miss.pheno.gid <- which(is.na(ix.pheno))
  if (length(miss.pheno.gid) > 0) { stop(paste("The following lines have phenotypes but no genotypes:", paste(not.miss.gid[miss.pheno.gid], collapse = " "))) }
  miss.gid <- setdiff(gid, not.miss.gid)
  ix <- c(ix.pheno, match(miss.gid, gid))
  Z <- Z[ix,]
  if (!is.null(weights)) {
    sqrt.weights <- sqrt(weights)
    X <- X/sqrt.weights
    Y <- Y/sqrt.weights
    Z <- Z/sqrt.weights
  }
  Y <- as.matrix(Y)
  cat("Running Ridge Regression Model...\n")
  UVM <- RRBLUP(Y, X, Z)
  rownames(UVM$u) <- mid
  GEBV <- Z %*% UVM$u
  resid[not.miss] <- Y - X %*% UVM$beta - GEBV
  pred <- as.matrix(GEBV + as.numeric(colMeans(X)%*%UVM$beta))
  yHat <- as.matrix(Y - (X%*%UVM$beta))
  if(isTRUE(CV)){
    lambda <- UVM$Ve/UVM$Vu
    if(is.null(folds)){
      cat("Performing Leave-One-Out Cross-Validation Without Doing Cross-Validation...\n")
      yHat.CV <- LOOrrDG(yHat, Z, as.matrix(UVM$u), lambda)
    } else {
      stop("Leave-D-Out CV Method is Not Implemented yet...")
      #cat("Performing Leave-d-Out Cross-Validation Without Doing Cross-Validation...\n")
      #yHat.CV <- LOOrrDG(yHat, Z, lambda, folds, as.matrix(UVM$u))
    }
    return(list(Vg = UVM$Vu,
                Ve = UVM$Ve, 
                beta = UVM$beta,
                LL = UVM$LL,
                AIC = as.vector((-2 * UVM$LL ) + (2 * dim(X)[1])),
                BIC = as.vector((-2 * UVM$LL ) + (log(dim(as.matrix(Y))[2]) * dim(X)[1])),
                SNP_Effect = UVM$u,
                g = GEBV, 
                resid = resid,
                pred = pred,
                g.CV=yHat.CV))
  } else {
    return(list(Vg = UVM$Vu,
                Ve = UVM$Ve, 
                beta = UVM$beta,
                LL = UVM$LL,
                AIC = as.vector((-2 * UVM$LL ) + ( 2 * dim(X)[1])),
                BIC = as.vector((-2 * UVM$LL ) + ( log(dim(as.matrix(Y))[2]) * dim(X)[1])),
                SNP_Effect = UVM$u,
                g = GEBV, 
                resid = resid,
                pred = pred))
  }
}
