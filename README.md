# `GBLUP_CV`: GBLUP and Leave-one-out Cross-validation.

## Description


 This function solves univariate linear mixed models by likelihood methods. The optimization methods is based on Efficient Mixed Model Association (Kang et al. 2008).


## Usage

```r
GBLUP_CV(formula, data, K, CV = TRUE, folds = NULL, weights = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
```formula```     |     an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
```data```     |     an optional data frame containing the variables in the model.
```K```     |     Covariance matrix for random effects.
```CV```     |     (logical) if TRUE the Leave-one-out Cross-validation is performed, default TRUE.
```folds```     |     is a vector with validation sets used to perform k-folds cross-validation. Should be NULL or a numeric vector. If it is NULL, LOOCV will be performed.
```weights```     |     is an optional vector to allow for heterogeneous error variance: $Var[\varepsilon_i] = R_i \sigma^2_e$ . Should be NULL or a numeric vector. The length of weights must be equal to the number of individuals.

## Value


 A list with genetic and residual variances; a vector with BLUPs for random effects; a vector for BLUEs of fixed effects; the log-likelihood. AIC and BIC statistics. If CV is TRUE, a vector with BLUPs for random effects using LOOCV method.


## References


 Hyun Min Kang, Noah A. Zaitlen, Claire M. Wade, Andrew Kirby, David Heckerman, Mark J. Daly and Eleazar Eskin, 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723. doi:10.1534/genetics.107.080101.
 
 Gianola D, Schon C-C, 2016. Cross-Validation Without Doing Cross-Validation in Genome-Enabled Prediction. G3: Genes|Genomes|Genetics. 6(10):3107-3128. doi:10.1534/g3.116.033381.


## Examples

```r 
 ## Not to run ##
 
 ## GBLUP_CV(Phen ~ Effect, data=Data, K=G)
 
 ## End(Not run)
 
 ``` 

# `cvBGBLUP`: Cross-Validation Without Doing Cross-Validation on Gibbs Sampling

## Description


 Performs a Leave-One-Out Cross-Validation Without Doing Cross-Validation using Gibbs Sampling results.


## Usage

```r
cvBGBLUP(n, Y, g, K, Vg, Ve)
```


## Arguments

Argument      |Description
------------- |----------------
```n```     |     is the number of observations. If length(n) > 1, the length is taken to be the number required.
```Y```     |     is a phenotypic matrix with n rows and 1 column.
```g```     |     is vector with the BLUP solution for the genetic values.
```K```     |     is a relationship matrix with m rows and m columns.
```Vg```     |     is the posterior mean of the genetic variance.
```Ve```     |     is the posterior mean of the residual variance.

## Examples

```r 
 ## Not to run ##
 
 ## cvBGBLUP(n=nSamp, Y=y, g=gHat, K=G, Vg=gVAR, Ve=eVAR)
 
 ## End(Not run)
 
 ``` 

# `RRBLUP_CV`: Ridge Regression, Leave-one-out and Leave-D-out Cross-validation

## Description


 This function solves univariate linear mixed models by likelihood methods. The optimization methods is based on Efficient Mixed Model Association (Kang et al. 2008).


## Usage

```r
RRBLUP_CV(formula, data, Z, CV = TRUE, folds = NULL, weights = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
```formula```     |     an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
```data```     |     an optional data frame containing the variables in the model.
```Z```     |     SNP markers matrix with animals in rows and markers in columns.
```CV```     |     (logical) if TRUE the Leave-one-out Cross-validation is performed, default TRUE.
```folds```     |     is a vector with validation sets used to perform k-folds cross-validation. Should be NULL or a numeric vector. If it is NULL, LOOCV will be performed.
```weights```     |     is an optional vector to allow for heterogeneous error variance: $Var[\varepsilon_i] = R_i \sigma^2_e$ . Should be NULL or a numeric vector. The length of weights must be equal to the number of individuals.

## Value


 A list with genetic and residual variances; a vector with BLUPs for random effects; a vector for BLUEs of fixed effects; the log-likelihood. AIC and BIC statistics. If CV is TRUE, a vector with BLUPs for random effects using LOOCV method.


## References


 Hyun Min Kang, Noah A. Zaitlen, Claire M. Wade, Andrew Kirby, David Heckerman, Mark J. Daly and Eleazar Eskin, 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723. doi:10.1534/genetics.107.080101.
 
 Gianola D, Schon C-C, 2016. Cross-Validation Without Doing Cross-Validation in Genome-Enabled Prediction. G3: Genes|Genomes|Genetics. 6(10):3107-3128. doi:10.1534/g3.116.033381.


## Examples

```r 
 ## Not to run ##
 
 ## RRBLUP_CV(Phen ~ Effect, data=Data, K=M)
 
 ## End(Not run)
 
 ``` 

# `cvBayes`: Cross-Validation Without Doing Cross-Validation on Gibbs Sampling

## Description


 Performs a Leave-One-Out Cross-Validation Without Doing Cross-Validation using Gibbs Sampling results.


## Usage

```r
cvBayes(Y, B, X, varE)
```


## Arguments

Argument      |Description
------------- |----------------
```Y```     |     is a phenotypic matrix with n rows and 1 column.
```B```     |     is a matrix with the posterior distribution of markers effect.
```X```     |     is a matrix with markers in columns and animals in rows.
```varE```     |     is a vector with the posterior distribution of residual variance.

## Examples

```r 
 ## Not to run ##
 
 ## cvBayes(Y=y, B=Effect, X=X, varE=BA_varE)
 
 ## End(Not run)
 
 ``` 

# `MME`: Mixed Model Equation

## Description


 Solves a univariate mixed model of form $y=X\beta+Zu+e$ .


## Usage

```r
MME(y, X, Z, K)
```


## Arguments

Argument      |Description
------------- |----------------
```y```     |     a matrix with n rows and 1 column.
```X```     |     a matrix with n rows and x columns.
```Z```     |     a matrix with n rows and m columns.
```K```     |     a matrix with m rows and m columns.

## Examples

```r 
 ## Not to run ##
 
 ## MME(Y, X, Z, K)
 
 ## End(Not run)
 
 ``` 

# `LOOCV_DG`: Cross-Validation Without Doing Cross-Validation

## Description


 Performs a Leave-One-Out Cross-Validation Without Doing Cross-Validation using a relationship matrix.


## Usage

```r
LOOCV_DG(y, K, lambda, g)
```


## Arguments

Argument      |Description
------------- |----------------
```y```     |     is a phenotypic matrix with n rows and 1 column.
```K```     |     is a relationship matrix with m rows and m columns.
```lambda```     |     is the ratio bwtween residual and genetic variance.
```g```     |     is the BLUP solution for the genetic values.

## Examples

```r 
 ## Not to run ##
 
 ## LOOCV_DG(Y, G, lambda)
 
 ## End(Not run)
 
 ``` 

# `kCV_DG`: K-folds Cross-Validation Without Doing Cross-Validation

## Description


 Performs a K-folds Cross-Validation Without Doing Cross-Validation using a relationship matrix.


## Usage

```r
kCV_DG(Y, K, lambda, folds, g)
```


## Arguments

Argument      |Description
------------- |----------------
```Y```     |     is a phenotypic matrix with n rows and 1 column.
```K```     |     is a relationship matrix with m rows and m columns.
```lambda```     |     is the ratio bwtween residual and genetic variance.
```folds```     |     is a vector with validation sets used to perform k-folds cross-validation.
```g```     |     is the BLUP solution for the genetic values.

## Examples

```r 
 ## Not to run ##
 
 ## kCV_DG(Y, K, lambda, folds)
 
 ## End(Not run)
 
 ``` 

# `RRBLUP`: Ridge Regression BLUP

## Description


 Solves a univariate mixed model of form $y=X\beta+Mu+e$ 


## Usage

```r
RRBLUP(y, X, M)
```


## Arguments

Argument      |Description
------------- |----------------
```y```     |     a matrix with n rows and 1 column
```X```     |     a matrix with n rows and x columns
```M```     |     a matrix with n rows and m columns

## Examples

```r 
 ## Not to run ##
 
 ## RRBLUP(y, X, M)
 
 ## End(Not run)
 
 ``` 

# `LOOrrDG`: rrBLUP Cross-Validation Without Doing Cross-Validation

## Description


 Performs a Leave-One-Out Cross-Validation Without Doing Cross-Validation using a relationship matrix.


## Usage

```r
LOOrrDG(Y, X, B, lambda)
```


## Arguments

Argument      |Description
------------- |----------------
```Y```     |     is a phenotypic matrix with n rows and 1 column.
```X```     |     is a matrix with markers in columns and animals in rows.
```B```     |     is matrix with markers effect.
```lambda```     |     is the ratio bwtween residual and genetic variance.

## Examples

```r 
 ## Not to run ##
 
 ## LOOrrDG(Y, X, B, lambda)
 
 ## End(Not run)
 
 ``` 

