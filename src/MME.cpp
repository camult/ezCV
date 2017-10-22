// Functions from AlphaSimR
#include <RcppArmadillo.h>
#include <string>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Note: Fortran compiler appends '_' to subroutine name
// See http://www.netlib.org/lapack/explore-html/ for description of args
extern "C" void dsyevr_(char* JOBZ, char* RANGE, char* UPLO, long long int* N, double* A, long long int* LDA, double* VL,
                       double* VU, long long int* IL, long long int* IU, double* ABSTOL, long long int* M, double* W, double* Z,
                       long long int* LDZ, long long int* ISUPPZ, double* WORK, long long int* LWORK, long long int* IWORK, 
                       long long int* LIWORK, long long int* INFO);

// Replacement for Armadillo's eig_sym
// Fixes an error with decompisition of large matrices on Eddie
// If calcVec = false, eigvec is not used
// It would be better to template this function
int eigen2(arma::vec& eigval, arma::mat& eigvec, arma::mat X, 
           bool calcVec = true){
  char JOBZ;
  if(calcVec){
    JOBZ = 'V';
  }else{
    JOBZ = 'N';
  }
  char RANGE = 'A';
  char UPLO = 'L';
  long long int N = X.n_rows;
  // A = X
  long long int LDA = N;
  double VL = 0.0;
  double VU = 0.0;
  long long int IL;
  long long int IU;
  double ABSTOL = 0.0;
  long long int M = N;
  // W=eigval
  // Z=eigvec
  long long int LDZ = N;
  arma::Col<long long int> ISUPPZ(2*M);
  // WORK length to be determined
  double tmpWORK;
  long long int LWORK = -1; // To be calculated
  // IWORK length to be determined
  long long int tmpIWORK;
  long long int LIWORK = -1; // To be calculated
  long long int INFO;
  // Calculate LWORK and LIWORK
  dsyevr_(&JOBZ,&RANGE,&UPLO,&N,&*X.begin(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,&*eigval.begin(),
          &*eigvec.begin(),&LDZ,&*ISUPPZ.begin(),&tmpWORK,&LWORK,&tmpIWORK,&LIWORK,&INFO);
  LWORK = (long long int) tmpWORK;
  LIWORK = tmpIWORK;
  // Allocate WORK and IWORK
  arma::vec WORK(LWORK);
  arma::Col<long long int> IWORK(LIWORK);
  // Perform decomposition
  dsyevr_(&JOBZ,&RANGE,&UPLO,&N,&*X.begin(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,&*eigval.begin(),
          &*eigvec.begin(),&LDZ,&*ISUPPZ.begin(),&*WORK.begin(),&LWORK,&*IWORK.begin(),&LIWORK,&INFO);
  return INFO; // Return error code
}

// Objective function for REML using the EMMA algorithm
Rcpp::List objREML(double param, Rcpp::List args){
  double df = args["df"];
  arma::vec eta = args["eta"];
  arma::vec lambda = args["lambda"];
  double value = df * log(sum(eta%eta/(lambda+param)));
  value += sum(log(lambda+param));
  return Rcpp::List::create(Rcpp::Named("objective") = value,
                            Rcpp::Named("output") = 0);
}

/*
* Modified version of Brent's method for numerical optimization
* objective is a pointer to a function returning a List
*   objective takes a double and a List as arguments
*   optimization is for the double
*   the returned list includes "objective" and list of additional output
* args is a list which is passed to the objective
* l is the lower bound of optimization
* u is the upper bound of optimization
* maxIter is the maximum number of iterations
*   the best obtained value is returned when maxIter is reached
* if firstLast is true, l and u are evaluated if convergence isn't reached
* maximize searches for a maximum value
* eps controls the tolerance for convergence
*/
Rcpp::List optimize(Rcpp::List (*objective)(double, Rcpp::List), 
                    Rcpp::List args, double l, double u, int maxIter=1000, 
                    bool maximize=false, bool evalU = false, bool evalL = false, 
                    double eps=1.0e-9){
  double MACHEPS_SQRT = sqrt(2.2204460492503131e-016);
  double c = (3.0 - sqrt(5.0)) / 2.0;
  double x = l+c*(u - l);
  double v = x;
  double w = x;
  double e = 0.0;
  double lInt = l;
  double uInt = u;
  List fOut = objective(x, args);
  List output = fOut["output"];
  double fx = fOut["objective"];
  if(maximize) fx = -fx;
  double fv = fx;
  double fw = fx;
  
  int numiter = 0;
  
  while(true){
    double m = 0.5*(l+u);
    double tol = MACHEPS_SQRT*std::abs(x)+eps;
    double tol2 = 2.0*tol;
    
    // Check the stopping criterion
    if(std::abs(x-m) <= (tol2-(0.5*(u-l)))){
      break;
    }
    // Check maximum iterations
    if (++numiter > maxIter){
      break;
    }
    
    double p = 0.0, q = 0.0, r = 0.0, d = 0.0, z = 0.0;
    
    if (std::abs(e) > tol){
      // Fit parabola
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q-(x-w)*r;
      q = 2.0*(q-r);
      if (q>0.0){
        p = -p;
      }else{
        q = -q;
      }
      r = e;
      e = d;
    }
    
    if((std::abs(p) < std::abs(0.5*q*r))&&
       (p < (q*(l-x)))&&
       (p < (q*(u-x)))){
      // Parabolic interpolation step
      d = p/q;
      z = x+d;
      // objective must not be evaluated too close to l or u
      if(((z-l) < tol2) || ((u-z) < tol2)){
        d = (x < m) ? tol : -tol;
      }
    }else{
      // Golden section step
      e = (x < m) ? (u - x) : (l - x);
      d = c*e;
    }
    
    // objective must not be evaluated too close to x
    if(std::abs(d) >= tol){
      z = x+d;
    }else if(d > 0.0){
      z = x+tol;
    }else{
      z = x-tol;
    }
    fOut = objective(z, args);
    double funcu = fOut["objective"];
    if(maximize) funcu = -funcu;
    
    // Update
    if(funcu <= fx){
      if(z < x){
        u = x;
      }else{
        l = x;
      }
      
      v = w; 
      fv = fw;
      w = x; 
      fw = fx;
      x = z;
      output = fOut["output"];
      fx = funcu;
    }else{
      if (z < x){
        l = z;
      }else{
        u = z;
      }
      
      if ((funcu <= fw) || (w == x)){
        v = w;
        fv = fw;
        w = z;
        fw = funcu;
      }else if((funcu <= fv) || (v == x) || (v == w)){
        v = z;
        fv = funcu;
      }
    }
  }
  if(evalU){
    fOut = objective(uInt, args);
    double funcu = fOut["objective"];
    if(maximize) funcu = -funcu;
    if(funcu<fx){
      fx = funcu;
      output = fOut["output"];
      x = uInt;
    }
  }
  if(evalL){
    fOut = objective(lInt, args);
    double funcu = fOut["objective"];
    if(maximize) funcu = -funcu;
    if(funcu<fx){
      fx = funcu;
      output = fOut["output"];
      x = uInt;
    }
  }
  if(maximize) fx = -fx;
  return List::create(Named("parameter") = x,
                      Named("objective") = fx,
                      Named("output") = output);
}

//' @title Mixed Model Equation
//' 
//' @description Solves a univariate mixed model of form \eqn{y=X\beta+Zu+e}.
//'
//' @param y a matrix with n rows and 1 column.
//' @param X a matrix with n rows and x columns.
//' @param Z a matrix with n rows and m columns.
//' @param K a matrix with m rows and m columns.
//' 
//' @examples
//' ## Not to run ##
//' 
//' ## MME(Y, X, Z, K)
//'
//' ## End(Not run)
//' 
//' @useDynLib ezCV
//' 
//' @export MME
// [[Rcpp::export]]
Rcpp::List MME(const arma::mat& y, const arma::mat& X, 
               const arma::mat& Z, const arma::mat& K){
  int n = y.n_rows;
  int q = X.n_cols;
  double df = double(n)-double(q);
  double offset = log(double(n));
  bool invPass;
  // Construct system of equations for eigendecomposition
  arma::mat S = arma::eye(n,n) - X*inv_sympd(X.t()*X)*X.t();
  arma::mat ZK = Z*K;
  arma::mat ZKZ = ZK*Z.t();
  S = S*(ZKZ+offset*arma::eye(n,n))*S;
  // Compute eigendecomposition
  arma::vec eigval(n);
  arma::mat eigvec(n,n);
  eigen2(eigval, eigvec, S);
  // Drop eigenvalues
  eigval = eigval(arma::span(q,eigvec.n_cols-1)) - offset;
  eigvec = eigvec(arma::span(0,eigvec.n_rows-1),
                  arma::span(q,eigvec.n_cols-1));
  // Estimate variances and solve equations
  arma::vec eta = eigvec.t()*y;
  Rcpp::List optRes = optimize(*objREML,
                               Rcpp::List::create(
                                 Rcpp::Named("df")=df,
                                 Rcpp::Named("eta")=eta,
                                 Rcpp::Named("lambda")=eigval), 
                                 1.0e-10, 1.0e10);
  double delta = optRes["parameter"];
  arma::mat Hinv; 
  invPass = inv_sympd(Hinv,ZKZ+delta*arma::eye(n,n));
  if(!invPass){
    Hinv = pinv(ZKZ+delta*arma::eye(n,n));
  }
  arma::mat XHinv = X.t()*Hinv;
  arma::mat beta = solve(XHinv*X,XHinv*y);
  arma::mat u = ZK.t()*(Hinv*(y-X*beta));
  double Vu = sum(eta%eta/(eigval+delta))/df;
  double Ve = delta*Vu;
  arma::colvec PEV = arma::diagvec(inv_sympd(Z.t()*Z + ((Ve/Vu)*inv_sympd(K))))*Ve;
  double ll = -0.5*(double(optRes["objective"])+df+df*log(2*PI/df));
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("PEV")=PEV,
                            Rcpp::Named("LL")=ll);
}

//' @title Cross-Validation Without Doing Cross-Validation
//' 
//' @description Performs a Leave-One-Out Cross-Validation Without Doing Cross-Validation using a relationship matrix.
//'
//' @param y is a phenotypic matrix with n rows and 1 column.
//' @param K is a relationship matrix with m rows and m columns.
//' @param lambda is the ratio bwtween residual and genetic variance.
//' @param g is the BLUP solution for the genetic values.
//' 
//' @examples
//' ## Not to run ##
//' 
//' ## LOOCV_DG(Y, G, lambda)
//'
//' ## End(Not run)
//' 
//' @useDynLib ezCV
//'  
//' @export LOOCV_DG
// [[Rcpp::export]]
arma::vec LOOCV_DG(arma::mat& y, arma::mat& K, double lambda, arma::vec& g){
  int nIDs = K.n_rows;
  arma::vec yHat(nIDs);
  arma::mat invC = arma::inv_sympd(arma::eye(nIDs, nIDs) + arma::inv_sympd(K)*lambda);
  //arma::vec g = invC*y;
  for (int nAn = 0; nAn < nIDs; nAn++) {
    yHat[nAn] = (g[nAn]-invC(nAn,nAn)*y[nAn])/(1-invC(nAn,nAn));
  }
  return yHat;
}

//' @title K-folds Cross-Validation Without Doing Cross-Validation
//' 
//' @description Performs a K-folds Cross-Validation Without Doing Cross-Validation using a relationship matrix.
//'
//' @param Y is a phenotypic matrix with n rows and 1 column.
//' @param K is a relationship matrix with m rows and m columns.
//' @param lambda is the ratio bwtween residual and genetic variance.
//' @param folds is a vector with validation sets used to perform k-folds cross-validation.
//' @param g is the BLUP solution for the genetic values.
//' 
//' @examples
//' ## Not to run ##
//' 
//' ## kCV_DG(Y, K, lambda, folds)
//'
//' ## End(Not run)
//' 
//' @useDynLib ezCV
//'  
//' @export kCV_DG
// [[Rcpp::export]]
arma::vec kCV_DG(arma::vec& Y, arma::mat& K, double lambda, arma::vec& folds, arma::vec& g){
  int k = arma::max(folds);
  int nIDs = Y.n_rows;
  arma::vec yHat(nIDs);
  arma::mat invC = arma::inv_sympd(arma::eye(nIDs, nIDs) + arma::inv_sympd(K)*lambda);
  //arma::vec g = invC*Y;
  for(int fold = 1; fold <= k; fold++) {
    arma::uvec idx = find(folds == fold);
    arma::mat subInvC  = invC.submat(idx,idx);
    int nIdx = subInvC.n_rows;
    yHat.elem(idx) = arma::inv_sympd(arma::eye(nIdx, nIdx)-subInvC)*(g.elem(idx)-(subInvC*Y.elem(idx)));
  }
  return yHat;
}


//' @title Ridge Regression BLUP
//'
//' @description Solves a univariate mixed model of form \eqn{y=X\beta+Mu+e}
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param M a matrix with n rows and m columns
//' 
//' @examples
//' ## Not to run ##
//' 
//' ## RRBLUP(y, X, M)
//'
//' ## End(Not run)
//' 
//' @useDynLib ezCV
//' 
//' @export RRBLUP
// [[Rcpp::export]]
Rcpp::List RRBLUP(const arma::mat& y, const arma::mat& X, const arma::mat& M){
  int n = y.n_rows;
  int q = X.n_cols;
  double df = double(n)-double(q);
  double offset = log(double(n));
  bool invPass;
  arma::mat S = arma::eye(n,n) - X*inv_sympd(X.t()*X)*X.t();
  S = S*((M*M.t())+offset*arma::eye(n,n))*S;
  arma::vec eigval(n);
  arma::mat eigvec(n,n);
  eigen2(eigval, eigvec, S);
  eigval = eigval(arma::span(q,eigvec.n_cols-1)) - offset;
  eigvec = eigvec(arma::span(0,eigvec.n_rows-1),
                  arma::span(q,eigvec.n_cols-1));
  arma::vec eta = eigvec.t()*y;
  Rcpp::List optRes = optimize(*objREML,
                               Rcpp::List::create(
                                 Rcpp::Named("df")=df,
                                 Rcpp::Named("eta")=eta,
                                 Rcpp::Named("lambda")=eigval),
                                 1.0e-10, 1.0e10);
  double delta = optRes["parameter"];
  arma::mat Hinv;
  invPass = inv_sympd(Hinv,M*M.t()+delta*arma::eye(n,n));
  if(!invPass){
    Hinv = pinv(M*M.t()+delta*arma::eye(n,n));
  }
  arma::mat XHinv = X.t()*Hinv;
  arma::mat beta = solve(XHinv*X,XHinv*y);
  arma::mat u = M.t()*(Hinv*(y-X*beta));
  double Vu = sum(eta%eta/(eigval+delta))/df;
  double Ve = delta*Vu;
  double ll = -0.5*(double(optRes["objective"])+df+df*log(2*PI/df));
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=ll);
}

//' @title rrBLUP Cross-Validation Without Doing Cross-Validation
//' 
//' @description Performs a Leave-One-Out Cross-Validation Without Doing Cross-Validation using a relationship matrix.
//'
//' @param Y is a phenotypic matrix with n rows and 1 column.
//' @param X is a matrix with markers in columns and animals in rows.
//' @param B is matrix with markers effect.
//' @param lambda is the ratio bwtween residual and genetic variance.
//' 
//' @examples
//' ## Not to run ##
//' 
//' ## LOOrrDG(Y, X, B, lambda)
//'
//' ## End(Not run)
//' 
//' @useDynLib ezCV
//'  
//' @export LOOrrDG
// [[Rcpp::export]]
arma::vec LOOrrDG(arma::mat& Y, arma::mat& X, arma::mat& B, double lambda){
  int n = X.n_rows;
  int m = X.n_cols;
  arma::vec yHat(n);
  arma::vec yHati(n);
  arma::mat invC = arma::inv_sympd(X.t()*X + lambda*arma::eye(m, m));
  arma::vec ei = Y - X*B;
  arma::mat h = X*invC*X.t();
  for (int nAn = 0; nAn < n; nAn++) {
    arma::vec betai = B - (invC*trans(X.row(nAn))*ei[nAn] / (1-h(nAn,nAn)));
    yHat = X*betai;
    yHati[nAn] = yHat[nAn];
  }
  return yHati;
}


//' @title Cross-Validation Without Doing Cross-Validation on Gibbs Sampling
//' 
//' @description Performs a Leave-One-Out Cross-Validation Without Doing Cross-Validation using Gibbs Sampling results.
//'
//' @param Y is a phenotypic matrix with n rows and 1 column.
//' @param B is a matrix with the posterior distribution of markers effect.
//' @param X is a matrix with markers in columns and animals in rows.
//' @param varE is a vector with the posterior distribution of residual variance.
//' 
//' @examples
//' ## Not to run ##
//' 
//' ## cvBayes(Y=y, B=Effect, X=X, varE=BA_varE)
//'
//' ## End(Not run)
//' 
//' @useDynLib ezCV
//'  
//' @export cvBayes
// [[Rcpp::export]]
arma::mat cvBayes(arma::mat Y, arma::mat B, arma::mat X, arma::vec varE){
  int S = B.n_rows;
  int n = Y.n_rows;
  arma::mat yMAT(n,S);
  arma::mat eMAT(S,n);
  for(int s = 0; s < S; s++){
    yMAT.col(s) = Y;
  }
  for(int i = 0; i < n; i++){
    eMAT.col(i) = varE;
  }
  arma::mat ssMAT = pow(trans(yMAT - (X*trans(B))),2);
  arma::mat exMat = exp(ssMAT/(2.0*eMAT));
  arma::mat sumLOOP(S,n);
  for(int s = 0; s < S; s++){
    sumLOOP.row(s) = sum(exp(ssMAT/(2.0*varE[s])),0);
  }
  arma::mat w = exMat/sumLOOP;
  arma::mat yHat(n,1);
  for(int i = 0; i < n; i++){
    yHat.row(i) = X.row(i) * trans(sum(B.each_col()%w.col(i),0));
  }
  return yHat;
}



//' @title Cross-Validation Without Doing Cross-Validation on Bayesian GBLUP
//' 
//' @description Performs a Leave-One-Out Cross-Validation Without Doing Cross-Validation using Gibbs Sampling results.
//'
//' @param n is the number of observations. If length(n) > 1, the length is taken to be the number required.
//' @param Y is a phenotypic matrix with n rows and 1 column.
//' @param g is vector with the BLUP solution for the genetic values.
//' @param K is a relationship matrix with m rows and m columns.
//' @param Vg is the posterior mean of the genetic variance.
//' @param Ve is the posterior mean of the residual variance.
//' 
//' @examples
//' ## Not to run ##
//' 
//' ## cvBGBLUP(n=nSamp, Y=y, g=gHat, K=G, Vg=gVAR, Ve=eVAR)
//'
//' ## End(Not run)
//' 
//' @useDynLib ezCV
//'  
//' @export cvBGBLUP
// [[Rcpp::export]]
arma::mat cvBGBLUP(int n, arma::vec Y, arma::vec g, arma::mat K, double Vg, double Ve) {
  double lambda = Ve/Vg;
  int nIDs = K.n_rows;
  arma::mat invC = arma::inv_sympd(arma::eye(nIDs, nIDs) + arma::inv_sympd(K)*lambda);
  arma::mat sig = invC*Ve;
  arma::mat rng = arma::randn(n, nIDs);
  arma::mat gs = arma::repmat(g, 1, n).t() + rng * arma::chol(sig);
  arma::mat wis_G(n, nIDs);
  for(int i = 0; i < nIDs; i++) {
    arma::colvec wi = exp(pow((Y[i]-gs.col(i)),2)/(2.0*Ve))/accu(exp(pow((Y[i]-gs.col(i)),2)/(2.0*Ve)));
    wis_G.col(i) = wi;
  }
  arma::mat yHat = sum(wis_G%gs, 0).t();
  return yHat;
}



 
 