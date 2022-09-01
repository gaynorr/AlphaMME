// solveUVM and solveMVM are based on R/EMMREML functions
// solveMKM is based on the mmer function in R/sommer
#include "alphamme.h"
#include <iostream>
#include <string>
#ifdef ARMA_USE_LAPACK

#if !defined(ARMA_BLAS_CAPITALS)
#define arma_dsyevr dsyevr
#else
#define arma_dsyevr DSYEVR
#endif

extern "C"
void arma_fortran(arma_dsyevr)(char* JOBZ, char* RANGE, char* UPLO, long long int* N, double* A, long long int* LDA, double* VL,
                                          double* VU, long long int* IL, long long int* IU, double* ABSTOL, long long int* M, double* W, double* Z,
                                          long long int* LDZ, long long int* ISUPPZ, double* WORK, long long int* LWORK, long long int* IWORK,
                                          long long int* LIWORK, long long int* INFO);
#endif

const double PI = 3.14159265358979323846;

// // Note: Fortran compiler appends '_' to subroutine name
// // See http://www.netlib.org/lapack/explore-html/ for description of args
// extern "C" void dsyevr_(char* JOBZ, char* RANGE, char* UPLO, long long int* N, double* A, long long int* LDA, double* VL,
//                        double* VU, long long int* IL, long long int* IU, double* ABSTOL, long long int* M, double* W, double* Z,
//                        long long int* LDZ, long long int* ISUPPZ, double* WORK, long long int* LWORK, long long int* IWORK,
//                        long long int* LIWORK, long long int* INFO);

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
  long long int IL = 0;
  long long int IU = 0;
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
  long long int tmpIWORK = 0;
  long long int LIWORK = -1; // To be calculated
  long long int INFO = 0;
  // Calculate LWORK and LIWORK
  F77_CALL(dsyevr)(&JOBZ,&RANGE,&UPLO,&N,&*X.begin(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,&*eigval.begin(),
           &*eigvec.begin(),&LDZ,&*ISUPPZ.begin(),&tmpWORK,&LWORK,&tmpIWORK,&LIWORK,&INFO);
  LWORK = (long long int) tmpWORK;
  LIWORK = tmpIWORK;
  // Allocate WORK and IWORK
  arma::vec WORK(LWORK);
  arma::Col<long long int> IWORK(LIWORK);
  // Perform decomposition
  F77_CALL(dsyevr)(&JOBZ,&RANGE,&UPLO,&N,&*X.begin(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,&*eigval.begin(),
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

//' @title Read Matrix
//'
//' @description
//' Uses C++ to quickly read a matrix from a text
//' file. Requires knowledge of the number of rows
//' and columns in the file.
//'
//' @param fileName path to the file being read
//' @param rows number of rows to read in
//' @param cols number of columns to read in
//' @param sep a single character delimiter seperating data entries
//' @param skipRows number of rows to skip
//' @param skipCols number of columns to skip
//'
//' @return a numeric matrix
//'
//' @export
// [[Rcpp::export]]
arma::mat readMat(std::string fileName, int rows, int cols,
                  char sep=' ', int skipRows=0, int skipCols=0){
  arma::mat output(rows,cols);
  std::ifstream file(fileName.c_str());
  std::string line;
  //Skip rows
  for(arma::uword i=0; i<skipRows; ++i){
    std::getline(file,line);
  }
  //Read rows
  for(arma::uword i=0; i<rows; ++i){
    std::getline(file,line);
    std::stringstream lineStream(line);
    std::string cell;
    //Skip columns
    for(arma::uword j=0; j<skipCols; ++j){
      std::getline(lineStream,cell,sep);
    }
    //Read columns
    for(arma::uword j=0; j<cols; ++j){
      std::getline(lineStream,cell,sep);
      output(i,j) = std::atof(cell.c_str());
    }
  }
  file.close();
  return output;
}

//' @title Solve Univariate Model
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+Zu+e}
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param Z a matrix with n rows and m columns
//' @param K a matrix with m rows and m columns
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveUVM(const arma::mat& y, const arma::mat& X,
                    const arma::mat& Z, const arma::mat& K){
  int n = y.n_rows;
  int q = X.n_cols;
  double df = double(n)-double(q);
  double offset = log(double(n));

  // Construct system of equations for eigendecomposition
  arma::mat S = -(X*inv_sympd(X.t()*X)*X.t());
  S.diag() += 1;
  arma::mat ZK = Z*K;
  arma::mat H = ZK*Z.t(); // Used later
  H.diag() += offset;
  S = S*H*S;

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
  H.diag() += (delta-offset);
  H = inv_sympd(H);
  arma::mat XH = X.t()*H;
  arma::mat beta = solve(XH*X,XH*y);
  arma::mat u = ZK.t()*(H*(y-X*beta));
  double Vu = sum(eta%eta/(eigval+delta))/df;
  double Ve = delta*Vu;
  double ll = -0.5*(double(optRes["objective"])+df+df*log(2*PI/df));
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=ll);
}

//' @title Solve animal model
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+u+e}
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param K the numeric relationship matrix
//' with n rows and n columns
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveAniModel(const arma::mat& y,
                         const arma::mat& X,
                         const arma::mat& K){
  int n = y.n_rows;
  int q = X.n_cols;
  double df = double(n)-double(q);
  double offset = log(double(n));

  // Construct system of equations for eigendecomposition
  arma::mat S = -(X*inv_sympd(X.t()*X)*X.t());
  S.diag() += 1;
  arma::mat H = K; //Used later
  H.diag() += offset;
  S = S*H*S;

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
  H.diag() += (delta-offset);
  H = inv_sympd(H);
  arma::mat XH = X.t()*H;
  arma::mat beta = solve(XH*X,XH*y);
  arma::mat u = K.t()*(H*(y-X*beta));
  double Vu = sum(eta%eta/(eigval+delta))/df;
  double Ve = delta*Vu;
  double ll = -0.5*(double(optRes["objective"])+df+df*log(2*PI/df));
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=ll);
}

//' @title Solve RR-BLUP
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+Mu+e}
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param M a matrix with n rows and m columns
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUP(const arma::mat& y, const arma::mat& X,
                       const arma::mat& M){
  int n = y.n_rows;
  int q = X.n_cols;
  double df = double(n)-double(q);
  double offset = log(double(n));

  // Construct system of equations for eigendecomposition
  arma::mat S = -(X*inv_sympd(X.t()*X)*X.t());
  S.diag() += 1;
  arma::mat H = M*M.t(); // Used later
  H.diag() += offset;
  S = S*H*S;

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
  H.diag() += (delta-offset);
  H = inv_sympd(H);
  arma::mat XH = X.t()*H;
  arma::mat beta = solve(XH*X,XH*y);
  arma::mat u = M.t()*(H*(y-X*beta));
  double Vu = sum(eta%eta/(eigval+delta))/df;
  double Ve = delta*Vu;
  double ll = -0.5*(double(optRes["objective"])+df+df*log(2*PI/df));
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=ll);
}

//' @title Solve Multivariate Model
//'
//' @description
//' Solves a multivariate mixed model of form \eqn{Y=X\beta+Zu+e}
//'
//' @param Y a matrix with n rows and q columns
//' @param X a matrix with n rows and x columns
//' @param Z a matrix with n rows and m columns
//' @param K a matrix with m rows and m columns
//' @param tol tolerance for convergence
//' @param maxIter maximum number of iteration
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveMVM(const arma::mat& Y, const arma::mat& X,
                    const arma::mat& Z, const arma::mat& K,
                    double tol=1e-6, int maxIter=1000){
  int n = Y.n_rows;
  int m = Y.n_cols;
  arma::mat ZK = Z*K;
  arma::mat ZKZ = ZK*Z.t();
  arma::vec eigval(n);
  arma::mat eigvec(n,n);
  eigen2(eigval, eigvec, ZKZ);
  arma::mat Yt = Y.t()*eigvec;
  arma::mat Xt = X.t()*eigvec;
  arma::mat Vu = cov(Y)/2;
  arma::mat Ve = Vu;
  arma::mat W = Xt.t()*inv_sympd(Xt*Xt.t());
  arma::mat B = Yt*W; //BLUEs
  arma::mat Gt(m,n), sigma(m,m), BNew,
  VeNew(m,m), VuNew(m,m);
  double denom, numer;
  bool converging=true;
  int iter=0;
  while(converging){
    ++iter;
    VeNew.fill(0.0);
    VuNew.fill(0.0);
    for(arma::uword i=0; i<n; ++i){
      Gt.col(i) = eigval(i)*Vu*inv_sympd(eigval(i)*Vu+
        Ve+tol*arma::eye(m,m))*(Yt.col(i)-B*Xt.col(i));
    }
    BNew = (Yt - Gt)*W;
    for(arma::uword i=0; i<n; ++i){
      sigma = eigval(i)*Vu-(eigval(i)*Vu)*inv_sympd(eigval(i)*Vu+
        Ve+tol*arma::eye(m,m))*(eigval(i)*Vu);
      VuNew += 1.0/(double(n)*eigval(i))*(Gt.col(i)*Gt.col(i).t()+sigma);
      VeNew += 1.0/double(n)*((Yt.col(i)-BNew*Xt.col(i)-Gt.col(i))*
        (Yt.col(i)-BNew*Xt.col(i)-Gt.col(i)).t()+sigma);
    }
    denom = fabs(sum(Ve.diag()));
    if(denom>0.0){
      numer = fabs(sum(VeNew.diag()-Ve.diag()));
      if((numer/denom)<tol) converging=false;
    }
    Ve = VeNew;
    Vu = VuNew;
    B = BNew;
    if(iter>=maxIter){
      Rf_warning("Reached maxIter without converging");
      break;
    }
  }
  arma::mat HI = inv_sympd(kron(ZKZ, Vu)+
    kron(arma::eye(n,n), Ve)+
    tol*arma::eye(n*m,n*m));
  arma::mat E = Y.t() - B*X.t();
  arma::mat U = kron(K, Vu)*kron(Z.t(),
                     arma::eye(m,m))*(HI*vectorise(E)); //BLUPs
  U.reshape(m,U.n_elem/m);
  //Log Likelihood calculation
  arma::mat ll = -0.5*arma::vectorise(E).t()*HI*vectorise(E);
  ll -= double(n*m)/2.0*log(2*PI);
  double value;
  double sign;
  log_det(value, sign, kron(ZKZ, Vu)+kron(arma::eye(n,n), Ve));
  ll -= 0.5*value*sign;
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=B.t(),
                            Rcpp::Named("u")=U.t(),
                            Rcpp::Named("LL")=arma::as_scalar(ll),
                            Rcpp::Named("iter")=iter);
}

//' @title Solve Multivariate RR-BLUP
//'
//' @description
//' Solves a multivariate mixed model of form \eqn{Y=X\beta+Mu+e}
//'
//' @param Y a matrix with n rows and q columns
//' @param X a matrix with n rows and x columns
//' @param M a matrix with n rows and m columns
//' @param tol tolerance for convergence
//' @param maxIter maximum number of iteration
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUPMV(const arma::mat& Y, const arma::mat& X,
                         const arma::mat& M, double tol=1e-6,
                         int maxIter=1000){
  int n = Y.n_rows;
  int m = Y.n_cols;
  arma::vec eigval(n);
  arma::mat eigvec(n,n);
  eigen2(eigval, eigvec, M*M.t());
  arma::mat Yt = Y.t()*eigvec;
  arma::mat Xt = X.t()*eigvec;
  arma::mat Vu = cov(Y)/2;
  arma::mat Ve = Vu;
  arma::mat W = Xt.t()*inv_sympd(Xt*Xt.t());
  arma::mat B = Yt*W; //BLUEs
  arma::mat Gt(m,n), sigma(m,m), BNew,
  VeNew(m,m), VuNew(m,m);
  double denom, numer;
  bool converging=true;
  int iter=0;
  while(converging){
    ++iter;
    VeNew.fill(0.0);
    VuNew.fill(0.0);
    for(arma::uword i=0; i<n; ++i){
      Gt.col(i) = eigval(i)*Vu*inv_sympd(eigval(i)*Vu+
        Ve+tol*arma::eye(m,m))*(Yt.col(i)-B*Xt.col(i));
    }
    BNew = (Yt - Gt)*W;
    for(arma::uword i=0; i<n; ++i){
      sigma = eigval(i)*Vu-(eigval(i)*Vu)*inv_sympd(eigval(i)*Vu+
        Ve+tol*arma::eye(m,m))*(eigval(i)*Vu);
      VuNew += 1.0/(double(n)*eigval(i))*(Gt.col(i)*Gt.col(i).t()+sigma);
      VeNew += 1.0/double(n)*((Yt.col(i)-BNew*Xt.col(i)-Gt.col(i))*
        (Yt.col(i)-BNew*Xt.col(i)-Gt.col(i)).t()+sigma);
    }
    denom = fabs(sum(Ve.diag()));
    if(denom>0.0){
      numer = fabs(sum(VeNew.diag()-Ve.diag()));
      if((numer/denom)<tol) converging=false;
    }
    Ve = VeNew;
    Vu = VuNew;
    B = BNew;
    if(iter>=maxIter){
      Rcpp::Rcerr<<"Warning: did not converge, reached maxIter\n";
      break;
    }
  }
  arma::mat HI = inv_sympd(kron(M*M.t(), Vu)+
    kron(arma::eye(n,n), Ve)+
    tol*arma::eye(n*m,n*m));
  arma::mat E = Y.t() - B*X.t();
  arma::mat U = kron(arma::eye(M.n_cols,M.n_cols), Vu)*kron(M.t(),
                     arma::eye(m,m))*(HI*vectorise(E)); //BLUPs
  U.reshape(m,U.n_elem/m);
  //Log Likelihood calculation
  arma::mat ll = -0.5*arma::vectorise(E).t()*HI*vectorise(E);
  ll -= double(n*m)/2.0*log(2*PI);
  double value;
  double sign;
  log_det(value, sign, kron(M*M.t(), Vu)+kron(arma::eye(n,n), Ve));
  ll -= 0.5*value*sign;
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=B.t(),
                            Rcpp::Named("u")=U.t(),
                            Rcpp::Named("LL")=arma::as_scalar(ll),
                            Rcpp::Named("iter")=iter);
}

//' @title Solve Multikernel Model
//'
//' @description
//' Solves a univariate mixed model with multiple random effects.
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param Zlist a list of Z matrices
//' @param Klist a list of K matrices
//' @param maxIter maximum number of iteration
//' @param tol tolerance for convergence
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveMKM(arma::mat& y, arma::mat& X,
                     arma::field<arma::mat>& Zlist,
                     arma::field<arma::mat>& Klist,
                     int maxIter=40, double tol=1e-4){
  int k = Klist.n_elem;
  int n = y.n_rows;
  int q = X.n_cols;
  double df = double(n)-double(q);
  arma::field<arma::mat> V(k);
  for(arma::uword i=0; i<k; ++i){
    V(i) = Zlist(i)*Klist(i)*Zlist(i).t();
  }
  arma::mat A(k+1,k+1), W0(n,n), W(n,n), WX(n,q), WQX(n,n);
  arma::vec qvec(k+1), sigma(k+1);
  double rss, ldet, llik, llik0, deltaLlik, taper,
  value, sign;
  bool invPass;
  arma::field<arma::mat> T(k);
  sigma.fill(var(y.col(0)));
  int iter = 0;
  while(true){
    ++iter;
    W0 = V(0)*sigma(0);
    W0.diag() += sigma(k);
    for(arma::uword i=1; i<k; ++i){
      W0 += V(i)*sigma(i);
    }
    invPass = inv_sympd(W,W0);
    if(!invPass){
      W = pinv(W0);
    }
    WX = W*X;
    WQX = W - WX*solve(X.t()*WX, WX.t());
    rss = as_scalar(y.t()*WQX*y);
    sigma = sigma*(rss/df);
    WQX = WQX*(df/rss);
    log_det(value, sign, WQX);
    ldet = value*sign;
    llik = ldet/2 - df/2;
    if(iter == 1) llik0 = llik;
    deltaLlik = llik - llik0;
    llik0 = llik;
    for(arma::uword i=0; i<k; ++i){
      T(i) = WQX*V(i);
    }
    for(arma::uword i=0; i<k; ++i){
      qvec(i) = as_scalar(y.t()*T(i)*WQX*y - sum(T(i).diag()));
      for(arma::uword j=0; j<k; ++j){
        A(i,j) = accu(T(i)%T(j).t());
      }
      A(i,k) = accu(T(i)%WQX.t());
    }
    for(arma::uword j=0; j<k; ++j){
      A(k,j) = accu(WQX%T(j).t());
    }
    A(k,k) = accu(WQX%WQX.t());
    qvec(k) = as_scalar(y.t()*WQX*WQX*y - sum(WQX.diag()));
    A = pinv(A);
    qvec = A*qvec;
    if(iter == 1){
      taper = 0.5;
    }else if(iter == 2){
      taper = 0.7;
    }else{
      taper = 0.9;
    }
    sigma += taper*qvec;
    while(sigma.min() < -(1e-6)){
      sigma(sigma.index_min()) = -(1e-6);
    }
    if(iter > 1 & fabs(deltaLlik) < tol*10){
      break;
    }
    if(max(abs(qvec)) < tol){
      break;
    }
    if(iter >= maxIter){
      Rf_warning("Reached maxIter without converging");
      break;
    }
  }
  while(sigma.min() < 0.0){
    sigma(sigma.index_min()) = 0.0;
  }
  arma::mat beta(q,1), ee(n,1);
  arma::field<arma::mat> u(k);
  beta = solve(X.t()*W*X,X.t()*W*y);
  ee = y - X*beta;
  for(arma::uword i=0; i<k; ++i){
    u(i) = (Klist(i)*sigma(i))*Zlist(i).t()*W*ee;
  }
  arma::vec Vu(k), Ve(1);
  Vu = sigma(arma::span(0,k-1));
  Ve = sigma(k);
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=llik,
                            Rcpp::Named("iter")=iter);
}

//' @title Solve Multikernel RR-BLUP
//'
//' @description
//' Solves a univariate mixed model with multiple random effects.
//'
//' @param y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param Mlist a list of M matrices
//' @param maxIter maximum number of iteration
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUPMK(arma::mat& y, arma::mat& X,
                         arma::field<arma::mat>& Mlist,
                         int maxIter=40){
  double tol = 1e-4;
  int k = Mlist.n_elem;
  int n = y.n_rows;
  int q = X.n_cols;
  double df = double(n)-double(q);
  arma::field<arma::mat> V(k);
  for(arma::uword i=0; i<k; ++i){
    V(i) = Mlist(i)*Mlist(i).t();
  }
  arma::mat A(k+1,k+1), W0(n,n), W(n,n), WX(n,q), WQX(n,n);
  arma::vec qvec(k+1), sigma(k+1);
  double rss, ldet, llik, llik0, deltaLlik, taper,
  value, sign;
  bool invPass;
  arma::field<arma::mat> T(k);
  sigma.fill(var(y.col(0)));
  int iter = 0;
  while(true){
    ++iter;
    W0 = V(0)*sigma(0);
    W0.diag() += sigma(k);
    for(arma::uword i=1; i<k; ++i){
      W0 += V(i)*sigma(i);
    }
    invPass = inv_sympd(W,W0);
    if(!invPass){
      W = pinv(W0);
    }
    WX = W*X;
    WQX = W - WX*solve(X.t()*WX, WX.t());
    rss = as_scalar(y.t()*WQX*y);
    sigma = sigma*(rss/df);
    WQX = WQX*(df/rss);
    log_det(value, sign, WQX);
    ldet = value*sign;
    llik = ldet/2 - df/2;
    if(iter == 1) llik0 = llik;
    deltaLlik = llik - llik0;
    llik0 = llik;
    for(arma::uword i=0; i<k; ++i){
      T(i) = WQX*V(i);
    }
    for(arma::uword i=0; i<k; ++i){
      qvec(i) = as_scalar(y.t()*T(i)*WQX*y - sum(T(i).diag()));
      for(arma::uword j=0; j<k; ++j){
        A(i,j) = accu(T(i)%T(j).t());
      }
      A(i,k) = accu(T(i)%WQX.t());
    }
    for(arma::uword j=0; j<k; ++j){
      A(k,j) = accu(WQX%T(j).t());
    }
    A(k,k) = accu(WQX%WQX.t());
    qvec(k) = as_scalar(y.t()*WQX*WQX*y - sum(WQX.diag()));
    A = pinv(A);
    qvec = A*qvec;
    if(iter == 1){
      taper = 0.5;
    }else if(iter == 2){
      taper = 0.7;
    }else{
      taper = 0.9;
    }
    sigma += taper*qvec;
    while(sigma.min() < -(1e-6)){
      sigma(sigma.index_min()) = -(1e-6);
    }
    if(iter > 1 & fabs(deltaLlik) < tol*10){
      break;
    }
    if(max(abs(qvec)) < tol){
      break;
    }
    if(iter >= maxIter){
      Rf_warning("Reached maxIter without converging");
      break;
    }
  }
  while(sigma.min() < 0.0){
    sigma(sigma.index_min()) = 0.0;
  }
  arma::mat beta(q,1), ee(n,1);
  arma::field<arma::mat> u(k);
  beta = solve(X.t()*W*X,X.t()*W*y);
  ee = y - X*beta;
  for(arma::uword i=0; i<k; ++i){
    u(i) = sigma(i)*Mlist(i).t()*W*ee;
  }
  arma::vec Vu(k), Ve(1);
  Vu = sigma(arma::span(0,k-1));
  Ve = sigma(k);
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=beta,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("LL")=llik,
                            Rcpp::Named("iter")=iter);
}

//' @title Solve RR-BLUP with EM
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+Mu+e} using
//' the Expectation-Maximization algorithm.
//'
//' @param Y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param M a matrix with n rows and m columns
//' @param Vu initial guess for variance of marker effects
//' @param Ve initial guess for error variance
//' @param tol tolerance for declaring convergence
//' @param maxIter maximum iteration for attempting convergence
//' @param useEM should EM algorithm be used. If false, no estimation of
//' variance components is performed. The initial values are treated as true.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUP_EM(const arma::mat& Y, const arma::mat& X,
                          const arma::mat& M, double Vu, double Ve,
                          double tol=1e-6, int maxIter=100,
                          bool useEM=true){
  double lambda = Ve/Vu;
  double delta=0;
  int iter=0;
  arma::uword n=Y.n_rows,m=M.n_cols,q=X.n_cols;
  arma::mat RHS(q+m,q+m),LHS(q+m,1),Rvec(q+m,1);
  RHS(arma::span(0,q-1),arma::span(0,q-1)) = X.t()*X;
  RHS(arma::span(0,q-1),arma::span(q,q+m-1)) = X.t()*M;
  RHS(arma::span(q,q+m-1),arma::span(0,q-1)) = M.t()*X;
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)) = M.t()*M;
  RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += lambda;
  Rvec(arma::span(0,q-1),0) = X.t()*Y;
  Rvec(arma::span(q,q+m-1),0) = M.t()*Y;
  arma::mat RHSinv = pinv(RHS);
  LHS = RHSinv*Rvec;
  if(useEM){
    Ve = as_scalar(Y.t()*Y-LHS.t()*Rvec)/(n-q);
    Vu = as_scalar(
      LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
        Ve*sum(RHSinv(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag())
    )/m;
    delta = Ve/Vu-lambda;
    while(fabs(delta)>tol){
      RHS(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag() += delta;
      lambda += delta;
      RHSinv = pinv(RHS);
      LHS = RHSinv*Rvec;
      iter++;
      if(iter>=maxIter){
        Rcpp::Rcerr<<"Warning: did not converge, reached maxIter\n";
        break;
      }
      Ve = as_scalar(Y.t()*Y-LHS.t()*Rvec)/(n-q);
      Vu = as_scalar(
        LHS(arma::span(q,q+m-1),0).t()*LHS(arma::span(q,q+m-1),0)+
          Ve*sum(RHSinv(arma::span(q,q+m-1),arma::span(q,q+m-1)).diag())
      )/m;
      delta = Ve/Vu-lambda;
    }
  }
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=LHS.rows(arma::span(0,q-1)),
                            Rcpp::Named("u")=LHS.rows(arma::span(q,q+m-1)),
                            Rcpp::Named("iter")=iter);
}

//' @title Solve RR-BLUP with EM and 2 random effects
//'
//' @description
//' Solves a univariate mixed model of form \eqn{y=X\beta+M_1u_1+M_2u_2+e} using
//' the Expectation-Maximization algorithm.
//'
//' @param Y a matrix with n rows and 1 column
//' @param X a matrix with n rows and x columns
//' @param M1 a matrix with n rows and m1 columns
//' @param M2 a matrix with n rows and m2 columns
//' @param Vu1 initial guess for variance of the first marker effects
//' @param Vu2 initial guess for variance of the second marker effects
//' @param Ve initial guess for error variance
//' @param tol tolerance for declaring convergence
//' @param maxIter maximum iteration for attempting convergence
//' @param useEM should EM algorithm be used. If false, no estimation of
//' variance components is performed. The initial values are treated as true.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List solveRRBLUP_EM2(const arma::mat& Y, const arma::mat& X,
                           const arma::mat& M1, const arma::mat& M2,
                           double Vu1, double Vu2, double Ve,
                           double tol=1e-6, int maxIter=100,
                           bool useEM=true){
  double lambda1 = Ve/Vu1;
  double lambda2 = Ve/Vu2;
  double delta1=0,delta2=0,VeN=0,Vu1N=0,Vu2N=0;
  int iter=0;
  arma::uword n=Y.n_rows,m1=M1.n_cols,m2=M2.n_cols,q=X.n_cols;
  arma::mat RHS(q+m1+m2,q+m1+m2),LHS(q+m1+m2,1),Rvec(q+m1+m2,1);
  // Top row
  RHS(arma::span(0,q-1),arma::span(0,q-1)) = X.t()*X;
  RHS(arma::span(0,q-1),arma::span(q,q+m1-1)) = X.t()*M1;
  RHS(arma::span(0,q-1),arma::span(q+m1,q+m1+m2-1)) = X.t()*M2;
  // Second row
  RHS(arma::span(q,q+m1-1),arma::span(0,q-1)) = M1.t()*X;
  RHS(arma::span(q,q+m1-1),arma::span(q,q+m1-1)) = M1.t()*M1;
  RHS(arma::span(q,q+m1-1),arma::span(q+m1,q+m1+m2-1)) = M1.t()*M2;
  // Third row
  RHS(arma::span(q+m1,q+m1+m2-1),arma::span(0,q-1)) = M2.t()*X;
  RHS(arma::span(q+m1,q+m1+m2-1),arma::span(q,q+m1-1)) = M2.t()*M1;
  RHS(arma::span(q+m1,q+m1+m2-1),arma::span(q+m1,q+m1+m2-1)) = M2.t()*M2;
  // Add to diagonal
  RHS(arma::span(q,q+m1-1),arma::span(q,q+m1-1)).diag() += lambda1;
  RHS(arma::span(q+m1,q+m1+m2-1),arma::span(q+m1,q+m1+m2-1)).diag() += lambda2;
  Rvec(arma::span(0,q-1),0) = X.t()*Y;
  Rvec(arma::span(q,q+m1-1),0) = M1.t()*Y;
  Rvec(arma::span(q+m1,q+m1+m2-1),0) = M2.t()*Y;
  arma::mat RHSinv = inv(RHS);
  LHS = RHSinv*Rvec;
  if(useEM){
    VeN = as_scalar(Y.t()*Y-LHS.t()*Rvec)/(n-q);
    Vu1N = as_scalar(
      LHS(arma::span(q,q+m1-1),0).t()*LHS(arma::span(q,q+m1-1),0)+
        Ve*sum(RHSinv(arma::span(q,q+m1-1),arma::span(q,q+m1-1)).diag())
    )/m1;
    Vu2N = as_scalar(
      LHS(arma::span(q+m1,q+m1+m2-1),0).t()*LHS(arma::span(q+m1,q+m1+m2-1),0)+
        Ve*sum(RHSinv(arma::span(q+m1,q+m1+m2-1),arma::span(q+m1,q+m1+m2-1)).diag())
    )/m2;
    delta1 = VeN/Vu1N-lambda1;
    delta2 = VeN/Vu2N-lambda2;
    while((fabs(delta1)>tol) || (fabs(delta2)>tol)){
      Ve = VeN;
      Vu1 = Vu1N;
      Vu2 = Vu2N;
      RHS(arma::span(q,q+m1-1),arma::span(q,q+m1-1)).diag() += delta1;
      RHS(arma::span(q+m1,q+m1+m2-1),arma::span(q+m1,q+m1+m2-1)).diag() += delta2;
      lambda1 += delta1;
      lambda2 += delta2;
      RHSinv = inv(RHS);
      LHS = RHSinv*Rvec;
      iter++;
      if(iter>=maxIter){
        Rcpp::Rcerr<<"Warning: did not converge, reached maxIter\n";
        break;
      }
      VeN = as_scalar(Y.t()*Y-LHS.t()*Rvec)/(n-q);
      Vu1N = as_scalar(
        LHS(arma::span(q,q+m1-1),0).t()*LHS(arma::span(q,q+m1-1),0)+
          Ve*sum(RHSinv(arma::span(q,q+m1-1),arma::span(q,q+m1-1)).diag())
      )/m1;
      Vu2N = as_scalar(
        LHS(arma::span(q+m1,q+m1+m2-1),0).t()*LHS(arma::span(q+m1,q+m1+m2-1),0)+
          Ve*sum(RHSinv(arma::span(q+m1,q+m1+m2-1),arma::span(q+m1,q+m1+m2-1)).diag())
      )/m2;
      delta1 = VeN/Vu1N-lambda1;
      delta2 = VeN/Vu2N-lambda2;
    }
  }
  arma::vec Vu(2);
  Vu(0) = Vu1;
  Vu(1) = Vu2;
  arma::field<arma::mat> u(2);
  u(0).set_size(m1,1);
  u(1).set_size(m2,1);
  u(0) = LHS.rows(arma::span(q,q+m1-1));
  u(1) = LHS.rows(arma::span(q+m1,q+m1+m2-1));
  return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("beta")=LHS.rows(arma::span(0,q-1)),
                            Rcpp::Named("u")=u,
                            Rcpp::Named("iter")=iter);
}


//' @title Calculate G Matrix
//'
//' @description
//' Calculates the genomic relationship matrix.
//'
//' @param X a matrix of marker genotypes scored as 0,1,2
//'
//' @return a matrix of the realized genomic relationships
//'
//' @export
// [[Rcpp::export]]
arma::mat calcG(arma::mat X){
  arma::rowvec p = mean(X,0)/2.0;
  X.each_row() -= 2*p;
  arma::mat G = X*X.t();
  G = G/(2.0*sum(p%(1-p)));
  return G;
}

//' @title Calculate Dominance Matrix
//'
//' @description
//' Calculates the dominance relationship matrix.
//'
//' @param X a matrix of marker genotypes scored as 0,1,2
//'
//' @references
//' \cite{Nishio, M, and M. Satoh. 2014. Including Dominance Effects in the Genomic BLUP Method for Genomic Evaluation. PLOS ONE 9(1): e85792.}
//'
//' @return a matrix of the realized dominance relationships
//'
//' @export
// [[Rcpp::export]]
arma::mat calcD(arma::mat X){
  arma::rowvec p = mean(X,0)/2.0;
  arma::rowvec q = 1-p;
  double zero,one,two;
  int x;
  for(arma::uword i=0; i<X.n_cols; ++i){
    zero = -2.0*(p(i)*p(i));
    one = 2.0*p(i)*q(i);
    two = -2.0*(q(i)*q(i));
    for(arma::uword j=0; j<X.n_rows; ++j){
      x = int(X(j,i));
      if(x==2){
        X(j,i) = two;
      }else if(x==1){
        X(j,i) = one;
      }else{
        X(j,i) = zero;
      }
    }
  }
  arma::mat D = X*X.t();
  D = D/(4.0*sum(p%p%q%q));
  return D;
}


//' @title Calculate IBS G Matrix
//'
//' @description
//' Calculates an identity-by-state genomic relationship matrix
//' based on simple matching.
//'
//' @param X a matrix of marker genotypes scored as 0,1,2
//'
//' @return a matrix of genomic relationships
//'
//' @export
// [[Rcpp::export]]
arma::mat calcGIbs(arma::mat X){
  X -= 1.0;
  arma::mat G = (X*X.t())/X.n_cols + 1.0;
  return G;
}

// Calculates a distance matrix from a marker matrix
// Uses a binomial theorem trick
// Inspired by code from:
// http://blog.felixriedel.com/2013/05/pairwise-distances-in-r/
// First described here:
// http://blog.smola.org/post/969195661/in-praise-of-the-second-binomial-formula
//' @title Calculate Euclidean distance
//'
//' @description
//' Calculates a Euclidean distance matrix using a binomial
//' theorem trick. Results in much faster computation than the
//' \code{dist} function in package \code{stats}.
//'
//' @param X a numeric matrix
//'
//' @return a matrix of columnwise distances
//'
//' @export
// [[Rcpp::export]]
arma::mat fastDist(const arma::mat& X){
  arma::colvec Xn = sum(square(X),1);
  arma::mat D = -2*(X*X.t());
  D.each_col() += Xn;
  D.each_row() += Xn.t();
  D = sqrt(D);
  D.diag().zeros(); //Removes NaN values
  if(D.has_nan()){
    D.elem(find_nonfinite(D)).fill(0.0);
  }
  return D;
}

//' @title Calculate Paired Euclidean distance
//'
//' @description
//' Calculates a Euclidean distance between two matrices using
//' a binomial theorem trick.
//'
//' @param X a numeric matrix
//' @param Y a numeric matrix
//'
//' @return a matrix of columnwise distances between matrices
//' X and Y
//'
//' @export
// [[Rcpp::export]]
arma::mat fastPairDist(const arma::mat& X, const arma::mat& Y){
  arma::colvec Xn = sum(square(X),1);
  arma::colvec Yn = sum(square(Y),1);
  arma::mat D = -2*(X*Y.t());
  D.each_col() += Xn;
  D.each_row() += Yn.t();
  D = sqrt(D);
  if(D.has_nan()){
    D.elem(find_nonfinite(D)).fill(0.0);
  }
  return D;
}

//' @title Calculate Gaussian Kernel
//'
//' @description
//' Calculates a Gaussian kernel using a Euclidean distance
//' matrix.
//'
//' @param D a matrix of Euclidean distances,
//' see \code{\link{fastDist}}
//' @param theta the tuning parameter
//'
//' @return a numeric matrix
//'
//' @export
// [[Rcpp::export]]
arma::mat gaussKernel(arma::mat& D, double theta){
  return exp(-1.0*square(D/theta));
}

// // Objective function for Gaussian kernel method
// Rcpp::List objRKHS(double theta, Rcpp::List args){
//   Rcpp::List output;
//   arma::mat D = args["D"];
//   output = solveAniModel(args["y"],args["X"],
//                          gaussKernel(D,theta));
//   return Rcpp::List::create(Rcpp::Named("objective")=output["LL"],
//                             Rcpp::Named("output")=output);
// }
//
// //' @title Solve RKHS
// //'
// //' @description
// //' Solves a Reproducing Kernel Hilbert Space regression
// //' using a Gaussian Kernel.
// //'
// //' @param y a matrix with n rows and 1 column
// //' @param X a matrix with n rows and x columns
// //' @param M a matrix with n rows and m columns
// //'
// //' @export
// // [[Rcpp::export]]
// Rcpp::List solveRKHS(const arma::mat& y, const arma::mat& X,
//                      const arma::mat& M){
//
//
//   return Rcpp::List::create(Rcpp::Named("Vu")=Vu,
//                             Rcpp::Named("Ve")=Ve,
//                             Rcpp::Named("beta")=beta,
//                             Rcpp::Named("u")=u,
//                             Rcpp::Named("LL")=ll);
// }
