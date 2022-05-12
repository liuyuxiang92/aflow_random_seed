// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_XFIT_CPP_
#define _AUROSTD_XFIT_CPP_

#ifndef _AUROSTD_XFIT_H_
#include "aurostd_xfit.h"
#endif

#define DEBUG_XFIT false // toggles debug output for xfit

//********************************************************************************
//              Functions to work with polynomials
namespace aurostd{
  //Evaluates the value of the polynomial with coefficients p at the value x.
  //
  //The polynomial is represented in the following form:
  //p(x) = p0+x*(p1+x*(p2+x*(p3+...x*(p(n-1)+x*(pn)))))
  template<class utype> utype evalPolynomial(const utype x, const xvector<utype>& p) {
    utype res = p[p.urows];
    for (int i=p.urows-1; i>=p.lrows; i--) res = res*x + p[i];
    return res;
  }

  //SD20220422
  template<class utype> xvector<utype> evalPolynomial_xv(const xvector<utype>& x, const xvector<utype>& p) {
    xvector<utype> res(x.rows);
    for (int i = 1; i <= res.rows; i++) {
      res(i) = evalPolynomial(x(i), p);
    }
    return res;
  }

  //SD20220505
  template<class utype> xmatrix<utype> evalPolynomial_xm(const xmatrix<utype>& x, const xvector<utype>& p) {
    xmatrix<utype> res(x.rows, x.cols);
    for (int i = 1; i <= res.rows; i++) {
      for (int j = 1; j <= res.cols; j++) {
        res(i, j) = evalPolynomial(x(i, j), p);
      }
    }
    return res;
  }

  //Evaluates the value and the derivatives of the polynomial with coefficients p at
  //the value x and outputs the derivatives in the dp array.
  //
  //The highest derivative is determined by the size of the dp array and the result
  //is stored in ascending order starting from the zero's derivative (function itself).
  template<class utype> void evalPolynomialDeriv(const utype x, const xvector<utype>& p, xvector<utype>& dp) {
    for (int i=dp.lrows; i<=dp.urows; i++) dp[i] = 0.0;
    dp[dp.lrows] = p[p.urows];
    for (int i=p.urows-1; i>=p.lrows; i--){
      for (int j=dp.urows; j>=dp.lrows+1; j--) dp[j] = dp[j]*x + dp[j-1];
      dp[dp.lrows] = dp[dp.lrows]*x + p[i];
    }

    int factor = 1;
    for (int i=dp.lrows+1; i<dp.urows; i++){
      factor *= i;
      dp[i+1] *= factor;
    }
  }

  //Evaluates the value and the derivatives of the polynomial with coefficients p at
  //the value x and returns the result as an array (dp).
  //
  //The result is stored in ascending order starting from the zero's derivative
  //(function itself).
  //@param n is the highest order of the derivative to be calculated.
  template<class utype> xvector<utype> evalPolynomialDeriv(const utype x, const xvector<utype>& p, const uint n) {
    xvector<utype> dp(n+1);
    evalPolynomialDeriv(x, p, dp);
    return dp;
  }

  //SD20220422
  //Returns the coefficients of the nth derivative of the polynomial
  template<class utype> xvector<utype> evalPolynomialCoeff(const xvector<utype>& p, const uint n) {
    if ((int)n >= p.rows) {return 0.0 * ones_xv<utype>(1);}
    xvector<utype> dp(p.rows - n);
    for (int i = p.rows; i > (int)n; i--) {
      dp(i - n) = p(i) * aurostd::factorial<utype>(i - 1) / aurostd::factorial<utype>(i - n - 1);
    }
    return dp;
  }

  //Constructs the Vandermonde matrix with columns up to a given number of columns n,
  //which corresponds to the polynomial of degree n-1.
  //Vandermonde matrix:
  // 1 x1^1 x1^2 .. x1^(n-1)
  // 1 x2^1 x2^2 .. x2^(n-1)
  // ...
  // 1 xm^1 xm^2 .. xm^(n-1)
  //where m is the size of x array
  template<class utype> xmatrix<utype> Vandermonde_matrix(const xvector<utype>& x, const int n) {
    xmatrix<utype> VM(x.rows, n);

    for (int i=1; i<=x.rows; i++){
      VM[i][1] = 1.0;
      for (int j=2; j<=n; j++) VM[i][j] = std::pow(x[i], j-1);
    }

    return VM;
  }

  //Calculates the extremum of the polynomial bounded by the region [xmin,xmax] by searching
  //for the value x where the derivative of polynomial equals to zero using the
  //bisection method.
  //
  //It is assumed that only one extremum exists between xmin and xmax.
  //The specific use of this function is to search for the equilibrium volume of a given
  //polynomial E(V) dependency.
  //
  //The main idea: if a continuous function f(x) has the root f(x0)=0 in the interval [a,b], 
  //then its sign in the interval [a,x0) is opposite to its sign in the interval (x0,b].
  //The root searching algorithm is iterative: one bisects the interval [a,b] into two 
  //subintervals and for the next iteration picks a subinterval where f(x) has 
  //opposite sign on its ends.
  //This process continues until the length of the interval is less than some negligibly
  //small number and one picks the middle of the interval as x0.
  //In this situation, f(x0)~0 and this condition is used as a criterion to stop the
  //iteration procedure.
  //SD20220427 - Throws error if zero is not found or bisection method is not converged
  template<class utype> utype polynomialFindExtremum(const xvector<utype>& p, const utype xmin, const utype xmax, const utype tol) {
    bool LDEBUG = (FALSE || DEBUG_XFIT || XHOST.DEBUG);
    if (LDEBUG) cerr << __AFLOW_FUNC__ << " begin" << std::endl;

    utype left_end = xmin; utype right_end = xmax;
    utype middle = 0.5*(left_end + right_end);

    xvector<utype> dp(2); // index 1 - value, index 2 - first derivative
    evalPolynomialDeriv(left_end, p, dp);
    utype f_at_left_end  = dp[2];
    evalPolynomialDeriv(right_end, p, dp);
    utype f_at_right_end = dp[2];
    evalPolynomialDeriv(middle, p, dp);
    utype f_at_middle    = dp[2];

    // no root within a given interval
    if (sign(f_at_left_end) == sign(f_at_right_end)) {
      string message = "No root within the given interval";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    if (LDEBUG){
      cerr << __AFLOW_FUNC__ << " left_end= "  << left_end  << "middle= ";
      cerr << middle  << "right_end= "  << right_end  << std::endl;
      cerr << __AFLOW_FUNC__ << " f_left= " << f_at_left_end;
      cerr << "f_middle= " << f_at_middle << "f_right= ";
      cerr << f_at_right_end << std::endl;
    }

    // Iterate until the convergence criterion is reached:
    // f(middle) is sufficiently close to zero.
    // Meanwhile, do a sanity check that the function has opposite signs at the interval ends.
    while ((sign(f_at_left_end) != sign(f_at_right_end)) && (std::abs(f_at_middle)>tol)){
      if (sign(f_at_left_end) == sign(f_at_middle)){
        std::swap(left_end, middle);
        std::swap(f_at_left_end, f_at_middle);
      }

      if (sign(f_at_right_end) == sign(f_at_middle)){
        std::swap(right_end, middle);
        std::swap(f_at_right_end, f_at_middle);
      }

      middle = 0.5*(left_end + right_end);
      evalPolynomialDeriv(middle, p, dp);
      f_at_middle = dp[2];

      if (LDEBUG){
        cerr << __AFLOW_FUNC__ << " left_end= "  << left_end  << "middle= ";
        cerr << middle  << "right_end= "  << right_end  << std::endl;
        cerr << __AFLOW_FUNC__ << " f_left= " << f_at_left_end;
        cerr << "f_middle= " << f_at_middle << "f_right= ";
        cerr << f_at_right_end << std::endl;
      }
    }

    // utype-check that the convergence criterion was reached
    if (std::abs(f_at_middle) > tol) {
      string message = "Bisection method did not converge";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    if (LDEBUG) cerr << __AFLOW_FUNC__ << " end" << std::endl;
    return middle;
  }

  //SD20220422
  //Calculate the coefficients for a polynomial p(x) of degree n that fits the y data with weights w
  template<class utype> xvector<utype> polynomialCurveFit(const xvector<utype>& x, const xvector<utype>& _y, const int n, const xvector<utype> _w) {
    bool LDEBUG = (FALSE || DEBUG_XFIT || XHOST.DEBUG);
    xvector<utype> p;
    utype wtot = 0.0;
    for (int i = 1; i <= _w.rows; i++) {
      if (_w(i) < (utype)0.0) {
        stringstream message;
        message << "Negative weight i=" << i;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
      else if (std::isnan(_w(i))) {
        stringstream message;
        message << "NaN weight i=" << i;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
      else if (aurostd::abs(_w(i)) == (utype)INFINITY) {
        stringstream message;
        message << "Inf weight i=" << i;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
      wtot += _w(i);
    }
    if (wtot == 0.0) {
      string message = "Weights cannot be zero";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    xmatrix<utype> VM = Vandermonde_matrix(x, n + 1);
    if (LDEBUG) {cerr << "VM old=" << VM << endl;}
    xvector<utype> w = _w / wtot; // normalize weights
    for (int i = 1; i <= VM.rows; i++) {
      VM.setmat(w(i) * VM.getmat(i, i, 1, VM.cols), i, 1);
    }
    if (LDEBUG) {cerr << "VM new=" << VM << endl;}
    xvector<utype> y = aurostd::elementwise_product(w, _y);
    return aurostd::LinearLeastSquares(VM, y);
  }

  //SD20220318
  //Calculate the (tranpose) n-by-n companion matrix of a polynomial, where n is the polynomial degree
  //transposed companion matrix:
  //  0   1   0  ..  0
  //  0   0   1  ..  0
  // ...
  //  0   0   0  ..  1
  // -p0 -p1 -p2 .. -p(n-1)
  template<class utype> xmatrix<utype> companion_matrix(const xvector<utype>& p) {
    if (p(p.rows) == 0.0) {
      string message = "Leading polynomial coefficient is zero";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    int n = p.rows - 1;
    xmatrix<utype> CM(n, n);
    for (int i = 1; i < n; i++) {
      for (int j = 1; j <= n; j++) {
        if (i == j - 1) {CM(i, j) = 1.0;}
      } 
    }
    for (int j = 1; j <= n; j++) {CM(n, j) = -p(j) / p(n + 1);}
    return CM;
  }

  //SD20220318
  //Calculates the roots of a polynomial as the eigenvalues of the companion matrix
  //The real and imaginary parts of the roots are returned in rr and ri, respectively
  //DOI: 10.1090/S0025-5718-1995-1262279-2
  template<class utype> void polynomialFindRoots(const xvector<utype>& p, xvector<utype>& rr, xvector<utype>& ri) {
    for (int i = 1; i <= p.rows; i++) {
      if (std::isnan(p(i))) {
        stringstream message;
        message << "NaN in polynomial coefficient i=" << i;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
      else if (aurostd::abs(p(i)) == (utype)INFINITY) {
        stringstream message;
        message << "Inf in polynomial coefficient i=" << i;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
    }
    aurostd::eigen(companion_matrix(p), rr, ri);
  }
}

//********************************************************************************
//              Definitions of the NonlinearFit class members
namespace aurostd{
  NonlinearFit::NonlinearFit()
  {
    free();
  }

  NonlinearFit::NonlinearFit(const NonlinearFit &nlf)
  {
    if (this==&nlf) return;
    free(); copy(nlf);
  }

  NonlinearFit::NonlinearFit(xvector<double>& x_in, xvector<double>& y_in,
      xvector<double>& guess_in,
      double foo(const double x, const xvector<double>& p, xvector<double>& dydp),
      double tol_in, double tau_in, int max_iter_in)
  {
    free();
    Npoints = x_in.rows;
    Nparams = guess_in.rows;
    tol = tol_in;
    tau = tau_in;
    max_iter = max_iter_in;
    x = x_in;
    y = y_in;
    residuals = xvector<double> (Npoints);
    guess = guess_in;
    p = xvector<double> (Nparams);
    dydp = xvector<double> (Nparams);
    A = xmatrix<double> (Nparams, Nparams);
    J = xmatrix<double> (Npoints, Nparams);
    f = foo;
  }

  NonlinearFit::~NonlinearFit() { free(); }

  const NonlinearFit& NonlinearFit::operator=(const NonlinearFit &nlf)
  {
    copy(nlf);
    return *this;
  }

  void NonlinearFit::clear() { free(); }
  void NonlinearFit::free()
  {
    Npoints = 0;
    Nparams = 0;
    tol = 0.0;
    tau = 0.0;
    max_iter = 0;
    x.clear();
    y.clear();
    residuals.clear();
    guess.clear();
    p.clear();
    dydp.clear();
    A.clear();
    J.clear();
    f = NULL;
  }

  void NonlinearFit::copy(const NonlinearFit &nlf)
  {
    if (this==&nlf) return;

    Npoints = nlf.Npoints;
    Nparams = nlf.Nparams;
    tol = nlf.tol;
    tau = nlf.tau;
    max_iter = nlf.max_iter;
    x = nlf.x;
    y = nlf.y;
    residuals = nlf.residuals;
    guess = nlf.guess;
    p = nlf.p;
    dydp = nlf.dydp;
    A = nlf.A;
    J = nlf.J;
    f = nlf.f;
  }

  /// Calculates the residual sum of squares of a model function w.r.t the given data
  ///
  double NonlinearFit::calculateResidualSquareSum(const xvector<double> &params)
  {
    double chi_sqr = 0.0;

    calculateResiduals(params); // calculate residuals and save them to residuals
    for (int i=1; i<= Npoints; i++) chi_sqr += std::pow(residuals[i], 2);

    return chi_sqr;
  }

  /// Calculates residuals of a given model function and stores the results
  ///
  void NonlinearFit::calculateResiduals(const xvector<double> &params)
  {
    for (int i=1; i<=Npoints; i++) residuals[i] = y[i] - f(x[i], params, dydp);
  }

  /// Calculates the Jacobian of a given model function for the given "guess" parameters
  ///
  void NonlinearFit::Jacobian(const xvector<double> &guess)
  {
    for (int i=1; i<=Npoints; i++){
      f(x[i], guess, dydp); 
      for (int j=1; j<=Nparams; j++){
        J[i][j] = -dydp[j]; // derivative of a given function w.r.t fit parameters
      }
    }
  }

  /// Nonlinear fit using the Levenberg-Marquardt algorithm.
  /// The implementation here is based on the ideas from Numerical Recipes and
  /// K.Madsen et al. Methods For Non-linear Least Squares Problems
  /// http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
  ///
  bool NonlinearFit::fitLevenbergMarquardt()
  {
    bool LDEBUG = (FALSE || DEBUG_XFIT || XHOST.DEBUG);
    if (LDEBUG) cerr << __AFLOW_FUNC__ << " begin"  << std::endl;

    static const double step_scaling_factor = 10.0;

    xvector<double> pnew(Nparams);
    xmatrix<double> G(Nparams, 1), M(Nparams, Nparams);

    int iter = 0;
    double chi_sqr = 0.0, new_chi_sqr =0.0; // old and new residual square sums

    p = guess; Jacobian(p); calculateResiduals(p); // p are the fitted parameters

    xmatrix<double> A = trasp(J)*J;
    xvector<double> g = -trasp(J)*residuals;

    // determine initial step
    double lambda = tau*maxDiagonalElement(A);
    if (LDEBUG) cerr << __AFLOW_FUNC__ << " lambda = " << lambda << std::endl;

    while (iter<max_iter){
      iter++;
      if (LDEBUG) cerr << __AFLOW_FUNC__ << " iteration: " << iter << std::endl;

      // M = A + lambda*diag(A)
      M = A; for (int i=1; i<=Nparams; i++) M[i][i] += lambda*A[i][i];

      if (LDEBUG) cerr << __AFLOW_FUNC__ << " M:" << std::endl << M << std::endl;

      // transformation: xvector(N) => xmatrix(N,1) to use GaussJordan function
      for (int i=1; i<=Nparams; i++) G[i][1] = g[i];
      aurostd::GaussJordan(M,G);

      pnew = p + M*g;

      chi_sqr = calculateResidualSquareSum(p);
      new_chi_sqr = calculateResidualSquareSum(pnew);

      if (std::abs(new_chi_sqr - chi_sqr) < tol){
        p = pnew;
        break;
      }

      // update only if step leads to a smaller residual sum of squares
      if (new_chi_sqr < chi_sqr){
        p = pnew; Jacobian(p); calculateResiduals(p);

        A = trasp(J)*J;
        g = -trasp(J)*residuals;

        lambda /= step_scaling_factor;
      }
      else{
        lambda *= step_scaling_factor;
      }
      if (LDEBUG) cerr << __AFLOW_FUNC__ << " pnew = " << p << std::endl;
    }

    if (LDEBUG) cerr << __AFLOW_FUNC__ << " end"  << std::endl;
    return iter < max_iter;
  }
}

#endif
// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2021              *
// *                                                                        *
// **************************************************************************
