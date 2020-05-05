// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
// Written by Andriy Smolyanyuk (andriy.smolyanyuk@duke.edu)
//
//    This file is part of Aflow.
//
//    Aflow is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Aflow is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Aflow.  If not, see <https://www.gnu.org/licenses/>.

#include "aflow_apl.h"

#if __cplusplus >= 201103L
template <typename T>
using auto_ptr = std::unique_ptr<T>;
#else
using std::auto_ptr;
#endif

#define SW 5  // width of separator when text is formatted
#define TW 15 // width of label/number when text is formatted
#define DCOEFF 1e-2 // coefficient used in numerical differentiation (should not
                    // be too small). Usage dT = DCOEFF*T

#define EOS_MURNAGHAN       0
#define EOS_POLYNOMIAL      1
#define EOS_BIRCH_MURNAGHAN 2

#define QHA_CALC   0
#define QHA3P_CALC 1
#define SCQHA_CALC 2

//================================================================================

// [TODO:] move this class to aurostd_xmatrix?

/** Fit to a nonlinear model using Levenberg-Marquardt algorithm.
 *
 *  The implementation here is based on the ideas from Numerical Recipes and
 *  K.Madsen et. al. Methods For Non-linear Least Squares Problems
 *  http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
 *
 *  Caution: default value for parameter tau (a scaling factor for initial step size)
 *  was picked to yield a correct fit to Murnaghan equation of state. 
 *  If you observe that this is not a good choice for your function,
 *  try 1e-6 if initial guess is believed to be a good approximation to a true
 *  parameters. Otherwise 1e-3 or even 1 might be a better choice.
 */
class NonlinearFit{
  static const int max_iter=1000;
  public:
    int Npoints, Nparams;
    double tol; /// convergence tolerance criterion
    double tau; /// scaling parameter for initial step size
    xvector<double> x,y;   // data points
    xvector<double> yvec;  // residuals of a given model function
    xvector<double> guess; // initial guess for fit parameters
    xvector<double> p;     // parameters obtained by fit
    xmatrix<double> A;     // J^T.J matrix
    xmatrix<double> J;     // Jacobian of a model function w.r.t parameters
    void (*f)(const double x, xvector<double> &a, double &y, xvector<double> &dyda);
    NonlinearFit(xvector<double> &x, xvector<double> &y, xvector<double> &guess,
        void foo(const double x, xvector<double> &a, double &y, xvector<double> &dyda),
        double tol=1e-6, double tau=1e-12)
        : Npoints(x.rows), Nparams(guess.rows), tol(tol), tau(tau), x(x), y(y),
        yvec(Npoints), guess(guess), p(Nparams), A(Nparams, Nparams),
        J(Npoints, Nparams), f(foo) {}
    void fit();
    void Jacobian(xvector<double> &guess);
    void fvec(xvector<double> &params);
    double err(xvector<double> &params);
    double inf_norm(xvector<double> &x);
};


/** Calculates residual sum of squares of a model function w.r.t given data
 */
double NonlinearFit::err(xvector<double> &params)
{
  double err = 0.0;

  fvec(params); // calculate residuals and save them to yvec
  for (int i=1; i<= Npoints; i++) err += pow(yvec[i], 2);

  return err;
}

/** Calculates residuals of a given model function and stores the result in 
 *  "yvec" array
 */
void NonlinearFit::fvec(xvector<double> &params)
{
  static xvector<double> dyda(Nparams);
  double yval = 0.0;
  for (int i=1; i<=Npoints; i++){
    f(x[i], params, yval, dyda);
    yvec[i] = y[i] - yval;
  }
}

/** Calculates Jacobian of a given model function for "guess" parameters
 */
void NonlinearFit::Jacobian(xvector<double> &guess)
{
  double yval = 0.0;
  xvector<double> dyda(Nparams);
  for (int i=1; i<=Npoints; i++){
    f(x[i], guess, yval, dyda); 
    for (int j=1; j<=Nparams; j++){
      J[i][j] = -dyda[j];
    }
  }
}

/** Nonlinear fit using Levenberg-Marquardt algorithm.
 *  The implementation here is based on the ideas from Numerical Recipes and
 *  K.Madsen et. al. Methods For Non-linear Least Squares Problems
 *  http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
 */
void NonlinearFit::fit()
{
  bool LDEBUG = false;
  if (LDEBUG) cout << "NonlinearFit::fit()\n";

  xvector<double> pnew(Nparams);
  xmatrix<double> G(Nparams, 1), M(Nparams, Nparams);

  int iter = 0;
  double khi = 0.0, new_khi =0.0; // old and new residual square sums

  p = guess; Jacobian(p); fvec(p);
  xmatrix<double> A = trasp(J)*J;
  xvector<double> g = -trasp(J)*yvec;

  bool found = false;

  // find maximum diagonal element
  double maxA = A[1][1];
  for (int i=2; i<=Nparams; i++) if (A[i][i] > maxA) maxA = A[i][i];
  
  double lambda = tau*maxA; // determine initial step
  if (LDEBUG) cout << "lambda = " << lambda << '\n';

  while (!found && iter<max_iter){
    iter++;
    if (LDEBUG) cout << "Iteration: " << iter << '\n';

    // M = A + lambda*diag(A)
    M = A; for (int i=1; i<=Nparams; i++) M[i][i] += lambda*A[i][i];

    if (LDEBUG) cout << "M: \n" << M << '\n';

    // transformation: xvector(N) => xmatrix(N,1) to use GaussJordan function
    for (int i=1; i<=Nparams; i++) G[i][1] = g[i];
    aurostd::GaussJordan(M,G);

    pnew = p + M*g;

    khi = err(p); new_khi = err(pnew);

    if (abs(new_khi-khi) < tol) found = true;

    // update only if step leads to a smaller residual sum of squares
    if (new_khi < khi){
      p = pnew; Jacobian(p); fvec(p);

      A = trasp(J)*J;
      g = -trasp(J)*yvec;

      lambda /= 10.0;
    }
    else{
      lambda *= 10.0;
    }
    if (LDEBUG) cout << "pnew = " << p << '\n';
  }

  if (iter==max_iter){
    pflow::logger("NonlinearFit", "fit()","Maximum number of iterations was reached.",
        cout, _LOGGER_WARNING_);
  }
}

//================================================================================
//                    EOS related

/** Murnaghan equations of state
 * https://en.wikipedia.org/wiki/Murnaghan_equation_of_state
 * Ref: Proceedings of the National Academy of Sciences of the United States of 
 * America, 30 (9): 244–247
 */
void murnaghan(const double x, xvector<double> &a, double &y, xvector<double> &dyda)
{
  double Eeq = a[1];
  double Veq = a[2];
  double B   = a[3];
  double Bp  = a[4];
  double V = x;

  y = Eeq-(B*Veq)/(Bp-1.0)+(B*V*(1+pow(Veq/V,Bp)/(Bp-1)))/Bp;

  dyda[1] = 1.0;
  dyda[2] = (B*(Veq-V*pow(Veq/V,Bp)))/(Veq-Bp*Veq);
  dyda[3] = -(Veq/(Bp-1.0))+(V*(1+pow(Veq/V,Bp)/(Bp-1.0)))/Bp;
  dyda[4] = (B*(-V+2.0*Bp*V-pow(Bp,2)*V+pow(Bp,2)*Veq+V*pow(Veq/V,Bp)*(1.0-2.0*Bp
            +(Bp-1.0)*Bp*log(Veq/V))))/(pow(Bp-1.0,2)*pow(Bp,2));
}

/** Birch-Murnaghan equation of state
 * https://en.wikipedia.org/wiki/Birch%E2%80%93Murnaghan_equation_of_state
 * Ref: Physical Review. 71 (11): 809–824
 */
void birch_murnaghan(const double x, xvector<double> &a, double &y, xvector<double> &dyda)
{
  double Eeq = a[1];
  double Veq = a[2];
  double B   = a[3];
  double Bp  = a[4];
  double V = x;

  double kappa = pow(Veq/V, 2.0/3.0);

  y = Eeq + 9.0/16.0*B*Veq*(pow(kappa-1.0,3.0)*Bp + pow(kappa-1.0,2.0)*(6.0-4.0*kappa));

  dyda[1] = 1.0;
  dyda[2] = 3.0*B*kappa*(kappa-1)*(3.0*(Bp-6.0)*V*pow(Veq/V,1.0/3.0)+
      Veq*(62.0-36.0*kappa+3*Bp*(3*kappa-4)))/(16.0*Veq);
  dyda[3] = 9.0/16.0*Veq*pow(kappa-1,2.0)*(6.0-4.0*kappa+Bp*(kappa-1));
  dyda[4] = 9.0/16.0*B*Veq*pow(kappa-1,3.0);
}


/** Check if minimum exists within a given data set 
 *  (at least one internal point should be lower than edge points) */
bool isMinimumWithinBounds(xvector<double> y){
  bool min_within_bounds = false;
  for (int i=y.lrows+1; i<y.urows; i++){
    if (y[y.lrows] > y[i] && y[y.urows] > y[i]){
      min_within_bounds = true;
      break;
    }
  }

  return min_within_bounds;
}

/** Polynomial model for equation of state fit
 *
 *  E = a + b*V^(-2/3) + c*V^(-4/3) + d*V^(-6/3) + e*V^(-8/3)
 *  Ref: (eq (5.2) page 106)  https://doi.org/10.1017/CBO9781139018265
 */
double EOS_poly(double V, xvector<double> p){
  if (p.rows != 5) return 0;
  return p[1] + p[2]*pow(V,-2/3.0) + p[3]*pow(V,-4/3.0)
              + p[4]*pow(V,-6/3.0) + p[5]*pow(V,-8/3.0);
}

/** First derivative w.r.t volume of polynomial model */
double dEOS_poly(double x, xvector<double> p){
  if (p.rows != 5) return 0;
  return -2/3.0*p[2]*pow(x,-5/3.0) - 4/3.0*p[3]*pow(x, -7/3.0)
         -6/3.0*p[4]*pow(x,-9/3.0) - 8/3.0*p[5]*pow(x,-11/3.0);
}

/** Second derivative w.r.t volume of polynomial model */
double d2EOS_poly(double x, xvector<double> p){
  if (p.rows != 5) return 0;
  return  10/9.0*p[2]*pow(x, -8/3.0) + 28/9.0*p[3]*pow(x,-10/3.0)
         +     6*p[4]*pow(x,-12/3.0) + 88/9.0*p[5]*pow(x,-14/3.0);
}

/** Third derivative w.r.t volume of polynomial model */
double d3EOS_poly(double x, xvector<double> p){
  if (p.rows != 5) return 0;
  return -(80.0/27.0*p[2]*pow(x, -11.0/3.0) + 280.0/27.0*p[3]*pow(x, -13.0/3.0) +
                24.0*p[4]*pow(x, -15.0/3.0) + 88*14/27.0*p[5]*pow(x, -17.0/3.0));
}

/** Calculate equilibrium volume for polynomial model
 *
 *  Veq is determined by searching for a solution of dE/dV=0 equation.
 *
 *  It is assumed that only one minimum exists between Vmin and Vmax,
 *  since it is the only physically meaningful situation
 */
double getEqVolume(xvector<double> p, double Vmin, double Vmax, double tol=_mm_epsilon) {
  bool LDEBUG = false;

  double a = Vmin; double c = Vmax; double b = (a+c)/2;

  double fa = dEOS_poly(a, p);
  double fb = dEOS_poly(b, p);
  double fc = dEOS_poly(c, p);

  if (sign(fa) == sign(fc)){
    pflow::logger("EOS", "getEqVolume","ERROR: no root within the given region.\n",
        cout, _LOGGER_ERROR_);
    return -1;
  }

  if (LDEBUG){
    cout << "a= "  << a  << "b= "  << b  << "c= "  << c  << '\n';
    cout << "fa= " << fa << "fb= " << fb << "fc= " << fc << '\n';
  }

  while ((sign(fa) != sign(fc)) && (abs(fb)>tol)){
    if (sign(fa) == sign(fb)){
      std::swap(a,b);
      std::swap(fa,fb);
    }

    if (sign(fc) == sign(fb)){
      std::swap(c,b);
      std::swap(fc,fb);
    }

    b = (a+c)/2;
    fb = dEOS_poly(b, p);

    if (LDEBUG){
      cout << "a= "  << a  << "b= "  << b  << "c= "  << c  << '\n';
      cout << "fa= " << fa << "fb= " << fb << "fc= " << fc << '\n';
    }
  }

  return b;
}

/** Perform a fit for polynomial EOS model for a given set of volumes and energies
 *  using linear least squares method
 */
xvector<double> fitEOS_poly(xvector<double> &Volumes, xvector<double> &E){
  xmatrix<double> M(Volumes.rows,5);
  xvector<double> sigma(Volumes.rows);
  for (int i=sigma.lrows;i<=sigma.urows;i++) sigma[i]=1.0;

  for (int i=Volumes.lrows; i<=Volumes.urows; i++) {
    for (int j=0; j<=4; j++) M[i][j+1] = pow(Volumes[i],-2*j/3.0);
  }

  aurostd::cematrix M_ce(M);
  M_ce.LeastSquare(E, sigma);
  return M_ce.GetFitVector();
}

/** Calculate bulk modulus for polynomial EOS model */
double BulkModulus(double V0, xvector<double> fit_parameters){
  return eV2GPa*V0*d2EOS_poly(V0,fit_parameters);
}

/** Calculate pressure derivative of bulk modulus for polynomial EOS model */
double Bprime(double V0, xvector<double> fit_parameters){
  return -(V0*d3EOS_poly(V0, fit_parameters)/d2EOS_poly(V0, fit_parameters) + 1);
}

apl::EOSfit::EOSfit()
{
  guess = xvector<double>(4);
  p     = xvector<double>(4);
}

/** Fit equation of state to one of the following models:
 * polynomial: (eq (5.2) page 106)  https://doi.org/10.1017/CBO9781139018265
 * Murnaghan: Proceedings of the National Academy of Sciences of the United States of 
 * America, 30 (9): 244–247
 * Birch-Murnaghan: Physical Review. 71 (11): 809–824
 */
void apl::EOSfit::fit()
{
  switch(method){
    case(EOS_POLYNOMIAL):
      p   = fitEOS_poly(V, E);
      Veq = getEqVolume(p, min(V), max(V));
      Eeq = EOS_poly(Veq, p);
      B   = BulkModulus(Veq, p);
      Bp  = Bprime(Veq, p);
      break;
    case(EOS_BIRCH_MURNAGHAN):
      {
      guess[1] = min(E);
      guess[2] = (max(V)+min(V))/2;
      guess[3] = V[1]*(E[3]-2*E[2]+E[1])/pow(V[2]-V[1],2); // B from central differences
      guess[4] = 3.5; // a reasonable initial value for most materials

      NonlinearFit bmfit(V,E,guess,birch_murnaghan);
      bmfit.fit();
      p   = bmfit.p;
      Eeq = p[1];
      Veq = p[2];
      B   = p[3]*eV2GPa;
      Bp  = p[4];
      }
      break;
    case(EOS_MURNAGHAN):
    default:
      guess[1] = min(E);
      guess[2] = (max(V)+min(V))/2;
      guess[3] = V[1]*(E[3]-2*E[2]+E[1])/pow(V[2]-V[1],2); // B from central differences 
      guess[4] = 3.5; // a reasonable intial value for most materials

      NonlinearFit nlfit(V,E,guess,murnaghan);
      nlfit.fit();
      p = nlfit.p;

      Eeq = p[1];
      Veq = p[2];
      B   = p[3]*eV2GPa;
      Bp  = p[4];

      break;
  }
}

/** Returns energy at a given volume for a chosen EOS model */
double apl::EOSfit::eval(double Vin)
{
  double res = 0;

  static xvector<double> dyda(4);
  switch(method){
    case(EOS_POLYNOMIAL):
      res = EOS_poly(Vin, p);
      break;
    case(EOS_BIRCH_MURNAGHAN):
      birch_murnaghan(Vin, p, res, dyda);
      break;
    case(EOS_MURNAGHAN):
    default:
      murnaghan(Vin, p, res, dyda);
      break;
  }

  return res;
}


//=============================================================================
//                         Definitions of QHA class members
namespace apl
{
  QHAN::QHAN() { free(); }
  QHAN::QHAN(const QHAN &qha){
    free();
    copy(qha);
  }
  QHAN::~QHAN() { free(); }
  const QHAN& QHAN::operator=(const QHAN &qha){
    copy(qha);
    return *this;
  }

  void QHAN::clear()
  {
    QHAN qha;
    copy(qha);
  }

  /** Initializes all values to "zero" and attempts to clear all containers.
   */
  void QHAN::free()
  {
    messageFile = NULL;
    oss = &cout;
    supercellopts.clear();
    isEOS = false; isGP_FD = false;
    ignore_imaginary = false;
    runQHA   = false; runQHA3P = false; runSCQHA = false;
    includeElectronicContribution = false;
    Ntemperatures = 0;
    N_GPvolumes = 3;
    N_EOSvolumes = 0;
    Nbranches = 0;
    NatomsOrigCell = 0;
    origStructure.clear();
    Temperatures.clear();
    GPvolumes.clear();
    EOSvolumes.clear();
    coefGPVolumes.clear();
    coefEOSVolumes.clear();
    DOS_Ef.clear();
    Efermi_V.clear();
    E0_V.clear();
    static_eigvals.clear();
    static_ibzkpts.clear();
    energies_V.clear();
    edos_V.clear();
    frequencies_V.clear();
    pdos_V.clear();
    qpWeights.clear();
    qPoints.clear();
    gp_fit_matrix.clear();
    omegaV_mesh.clear();
    omegaV_mesh_EOS.clear();
    gp_ph_dispersions.clear();
    eos_vib_thermal_properties.clear();
    subdirectories_apl_gp.clear();
    subdirectories_apl_eos.clear();
    subdirectories_static.clear();
    arun_runnames_apl_gp.clear();
    arun_runnames_apl_eos.clear();
    arun_runnames_static.clear();
    xinput.clear();
    xflags.clear();
    aflags.clear();
    kflags.clear();
    currentDirectory = ".";
  }

  void QHAN::copy(const QHAN &qha)
  {
    if (this==&qha) return;

    apl_options       = qha.apl_options;
    eos               = qha.eos;
    messageFile       = qha.messageFile;
    oss               = qha.oss;
    supercellopts     = qha.supercellopts;
    isEOS             = qha.isEOS;
    isGP_FD           = qha.isGP_FD;
    ignore_imaginary  = qha.ignore_imaginary;
    runQHA            = qha.runQHA;
    runQHA3P          = qha.runQHA3P;
    runSCQHA          = qha.runSCQHA;
    includeElectronicContribution = qha.includeElectronicContribution;
    Ntemperatures     = qha.Ntemperatures;
    N_GPvolumes       = qha.N_GPvolumes;
    N_EOSvolumes      = qha.N_EOSvolumes;
    Nbranches         = qha.Nbranches;
    NatomsOrigCell    = qha.NatomsOrigCell;
    origStructure     = qha.origStructure;
    Temperatures      = qha.Temperatures;
    GPvolumes         = qha.GPvolumes;
    EOSvolumes        = qha.EOSvolumes;
    coefGPVolumes     = qha.coefGPVolumes;
    coefEOSVolumes    = qha.coefEOSVolumes;
    DOS_Ef            = qha.DOS_Ef;
    Efermi_V          = qha.Efermi_V;
    E0_V              = qha.E0_V;
    static_eigvals    = qha.static_eigvals;
    static_ibzkpts    = qha.static_ibzkpts;
    energies_V        = qha.energies_V;
    edos_V            = qha.edos_V;
    frequencies_V     = qha.frequencies_V;
    pdos_V            = qha.pdos_V;
    qpWeights         = qha.qpWeights;
    qPoints           = qha.qPoints;
    gp_fit_matrix     = qha.gp_fit_matrix;
    omegaV_mesh       = qha.omegaV_mesh;
    omegaV_mesh_EOS   = qha.omegaV_mesh_EOS;
    gp_ph_dispersions = qha.gp_ph_dispersions;
    eos_vib_thermal_properties = qha.eos_vib_thermal_properties;
    subdirectories_apl_eos  = qha.subdirectories_apl_eos;
    subdirectories_apl_gp   = qha.subdirectories_apl_gp;
    subdirectories_static   = qha.subdirectories_static;
    arun_runnames_apl_eos   = qha.arun_runnames_apl_eos;
    arun_runnames_apl_gp    = qha.arun_runnames_apl_gp;
    arun_runnames_static    = qha.arun_runnames_static;
    xinput            = qha.xinput;
    xflags            = qha.xflags;
    aflags            = qha.aflags;
    kflags            = qha.kflags;
    currentDirectory  = qha.currentDirectory;
  }

///////////////////////////////////////////////////////////////////////////////

  QHAN::QHAN(string &tpt, _xinput &xinput, _aflags &aflags, _kflags &kflags,
      _xflags &xflags, string *aflowin, xoption &supercellopts, ofstream &messageFile, 
              ostream &oss)
  {
    initialize(tpt, xinput, aflags, kflags, xflags, aflowin, supercellopts, messageFile,
        oss);
  }

  /** Used to initialize QHAN class with all the necessary data.
   */
  void QHAN::initialize(string &tpt, _xinput &xinput, _aflags &aflags, _kflags &kflags,
          _xflags &xflags, string *aflowin, xoption &supercellopts,
          ofstream &messageFile, ostream &oss)
  {
    string function = "QHAN::initialize()";
    string msg  = "Initializing QHA.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, messageFile, oss, 
        _LOGGER_MESSAGE_);

    free();

    this->oss = &oss; this->messageFile = &messageFile;

    this->xinput = xinput; this->xflags = xflags;
    this->kflags = kflags; this->aflags = aflags;

    this->supercellopts = supercellopts;
    this->aflowin = aflowin;

    currentDirectory = xinput.xvasp.Directory; // remember current directory

    // initialization of data that is derived from the structure info
    origStructure = xinput.xvasp.str;
    origStructure.ReScale(1.0);
    NatomsOrigCell = origStructure.atoms.size();
    Nbranches = NatomsOrigCell * 3;
    double Volume = origStructure.GetVolume();

    // parse QHA-related aflow.in options
    string dirname;
    double gp_distortion;
    vector<double> eosrange(3);
    vector<string> tokens;
    vector<xoption>::iterator option;
    for (option  = kflags.KBIN_MODULE_OPTIONS.qhaflags.begin();
         option != kflags.KBIN_MODULE_OPTIONS.qhaflags.end(); ++option){
      if (option->keyword=="EOS") isEOS = option->option;
      if (option->keyword=="INCLUDE_ELE") includeElectronicContribution = option->option;
      if (option->keyword=="GP_FINITE_DIFF") isGP_FD = option->option;
      if (option->keyword=="IGNORE_IMAGINARY") ignore_imaginary = option->option;

      if (option->keyword=="EOS_DISTORTION_RANGE"){
        apl::tokenize(option->content_string, tokens, string(" :"));

        if (tokens.size() != 3) {
          string message = "Wrong setting in the ";
          message += "[AFLOW_QHA]EOS_DISTORTION_RANGE.\n";
          message += "Specify as EOS_DISTORTION_RANGE=-3:6:1.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, message,
              _INPUT_NUMBER_);
        }

        for (uint j=0; j<3; j++){
          eosrange[j] = aurostd::string2utype<double>(tokens[j]);
        }
      } // eos_distortion_range

      if (option->keyword=="GP_DISTORTION"){
        gp_distortion = option->content_double/100.0;
      }

      if (option->keyword=="MODE"){
        tokens.clear();
        aurostd::string2tokens(option->content_string, tokens, ",");

        if ((tokens.size()<1) || (tokens.size() > 3)){
          string message = "Wrong setting in the ";
          message += "[AFLOW_QHA]MODE.\n";
          message += "Specify as MODE=QHA,QHA3P,SCQHA";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, message,
              _INPUT_NUMBER_);
        }

        for (uint i=0; i<tokens.size(); i++){
          if (tokens[i].length()==3){
            if (tokens[i].find("QHA")!=std::string::npos) runQHA = true;
          }
          else if (tokens[i].length()==5){
            if (tokens[i].find("QHA3P")!=std::string::npos) runQHA3P = true;
            if (tokens[i].find("SCQHA")!=std::string::npos) runSCQHA = true;
          }
        }
      }
    }

    // determine names for directories for the calculation of Grueneisen parameter
    // (calculated using finite differences)
    vector<double> gprange(3);
    gprange[0] = 1.0-gp_distortion; gprange[1] = 1.0; gprange[2] = 1.0+gp_distortion;
    N_GPvolumes = 3;
    for (int j=0; j<N_GPvolumes; j++){
      arun_runnames_apl_gp.push_back("PHONON_" + apl::stringify(gprange[j]));
      dirname = ARUN_DIRECTORY_PREFIX + QHA_ARUN_MODE + '_' +
                arun_runnames_apl_gp.back();
      subdirectories_apl_gp.push_back(dirname);

      coefGPVolumes.push_back(gprange[j]);
      GPvolumes.push_back(gprange[j]*Volume/NatomsOrigCell);
    }

    // determine names for directories for EOS calculation
    if (isEOS){
      eosrange[0] = 1.0 + eosrange[0]/100.0;
      eosrange[1] = 1.0 + eosrange[1]/100.0;
      eosrange[2] = eosrange[2]/100.0;

      N_EOSvolumes = floor((eosrange[1]-eosrange[0])/eosrange[2])+1;

      // get a set of volumes that would be used for QHA-EOS calculation
      for (double i=eosrange[0]; i<=eosrange[1]; i+=eosrange[2]){
        arun_runnames_apl_eos.push_back("PHONON_" + apl::stringify(i));
        dirname = ARUN_DIRECTORY_PREFIX + QHA_ARUN_MODE + '_' +
          arun_runnames_apl_eos.back();
        subdirectories_apl_eos.push_back(dirname);

        arun_runnames_static.push_back("STATIC_" + apl::stringify(i));
        dirname = ARUN_DIRECTORY_PREFIX + QHA_ARUN_MODE + '_' +
          arun_runnames_static.back();
        subdirectories_static.push_back(dirname);

        coefEOSVolumes.push_back(i);
        EOSvolumes.push_back(i*Volume/NatomsOrigCell);
      }
      eos.V = aurostd::vector2xvector(EOSvolumes);
      eos.E = xvector<double>(N_EOSvolumes);
    }

    tokens.clear();
    apl::tokenize(tpt, tokens, string (" :"));
    if (tokens.size() != 3){
      msg = "Wrong setting in the ";
      msg += "[AFLOW_APL]TPT.\n";
      msg += "Specify as TPT=0:2000:10.\n";
      msg += "See README_AFLOW_APL.TXT for the details.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg,
          _INPUT_NUMBER_);
    }
    double tp_start = aurostd::string2utype<double>(tokens[0]);
    double tp_end   = aurostd::string2utype<double>(tokens[1]);
    double tp_step  = aurostd::string2utype<double>(tokens[2]);


    // define a set of temperatures for thermodynamic calculations
    Ntemperatures = floor((tp_end - tp_start)/tp_step) + 1;

    Temperatures = vector<double> (Ntemperatures);
    uint Tcount = 0;
    for (int T=tp_start; T<=tp_end; T+=tp_step){
      Temperatures[Tcount] = T;
      Tcount++;
    }

    // this matrix will be used in the fit of frequency-volume dependency
    // w(V) = a + b*V + c*V^2 + d*V^3
    gp_fit_matrix = xmatrix<double> (N_EOSvolumes, 4);
    for (int Vid=1; Vid<=N_EOSvolumes; Vid++){
      gp_fit_matrix[Vid][1] = 1;
      gp_fit_matrix[Vid][2] = EOSvolumes[Vid-1];
      gp_fit_matrix[Vid][3] = pow(EOSvolumes[Vid-1],2);
      gp_fit_matrix[Vid][4] = pow(EOSvolumes[Vid-1],3);
    }
  }

  /** Perform QHA calculation
   *
   * For a regular QHA calculation there is a choice between two possibilities:
   * 1) calculate Grueneisen parameter using finite difference method.
   *    This calculation requires 3 phonon calculations (for volumes V, V-dV and V+dV);
   * 2) calculate temperature-dependent parameters (such as equilibrium volume,
   * free energy, bulk modulus, thermal expansion, isochoric and isobaric specific heat,
   * average Grueneisen parameter) employing a set of phonon calculations and
   * making a fit to some model equation of state.
   * It requires N static and N phonon calculations.
   *
   * QHA3P and SCQHA calculations require 3 phonon calculations and N static
   * calculations.
   */
  void QHAN::run()
  {
    bool LDEBUG = false | XHOST.DEBUG;

    string function = "QHAN::run()";
    string msg = "Performing QHA calculation.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);

    try{
      bool eos_data_available = false;
      bool gp_data_available = false;

      // QHA3P and SCQHA require static EOS calculation
      // For QHA calculation EOS calculation is to be set with EOS flag
      if (isEOS || runQHA3P || runSCQHA){
        if (check_static() == N_EOSvolumes){
          read_static();
        }
        else{
          createSubdirectoriesStaticRun();
        }
      }

      // In QHA calculation EOS flag enables APL calculation for a set of volumes.
      // This flag is used when one is interested in T-dependent properties.
      if (isEOS){
        eos_data_available = run_apl(subdirectories_apl_eos,
             arun_runnames_apl_eos, coefEOSVolumes, false);

        if (eos_data_available){
          if (includeElectronicContribution) DOSatEf();
          if (LDEBUG) writeFrequencies();
          writeFVT();

          writeThermalProperties(EOS_POLYNOMIAL, QHA_CALC);
          writeThermalProperties(EOS_MURNAGHAN, QHA_CALC);
          writeThermalProperties(EOS_BIRCH_MURNAGHAN, QHA_CALC);
        }
      }

      // QHA3P and SCQHA require that Grueneisen parameter is calculated via
      // finite difference numerical derivative.
      // For QHA calculation use GP_FD flag to turn on this type of calculation
      if (isGP_FD || runQHA3P || runSCQHA){
        gp_data_available = run_apl(subdirectories_apl_gp, arun_runnames_apl_gp,
            coefGPVolumes);

        if (gp_data_available){
          if (runQHA || runQHA3P){
            double V0 = origStructure.GetVolume()/NatomsOrigCell;
            writeGPpath(V0);
            writeGPmeshFD();
            writeAverageGP_FD();

            if (runQHA3P){
              writeThermalProperties(EOS_POLYNOMIAL, QHA3P_CALC);
              writeThermalProperties(EOS_MURNAGHAN, QHA3P_CALC);
              writeThermalProperties(EOS_BIRCH_MURNAGHAN, QHA3P_CALC);
            }
          }

          if (runSCQHA){
            run_scqha(EOS_POLYNOMIAL, true);
            run_scqha(EOS_POLYNOMIAL, false);
          }
        }
      }
    }
    catch(aurostd::xerror &e){
      pflow::logger(e.whereFileName(),e.whereFunction(), e.what(), currentDirectory,
          *messageFile, *oss, _LOGGER_ERROR_);
      return;
    }
  }

  /** Creates subdirectories with aflow.in for static DFT calculations.
   */
  void QHAN::createSubdirectoriesStaticRun()
  {
    // prepare static_bands calculations to get reasonable electronic DOS
    xinput.xvasp.AVASP_flag_RUN_STATIC_BANDS       = true;
    xinput.xvasp.AVASP_flag_RUN_STATIC             = false;
    xinput.xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS = false;
    xinput.xvasp.AVASP_flag_RUN_RELAX_STATIC       = false;
    xinput.xvasp.AVASP_flag_RUN_RELAX              = false;
    xinput.xvasp.AVASP_flag_GENERATE               = false;
    xinput.xvasp.aopts.pop_attached("AFLOWIN_FLAG::MODULE");

    stringstream aflowin;

    for (uint i=0; i<subdirectories_static.size(); i++){
      if (aurostd::FileExist(subdirectories_static[i]+"/aflow.in")) continue;

      xinput.xvasp.str = origStructure;
      xinput.xvasp.str.InflateVolume(coefEOSVolumes[i]);

      xinput.setDirectory(currentDirectory);
      xinput.xvasp.AVASP_arun_mode = QHA_ARUN_MODE;
      xinput.xvasp.AVASP_arun_runname = arun_runnames_static[i];

      if (xinput.xvasp.AVASP_path_BANDS.empty()) xinput.xvasp.AVASP_path_BANDS = "AUTO";
      if (!xinput.xvasp.AVASP_value_BANDS_GRID)
        xinput.xvasp.AVASP_value_BANDS_GRID=DEFAULT_BANDS_GRID;
      AVASP_populateXVASP(aflags,kflags,xflags.vflags, xinput.xvasp);
      AVASP_MakeSingleAFLOWIN(xinput.xvasp, aflowin, true);
    }
  }

  /** Check if all static calculations needed to proceed with QHA calculation
   * exist and returns the number of finished calculations.
   *
   * If not all calculations are finished returns 0.
   */
  int QHAN::check_static()
  {
    string function = "QHAN::check_static()";
    int count = 0;
    bool all_calculations_done = true;
    string dosfile, outcarfile;
    for (uint i=0; i<subdirectories_static.size(); i++){
      if (aurostd::IsDirectory(subdirectories_static[i])){
        count++;
      }
      else{
        string msg = "QHA is not able to proceed: ";
        msg += subdirectories_static[i] + " is missing.\n";
        pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
           _LOGGER_WARNING_);
        return 0;
      }
    }

    if (!all_calculations_done) count = 0;

    return count;
  }

  /** Create APL subdirectory or gather and process data from finished
   *  calculation.
   */
  double QHAN::run_apl(vector<string> &subdirectories,
     vector<string> &arun_subdirectories,
     vector<double> &coefVolumes, bool gp)
  {
    string function = "QHAN::run_apl()";
    string msg = "Reading phonon DOS and dispersion relations.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);

    int Nqpoints = 0;
    bool apl_data_calculated = true;

    for (uint i=0; i<subdirectories.size(); i++){
      xinput.xvasp.str = origStructure;
      xinput.xvasp.str.InflateVolume(coefVolumes[i]);
      // update cartesian coordinates
      for (uint at=0; at<xinput.xvasp.str.atoms.size(); at++){
        xinput.xvasp.str.atoms.at(at).cpos = F2C(xinput.xvasp.str.scale,
            xinput.xvasp.str.lattice, xinput.xvasp.str.atoms.at(at).fpos);
      }

      apl::PhononCalculator phcalc(*messageFile, *oss);
      phcalc.initialize_supercell(xinput.getXStr());
      phcalc.getSupercell().build(supercellopts);

      auto_ptr<apl::ForceConstantCalculator> fccalc;
      apl::DirectMethodPC* dmPC = new apl::DirectMethodPC(phcalc.getSupercell(),
          *messageFile, *oss);

      apl::PhononDispersionCalculator pdisc(phcalc);

      fccalc.reset(dmPC);

      // set options for direct method phonon calculations
      dmPC->setGeneratePlusMinus(apl_options.flag("AUTO_DIST"),apl_options.flag("DPM"));
      dmPC->setGenerateOnlyXYZ(apl_options.flag("XYZONLY"));
      dmPC->setDistortionSYMMETRIZE(apl_options.flag("SYMMETRIZE"));
      dmPC->setDistortionINEQUIVONLY(apl_options.flag("INEQUIVONLY"));
      dmPC->setDistortionMagnitude(aurostd::string2utype<double>(
          apl_options.getattachedscheme("DIST_MAGNITUDE")));
      dmPC->setCalculateZeroStateForces(apl_options.flag("ZEROSTATE"));
      fccalc->setPolarMaterial(apl_options.flag("POLAR"));

      // set name of subdirectory
      xinput.setDirectory(currentDirectory);
      xinput.xvasp.AVASP_arun_mode = QHA_ARUN_MODE;
      xinput.xvasp.AVASP_arun_runname = arun_subdirectories[i];

      AVASP_populateXVASP(aflags,kflags,xflags.vflags, xinput.xvasp);

      // determine APL subdirectories: if it is first run skip to the next
      // volume calculation
      if (fccalc->runVASPCalculations(xinput, aflags, kflags, xflags, *aflowin, false)){
        apl_data_calculated = false;
        continue;
      }

      // check if APL calculation is ready and calculate force constants
      if (!fccalc->run()){
        apl_data_calculated = false;
        continue;
      }
      phcalc.setHarmonicForceConstants(*fccalc);

      // calculate all phonon-related data: DOS, frequencies at q-mesh and
      // phonon dispersions
      vector<xvector<double> > dummy_dos_projections; // do not need in QHA
      vector<int> dos_mesh(3);

      vector<string> tokens;

      apl::tokenize(apl_options.getattachedscheme("DOS_MESH"), tokens, string(" xX"));
      for (uint j=0; j<tokens.size(); j++){
        dos_mesh[j] = aurostd::string2utype<int>(tokens[j]);
      }

      phcalc.initialize_qmesh(dos_mesh);

      apl::DOSCalculator dosc(phcalc, apl_options.getattachedscheme("DOS_METHOD"),
           dummy_dos_projections);
      dosc.calc(aurostd::string2utype<double>(apl_options.getattachedscheme("DOS_NPOINTS")),
          aurostd::string2utype<double>(apl_options.getattachedscheme("DOS_SMEAR")));
      dosc.writePHDOSCAR(subdirectories[i]);

      if (dosc.hasNegativeFrequencies()){
        msg = "Imaginary frequencies for APL calculation in ";
        msg += subdirectories[i] + "folder were detected.\n";
        if (ignore_imaginary){
          msg += "They will be ignored. Check if results are still meaningful!\n";
          pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile,
              *oss, _LOGGER_WARNING_);
        }
        else{
          msg += "QHA calculation will be stopped after checking all available APL calculations.\n";
          msg += "Please check phonon DOS and band structure.\n";
          msg += "Workaround: adjust EOS_DISTORTION_RANGE to exclude problematic calculations or\n";
          msg += "ignore this error by setting IGNORE_IMAGINARY=ON.\n";

          pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile,
              *oss, _LOGGER_ERROR_);
          apl_data_calculated = false;
        }
      }

      string USER_DC_INITLATTICE="";
      int USER_DC_NPOINTS = aurostd::string2utype<int>(
          apl_options.getattachedscheme("BAND_NPOINTS"));
      pdisc.initPathLattice(USER_DC_INITLATTICE, USER_DC_NPOINTS);
      pdisc.calc(apl::THZ | apl::ALLOW_NEGATIVE);
      pdisc.writePHEIGENVAL(subdirectories[i]);

      // save all the data that is necessary for QHA calculations
      if (i==0){// qmesh data is the same for all volumes: need to save it only once
        qpWeights = phcalc.getQMesh().getWeights();
        qPoints   = phcalc.getQMesh().getIrredQPointsFPOS();
        Nqpoints  = qPoints.size();
      }

      // we have two different sets of data: one for finite difference Grueneisen
      // parameter calculation and one for EOS APL calculation
      if (gp){
        gp_ph_dispersions.push_back(pdisc.createEIGENVAL());

        // allocate memory at first run (entire chunk of memory is needed because
        // data will be reordered)
        if (i==0){
          omegaV_mesh = vector<vector<vector<double> > > (Nqpoints,
            vector<vector<double> >(Nbranches, vector<double>(N_GPvolumes)));
        }

        vector<xvector<double> > freqs = dosc.getFreqs();

        for (int q=0; q<Nqpoints; q++){
          for (int branch=0; branch<Nbranches; branch++){
            // clean imaginary frequencies in case if user wants to ignore them
            omegaV_mesh[q][branch][i] = freqs[q][branch+1] < 0 ? 0 : freqs[q][branch+1];
          }
        }
      }
      else {
        eos_vib_thermal_properties.push_back(ThermalPropertiesCalculator(dosc,
              *messageFile));

        // allocate memory at first run (entire chunk of memory is needed because
        // data will be reordered)
        if (i==0){
          omegaV_mesh_EOS = vector<vector<vector<double> > > (Nqpoints,
            vector<vector<double> >(Nbranches, vector<double>(N_EOSvolumes)));
        }

        vector<xvector<double> > freqs = dosc.getFreqs();

        for (int q=0; q<Nqpoints; q++){
          for (int branch=0; branch<Nbranches; branch++){
            // clean imaginary frequencies in case if user wants to ignore them
            omegaV_mesh_EOS[q][branch][i] = freqs[q][branch+1] < 0 ? 0 : freqs[q][branch+1];
          }
        }
      }
    }

    return apl_data_calculated;
  }

  /** Read data from static calculations. */
  void QHAN::read_static()
  {
    string function = "QHAN::read_static()";
    string msg = "";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);

    xDOSCAR doscar;
    xOUTCAR outcar;
    string  outcarfile, dosfile;
    for (uint i=0; i<subdirectories_static.size(); i++){
      msg = "Reading static data from " + subdirectories_static[i];
      pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);

      outcarfile = subdirectories_static[i]+'/'+"OUTCAR.static";
      outcar.GetPropertiesFile(outcarfile);

      dosfile = subdirectories_static[i]+'/'+"DOSCAR.static";
      doscar.GetPropertiesFile(dosfile);

      vector<double> edos(doscar.number_energies);
      vector<double> energies = aurostd::deque2vector(doscar.venergy);

      for (uint j=0; j<doscar.number_energies; j++){
        // sum spin contributions if we deal with magnetic calculation
        for (uint s=0; s<doscar.vDOS.at(0).at(0).size(); s++){
          // if Methfessel-Paxton method was used some of the DOS values might
          // be negative => fix them to be zero
          if (doscar.vDOS[0][0][s][j]>=0){
            edos[j] += doscar.vDOS[0][0][s][j]/outcar.natoms;
          }
        }
      }

      edos_V.push_back(edos);
      energies_V.push_back(energies);
      Efermi_V.push_back(doscar.Efermi);
      E0_V.push_back(outcar.energy_cell/outcar.natoms);

      static_eigvals.push_back(xEIGENVAL(subdirectories_static[i]+'/'+"EIGENVAL.static"));
      static_ibzkpts.push_back(xIBZKPT(subdirectories_static[i]+'/'+"IBZKPT.static"));
    }
  }

  /** Calculates Grueneisen parameter of an individual vibrational mode for a
   * given volume.
   *
   * gamma_i = -V/w*dw/dV,
   * where w is frequency of a given mode a at given volume
   *
   * Volume dependence of w is approximated by the polynomial
   * w = a + b*V  + c*V**2 + d*V**3
   *
   * @param V is a volume at which Grueneisen parameter is calculated.
   * @param xomega is an array of frequency-volume dependence for a given
   * phonon branch.
   */
  double QHAN::calcGamma(double V, xvector<double> &xomega, double &w)
  {
    bool LDEBUG = false | XHOST.DEBUG;
    // no weights in fit
    xvector<double> s(xomega.rows); for (int i=s.lrows;i<=s.urows;i++) s[i]=1;
    aurostd::cematrix lsfit(gp_fit_matrix);

    double a,b,c,d; //  fit parameters
    double gamma = 0;
    w = 0;

    /** Workaround: for w=0, gamma_i is assumed to be zero, since
    / the actual value depends on the path used to approach w->0. */
    bool well_behavied = true;
    for (int i=xomega.lrows; i<=xomega.urows; i++){
      if (!(xomega[i] > _mm_epsilon)){
        well_behavied = false;
        break;
      }
    }

    double err = 0;

    // calculate gamma_i only if w is nonzero
    if (well_behavied){
      lsfit.LeastSquare(xomega,s);
      a = lsfit.AVec()[0];
      b = lsfit.AVec()[1];
      c = lsfit.AVec()[2];
      d = lsfit.AVec()[3];

      if (LDEBUG){
        xvector<double> tmp(4);
        tmp[1] = a; tmp[2] = b; tmp[3] = c; tmp[4] = d;

        // check fit error
        for (int Vid=0; Vid<N_EOSvolumes; Vid++){
          err = abs(a+b*EOSvolumes[Vid]+c*pow(EOSvolumes[Vid],2)+d*pow(EOSvolumes[Vid],3)
                -xomega[Vid+1])/xomega[Vid+1];
          err *= 100.0;
          if (err>=10.0){
            string msg="Relative error of log(w)=f(V) fit (used to ";
            msg+="determine mode-decomposed Grueneisen parameter) is larger than ";
            msg+="10\% for V="+apl::stringify(EOSvolumes[Vid])+'\n';

            pflow::logger(QHA_ARUN_MODE, "calcGamma()", msg, currentDirectory,
                *messageFile, *oss, _LOGGER_MESSAGE_);
            cout << xomega << '\n';
            cout << gp_fit_matrix * tmp << '\n';
            break;
          }
        }
      }

      w = a + b*V + c*pow(V,2) + d*pow(V,3);
      gamma = -V/w * (b + 2*c*V + 3*d*pow(V,2));
    }

    return gamma;
  }

  /** Calculates Grueneisen parameter of an individual vibrational mode using
   * central finite difference method.
   *
   * gamma = -V0/w0*dw/dV |V->V0 ~ -V0/w0*(w(V0+dV)-w(V0-dV))/(2*dV)
   */
  double QHAN::calcGammaFD(xvector<double> &xomega)
  {
    if (xomega.rows != 3) return 0.0;

    if (xomega[2] > 0){
      return -GPvolumes[1]/xomega[2]*(xomega[3]-xomega[1])/
        (GPvolumes[2]-GPvolumes[0]);
    }
    else{
      return 0.0;
    }
  }


  /** Calculates free energy (without electronic contribution) as a function of volume
   *  and temperature.
   *  Volume dependency is obtained via fit to QHA model.
   *  @param qha_method defines what type of QHA calculation is performed.
   */
  double QHAN::FreeEnergyFit(double T, double V, int qha_method)
  {
    if (T<0) return 0;

    switch(qha_method){
      case (QHA3P_CALC):
        for (int i=eos.E.lrows; i<=eos.E.urows; i++){
          eos.E[i] = FreeEnergyTE(T, i-1);
        }
        break;
      case(QHA_CALC):
      default:
        for (int i=eos.E.lrows; i<=eos.E.urows; i++){
          eos.E[i] = FreeEnergy(T, i-1);
        }
        break;
    }
    if (includeElectronicContribution) eos.E += electronicFreeEnergySommerfeld(T);
    eos.fit();
    return eos.eval(V);
  }

  /** Calculates internal energy as a function of volume and temperature.
   */
  double QHAN::InternalEnergyFit(double T, double V)
  {
    static xvector<double> U(N_EOSvolumes);
    static xvector<double> xvolumes = aurostd::vector2xvector(EOSvolumes);

    if (T<0) return 0;

    for (int i=0; i<N_EOSvolumes; i++){
      U[i+1] = eos_vib_thermal_properties[i].getInternalEnergy(T, apl::eV);
      U[i+1] /= NatomsOrigCell;
    }

    return EOS_poly(V, fitEOS_poly(xvolumes, U));
  }

  /** Free energy (without electronic contribution) for a calculation at specific
   *  volume (given by id).
   */
  double QHAN::FreeEnergy(double T, int id)
  {
    return eos_vib_thermal_properties[id].getVibrationalFreeEnergy(T, apl::eV)
           /NatomsOrigCell + E0_V[id];
  }

  double FermiDirac(double E, double mu, double T){
    if (T<0) return 0;

    if (T>0){
      return 1/(1+exp((E-mu)/(KBOLTZEV*T)));
    }
    else
      // At T=0 FD transforms to Heaviside step function
      if (E<mu){
        return 1;
      }
      else if (E>mu){
        return 0;
      }
      else return 0.5;
  }

  /** Calculates electronic free energy using integration over DOS.
   *  Currently this function does not work (to be fixed in the future)
   *
   *  @param id selects the volume at which corresponding data was obtained 
   */
  double QHAN::electronicFreeEnergy(double T, int id)
  {
    double EelecT = 0; double SelecT = 0; double FelecT = 0; double Eelec0K = 0;
    double f, f0;

    double dE = (max(energies_V[id])-min(energies_V[id]))/(energies_V[id].size()-1);
    for (uint i=0; i<edos_V[id].size(); i++){
      f0 = FermiDirac(energies_V[id][i], Efermi_V[id], 0);
      f  = FermiDirac(energies_V[id][i], Efermi_V[id], T);

      Eelec0K += energies_V[id][i] * edos_V[id][i] * f0;
      EelecT  += energies_V[id][i] * edos_V[id][i] * f;
      // limit of values for f->0 or f->1 is 0
      if (f>0 && f<1) SelecT -= (f*log(f) + (1-f)*log(1-f)) * edos_V[id][i];
    }

    EelecT  *= dE;
    Eelec0K *= dE;
    SelecT  *= KBOLTZEV*dE;

    FelecT = (EelecT-Eelec0K) - T*SelecT;

    return FelecT;
  }

  /** Calculate equilibrium volume for a given temperature.
   * @param qha_method defines what type of QHA calculation is performed.
   */
  double QHAN::getEqVolumeT(double T, int qha_method)
  {
    switch(qha_method){
      case (QHA3P_CALC):
        for (int i=eos.E.lrows; i<=eos.E.urows; i++){
          eos.E[i] = FreeEnergyTE(T, i-1);
        }
        break;
      case(QHA_CALC):
      default:
        for (int i=eos.E.lrows; i<=eos.E.urows; i++){
          eos.E[i] = FreeEnergy(T, i-1);
        }
        break;
    }

    if (includeElectronicContribution) eos.E += electronicFreeEnergySommerfeld(T);
    eos.fit();

    return eos.Veq;
  }


  /** Calculates average Grueneisen parameter (GP) and isochoric specific heat (CV)
   *  using weighted sum over k-points mesh.
   *  This function is used when Grueneisen parameter is calculated using finite
   *  difference method.
   */
  void QHAN::calcCV_GP(double T, double &CV, double &GP)
  {
    CV = 0; GP = 0;
    if (!(T>0)) return;

    xvector<double> xomega;

    uint NIrQpoints = omegaV_mesh.size();
    int NQpoints = 0; // the total number of q-points

    double Cvi  = 0;   // mode-dependent specific heat at V=const
    double w    = 0;   // frequency for a given volume
    double expx = 0;   // temperature-dependent exponential factor

    for (uint q=0; q<NIrQpoints; q++){
      for (uint branch=0; branch<omegaV_mesh[q].size(); branch++){
          xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);

          w = xomega[2]*THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]
          if (w>0){
            expx = exp(w/KBOLTZEV/T);

            Cvi = pow(w,2)*expx/pow(expx-1.0,2) * qpWeights[q];

            GP += calcGammaFD(xomega) * Cvi;
            CV += Cvi;
          }
        }
        NQpoints += qpWeights[q];
      }

    GP /= CV;
    CV /= NQpoints; CV /= Nbranches;
    CV *= 3*pow(1/KBOLTZEV/T,2);
  }

  /** Calculates average Grueneisen parameter (GP) and isochoric specific heat (CV)
   *  This function requires to input the volume at which those parameters are
   *  calculated.
   *  This function is used when Grueneisen parameter is calculated using fit to
   *  w(V) relation obtained from EOS calculations.
   */
  void QHAN::calcCV_GP_fit(double T, double V, double &CV, double &GP)
  {
    CV = 0; GP = 0;
    if (!(T>0)) return;

    xvector<double> xomega;

    uint NIrQpoints = omegaV_mesh_EOS.size();
    int NQpoints = 0; // the total number of q-points

    double Cvi  = 0;  // mode-dependent specific heat at V=const
    double w    = 0;  // frequency for a given volume
    double expx = 0;  // temperature-dependent exponential factor
    double gamma = 0;

    for (uint q=0; q<NIrQpoints; q++){
      for (int branch=0; branch<Nbranches; branch++){
          xomega = aurostd::vector2xvector(omegaV_mesh_EOS[q][branch]);

          gamma = calcGamma(V, xomega, w);
          w *= THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]
          if (w > _mm_epsilon){
            expx = exp(w/KBOLTZEV/T);

            Cvi = pow(w,2)*expx/pow(expx-1.0,2) * qpWeights[q];

            GP += gamma * Cvi;
            CV += Cvi;
          }
        }
        NQpoints += qpWeights[q];
      }

    GP /= CV;
    CV /= NQpoints; CV /= NatomsOrigCell;
    CV *= pow(1/KBOLTZEV/T,2);
  }

  /** Calculates volumetric thermal expansion coefficient (beta).
   *
   * Central finite differences are used to calculate beta = 1/V dV/dT.
   * Values of V(T+dT) and V(T-dT) are taken from the EOS fit to the corresponding
   * free energies.
   *
   * For T=0K function returns 0, negative temperatures are ignored.
   *
   * @param T is the temperature at which the calculation is done.
   * @param eos_method defines which model is used for EOS fit
   * @return volumetric thermal expansion coefficient.
   */
  double QHAN::ThermalExpansion(double T, int eos_method)
  {
    if (!(T>0)) return 0;

    double dT = DCOEFF*T;
    return 0.5*(getEqVolumeT(T+dT,eos_method)-getEqVolumeT(T-dT,eos_method))/dT/
           getEqVolumeT(T,eos_method);
  }

  /** Calculates isochoric specific heat as a temperature derivative of free
   * energy using central finite differences method.
   *
   * @param eos_method defines which model is used for EOS fit
   */
  double QHAN::IsochoricSpecificHeat(double T, double V, int eos_method)
  {
    double dT = DCOEFF*T;
    double CV = 0;
    if (T>0){
      CV = -(FreeEnergyFit(T+dT,V,eos_method)-2*FreeEnergyFit(T,V,eos_method)
             +FreeEnergyFit(T-dT,V,eos_method));
      CV *= T/pow(dT,2);
    }

    return CV;
  }

  /** Calculates entropy as a temperature derivative of free energy using
   * central finite differences method.
   *
   * @param eos_method defines which model is used for EOS fit
   */
  double QHAN::Entropy(double T, double V, int eos_method)
  {
    double dT = DCOEFF*T;
    return -0.5*(FreeEnergyFit(T+dT,V,eos_method)-FreeEnergyFit(T-dT,V,eos_method))/dT;
  }

  /** Calculates DOS value at Fermi level using linear tetrahedron method. 
   *  For the details check: https://doi.org/10.1103/PhysRevB.49.16223
   */
  xvector<double> QHAN::DOSatEf()
  {
    xvector<double> ET(4); // energy eigenvalues at corners of tetrahedra
    double VOLWGT = 0; // weighted volume of tetrahedron
    double DEf; // DOS at Fermi level
    double Ef;  // Fermi energy
    double E21, E31, E41, E32, E42, E43; // Eij = ET[i] - ET[j]

    DOS_Ef =  xvector<double>(N_EOSvolumes);

    for (int id=0; id<N_EOSvolumes; id++){
      DEf = 0;
      for (uint band=0; band<static_eigvals[id].number_bands; band++){
        for (uint i=0; i<static_ibzkpts[id].ntetrahedra; i++){
          VOLWGT = static_ibzkpts[id].wtetrahedra * static_ibzkpts[id].vtetrahedra[i][1];
       
          for (uint s=0; s<static_eigvals[id].spin+1; s++){
            for (int j=1; j<=4; j++)
              ET[j] = static_eigvals[id].venergy[static_ibzkpts[id].vtetrahedra[i][j+1]-1][band][s];

            ET = aurostd::sort(ET);
            Ef = Efermi_V[id];

            if (Ef < ET[1] && ET[4] < Ef) DEf += 0;

            if (ET[1] < Ef && Ef <= ET[2]){
              E21 = ET[2] - ET[1];
              E31 = ET[3] - ET[1];
              E41 = ET[4] - ET[1];

              DEf += VOLWGT*3*pow(Ef-ET[1],2)/(E21*E31*E41);
            }
            else if (ET[2] < Efermi_V[id] && Efermi_V[id] <= ET[3]){
              E21 = ET[2] - ET[1];
              E31 = ET[3] - ET[1];
              E41 = ET[4] - ET[1];

              E32 = ET[3] - ET[2];
              E42 = ET[4] - ET[2];

              DEf += VOLWGT/(E31*E41)*(3*E21+6*(Ef-ET[2])-3*(E31+E42)*pow(Ef-ET[2],2)/(E32*E42));
            }
            else if (ET[3] < Ef && Ef <= ET[4]){
              E41 = ET[4] - ET[1];
              E42 = ET[4] - ET[2];
              E43 = ET[4] - ET[3];

              DEf += 3*VOLWGT*pow(ET[4]-Ef,2)/(E41*E42*E43);
            }
          }
        }
      }
      DOS_Ef[id+1] = (2-static_eigvals[id].spin)*DEf; // factor 2 if non-magnetic and 1 otherwise
    }
    return DOS_Ef;
  }

  /** Calculates electronic free energy using Sommerfeld expansion. */
  xvector<double> QHAN::electronicFreeEnergySommerfeld(double T)
  {
    xvector<double> F_Som = DOS_Ef;
    for (int i=F_Som.lrows; i<=F_Som.urows; i++){
      F_Som[i] *= -pow(M_PI*KBOLTZEV*T,2)/6.0;
    }
    return F_Som;
  }

  // QHA3P-related functions

  /** Returns frequency extrapolated by Taylor expansion
   */
  double QHAN::extrapolateFrequency(double V, xvector<double> &xomega)
  {
    double dwdV = (xomega[1]-xomega[3])/(GPvolumes[0]-GPvolumes[2]);
    double d2wdV2 = (xomega[1]+xomega[3]-2.0*xomega[2])/
      pow(0.5*(GPvolumes[0]-GPvolumes[2]),2);

    double dV = (V-GPvolumes[1]);
    return xomega[2] + dwdV*dV + 0.5*d2wdV2*pow(dV,2);
  }

  /** Returns Grueneisen parameter extrapolated by Taylor expansion
   */
  double QHAN::extrapolateGamma(double V, xvector<double> &xomega)
  {
    double gamma = 0;
    double V0 = GPvolumes[1];
    double w  = extrapolateFrequency(V, xomega);

    double dwdV = 0, d2wdV2 = 0;

    if (w>1e-6){
      dwdV = (xomega[1]-xomega[3])/(GPvolumes[0]-GPvolumes[2]);
      d2wdV2 = (xomega[1]+xomega[3]-2.0*xomega[2])/
        pow(0.5*(GPvolumes[0]-GPvolumes[2]),2);

      gamma = -V/w*(dwdV + d2wdV2 * (V-V0));
    }

    return gamma;
  }

  /* Free energy (total energy + vibrational contribution) obtained using Taylor 
   * expansion of frequencies.
   * This function is related to QHA3P and SCQHA methods.
   */
  double QHAN::FreeEnergyTE(double T, int Vid)
  {
    double w = 0.0; // extrapolated frequency at V_id volume
    double F = 0.0; // Free energy
    double fi = 0.0;
    int NQpoints = 0;
    xvector<double> xomega;
    for (uint q=0; q<omegaV_mesh.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);

        w = extrapolateFrequency(EOSvolumes[Vid], xomega) * THz2Hz *
            PLANCKSCONSTANTEV_h;
        fi = 0.5*w;

        if (w>_mm_epsilon && T>_mm_epsilon) fi += KBOLTZEV*T*log(1-exp(-w/KBOLTZEV/T));

        fi *= qpWeights[q];
        F += fi;
      }
      NQpoints += qpWeights[q];
    }
 
    F /= NQpoints;
    F /= NatomsOrigCell;

    return F + E0_V[Vid];
  }

  /* Vibrational internal energy obtained using Taylor expansion of frequencies.
   * This function is related to QHA3P and SCQHA methods. */
  double QHAN::InternalEnergyTE(double T, double V)
  {
    double U = 0.0, ui = 0.0,  w = 0.0; 
    int NQpoints = 0;
    xvector<double> xomega;
    for (uint q=0; q<omegaV_mesh.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);

        w = extrapolateFrequency(V, xomega) * THz2Hz * PLANCKSCONSTANTEV_h;
        ui = 0.5*w;

        if (w>_mm_epsilon && T>_mm_epsilon) ui += w/(exp(w/KBOLTZEV/T)-1.0);

        ui *= qpWeights[q];
        U += ui;
      }
      NQpoints += qpWeights[q];
    }
 
    U /= NQpoints;
    U /= NatomsOrigCell;

    return U;
  }

  // SCQHA-related functions
  // Implementation is based on http://dx.doi.org/10.1103/PhysRevMaterials.3.073801
  // and https://doi.org/10.1016/j.commatsci.2016.04.012

  /** Phononic pressure multiplied by volume.
   * Check http://dx.doi.org/10.1103/PhysRevMaterials.3.073801
   * for more details
   */
  double QHAN::VPgamma(double T, double V)
  {
    double VPgamma = 0.0, ui = 0.0,  w = 0.0; 
    int NQpoints = 0;
    xvector<double> xomega;
    for (uint q=0; q<omegaV_mesh.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);

        w = extrapolateFrequency(V, xomega) * THz2Hz * PLANCKSCONSTANTEV_h;
        ui = 0.5*w;

        if (w>_mm_epsilon && T>_mm_epsilon) ui += w/(exp(w/KBOLTZEV/T)-1.0);

        ui *= qpWeights[q];
        VPgamma += ui * extrapolateGamma(V, xomega);
      }
      NQpoints += qpWeights[q];
    }

    VPgamma /= NQpoints;
    VPgamma /= NatomsOrigCell;

    return VPgamma;
  }

  /** Performs SCQHA calculation.
   *  Here there are two implementations:
   *  1) perform self-consistent loop for initial nonzero temperature and extrapolate
   *  the volume at next temperature step using V(T+dT) ~ (1+beta dT)*V.
   *  Expect it to be inaccurate at high temperatures.
   *  2) perform self-consistent loop for each temperature
   */
  void QHAN::run_scqha(int method, bool all_iterations_self_consistent)
  {
    int max_scqha_iteration = 10000;
    double dV = 1e-3;
    double Vtol = 1e-5;

    string function = "QHAN::run_scqha()";
    pflow::logger(QHA_ARUN_MODE, function, "", currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);

    // get equilibrium volume from the fit to EOS based on energies from static
    // calculations
    eos.method = method;
    eos.E = aurostd::vector2xvector<double>(E0_V);
    eos.fit();
    double Pe, VPg; // electronic pressure and volume multiplied by phononic pressure
    double V = 1.1*eos.Veq; // to avoid division by zero in self-consistent loop
    // initial volume is taken to be 10% bigger
    double Vnew = 0;

    // self-consistent loop for 0K temperature
    double V0K = 1.1*eos.Veq;
    for (int i=0; i<max_scqha_iteration; i++){
      Pe   = dEOS_poly(V, eos.p);
      VPg  = VPgamma(0.0,V);
      Vnew = VPg/Pe;
      if (abs(V0K - Vnew)/V > Vtol) V0K += (Vnew - V0K) * dV; else break;
    }

    // output file name depends on the used EOS fit method
    ofstream file;
    string filename;
    string sc = all_iterations_self_consistent ? "sc." : "";
    switch(method){
      case(EOS_POLYNOMIAL):
        filename = "aflow.scqha."+sc+"polynomial.thermo.out";
        break;
      case(EOS_BIRCH_MURNAGHAN):
        filename = "aflow.scqha."+sc+"birch-murnaghan.thermo.out";
        break;
      case(EOS_MURNAGHAN):
      default:
        filename = "aflow.scqha."+sc+"murnaghan.thermo.out";
        break;
    }

    file.open(filename.c_str());
    file.precision(10);

    // print header
    file << setw(5)  << "# T[K]"               << setw(SW) << ' ' <<
            setw(TW) << "Veq[ev/atom]"         << setw(SW) << ' ' <<
            setw(TW) << "F(Veq)[eV/atom]"      << setw(SW) << ' ' <<
            setw(TW) << "B[GPa]"               << setw(SW) << ' ' <<
            setw(TW) << "beta[10^-6/K]"        << setw(SW) << ' ' <<
            setw(TW) << "Cv[kB/atoms]"         << setw(SW) << ' ' <<
            setw(TW) << "Cp[kB/atoms]"         << setw(SW) << ' ' <<
            setw(TW) << "gamma"                << setw(SW) << ' ' <<
            '\n';

    double T = Temperatures[0];

    // self-consistent loop for initial temperature (defined by user)
    int iter=1;
    while (iter <= max_scqha_iteration){
      Pe   = dEOS_poly(V, eos.p);
      VPg  = VPgamma(T, V);
      Vnew = VPg/Pe;
      if (abs(V - Vnew)/V > Vtol) V += (Vnew - V) * dV; else break;
      iter++;
    }
    if (iter == max_scqha_iteration){
      string msg="Maximum number of iterations in self consistent loop is reached";
      msg += " at T="+apl::stringify(T)+"K.";
      pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);
    }

    double dT = (Temperatures[Ntemperatures-1]-Temperatures[0])/(Ntemperatures-1);

    // a set of variables predefined here
    double Cvi     = 0.0;  // mode-dependent specific heat at V=const
    double w       = 0.0;  // frequency for a given volume
    double w0K     = 0.0;  // frequency for a volume at T=0K
    double expx    = 0.0;  // temperature-dependent exponential factor
    double expx0   = 0.0;  // temperature-dependent exponential factor
    double ui      = 0.0;  // mode-dependent internal energy
    double Bgamma  = 0.0;  // contribution to bulk modulus of 2nd order gamma component
    double Bdgamma = 0.0;  // check https://doi.org/10.1016/j.commatsci.2016.04.012
    double Belec   = 0.0;  // "electronic" bulk modulus
    double B       = 0.0;  // total bulk modulus
    double Pgamma  = 0.0;  // "phononic" pressure
    double gamma   = 0.0;  // Grueneisen parameter
    double fi      = 0.0;  // mode-dependet free energy
    double Feq     = 0.0;  // total free energy for equilibrium volume at given T
    double CP      = 0.0;  // isobaric specific heat
    double CV      = 0.0;  // isochoric specific heat
    double GP      = 0.0;  // average Grueneisen parameter
    double beta    = 0.0;  // coefficient of thermal expansion
    double d2wdV2  = 0.0;  // second derivative of frequency w.r.t volume

    uint NIrQpoints = omegaV_mesh.size(); // number of irreducible q-points
    int NQpoints = 0; // the total number of q-points (to be determined in loop)

    for (int Tid=0; Tid<Ntemperatures; Tid++){
      T = Temperatures[Tid];

      // calculate next equilibrium volume using self-consistent loop if beta=0 or
      // if user wants so
      if (all_iterations_self_consistent || !(abs(beta)>0)){
        iter=1;
        while (iter <= max_scqha_iteration){
          Pe   = dEOS_poly(V, eos.p);
          VPg  = VPgamma(T, V);
          Vnew = VPg/Pe;
          if (abs(V - Vnew)/V > Vtol) V += (Vnew - V) * dV; else break;
          iter++;
        }
        if (iter == max_scqha_iteration){
          string msg="Maximum number of iterations in self consistent loop is reached";
          msg += " at T="+apl::stringify(T)+"K";
          pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, 
              *oss, _LOGGER_MESSAGE_);
        }
      }
      else{
        V *= (1 + beta*dT);
      }

      // calculate thermodynamic properties
      Cvi = 0.0; CP = 0; CV = 0.0; 
      GP  = 0.0; gamma   = 0.0;
      w       = 0.0; w0K     = 0.0; expx    = 0.0; expx0   = 0.0;
      ui      = 0.0;
      Bgamma  = 0.0; Bdgamma = 0.0; Belec   = 0.0; Pgamma  = 0.0; B = 0.0;
      fi = 0.0; Feq = 0.0;

      for (uint q=0; q<NIrQpoints; q++){
        for (uint branch=0; branch<omegaV_mesh[q].size(); branch++){
          xvector<double> xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);
          w = extrapolateFrequency(V, xomega);
          w *= THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]

          w0K = extrapolateFrequency(V0K, xomega);
          w0K *= THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]

          ui = 0.5*w;
          fi = 0.5*w;
          if (w > 1e-12 && T > 0){
            expx = exp(w/KBOLTZEV/T);

            // use volume at T=0K for calculation of isochoric specific heat
            expx0 = exp(w0K/KBOLTZEV/T);
            Cvi = pow(w,2)*expx0/pow(expx0-1.0,2) * qpWeights[q];

            ui  += w/(expx - 1.0);
            ui  *= qpWeights[q];

            fi += KBOLTZEV*T*log(1-exp(-w/KBOLTZEV/T));
            fi *= qpWeights[q];

            d2wdV2 = (xomega[1]+xomega[3]-2.0*xomega[2])/
                      pow(0.5*(GPvolumes[0]-GPvolumes[2]),2);
            d2wdV2 *= THz2Hz*PLANCKSCONSTANTEV_h;

            gamma = extrapolateGamma(V, xomega);

            Bgamma  += (ui - T*Cvi*pow(1/KBOLTZEV/T,2)*KBOLTZEV)*pow(gamma, 2);
            Bdgamma -= ui*((1+gamma)*gamma - pow(V,2)/w*d2wdV2);

            GP += extrapolateGamma(V, xomega) * Cvi;
            CV += Cvi;
            Feq  += fi;
          }
        }
        NQpoints += qpWeights[q];
      }

      Feq /= NQpoints;
      Feq /= NatomsOrigCell;
      Feq += eos.eval(V);

      Bgamma /= NQpoints; Bdgamma /= NQpoints;
      Bgamma /= NatomsOrigCell; Bdgamma /= NatomsOrigCell;
      Bgamma /= V; Bdgamma /= V;

      Belec  = BulkModulus(V, eos.p);
      Pgamma = VPgamma(T, V)/V;

      B = Belec + (Bgamma + Bdgamma + Pgamma)*eV2GPa;

      GP /= CV;
      CV /= NQpoints; CV /= Nbranches;
      CV *= 3*pow(1/KBOLTZEV/T,2); // [kB/atom]

      beta = KBOLTZEV*CV*GP/V/(B/eV2GPa); // [K^-1]

      CP = CV + V*T*B*pow(beta,2)/eV2GPa/KBOLTZEV; // [kB/atom]
      file << setw(5)  << T                   << setw(SW) << ' ' <<
              setw(TW) << V                   << setw(SW) << ' ' <<
              setw(TW) << Feq                 << setw(SW) << ' ' <<
              setw(TW) << B                   << setw(SW) << ' ' <<
              setw(TW) << beta * 1e6          << setw(SW) << ' ' << //[10^-6/K]
              setw(TW) << CV                  << setw(SW) << ' ' <<
              setw(TW) << CP                  << setw(SW) << ' ' <<
              setw(TW) << GP                  << setw(SW) << ' ' <<
              '\n';
    }
    file.close();
  }

  ////////////////////////////////////////////////////////////////////////////////
  // output-related functions

  /** Writes temperature-dependent properties to file.
   *
   * List of properties: temperature, volume, free energy, bulk modulus,
   * thermal expansion coefficient, specific heat at const V, specific heat at
   * const P, average Grueneisen parameter calculated from mode-dependent
   * Grueneisen parameters and average Grueneisen calculated from thermal
   * expansion coefficient.
   *
   */
  void QHAN::writeThermalProperties(int eos_method, int qha_method)
  {
    eos.method = eos_method;

    // type of qha calculation
    string qha;
    switch(qha_method){
      case (QHA3P_CALC):
        qha = "qha3p";
        break;
      case(QHA_CALC):
      default:
        qha = "qha";
        break;
    }

    // the name of the output file depends on the used EOS fit method and on type of
    // QHA calculation
    ofstream file;
    string filename;
    switch(eos_method){
      case(EOS_POLYNOMIAL):
        filename = "aflow." + qha + ".polynomial.thermo.out";
        break;
      case(EOS_BIRCH_MURNAGHAN):
        filename = "aflow." + qha + ".birch-murnaghan.thermo.out";
        break;
      case(EOS_MURNAGHAN):
      default:
        filename = "aflow." + qha + ".murnaghan.thermo.out";
        break;
    }

    string function = "QHAN::writeThermalProperties()";
    string msg = "Writing T-dependent properties to "+filename;
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);

    file.open(filename.c_str());
    file.precision(10);

    // write header
    file << setw(5)  << "# T[K]"               << setw(SW) << ' ' <<
            setw(TW) << "Veq[ev/atom]"         << setw(SW) << ' ' <<
            setw(TW) << "F(V0)[eV/atom]"       << setw(SW) << ' ' <<
            setw(TW) << "B[GPa]"               << setw(SW) << ' ' <<
            setw(TW) << "beta[10^-6/K]"        << setw(SW) << ' ' <<
            setw(TW) << "Cv[kB/atoms]"         << setw(SW) << ' ' <<
            setw(TW) << "Cp[kB/atoms]"         << setw(SW) << ' ' <<
            setw(TW) << "gamma(beta,B,Cv)"     << setw(SW) << ' ' <<
            setw(TW) << "Bprime"               << setw(SW) << ' ' <<
            setw(TW) << "beta_mesh[10^-6/K]"   << setw(SW) << ' ' <<
            setw(TW) << "Cv_mesh[kB/atoms]"    << setw(SW) << ' ' <<
            setw(TW) << "Cp_mesh[kB/atoms]"    << setw(SW) << ' ' <<
            setw(TW) << "gamma_mesh"           << setw(SW) << ' ' <<
            '\n';

    xvector<double> F(N_EOSvolumes); // free energy
    xvector<double> xvolumes = aurostd::vector2xvector(EOSvolumes);

    switch(qha_method){
      case(QHA3P_CALC):
        for (int Vid=0; Vid<N_EOSvolumes; Vid++) eos.E[Vid+1] = FreeEnergyTE(0, Vid);
        break;
      case(QHA_CALC):
      default:
        for (int Vid=0; Vid<N_EOSvolumes; Vid++) eos.E[Vid+1] = FreeEnergy(0, Vid);
        break;
    }

    eos.fit();
    double V0K = eos.Veq; // equilibrium volume at 0K

    double T, Veq, Feq, B, Bp, beta, CV, CP, GP;
    double CV_mesh, GP_mesh, CP_mesh, beta_mesh; // these properties are calculated 
    // by weighted sum over q-points mesh
    for (int Tid=0; Tid<Ntemperatures; Tid++){
      T = Temperatures[Tid];

      switch(qha_method){
        case(QHA3P_CALC):
          for (int Vid=0; Vid<N_EOSvolumes; Vid++) eos.E[Vid+1] = FreeEnergyTE(T, Vid);
          break;
        case(QHA_CALC):
        default:
          for (int Vid=0; Vid<N_EOSvolumes; Vid++) eos.E[Vid+1] = FreeEnergy(T, Vid);
          break;
      }

      if (includeElectronicContribution) eos.E += electronicFreeEnergySommerfeld(T);

      // stop if energy minimum is no longer within a given set of volumes
      if (!isMinimumWithinBounds(eos.E)){
        msg = "Stopping at T=" + apl::stringify(T) + " [K]";
        msg+= "since there is no free energy minimum within a given volume range.";
        pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);
        break;
      }

      eos.fit();
      // note that state of EOSfit is changed by ThermalExpansion and/or
      // IsochoricSpecificHeat functions, so save Veq, Feq and B for future use
      Veq = eos.Veq;
      Feq = eos.Eeq;
      B   = eos.B;  // [GPa]
      Bp  = eos.Bp;
      beta = ThermalExpansion(T, qha_method); // [K^-1]
      CV   = IsochoricSpecificHeat(T, V0K, qha_method)/KBOLTZEV; // [kB/atom]
      CP   = CV + Veq*T*B*pow(beta,2)/eV2GPa/KBOLTZEV; // [kB/atom]
      GP   = (beta/CV)*B*Veq/eV2GPa/KBOLTZEV;

      calcCV_GP_fit(T, Veq, CV_mesh, GP_mesh);
      beta_mesh = KBOLTZEV*CV_mesh*GP_mesh/Veq/(B/eV2GPa); // [K^-1]
      CP_mesh   = CV_mesh + Veq*T*B*pow(beta_mesh,2)/eV2GPa/KBOLTZEV; // [kB/atom]

      // write values to file
      file << setw(5)  << T                   << setw(SW) << ' ' <<
              setw(TW) << Veq                 << setw(SW) << ' ' <<
              setw(TW) << Feq                 << setw(SW) << ' ' <<
              setw(TW) << B                   << setw(SW) << ' ' <<
              setw(TW) << beta * 1e6          << setw(SW) << ' ' << //[10^-6/K]
              setw(TW) << CV                  << setw(SW) << ' ' <<
              setw(TW) << CP                  << setw(SW) << ' ' <<
              setw(TW) << GP                  << setw(SW) << ' ' <<
              setw(TW) << Bp                  << setw(SW) << ' ' <<
              setw(TW) << beta_mesh * 1e6     << setw(SW) << ' ' << //[10^-6/K]
              setw(TW) << CV_mesh             << setw(SW) << ' ' <<
              setw(TW) << CP_mesh             << setw(SW) << ' ' <<
              setw(TW) << GP_mesh             << setw(SW) << ' ' <<
              '\n';
    }

    file.close();
  }

  /** Writes F(V,T) data to aflow.qha.FVT.out file
   */
  void QHAN::writeFVT()
  {
    string function = "QHA::writeFVT()";
    string msg = "Writing F(V,T) relations to file.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);

    ofstream file;
    file.open("aflow.qha.FVT.out");
    file.precision(10);

    xvector<double> Felec(N_EOSvolumes);
    double T;

    for (int Tid = 0; Tid < Ntemperatures; Tid++){
      for (int Vid = 0; Vid < N_EOSvolumes; Vid++){
        T = Temperatures[Tid];
        if (includeElectronicContribution) Felec = electronicFreeEnergySommerfeld(T);
        file << setw(TW) << EOSvolumes[Vid]    << setw(SW) << ' '
             << setw(TW) << FreeEnergy(T, Vid) << setw(SW) << ' '
             << setw(TW) << Felec[Vid+1]       << setw(SW) << ' '
             << setw(TW) << E0_V[Vid] <<
	     '\n';
      }
      file <<  "\n\n";
    }
    file.close();
  }

  /** Writes mode-dependent Grueneisen parameter along a given path in k-space.
   */
  void QHAN::writeGPpath(double V, const string &directory)
  {
    string function = "QHAN::writeGPpath()";
    string msg = "Calculate and save Grueneisen parameters along a path.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);
    // we will save bands-projected Grueneisen parameter in xEIGENVAL
    xEIGENVAL GPpath(gp_ph_dispersions.front());
    GPpath.Vol = V;
    GPpath.temperature = 0;

    xstructure struc = origStructure;
    struc.SetVolume(V);
    xvector<double> lattice(3);
    lattice[1] = struc.a * 1E-10;
    lattice[2] = struc.b * 1E-10;
    lattice[3] = struc.c * 1E-10;
    GPpath.lattice = lattice;
    GPpath.carstring = "QHA";

    xvector<double> xomega(N_GPvolumes);
    //venergy.at(kpoint number).at(band number).at(spin number)
    for (uint q=0; q<GPpath.venergy.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        for (int Vid=0; Vid<N_GPvolumes; Vid++){
          xomega[Vid+1] = gp_ph_dispersions[Vid].venergy.at(q).at(branch).at(0);
        }

        //GPpath.venergy.at(q).at(branch).at(0) = calcGamma(V, xomega);
        GPpath.venergy.at(q).at(branch).at(0) = calcGammaFD(xomega);
      }
    }

    stringstream eigenval;
    eigenval << GPpath;

    string filename = directory + '/' + DEFAULT_QHA_GP_PATH_FILE;
    aurostd::stringstream2file(eigenval, filename);
    if (!aurostd::FileExist(filename)){
      msg = "Cannot open "+filename+" file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
    }
  }


  /** Writes average Grueneisen parameter, which is calculated using finite
   *  difference method and isochoric specific heat
   *  Output file: "aflow.qha.gp.avg.out"
   */
  void QHAN::writeAverageGP_FD()
  {
    string function = "QHAN::writeAverageGP_FD";
    string msg = "Writing T-dependence of average Grueneisen parameter.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);

    ofstream file;

    file.open("aflow.qha.gp.avg.out");
    file.precision(10);

    file << setw(5)  << "# T[K]"           << setw(SW) << ' ' <<
            setw(TW) << "gamma"            << setw(SW) << ' ' <<
            setw(TW) << "CV"               << setw(SW) << ' ' <<
            '\n';

    double T, CV, GP;
    for (int Tid=0; Tid<Ntemperatures; Tid++){
      T = Temperatures[Tid];

      calcCV_GP(T, CV, GP);

      file << setw(5)  << T                   << setw(SW) << ' ' <<
              setw(TW) << GP                  << setw(SW) << ' ' <<
              setw(TW) << CV                  << setw(SW) << ' ' <<
              '\n';
    }

    file.close();
  }

  /** Writes mode-dependent Grueneisen parameter calculated at each q-point in IBZ.
   */
  void QHAN::writeGPmeshFD()
  {
    string function = "QHAN::writeGPmeshFD()";
    string msg = "Writing Grueneisen parameter calculated on mesh of q-points.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *messageFile, *oss,
        _LOGGER_MESSAGE_);

    ofstream file;
    file.open((DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_GP_MESH_FILE).c_str());
    file.precision(6);

    // write header
    file << "#" << setw(TW) << "Grueneisen" << setw(SW) << ' ' <<
                   setw(TW) << "Freq(THz)"  << setw(SW) << ' ' <<
                   setw(TW) << "phonon_mode\n";

    xvector<double> xomega;

    uint NIrQpoints = omegaV_mesh.size();

    for (uint q=0; q<NIrQpoints; q++){
      file << "#k-point: " << qPoints[q] << " weight: " << qpWeights[q] << '\n';
      for (uint branch=0; branch<omegaV_mesh[q].size(); branch++){
          xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);

          file << setw(TW) << calcGammaFD(xomega) << setw(SW) << " " <<
                  setw(TW) << xomega(2)           << setw(SW) << " " << 
                  setw(TW) << branch << '\n';
        }
      file << '\n';
      }
    file.close();
  }

  /** Writes volume-dependent phonon frequencies obtained from series of EOS APL
   * calculations
   */
  void QHAN::writeFrequencies()
  {
    ofstream file;
    file.open("aflow.qha.frequencies.out");
    for (uint q=0; q<omegaV_mesh_EOS.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        for (int Vid=0; Vid<N_EOSvolumes; Vid++){
          file << EOSvolumes[Vid] << " " << omegaV_mesh_EOS[q][branch][Vid] << '\n';
        }
        file << "\n\n";
      }
    }
    file.close();
  }
}
