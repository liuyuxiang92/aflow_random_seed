// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************
// Marco Esters
// Contains thread manager class xThread

#ifdef AFLOW_MULTITHREADS_ENABLE

#include "aflow.h"
#include "APL/aflow_apl.h"

namespace xthread {

  xThread::xThread(int nmax, int nmin) {
    free();
    setCPUs(nmax, nmin);
  }

  xThread::xThread(const xThread& xt) {
    copy(xt);
  }

  const xThread& xThread::operator=(const xThread& xt) {
    copy(xt);
    return *this;
  }

  void xThread::copy(const xThread& xt) {
    if (this == &xt) return;
    // std::mutex should not be copied because
    // it needs to stay immutable
    ncpus_max = xt.ncpus_max;
    ncpus_max = xt.ncpus_min;
    progress_bar = xt.progress_bar;
    progress_bar_set = progress_bar_set;
  }

  xThread::~xThread() {
    free();
  }

  void xThread::free() {
    ncpus_max = 0;
    ncpus_min = 0;
    progress_bar_set = false;
  }

  void xThread::setCPUs(int nmax, int nmin) {
    if (nmax < nmin) std::swap(nmax, nmin);
    ncpus_max = nmax;
    if (nmin == 0) ncpus_min = nmax;
  }

  void xThread::setProgressBar(ostream& oss) {
    progress_bar = &oss;
    progress_bar_set = true;
  }

  void xThread::unsetProgressBar() {
    progress_bar_set = false;
  }

  template <typename F, typename... A>
  void xThread::spawnWorker(int& task_index, int nbins,
                            F& func, A&... args) {
    int icurr = -1;
    if (task_index < nbins) {
      std::unique_lock<std::mutex> lk(mtx);
      icurr = task_index++;
    } else {
      return;
    }

    while (icurr < nbins) {
      func(icurr, args...);
      std::unique_lock<std::mutex> lk(mtx);
      icurr = task_index++;
      if (progress_bar_set && (task_index <= nbins)) pflow::updateProgressBar(task_index, nbins, *progress_bar);
    }
  }

  template <typename F, typename... A>
  void xThread::run(int nbins, F& func, A&... args) {
    uint sleep_second = 10;
    int ncpus_max_available = init::GetCPUCores();
    int ncpus_available = ncpus_max_available - XHOST.CPU_active;
    while (ncpus_available < ncpus_min) {
      ncpus_available = ncpus_max_available - XHOST.CPU_active;
      aurostd::Sleep(sleep_second);
    }
    int ncpus = (ncpus_available > ncpus_max)?ncpus_max:ncpus_available;
    XHOST.CPU_active += ncpus;

    int task_index = 0;
    if (progress_bar_set) pflow::updateProgressBar(0, nbins, *progress_bar);

    vector<std::thread*> threads;
    for (int i = 0; i < ncpus; i++) {
      threads.push_back(new std::thread(&xThread::spawnWorker<F, A...>, this,
                                        std::ref(task_index), nbins,
                                        std::ref(func), std::ref(args)...)
      );
    }

    for (std::thread* t : threads) {
      t->join();
      delete t;
    }
    XHOST.CPU_active -= ncpus;
  }

}

namespace xthread {

  void initializeXThread() {
    xThread xt;
    int i = 0;

    // Pure vector types
    vector<uint> vuint;

    vector<vector<double> > vvdbl;
    const vector<vector<double> > cvvdbl = vvdbl;
    vector<vector<vector<vector<double> > > > vvvvdbl;
    const vector<vector<vector<vector<double> > > > cvvvvdbl = vvvvdbl;

    vector<string> vs;
    const vector<string> cvs;
    vector<vector<string> > vvs;

    vector<vector<vector<xcomplex<double> > > > vvvxcdbl;


    // xvector types
    vector<vector<xvector<double> > > vvxvdbl;
    const vector<vector<xvector<double> > > cvvxvdbl = vvxvdbl;

    // xmatrix types
    vector<xmatrix<xcomplex<double> > > vxmxcdbl;

    // special types
    vector<aflowlib::DBStats> vdbs;
    aurostd::xoption xopt;
    vector<apl::DOSCalculator> vdoscalc;
    vector<xDOSCAR> vxdos;


    //Instantiate function with no args
    //apl::AtomicDisplacements::calculateEigenvectorsInThread
    //apl::DOSCalculator::calculateInOneThread
    //apl::PhononDispersionCalculator::calculateInOneThread
    //apl::TCONDCalculator::calculateTransitionProbabilitiesIsotope
    std::function<void(int)> fnnoargs;
    xt.run(i, fnnoargs);

    //AAPL
    //apl::TCONDCalculator::calculateTransitionProbabilitiesPhonon
    std::function<void(int, vector<vector<vector<vector<double> > > >&, const vector<vector<vector<xcomplex<double> > > >&)> fn_vvvvdbl_vvvxcdbl;
    xt.run(i, fn_vvvvdbl_vvvxcdbl, vvvvdbl, vvvxcdbl);
    //apl::TCONDCalculator::calculateAnharmonicRates
    using ft_vvdbl_vvdbl = std::function<void(int, const vector<vector<double> >&, vector<vector<double> >&)>;
    ft_vvdbl_vvdbl fn_vvdbl_vvdbl;
    xt.run(i, fn_vvdbl_vvdbl, cvvdbl, vvdbl);

    //apl::TCONDCalculator::calculateDelta
    std::function<void(int, const vector<vector<double> >&, const vector<vector<xvector<double> > >&, vector<vector<xvector<double> > >&)> fn_vvdbl_vvxvdbl_vvxvdbl;
    xt.run(i, fn_vvdbl_vvxvdbl_vvxvdbl, cvvdbl, vvxvdbl, vvxvdbl);

    //APL
    //apl::PhononCalculator::calculateGroupVelocitiesThread
    std::function<void(int, vector<vector<double> >&, vector<xmatrix<xcomplex<double> > >&, vector<vector<xvector<double> > >&)> fn_vvdbl_vxmxcdbl_vvxvdbl;
    xt.run(i, fn_vvdbl_vxmxcdbl_vvxvdbl, vvdbl, vxmxcdbl, vvxvdbl);

    //AflowDB
    //aflowlib::AflowDB::createTable
    std::function<void(int, const vector<string>&, const vector<string>&)> fn_vs_vs;
    xt.run(i, fn_vs_vs, vs, vs);
    //aflowlib::AflowDB::getColStats
    std::function<void(int, const vector<string>&, vector<aflowlib::DBStats>&)> fn_vs_vdbs;
    xt.run(i, fn_vs_vdbs, cvs, vdbs);

    //POCC
    //POccCalculator::calculatePhononDOSThread
    std::function<void(int, const vector<uint>&, const aurostd::xoption&, vector<apl::DOSCalculator>&, vector<xDOSCAR>&)> fn;
    xt.run(i, fn, vuint, xopt, vdoscalc, vxdos);

  }

}

#endif
