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

// Global mutex that prevents two xThread instances from checking the number
// of available CPUs at the same time.
static std::mutex xthread_cpu_check;

namespace xthread {

  /// @brief Constructur for xThread
  ///
  /// @param nmax Maximum number of CPUs used by xThread (default: 0 for all available CPUs)
  /// @param nmin Mininum number of CPUs required to spawn thread workers default: 0 for nmin = nmax)
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

  /// @brief Sets the minimum and maximum number of CPUs used for threading
  ///
  /// @param nmax Maximum number of CPUs used by xThread (default: 0 for all available CPUs)
  /// @param nmin Mininum number of CPUs required to spawn thread workers default: 0 for nmin = nmax)
  void xThread::setCPUs(int nmax, int nmin) {
    if (nmax < nmin) std::swap(nmax, nmin);
    if (nmax <= 0) nmax = init::GetCPUCores();
    ncpus_max = nmax;
    if (nmin <= 0) ncpus_min = nmax;
  }

  /// @brief Tells xThread that it has to update a progress bar when running
  ///
  /// @param oss The ostream for the progress bar
  void xThread::setProgressBar(ostream& oss) {
    progress_bar = &oss;
    progress_bar_set = true;
  }

  void xThread::unsetProgressBar() {
    progress_bar_set = false;
  }

  /// @brief Worker called by the threads
  ///
  /// @param task_index The index of the next bin to run, shared across all threads
  /// @param ntasks The number of task over which the function is parallelized
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  template <typename F, typename... A>
  void xThread::spawnWorker(int& task_index, int ntasks,
                            F& func, A&... args) {
    // Initial distribution of the tasks over all threads
    int icurr = -1;
    if (task_index < ntasks) {
      mtx.lock();
      icurr = task_index++;
      mtx.unlock();
    } else {
      return;
    }

    while (icurr < ntasks) {
      func(icurr, args...); // Call function
      // Update the task index to the next available task
      mtx.lock();
      icurr = task_index++;
      if (progress_bar_set && (task_index <= ntasks)) pflow::updateProgressBar(task_index, ntasks, *progress_bar);
      mtx.unlock();
    }
  }

  /// @brief Executes a function over at least ncpus_min and at most ncus_max threads
  ///
  /// @param ntasks The number of bins over which the function is parallelized
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  template <typename F, typename... A>
  void xThread::run(int ntasks, F& func, A&... args) {
    // First check if enough threads are available using XHOST.CPU_active,
    // which every multi-threaded call should update.
    // This prevents threaded functions that spawn other multi-threaded
    // processes from allocating more threads than the machine can afford.
    // This only works within a single AFLOW run
    uint sleep_second = 10;
    int ncpus_max_available = init::GetCPUCores();
    int ncpus_available = ncpus_max_available - XHOST.CPU_active;
    int ncpus = 0;
    while (ncpus_available < ncpus_min) {
      xthread_cpu_check.lock();
      ncpus_available = ncpus_max_available - XHOST.CPU_active;
      if (ncpus_available >= ncpus_min) {
        ncpus = (ncpus_available > ncpus_max)?ncpus_max:ncpus_available;
        // "reserve" threads globally
        XHOST.CPU_active += ncpus;
      }
      xthread_cpu_check.unlock();
      aurostd::Sleep(sleep_second);
    }

    // Initialize progress bar
    if (progress_bar_set) pflow::updateProgressBar(0, ntasks, *progress_bar);

    int task_index = 0;
    if (ncpus > 1) {
      vector<std::thread*> threads;
      for (int i = 0; i < ncpus; i++) {
        threads.push_back(new std::thread(&xThread::spawnWorker<F, A...>, this,
                                          std::ref(task_index), ntasks,
                                          std::ref(func), std::ref(args)...)
        );
      }
      for (std::thread* t : threads) {
        t->join();
        delete t;
      }
    } else {
      // No need for thread overhead when ncpus == 1
      while (task_index < ntasks) {
        func(task_index, args...);
        if (progress_bar_set) pflow::updateProgressBar(++task_index, ntasks, *progress_bar);
      }
    }

    // "free" threads globally
    xthread_cpu_check.lock();
    XHOST.CPU_active -= ncpus;
    xthread_cpu_check.unlock();
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
    std::function<void(int, const vector<uint>&, const aurostd::xoption&, vector<apl::DOSCalculator>&, vector<xDOSCAR>&)> fnpoccapl;
    xt.run(i, fnpoccapl, vuint, xopt, vdoscalc, vxdos);
  }

}

#endif
