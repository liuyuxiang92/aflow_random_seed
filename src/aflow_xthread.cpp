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
    progress_bar_counter = xt.progress_bar_counter;
    progress_bar_set = progress_bar_set;
  }

  xThread::~xThread() {
    free();
  }

  void xThread::free() {
    ncpus_max = 0;
    ncpus_min = 0;
    progress_bar = nullptr;
    progress_bar_set = false;
    progress_bar_counter = 0;
  }

  /// @brief Sets the minimum and maximum number of CPUs used for threading
  ///
  /// @param nmax Maximum number of CPUs used by xThread (default: 0 for all available CPUs)
  /// @param nmin Mininum number of CPUs required to spawn thread workers default: 0 for nmin = nmax)
  void xThread::setCPUs(int nmax, int nmin) {
    if (nmax < nmin) std::swap(nmax, nmin);
    ncpus_max = (nmax > 0)?nmax:(init::GetCPUCores());
    ncpus_min = (nmin > 0)?nmin:nmax;
  }

  /// @brief Tells xThread that it has to update a progress bar when running
  ///
  /// @param oss The ostream for the progress bar
  void xThread::setProgressBar(ostream& oss) {
    progress_bar = &oss;
    progress_bar_set = true;
    progress_bar_counter = 0;
  }

  void xThread::unsetProgressBar() {
    progress_bar_set = false;
  }

  void xThread::initializeProgressBar(unsigned long long int ntasks) {
    progress_bar_counter = 0;
    pflow::updateProgressBar(0, ntasks, *progress_bar);
  }

  /// @brief Worker called by the threads
  ///
  /// @param it The iterator or task index
  /// @param end The end() iterator or the last task index
  /// @param ntasks Number of tasks (for the progress bar)
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  template <typename I, typename F, typename... A>
  void xThread::spawnWorker(I& it, I& end,
                            unsigned long long int ntasks,
                            F& func, A&... args) {
    I icurr = advance(it, end, ntasks);

    while (icurr != end) {
      func(icurr, args...); // Call function
      icurr = advance(it, end, ntasks, progress_bar_set);
    }
  }

  /// @brief Advances to the next task
  ///
  /// @param it The iterator or task index
  /// @param end The end() iterator or the last task index
  /// @param ntasks Number of tasks (for the progress bar)
  /// @param update_progress_bar Update progress bar if true (default: false)
  ///
  /// @return The next task index or iterator position
  template <typename I>
  I xThread::advance(I& it, I& end,
                     unsigned long long int ntasks,
                     bool update_progress_bar) {
    std::lock_guard<std::mutex> lk(mtx);
    if (update_progress_bar && (progress_bar_counter <= ntasks)) {
      pflow::updateProgressBar(++progress_bar_counter, ntasks, *progress_bar);
    }
    return (it == end)?end:(it++);
  }

  /// @brief Overload for running threads using indices to keep track of tasks
  ///
  /// @param ntasks The number of bins over which the function is parallelized
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  template <typename F, typename... A>
  void xThread::run(int ntasks, F& func, A&... args) {
    if (ntasks <= 0) return;
    int task_index = 0;
    run<int, F, A...>(task_index, ntasks, (unsigned long long int) ntasks, func, args...);
  }

  template <typename F, typename... A>
  void xThread::run(uint ntasks, F& func, A&... args) {
    if (ntasks == 0) return;
    uint task_index = 0;
    run<uint, F, A...>(task_index, ntasks, (unsigned long long int) ntasks, func, args...);
  }

  template <typename F, typename... A>
  void xThread::run(unsigned long long int ntasks, F& func, A&... args) {
    if (ntasks == 0) return;
    unsigned long long int task_index = 0;
    run<unsigned long long int, F, A...>(task_index, ntasks, ntasks, func, args...);
  }

  /// @brief Overload for running threads using iterables
  ///
  /// @param it the iterable over which to parallelize
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  template <typename IT, typename F, typename... A>
  void xThread::run(IT& it, F& func, A&... args) {
    unsigned long long int ntasks = (unsigned long long int) std::distance(it.begin(), it.end());
    typename IT::iterator start = it.begin();
    typename IT::iterator end = it.end();
    run(start, end, ntasks, func, args...);
  }

  /// @brief Executes a function over at least ncpus_min and at most ncus_max threads
  ///
  /// @param it The iterator or task index
  /// @param end The end() iterator or the last task index
  /// @param ntasks Number of tasks (for the progress bar)
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  template <typename I, typename F, typename... A>
  void xThread::run(I& it, I& end, unsigned long long int ntasks, F& func, A&... args) {
    if (ntasks == 0) return;
    // First check if enough threads are available using XHOST.CPU_active,
    // which every multi-threaded call should update.
    // This prevents threaded functions that spawn other multi-threaded
    // processes from allocating more threads than the machine can afford.
    // This only works within a single AFLOW run
    uint sleep_second = 10;
    int ncpus_max_available = init::GetCPUCores();
    int ncpus_available = ncpus_max_available - XHOST.CPU_active;
    int ncpus = 0;
    do {
      xthread_cpu_check.lock();
      ncpus_available = ncpus_max_available - XHOST.CPU_active;
      if (ncpus_available >= ncpus_min) {
        ncpus = (ncpus_available > ncpus_max)?ncpus_max:ncpus_available;
        // "reserve" threads globally
        XHOST.CPU_active += ncpus;
        xthread_cpu_check.unlock();
      } else {
        xthread_cpu_check.unlock();
        aurostd::Sleep(sleep_second);
      }
    } while (ncpus_available < ncpus_min);

    // Initialize progress bar
    if (progress_bar_set) initializeProgressBar(ntasks);

    if (ncpus > 1) {
      vector<std::thread*> threads;
      for (int i = 0; i < ncpus; i++) {
        threads.push_back(new std::thread(&xThread::spawnWorker<I, F, A...>, this,
                                          std::ref(it), std::ref(end), ntasks,
                                          std::ref(func), std::ref(args)...)
        );
      }
      for (std::thread* t : threads) {
        t->join();
        delete t;
      }
    } else {
      // No need for thread overhead when ncpus == 1
      for (; it != end; ++it) {
        func(it, args...);
        if (progress_bar_set) pflow::updateProgressBar(++progress_bar_counter, ntasks, *progress_bar);
      }
    }

    // "free" threads globally
    std::lock_guard<std::mutex> lk(xthread_cpu_check);
    XHOST.CPU_active -= ncpus;
  }

}

namespace xthread {

  //apl::AtomicDisplacements::calculateEigenvectorsInThread
  //apl::DOSCalculator::calculateInOneThread
  //apl::PhononDispersionCalculator::calculateInOneThread
  //apl::TCONDCalculator::calculateTransitionProbabilitiesIsotope
  template void xThread::run<
    std::function<void(int)>
  >(int, std::function<void(int) >&
  );

  //apl::PhononCalculator::calculateGroupVelocitiesThread
  template void xThread::run<
    std::function<void(int, vector<vector<double> >&, vector<xmatrix<xcomplex<double> > >&, vector<vector<xvector<double> > >&)>,
    vector<vector<double> >,
    vector<xmatrix<xcomplex<double> > >,
    vector<vector<xvector<double> > >
  >(int, std::function<void(int, vector<vector<double> >&, vector<xmatrix<xcomplex<double> > >&, vector<vector<xvector<double> > >&)>&,
    vector<vector<double> >&,
    vector<xmatrix<xcomplex<double> > >&,
    vector<vector<xvector<double> > >&
  );

  //apl::TCONDCalculator::calculateTransitionProbabilitiesPhonon
  template void xThread::run<
    std::function<void(int, vector<vector<vector<vector<double> > > >&, const vector<vector<vector<xcomplex<double> > > >&)>,
    vector<vector<vector<vector<double> > > >,
    vector<vector<vector<xcomplex<double> > > >
  >(int, std::function<void(int, vector<vector<vector<vector<double> > > >&, const vector<vector<vector<xcomplex<double> > > >&)>&,
    vector<vector<vector<vector<double> > > >&,
    vector<vector<vector<xcomplex<double> > > >&
  );

  //POccCalculator::calculatePhononDOSThread
  template void xThread::run<
    std::function<void(int, const vector<uint>&, const aurostd::xoption&, vector<apl::DOSCalculator>&, vector<xDOSCAR>&)>,
    vector<uint>,
    aurostd::xoption,
    vector<apl::DOSCalculator>,
    vector<xDOSCAR>
  >(int, std::function<void(int, const vector<uint>&, const aurostd::xoption&, vector<apl::DOSCalculator>&, vector<xDOSCAR>&)>&,
    vector<uint>&,
    aurostd::xoption&,
    vector<apl::DOSCalculator>&,
    vector<xDOSCAR>&
  );

  //apl::TCONDCalculator::calculateAnharmonicRates
  template void xThread::run<std::function<void(int, const vector<vector<double> >&, vector<vector<double> >&)>,
    const vector<vector<double> >,
    vector<vector<double> >
  >(int, std::function<void(int, const vector<vector<double> >&, vector<vector<double> >&)>&,
    const vector<vector<double> >&,
    vector<vector<double> >&
  );

  //apl::TCONDCalculator::calculateDelta
  template void xThread::run<
    std::function<void(int, const vector<vector<double> >&, const vector<vector<xvector<double> > >&, vector<vector<xvector<double> > >&)>,
    const vector<vector<double> >,
    vector<vector<xvector<double> > >,
    vector<vector<xvector<double> > >
  >(int, std::function<void(int, const vector<vector<double> >&, const vector<vector<xvector<double> > >&, vector<vector<xvector<double> > >&)>&,
    const vector<vector<double> >&,
    vector<vector<xvector<double> > >&,
    vector<vector<xvector<double> > >&
  );


  //aflowlib::AflowDB::createTable
  template void xThread::run<
    std::function<void(int, const vector<string>&, const vector<string>&)>,
    vector<string>,
    vector<string>
  >(int, std::function<void(int, const vector<string>&, const vector<string>&)>&,
    vector<string>&,
    vector<string>&
  );

  //aflowlib::AflowDB::getColStats
  template void xThread::run<
    std::function<void(int, const vector<string>&, vector<aflowlib::DBStats>&)>,
    const vector<string>,
    vector<aflowlib::DBStats>
  >(int, std::function<void(int, const vector<string>&, vector<aflowlib::DBStats>&)>&,
    const vector<string>&,
    vector<aflowlib::DBStats>&
  );

  template void xThread::run<
    void(int, vector<int>&),
    vector<int>
  >(int, void(&) (int, vector<int>&), vector<int>&);

  template void xThread::run<
    std::vector<int>,
    void(std::vector<int>::iterator&, int),
    int
  >(std::vector<int>&, void (&) (std::vector<int>::iterator&, int), int&);
}

#endif
