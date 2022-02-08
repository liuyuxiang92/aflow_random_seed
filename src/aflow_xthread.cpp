// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************
// Written by Marco Esters
//
// Contains thread manager class xThread, which executes functions in parallel
// and handles all index/iterator management and progress bar updates.
//
// It supports two modes: on-the-fly task distribution (default) or with
// pre-distributed task lists. The latter is useful for many small tasks or
// if the tasks involve expensive pre-computation steps.
//
// ----------
//
// Usage notes:
//
// Requirements for the functions that can be run:
// For on-the-fly redistribution, the function must have either an index
// (integer type) or an iterable to parallelize over as its first parameter.
//
// For predistributed execution, the function must have two integer types
// (start and end points) as first paramters. The pre-distributed scheme
// cannot be used with iterators and currently does not support progress bars.
//
// Additionally, function inputs cannot be prvalues.
//
// ----------
//
// Running functions with xThread:
//
// A function can be run in parallel with the on-the-fly scheme using the
// run() function.
//
// For functions with an index: run(max_index, function, args...)
// For functions over an iterable: run(iterable, function, args...)
//
// For the pre-distributed scheme: runPredistributed(max_index, function, args...)
//
// Every function handled by this class needs to be instantiated here to avoid
// linker errors. There is a section on how to do this at the end of this file.
//
// Static functions are called differently than non-static member functions.
// The following examples show the difference.
//
// Example functions:
// void f1(uint i, const vector<int>& vint, vector<double>& vdbl)
//   Parallelize over vdbl
// void f2(vector<int>::iterator& i, const xmatrix<double>& mdbl)
//   Parallelize over vector<int> v that i will iterate over
// void f3(int i)
//   With ntasks number of tasks.
//
// If functions are static functions:
// Static member functions can directly be plugged into run()
//
// f1:
//   uint ntasks = vdbl.size();
//   xThread xt;
//   xt.run(ntasks, f1, vint, vdbl);
//
// f2:
//   xThread xt;
//   xt.run(v, mdbl);
//
// f3:
//   xThread xt;
//   xt.run(ntasks, f3);
//
// Non-static member functions with xThread:
//
// Member functions of a class cannot be directly plugged into run because they
// have to be bound to an instance of the class. For example, let all functions
// belong to class C instantiated as cls. Let f1 and f2 be called inside another
// class function of C and let f3 be called outside.
//
// f1:
//   uint ntasks = vdbl.size();
//   xThread xt;
//   std::function<void(uint, const vector<int>&, const vector<double>&)> fn1 =
//     std::bind(&C::f1, this,
//       std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
//   xt.run(ntasks, fn1, vint, vdbl);
//
// f2:
//   xThread xt;
//   std::function<void(vector<int>::iterator&, const xmatrix<double>&)> fn2 =
//     std::bind(&C::f2, this, std::placeholders::_1, std::placeholders::_2);
//   xt.run(v, fn2, mdbl);
//
// f3:
//   xThread xt;
//   std::function<void(int)> fn3 = std::bind(&C::f3, cls*, std::placeholders::_1);
//   xt.run(ntasks, fn3);
//
// The template parameter in std::function takes the function return type and
// the type of the function inputs exactly as written in the function defintion.
// std::bind takes an address to the function, a pointer to the class instance
// and one std::placeholders::_N for each argument of the function.
// Using namespace std::placeholders can make std::bind more legible.
//
// ----------
//
// Setting the number of CPUs:
// The default constructor uses as many threads as possible, but the number of
// CPUs can be passed into the constructor or changes via setCPUs(). A minimum
// number of threads can be set as well. In that case, xThread will wait until
// that minimum number of threads is available.
//
// ----------
//
// Progress bars:
// To use a progress bar, pass an ostream into the setProgressBar() function
// before calling run(). unsetProgressBar() removes the progress bar.
//
// ----------
//
// Thread safety:
// xThread guarantees only that no two instances of the called function writes
// to the same index or iterator using mutexes. It cannot guarantee that the
// passed function itself is thread-safe. Functions that perform actions that
// are not thread-safe should pass their own mutex as a parameter.
//
// ----------

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
  /// @param nmin Mininum number of CPUs required to spawn thread workers default: 1)
  xThread::xThread(int nmax, int nmin) {
    free();
    setCPUs(nmax, nmin);
  }

  xThread::xThread(const xThread& that) {
    copy(that);
  }

  const xThread& xThread::operator=(const xThread& that) {
    copy(that);
    return *this;
  }

  void xThread::copy(const xThread& that) {
    if (this == &that) return;
    // std::mutex should not be copied because
    // it needs to stay immutable
    ncpus_max = that.ncpus_max;
    ncpus_max = that.ncpus_min;
    progress_bar = that.progress_bar;
    progress_bar_counter = that.progress_bar_counter;
    progress_bar_set = that.progress_bar_set;
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

  // Progress bar -------------------------------------------------------------

  /// @brief Tells xThread that it has to update a progress bar when running
  ///
  /// @param oss The ostream for the progress bar
  void xThread::setProgressBar(ostream& oss) {
    progress_bar = &oss;
    progress_bar_set = true;
    progress_bar_counter = 0;
  }

  void xThread::unsetProgressBar() {
    progress_bar = nullptr;
    progress_bar_set = false;
    progress_bar_counter = 0;
  }

  void xThread::initializeProgressBar(unsigned long long int ntasks) {
    progress_bar_counter = 0;
    pflow::updateProgressBar(0, ntasks, *progress_bar);
  }

  // CPU functions ------------------------------------------------------------

  /// @brief Sets the minimum and maximum number of CPUs used for threading
  ///
  /// @param nmax Maximum number of CPUs used by xThread (default: 0 for all available CPUs)
  /// @param nmin Mininum number of CPUs required to spawn thread workers default: 1)
  void xThread::setCPUs(int nmax, int nmin) {
    if (nmax < nmin) std::swap(nmax, nmin);
    ncpus_max = (nmax > 0)?nmax:(KBIN::get_NCPUS());
    ncpus_min = (nmin > 0)?nmin:nmax;
  }

  /// @brief Checks if enough threads are available and "reserves" them
  /// using XHOST.CPU_active
  ///
  /// @param ntasks Number of tasks
  ///
  /// @return Number of threads to use for the function
  ///
  /// This prevents threaded functions that spawn other multi-threaded
  /// processes from allocating more threads than the machine can afford.
  /// This only works within a single AFLOW run
  int xThread::reserveThreads(unsigned long long int ntasks) {
    uint sleep_second = 10;
    int ncpus_max_available = KBIN::get_NCPUS();
    int nmax = ncpus_max;
    int nmin = ncpus_min;
    // Adjust max. and min. number of CPUs to the number of tasks
    if ((ntasks <= (unsigned long long int) AUROSTD_MAX_INT)) {
      int n = (unsigned long long int) ntasks;
      if (n < nmin) nmin = n;
      if (n < nmax) nmax = n;
    }
    int ncpus = 0, ncpus_available = 0;
    do {
      xthread_cpu_check.lock();
      ncpus_available = ncpus_max_available - XHOST.CPU_active;
      if (ncpus_available >= nmin) {
        ncpus = std::min(ncpus_available, nmax);
        // "reserve" threads globally
        XHOST.CPU_active += ncpus;
        xthread_cpu_check.unlock();
      } else {
        xthread_cpu_check.unlock();
        aurostd::Sleep(sleep_second);
      }
    } while (ncpus_available < nmin);
    return ncpus;
  }

  /// @brief "Frees" threads in XHOST.CPU_active
  ///
  /// @param ncpus Number of threads to free
  void xThread::freeThreads(int ncpus) {
    std::lock_guard<std::mutex> lk(xthread_cpu_check);
    XHOST.CPU_active -= ncpus;
  }

  // On-the-fly scheme --------------------------------------------------------

  /// @brief Overload for running threads using indices to keep track of tasks
  ///
  /// @param ntasks The number of tasks over which the function is parallelized
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
  ///
  /// ntasks will be converted to unsigned long long int because
  /// the progress bar counter in pflow is of that type.
  template <typename IT, typename F, typename... A>
  void xThread::run(IT& it, F& func, A&... args) {
    int dist = (int) std::distance(it.begin(), it.end());
    if (dist <= 0) return; // Cannot iterate backwards (yet)
    unsigned long long int ntasks = (unsigned long long int) dist;
    typename IT::iterator start = it.begin();
    typename IT::iterator end = it.end();
    run(start, end, ntasks, func, args...);
  }

  /// @brief Executes a function over at least ncpus_min and at most ncus_max threads
  ///
  /// @param it The iterator or task index
  /// @param end The end() iterator or the last task index
  /// @param ntasks Number of tasks
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  template <typename I, typename F, typename... A>
  void xThread::run(I& it, I& end, unsigned long long int ntasks, F& func, A&... args) {
    if (ntasks == 0) return;
    int ncpus = reserveThreads(ntasks);
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

    freeThreads(ncpus);
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

  // Pre-distributed scheme ---------------------------------------------------

  /// @brief Overload for running threads using indices to keep track of tasks
  ///
  /// @param ntasks The number of tasks over which the function is parallelized
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  ///
  /// The overloads are not necessary, but they allow for the same instantiation
  /// scheme as the on-the-fly methods.
  template <typename F, typename... A>
  void xThread::runPredistributed(int ntasks, F& func, A&... args) {
    if (ntasks <= 0) return;
    runPredistributed<int, F, A...>(ntasks, func, args...);
  }

  template <typename F, typename... A>
  void xThread::runPredistributed(uint ntasks, F& func, A&... args) {
    if (ntasks == 0) return;
    runPredistributed<uint, F, A...>(ntasks, func, args...);
  }

  template <typename F, typename... A>
  void xThread::runPredistributed(unsigned long long int ntasks, F& func, A&... args) {
    if (ntasks == 0) return;
    runPredistributed<unsigned long long int, F, A...>(ntasks, func, args...);
  }

  template <typename I, typename F, typename... A>
  void xThread::runPredistributed(I ntasks, F& func, A&... args) {
    if (ntasks <= 0) return;
    int ncpus = reserveThreads(ntasks);

    if (ncpus > 1) {
      vector<std::thread*> threads;
      I n = (I) ncpus; // convert ncpus to same type
      I tasks_per_thread = ntasks/n;
      I remainder = ntasks % n;
      I startIndex = 0, endIndex = 0;
      for (I t = 0; t < n; t++) {
        if (t < remainder) {
          startIndex = (tasks_per_thread + 1) * t;
          endIndex = startIndex + tasks_per_thread + 1;
        } else {
          startIndex = tasks_per_thread * t + remainder;
          endIndex = startIndex + tasks_per_thread;
        }
        threads.push_back(new std::thread(std::ref(func), startIndex, endIndex, std::ref(args)...));
      }
      for (std::thread* t : threads) {
        t->join();
        delete t;
      }
    } else {
      // No need for thread overhead when ncpus == 1
      func(0, ntasks, args...);
    }

    freeThreads(ncpus);
  }

}

// Template Instantiation
//
// Functions run with xThread need to be instantiated to avoid linker errors.
// The way the instantiation needs to be done depends on whether an index or
// an iterator is used, on whether the function is static or a std::function,
// and on the origin of the arguments passed into run().
//
// On-the-fly and pre-distributed schemes are instantiated the same way though.
//
// ----------
//
// Instantiating functions using an index:
//
// It is best to learn by example. Take the following function f1:
//
// void f1(int, const vector<int>&, const vector<double>&, vector<int>&, vector<double>&);
//
// Let if be called inside function f2:
//
// bool f2(const vector<int>& vint1, vector<double>& vdbl2) {
//   vector<int> vint2;
//   vector<double> vdbl1;
//   uint ntasks = vdbl.size();
//   xthread::xThread xt;
//   xt.run(ntasks, f1, vint1, vdbl1, vint2, vdbl2);
// }
//
// The instantiation is:
//
// template void xThread::run<
//   void(int, const vector<int>&, const vector<double>&, vector<int>&, vector<double>&),
//   const vector<int>,
//   vector<double>,
//   vector<int>,
//   vector<double>
// >(uint, void(&) (int, const vector<int>&, const vector<double>&, vector<int>&, vector<double>&),
//   const vector<int>&,
//   vector<double>&,
//   vector<int>&,
//   vector<double>&
// );
//
// The first parameter inside <> is function type with the argument types in
// parentheses. The remaining parameters are the argument types but without &.
// The usage of const is different (see explanation later).
//
// The first parameter inside () is the index type, then the function type with
// the argument types. This time, the function type has an additional (&). The
// parentheses are mandatory. The remaining parameters are the argument types,
// but with & this time.
//
// The remaining parameters do not follow the function definition since the
// first vector<double> does not get const, even though one is const & in the
// function definition. On the other hand, the first vector<int> is const &.
// The important part is how it passed into xt.run(), not how it is passed into
// f1. vdbl1 is created inside f2 as vector<double> and is thus passed into
// run() as vector<double>&. vint1 on the other hand is passed into f2 as
// const vector<int>& and will thus be passed into run() as const vector<int>&.
//
// So, if an argument is created in the same function that calls run(), do not
// use const inside (). If it is passed down by another function as const &,
// use const &.
//
//
// If f2 was converted to a std::function, e.g. via std::bind for non-static
// member function, the instantiation is:
//
// template void xThread::run<
//   uint, std::function<void(int, const vector<int>&, const vector<double>&, vector<int>&, vector<double>&)>,
//   const vector<int>,
//   vector<double>,
//   vector<int>,
//   vector<double>
// >(uint, std::function<(int, const vector<int>&, const vector<double>&, vector<int>&, vector<double>&)>&,
//   const vector<int>&,
//   vector<double>&,
//   vector<int>&,
//   vector<double>&
// );
//
// Note: "template void" will always be "template void" since it is the type of
// xThread::run(). It does not depend on the types of f1 and f2.
//
// ----------
//
// Instantiating functions using an iterator:
//
// Take the function f1:
//
// void f1(vector<int>::iterator& it, const vector<double>&);
//
// Let it be called as:
//
//   vector<int> v1;
//   vector<double> v2;
//   xthread::xThread xt;
//   xt.run(v1, v2);
//
// The instantiation is then:
//
//   template void xThread::run<
//     vector<int>,
//     void(vector<int>::iterator&, const vector<double>&),
//     vector<double>
//   >(
//     vector<int>&,
//     void (&) (vector<int>::iterator&, const vector<double>&),
//     vector<double>&
//   );
//
// Note that the iterable has to appear as the first template parameter.
// The rules and caveats for member functions and const references are the same
// as for functions using an index.
//
// ----------
//
// All instantiations should be added below with a comment on which function is
// instantiated. One instantiation can cover multiple functions.
//
// ----------
//
// Common mistakes
//
// When encountering linker errors, check if:
//
// - the index variable is exactly the same type (no implicit conversions)
// - all arguments inside () have an & (except for the index variable)
// - no parameter inside <> has an &
// - this includes the function inside () - double-check
// - the parameters inside the function are identical to the function definition
// - only use const when passed in by parent function or when defined as const
// - const inside () and <> are used the same way
// - the iterable is added before functions inside <>
// - on the other hand, the index variable type is not to be added inside <>
//
// ----------

namespace xthread {

  // run ----------------------------------------------------------------------

  //apl::AtomicDisplacements::calculateEigenvectorsInThread
  //apl::DOSCalculator::calculateInOneThread
  //apl::PhononDispersionCalculator::calculateInOneThread
  //apl::TCONDCalculator::calculateTransitionProbabilitiesIsotope
  template void xThread::run<
    std::function<void(int)>
  >(int, std::function<void(int)>&
  );

  //lambda function inside aurostd::multithread_execute
  template void xThread::run<
    vector<string>,
    std::function<void(vector<string>::iterator&)>
  >(vector<string>&, std::function<void(vector<string>::iterator&)>&);

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
    std::function<void(int, const vector<uint>&, const aurostd::xoption&, vector<apl::DOSCalculator>&, vector<xDOSCAR>&, std::mutex&)>,
    vector<uint>,
    aurostd::xoption,
    vector<apl::DOSCalculator>,
    vector<xDOSCAR>,
    std::mutex
  >(int, std::function<void(int, const vector<uint>&, const aurostd::xoption&, vector<apl::DOSCalculator>&, vector<xDOSCAR>&, std::mutex&)>&,
    vector<uint>&,
    aurostd::xoption&,
    vector<apl::DOSCalculator>&,
    vector<xDOSCAR>&,
    std::mutex&
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

  //XtalFinderCalculator::performStructureConversions
  template void xThread::run<
    std::function<void(uint, const vector<bool>&, const vector<bool>&, const vector<bool>&)>,
    vector<bool>,
    vector<bool>,
    vector<bool>
  >(uint, std::function<void(uint, const vector<bool>&, const vector<bool>&, const vector<bool>&)>&,
    vector<bool>&,
    vector<bool>&,
    vector<bool>&
  );

  //XtalFinderCalculator::runComparisonThreads
  template void xThread::run<
    std::function<void(uint,
      vector<StructurePrototype>&,
      const vector<std::pair<uint, uint> >&,
      const vector<std::pair<uint, uint> >&,
      bool, bool, bool)>,
    vector<StructurePrototype>,
    vector<std::pair<uint, uint> >,
    vector<std::pair<uint, uint> >,
    bool, bool, bool
  >(
    uint,
    std::function<void(uint,
      vector<StructurePrototype>&,
      const vector<std::pair<uint, uint> >&,
      const vector<std::pair<uint, uint> >&,
      bool, bool, bool)>&,
    vector<StructurePrototype>&,
    vector<std::pair<uint, uint> >&,
    vector<std::pair<uint, uint> >&,
    bool&, bool&, bool&
  );

  //XtalFinderCalculator::getPrototypeDesignations
  template void xThread::run<
    vector<StructurePrototype>,
    std::function<void(vector<StructurePrototype>::iterator&)>
  >(
    vector<StructurePrototype>&,
    std::function<void(vector<StructurePrototype>::iterator&)>&
  );

  //XtalFinderCalculator::getMatchingAFLOWPrototypes
  template void xThread::run<
    std::function<void(uint, vector<StructurePrototype>&, const aurostd::xoption&)>,
    vector<StructurePrototype>,
    aurostd::xoption
  >(uint, std::function<void(uint, vector<StructurePrototype>&, const aurostd::xoption&)>&,
    vector<StructurePrototype>&,
    aurostd::xoption&
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

  // runPredistributed --------------------------------------------------------

  //XTalFinderCalculator::calculateSpaceGroups
  template void xThread::runPredistributed<
    std::function<void(uint, uint, uint)>,
    uint
  >(
    uint, std::function<void(uint, uint, uint)>&, uint&
  );

  //XtalFinderCalculator::computeLFAEnvironments
  //XtalFinderCalculator::calculateNearestNeighbors
  template void xThread::runPredistributed<
    std::function<void(uint, uint)>
  >(
    uint, std::function<void(uint, uint)>&
  );

  //XtalFinderCalculator::searchAtomMappings
  template void xThread::runPredistributed<
    std::function<bool(
      uint, uint, const xstructure&,
      const vector<double>&, const xstructure&, const string&,
      vector<xmatrix<double> >&, vector<structure_mapping_info>&, bool, bool
    )>,
    xstructure,
    vector<double>,
    xstructure,
    string,
    vector<xmatrix<double> >,
    vector<structure_mapping_info>,
    bool, bool
  >(
    uint, std::function<bool(uint, uint, const xstructure&,
      const vector<double>&, const xstructure&, const string&,
      vector<xmatrix<double> >&, vector<structure_mapping_info>&, bool, bool)>&,
    xstructure&,
    vector<double>&,
    xstructure&,
    string&,
    vector<xmatrix<double> >&,
    vector<structure_mapping_info>&,
    bool&, bool&
  );

  //aflowlib::AflowDB::getColStats
  template void xThread::runPredistributed<
    std::function<void(
      int, int, const vector<string>&, vector<aflowlib::DBStats>&
    )>,
    const vector<string>,
    vector<aflowlib::DBStats>
  >(int,
    std::function<void(
      int, int, const vector<string>&, vector<aflowlib::DBStats>&
    )>&,
    const vector<string>&,
    vector<aflowlib::DBStats>&
  );

}

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************
