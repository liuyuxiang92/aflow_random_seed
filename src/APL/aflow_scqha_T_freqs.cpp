// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                Aflow PINKU NATH - Duke University 2014-2016             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu

#include "aflow_apl.h"
#include <iterator>

#define _isnegative(a) (a<MIN_EIGEN_TRESHOLD) ? true : false

#undef AFLOW_APL_MULTITHREADS_ENABLE

#if GCC_VERSION >= 40400   // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

//This class computes temperature dependent PDIS curves using SCQHA temperature and volume data.

namespace apl
{
  // ***************************************************************************************
  T_spectra_SCQHA_QHA3P::T_spectra_SCQHA_QHA3P( PhononCalculator& pc, QHA_AFLOWIN_CREATOR& runeos, Logger& l): _pc(pc),_runeos(runeos),_logger(l)
  {
    //clear all memories before using
    clear();
  }
  // ***************************************************************************************
  T_spectra_SCQHA_QHA3P::~T_spectra_SCQHA_QHA3P()
  {
    this->clear();
  }
  // ***************************************************************************************
  void T_spectra_SCQHA_QHA3P::clear()
  {
    _is_negative_freq=false;
    _DMp.clear();
    _DMm.clear();
    _nBranches=0;
    _kpoints.clear();
    _gp_path_test.clear();
    _freqs_path.clear();
    _freqs_pathM.clear();
    _freqs_pathP.clear();
    _freqs_T.clear();
    _d1fdv1.clear();
    _d2fdv2.clear();
    _qha_gpdir.clear();
    _qha_gpvol.clear();
    _TV.clear();
    _delta_V=0.0;
    _V0=0.0;
    _tmp_dir="";
    _cutoff_freq=0.0;
  }
  // ***************************************************************************************
  void T_spectra_SCQHA_QHA3P::get_tmp_dir_name(const string dir)
  {
    _tmp_dir=dir;
  }
  // ***************************************************************************************
  bool T_spectra_SCQHA_QHA3P::calculation_freqs(const std::vector< aurostd::xvector<double> > &kpoints)
  {
    _logger<<"Calculating SCQHA temperature dependent phonon spectra, "<<apl::endl;
    _kpoints=kpoints;

    if(_kpoints.size()==0)
    {
      _logger<< apl::error <<"_kpoints.size()==0"<<apl::endl;
      return false;
    }

    if(!calculation_freqs()) return false;
    return true;
  }
  // ***************************************************************************************
  bool T_spectra_SCQHA_QHA3P::calculation_freqs()
  {
    if(_kpoints.size()==0)
    {
      _logger<<apl::error<<"kpoints.size()=0"<<apl::endl;
      return false;
    }
    return get_dynamicalmatrices_along_path();
  }
  // ***************************************************************************************
  bool T_spectra_SCQHA_QHA3P::set_imported_variables()
  {
    _logger<<"Importing variables"<<apl::endl;
    _nBranches=_pc.getNumberOfBranches();
    _qha_gpdir = _runeos.get_scqha_dir_names();
    _qha_gpvol = _runeos.get_scqha_volumes();

    if( (_qha_gpdir.size()==0))
    {
      _logger<<  apl::warning <<"directory size is zero remove LOCK and apl.xml and run again"<<apl::endl;
      return false;
    }
    if((_qha_gpvol.size()==0))
    {
      _logger<<  apl::warning<<"volume size is zero remove LOCK and apl.xml and run again"<<apl::endl;
      return false;
    }
    if((_qha_gpvol.size()!=3))
    {
      _is_vol_err=true;
      _logger<<  apl::warning<< "qha_gpvol.size()!=3, SCQHA calculations skipped"<<apl::endl;
      return false;
    }

    _delta_V=0.5*((_qha_gpvol[2]-_qha_gpvol[1])+(_qha_gpvol[1]-_qha_gpvol[0]));
    _V0=_qha_gpvol[1];

    return true;
  }
  // ***************************************************************************************
  //read dynamical matrices generated by other directories
  bool T_spectra_SCQHA_QHA3P::get_dynamicalmatrices_along_path()
  {
    xmatrix< xcomplex<double> > dd(_nBranches,_nBranches,1,1);
    _DMp.resize(_kpoints.size(), dd);
    _DMm.resize(_kpoints.size(), dd);

    //dynamic matrix file names
    string DMmfile=string(_tmp_dir)+string("/")+string("DM_")+string("path.")+string(_qha_gpdir[0]);
    string DMpfile=string(_tmp_dir)+string("/")+string("DM_")+string("path.")+string(_qha_gpdir[1]);

    _logger <<"Reading dynamical matrices "<<apl::endl;
    //else if(_is_scqha)_logger<<_SCQHA_MSG <<"Reading dynamical matrices "<<apl::endl;

    if(!read_matrix(_DMp,DMpfile)) return false;
    if(!read_matrix(_DMm,DMmfile)) return false;

    xvector<double> d(_nBranches,1);
    _freqs_path.resize(_kpoints.size(),d);
    _freqs_pathM.resize(_kpoints.size(),d);
    _freqs_pathP.resize(_kpoints.size(),d);
    _freqs_T.resize(_kpoints.size(),d);
    _d1fdv1.resize(_kpoints.size(),d);
    _d2fdv2.resize(_kpoints.size(),d);

    bool gppass=false;

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    // Get the number of CPUS
    int ncpus = sysconf(_SC_NPROCESSORS_ONLN);// AFLOW_MachineNCPUs;
    if(ncpus<1) ncpus=1;
    //    int qpointsPerCPU = _kpoints.size() / ncpus;  OBSOLETE ME180801
    _gp_path_test.resize(ncpus, false);
    // Show info 
    string msg="";
    if( ncpus == 1 )
    {
      msg="Calculating dynamical matrices along path";
      _logger.initProgressBar(msg);
    }
    else
    {
      msg="Calculating dynamical matrices along path (" + stringify(ncpus) + " threads)";
      _logger.initProgressBar(msg);
    }

    // Distribute the calculation
    int startIndex, endIndex;
    std::vector< std::thread* > threads;
    vector<vector<int> > thread_dist = getThreadDistribution((int) _kpoints.size(), ncpus);
    for (int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = thread_dist[icpu][0];
      endIndex = thread_dist[icpu][1];
      threads.push_back( new std::thread(&T_spectra_SCQHA_QHA3P::get_freqs_in_threads,this,startIndex,endIndex, icpu) );
    }

    /* OBSOLETE ME180801
       for(int icpu = 0; icpu < ncpus; icpu++) {
       startIndex = icpu * qpointsPerCPU;
       endIndex = startIndex + qpointsPerCPU;
       if( ( (uint)endIndex > _kpoints.size() ) ||
       ( ( icpu == ncpus-1 ) && ( (uint)endIndex < _kpoints.size() ) ) )
       endIndex = _kpoints.size();
       threads.push_back( new std::thread(&T_spectra_SCQHA_QHA3P::get_freqs_in_threads,this,startIndex,endIndex, icpu) );
       }
    // Wait to finish all threads here!
    for(uint i = 0; i < threads.size(); i++) {
    threads[i]->join();
    delete threads[i];
    }
    */

    // Done
    _logger.finishProgressBar();


    uint bool_size=0;
    for(uint id=0; id!=_gp_path_test.size(); id++)
      bool_size+=_gp_path_test[id];
    if(bool_size==(uint)ncpus)gppass=true;
    _gp_path_test.clear();
#else
    _gp_path_test.resize(1, false);
    get_freqs_in_threads(0, (int)_kpoints.size(), 0);
    uint bool_size=0;
    for(uint id=0; id!=_gp_path_test.size(); id++)
      bool_size+=_gp_path_test[id];
    if(bool_size==1)gppass=true;  
#endif 
    return gppass;
  }
  // ***************************************************************************************
  void T_spectra_SCQHA_QHA3P::get_freqs_in_threads(int startIndex, int endIndex, int cpuid)
  {
    double dv2=0.25*_delta_V*_delta_V;

    for(int In=startIndex;In<endIndex;In++)
    {
      xvector<double> qpoint;
      xmatrix<xcomplex<double> >  DM0(_nBranches,_nBranches,1,1);//at volume 0
      qpoint=_kpoints[In];
      DM0=_pc.getDynamicalMatrix(qpoint);

      xvector<double> eigenvalues(_nBranches, 1);
      xmatrix<xcomplex<double> > eigenvectors(_nBranches, _nBranches, 1,1);

      xvector<double> eigenvaluesM(_nBranches, 1);
      xmatrix<xcomplex<double> > eigenvectorsM(_nBranches, _nBranches, 1,1);

      xvector<double> eigenvaluesP(_nBranches, 1);
      xmatrix<xcomplex<double> > eigenvectorsP(_nBranches, _nBranches, 1,1);

      apl::aplEigensystems e;
      e.eigen_calculation(DM0, eigenvalues, eigenvectors, APL_MV_EIGEN_SORT_VAL_ASC);


      apl::aplEigensystems eM;
      eM.eigen_calculation(_DMm[In], eigenvaluesM, eigenvectorsM, APL_MV_EIGEN_SORT_VAL_ASC);

      apl::aplEigensystems eP;
      eP.eigen_calculation(_DMp[In], eigenvaluesP, eigenvectorsP, APL_MV_EIGEN_SORT_VAL_ASC);

      for(uint j=1; j<=_nBranches; j++)
      {
        if(_iszero(eigenvalues[j]))eigenvalues[j]=0.0;
        if(eigenvalues[j]<0){
          if(_isnegative(eigenvalues[j]))
          {
            _logger<<  apl::warning <<"negative eigenvalue: = " <<eigenvalues[j]<<apl::endl;
            _is_negative_freq=true;
            return;
          } else eigenvalues[j]=0.00;
        }
        _freqs_path[In][j]=sqrt(eigenvalues[j])*RAW2Hz;
      }// nBranch loop end

      for(uint j=1; j<=_nBranches; j++)
      {
        if(_iszero(eigenvaluesM[j]))eigenvaluesM[j]=0.0;
        if(eigenvaluesM[j]<0){
          if(_isnegative(eigenvaluesM[j]))
          {
            _logger<<  apl::warning <<"negative eigenvalueM: = " <<eigenvaluesM[j]<<apl::endl;
            _is_negative_freq=true;
            return;
          } else eigenvaluesM[j]=0.00;
        }
        _freqs_pathM[In][j]=sqrt(eigenvaluesM[j])*RAW2Hz;
      }// nBranch loop end

      for(uint j=1; j<=_nBranches; j++)
      {
        if(_iszero(eigenvaluesP[j]))eigenvaluesP[j]=0.0;
        if(eigenvaluesP[j]<0){
          if(_isnegative(eigenvaluesP[j]))
          {
            _logger<<  apl::warning <<"negative eigenvalueP: = " <<eigenvaluesP[j]<<apl::endl;
            _is_negative_freq=true;
            return;
          } else eigenvaluesP[j]=0.00;
        }
        _freqs_pathP[In][j]=sqrt(eigenvaluesP[j])*RAW2Hz;
      }// nBranch loop end

      xvector<double> d1fdv1(_nBranches, 1);
      xvector<double> d2fdv2(_nBranches, 1);
      for(uint j=1; j<=_nBranches; j++)
      {
        _d1fdv1[In][j]= (-_freqs_pathM[In][j] + _freqs_pathP[In][j])/_delta_V;
        _d2fdv2[In][j]= (_freqs_pathM[In][j] + _freqs_pathP[In][j] - 2.0* _freqs_path[In][j])/dv2;
      }
    }

    _gp_path_test[cpuid]=true;
  }//fn end
  // ***************************************************************************************
  void T_spectra_SCQHA_QHA3P::write_T_dispersion(double V_T, double T, const vector <double> &path, const vector <int> &path_seg)
  { 
    double dv=V_T-_qha_gpvol[1];

    for(uint In=0; In<_kpoints.size(); In++)
    {
      for(uint j=1; j<=_nBranches; j++)
      {
        _freqs_T[In][j]=_freqs_path[In][j]+_d1fdv1[In][j]*dv+0.5*_d2fdv2[In][j]*dv*dv;
      }
    }

    if(_freqs_T.size()==0)
    {
      _logger<<apl::error <<"_freqs_T size is zero "<<apl::endl; return;
    }
    else if(path.size()==0)
    {
      _logger<<apl::error <<"path size is zero "<<apl::endl; return;
    }
    else if(path_seg.size()==0)
    {
      _logger<<apl::error <<"path_seg size is zero "<<apl::endl; return;
    }

    string outfile =  "aflow.scqha_pdis_"+NumToStr<double>(T);

    string msg="Writing "+outfile+" ,";

    _logger <<msg << apl::endl;

    vector<string> hash_lines;
    if(!read_PDIS(hash_lines))return ;

    //writing to a file
    stringstream os_gp;
    uint isize=_freqs_T.size();
    uint jsize=_freqs_T[0].rows;

    os_gp<<"#frequencies along high-symmetry q-points. File format is same as PDIS file"<<"\n";
    for(uint i=1; i!=hash_lines.size(); i++)
      os_gp<<hash_lines[i]<<"\n";

    for(uint j=1; j<=jsize; j++){
      if(j==1) os_gp<<"#"<<setw(10)<<"        "<<setw(10) <<"q-path"<<setw(10)<<string("Br")+aurostd::utype2string<uint>(j);
      else os_gp<< setw(15)<<string("Br")+aurostd::utype2string<uint>(j);
    }
    os_gp<<"\n";

    for(uint i=0; i<isize; i++){
      os_gp<< setw(5) << path_seg[i];
      os_gp<< setprecision(6)<<std::fixed << std::showpoint
        << setw(15) << path[i];

      for(uint j=1; j<=jsize; j++){
        os_gp<<setprecision(6)<<std::fixed << std::showpoint
          <<setw(15)<<_freqs_T[i][j];
      }os_gp<<"\n";}

    if(!aurostd::stringstream2file(os_gp, outfile, "WRITE")) {
      // ME191031 - use xerror
      //throw APLRuntimeError("Cannot write aflow.scqha_pdis_T");
      string function = "T_spectra_SCQHA_QHA3P::write_T_dispersion()";
      string message = "Cannot write aflow.scqha_pdis_T";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }
    aurostd::StringstreamClean(os_gp);
  }//fn end
  // ***************************************************************************************
  bool T_spectra_SCQHA_QHA3P::read_matrix(vector<xmatrix<xcomplex<double> > >&A, const string file)
  {
    if (!exists_test0(file) && !aurostd::EFileExist(file)) {
      // ME191031 - use xerror
      //throw apl::APLRuntimeError("T_spectra_SCQHA_QHA3P:: Missing file: "+file);
      string function = "T_spectra_SCQHA_QHA3P::read_matrix()";
      string message = "Missing file: " + file;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
    }

    ifstream in;
    in.open(file.c_str(),ios_base::binary);
    if (!in.is_open()){_logger<<apl::error <<file<<" not able to open "<<apl::endl; return false;}
    for(uint I=0;I<A.size();I++){
      for(uint J=1;J<=_nBranches;J++){
        in.read( (char *)(&A[I][J][1]), _nBranches*sizeof(xcomplex<double>) );
      }}
    in.clear();
    in.close();
    return true;
  }
  // ***************************************************************************************
  bool T_spectra_SCQHA_QHA3P::calculate_pdis_T(const vector <double> &path, const vector <int> &path_seg)
  {
    _logger<<"Calculating temperature dependent SCQHA PDIS "<<apl::endl;

    if(_TV.size()==0){
      _logger<<apl::error <<" _TV.size()==0, no matching temperature to calculate SCQHA PDIS, "<<apl::endl;
      return false;
    }

    for(uint i=0; i!=_TV.size(); i++){
      write_T_dispersion(_TV[i][1], _TV[i][0], path, path_seg);
    }
    return true;
  }
  // ***************************************************************************************
  bool T_spectra_SCQHA_QHA3P::read_PDIS(vector<string> &hash_lines)
  {
    string file = DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_PDIS_FILE;  // ME190428
    if (!exists_test0(file) && !aurostd::EFileExist(file)) {
      // ME191031
      //throw apl::APLRuntimeError("T_spectra_SCQHA_QHA3P::read_PDIS() Missing file: "+file);
      string function = "T_spectra_SCQHA_QHA3P::read_PDIS()";
      string message = "Missing file: " + file;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
    }
    vector<string> vlines;
    aurostd::efile2vectorstring(file, vlines);
    if (!vlines.size()) {
      // ME191031 - use xerror
      //throw apl::APLRuntimeError("T_spectra_SCQHA_QHA3P::read_PDIS() Missing file: "+file);
      string function = "T_spectra_SCQHA_QHA3P::read_PDIS()";
      string message = "Missing file: " + file;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
    } 
    uint line_count = 0;  
    hash_lines.clear();
    string line;
    while (line_count < vlines.size()) {
      line = vlines[line_count++];
      if(line=="")continue;
      if(line[0]=='#')hash_lines.push_back(line);
      else break;
    }
    return true;
  }
  // ***************************************************************************************
  void T_spectra_SCQHA_QHA3P::get_input_data(const std::vector<std::vector<double> > &TV)
  {
    _TV=TV;
  }
  // ***************************************************************************************
  template <typename T>
    string T_spectra_SCQHA_QHA3P::NumToStr ( T Number )
    {
      ostringstream ss;
      ss << Number;
      return ss.str();
    }
  // ***************************************************************************************
  bool T_spectra_SCQHA_QHA3P::exists_test0 (const std::string& name)
  {
    ifstream f(name.c_str());
    if (f.good()) {
      f.close();
      return true;
    } else {
      f.close();
      return false;
    }
  }
  // ***************************************************************************************
}//apl end
