// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                Aflow CORMAC TOHER - Duke University 2013-2020           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AGL_PLOT_WRITE_CPP
#define _AFLOW_AGL_PLOT_WRITE_CPP
#include "aflow.h"
#include "aflow_agl_debye.h"


// ###############################################################################
//                  AFLOW Automatic GIBBS Library (AGL) (2013-2018)
// ###############################################################################
//
// Uses quasi-harmonic Debye model to obtain thermodynamic properties of materials
// Based on original Fortran program written by M. A. Blanco et al.
// See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details of original GIBBS program
// See C. Toher et al., Phys. Rev. B 90, 174107 (2014), Phys. Rev. Materials 1, 015401 (2017) and references therein for description of this AGL implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AGL
//


// **************************************************************************************
//  This set of functions plot data or write input files
// **************************************************************************************

// ***************************************************************************
// AGL_functions::plotgibbsresults
// ***************************************************************************
namespace AGL_functions {
  //
  // plotgibbsresults: uses gnuplot generate x-y plots of data arrays xdata and ydata, in the range xmin to xmax and ymin to ymax
  // The user can choose to save the plot in png, pdf, jpg, gif, or eps format
  //
  uint plotaglresults(vector<double>& xdata, vector<double>& ydata, double& xmin, double& xmax, double& ymin, double& ymax, int& nplotpoints, string& datalabel, string& xlabel, string& ylabel, string& plotname, string& plotfilename, string& sysname, ofstream& FileMESSAGE) {
    ostringstream aus;
    stringstream fdatainss;
    string datafilename = datalabel + ".dat";
    string plotfilenameeps = plotfilename + ".eps";
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Plotting data for " << plotname.c_str() << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    aurostd::StringstreamClean(fdatainss);
    for (int i = 0; i < nplotpoints; i++) {
      fdatainss << xdata.at(i) << "\t" << ydata.at(i) << endl;
    }		  
    if(!aurostd::stringstream2file(fdatainss, datafilename, "WRITE")) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Unable to open file " << datafilename.c_str() <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    
    // Create plot title by replacing "_" by "\_" to avoid it being interpreted as a subscript command
    string plot_title;
    char chwork;
    for (uint i = 0; i < sysname.size(); i++) {
      chwork = sysname.at(i);
      if(chwork != '_') {
	plot_title.push_back(chwork);
      } else {
	plot_title.push_back('\\');
	// OBSOLETE plot_title.push_back('\\');
	plot_title.push_back(chwork);
      }
    }
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ + "plot_title = " << plot_title <<  endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	     

    //***********************************GENERATING GNUPLOT SCRIPT**************************************************************
    // Writing Gnuplot Script
    //
    string scriptname;
    istringstream iss(sysname);
    iss >> scriptname;
    string gnuplotscript = "GNUPLOT_" + scriptname + "_" + datalabel  + ".gp";
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Generating gnuplot script for " << plotname.c_str() << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    stringstream finss;
    aurostd::StringstreamClean(finss);
    finss << "#Generated by AFLOW (Cormac Toher [cormac.toher@duke.edu], 2013, Duke)" << endl;
    finss << "set term postscript eps enhanced color font \"Times-Roman, 40\" size 18, 10.125" << endl;
    finss << "set output " << "\"" << plotfilename.c_str() << ".eps" << "\"" << endl;
    finss << "set label '" << AFLOWLIB_CONSORTIUM_STRING << "' at screen 0.75, 0.02 font \"Times-Roman, 32\"" << endl;
    finss << endl;

    finss << "# " << plotname.c_str() << " PLOT" << endl;
    finss << "set title '" << plotname << "\t" << plot_title << "'" << endl;
    finss << "set xtics " << endl;
    finss << "set ytics" << endl;
    finss << "set yrange [" << ymin << ":" << ymax << "]" << endl;
    finss << "set xrange [" << xmin << ":" << xmax << "]" << endl;
    finss << endl;

    finss << "set xlabel '" << xlabel << "' offset graph 0.00" << endl;
    finss << "set ylabel '" << ylabel << "' offset graph 0.00" << endl;
    finss << "set arrow from 0, 0 to graph 1, first 0 nohead lt 3 lw 1.5" << endl;
    finss << "set arrow from 0, 0 to graph 0, first 0 nohead lt 3 lw 1.5" << endl;
    finss << "set key font \"Times-Roman, 40\"" << endl;       
    finss << endl;
    finss << "plot[][] \\" << endl;
    finss << "\"" << datafilename << "\""  << endl;

    finss << endl;
  
    if(!aurostd::stringstream2file(finss, gnuplotscript, "WRITE")) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Unable to open file " << gnuplotscript.c_str() <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }	
    // OBSOLETE aurostd::StringstreamClean(finss);
    finss.clear();
    finss.str(std::string());

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Plotting data for " << plotname.c_str() << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Call gnuplot to plot the data	
    aurostd::execute(XHOST.command("gnuplot")+" " + gnuplotscript);
 
    // Convert the plot file to the format requested by the user
    aurostd::StringstreamClean(aus);
    // OBSOLETE if(XHOST.vflag_control.isflag("KEEP::PNG")) {
    if(XHOST.vflag_control.flag("KEEP::PNG")) {
      aus << XHOST.command("convert") << " " << plotfilename << ".eps" << " " << plotfilename << ".png" << endl;
      // OBSOLETE }  else if(XHOST.vflag_control.isflag("KEEP::PDF")) {
    }  else if(XHOST.vflag_control.flag("KEEP::PDF")) {
      aus << XHOST.command("convert") << " " << plotfilename << ".eps" << " " << plotfilename << ".pdf" << endl;
      // OBSOLETE }  else if(XHOST.vflag_control.isflag("KEEP::JPG")) 
    }  else if(XHOST.vflag_control.flag("KEEP::JPG")) {
      aus << XHOST.command("convert") << " " << plotfilename << ".eps" << " " << plotfilename << ".jpg" << endl;
      // OBSOLETE }  else if(XHOST.vflag_control.isflag("KEEP::GIF")) {
    }  else if(XHOST.vflag_control.flag("KEEP::GIF")) {
      aus << XHOST.command("convert") << " " << plotfilename << ".eps" << " " << plotfilename << ".gif" << endl;
    }
    aurostd::execute(aus);
    aurostd::StringstreamClean(aus);
	
    // Postprocess clean up
    if(!XHOST.vflag_control.flag("KEEP::GPL")) {
      aurostd::RemoveFile(datafilename);
      aurostd::RemoveFile(gnuplotscript);
    }
    if(!XHOST.vflag_control.flag("KEEP::EPS")) {
      aurostd::RemoveFile(plotfilenameeps);
    }

    return 0;
  } 
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::gibbsinpwrite
// ***************************************************************************
namespace AGL_functions {
  // gibbsinpwrite: writes input file for the original Fortran version of GIBBS
  // Useful for testing and comparison
  // Activated by setting [AFLOW_GIBBS]GIBBSWRITEINPUT=ON in the _AFLOWIN_ file
  uint gibbsinpwrite(_AGL_data& AGL_data, ofstream& FileMESSAGE) {
    ostringstream aus;
    // OBSOLETE double amu2au = physconstamu/physconstme;
    // OBSOLETE double cmamu = AGL_data.cellmass / amu2au;
    string filename, inpname;
    vector<double> volumeinputBohr, energyinputHartree;
    inpname = aurostd::CleanStringASCII(AGL_data.sysname);
    inpname = aurostd::RemoveWhiteSpaces(inpname);
    string outputfilename = AGL_data.dirpathname + "/" + AGL_data.sysname + ".out";
    filename = AGL_data.dirpathname + "/" + AGL_data.sysname + ".inp";
    stringstream outfiless;
    aurostd::StringstreamClean(outfiless);

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Preparing input file for original, non-AFLOW version of GIBBS" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Converts units of energy values calculated by VASP from eV to Hartree
    // Converts units of volume from cubic Angstrom to cubic Bohr
    for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
      volumeinputBohr.push_back(AGL_data.volumeinput.at(i) * pow(angstrom2bohr, 3.0));
      energyinputHartree.push_back(AGL_data.energyinput.at(i) / hart2ev);	
    }
 
    outfiless << AGL_data.sysname << endl;
    outfiless << outputfilename << endl;
    outfiless << AGL_data.natoms << endl;
    // OBSOLETE outfiless << cmamu << endl;
    outfiless << AGL_data.cellmass << endl;
    outfiless << AGL_data.energy_infinity << endl;
    outfiless << AGL_data.i_eqn_of_state << endl;
    outfiless << AGL_data.i_debye << "\t" << AGL_data.poissonratio << endl;
    outfiless << AGL_data.pressure_external.size() << "\t";
    for (uint i = 0; i < AGL_data.pressure_external.size(); i++) {
      outfiless << AGL_data.pressure_external.at(i) << "\t";
    }
    outfiless << endl;
    outfiless << AGL_data.temperature_external.size() << "\t";
    for (uint i = 0; i < AGL_data.temperature_external.size(); i++) {
      outfiless << AGL_data.temperature_external.at(i) << "\t";
    }
    outfiless << endl;
    // OBSOLETE outfiless << AGL_data.volumeinput.size() << endl;
    outfiless << volumeinputBohr.size() << endl;
    if(AGL_data.i_debye == 1) {
      // OBSOLETE for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
      for (uint i = 0; i < volumeinputBohr.size(); i++) {
	// OBSOLETE outfiless << AGL_data.volumeinput.at(i) << "\t" << AGL_data.energyinput.at(i) << "\t" << AGL_data.tdebye.at(i) << endl;
	outfiless << volumeinputBohr.at(i) << "\t" << energyinputHartree.at(i) << "\t" << AGL_data.tdebye.at(i) << endl;
      }
    } else {
      // OBSOLETE for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
      for (uint i = 0; i < volumeinputBohr.size(); i++) {
	// OBSOLETE outfiless << AGL_data.volumeinput.at(i) << "\t" << AGL_data.energyinput.at(i) << endl;
	outfiless << volumeinputBohr.at(i) << "\t" << energyinputHartree.at(i) << endl;
      }
    }

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Writing file" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    if(!aurostd::stringstream2file(outfiless, filename, "WRITE")) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Unable to open file " << filename.c_str() <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }	
    aurostd::StringstreamClean(outfiless);

    return 0;
  }
} // namespace AGL_functions

// **************************************************************************
//  End of AFLOW AGL plot and write
// **************************************************************************

#endif // _AFLOW_AGL_PLOT_WRITE_CPP
