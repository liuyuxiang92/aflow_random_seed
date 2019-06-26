// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
// ***************************************************************************
// 
// A namespace for functions to plot electronic structures, phonon properties,
// and physical properties from AFLOW calculations. The default plotting
// engine is gnuplot with the epslatex terminal.
//
// Each plot type such as a DOS plot has a main function that takes the plot
// options and a stringstringstream to output the plots into the desired
// format (e.g. a gnuplot script). Each format should have its own wrapper
// function (only overload when using gnuplot as the output format).
// 
// Plot options should be handled with xoptions only to allow for as much
// flexibility and customizability as possible.

#include "aflow.h"

using std::deque;
using std::string;
using std::stringstream;
using std::vector;
using aurostd::xoption;

static const string BANDDOS_SIZE = "8, 4.5";

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                              PLOT FUNCTIONS                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace plotter {

// Plot options --------------------------------------------------------------

//getPlotOptions//////////////////////////////////////////////////////////////
// Sets plot options for all plots based on command line arguments.
xoption getPlotOptions(const aurostd::xoption& xopt, const string& key, bool datasets) {
  xoption plotoptions;

  // Get options
  vector<string> tokens;
  string scheme = xopt.getattachedscheme(key);
  uint ntokens = aurostd::string2tokens(scheme, tokens, ",", true); // Keep empty arguments

  if (ntokens >= 1) {
    if (tokens[0][tokens[0].size() - 1] != '/') tokens[0] += '/';
    plotoptions.push_attached("DIRECTORY", tokens[0]);
  }
  // Need to shift options if partial DOS are plotted
  uint shift = 0;
  if (datasets) { // DATASET = -1: all data sets
    if ((ntokens >= 2) && !tokens[1].empty()) {
      plotoptions.push_attached("DATASET", tokens[1]);
    } else {
      plotoptions.push_attached("DATASET", "-1");
    }
    shift = 1;
  } else {  // DATASET = 0: no data set specified 
    plotoptions.push_attached("DATASET", "0");
  }
  if (ntokens >= 2 + shift) plotoptions.push_attached("XMIN", tokens[1 + shift]);
  if (ntokens >= 3 + shift) plotoptions.push_attached("XMAX", tokens[2 + shift]);
  if (ntokens >= 4 + shift) plotoptions.push_attached("YSCALE", tokens[3 + shift]);

  // Get title if present
  scheme = xopt.getattachedscheme("PLOTTER::TITLE");
  if (!scheme.empty()) plotoptions.push_attached("TITLE", scheme);

  // Get image format
  scheme = xopt.getattachedscheme("PLOTTER::PRINT");
  if (!scheme.empty()) plotoptions.push_attached("IMAGE_FORMAT", aurostd::tolower(scheme));

  // Set standard background color
  plotoptions.push_attached("BACKGROUND_COLOR", "#FFFFFF");

  // Set standard grid options
  plotoptions.push_attached("GRID_COLOR", "#808080");
  plotoptions.push_attached("GRID_WIDTH", "1");
  plotoptions.push_attached("GRID_LINE_TYPE", "0");

  // Set standard size
  plotoptions.push_attached("PLOT_SIZE", "5.333, 3");
  
  plotoptions.flag("NOWATERMARK", xopt.flag("PLOTTER::NOWATERMARK"));
  return plotoptions;
}

// Electronic structure plots ------------------------------------------------

//getPlotOptionsEStructure////////////////////////////////////////////////////
// Sets the plot options that are specific to electronic structure plots.
xoption getPlotOptionsEStructure(const aurostd::xoption& xopt, const string& key, bool datasets) {
  xoption plotoptions = getPlotOptions(xopt, key, datasets);

  // Projection
  string scheme = xopt.getattachedscheme("PLOTTER::PROJECTION");
  if (scheme.empty()) {
    plotoptions.push_attached("PROJECTION", "ORBITALS");
  } else {
    plotoptions.push_attached("PROJECTION", aurostd::toupper(scheme));
  }

  // No border
  plotoptions.flag("NOBORDER", true);

  // Set to true to not set Efermi to zero
  plotoptions.flag("NOSHIFT", xopt.flag("PLOTTER::NOSHIFT"));

  // Set gray background
  plotoptions.pop_attached("BACKGROUND_COLOR");
  plotoptions.push_attached("BACKGROUND_COLOR", "#E4E4E4");

  // Change grid options
  plotoptions.pop_attached("GRID_COLOR");
  plotoptions.push_attached("GRID_COLOR", "#FFFFFF");
  plotoptions.pop_attached("GRID_WIDTH");
  plotoptions.push_attached("GRID_WIDTH", "2");
  plotoptions.pop_attached("GRID_LINE_TYPE");
  plotoptions.push_attached("GRID_LINE_TYPE", "1");

  return plotoptions;
}

//getPlotOptionsPhonons///////////////////////////////////////////////////////
// Sets the plot options that are specific to phonon dispersions and DOS.
xoption getPlotOptionsPhonons(const aurostd::xoption& xopt, const string& key) {
  xoption plotoptions = getPlotOptionsEStructure(xopt, key);
  string scheme = xopt.getattachedscheme("PLOTTER::UNIT");
  if (scheme.empty()) {
    plotoptions.push_attached("UNIT", "THZ");
  } else {
    plotoptions.push_attached("UNIT", aurostd::toupper(scheme));
  }
  // There is no Fermi level for phonons, so do not shift
  plotoptions.flag("NOSHIFT", true);

  // Orbital-projections do not exist for phonons
  scheme = plotoptions.getattachedscheme("PROJECTION");
  if ((scheme == "ORBITALS") || (scheme == "LM")) {
    plotoptions.pop_attached("PROJECTION");
    plotoptions.push_attached("PROJECTION", "NONE");
  }
  return plotoptions;
}

// Plot functions ------------------------------------------------------------

//generateHeader//////////////////////////////////////////////////////////////
// Creates the header in the desired output format.
void generateHeader(stringstream& out, const aurostd::xoption& plotoptions, bool multiplot) {
  string plottitle = plotoptions.getattachedscheme("PLOT_TITLE");
  string outformat = plotoptions.getattachedscheme("OUTPUT_FORMAT");
  if (outformat == "GNUPLOT") {
    out << "# Generated by AFLOW" << std::endl;
    out << "set terminal epslatex standalone color"
        << " size " << plotoptions.getattachedscheme("PLOT_SIZE") << " linewidth 2" << std::endl;
    out << "set output " << "'" << plotoptions.getattachedscheme("FILE_NAME_LATEX") << ".tex'" << std::endl;
    if (!plottitle.empty())
      out << "set " << (multiplot?"multiplot ":"")
          << "title '" << plottitle << "' offset 0, -0.5" << std::endl;
    if (plotoptions.flag("NOBORDER")) out << "unset border" << std::endl;
    out << "set object 1 rectangle from graph 0,0 to graph 1,1 fc"
        << " rgb '" << plotoptions.getattachedscheme("BACKGROUND_COLOR") << "' behind fs noborder" << std::endl;
    out << "set grid back lt " << plotoptions.getattachedscheme("GRID_LINE_TYPE")
        << " lc rgb '" << plotoptions.getattachedscheme("GRID_COLOR") << "'"
        << " lw " << plotoptions.getattachedscheme("GRID_WIDTH") << std::endl;
    if (!plotoptions.flag("NOWATERMARK")) {
      out << "set label right '"
          << (plotoptions.flag("BANDDOS")?"":"\\scriptsize") << " "
          << AFLOWLIB_CONSORTIUM_STRING << "' at screen 0.98, 0.025" << std::endl;
    }
    out << std::endl;
  }
}

//savePlotGNUPLOT/////////////////////////////////////////////////////////////
// Executes the gnuplot script and converts into the desired image format.
void savePlotGNUPLOT(const xoption& plotoptions, const stringstream& gpfile) {
  string directory = plotoptions.getattachedscheme("DIRECTORY");
  string filename = plotoptions.getattachedscheme("FILE_NAME");
  string filename_latex = plotoptions.getattachedscheme("FILE_NAME_LATEX");
  // PDF is default since we use pdflatex to compile
  string format = plotoptions.getattachedscheme("IMAGE_FORMAT");
  if (format.empty()) format = "pdf";
  string current_dir = aurostd::execute2string("pwd");
  // Create temp directory
  string tmp = aurostd::TmpDirectoryCreate("plotLATEX") + "/";
  chdir(tmp.c_str());
  // Execute gnuplot and pdflatex
  aurostd::stringstream2file(gpfile, filename + ".plt");
  aurostd::execute(XHOST.command("gnuplot") + " " + filename + ".plt");
  aurostd::execute(XHOST.command("pdflatex") + " -interaction=nonstopmode -halt-on-error " + filename_latex + ".tex > /dev/null");
  // Convert to the desired format if not pdf
  if (format != "pdf") {
    aurostd::execute(XHOST.command("convert") + " -quiet -density 300 -background white " + filename_latex + ".pdf " + filename_latex  + "." + format);
  }
  chdir(current_dir.c_str());
  aurostd::CopyFile(tmp + filename_latex + "." + format, directory + filename + "." + format);
  // Keep gnuplot file if aflow was called with --keep=gpl
  if (XHOST.vflag_control.flag("KEEP::GPL")) {
    aurostd::CopyFile(tmp + filename + ".plt", directory);
  }
  // Clean up
  aurostd::RemoveDirectory(tmp);
  if (!aurostd::FileExist(directory + filename + "." + format)) {
    string function = "plotter::savePlotGNUPLOT()";
    string message = "Error while generating plot.";
    throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
  }
}

//setFileName/////////////////////////////////////////////////////////////////
// Sets the file name of the final plot. FILE_NAME_LATEX is the name of the
// tex file that is generated by gnuplot, which has different limitations than
// the output image.
void setFileName(xoption& plotoptions, string filename) {
  if (filename.empty()) {
    string default_title = plotoptions.getattachedscheme("DEFAULT_TITLE");
    filename = default_title;
    // Get filename
    string ext = plotoptions.getattachedscheme("EXTENSION");
    if (!ext.empty()) {
      if (filename.empty()) filename = ext;
      else filename += "_" + ext;
    }
    filename = aurostd::StringSubst(filename, " ", "_");
    string set = plotoptions.getattachedscheme("DATASET");
    if (aurostd::string2utype<int>(set) > 0) {
      filename += "_" + plotoptions.getattachedscheme("DATALABEL");
      filename += "_" + set;
    }
  }
  plotoptions.push_attached("FILE_NAME", filename);
  // The .tex file created by gnuplot cannot have . or includegraphics will break
  plotoptions.push_attached("FILE_NAME_LATEX", aurostd::StringSubst(filename, ".", "_"));
}

//setTitle////////////////////////////////////////////////////////////////////
// Sets the plot title.
void setTitle(xoption& plotoptions) {
  string title = plotoptions.getattachedscheme("TITLE");
  // Format title
  if (title.empty()) title = formatDefaultPlotTitle(plotoptions);
  plotoptions.push_attached("PLOT_TITLE", title);
}

//formatDefaultPlotTitle//////////////////////////////////////////////////////
// Checks if the default title is in a known AFLOW format and formats it
// appropriately.
string formatDefaultPlotTitle(const xoption& plotoptions) {
  string default_title = plotoptions.getattachedscheme("DEFAULT_TITLE");
  if (default_title.empty()) return default_title;
  string title;
  if (aurostd::substring2bool(default_title, "_ICSD_")) {  // Check if AFLOW ICSD format
    vector<string> tokens;
    aurostd::string2tokens(default_title, tokens, "_");
    if (tokens.size() == 3) {
      title = formatCompoundLATEX(tokens[0]) + " (ICSD \\#" + tokens[2];
      string lattice = plotoptions.getattachedscheme("LATTICE");
      if (lattice.empty()) title += ")";
      else title += ", " + lattice + ")";
    } else { // Title not in ICSD format
      return default_title;
    }
  } else if (aurostd::substring2bool(default_title, ".")) {  // Check if AFLOW prototype format
    vector<string> tokens;
    aurostd::string2tokens(default_title, tokens, ".");
    if (tokens.size() == 2) {
      tokens[1] = aurostd::fixStringLatex(tokens[1], false, false);  // Compound should not habe special characters after formatting
      title = formatCompoundLATEX(tokens[0]) + " (" + tokens[1];
      string lattice = plotoptions.getattachedscheme("LATTICE");
      if (lattice.empty()) title += ")";
      else title += ", " + lattice + ")";
    } else {  // Title not in proptype format
      return default_title;
    }
  } else {  // Not an AFLOW-formatted default
    return default_title;
  }
  // Code only gets here if the title is AFLOW-formatted
  string set = plotoptions.getattachedscheme("DATASET");
  if (aurostd::string2utype<int>(set) > 0) {
    title += " " + plotoptions.getattachedscheme("SETLABEL");
    title += " " + plotoptions.getattachedscheme("DATALABEL");
    title += " (" + set + ")";
  }
  return title;
}

string formatCompoundLATEX(string compound) {
  // First, remove any kind of pseudopotential information
  if (aurostd::substring2bool(compound, ":")) {
    string::size_type t;
    t = compound.find_first_of(':');
    compound = compound.substr(0, t);
  }
  if (aurostd::substring2bool(compound, "_")) {
    string pps = "_d,_h,_pv,_s,_sv,_3,_2,_2_n";
    vector<string> tokens;
    aurostd::string2tokens(pps, tokens, ",");
    for (uint i = 0; i < tokens.size(); i++) {
      compound = aurostd::RemoveSubString(compound, tokens[i]);
    }
  }

  string latex;
  bool sub = false;  // $_{ open?
  for (uint i = 0; i < compound.size(); i++) {
    if (isdigit(compound[i]) && !sub) {
      latex += "$_{";
      sub = true;
    } else if (isalpha(compound[i]) && sub) {
      latex += "}$";
      sub = false;
    }
    latex += compound[i];
  }
  // Close possibly open curly bracket
  if (sub) latex += "}$";
  // Remove sub-1 from chemical formula
  latex = aurostd::RemoveSubString(latex, "$_{1}$");
  return latex;
}

} // namespace plotter

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                           ELECTRONIC STRUCTURE                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

static const string EFERMI_COLOR = "#0000FF";
static const int ESTRUCTURE_NCOLORS = 11;
static const string ESTRUCTURE_COLORS[ESTRUCTURE_NCOLORS] = {
    "#000000",  // black
    "#4C72B0",  // blue
    "#55A868",  // green
    "#C44E52",  // red
    "#CCB974",  // yellow
    "#8172B2",  // purple
    "#64B5CD",  // light blue
    "#E08000",  // orange
    "#006060",  // blue-green
    "#A06000",  // brown
    "#BE80FF"   // light purple
};
static const string ISPIN_COLORS[2] = {"#000000", "#C44E52"};
static const string ORBITALS[4] = {"s", "p", "d", "f"};
static const string LM_ORBITALS[16] = {"s", "p_y", "p_z", "p_x",
                                       "d_{xy}", "d_{yz}", "d_{z^2}", "d_{xz}", "d_{x^2-y^2}",
                                       "f_{y(3x^2-y^2)}", "f_{xyz}", "f_{yz^2}", "f_{z^3}",
                                       "f_{xz^2}", "f_{z(x^2-y^2)}", "f_{x(x^2-3y^2)}"};

namespace plotter {

// Plot functions -------------------------------------------------------------

//PLOT_DOS////////////////////////////////////////////////////////////////////
// PLots electronic densities of states.
void PLOT_DOS(xoption& plotoptions) {
  // Set output format to gnuplot
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

  stringstream out;
  PLOT_DOS(plotoptions, out);
  savePlotGNUPLOT(plotoptions, out);
}

void PLOT_DOS(xoption& plotoptions, stringstream& out) {
  plotoptions.push_attached("EXTENSION", "dos");
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
  // Read files
  string directory = plotoptions.getattachedscheme("DIRECTORY");
  xDOSCAR xdos;
  xdos.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "DOSCAR"));
  // Make sure the projections are consistent with the DOSCAR file
  if ((plotoptions.getattachedscheme("PROJECTION") == "LM") && !(xdos.lmResolved)) {
    std::cerr << "Found --projection=lm, but DOSCAR is not lm-resolved."
              << " Will choose --projection=orbitals instead." << std::endl;
    plotoptions.pop_attached("PROJECTION");
    plotoptions.push_attached("PROJECTION", "ORBITALS");
  }

  plotoptions.push_attached("DEFAULT_TITLE", xdos.title);
  setFileName(plotoptions);
  setTitle(plotoptions);

  plotoptions.flag("BANDDOS", false);
  // Set Fermi energy to zero
  if (plotoptions.flag("NOSHIFT")) {
    plotoptions.push_attached("EFERMI", aurostd::utype2string<double>(xdos.Efermi));
  } else {
    plotoptions.push_attached("EFERMI", "0.0");
  }

  // Get Emin and Emax
  setEMinMax(plotoptions, xdos.energy_min, xdos.energy_max);

  // Plot
  generateHeader(out, plotoptions, false);
  generateDosPlot(out, xdos, plotoptions);
}

//PLOT_PDOS///////////////////////////////////////////////////////////////////
// Plots projected density of states. If PDOS == -1, the projected DOS of all
// atoms will be plotted into separate files.
void PLOT_PDOS(xoption& plotoptions) {
  // Set output format to gnuplot
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

  stringstream out;
  PLOT_PDOS(plotoptions, out);
}

void PLOT_PDOS(xoption& plotoptions, stringstream& out) {
  string set = plotoptions.getattachedscheme("DATASET");
  string outformat = plotoptions.getattachedscheme("OUTPUT_FORMAT");
  int pdos = -1;
  if (!set.empty()) pdos = aurostd::string2utype<int>(set);
  plotoptions.push_attached("SETLABEL", "PDOS");
  if (pdos == 0) {  // Plot total DOS
    PLOT_DOS(plotoptions);
    if (outformat == "GNUPLOT") savePlotGNUPLOT(plotoptions, out);
  } else {
    xstructure xstr = getStructureWithNames(plotoptions);
    if (pdos == -1) {  // Plot partial DOS of all inequivalent atoms
      pflow::PerformFullSymmetry(xstr);
      int iat;
      for (uint i = 0; i < xstr.iatoms.size(); i++) {
        iat = xstr.iatoms[i][0];
        xstr.atoms[iat].CleanName();
        plotoptions.pop_attached("DATASET");
        plotoptions.pop_attached("DATALABEL");
        plotoptions.push_attached("DATASET", aurostd::utype2string<int>(iat + 1));
        plotoptions.push_attached("DATALABEL", xstr.atoms[iat].cleanname);
        PLOT_DOS(plotoptions, out);
        if (outformat == "GNUPLOT") {
          savePlotGNUPLOT(plotoptions, out);
          out.str(string());
          out.clear();
        }
      }
    } else {
      xstr.atoms[pdos - 1].CleanName();
      plotoptions.push_attached("DATALABEL", xstr.atoms[pdos - 1].cleanname);
      PLOT_DOS(plotoptions, out);
      if (outformat == "GNUPLOT") savePlotGNUPLOT(plotoptions, out);
    }
  }
}

//PLOT_BAND///////////////////////////////////////////////////////////////////
// Plots band structures.
void PLOT_BAND(xoption& plotoptions) {
  // Set k-points format to LaTeX
  plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
  // Set output format to gnuplot
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

  stringstream out;
  PLOT_BAND(plotoptions, out);
  savePlotGNUPLOT(plotoptions, out);
}

void PLOT_BAND(xoption& plotoptions, stringstream& out) {
  plotoptions.push_attached("EXTENSION", "band");
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
  // Read files
  string directory = plotoptions.getattachedscheme("DIRECTORY");
  xDOSCAR xdos;
  xdos.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "DOSCAR"));
  xEIGENVAL xeigen;
  xeigen.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "EIGENVAL"));
  xKPOINTS xkpts;
  xkpts.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "KPOINTS"));
  std::stringstream poscar;
  aflowlib::vaspfile2stringstream(directory, "POSCAR", poscar);
  xstructure xstr(poscar);

  plotoptions.push_attached("DEFAULT_TITLE", xeigen.title);
  plotoptions.push_attached("LATTICE", getLatticeFromKpointsTitle(xkpts.title));
  setFileName(plotoptions);
  setTitle(plotoptions);

  plotoptions.flag("BANDDOS", false);
  // Set Fermi energy to zero
  if (plotoptions.flag("NOSHIFT")) {
    plotoptions.push_attached("EFERMI", aurostd::utype2string<double>(xdos.Efermi));
  } else {
    shiftEfermiToZero(xeigen, xdos.Efermi);
    plotoptions.push_attached("EFERMI", "0.0");
  }

  // Get Emin and Emax
  setEMinMax(plotoptions, xeigen.energy_min, xeigen.energy_max);

  // Plot
  generateHeader(out, plotoptions, false);
  generateBandPlot(out, xeigen, xkpts, xstr, plotoptions);
}

//PLOT_BANDDOS////////////////////////////////////////////////////////////////
// Plots combined band structure + DOS plots.
void PLOT_BANDDOS(xoption& plotoptions) {
  // Set k-points format to LaTeX
  plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
  // Set output format to gnuplot
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

  stringstream out;
  PLOT_BANDDOS(plotoptions, out);
  savePlotGNUPLOT(plotoptions, out);
}

void PLOT_BANDDOS(xoption& plotoptions, stringstream& out) {
  plotoptions.push_attached("EXTENSION", "banddos");
  // Increase plot size
  plotoptions.pop_attached("PLOT_SIZE");
  plotoptions.push_attached("PLOT_SIZE", BANDDOS_SIZE);

  // Read files
  string directory = plotoptions.getattachedscheme("DIRECTORY");
  xDOSCAR xdos;
  xdos.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "DOSCAR"));
  xEIGENVAL xeigen;
  xeigen.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "EIGENVAL"));
  xKPOINTS xkpts;
  xkpts.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "KPOINTS"));
  std::stringstream poscar;
  aflowlib::vaspfile2stringstream(directory, "POSCAR", poscar);
  xstructure xstr(poscar);

  plotoptions.push_attached("DEFAULT_TITLE", xeigen.title);
  plotoptions.push_attached("LATTICE", getLatticeFromKpointsTitle(xkpts.title));
  setFileName(plotoptions);
  setTitle(plotoptions);

  plotoptions.flag("BANDDOS", true);
  // Set Fermi energy to zero
  if (plotoptions.flag("NOSHIFT")) {
    plotoptions.push_attached("EFERMI", aurostd::utype2string<double>(xdos.Efermi));
  } else {
    shiftEfermiToZero(xeigen, xdos.Efermi);
    plotoptions.push_attached("EFERMI", "0.0");
  }

  // Get Emin and Emax
  setEMinMax(plotoptions, xdos.energy_min, xdos.energy_max);

  // Plot
  generateHeader(out, plotoptions, true);
  generateBandPlot(out, xeigen, xkpts, xstr, plotoptions);
  generateDosPlot(out, xdos, plotoptions);
}

// Helper functions ----------------------------------------------------------

//getStructureWithNames///////////////////////////////////////////////////////
// Extracts the structure from VASP input files, including species names.
xstructure getStructureWithNames(const xoption& plotoptions) {
  string directory = plotoptions.getattachedscheme("DIRECTORY");
  std::stringstream poscar;
  aflowlib::vaspfile2stringstream(directory, "POSCAR", poscar);
  xstructure xstr(poscar);
  if (xstr.is_vasp4_poscar_format) {
    vector<string> atoms = KBIN::ExtractAtomicSpecies(directory);
    for (uint i = 0; i < atoms.size(); i++) {
      xstr.atoms[i].name = atoms[i];
    }
  }
  return xstr;
}

//getLatticeFromKpointsTitle//////////////////////////////////////////////////
// If the KPOINTS file is formatted according to the AFLOW standard, return
// the lattice of the system. Otherwise, return nothing.
string getLatticeFromKpointsTitle(const string& title) {
  vector<string> tokens;
  aurostd::string2tokens(title, tokens, " ");
  if (tokens.size() >= 3) {  // AFLOW-formatted KPOINTS titles have at least three columns
    if (tokens[0] == "CUB") return tokens[0];
    else if (tokens[0] == "FCC") return tokens[0];
    else if (tokens[0] == "BCC") return tokens[0];
    else if (tokens[0] == "TET") return tokens[0];
    else if (tokens[0] == "BCT") return tokens[0];
    else if (tokens[0] == "BCT1") return tokens[0];
    else if (tokens[0] == "BCT2") return tokens[0];
    else if (tokens[0] == "ORC") return tokens[0];
    else if (tokens[0] == "ORCF") return tokens[0];
    else if (tokens[0] == "ORCF1") return tokens[0];
    else if (tokens[0] == "ORCF2") return tokens[0];
    else if (tokens[0] == "ORCF3") return tokens[0];
    else if (tokens[0] == "ORCI") return tokens[0];
    else if (tokens[0] == "ORCC") return tokens[0];
    else if (tokens[0] == "HEX") return tokens[0];
    else if (tokens[0] == "RHL") return tokens[0];
    else if (tokens[0] == "RHL1") return tokens[0];
    else if (tokens[0] == "RHL2") return tokens[0];
    else if (tokens[0] == "MCL") return tokens[0];
    else if (tokens[0] == "MCLC") return tokens[0];
    else if (tokens[0] == "MCLC1") return tokens[0];
    else if (tokens[0] == "MCLC2") return tokens[0];
    else if (tokens[0] == "MCLC3") return tokens[0];
    else if (tokens[0] == "MCLC4") return tokens[0];
    else if (tokens[0] == "MCLC5") return tokens[0];
    else if (tokens[0] == "TRI") return tokens[0];
    else if (aurostd::toupper(tokens[0]) == "TRI1A") return tokens[0];
    else if (aurostd::toupper(tokens[0]) == "TRI1B") return tokens[0];
    else if (aurostd::toupper(tokens[0]) == "TRI2A") return tokens[0];
    else if (aurostd::toupper(tokens[0]) == "TRI2B") return tokens[0];
    else return "";
  } else {
    return "";
  }
}


//shiftEfermiToZero///////////////////////////////////////////////////////////
// Shift the energies in an xEIGENVAL object so that the Fermi energy is at
// zero. This is not necessary for xDOSCAR because it has a separate vector
// for that purpose.
void shiftEfermiToZero(xEIGENVAL& xeigen, double Efermi) {
  for (uint k = 0; k < xeigen.number_kpoints; k++) {
    for (uint b = 0; b < xeigen.number_bands; b++) {
      for (uint s = 0; s < xeigen.spin + 1; s++) {
        xeigen.venergy[k][b][s] -= Efermi;
      }
    }
  }
}

//setEMinMax//////////////////////////////////////////////////////////////////
// Sets the minimum and maximum energy values for electronic structure plots.
void setEMinMax(xoption& plotoptions, double Emin, double Emax) {
  if (plotoptions.getattachedscheme("XMIN").empty()) {
    if (plotoptions.flag("NOSHIFT")) {
      plotoptions.push_attached("XMIN", aurostd::utype2string<double>(Emin));
    } else {
      plotoptions.push_attached("XMIN", aurostd::utype2string<double>(DEFAULT_DOS_EMIN));
    }
  }
  if (plotoptions.getattachedscheme("XMAX").empty()) {
    if (plotoptions.flag("NOSHIFT")) {
      plotoptions.push_attached("XMAX", aurostd::utype2string<double>(Emax));
    } else {
      plotoptions.push_attached("XMAX", aurostd::utype2string<double>(DEFAULT_DOS_EMAX));
    }
  }
}

// DOS -----------------------------------------------------------------------

//generateDosPlot/////////////////////////////////////////////////////////////
// Generates the data for a DOS plot. 
void generateDosPlot(stringstream& out, const xDOSCAR& xdos, const xoption& plotoptions) {
  deque<deque<deque<double> > > dos;
  vector<string> labels;
  labels.push_back("total");  // There is always a total DOS
  string projection = plotoptions.getattachedscheme("PROJECTION");
  int pdos = aurostd::string2utype<int>(plotoptions.getattachedscheme("DATASET"));
  if (projection == "ORBITALS") {
    // If the DOSCAR is lm-resolved, the orbital projection is the sum of all individual
    // orbitals with the same quantum number
    int norbitals;
    if (xdos.lmResolved) norbitals = (int) std::sqrt(xdos.vDOS[pdos].size());  // size is either 17 or 10
    else norbitals = (int) xdos.vDOS[pdos].size() - 1;
    for (int i = 0; i < norbitals; i++) {
      labels.push_back("$" + ORBITALS[i] + "$");
    }
    if (xdos.lmResolved) {
      // Total DOS and s-orbitals
      dos.push_back(xdos.vDOS[pdos][0]);
      dos.push_back(xdos.vDOS[pdos][1]);
      // p-, d-, and maybe f-orbitals
      deque<deque<deque<double> > > dospart(norbitals - 1, deque<deque<double> >(xdos.spin + 1, deque<double>(xdos.number_energies)));
      for (uint e = 0; e < xdos.number_energies; e++) {
        for (int d = 0; d < norbitals - 1; d++) {
          for (int i = d * (2 + d) + 2; i < d * (d + 4) + 5; i++) {
            for (uint s = 0; s < xdos.spin + 1; s++) {
              dospart[d][s][e] += xdos.vDOS[pdos][i][s][e];
            }
          }
        }
      }
      for (int d = 0; d < norbitals - 1; d++) dos.push_back(dospart[d]);
    } else {
      dos = xdos.vDOS[pdos];
    }
  } else if (projection == "LM") {
    // Safety check
    if (!xdos.lmResolved) {
      string function = "plotter::generateDosPlot()";
      string message = "Projection scheme LM chosen, but DOSCAR is not lm-resolved.";
      throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
    }
    for (uint i = 1; i < xdos.vDOS[pdos].size(); i++) {
      labels.push_back("$" + LM_ORBITALS[i-1] + "$");
    }
    dos = xdos.vDOS[pdos];
  } else if (projection == "ATOMS") {
    xstructure xstr = getStructureWithNames(plotoptions);
    if (pdos == 0) {
      dos.push_back(xdos.vDOS[0][0]);
      if (xdos.vDOS.size() > 1) {
        pflow::PerformFullSymmetry(xstr);
        int iat;
        for (uint i = 0; i < xstr.iatoms.size(); i++) {
          iat = xstr.iatoms[i][0];
          dos.push_back(xdos.vDOS[iat + 1][0]);
          xstr.atoms[iat].CleanName();
          labels.push_back(xstr.atoms[iat].cleanname + "(" + aurostd::utype2string<int>(iat + 1) + ")");
        }
      }
    } else { // In case someone uses pdos option and projection=atoms
      xstr.atoms[pdos - 1].CleanName();
      labels.push_back(xstr.atoms[pdos - 1].cleanname + "(" + aurostd::utype2string<int>(pdos) + ")");
      dos.push_back(xdos.vDOS[pdos][0]);
    }
  } else if (projection == "NONE") {  // Total DOS only without projections
    dos.push_back(xdos.vDOS[0][0]);
  } else {
    string function = "plotter::genertateDosPlot()";
    string message = "Unknown projection scheme " + projection + ".";
    throw aurostd::xerror(function, message, _INPUT_ILLEGAL_);
  }
  string outformat = plotoptions.getattachedscheme("OUTPUT_FORMAT");
  if (outformat == "GNUPLOT") {
    if (plotoptions.flag("NOSHIFT")) {
      generateDosPlotGNUPLOT(out, xdos, xdos.venergy, dos, labels, plotoptions);
    } else {
      generateDosPlotGNUPLOT(out, xdos, xdos.venergyEf, dos, labels, plotoptions);
    }
  }
}

// Bands ---------------------------------------------------------------------

//generateBandPlot////////////////////////////////////////////////////////////
// Generates the data for a band structure plot.
void generateBandPlot(stringstream& out, const xEIGENVAL& xeigen, const xKPOINTS& xkpts,
                      const xstructure& xstr, const xoption& plotoptions) {
  // Create segments
  uint nsegments = xkpts.vpath.size()/2;
  // Make sure that the number of k-points is consistent with EIGENVAL
  if (xeigen.number_kpoints != nsegments * xkpts.path_grid) {
    string function = "plotter::generateBandPlot()";
    string message = "Number of k-points in EIGENVAL and KPOINTS files do not match.";
    throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
  }

  // Labels
  string format = plotoptions.getattachedscheme("KPOINT_FORMAT");
  vector<string> labels(nsegments + 1);
  labels[0] = convertKPointLabel(xkpts.vpath[0], format);
  for (uint i = 2; i < 2 * nsegments; i += 2) {
    labels[i/2] = convertKPointLabel(xkpts.vpath[i - 1], format);
    if (xkpts.vpath[i-1] != xkpts.vpath[i]) {
      labels[i/2] += "|" + convertKPointLabel(xkpts.vpath[i], format);
    }
  }
  labels.back() = convertKPointLabel(xkpts.vpath.back(), format);

  // k-points
  xmatrix<double> f2c = trasp(ReciprocalLattice(xstr.lattice));
  double total_length = 0.0;
  vector<double> segment_points(nsegments + 1, 0.0);
  double dk;
  for (uint i = 0; i < nsegments; i++) {
    dk = aurostd::modulus(f2c * (xkpts.vkpoints[2*i] - xkpts.vkpoints[2*i+1]));
    segment_points[i+1] = segment_points[i] + dk;
    total_length += dk;
  }
  for (uint i = 0; i < nsegments + 1; i++) segment_points[i] /= total_length;

  // Project to path lengths
  vector<double> distances(xeigen.number_kpoints, 0.0);
  int k = 0;
  double ds;
  for (uint i = 0; i < nsegments; i++) {
    ds = (segment_points[i+1] - segment_points[i])/(xkpts.path_grid - 1);
    distances[k] = segment_points[i];
    k++;
    for (int j = 1; j < xkpts.path_grid; j++) {
      distances[k] = distances[k -1] + ds;
      k++;
    }
  }
  string outformat = plotoptions.getattachedscheme("OUTPUT_FORMAT");
  if (outformat == "GNUPLOT") {
    generateBandPlotGNUPLOT(out, xeigen, distances, segment_points, labels, plotoptions);
  }
}

//convertKPointLabel//////////////////////////////////////////////////////////
// Converts a raw k-point string a k-point label into the desired format.
string convertKPointLabel(const string& kpoint, const string& format) {
  vector<string> parts;
  string formatted_label;
  aurostd::string2tokens(kpoint, parts, "_");
  if (parts.size() > 2) return kpoint;
  formatted_label = convertKPointLetter(parts[0], format);
  if (parts.size() == 2) {
    if (format == "LATEX") {
      formatted_label += "$_{" + parts[1] + "}$";
      aurostd::RemoveSubStringFirst(kpoint, "$$");
    } else if (format == "HTML") {
      formatted_label += "<sub>" + parts[1] + "</sub>";
    }
  }
  return formatted_label;
}

//convertKPointLetter/////////////////////////////////////////////////////////
// Converts a raw k-point letter string into the desired format.
string convertKPointLetter(string letter, const string& format) {
  if (format == "LATEX") {
    if (aurostd::substring2bool(letter, "\\") && !aurostd::substring2bool(letter, "\\Gamma")) {
      letter = "$\\mathit{" + letter + "}$";
    } else {
      letter = "$" + letter + "$";
    }
  } else if (format == "HTML") {
    if (aurostd::substring2bool(letter, "\\")) {
      aurostd::StringSubst(letter, "\\", "&");
      letter += ";";
    }
  }
  return letter;
}

// Gnuplot -------------------------------------------------------------------

//generateDosPlotGNUPLOT//////////////////////////////////////////////////////
// Generates the gnuplot script for DOS plots.
void generateDosPlotGNUPLOT(stringstream& out, const xDOSCAR& xdos, const deque<double>& energies,
                            const deque<deque<deque<double> > >& dos, const vector<string>& labels,
                            const xoption& plotoptions) {
  // Initialize variables
  double Efermi = aurostd::string2utype<double>(plotoptions.getattachedscheme("EFERMI"));
  double Emin = aurostd::string2utype<double>(plotoptions.getattachedscheme("XMIN"));
  double Emax = aurostd::string2utype<double>(plotoptions.getattachedscheme("XMAX"));
  bool banddos = plotoptions.flag("BANDDOS");
  bool swap = (banddos || plotoptions.flag("SWAP_AXES"));
  uint ndos = dos.size();

  string yscale, maxdos, mindos;
  yscale = plotoptions.getattachedscheme("YSCALE");
  if (!yscale.empty()) {
    double dosscale = aurostd::string2utype<double>(yscale);
    if (dosscale < _ZERO_TOL_) dosscale = DEFAULT_DOS_SCALE;
    double dosmax = 0.0;
    // Determine the maximum displayed DOS for the xrange
    // Find energy range that is displayed first
    int e1 = -1, e2 = -1;
    for (uint e = 0; e < xdos.number_energies; e++) {
      if ((e1 < 0) && (energies[e] >= Emin)) e1 = (int) e;
      if ((e2 < 0) && (energies[e] > Emax)) e2 = (int) e;
      if ((e1 > 0) && (e2 > 0)) break;
    }
    if (e1 < 0) e1 = 0;  // Emin not found
    if (e2 < 0) e2 = (int) xdos.number_energies; // Emax not found

    for (uint d = 0; d < ndos; d++) {
      for (uint s = 0; s < xdos.spin + 1; s++) {
        for (int e = e1; e < e2; e++) {
          if (dos[d][s][e] > dosmax) dosmax = dos[d][s][e];
        }
      }
    }
    dosmax *= dosscale;
    maxdos = aurostd::utype2string<double>(dosmax);
    mindos = aurostd::utype2string<double>(-(xdos.spin * dosmax));
  }

  string unit = plotoptions.getattachedscheme("UNIT");
  if (unit.empty()) unit = "EV";
  string energyLabel;
  if (aurostd::substring2bool(unit, "EV")) energyLabel = "energy";
  else energyLabel = "frequency";
  unit = getFormattedUnit(unit);

  out << std::endl << "# DOS plot" << std::endl;

  // Create data block
  out << std::endl << "$dos_data << EOD" << std::endl;
  for (uint e = 0; e < xdos.number_energies; e++) {
    out << "  " << energies[e];
    for (uint d = 0; d < ndos; d++) {
      out << " " << dos[d][0][e];
    }
    if (xdos.spin == 1) {
      for (uint d = 0; d < ndos; d++) {
        out << " " << -dos[d][1][e];
      }
    }
    out << std::endl;
  }
  out << "EOD" << std::endl;

  // Margins
  out << "# Margins" << std::endl;
  if (banddos) {
    out << " set lmargin at screen 0.73" << std::endl;
    out << " set rmargin at screen 0.98" << std::endl;
    out << " set tmargin at screen 0.9" << std::endl;
    out << " set bmargin at screen 0.12" << std::endl;
  } else {
    out << " set tmargin at screen 0.9" << std::endl;
    out << " set bmargin at screen 0.2" << std::endl;
  }

  // Key
  out << std::endl << "# Key" << std::endl;
  if (dos.size() == xdos.spin + 1) { // no need for key when only total DOS is plotted
    out << " unset key" << std::endl;
  } else {
    out << " set key samplen 2.5" << std::endl;  // Shorter lines to fit key into image
  }
  
  // Axes
  out << std::endl << "# Axes" << std::endl;
  if (banddos) {
    out << " unset xtics" << std::endl;
  }
  out << " set xtics" << std::endl;
  out << " set ytics" << (banddos?" format \"\"":"") << std::endl;
  out << " set tic scale 0" << std::endl;
  out << " set " << (swap?"y":"x") << "range [" << Emin << ":" << Emax << "]" << std::endl;
  out << " set " << (swap?"x":"y") << "range [" << mindos << ":" << maxdos << "]" << std::endl;
  if (banddos) {
    out << " unset ylabel" << std::endl;
    out << " set title 'DOS (states/" << unit << ")' offset 0,-0.7" << std::endl;
  } else {
    out << " set " << (swap?"y":"x") << "label '" << energyLabel << " (" << unit << ")' offset graph 0.00" << std::endl;
    out << " set " << (swap?"x":"y") << "label 'DOS (states/" << unit << ")' offset graph 0.00" << std::endl;
  }

  // Fermi level
  if (Efermi > Emin) {
    out << std::endl << "# Fermi level" << std::endl;
    if (swap) {
      out << "plot " << Efermi << " lt 1 lc rgb '" << EFERMI_COLOR << "' lw 3 notitle" << std::endl;
    } else {
      out << " set arrow from " << Efermi << ", graph 0 to " << Efermi << ", graph 1"
          << " nohead lt 1 lc rgb '" << EFERMI_COLOR << "' lw 3" << std::endl;
    }
  }

  // Plot data
  out << std::endl << "# Data" << std::endl;
  out << " plot ";
  int xcol, ycol;
  if (swap) ycol = 1;
  else xcol = 1;

  // Majority spin
  for (uint i = 0; i < ndos; i++) {
    if (swap) xcol = i + 2;
    else ycol = i + 2;
    if (i > 0) out << "      ";
    out << "'$dos_data' u " << xcol << ":" << ycol << " w l lt -1 "
        << "lc rgb '" << ESTRUCTURE_COLORS[i % ESTRUCTURE_NCOLORS] << "' lw 2 title '" << labels[i] << "'"
        << (((xdos.spin == 1) || (i < ndos - 1))?",\\":"") << std::endl;
  }
  // Minority spin
  if (xdos.spin == 1) {
    for (uint i = 0; i < ndos; i++) {
      if (swap) xcol = i + ndos + 2;
      else ycol = i + ndos + 2;
      out << "      ";
      out << "'$dos_data' u "  << xcol << ":" << ycol << " w l lt -1 "
          << "lc rgb '" << ESTRUCTURE_COLORS[i % ESTRUCTURE_NCOLORS] << "' lw 2 notitle"  // No title to prevent redundant key entries
          << ((i < ndos - 1)?",\\":"") << std::endl;
    }
  }
}

//generateBandPlotGNUPLOT/////////////////////////////////////////////////////
// Generates the gnuplot script for band structure plots.
void generateBandPlotGNUPLOT(stringstream& out, const xEIGENVAL& xeigen,
                             const vector<double>& xvals, const vector<double>& ticvals,
                             const vector<string>& ticlabels, const xoption& plotoptions) {
  // Initialize variables
  double Efermi = aurostd::string2utype<double>(plotoptions.getattachedscheme("EFERMI"));
  double Emin = aurostd::string2utype<double>(plotoptions.getattachedscheme("XMIN"));
  double Emax = aurostd::string2utype<double>(plotoptions.getattachedscheme("XMAX"));
  bool banddos = plotoptions.flag("BANDDOS");

  string unit = plotoptions.getattachedscheme("UNIT");
  if (unit.empty()) unit = "EV";
  string energyLabel;
  if (aurostd::substring2bool(unit, "EV")) energyLabel = "energy";
  else energyLabel = "frequency";
  unit = getFormattedUnit(unit);

  out << "# Band structure plot" << std::endl;

  uint kpts_per_segment = xeigen.number_kpoints/(ticvals.size() - 1);
  // Create data block
  out << std::endl << "$band_data << EOD" << std::endl;
  for (uint k = 0; k < xeigen.number_kpoints; k++) {
    out << "  " << xvals[k];
    for (uint s = 0; s < xeigen.spin + 1; s++) {
      for (uint b = 0; b < xeigen.number_bands; b++) {
        out << " " << xeigen.venergy[k][b][s];
      }
    }
    out << std::endl;

    // Put each segment into one block. This prevents discontinuities
    // in the band structure from being connected in the plot
    if (((k + 1) % kpts_per_segment == 0) && (k + 1 < xeigen.number_kpoints)) {
      out << std::endl << std::endl;
    }
  }
  out << "EOD" << std::endl << std::endl;

  // Margins
  out << "# Margins" << std::endl;
  if (banddos) {
    out << " set lmargin at screen 0.08" << std::endl;
    out << " set rmargin at screen 0.70" << std::endl;
    out << " set tmargin at screen 0.9" << std::endl;
    out << " set bmargin at screen 0.12" << std::endl;
  } else {
    out << " set tmargin at screen 0.9" << std::endl;
    out << " set bmargin at screen 0.15" << std::endl;
  }

  // Key
  out << std::endl << "# Key" << std::endl;
  out << " unset key" << std::endl;
  
  // Axes
  out << std::endl << "# Axes" << std::endl;
  out << " unset xtics" << std::endl;
  out << " set xtics(";
  uint ntics = ticvals.size();
  for (uint i = 0; i < ntics; i++) {
    out << "'" << ticlabels[i] << "' " << ticvals[i];
    if (i < ntics - 1) out << ", ";
  }
  out << ")" << std::endl;
  out << " set tic scale 0" << std::endl;
  out << " set xrange [0:1]" << std::endl;
  out << " set yrange [" << Emin << ":" << Emax << "]" << std::endl;
  out << " set ylabel '" << energyLabel << " (" << unit << ")'" << std::endl;

  // Fermi level
  if (Efermi > Emin) {
    out << std::endl << "# Fermi level" << std::endl;
    out << " set arrow from 0, " << Efermi << " to graph 1, first " << Efermi
        << " nohead lt 1 lc rgb '" << EFERMI_COLOR << "' lw 3" << std::endl;
  }

  // Plot data
  out << std::endl << "# Data" << std::endl;
  out << " plot ";
  // Majority spin
  for (uint b = 0; b < xeigen.number_bands; b++) {
    if (b > 0) out << "      ";
    out << "'$band_data' u 1:" << (b + 2)
        << " w l lt -1 lc rgb '" <<  ISPIN_COLORS[0] << "' lw 2"
        << (((xeigen.spin == 1) || (b < xeigen.number_bands - 1))?",\\":"") << std::endl;
  }
  
  // Minority spin
  if (xeigen.spin == 1) {
    for (uint b = 0; b < xeigen.number_bands; b++) {
      out << "      ";
      out << "'$band_data' u 1:" << (b + xeigen.number_bands + 2)
          << " w l lt 1 lc rgb '" << ISPIN_COLORS[1] << "' lw 2"
          << ((b < xeigen.number_bands - 1)?",\\":"") << std::endl;
    }
  }
}

//getFormattedUnit////////////////////////////////////////////////////////////
// Formats the energy/frequency unit for band structures. This is especially
// useful for phonons.
string getFormattedUnit(const string& unit) {
  if (unit == "EV") return "eV";
  if (unit == "MEV") return "meV";
  if (unit == "THZ") return "THz";
  if (unit == "HZ") return "Hz";
  if ((unit == "CM-1") || (unit == "RCM")) return "cm$^{-1}$";
  return unit;
}

}  // namespace plotter

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                 PHONONS                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace plotter {

//PLOT_PHDOS//////////////////////////////////////////////////////////////////
// Plots phonon DOS.
void PLOT_PHDOS(xoption& plotoptions) {
  // Set k-points format to LaTeX
  plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
  // Set output format to gnuplot
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

  stringstream out;
  PLOT_PHDOS(plotoptions, out);
  savePlotGNUPLOT(plotoptions, out);
}

void PLOT_PHDOS(xoption& plotoptions, stringstream& out) {
  plotoptions.push_attached("EXTENSION", "phdos");
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
  // Read files
  string directory = plotoptions.getattachedscheme("DIRECTORY");
  xDOSCAR xdos;
  xdos.GetPropertiesFile(directory + DEFAULT_APL_PHDOSCAR_FILE);
  
  plotoptions.push_attached("DEFAULT_TITLE", xdos.title);
  setFileName(plotoptions);
  setTitle(plotoptions);
  
  plotoptions.flag("BANDDOS", false);

  // Convert energies if necessary
  string unit = plotoptions.getattachedscheme("UNIT");
  if (!unit.empty() && (unit != "EV")) {
    convertEnergies(xdos, unit);
  }

  // Get Emin and Emax
  setEMinMax(plotoptions, xdos.energy_min, xdos.energy_max);

  generateHeader(out, plotoptions, false);
  generateDosPlot(out, xdos, plotoptions);
}

//PLOT_PHDISP/////////////////////////////////////////////////////////////////
// Plots phonon dispersion curves.
void PLOT_PHDISP(xoption& plotoptions) {
  // Set k-points format to LaTeX
  plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
  // Set output format to gnuplot
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

  stringstream out;
  PLOT_PHDISP(plotoptions, out);
  savePlotGNUPLOT(plotoptions, out);
}

void PLOT_PHDISP(xoption& plotoptions, stringstream& out) {
  plotoptions.push_attached("EXTENSION", "phdisp");
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
  // Read files
  string directory = plotoptions.getattachedscheme("DIRECTORY");
  xEIGENVAL xeigen;
  xeigen.GetPropertiesFile(directory + DEFAULT_APL_PHEIGENVAL_FILE);
  xKPOINTS xkpts;
  xkpts.GetPropertiesFile(directory + DEFAULT_APL_PHKPOINTS_FILE);
  stringstream poscar;
  aurostd::efile2stringstream(DEFAULT_APL_PHPOSCAR_FILE, poscar);
  xstructure xstr(poscar);

  plotoptions.push_attached("DEFAULT_TITLE", xeigen.title);
  plotoptions.push_attached("LATTICE", getLatticeFromKpointsTitle(xkpts.title));
  setFileName(plotoptions);
  setTitle(plotoptions);

  plotoptions.flag("BANDDOS", false);

  // Convert energies if necessary
  string unit = plotoptions.getattachedscheme("UNIT");
  if (!unit.empty() && (unit != "EV")) {
    convertEnergies(xeigen, unit);
  }

  // Get Emin and Emax
  setEMinMax(plotoptions, xeigen.energy_min, xeigen.energy_max);

  generateHeader(out, plotoptions, false);
  generateBandPlot(out, xeigen, xkpts, xstr, plotoptions);
}

//PLOT_PHDISPDOS//////////////////////////////////////////////////////////////
// Plots combined phonon band structure + DOS plots.
void PLOT_PHDISPDOS(xoption& plotoptions) {
  // Set k-points format to LaTeX
  plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
  // Set output format to gnuplot
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

  stringstream out;
  PLOT_PHDISPDOS(plotoptions, out);
  savePlotGNUPLOT(plotoptions, out);
}

void PLOT_PHDISPDOS(xoption& plotoptions, stringstream& out) {
  plotoptions.push_attached("EXTENSION", "phdispdos");
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
  plotoptions.push_attached("PLOT_SIZE", BANDDOS_SIZE);
  // Read files
  string directory = plotoptions.getattachedscheme("DIRECTORY");
  xDOSCAR xdos;
  xdos.GetPropertiesFile(directory + DEFAULT_APL_PHDOSCAR_FILE);
  xEIGENVAL xeigen;
  xeigen.GetPropertiesFile(directory + DEFAULT_APL_PHEIGENVAL_FILE);
  xKPOINTS xkpts;
  xkpts.GetPropertiesFile(directory + DEFAULT_APL_PHKPOINTS_FILE);
  stringstream poscar;
  aurostd::efile2stringstream(DEFAULT_APL_PHPOSCAR_FILE, poscar);
  xstructure xstr(poscar);

  plotoptions.push_attached("DEFAULT_TITLE", xeigen.title);
  plotoptions.push_attached("LATTICE", getLatticeFromKpointsTitle(xkpts.title));
  setFileName(plotoptions);
  setTitle(plotoptions);

  plotoptions.flag("BANDDOS", true);

  // Convert energies if necessary
  string unit = plotoptions.getattachedscheme("UNIT");
  if (!unit.empty() && (unit != "EV")) {
    convertEnergies(xdos, unit);
    convertEnergies(xeigen, unit);
  }

  // Get Emin and Emax
  setEMinMax(plotoptions, xdos.energy_min, xdos.energy_max);

  generateHeader(out, plotoptions, true);
  generateBandPlot(out, xeigen, xkpts, xstr, plotoptions);
  generateDosPlot(out, xdos, plotoptions);
}

//convertEnergies/////////////////////////////////////////////////////////////
// Converts the energies in an electronic structure object (xDOSCAR/xEIGENVAL)
// into the desired energy/frequency unit.
void convertEnergies(xEIGENVAL& xeigen, const string& unit) {
  double conversion_factor = getEnergyConversionFactor(unit);
  for (uint k = 0; k < xeigen.number_kpoints; k++) {
    for (uint b = 0; b < xeigen.number_bands; b++) {
      for (uint s = 0; s < xeigen.spin + 1; s++) {
        xeigen.venergy[k][b][s] *= conversion_factor;
      }
    }
  }
  xeigen.energy_min *= conversion_factor;
  xeigen.energy_max *= conversion_factor;
}

void convertEnergies(xDOSCAR& xdos, const string& unit) {
  double conversion_factor = getEnergyConversionFactor(unit);
  for (uint i = 0; i < xdos.number_energies; i++) {
    xdos.venergy[i] *= conversion_factor;
    xdos.venergyEf[i] *= conversion_factor;
  }
  xdos.energy_min *= conversion_factor;
  xdos.energy_max *= conversion_factor;
}

//getEnergyConversionFactor///////////////////////////////////////////////////
// Returns the factor to convert eV into the desired energy/frequency unit.
// Supported units are meV, THz, Hz, and reciprocal cm (CM-1/RCM).
double getEnergyConversionFactor(const string& unit) {
  if (unit == "MEV") return 1000.0;
  if (unit == "THZ") return 241.7990504024;
  if (unit == "HZ") return 2.417990504024E14;
  if ((unit == "CM-1") || (unit == "RCM")) return 8065.548153549;
  return 1.0;
}

}  // namespace plotter

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                           PROPERTIES PLOTTERS                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace plotter {

//PLOT_THERMO/////////////////////////////////////////////////////////////////
// Plots APL thermal properties.
void PLOT_THERMO(xoption& plotoptions) {
  stringstream out;
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
  plotoptions.push_attached("COLOR", "$000000");
  plotoptions.push_attached("LINETYPES", "-1");
  PLOT_THERMO(plotoptions, out);
}

void PLOT_THERMO(xoption& plotoptions, stringstream& out) {
  // Set labels
  static const int nprops = 4;
  string ylabels[nprops] = {"U", "F", "S", "c_V"};
  string extensions[nprops] = {"vib_internal_energy", "vib_free_energy", "vib_entropy", "apl_cV"};
  string yunits[nprops] = {"meV/cell", "meV/cell", "$k_B$/cell", "$k_B$/cell"};
  string ymin[nprops] = {"", "", "0", "0"};

  // Get data
  string thermo_file = DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_THERMO_FILE;
  if (aurostd::EFileExist(thermo_file)) {
    string outformat = plotoptions.getattachedscheme("OUTPUT_FORMAT");
    plotoptions.push_attached("DATA_FILE", thermo_file);
    plotoptions.push_attached("KEYWORD", "APL_THERMO");
    vector<vector<double> > data = readAflowDataFile(plotoptions);
    for (int i = 0; i < nprops; i++) {
      plotoptions.pop_attached("YMIN");
      if (!ymin[i].empty()) plotoptions.push_attached("YMIN", ymin[i]);
      plotoptions.push_attached("EXTENSION", extensions[i]);
      setPlotLabels(plotoptions, "T", "K", ylabels[i], yunits[i]);
      plotSingleFromSet(plotoptions, out, data, i + 2);
      if (outformat == "GNUPLOT") {
        savePlotGNUPLOT(plotoptions, out);
      }
    }
  } else {
    string function = "plotter::PLOT_THERMO()";
    string message = "Could not find file " + thermo_file + ".";
    throw aurostd::xerror(function, message, _FILE_NOT_FOUND_);
  }
}

//PLOT_TCOND//////////////////////////////////////////////////////////////////
// Plots AAPL thermal conductivity tensors
void PLOT_TCOND(xoption& plotoptions) {
  stringstream out;
  plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
  plotoptions.flag("LINESPOINTS", true);
  PLOT_TCOND(plotoptions, out);
  savePlotGNUPLOT(plotoptions, out);
}

void PLOT_TCOND(xoption& plotoptions, stringstream& out) {
  plotoptions.push_attached("EXTENSION", "thermal_conductivity");
  string tcond_file = DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_TCOND_FILE;
  if (aurostd::EFileExist(tcond_file)) {
    plotoptions.push_attached("DATA_FILE", tcond_file);
    plotoptions.push_attached("KEYWORD", "AAPL_THERMAL_CONDUCTIVITY");
    plotoptions.flag("CONTRAVARIANT", true);
    plotoptions.push_attached("YMIN", "0");
    setPlotLabels(plotoptions, "T", "K", "\\kappa", "W/m K");
    plotMatrix(plotoptions, out);
  } else {
    string function = "plotter::PLOT_TCOND()";
    string message = "Could not find file " + tcond_file + ".";
    throw aurostd::xerror(function, message, _FILE_NOT_FOUND_);
  }
}

}  // namespace plotter

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                              GENERAL PLOTS                               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace plotter {

// Color palettes ------------------------------------------------------------

// The color palette was designed to be accessible for people with color vision
// deficiencies. When changing the colors, please make sure that they are still
// distinguishable for everyone!
static const string MATRIX_COLORS = "#000000,#004949,#009292,#490092,#B66DFF,#6DB6FF,#924900,#D55E00,#EDB120";

// Point styles --------------------------------------------------------------

static const string MATRIX_POINT_STYLES = "17,35,51,44,18,9,60,11,20";

// Line types ----------------------------------------------------------------

static const string MATRIX_LINE_TYPES = "-1";

// Pre-set labels -----------------------------------------------------------

static const string MATRIX_LABELS[9] = {"xx", "yx", "zx",
                                        "xy", "yy", "zy",
                                        "xz", "yz", "zz"};


//plotSingleFromSet///////////////////////////////////////////////////////////
// Plots a single column from a dataset.
void plotSingleFromSet(xoption& plotoptions, stringstream& out,
                       const vector<vector<double> >& data_set, int col) {
  setFileName(plotoptions);
  setTitle(plotoptions);
  
  vector<vector<double> > data(data_set.size(), vector<double>(2));
  for (uint i = 0; i < data.size(); i++) {
    data[i][0] = data_set[i][0];
    data[i][1] = data_set[i][col];
  }

  generateHeader(out, plotoptions);
  generatePlotGNUPLOT(out, plotoptions, data);
}

//plotMatrix//////////////////////////////////////////////////////////////////
// Plots a 3 x 3 matrix.
void plotMatrix(xoption& plotoptions, stringstream& out) {
  vector<vector<double> > data = readAflowDataFile(plotoptions);
  setFileName(plotoptions);
  setTitle(plotoptions);

  // Set additional plot options
  bool contravariant = plotoptions.flag("CONTRAVARIANT");
  string title_base = plotoptions.getattachedscheme("YLABEL");
  vector<string> titles(9);
  vector<string> colors(9);
  vector<int> point_styles(9);
  for (int i = 0; i < 9; i++) {
    if (contravariant) titles[i] = "$" + title_base + "^{" + MATRIX_LABELS[i] + "}$";
    else titles[i] = "$" + title_base + "_{" + MATRIX_LABELS[i] + "}$";
  }
  plotoptions.push_attached("TITLES", aurostd::joinWDelimiter(titles, ","));
  plotoptions.push_attached("COLORS", MATRIX_COLORS);
  plotoptions.push_attached("POINTSTYLES", MATRIX_POINT_STYLES);
  plotoptions.push_attached("LINETYPES", MATRIX_LINE_TYPES);

  generateHeader(out, plotoptions);
  generatePlotGNUPLOT(out, plotoptions, data);
}

//setPlotLabels///////////////////////////////////////////////////////////////
// Stores the labels and units for the plots.
void setPlotLabels(xoption& plotoptions,
                   const string& xlabel, const string& xunit,
                   const string& ylabel, const string& yunit) {
  plotoptions.pop_attached("XLABEL");
  plotoptions.push_attached("XLABEL", xlabel);
  plotoptions.pop_attached("XUNIT");
  plotoptions.push_attached("XUNIT", xunit);
  plotoptions.pop_attached("YLABEL");
  plotoptions.push_attached("YLABEL", ylabel);
  plotoptions.pop_attached("YUNIT");
  plotoptions.push_attached("YUNIT", yunit);
}

//readDataFile////////////////////////////////////////////////////////////////
// Reads data from an AFLOW data file. Requires a START and STOP string to be
// present so that it can skip headers and other data sets.
vector<vector<double> > readAflowDataFile(xoption& plotoptions) {
  vector<vector<double> > data;
  vector<double> row;
  vector<string> vcontent;
  string keyword = plotoptions.getattachedscheme("KEYWORD");
  string path_to_file = plotoptions.getattachedscheme("DIRECTORY") + plotoptions.getattachedscheme("DATA_FILE");
  string startstring = "[" + keyword + "]START";
  string stopstring = "[" + keyword + "]STOP";
  string systemstring = "[" + keyword + "]SYSTEM=";
  aurostd::efile2vectorstring(path_to_file, vcontent);
  string line;
  uint nlines = vcontent.size();
  for (uint iline = 0; iline < nlines; iline++) {
    line = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vcontent[iline]);
    if (aurostd::substring2bool(line, systemstring)) {
      vector<string> tokens;
      aurostd::string2tokens(line, tokens, "=");
      if (tokens.size() == 2) plotoptions.push_attached("DEFAULT_TITLE", tokens.back());
    }
    if (vcontent[iline] == startstring) {
      iline++;
      while ((vcontent[iline] != stopstring) && (iline < nlines)) {
        line = aurostd::RemoveWhiteSpacesFromTheFront(vcontent[iline]);
        if (line[0] != '#') {
          aurostd::string2tokens(vcontent[iline], row, " ");
          data.push_back(row);
        }
        iline++;
      }
    }
    if (vcontent[iline] == stopstring) break;
    if (iline == nlines) {
      string function = "plotter::readAflowDataFile()";
      string message = "Wrong file format. No STOP tag found.";
      throw aurostd::xerror(function, message, _FILE_WRONG_FORMAT_);
    }
  }
  if (data.size() == 0) {
    string function = "plotter::readAflowDataFile()";
    string message = "No data extracted from file " + path_to_file + ".";
    message += "File is either empty or has the wrong format.";
    throw aurostd::xerror(function, message, _FILE_WRONG_FORMAT_);
  }
  return data;
}

//generatePlotGNUPLOT/////////////////////////////////////////////////////////
// Generate the gnuplot scripts for general plots.
void generatePlotGNUPLOT(stringstream& out, const xoption& plotoptions,
                         const vector<vector<double> >& data) {
  uint ndata = data[0].size() - 1;

  // Axes settings
  string xmax = plotoptions.getattachedscheme("XMAX");
  string xmin = plotoptions.getattachedscheme("XMIN");
  string ymax = plotoptions.getattachedscheme("YMAX");
  string ymin = plotoptions.getattachedscheme("YMIN");
  string xlabel = plotoptions.getattachedscheme("XLABEL");
  string xunit = plotoptions.getattachedscheme("XUNIT");
  if (!xunit.empty()) xunit = " (" + xunit + ")";
  string ylabel = plotoptions.getattachedscheme("YLABEL");
  string yunit = plotoptions.getattachedscheme("YUNIT");
  if (!yunit.empty()) yunit = " (" + yunit + ")";
  
  // Plot types: lines (default), linespoints, or points only (nolines)
  bool linespoints = plotoptions.flag("LINESPOINTS");
  bool points = linespoints || plotoptions.flag("NOLINES");
  bool lines = (linespoints || !points);
  string plotstyle;
  if (lines) plotstyle += "l";
  if (points) plotstyle += "p";

  // Colors, point styles, line types
  vector<string> colors;
  vector<int> point_styles, line_types;
  aurostd::string2tokens(plotoptions.getattachedscheme("COLORS"), colors, ", ");
  aurostd::string2tokens(plotoptions.getattachedscheme("POINTSTYLES"), point_styles, ", ");
  aurostd::string2tokens(plotoptions.getattachedscheme("LINETYPES"), line_types, ", ");
  // Assume the same color, etc. for all plots if only one is given
  if ((ndata > 1) && (colors.size() == 1)) colors.assign(ndata, colors[0]);
  if ((ndata > 1) && (point_styles.size() == 1)) point_styles.assign(ndata, point_styles[0]);
  if ((ndata > 1) && (line_types.size() == 1)) line_types.assign(ndata, line_types[0]);
  bool colors_set = (colors.size() == ndata);
  bool points_set = ((linespoints || points) && (point_styles.size() == ndata));
  bool lines_set = (lines && (line_types.size() == ndata));

  // Titles
  vector<string> titles;
  aurostd::string2tokens(plotoptions.getattachedscheme("TITLES"), titles, ", ");
  if (titles.size() < ndata) {  // Fill up with empty titles if there aren't enough
    for (uint i = titles.size(); i < ndata; i++) titles.push_back("");
  }

  // Data block
  out << "$matrix_data << EOD" << std::endl;
  for (uint i = 0; i < data.size(); i++) {
    for (uint j = 0; j < ndata + 1; j++) {  // ndata + 1 to include x-values
      out << " " << data[i][j];
    }
    out << std::endl;
  }
  out << "EOD" << std::endl << std::endl;

  // Margins
  out << "# Margins" << std::endl;
  out << " set tmargin at screen 0.9" << std::endl;
  out << " set bmargin at screen 0.22" << std::endl;

  // Axes
  out << "# Axes" << std::endl;
  out << " set xrange [" << xmin << ":" << xmax << "]" << std::endl;
  out << " set yrange [" << ymin << ":" << ymax << "]" << std::endl;
  out << " set xlabel '$" << xlabel << "$" << xunit << "'" << std::endl;
  out << " set ylabel '$" << ylabel << "$" << yunit << "'" << std::endl;
  out << " set tics nomirror out" << std::endl;

  // Key
  out << std::endl << "# Key" << std::endl;
  if (ndata == 1) {  // No need for legend if only one set of data to plot
    out << " unset key" << std::endl;
  } else {
    out << " set key horizontal maxcols 3" << std::endl;
  }

  // Plot
  out << std::endl << "# Plot" << std::endl;
  out << " plot ";
  for (uint i = 0; i < ndata; i++) {
    if (i > 0) out << "      ";
    out << "'$matrix_data' u 1:" << (i + 2) << " w " << plotstyle;
    if (lines) out << " lw 2";
    if (lines_set) out << " lt " << line_types[i];
    if (points) out << " ps 1.5";
    if (colors_set) out << " lc rgb '" << colors[i] << "'";
    if (points_set) out << " pt " << point_styles[i];
    if (titles[i].empty()) {
      out << " notitle";
    } else {
      out << " title '" << titles[i] << "'";
    }
    out << ((i < ndata - 1)?",\\":"") << std::endl;
  }
}

}  // namespace plotter

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
// ***************************************************************************

