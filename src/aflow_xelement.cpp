// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2019
#ifndef _AFLOW_XELEMENT_CPP
#define _AFLOW_XELEMENT_CPP
#include "aflow.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// _XELEMENT_PROTOTYPES_

/* 
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// XELEMENT
// look into aflow.h for the definitions

namespace xelement {
class xelement { // simple class.. nothing fancy
public:
// constructor destructor                              // constructor/destructor
xelement();                                            // default, just allocate
xelement(uint);                                            // default, just allocate
xelement(string);                                      // look at it by symbol or name
~xelement();                                           // kill everything
const xelement& operator=(const xelement &b);          // copy
void clear();
// content                                             // content
bool verbose;
// [AFLOW]START=DECLARATION
int Z; 
string symbol;                          // http://periodictable.com      // DU 2019/05/17
string name;                            // http://periodictable.com      // DU 2019/05/17
double Period;                          // http://periodictable.com      // DU 2019/05/17
double Group;                           // http://periodictable.com      // DU 2019/05/17
string Series;                          // http://periodictable.com For Nh,Fl,Mc,Lv,Ts Value is a guess based on periodic table trend.      // DU 2019/05/17 
string Block;                           // http://periodictable.com      // DU 2019/05/17
//                                          
double mass;                            // (kg)
double MolarVolume;                     // (m^3/mol) http://periodictable.com      // DU 2019/05/17
double volume;                          // atomic volume in A^3 from the FCC vasp table and/or successive calculations
double Miedema_Vm;                      // (V_m^{2/3} in (cm^2)) Miedema Rule Table 1a Physica 100B (1980) 1-28
// for lanthines from J.A. Alonso and N.H. March. Electrons in Metals and Alloys, Academic Press, London (1989) (except La)
double valence_std;                     // http://en.wikipedia.org/wiki/Valence_(chemistry) standard: number electrons minus closed shell at leff (noble gas)
double valence_iupac;                   // http://en.wikipedia.org/wiki/Valence_(chemistry) IUPAC Maximum number of univalent atoms that may combine with an atom of the element under consideration, or with a fragment, or for which an atom of this element can be substituted.
double valence_PT;                      //           http://periodictable.com      // DU 2019/05/17
double Density_PT;                      // (g/cm^3)  http://periodictable.com      // DU 2019/05/17
string crystal;                     // Ashcroft-Mermin                                                                                                                   
string CrystalStr_PT;                   // http://periodictable.com      // DU 2019/05/17
string space_group;                     // http://periodictable.com      // DU 2019/05/17
uint space_group_number;                // http://periodictable.com      // DU 2019/05/17


double Pearson_coefficient;             // Pearson mass deviation coefficient // ME20181020
xvector<double> lattice_constant;       // (pm) http://periodictable.com      // DU 2019/05/17
xvector<double> lattice_angle;          // (rad) http://periodictable.com      // DU 2019/05/17
string phase;                           //      http://periodictable.com      // DU 2019/05/17
double radius;                      // Saxena (nm)
double radius_PT;                       // (pm)       http://periodictable.com      // DU 2019/05/17
double radius_covalent_PT;              // (pm)       http://periodictable.com      // DU 2019/05/17
double radius_covalent;             // (Angstrom) Dalton Trans. 2836, 2832-2838 (2008) // DX and CO -20170904
double radius_VanDerWaals_PT;           // (pm)       http://periodictable.com      // DU 2019/05/17
double radii_Ghosh08;                    // (Angstrom) Journal of Molecular Structure: THEOCHEM 865, 60–67 (2008)      // DU 2019/05/17
double radii_Slatter;                    // (Angstrom) J. of Chem. Phys. 41, 3199 (1964)      // DU 2019/05/17
double radii_Pyykko;                     // (pm) single bond covalent radii  Chem. Eur. J. 15, 186-197 (2009)      // DU 2019/05/17
//                                          
double electrical_conductivity;          // (S/m)  http://periodictable.com  Value given for graphite. Diamond electrical conductivity is approximately 0.001.      // DU 2019/05/17
double electronegativity_vec;           // Saxena
double hardness_Ghosh;                   // (eV) Int. J. Quantum Chem 110, 1206-1213 (2010) Table III       // DU 2019/05/17
double electronegativity_Pearson;                  // (eV) Inorg. Chem., 27(4), 734–740 (1988)      // DU 2019/05/17
double electronegativity_Ghosh;                    // (eV) Journal of Theoretical and Computational Chemistry, 4, 21-33 (2005)      // DU 2019/05/17
double electron_affinity_PT;             // (kJ/mol)  http://periodictable.com       // DU 2019/05/17
double Miedema_phi_star;                // (V)        (phi^\star   Miedema Rule Table 1a Physica 100B 1-28 (1980)
double Miedema_nws;                     // (d.u.)^1/3 n_{ws}^{1/3} Miedema Rule Table 1a Physica 100B 1-28 (1980)
double Miedema_gamma_s;                 // (mJ/m^2)   \gamma_s^0   Miedema Rule Table 1a Physica 100B 1-28 (1980)
double Pettifor_scale;                  // Chemical Scale Pettifor Solid State Communications 51 31-34 (1984)
//                                          
double boiling_point;                   // (Celsius), http://periodictable.com C:diamond, P:"YELLOW" Phosphorus, As:sublimates at this T.      // DU 2019/05/17
double melting_point;                   // (Celsius), http://periodictable.com He does not solidify at standard pressure,C: Value given for diamond form, P : Value given for "YELLOW" phosphorus form, S : Value given for monoclinic, beta form, Se: Value given for hexagonal, gray form, Bk: Value given for alpha form.           // DU 2019/05/17
double vaporization_heat_PT;             // (kJ/mol)   http://periodictable.com      // DU 2019/05/17
double specific_heat_PT;                 // (J/(kg.K)) http://periodictable.com Gas_Phase:H(H2),He,N(N2),O(O2),F(F2),Ne,Cl(Cl2),Ar,Kr,Tc,Xe,Rn,Ra,Pa -- Liquid_Phase:Br,Hg -- Solid Phase: B(rhombic),C(graphite),S(rhombic),P(phase of P.4),As(alpha),Se(hexagonal),Cd(gamma),Sn(gray),Li,In,Be,Na,Mg,Al,Si,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,Rb,Sr,Y,Zr,Nb,Mo,Ru,Rh,Pd,Ag,Sb,Te,I,Cs,Ba,La,Ce,Pr,Nd,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Tl,Pb,Bi,Ac,Th,U.      // DU 2019/05/17 
double critical_Pressure;                // (Atm)      http://periodictable.com Li,Na,K,Rb: Value estimated based on extrapolation.      // DU 2019/05/17
double critical_Temperature_PT;          // (K)        http://periodictable.com Li,Na,K,Rb: Value estimated based on extrapolation.      // DU 2019/05/17
double thermal_expansion;               // (K^{-1})   http://periodictable.com C:graphite      // DU 2019/05/17
double thermal_conductivity;            // (W/(mK))   http://periodictable.com      // DU 2019/05/17
//                                         
double Brinelll_hardness;               // (MPa)  http://periodictable.com For Ge value is converted from Mohs scale      // DU 2019/05/17
double Mohs_hardness;                   //        http://periodictable.com For C, value given for graphite. Diamond value is 10.0; For Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Lu converted from Vickers scale.      // DU 2019/05/17
double Vickers_hardness;                // (MPa)  http://periodictable.com For Si,Ge,As,Ru,Os converted from Brinell scale.      // DU 2019/05/17
double Hardness_Pearson;                // (eV)   Inorg. Chem. 27(4) 734-740 (1988).      // DU 2019/05/17
double Hardness_Putz;                   // (eV/atom) International Journal of Quantum Chemistry, Vol 106, 361–389 (2006), TABLE-V.      // DU 2019/05/17
double Hardness_RB;                     // (eV)   Robles and Bartolotti, J. Am. Chem. Soc. 106, 3723-3727 (1984).      // DU 2019/05/17
double shear_modulus;                    // (GPa)  http://periodictable.com      // DU 2019/05/17
double Young_modulus;                    // (GPa)  http://periodictable.com      // DU 2019/05/17
double bulk_modulus;                     // (GPa)  http://periodictable.com      // DU 2019/05/17
double Poisson_ratio_PT;                 // (--)   http://periodictable.com      // DU 2019/05/17
double Miedema_BVm;                     // (kJ/mole) BV_m Miedema Rule Table 1a Physica 100B 1-28 (1980) 
//
string Magnetic_Type_PT;                 //           http://periodictable.com       // DU 2019/05/17
double Mass_Magnetic_Susceptibility;      // (m^3/K)   http://periodictable.com       // DU 2019/05/17
double Volume_Magnetic_Susceptibility;    //           http://periodictable.com       // DU 2019/05/17
double Molar_Magnetic_Susceptibility;     // (m^3/mol) http://periodictable.com       // DU 2019/05/17
double Curie_point;                     // (K)       http://periodictable.com       // DU 2019/05/17
//
double refractive_index;                 // http://periodictable.com C:diamond      // DU 2019/05/17
string color_PT;                        // http://periodictable.com      // DU 2019/05/17
//
double HHIP;                            // Chem. Mater. 25(15), 2911–2920 (2013) Herfindahl–Hirschman Index (HHI), HHIP: for elemental production, Uncertinities in HHI_P: C,O,F,Cl,Sc,Ga,Rb,Ru,Rh,Cs,Hf,Os,Ir,Tl.      // DU 2019/05/17
double HHIR;                            // Chem. Mater. 25(15), 2911–2920 (2013) Herfindahl–Hirschman Index (HHI), HHIR: for elemental reserves,   Uncertinities in HHI_R: Be,C,N,O,F,Na,Mg,Al,Si,S,Cl,Ca,Sc,Ga,Ge,As,Rb,Sr,Ru,Rh,Pd,In,Cs,Hf,Os,Ir,Pt,Tl.      // DU 2019/05/17
double xray_scatt;                  // shift+1 // All data collected from the NIST online tables: http://physics.nist.gov/PhysRefData/FFast/html/form.html//

// Xray_scatt_vector All data collected from the NIST online tables
// http://physics.nist.gov/PhysRefData/FFast/html/form.html
// All data are ideally for f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
// These are for E=7.9026keV (Cu-alpha is wavelength=1.5418A, E=8.0416keV).

// All data collected from the online tables:
// http://www-cxro.lbl.gov/optical_constants/pert_form.html
// All data are f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV].

// [AFLOW]STOP=DECLARATION
// operators/functions                                    // operator/functions
friend ostream& operator<<(ostream &,const xelement&);    // print
xelement Iinitialize(uint Z);                      // function to clean up the name

private:                                                    //
void free();                                              // free space
};
}
*/


std::vector<xelement::xelement> velement(NUM_ELEMENTS);        // store starting from ONE

namespace pflow {
  void XelementPrint(string options,ostream& oss) {
    bool LDEBUG=0;//(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XHOST.sPID << "pflow::XelementPrint [BEGIN]" << endl;
    if(LDEBUG) cerr << "options=" << options << endl;
    if(LDEBUG) cerr << "velement.size()=" << velement.size() << endl;

    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(LDEBUG) cerr << "tokens.size()=" << tokens.size() << endl;
    if(tokens.size()==0) {
      init::ErrorOption(cout,options,"pflow::XelementPrint","aflow --element=Z|name|symbol[,property[,property]....]");
      exit(0);
    } 
    // move on
    string species=tokens.at(0);
    uint Z=0; // some defaults
    // try with number
    if(tokens.size()>=1) if(aurostd::string2utype<uint>(species)>0) Z=aurostd::string2utype<uint>(species);
    if(Z>103) {
      init::ErrorOption(cout,options,"pflow::XelementPrint",aurostd::liststring2string("aflow --element=Z|name|symbol[,property[,property]....]","Z outside [1,103] or name or symbol unrecognized"));
    }
    // try with symbol
    if(Z==0) {
      for(uint i=1;i<=103;i++)
	if(aurostd::toupper(species)==aurostd::toupper(xelement::xelement(i).symbol)) Z=i;
    }
    // try with name
    if(Z==0) {
      for(uint i=1;i<=103;i++)
	if(aurostd::toupper(species)==aurostd::toupper(xelement::xelement(i).name)) Z=i;
    }
    
    if(LDEBUG) cerr << "Z=" << Z << endl;
    oss << "AFLOW element property finder" << endl;
    if(Z>0) {
      oss << "Element Z=" << xelement::xelement(Z).Z << " - " << xelement::xelement(Z).symbol << " - " << xelement::xelement(Z).name << endl;
      string space="        ";
    
      // found
    
      //      cerr << xelement::xelement(3).name << endl;
      //      cerr << xelement::xelement("Li").name << endl;
      //      cerr << xelement::xelement("LiThIuM").name << endl;
      
      // now look at properties
      if(tokens.size()>=2) {
	for(uint i=1;i<tokens.size();i++) {
	  string c=aurostd::toupper(tokens.at(i));
	  vector<string> vs; uint len=52;
	  vs.clear();
	  int prec=10;
	  if(c=="ALL" || c==aurostd::toupper("name")) vs.push_back(aurostd::PaddedPOST("name="+xelement::xelement(Z).name,len));
	  if(c=="ALL" || c==aurostd::toupper("symbol")) vs.push_back(aurostd::PaddedPOST("symbol="+xelement::xelement(Z).symbol,len));
	  if(c=="ALL" || c==aurostd::toupper("Z")) vs.push_back(aurostd::PaddedPOST("Z="+aurostd::utype2string(xelement::xelement(Z).Z),len));
	  if(c=="ALL" || c==aurostd::toupper("Period")) vs.push_back(aurostd::PaddedPOST("Period="+aurostd::utype2string(xelement::xelement(Z).Period),len));
	  if(c=="ALL" || c==aurostd::toupper("Group")) vs.push_back(aurostd::PaddedPOST("Group="+aurostd::utype2string(xelement::xelement(Z).Group),len));
	  if(c=="ALL" || c==aurostd::toupper("Series")) vs.push_back(aurostd::PaddedPOST("Series="+xelement::xelement(Z).Series,len));
	  if(c=="ALL" || c==aurostd::toupper("Block")) vs.push_back(aurostd::PaddedPOST("Block="+xelement::xelement(Z).Block,len));
	  if(c=="ALL" || c==aurostd::toupper("mass")) vs.push_back(aurostd::PaddedPOST("mass="+aurostd::utype2string(xelement::xelement(Z).mass,prec),len)+"// (kg)");
	  if(c=="ALL" || c==aurostd::toupper("MolarVolume")) vs.push_back(aurostd::PaddedPOST("MolarVolume="+aurostd::utype2string(xelement::xelement(Z).MolarVolume,prec),len)+"// (m^3/mol)");
	  if(c=="ALL" || c==aurostd::toupper("volume")) vs.push_back(aurostd::PaddedPOST("volume="+aurostd::utype2string(xelement::xelement(Z).volume,prec),len)+"// A^3");
	  if(c=="ALL" || c==aurostd::toupper("Miedema_Vm")) vs.push_back(aurostd::PaddedPOST("Miedema_Vm="+aurostd::utype2string(xelement::xelement(Z).Miedema_Vm,prec),len)+"// (V_m^{2/3} in (cm^2))");
	  if(c=="ALL" || c==aurostd::toupper("valence_std")) vs.push_back(aurostd::PaddedPOST("valence_std="+aurostd::utype2string(xelement::xelement(Z).valence_std),len));
	  if(c=="ALL" || c==aurostd::toupper("valence_iupac")) vs.push_back(aurostd::PaddedPOST("valence_iupac="+aurostd::utype2string(xelement::xelement(Z).valence_iupac),len));
	  if(c=="ALL" || c==aurostd::toupper("valence_PT")) vs.push_back(aurostd::PaddedPOST("valence_PT="+aurostd::utype2string(xelement::xelement(Z).valence_PT),len));
	  if(c=="ALL" || c==aurostd::toupper("Density_PT")) vs.push_back(aurostd::PaddedPOST("Density_PT="+aurostd::utype2string(xelement::xelement(Z).Density_PT,prec),len)+"// (g/cm^3)");
	  if(c=="ALL" || c==aurostd::toupper("crystal")) vs.push_back(aurostd::PaddedPOST("crystal="+xelement::xelement(Z).crystal,len));
	  if(c=="ALL" || c==aurostd::toupper("CrystalStr_PT")) vs.push_back(aurostd::PaddedPOST("CrystalStr_PT="+xelement::xelement(Z).CrystalStr_PT,len));
	  if(c=="ALL" || c==aurostd::toupper("space_group")) vs.push_back(aurostd::PaddedPOST("space_group="+xelement::xelement(Z).space_group,len));
	  if(c=="ALL" || c==aurostd::toupper("space_group_number")) vs.push_back(aurostd::PaddedPOST("space_group_number="+aurostd::utype2string(xelement::xelement(Z).space_group_number),len));

	  if(c=="ALL" || c==aurostd::toupper("Pearson_coefficient")) vs.push_back(aurostd::PaddedPOST("Pearson_coefficient="+aurostd::utype2string(xelement::xelement(Z).Pearson_coefficient,prec),len));
	  if(c=="ALL" || c==aurostd::toupper("lattice_constant")) vs.push_back(aurostd::PaddedPOST("lattice_constant="+aurostd::utype2string(xelement::xelement(Z).lattice_constant[1],prec)+","+aurostd::utype2string(xelement::xelement(Z).lattice_constant[2],prec)+","+aurostd::utype2string(xelement::xelement(Z).lattice_constant[3],prec),len)+"// (pm)");
	  if(c=="ALL" || c==aurostd::toupper("lattice_angle")) vs.push_back(aurostd::PaddedPOST("lattice_angle="+aurostd::utype2string(xelement::xelement(Z).lattice_angle[1],prec)+","+aurostd::utype2string(xelement::xelement(Z).lattice_angle[2],prec)+","+aurostd::utype2string(xelement::xelement(Z).lattice_angle[3],prec),len)+"// (rad)");
	  if(c=="ALL" || c==aurostd::toupper("phase")) vs.push_back(aurostd::PaddedPOST("phase="+xelement::xelement(Z).phase,len));
	  if(c=="ALL" || c==aurostd::toupper("radius")) vs.push_back(aurostd::PaddedPOST("radius="+aurostd::utype2string(xelement::xelement(Z).radius,prec),len)+"// (nm)");
	  if(c=="ALL" || c==aurostd::toupper("radius_PT")) vs.push_back(aurostd::PaddedPOST("radius_PT="+aurostd::utype2string(xelement::xelement(Z).radius_PT,prec),len)+"// (pm)");
	  if(c=="ALL" || c==aurostd::toupper("radius_covalent_PT")) vs.push_back(aurostd::PaddedPOST("radius_covalent_PT="+aurostd::utype2string(xelement::xelement(Z).radius_covalent_PT,prec),len)+"// (pm)");
	  if(c=="ALL" || c==aurostd::toupper("radius_covalent")) vs.push_back(aurostd::PaddedPOST("radius_covalent="+aurostd::utype2string(xelement::xelement(Z).radius_covalent,prec),len)+"// (Angstrom)");
	  if(c=="ALL" || c==aurostd::toupper("radius_VanDerWaals_PT")) vs.push_back(aurostd::PaddedPOST("radius_VanDerWaals_PT="+aurostd::utype2string(xelement::xelement(Z).radius_VanDerWaals_PT,prec),len)+"// (pm)");
	  if(c=="ALL" || c==aurostd::toupper("radii_Ghosh08")) vs.push_back(aurostd::PaddedPOST("radii_Ghosh08="+aurostd::utype2string(xelement::xelement(Z).radii_Ghosh08,prec),len)+"// (Angstrom)");
	  if(c=="ALL" || c==aurostd::toupper("radii_Slatter")) vs.push_back(aurostd::PaddedPOST("radii_Slatter="+aurostd::utype2string(xelement::xelement(Z).radii_Slatter,prec),len)+"// (Angstrom)");
	  if(c=="ALL" || c==aurostd::toupper("radii_Pyykko")) vs.push_back(aurostd::PaddedPOST("radii_Pyykko="+aurostd::utype2string(xelement::xelement(Z).radii_Pyykko,prec),len)+"// (pm)");

	  if(c=="ALL" || c==aurostd::toupper("electrical_conductivity")) vs.push_back(aurostd::PaddedPOST("electrical_conductivity="+aurostd::utype2string(xelement::xelement(Z).electrical_conductivity,prec),len)+"// (S/m)");
	  if(c=="ALL" || c==aurostd::toupper("electronegativity_vec")) vs.push_back(aurostd::PaddedPOST("electronegativity_vec="+aurostd::utype2string(xelement::xelement(Z).electronegativity_vec,prec),len));
	  if(c=="ALL" || c==aurostd::toupper("hardness_Ghosh")) vs.push_back(aurostd::PaddedPOST("hardness_Ghosh="+aurostd::utype2string(xelement::xelement(Z).hardness_Ghosh,prec),len)+"// (eV)");
	    if(c=="ALL" || c==aurostd::toupper("electronegativity_Pearson")) vs.push_back(aurostd::PaddedPOST("electronegativity_Pearson="+aurostd::utype2string(xelement::xelement(Z).electronegativity_Pearson,prec),len)+"// (eV)");
	    if(c=="ALL" || c==aurostd::toupper("electronegativity_Ghosh")) vs.push_back(aurostd::PaddedPOST("electronegativity_Ghosh="+aurostd::utype2string(xelement::xelement(Z).electronegativity_Ghosh,prec),len)+"// (eV)");
	    if(c=="ALL" || c==aurostd::toupper("electron_affinity_PT")) vs.push_back(aurostd::PaddedPOST("electron_affinity_PT="+aurostd::utype2string(xelement::xelement(Z).electron_affinity_PT,prec),len)+"// (kJ/mol)");
	    if(c=="ALL" || c==aurostd::toupper("Miedema_phi_star")) vs.push_back(aurostd::PaddedPOST("Miedema_phi_star="+aurostd::utype2string(xelement::xelement(Z).Miedema_phi_star,prec),len)+"// (V) (phi^star)");
	    if(c=="ALL" || c==aurostd::toupper("Miedema_nws")) vs.push_back(aurostd::PaddedPOST("Miedema_nws="+aurostd::utype2string(xelement::xelement(Z).Miedema_nws,prec),len)+"// (d.u.)^1/3 n_{ws}^{1/3}");
	    if(c=="ALL" || c==aurostd::toupper("Miedema_gamma_s")) vs.push_back(aurostd::PaddedPOST("Miedema_gamma_s="+aurostd::utype2string(xelement::xelement(Z).Miedema_gamma_s,prec),len)+"// (mJ/m^2)");
	    if(c=="ALL" || c==aurostd::toupper("Pettifor_scale")) vs.push_back(aurostd::PaddedPOST("Pettifor_scale="+aurostd::utype2string(xelement::xelement(Z).Pettifor_scale,prec),len)); 
	    if(c=="ALL" || c==aurostd::toupper("boiling_point")) vs.push_back(aurostd::PaddedPOST("boiling_point="+aurostd::utype2string(xelement::xelement(Z).boiling_point,prec),len)+"// (Celsius)");
	    if(c=="ALL" || c==aurostd::toupper("melting_point")) vs.push_back(aurostd::PaddedPOST("melting_point="+aurostd::utype2string(xelement::xelement(Z).melting_point,prec),len)+"// (Celsius)");
	    if(c=="ALL" || c==aurostd::toupper("vaporization_heat_PT")) vs.push_back(aurostd::PaddedPOST("vaporization_heat_PT="+aurostd::utype2string(xelement::xelement(Z).vaporization_heat_PT,prec),len)+"// (kJ/mol)");
	    if(c=="ALL" || c==aurostd::toupper("specific_heat_PT")) vs.push_back(aurostd::PaddedPOST("specific_heat_PT="+aurostd::utype2string(xelement::xelement(Z).specific_heat_PT,prec),len)+"// (J/(kg.K))");
	    if(c=="ALL" || c==aurostd::toupper("critical_Pressure")) vs.push_back(aurostd::PaddedPOST("critical_Pressure="+aurostd::utype2string(xelement::xelement(Z).critical_Pressure,prec),len)+"// (Atm) "); 
	    if(c=="ALL" || c==aurostd::toupper("critical_Temperature_PT")) vs.push_back(aurostd::PaddedPOST("critical_Temperature_PT="+aurostd::utype2string(xelement::xelement(Z).critical_Temperature_PT,prec),len)+"// (K)"); 
	    if(c=="ALL" || c==aurostd::toupper("thermal_expansion")) vs.push_back(aurostd::PaddedPOST("thermal_expansion="+aurostd::utype2string(xelement::xelement(Z).thermal_expansion,prec),len)+"// (K^{-1})");
	    if(c=="ALL" || c==aurostd::toupper("thermal_conductivity")) vs.push_back(aurostd::PaddedPOST("thermal_conductivity="+aurostd::utype2string(xelement::xelement(Z).thermal_conductivity,prec),len)+"// (W/(mK))");
	    if(c=="ALL" || c==aurostd::toupper("Brinelll_hardness")) vs.push_back(aurostd::PaddedPOST("Brinelll_hardness="+aurostd::utype2string(xelement::xelement(Z).Brinelll_hardness,prec),len)+"// (MPa)");
	    if(c=="ALL" || c==aurostd::toupper("Mohs_hardness")) vs.push_back(aurostd::PaddedPOST("Mohs_hardness="+aurostd::utype2string(xelement::xelement(Z).Mohs_hardness,prec),len));
	    if(c=="ALL" || c==aurostd::toupper("Vickers_hardness")) vs.push_back(aurostd::PaddedPOST("Vickers_hardness="+aurostd::utype2string(xelement::xelement(Z).Vickers_hardness,prec),len)+"// (MPa)");
	    if(c=="ALL" || c==aurostd::toupper("Hardness_Pearson")) vs.push_back(aurostd::PaddedPOST("Hardness_Pearson="+aurostd::utype2string(xelement::xelement(Z).Hardness_Pearson,prec),len)+"// (eV)");
	    if(c=="ALL" || c==aurostd::toupper("Hardness_Putz")) vs.push_back(aurostd::PaddedPOST("Hardness_Putz="+aurostd::utype2string(xelement::xelement(Z).Hardness_Putz,prec),len)+"// (eV/atom)");
	    if(c=="ALL" || c==aurostd::toupper("Hardness_RB")) vs.push_back(aurostd::PaddedPOST("Hardness_RB="+aurostd::utype2string(xelement::xelement(Z).Hardness_RB,prec),len)+"// (eV)");
	    if(c=="ALL" || c==aurostd::toupper("shear_modulus")) vs.push_back(aurostd::PaddedPOST("shear_modulus="+aurostd::utype2string(xelement::xelement(Z).shear_modulus,prec),len)+"// (GPa)");
	    if(c=="ALL" || c==aurostd::toupper("Young_modulus")) vs.push_back(aurostd::PaddedPOST("Young_modulus="+aurostd::utype2string(xelement::xelement(Z).Young_modulus,prec),len)+"// (GPa)");
	    if(c=="ALL" || c==aurostd::toupper("bulk_modulus")) vs.push_back(aurostd::PaddedPOST("bulk_modulus="+aurostd::utype2string(xelement::xelement(Z).bulk_modulus,prec),len)+"// (GPa)");
	    if(c=="ALL" || c==aurostd::toupper("Poisson_ratio_PT")) vs.push_back(aurostd::PaddedPOST("Poisson_ratio_PT="+aurostd::utype2string(xelement::xelement(Z).Poisson_ratio_PT,prec),len));
	    if(c=="ALL" || c==aurostd::toupper("Miedema_BVm")) vs.push_back(aurostd::PaddedPOST("Miedema_BVm="+aurostd::utype2string(xelement::xelement(Z).Miedema_BVm,prec),len)+"// (kJ/mole)");

	    if(c=="ALL" || c==aurostd::toupper("Magnetic_Type_PT")) vs.push_back(aurostd::PaddedPOST("Magnetic_Type_PT="+xelement::xelement(Z).Magnetic_Type_PT,len));
	    if(c=="ALL" || c==aurostd::toupper("Mass_Magnetic_Susceptibility")) vs.push_back(aurostd::PaddedPOST("Mass_Magnetic_Susceptibility="+aurostd::utype2string(xelement::xelement(Z).Mass_Magnetic_Susceptibility,prec),len)+"// (m^3/K)");
	    if(c=="ALL" || c==aurostd::toupper("Volume_Magnetic_Susceptibility")) vs.push_back(aurostd::PaddedPOST("Volume_Magnetic_Susceptibility="+aurostd::utype2string(xelement::xelement(Z).Volume_Magnetic_Susceptibility,prec),len));
	    if(c=="ALL" || c==aurostd::toupper("Molar_Magnetic_Susceptibility")) vs.push_back(aurostd::PaddedPOST("Molar_Magnetic_Susceptibility="+aurostd::utype2string(xelement::xelement(Z).Molar_Magnetic_Susceptibility,prec),len)+"// (m^3/mol)");
	    if(c=="ALL" || c==aurostd::toupper("Curie_point")) vs.push_back(aurostd::PaddedPOST("Curie_point="+aurostd::utype2string(xelement::xelement(Z).Curie_point,prec),len)+"// (K)");
	
	    if(c=="ALL" || c==aurostd::toupper("refractive_index")) vs.push_back(aurostd::PaddedPOST("refractive_index="+aurostd::utype2string(xelement::xelement(Z).refractive_index,prec),len));
	    if(c=="ALL" || c==aurostd::toupper("color_PT")) vs.push_back(aurostd::PaddedPOST("color_PT="+xelement::xelement(Z).color_PT,len));
	    if(c=="ALL" || c==aurostd::toupper("HHIP")) vs.push_back(aurostd::PaddedPOST("HHIP="+aurostd::utype2string(xelement::xelement(Z).HHIP,prec),len));
	    if(c=="ALL" || c==aurostd::toupper("HHIR")) vs.push_back(aurostd::PaddedPOST("HHIR="+aurostd::utype2string(xelement::xelement(Z).HHIR,prec),len));
	    if(c=="ALL" || c==aurostd::toupper("xray_scatt")) vs.push_back(aurostd::PaddedPOST("xray_scatt="+aurostd::utype2string(xelement::xelement(Z).xray_scatt,prec),len)+"// shift+1");

	    if(vs.size())
	      for(uint j=0;j<vs.size();j++)
		oss << vs.at(j) << endl;
	}
      }
    }
    
    if(LDEBUG) cerr << XHOST.sPID << "pflow::XelementPrint [END]" << endl;
  }
}

/*
  std::vector<string> vatom_symbol(NUM_ELEMENTS);   // store starting from ONE // DONE
  std::vector<string> vatom_name(NUM_ELEMENTS);   // store starting from ONE // DONE
  std::vector<double> vatom_mass(NUM_ELEMENTS);     // store starting from ONE // DONE
  std::vector<double> vatom_volume(NUM_ELEMENTS);       // store starting from ONE // DONE
  std::vector<int> vatom_valence_iupac(NUM_ELEMENTS);   // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry) // DONE
  std::vector<int> vatom_valence_std(NUM_ELEMENTS);     // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry) // DONE
  std::vector<double> vatom_miedema_phi_star(NUM_ELEMENTS); // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28  
  std::vector<double> vatom_miedema_nws(NUM_ELEMENTS);      // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
  std::vector<double> vatom_miedema_Vm(NUM_ELEMENTS);       // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
  std::vector<double> vatom_miedema_gamma_s(NUM_ELEMENTS);  // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
  std::vector<double> vatom_miedema_BVm(NUM_ELEMENTS);      // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
  // for lanthines from J.A. Alonso and N.H. March. Electrons in Metals and Alloys, Academic Press, London (1989) (except La)
  std::vector<double> vatom_radius(NUM_ELEMENTS);       // store starting from ONE  // DONE
  std::vector<double> vatom_radius_covalent(NUM_ELEMENTS);// store starting from ONE// DX and CO -20170904 
  std::vector<double> vatom_electronegativity(NUM_ELEMENTS);       // store starting from ONE
  std::vector<string> vatom_crystal(NUM_ELEMENTS);       // store starting from ONE  // DONE
  std::vector<double> vatom_xray_scatt(NUM_ELEMENTS);        // store starting from ONE
  std::vector<double> vatom_pettifor_scale(NUM_ELEMENTS);        // store starting from ONE Chemical Scale Pettifor Solid State Communications 51 31-34 1984
  std::vector<double> vatom_pearson_coefficient(NUM_ELEMENTS);   // ME20181020 Pearson mass deviation coefficient
  
*/
  
namespace xelement {

    
  // initialize them all
  void Initialize(void) { for(uint Z=0;Z<NUM_ELEMENTS;Z++) { velement.at(Z)=xelement(Z); } }
  string symbol2name(const string& symbol) { return xelement(symbol).name; }
  string name2symbol(const string& name) { return xelement(name).symbol; }
  int symbol2Z(const string& symbol) { return xelement(symbol).Z; }
  string Z2symbol(const int& Z) { return xelement(Z).symbol; }
  string Z2name(const int& Z) { return xelement(Z).name; }
  int name2Z(const string& name) { return xelement(name).Z; }

  // constructors
  xelement::xelement() {
    // DEFAULT
    verbose=FALSE;
    // [AFLOW]START=CONSTRUCTOR
    Z=0;
    symbol="XX";//"UNDEFINED";
    name="UNDEFINED";
    Period=NNN;
    Group=NNN; 
    Series="UNDEFINED";
    Block="nnn";      
    //                                          
    mass=NNN;//  AMU2KILOGRAM goes inside.
    MolarVolume=NNN;  
    volume=NNN;      
    Miedema_Vm=NNN;      
    //
    valence_std=NNN;  
    valence_iupac=NNN;
    valence_PT=NNN;       
    Density_PT=NNN;       
    crystal="nnn";    
    CrystalStr_PT="UNDEFINED";
    space_group="nnn";     
    space_group_number=NNN;
    Pearson_coefficient=NNN;
    lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
    lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN; 
    phase="nnn";         
    radius=NNN;         
    radius_PT=NNN;          
    radius_covalent_PT=NNN;   
    radius_covalent=NNN;  
    radius_VanDerWaals_PT=NNN;
    radii_Ghosh08=NNN;         
    radii_Slatter=NNN;         
    radii_Pyykko=NNN;          
    //                                          
    electrical_conductivity=NNN;
    electronegativity_vec=NNN;    
    hardness_Ghosh=NNN;            
    electronegativity_Pearson=NNN;           
    electronegativity_Ghosh=NNN;             
    electron_affinity_PT=NNN;      
    Miedema_phi_star=NNN;         
    Miedema_nws=NNN;              
    Miedema_gamma_s=NNN;          
    //
    Pettifor_scale=NNN;          
    //
    boiling_point=NNN;         
    melting_point=NNN;         
    vaporization_heat_PT=NNN;     
    specific_heat_PT=NNN;         
    critical_Pressure=NNN;     
    critical_Temperature_PT=NNN;  
    thermal_expansion=NNN;     
    thermal_conductivity=NNN;  
    //                                         
    Brinelll_hardness=NNN;
    Mohs_hardness=NNN;    
    Vickers_hardness=NNN; 
    Hardness_Pearson=NNN;   
    Hardness_Putz=NNN;      
    Hardness_RB=NNN;        
    shear_modulus=NNN;    
    Young_modulus=NNN;    
    bulk_modulus=NNN;     
    Poisson_ratio_PT=NNN;    
    Miedema_BVm=NNN;        
    //
    Magnetic_Type_PT="UNDEFINED";     
    Mass_Magnetic_Susceptibility=NNN;
    Volume_Magnetic_Susceptibility=NNN;
    Molar_Magnetic_Susceptibility=NNN; 
    Curie_point=NNN;                  
    //
    refractive_index=NNN;             
    color_PT="UNDEFINED";               
    //
    HHIP=NNN;                           
    HHIR=NNN;                           
    xray_scatt=NNN;   
    // [AFLOW]STOP=CONSTRUCTOR
  }

  // destructor
  xelement::~xelement() {
    free();
  }

  void xelement::free() {
    // will populate
    // [AFLOW]START=FREE
    // [AFLOW]STOP=FREE
  }

  const xelement& xelement::operator=(const xelement& b) {      // operator=
    if(this != &b) {
      free();
      // will populate
      verbose=b.verbose;
      // [AFLOW]START=ASSIGNMENT
      Z=b.Z;
      symbol=b.symbol;
      name=b.name;
      Period=b.Period;
      Group=b.Group; 
      Series=b.Series;
      Block=b.Block;      
      //                                          
      mass=b.mass;
      MolarVolume=b.MolarVolume;  
      volume=b.volume;      
      Miedema_Vm=b.Miedema_Vm;      
      //
      valence_std=b.valence_std;  
      valence_iupac=b.valence_iupac;
      valence_PT=b.valence_PT;       
      Density_PT=b.Density_PT;       
      crystal=b.crystal;    
      CrystalStr_PT=b.CrystalStr_PT;
      space_group=b.space_group;
      space_group_number=b.space_group_number;    
      Pearson_coefficient=b.Pearson_coefficient;
      lattice_constant=b.lattice_constant; 
      lattice_angle=b.lattice_angle;   
      phase=b.phase;
      radius=b.radius;         
      radius_PT=b.radius_PT;          
      radius_covalent_PT=b.radius_covalent_PT;   
      radius_covalent=b.radius_covalent;  
      radius_VanDerWaals_PT=b.radius_VanDerWaals_PT;
      radii_Ghosh08=b.radii_Ghosh08;         
      radii_Slatter=b.radii_Slatter;         
      radii_Pyykko=b.radii_Pyykko;          
      //                                          
      electrical_conductivity=b.electrical_conductivity;
      electronegativity_vec=b.electronegativity_vec;    
      hardness_Ghosh=b.hardness_Ghosh;            
      electronegativity_Pearson=b.electronegativity_Pearson;           
      electronegativity_Ghosh=b.electronegativity_Ghosh;             
      electron_affinity_PT=b.electron_affinity_PT;      
      Miedema_phi_star=b.Miedema_phi_star;         
      Miedema_nws=b.Miedema_nws;              
      Miedema_gamma_s=b.Miedema_gamma_s;          
      //
      Pettifor_scale=b.Pettifor_scale;          
      //
      boiling_point=b.boiling_point;         
      melting_point=b.melting_point;         
      vaporization_heat_PT=b.vaporization_heat_PT;     
      specific_heat_PT=b.specific_heat_PT;         
      critical_Pressure=b.critical_Pressure;     
      critical_Temperature_PT=b.critical_Temperature_PT;  
      thermal_expansion=b.thermal_expansion;     
      thermal_conductivity=b.thermal_conductivity;  
      //                                         
      Brinelll_hardness=b.Brinelll_hardness;
      Mohs_hardness=b.Mohs_hardness;    
      Vickers_hardness=b.Vickers_hardness; 
      Hardness_Pearson=b.Hardness_Pearson;   
      Hardness_Putz=b.Hardness_Putz;      
      Hardness_RB=b.Hardness_RB;        
      shear_modulus=b.shear_modulus;    
      Young_modulus=b.Young_modulus;    
      bulk_modulus=b.bulk_modulus;     
      Poisson_ratio_PT=b.Poisson_ratio_PT;    
      Miedema_BVm=b.Miedema_BVm;        
      //
      Magnetic_Type_PT=b.Magnetic_Type_PT;
      Mass_Magnetic_Susceptibility=b.Mass_Magnetic_Susceptibility;
      Volume_Magnetic_Susceptibility=b.Volume_Magnetic_Susceptibility;
      Molar_Magnetic_Susceptibility=b.Molar_Magnetic_Susceptibility; 
      Curie_point=b.Curie_point;                  
      //
      refractive_index=b.refractive_index;             
      color_PT=b.color_PT;         
      //
      HHIP=b.HHIP;                           
      HHIR=b.HHIR;                           
      xray_scatt=b.xray_scatt;    
      // [AFLOW]STOP=ASSIGNMENT
    }
    return *this;
  }

  void xelement::clear(){
    xelement a;(*this)=a;
  }

  ostream& operator<<(ostream& oss,const xelement& element) {
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    // [AFLOW]START=COUT
    oss << "verbose=" << element.verbose << endl;
    // [AFLOW]STOP=COUT
    return oss;
  }
  
  
  // ********************************************************************************************************************************************************
  // constructor by name or symbol
  xelement::xelement(string element) {
    free();
    // DEFAULT
    verbose=FALSE;
    uint Z=0;
  
    // try with symbol
    if(Z==0) {
      for(uint i=1;i<=103;i++)
	if(aurostd::toupper(element)==aurostd::toupper(xelement(i).symbol)) Z=i;
    }
    // try with name
    if(Z==0) {
      for(uint i=1;i<=103;i++)
	if(aurostd::toupper(element)==aurostd::toupper(xelement(i).name)) Z=i;
    }
    if(Z!=0) (*this)=xelement(Z);

  }

  // ********************************************************************************************************************************************************
  // constructor by Z
  xelement::xelement(uint ZZ) {
    free();
    // DEFAULT
    verbose=FALSE;

    // OFFSET
    if(ZZ==0) {
      xelement a;
      (*this)=a;
      mass=0.0;// override
      valence_iupac=0;// override
      valence_std=0;// override
    }
    // ROW 1
    // s-electron systems

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Hydrogen
    // Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen
    if(ZZ==1) { // Hydrogen
      Z=ZZ;
      symbol="H";
      name="Hydrogen";
      Period=1;
      Group=1;
      Series="Nonmetal";
      Block="s";
      mass=AMU2KILOGRAM*1.0079;
      MolarVolume=0.01121;
      volume=0.75110;
      Miedema_Vm=NNN;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      Density_PT=0.899E-4;
      crystal="hex";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.00011460743;
      lattice_constant[1]=470;lattice_constant[2]=470;lattice_constant[3]=340;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Gas";
      radius=0.046;
      radius_PT=53;
      radius_covalent=0.31;
      radius_covalent_PT=31;
      radius_VanDerWaals_PT=120;
      radii_Ghosh08=0.5292;
      radii_Slatter=0.25;
      radii_Pyykko=0.32;
      electrical_conductivity=NNN;
      electronegativity_vec=2.10;
      hardness_Ghosh=6.4299;
      electronegativity_Pearson=7.18;
      electronegativity_Ghosh=7.178;
      electron_affinity_PT=72.8;
      Miedema_phi_star=5.2;
      Miedema_nws=1.5;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=-252.87;
      melting_point=-259.14;
      vaporization_heat_PT=0.452;
      specific_heat_PT=14300;
      critical_Pressure=12.76;
      critical_Temperature_PT=32.97;
      thermal_expansion=NNN;
      thermal_conductivity=0.1805;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=6.43;
      Hardness_Putz=6.45;
      Hardness_RB=6.83;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-2.48E-8;
      Volume_Magnetic_Susceptibility=-2.23E-9;
      Molar_Magnetic_Susceptibility=-4.999E-11;
      Curie_point=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000132;
      HHIP=NNN;
      HHIR=NNN;
      xray_scatt=1.000;
      // H volume wrong *dimer* MIEDEMA =PAUL VAN DER PUT book
    }
    // [AFLOW]STOP=Hydrogen
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Helium
    // Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium
    if(ZZ==2) { // Helium
      Z=ZZ;
      symbol="He";
      name="Helium";
      Period=1;
      Group=18;
      Series="NobleGas";
      Block="s";
      mass=AMU2KILOGRAM*4.0026;
      MolarVolume=0.022424;
      volume=-1.000;
      Miedema_Vm=NNN;
      valence_std=0;
      valence_iupac=0;
      valence_PT=0;
      Density_PT=1.785E-4;
      crystal="hcp";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=8.32328E-8;
      lattice_constant[1]=424.2;lattice_constant[2]=424.2;lattice_constant[3]=424.2;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Gas";
      radius=NNN;
      radius_PT=31;
      radius_covalent=0.28;
      radius_covalent_PT=28;
      radius_VanDerWaals_PT=140;
      radii_Ghosh08=0.3113;
      radii_Slatter=NNN;
      radii_Pyykko=0.46;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=12.5449;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=12.046;
      electron_affinity_PT=0;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=-268.93;
      melting_point=NNN;
      vaporization_heat_PT=0.083;
      specific_heat_PT=5193.1;
      critical_Pressure=2.24;
      critical_Temperature_PT=5.19;
      thermal_expansion=NNN;
      thermal_conductivity=0.1513;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=25.79;
      Hardness_RB=16.88;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-5.9E-9;
      Volume_Magnetic_Susceptibility=-1.05E-9;
      Molar_Magnetic_Susceptibility=-2.36E-11;
      Curie_point=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000035;
      HHIP=3200;
      HHIR=3900;
      xray_scatt=2.000;
      // He
    }
    // [AFLOW]STOP=Helium
    // ********************************************************************************************************************************************************

    // ROW2
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Lithium
    // Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium
    if(ZZ==3) { // Lithium
      Z=ZZ;
      symbol="Li";
      name="Lithium";
      Period=2;
      Group=1;
      Series="AlkaliMetal";
      Block="s";
      mass=AMU2KILOGRAM*6.941;
      MolarVolume=0.00001297;
      volume=20.24110;
      Miedema_Vm=5.5;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      Density_PT=0.535;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=0.0014588232;
      lattice_constant[1]=351;lattice_constant[2]=351;lattice_constant[3]=351;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.152;
      radius_PT=167;
      radius_covalent=1.28;
      radius_covalent_PT=128;
      radius_VanDerWaals_PT=182;
      radii_Ghosh08=1.6283;
      radii_Slatter=1.45;
      radii_Pyykko=1.33;
      electrical_conductivity=1.1E7;
      electronegativity_vec=0.98;
      hardness_Ghosh=2.3746;
      electronegativity_Pearson=3.01;
      electronegativity_Ghosh=2.860;
      electron_affinity_PT=59.6;
      Miedema_phi_star=2.85;
      Miedema_nws=0.98;
      Miedema_gamma_s=530;
      Pettifor_scale=0.45;
      boiling_point=1342;
      melting_point=180.54;
      vaporization_heat_PT=147;
      specific_heat_PT=3570;
      critical_Pressure=661.2;
      critical_Temperature_PT=3223;
      thermal_expansion=0.000046;
      thermal_conductivity=85;
      Brinelll_hardness=NNN;
      Mohs_hardness=0.6;
      Vickers_hardness=NNN;
      Hardness_Pearson=2.39;
      Hardness_Putz=0.65;
      Hardness_RB=3.06;
      shear_modulus=4.2;
      Young_modulus=4.9;
      bulk_modulus=11;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=1.5;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=2.56E-8;
      Volume_Magnetic_Susceptibility=0.0000137;
      Molar_Magnetic_Susceptibility=1.78E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=2900;
      HHIR=4200;
      xray_scatt=3.00145;
      // Li
    }
    // [AFLOW]STOP=Lithium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Beryllium
    // Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium
    if(ZZ==4) { // Beryllium
      Z=ZZ;
      symbol="Be";
      name="Beryllium";
      Period=2;
      Group=2;
      Series="AlkalineEarthMetal";
      Block="s";
      mass=AMU2KILOGRAM*9.0122;
      MolarVolume=4.8767E-6;
      volume=7.83290;
      Miedema_Vm=2.9;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      Density_PT=1.848;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=228.58;lattice_constant[2]=228.58;lattice_constant[3]=358.43;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.114;
      radius_PT=112;
      radius_covalent=0.96;
      radius_covalent_PT=96;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.0855;
      radii_Slatter=1.05;
      radii_Pyykko=1.02;
      electrical_conductivity=2.5E7;
      electronegativity_vec=1.57;
      hardness_Ghosh=3.4968;
      electronegativity_Pearson=4.90;
      electronegativity_Ghosh=3.945;
      electron_affinity_PT=0;
      Miedema_phi_star=4.20;
      Miedema_nws=1.60;
      Miedema_gamma_s=1900;
      Pettifor_scale=1.50;
      boiling_point=2470;
      melting_point=1287;
      vaporization_heat_PT=297;
      specific_heat_PT=1820;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000113;
      thermal_conductivity=190;
      Brinelll_hardness=600;
      Mohs_hardness=5.5;
      Vickers_hardness=1670;
      Hardness_Pearson=4.50;
      Hardness_Putz=1.69;
      Hardness_RB=5.16;
      shear_modulus=132;
      Young_modulus=287;
      bulk_modulus=130;
      Poisson_ratio_PT=0.032;
      Miedema_BVm=4.9;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-1.26E-8;
      Volume_Magnetic_Susceptibility=-0.00002328;
      Molar_Magnetic_Susceptibility=-1.136E-10;
      Curie_point=NNN;
      color_PT="SLATEGRAY";
      refractive_index=NNN;
      HHIP=8000;
      HHIR=4000;
      /*xray_scatt=NNN;*/
      // Be
    }
    // [AFLOW]STOP=Beryllium
    // ********************************************************************************************************************************************************

    // p-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Boron
    // Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron
    if(ZZ==5) { // Boron
      Z=ZZ;
      symbol="B";
      name="Boron";
      Period=2;
      Group=13;
      Series="Metalloid";
      Block="p";
      mass=AMU2KILOGRAM*10.81;
      MolarVolume=4.3943E-6;
      volume=5.88420;
      Miedema_Vm=2.8;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=2.46;
      crystal="tet";
      CrystalStr_PT="Simple_Trigonal";
      space_group="R_3m";
      space_group_number=166;
      Pearson_coefficient=0.00135391428;
      lattice_constant[1]=506;lattice_constant[2]=506;lattice_constant[3]=506;
      lattice_angle[1]=1.01334;lattice_angle[2]=1.01334;lattice_angle[3]=1.01334;
      phase="Solid";
      radius=0.097;
      radius_PT=87;
      radius_covalent=0.84;
      radius_covalent_PT=85;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8141;
      radii_Slatter=0.85;
      radii_Pyykko=0.85;
      electrical_conductivity=0.0001;
      electronegativity_vec=2.04;
      hardness_Ghosh=4.6190;
      electronegativity_Pearson=4.29;
      electronegativity_Ghosh=5.031;
      electron_affinity_PT=26.7;
      Miedema_phi_star=4.75;
      Miedema_nws=1.55;
      Miedema_gamma_s=NNN;
      Pettifor_scale=2.00;
      boiling_point=4000;
      melting_point=2075;
      vaporization_heat_PT=507;
      specific_heat_PT=1030;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6E-6;
      thermal_conductivity=27;
      Brinelll_hardness=NNN;
      Mohs_hardness=9.3;
      Vickers_hardness=49000;
      Hardness_Pearson=4.01;
      Hardness_Putz=3.46;
      Hardness_RB=4.39;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=320;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-8.7E-9;
      Volume_Magnetic_Susceptibility=-0.0000214;
      Molar_Magnetic_Susceptibility=-9.41E-11;
      Curie_point=NNN;
      color_PT="BLACK";
      refractive_index=NNN;
      HHIP=2900;
      HHIR=2000;
      /*xray_scatt=NNN;*/
      // B
    }
    // [AFLOW]STOP=Boron
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Carbon
    // Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon
    if(ZZ==6) { // Carbon
      Z=ZZ;
      symbol="C";
      name="Carbon";
      Period=2;
      Group=14;
      Series="Nonmetal";
      Block="p";
      mass=AMU2KILOGRAM*12.011;
      MolarVolume=5.3146E-6;
      volume=5.59490;
      Miedema_Vm=1.8;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=2.26;
      crystal="dia";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.00007387218;
      lattice_constant[1]=246.4;lattice_constant[2]=246.4;lattice_constant[3]=671.1;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.077;
      radius_PT=67;
      radius_covalent=0.76;
      radius_covalent_PT=76;
      radius_VanDerWaals_PT=170;
      radii_Ghosh08=0.6513;
      radii_Slatter=0.70;
      radii_Pyykko=0.75;
      electrical_conductivity=100000;
      electronegativity_vec=2.55;
      hardness_Ghosh=5.7410;
      electronegativity_Pearson=6.27;
      electronegativity_Ghosh=6.116;
      electron_affinity_PT=153.9;
      Miedema_phi_star=6.20;
      Miedema_nws=1.90;
      Miedema_gamma_s=NNN;
      Pettifor_scale=2.50;
      boiling_point=4027;
      melting_point=3550;
      vaporization_heat_PT=715;
      specific_heat_PT=710;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=7.1E-6;
      thermal_conductivity=140;
      Brinelll_hardness=NNN;
      Mohs_hardness=0.5;
      Vickers_hardness=NNN;
      Hardness_Pearson=5.00;
      Hardness_Putz=6.21;
      Hardness_RB=5.49;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=33;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-6.2E-9;
      Volume_Magnetic_Susceptibility=-0.000014;
      Molar_Magnetic_Susceptibility=-7.45E-11;
      Curie_point=NNN;
      color_PT="BLACK";
      refractive_index=2.417;
      HHIP=500;
      HHIR=500;
      xray_scatt=6.019;
      // C//DX and CO  20170904 radius_covalent uses sp3 hybridization (most common)
    }
    // [AFLOW]STOP=Carbon
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Nitrogen
    // Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen
    if(ZZ==7) { // Nitrogen
      Z=ZZ;
      symbol="N";
      name="Nitrogen";
      Period=2;
      Group=15;
      Series="Nonmetal";
      Block="p";
      mass=AMU2KILOGRAM*14.0067;
      MolarVolume=0.011197;
      volume=7.59940;
      Miedema_Vm=2.2;
      valence_std=5;
      valence_iupac=5;
      valence_PT=3;
      Density_PT=12.51E-4;
      crystal="hex";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.00001857771;
      lattice_constant[1]=386.1;lattice_constant[2]=386.1;lattice_constant[3]=626.5;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Gas";
      radius=0.071;
      radius_PT=56;
      radius_covalent=0.71;
      radius_covalent_PT=71;
      radius_VanDerWaals_PT=155;
      radii_Ghosh08=0.5428;
      radii_Slatter=0.65;
      radii_Pyykko=0.71;
      electrical_conductivity=NNN;
      electronegativity_vec=3.04;
      hardness_Ghosh=6.8625;
      electronegativity_Pearson=7.30;
      electronegativity_Ghosh=7.209;
      electron_affinity_PT=7;
      Miedema_phi_star=7.00;
      Miedema_nws=1.60;
      Miedema_gamma_s=NNN;
      Pettifor_scale=3.00;
      boiling_point=-195.79;
      melting_point=-210.1;
      vaporization_heat_PT=2.79;
      specific_heat_PT=1040;
      critical_Pressure=33.46;
      critical_Temperature_PT=126.21;
      thermal_expansion=NNN;
      thermal_conductivity=0.02583;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=7.23;
      Hardness_Putz=9.59;
      Hardness_RB=8.59;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-5.4E-9;
      Volume_Magnetic_Susceptibility=-6.8E-9;
      Molar_Magnetic_Susceptibility=-1.5E-10;
      Curie_point=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000298;
      HHIP=1300;
      HHIR=500;
      /*xray_scatt=NNN;*/
      //N JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Nitrogen
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Oxygen
    // Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen
    if(ZZ==8) { // Oxygen
      Z=ZZ;
      symbol="O";
      name="Oxygen";
      Period=2;
      Group=16;
      Series="Chalcogen";
      Block="p";
      mass=AMU2KILOGRAM*15.9994;
      MolarVolume=0.011196;
      volume=7.78230;
      Miedema_Vm=2.656;
      valence_std=6;
      valence_iupac=2;
      valence_PT=2;
      Density_PT=14.29E-4;
      crystal="cub";
      CrystalStr_PT="Base-centered_Monoclinic";
      space_group="C12/m1";
      space_group_number=12;
      Pearson_coefficient=0.00003358805;
      lattice_constant[1]=540.3;lattice_constant[2]=342.9;lattice_constant[3]=508.6;
      lattice_angle[1]=PI/2;lattice_angle[2]=2.313085;lattice_angle[3]=PI/2;
      phase="Gas";
      radius=0.060;
      radius_PT=48;
      radius_covalent=0.66;
      radius_covalent_PT=66;
      radius_VanDerWaals_PT=152;
      radii_Ghosh08=0.4652;
      radii_Slatter=0.60;
      radii_Pyykko=0.63;
      electrical_conductivity=NNN;
      electronegativity_vec=3.44;
      hardness_Ghosh=7.9854;
      electronegativity_Pearson=7.54;
      electronegativity_Ghosh=8.287;
      electron_affinity_PT=141;
      Miedema_phi_star=6.97;
      Miedema_nws=1.70;
      Miedema_gamma_s=NNN;
      Pettifor_scale=3.50;
      boiling_point=-182.9;
      melting_point=-218.3;
      vaporization_heat_PT=3.41;
      specific_heat_PT=919;
      critical_Pressure=49.77;
      critical_Temperature_PT=154.59;
      thermal_expansion=NNN;
      thermal_conductivity=0.02658;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=6.08;
      Hardness_Putz=13.27;
      Hardness_RB=6.42;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.335E-6;
      Volume_Magnetic_Susceptibility=1.90772E-6;
      Molar_Magnetic_Susceptibility=4.27184E-8;
      Curie_point=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000271;
      HHIP=500;
      HHIR=500;
      xray_scatt=8.052;
      // O Table 27 of JUNKAI
    }
    // [AFLOW]STOP=Oxygen
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Fluorine
    // Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine
    if(ZZ==9) { // Fluorine
      Z=ZZ;
      symbol="F";
      name="Fluorine";
      Period=2;
      Group=17;
      Series="Halogen";
      Block="p";
      mass=AMU2KILOGRAM*18.9984;
      MolarVolume=0.011202;
      volume=9.99090;
      Miedema_Vm=NNN;
      valence_std=7;
      valence_iupac=1;
      valence_PT=1;
      Density_PT=16.96E-4;
      crystal="mcl";
      CrystalStr_PT="Base-centered_Monoclinic";
      space_group="C12/c1";
      space_group_number=15;
      Pearson_coefficient=0.0;
      lattice_constant[1]=550;lattice_constant[2]=328;lattice_constant[3]=728;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Gas";
      radius=NNN;
      radius_PT=42;
      radius_covalent=0.57;
      radius_covalent_PT=57;
      radius_VanDerWaals_PT=147;
      radii_Ghosh08=0.4071;
      radii_Slatter=0.50;
      radii_Pyykko=0.64;
      electrical_conductivity=NNN;
      electronegativity_vec=3.98;
      hardness_Ghosh=9.1065;
      electronegativity_Pearson=10.41;
      electronegativity_Ghosh=9.372;
      electron_affinity_PT=328;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=4.00;
      boiling_point=-188.12;
      melting_point=-219.6;
      vaporization_heat_PT=3.27;
      specific_heat_PT=824;
      critical_Pressure=51.04;
      critical_Temperature_PT=144.13;
      thermal_expansion=NNN;
      thermal_conductivity=0.0277;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=7.01;
      Hardness_Putz=16.16;
      Hardness_RB=7.52;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000195;
      HHIP=1500;
      HHIR=1500;
      /*xray_scatt=NNN;*/
      //F
    }
    // [AFLOW]STOP=Fluorine
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Neon
    // Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon
    if(ZZ==10) { // Neon
      Z=ZZ;
      symbol="Ne";
      name="Neon";
      Period=2;
      Group=18;
      Series="NobleGas";
      Block="p";
      mass=AMU2KILOGRAM*20.179;
      MolarVolume=0.02242;
      volume=19.9052;
      Miedema_Vm=NNN;
      valence_std=0;
      valence_iupac=0;
      valence_PT=0;
      Density_PT=9E-4;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.00082783369;
      lattice_constant[1]=442.9;lattice_constant[2]=442.9;lattice_constant[3]=442.9;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Gas";
      radius=0.160;
      radius_PT=38;
      radius_covalent=0.58;
      radius_covalent_PT=58;
      radius_VanDerWaals_PT=154;
      radii_Ghosh08=0.3618;
      radii_Slatter=NNN;
      radii_Pyykko=0.67;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=10.2303;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=10.459;
      electron_affinity_PT=0;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=-246.08;
      melting_point=-248.59;
      vaporization_heat_PT=1.75;
      specific_heat_PT=1030;
      critical_Pressure=27.24;
      critical_Temperature_PT=44.4;
      thermal_expansion=NNN;
      thermal_conductivity=0.0491;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=17.87;
      Hardness_RB=15.45;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-4.1E-9;
      Volume_Magnetic_Susceptibility=-3.69E-9;
      Molar_Magnetic_Susceptibility=-8.27E-11;
      Curie_point=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000067;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Ne volume calculated with fcc-pawpbe
    }
    // [AFLOW]STOP=Neon
    // ********************************************************************************************************************************************************

    // ROW3
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Sodium
    // Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium
    if(ZZ==11) { // Sodium
      Z=ZZ;
      symbol="Na";
      name="Sodium";
      Period=3;
      Group=1;
      Series="AlkaliMetal";
      Block="s";
      mass=AMU2KILOGRAM*22.9898;
      MolarVolume=0.00002375;
      volume=36.9135;
      Miedema_Vm=8.3;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      Density_PT=0.968;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=0.0;
      lattice_constant[1]=429.06;lattice_constant[2]=429.06;lattice_constant[3]=429.06;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.186;
      radius_PT=190;
      radius_covalent=1.66;
      radius_covalent_PT=166;
      radius_VanDerWaals_PT=227;
      radii_Ghosh08=2.1650;
      radii_Slatter=1.80;
      radii_Pyykko=1.55;
      electrical_conductivity=2.1E7;
      electronegativity_vec=0.93;
      hardness_Ghosh=2.4441;
      electronegativity_Pearson=2.85;
      electronegativity_Ghosh=2.536;
      electron_affinity_PT=52.8;
      Miedema_phi_star=2.70;
      Miedema_nws=0.82;
      Miedema_gamma_s=260;
      Pettifor_scale=0.40;
      boiling_point=883;
      melting_point=97.72;
      vaporization_heat_PT=97.7;
      specific_heat_PT=1230;
      critical_Pressure=345.4;
      critical_Temperature_PT=2573;
      thermal_expansion=0.00007;
      thermal_conductivity=140;
      Brinelll_hardness=0.69;
      Mohs_hardness=0.5;
      Vickers_hardness=NNN;
      Hardness_Pearson=2.30;
      Hardness_Putz=0.66;
      Hardness_RB=2.91;
      shear_modulus=3.3;
      Young_modulus=10;
      bulk_modulus=6.3;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=1.6;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=8.8E-9;
      Volume_Magnetic_Susceptibility=8.6E-6;
      Molar_Magnetic_Susceptibility=2E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1100;
      HHIR=500;
      /*xray_scatt=NNN;*/
      // Na
    }
    // [AFLOW]STOP=Sodium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Magnesium
    // Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium
    if(ZZ==12) { // Magnesium
      Z=ZZ;
      symbol="Mg";
      name="Magnesium";
      Period=3;
      Group=2;
      Series="AlkalineEarthMetal";
      Block="s";
      mass=AMU2KILOGRAM*24.305;
      MolarVolume=0.000013984;
      volume=22.8178;
      Miedema_Vm=5.8;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      Density_PT=1.738;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.00073988271;
      lattice_constant[1]=320.94;lattice_constant[2]=320.94;lattice_constant[3]=521.08;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.160;
      radius_PT=145;
      radius_covalent=1.41;
      radius_covalent_PT=141;
      radius_VanDerWaals_PT=173;
      radii_Ghosh08=1.6711;
      radii_Slatter=1.50;
      radii_Pyykko=1.39;
      electrical_conductivity=2.3E7;
      electronegativity_vec=1.31;
      hardness_Ghosh=3.0146;
      electronegativity_Pearson=3.75;
      electronegativity_Ghosh=3.310;
      electron_affinity_PT=0;
      Miedema_phi_star=3.45;
      Miedema_nws=1.17;
      Miedema_gamma_s=790;
      Pettifor_scale=1.28;
      boiling_point=1090;
      melting_point=650;
      vaporization_heat_PT=128;
      specific_heat_PT=1020;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000248;
      thermal_conductivity=160;
      Brinelll_hardness=260;
      Mohs_hardness=2.5;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.90;
      Hardness_Putz=0.93;
      Hardness_RB=4.63;
      shear_modulus=17;
      Young_modulus=45;
      bulk_modulus=45;
      Poisson_ratio_PT=0.29;
      Miedema_BVm=5.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=6.9E-9;
      Volume_Magnetic_Susceptibility=0.000012;
      Molar_Magnetic_Susceptibility=1.68E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=5300;
      HHIR=500;
      /*xray_scatt=NNN;*/
      //Mg
    }
    // [AFLOW]STOP=Magnesium
    // ********************************************************************************************************************************************************

    // p-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Aluminium
    //Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium
    if(ZZ==13) { // Aluminium
      Z=ZZ;
      symbol="Al";
      name="Aluminium";
      Period=3;
      Group=13;
      Series="PoorMetal";
      Block="p";
      mass=AMU2KILOGRAM*26.9815;
      MolarVolume=9.99E-6;
      volume=16.4000;
      Miedema_Vm=4.6;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=2.7;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.0;
      lattice_constant[1]=404.95;lattice_constant[2]=404.95;lattice_constant[3]=404.95;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.143;
      radius_PT=118;
      radius_covalent=1.21;
      radius_covalent_PT=121;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.3608;
      radii_Slatter=1.25;
      radii_Pyykko=1.26;
      electrical_conductivity=3.8E7;
      electronegativity_vec=1.61;
      hardness_Ghosh=3.5849;
      electronegativity_Pearson=3.23;
      electronegativity_Ghosh=4.084;
      electron_affinity_PT=42.5;
      Miedema_phi_star=4.20;
      Miedema_nws=1.39;
      Miedema_gamma_s=1200;
      Pettifor_scale=1.66;
      boiling_point=2519;
      melting_point=660.32;
      vaporization_heat_PT=293;
      specific_heat_PT=904;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000231;
      thermal_conductivity=235;
      Brinelll_hardness=245;
      Mohs_hardness=2.75;
      Vickers_hardness=167;
      Hardness_Pearson=2.77;
      Hardness_Putz=1.42;
      Hardness_RB=2.94;
      shear_modulus=26;
      Young_modulus=70;
      bulk_modulus=76;
      Poisson_ratio_PT=0.35;
      Miedema_BVm=7.2;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=7.8E-9;
      Volume_Magnetic_Susceptibility=0.0000211;
      Molar_Magnetic_Susceptibility=2.1E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1600;
      HHIR=1000;
      /*xray_scatt=NNN;*/
      //Al
    }
    // [AFLOW]STOP=Aluminium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Silicon
    // Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon
    if(ZZ==14) { // Silicon
      Z=ZZ;
      symbol="Si";
      name="Silicon";
      Period=3;
      Group=14;
      Series="Metalloid";
      Block="p";
      mass=AMU2KILOGRAM*28.0855;
      MolarVolume=0.000012054;
      volume=14.3536;
      Miedema_Vm=4.2;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=2.33;
      crystal="dia";
      CrystalStr_PT="Tetrahedral_Packing";
      space_group="Fd_3m";
      space_group_number=227;
      Pearson_coefficient=0.00020046752;
      lattice_constant[1]=543.09;lattice_constant[2]=543.09;lattice_constant[3]=543.09;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.117;
      radius_PT=111;
      radius_covalent=1.11;
      radius_covalent_PT=111;
      radius_VanDerWaals_PT=210;
      radii_Ghosh08=1.1477;
      radii_Slatter=1.10;
      radii_Pyykko=1.16;
      electrical_conductivity=1000;
      electronegativity_vec=1.90;
      hardness_Ghosh=4.1551;
      electronegativity_Pearson=4.77;
      electronegativity_Ghosh=4.857;
      electron_affinity_PT=133.6;
      Miedema_phi_star=4.70;
      Miedema_nws=1.50;
      Miedema_gamma_s=1290;
      Pettifor_scale=1.92;
      boiling_point=2900;
      melting_point=1414;
      vaporization_heat_PT=359;
      specific_heat_PT=710;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=2.6E-6;
      thermal_conductivity=150;
      Brinelll_hardness=NNN;
      Mohs_hardness=6.5;
      Vickers_hardness=9630.1303;
      Hardness_Pearson=3.38;
      Hardness_Putz=2.10;
      Hardness_RB=3.61;
      shear_modulus=NNN;
      Young_modulus=47;
      bulk_modulus=100;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=11.9;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-1.6E-9;
      Volume_Magnetic_Susceptibility=-3.73E-6;
      Molar_Magnetic_Susceptibility=-4.49E-11;
      Curie_point=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=4700;
      HHIR=1000;
      xray_scatt=14.43;
      //Si ???
    }
    // [AFLOW]STOP=Silicon
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Phosphorus
    // Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus
    if(ZZ==15) { // Phosphorus
      Z=ZZ;
      symbol="P";
      name="Phosphorus";
      Period=3;
      Group=15;
      Series="Nonmetal";
      Block="p";
      mass=AMU2KILOGRAM*30.9738;
      MolarVolume=0.000016991;
      volume=14.1995;
      Miedema_Vm=NNN;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      Density_PT=1.823;
      crystal="cub";
      CrystalStr_PT="Simple_Triclinic";
      space_group="P-1";
      space_group_number=2;
      Pearson_coefficient=0.0;
      lattice_constant[1]=1145;lattice_constant[2]=550.3;lattice_constant[3]=1126.1;
      lattice_angle[1]=1.25384;lattice_angle[2]=1.57725;lattice_angle[3]=1.24896;
      phase="Solid";
      radius=0.109;
      radius_PT=98;
      radius_covalent=1.07;
      radius_covalent_PT=107;
      radius_VanDerWaals_PT=180;
      radii_Ghosh08=0.9922;
      radii_Slatter=1.00;
      radii_Pyykko=1.11;
      electrical_conductivity=1E7;
      electronegativity_vec=2.19;
      hardness_Ghosh=4.7258;
      electronegativity_Pearson=5.62;
      electronegativity_Ghosh=5.631;
      electron_affinity_PT=71;
      Miedema_phi_star=5.5;
      Miedema_nws=1.65;
      Miedema_gamma_s=NNN;
      Pettifor_scale=2.18;
      boiling_point=280.5;
      melting_point=44.2;
      vaporization_heat_PT=12.4;
      specific_heat_PT=769.7;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=0.236;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=4.88;
      Hardness_Putz=2.92;
      Hardness_RB=5.42;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=11;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-1.13E-8;
      Volume_Magnetic_Susceptibility=-0.0000206;
      Molar_Magnetic_Susceptibility=-3.5E-10;
      Curie_point=NNN;
      color_PT="COLORLESS";
      refractive_index=1.001212;
      HHIP=2000;
      HHIR=5100;
      xray_scatt=15.3133;
      //P MIEDEMA =PAUL VAN DER PUT book
    }
    // [AFLOW]STOP=Phosphorus
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Sulphur
    // Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur
    if(ZZ==16) { // Sulphur
      Z=ZZ;
      symbol="S";
      name="Sulphur";
      Period=3;
      Group=16;
      Series="Chalcogen";
      Block="p";
      mass=AMU2KILOGRAM*32.06;
      MolarVolume=0.000016357;
      volume=15.7301;
      Miedema_Vm=4.376;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      Density_PT=1.96;
      crystal="orc";
      CrystalStr_PT="Face-centered_Orthorhombic";
      space_group="Fddd";
      space_group_number=70;
      Pearson_coefficient=0.00016807795;
      lattice_constant[1]=1043.7;lattice_constant[2]=1284.5;lattice_constant[3]=2436.9;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.106;
      radius_PT=87;
      radius_covalent=1.05;
      radius_covalent_PT=105;
      radius_VanDerWaals_PT=180;
      radii_Ghosh08=0.8739;
      radii_Slatter=1.00;
      radii_Pyykko=1.03;
      electrical_conductivity=1E-15;
      electronegativity_vec=2.58;
      hardness_Ghosh=5.2960;
      electronegativity_Pearson=6.22;
      electronegativity_Ghosh=6.420;
      electron_affinity_PT=200;
      Miedema_phi_star=5.6;
      Miedema_nws=1.46;
      Miedema_gamma_s=NNN;
      Pettifor_scale=2.44;
      boiling_point=444.72;
      melting_point=115.21;
      vaporization_heat_PT=9.8;
      specific_heat_PT=705;
      critical_Pressure=204.3;
      critical_Temperature_PT=1314;
      thermal_expansion=NNN;
      thermal_conductivity=0.205;
      Brinelll_hardness=NNN;
      Mohs_hardness=2;
      Vickers_hardness=NNN;
      Hardness_Pearson=4.14;
      Hardness_Putz=3.82;
      Hardness_RB=4.28;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=7.7;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-6.2E-9;
      Volume_Magnetic_Susceptibility=-0.0000122;
      Molar_Magnetic_Susceptibility=-1.99E-10;
      Curie_point=NNN;
      color_PT="YELLOW";
      refractive_index=1.001111;
      HHIP=700;
      HHIR=1000;
      /*xray_scatt=NNN;*/
      //S Table 27 of JUNKAI
    }
    // [AFLOW]STOP=Sulphur
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Chlorine
    // Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine
    if(ZZ==17) { // Chlorine
      Z=ZZ;
      symbol="Cl";
      name="Chlorine";
      Period=3;
      Group=17;
      Series="Halogen";
      Block="p";
      mass=AMU2KILOGRAM*35.453;
      MolarVolume=0.01103;
      volume=21.2947;
      Miedema_Vm=6.71;
      valence_std=7;
      valence_iupac=7;
      valence_PT=5;
      Density_PT=32.14E-4;
      crystal="orc";
      CrystalStr_PT="Base_Orthorhombic";
      space_group="Cmca";
      space_group_number=64;
      Pearson_coefficient=0.00058238731;
      lattice_constant[1]=622.35;lattice_constant[2]=445.61;lattice_constant[3]=817.85;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Gas";
      radius=0.107;
      radius_PT=79;
      radius_covalent=1.02;
      radius_covalent_PT=102;
      radius_VanDerWaals_PT=175;
      radii_Ghosh08=0.7808;
      radii_Slatter=1.00;
      radii_Pyykko=0.99;
      electrical_conductivity=0.01;
      electronegativity_vec=3.16;
      hardness_Ghosh=5.8662;
      electronegativity_Pearson=8.30;
      electronegativity_Ghosh=7.178;
      electron_affinity_PT=349;
      Miedema_phi_star=5.32;
      Miedema_nws=0.34;
      Miedema_gamma_s=1013;
      Pettifor_scale=2.70;
      boiling_point=-34.04;
      melting_point=-101.5;
      vaporization_heat_PT=10.2;
      specific_heat_PT=478.2;
      critical_Pressure=78.87;
      critical_Temperature_PT=416.9;
      thermal_expansion=NNN;
      thermal_conductivity=0.0089;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=4.68;
      Hardness_Putz=5.01;
      Hardness_RB=4.91;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=1.1;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-7.2E-9;
      Volume_Magnetic_Susceptibility=-2.31E-8;
      Molar_Magnetic_Susceptibility=-5.11E-10;
      Curie_point=NNN;
      color_PT="YELLOW";
      refractive_index=1.000773;
      HHIP=1500;
      HHIR=1500;
      /*xray_scatt=NNN;*/
      //Cl interpolation phi_star, nws, Vm, gamma JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Chlorine
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Argon
    //Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon
    if(ZZ==18) { // Argon
      Z=ZZ;
      symbol="Ar";
      name="Argon";
      Period=3;
      Group=18;
      Series="NobleGas";
      Block="p";
      mass=AMU2KILOGRAM*39.948;
      MolarVolume=0.022392;
      volume=22.000;
      Miedema_Vm=NNN;
      valence_std=0;
      valence_iupac=2;
      valence_PT=0;
      Density_PT=17.84E-4;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.00003509919;
      lattice_constant[1]=525.6;lattice_constant[2]=525.6;lattice_constant[3]=525.6;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Gas";
      radius=0.192;
      radius_PT=71;
      radius_covalent=1.06;
      radius_covalent_PT=106;
      radius_VanDerWaals_PT=188;
      radii_Ghosh08=0.7056;
      radii_Slatter=NNN;
      radii_Pyykko=0.96;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=6.4366;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=7.951;
      electron_affinity_PT=0;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=-185.8;
      melting_point=-189.3;
      vaporization_heat_PT=6.5;
      specific_heat_PT=520.33;
      critical_Pressure=48.34;
      critical_Temperature_PT=150.87;
      thermal_expansion=NNN;
      thermal_conductivity=0.01772;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=6.16;
      Hardness_RB=10.69;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-6E-9;
      Volume_Magnetic_Susceptibility=-1.07E-8;
      Molar_Magnetic_Susceptibility=-2.4E-10;
      Curie_point=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000281;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Ar guessed volume, must double check from results JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Argon
    // ********************************************************************************************************************************************************

    // ROW4
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Potassium
    // Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium
    if(ZZ==19) { // Potassium
      Z=ZZ;
      symbol="K";
      name="Potassium";
      Period=4;
      Group=1;
      Series="AlkaliMetal";
      Block="s";
      mass=AMU2KILOGRAM*39.0983;
      MolarVolume=0.00004568;
      volume=73.9091;
      Miedema_Vm=12.8;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      Density_PT=0.856;
      crystal="fcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=0.000164;
      lattice_constant[1]=532.8;lattice_constant[2]=532.8;lattice_constant[3]=532.8;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.231;
      radius_PT=243;
      radius_covalent=2.03;
      radius_covalent_PT=203;
      radius_VanDerWaals_PT=275;
      radii_Ghosh08=3.2930;
      radii_Slatter=2.20;
      radii_Pyykko=1.96;
      electrical_conductivity=1.4E7;
      electronegativity_vec=0.82;
      hardness_Ghosh=2.3273;
      electronegativity_Pearson=2.42;
      electronegativity_Ghosh=2.672;
      electron_affinity_PT=48.4;
      Miedema_phi_star=2.25;
      Miedema_nws=0.65;
      Miedema_gamma_s=150;
      Pettifor_scale=0.35;
      boiling_point=759;
      melting_point=63.38;
      vaporization_heat_PT=76.9;
      specific_heat_PT=757;
      critical_Pressure=157.9;
      critical_Temperature_PT=2223;
      thermal_expansion=NNN;
      thermal_conductivity=100;
      Brinelll_hardness=0.363;
      Mohs_hardness=0.4;
      Vickers_hardness=NNN;
      Hardness_Pearson=1.92;
      Hardness_Putz=0.18;
      Hardness_RB=2.35;
      shear_modulus=1.3;
      Young_modulus=NNN;
      bulk_modulus=3.1;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=1.5;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=6.7E-9;
      Volume_Magnetic_Susceptibility=5.74E-6;
      Molar_Magnetic_Susceptibility=2.62E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1700;
      HHIR=7200;
      /*xray_scatt=NNN;*/
      //K
    }
    // [AFLOW]STOP=Potassium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Calcium
    // Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium
    if(ZZ==20) { // Calcium
      Z=ZZ;
      symbol="Ca";
      name="Calcium";
      Period=4;
      Group=2;
      Series="AlkalineEarthMetal";
      Block="s";
      mass=AMU2KILOGRAM*40.08;
      MolarVolume=0.000025857;
      volume=42.1927;
      Miedema_Vm=8.8;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      Density_PT=1.55;
      crystal="bcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.000297564;
      lattice_constant[1]=558.84;lattice_constant[2]=558.84;lattice_constant[3]=558.84;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.197;
      radius_PT=194;
      radius_covalent=1.76;
      radius_covalent_PT=176;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.5419;
      radii_Slatter=1.80;
      radii_Pyykko=1.71;
      electrical_conductivity=2.9E7;
      electronegativity_vec=1.00;
      hardness_Ghosh=2.7587;
      electronegativity_Pearson=2.2;
      electronegativity_Ghosh=3.140;
      electron_affinity_PT=2.37;
      Miedema_phi_star=2.55;
      Miedema_nws=0.91;
      Miedema_gamma_s=490;
      Pettifor_scale=0.60;
      boiling_point=1484;
      melting_point=842;
      vaporization_heat_PT=155;
      specific_heat_PT=631;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000223;
      thermal_conductivity=200;
      Brinelll_hardness=167;
      Mohs_hardness=1.75;
      Vickers_hardness=NNN;
      Hardness_Pearson=4.00;
      Hardness_Putz=0.25;
      Hardness_RB=3.07;
      shear_modulus=7.4;
      Young_modulus=20;
      bulk_modulus=17;
      Poisson_ratio_PT=0.31;
      Miedema_BVm=4.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.38E-8;
      Volume_Magnetic_Susceptibility=0.00002139;
      Molar_Magnetic_Susceptibility=5.531E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3900;
      HHIR=1500;
      /*xray_scatt=NNN;*/
      //Ca
    }
    // [AFLOW]STOP=Calcium
    // ********************************************************************************************************************************************************

    // d-electron systems: transition metals
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Scandium
    // Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium
    if(ZZ==21) { // Scandium
      Z=ZZ;
      symbol="Sc";
      name="Scandium";
      Period=4;
      Group=3;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*44.9559;
      MolarVolume=0.000015061;
      volume=24.6739;
      Miedema_Vm=6.1;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=2.985;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=330.9;lattice_constant[2]=330.9;lattice_constant[3]=527.33;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.160;
      radius_PT=184;
      radius_covalent=1.70;
      radius_covalent_PT=170;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.4149;
      radii_Slatter=1.60;
      radii_Pyykko=1.48;
      electrical_conductivity=1.8E6;
      electronegativity_vec=1.36;
      hardness_Ghosh=2.8582;
      electronegativity_Pearson=3.34;
      electronegativity_Ghosh=3.248;
      electron_affinity_PT=18.1;
      Miedema_phi_star=3.25;
      Miedema_nws=1.27;
      Miedema_gamma_s=1200;
      Pettifor_scale=0.74;
      boiling_point=2830;
      melting_point=1541;
      vaporization_heat_PT=318;
      specific_heat_PT=567;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000102;
      thermal_conductivity=16;
      Brinelll_hardness=750;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.20;
      Hardness_Putz=0.31;
      Hardness_RB=2.52;
      shear_modulus=29;
      Young_modulus=74;
      bulk_modulus=57;
      Poisson_ratio_PT=0.28;
      Miedema_BVm=6.6;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=8.8E-8;
      Volume_Magnetic_Susceptibility=0.0002627;
      Molar_Magnetic_Susceptibility=3.956E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=5500;
      HHIR=4500;
      xray_scatt=21.34;
      //Sc
    }
    // [AFLOW]STOP=Scandium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Titanium
    // Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium
    if(ZZ==22) { // Titanium
      Z=ZZ;
      symbol="Ti";
      name="Titanium";
      Period=4;
      Group=4;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*47.9;
      MolarVolume=0.000010621;
      volume=17.1035;
      Miedema_Vm=4.8;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=4.507;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.000286456;
      lattice_constant[1]=295.08;lattice_constant[2]=295.08;lattice_constant[3]=468.55;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.147;
      radius_PT=176;
      radius_covalent=1.60;
      radius_covalent_PT=160;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.2998;
      radii_Slatter=1.40;
      radii_Pyykko=1.36;
      electrical_conductivity=2.5E6;
      electronegativity_vec=1.54;
      hardness_Ghosh=2.9578;
      electronegativity_Pearson=3.45;
      electronegativity_Ghosh=3.357;
      electron_affinity_PT=7.6;
      Miedema_phi_star=3.65;
      Miedema_nws=1.47;
      Miedema_gamma_s=2050;
      Pettifor_scale=0.79;
      boiling_point=3287;
      melting_point=1668;
      vaporization_heat_PT=425;
      specific_heat_PT=520;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=8.6E-6;
      thermal_conductivity=22;
      Brinelll_hardness=715;
      Mohs_hardness=6;
      Vickers_hardness=970;
      Hardness_Pearson=3.37;
      Hardness_Putz=0.38;
      Hardness_RB=2.03;
      shear_modulus=44;
      Young_modulus=116;
      bulk_modulus=110;
      Poisson_ratio_PT=0.32;
      Miedema_BVm=11.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=4.01E-8;
      Volume_Magnetic_Susceptibility=0.0001807;
      Molar_Magnetic_Susceptibility=1.919E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1100;
      HHIR=1600;
      xray_scatt=22.24;
      //Ti
    }
    // [AFLOW]STOP=Titanium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Vanadium
    // Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium
    if(ZZ==23) { // Vanadium
      Z=ZZ;
      symbol="V";
      name="Vanadium";
      Period=4;
      Group=5;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*50.9415;
      MolarVolume=8.3374E-6;
      volume=13.2086;
      Miedema_Vm=4.1;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      Density_PT=6.11;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=9.54831E-07;
      lattice_constant[1]=303;lattice_constant[2]=303;lattice_constant[3]=303;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.132;
      radius_PT=171;
      radius_covalent=1.53;
      radius_covalent_PT=153;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.1953;
      radii_Slatter=1.35;
      radii_Pyykko=1.34;
      electrical_conductivity=5E6;
      electronegativity_vec=1.63;
      hardness_Ghosh=3.0573;
      electronegativity_Pearson=3.6;
      electronegativity_Ghosh=3.465;
      electron_affinity_PT=50.6;
      Miedema_phi_star=4.25;
      Miedema_nws=1.64;
      Miedema_gamma_s=2600;
      Pettifor_scale=0.84;
      boiling_point=3407;
      melting_point=1910;
      vaporization_heat_PT=453;
      specific_heat_PT=489;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=8.4E-6;
      thermal_conductivity=31;
      Brinelll_hardness=628;
      Mohs_hardness=7;
      Vickers_hardness=628;
      Hardness_Pearson=3.10;
      Hardness_Putz=0.45;
      Hardness_RB=NNN;
      shear_modulus=47;
      Young_modulus=128;
      bulk_modulus=160;
      Poisson_ratio_PT=0.37;
      Miedema_BVm=14.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=6.28E-8;
      Volume_Magnetic_Susceptibility=0.0003837;
      Molar_Magnetic_Susceptibility=3.199E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3300;
      HHIR=3400;
      /*xray_scatt=NNN;*/
      //V
    }
    // [AFLOW]STOP=Vanadium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Chromium
    // Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium
    if(ZZ==24) { // Chromium
      Z=ZZ;
      symbol="Cr";
      name="Chromium";
      Period=4;
      Group=6;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*51.996;
      MolarVolume=7.2317E-6;
      volume=11.4136;
      Miedema_Vm=3.7;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      Density_PT=7.19;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=0.00013287;
      lattice_constant[1]=291;lattice_constant[2]=291;lattice_constant[3]=291;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.125;
      radius_PT=166;
      radius_covalent=1.39;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.1000;
      radii_Slatter=1.40;
      radii_Pyykko=1.22;
      electrical_conductivity=7.9E6;
      electronegativity_vec=1.66;
      hardness_Ghosh=3.1567;
      electronegativity_Pearson=3.72;
      electronegativity_Ghosh=3.573;
      electron_affinity_PT=64.3;
      Miedema_phi_star=4.65;
      Miedema_nws=1.74;
      Miedema_gamma_s=2400;
      Pettifor_scale=0.89;
      boiling_point=2671;
      melting_point=1907;
      vaporization_heat_PT=339;
      specific_heat_PT=448;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=4.9E-6;
      thermal_conductivity=94;
      Brinelll_hardness=1120;
      Mohs_hardness=8.5;
      Vickers_hardness=1060;
      Hardness_Pearson=3.06;
      Hardness_Putz=0.54;
      Hardness_RB=4.06;
      shear_modulus=115;
      Young_modulus=279;
      bulk_modulus=160;
      Poisson_ratio_PT=0.21;
      Miedema_BVm=14.0;
      Magnetic_Type_PT="Antiferromagnetic";
      Mass_Magnetic_Susceptibility=4.45E-8;
      Volume_Magnetic_Susceptibility=0.0003177;
      Molar_Magnetic_Susceptibility=2.314E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3100;
      HHIR=4100;
      xray_scatt=23.84;
      //Cr
    }
    // [AFLOW]STOP=Chromium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Manganese
    // Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese
    if(ZZ==25) { // Manganese
      Z=ZZ;
      symbol="Mn";
      name="Manganese";
      Period=4;
      Group=7;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*54.93805;
      MolarVolume=7.3545E-6;
      volume=10.6487;
      Miedema_Vm=3.8;
      valence_std=7;
      valence_iupac=7;
      valence_PT=4;
      Density_PT=7.47;
      crystal="cub";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="I_43m";
      space_group_number=217;
      Pearson_coefficient=1.67276E-32;
      lattice_constant[1]=891.25;lattice_constant[2]=891.25;lattice_constant[3]=891.25;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.112;
      radius_PT=161;
      radius_covalent=1.61;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.0124;
      radii_Slatter=1.40;
      radii_Pyykko=1.19;
      electrical_conductivity=620000;
      electronegativity_vec=1.55;
      hardness_Ghosh=3.2564;
      electronegativity_Pearson=3.72;
      electronegativity_Ghosh=3.681;
      electron_affinity_PT=0;
      Miedema_phi_star=4.45;
      Miedema_nws=1.61;
      Miedema_gamma_s=1600;
      Pettifor_scale=0.94;
      boiling_point=2061;
      melting_point=1246;
      vaporization_heat_PT=220;
      specific_heat_PT=479;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000217;
      thermal_conductivity=7.7;
      Brinelll_hardness=196;
      Mohs_hardness=6;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.72;
      Hardness_Putz=0.64;
      Hardness_RB=2.88;
      shear_modulus=NNN;
      Young_modulus=198;
      bulk_modulus=120;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=4.4;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.21E-7;
      Volume_Magnetic_Susceptibility=0.00090387;
      Molar_Magnetic_Susceptibility=6.6475E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1600;
      HHIR=1800;
      xray_scatt=24.46;
      //xray_scatt=24.3589; Mn JUNKAI CHANGED VALENCE// DX and CO-20170904 radius_covalent[i] uses high spin configuration (most frequent)
    }
    // [AFLOW]STOP=Manganese
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Iron
    // Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron
    if(ZZ==26) { // Iron
      Z=ZZ;
      symbol="Fe";
      name="Iron";
      Period=4;
      Group=8;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*55.847;
      MolarVolume=7.0923E-6;
      volume=10.2315;
      Miedema_Vm=3.7;
      valence_std=8;
      valence_iupac=6;
      valence_PT=3;
      Density_PT=7.874;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=9.17912E-05;
      lattice_constant[1]=286.65;lattice_constant[2]=286.65;lattice_constant[3]=286.65;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.124;
      radius_PT=156;
      radius_covalent=1.52;
      radius_covalent_PT=132;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.9319;
      radii_Slatter=1.40;
      radii_Pyykko=1.16;
      electrical_conductivity=1E7;
      electronegativity_vec=1.83;
      hardness_Ghosh=3.3559;
      electronegativity_Pearson=4.06;
      electronegativity_Ghosh=3.789;
      electron_affinity_PT=15.7;
      Miedema_phi_star=4.93;
      Miedema_nws=1.77;
      Miedema_gamma_s=2550;
      Pettifor_scale=0.99;
      boiling_point=2861;
      melting_point=1538;
      vaporization_heat_PT=347;
      specific_heat_PT=449;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000118;
      thermal_conductivity=79;
      Brinelll_hardness=490;
      Mohs_hardness=4;
      Vickers_hardness=608;
      Hardness_Pearson=3.81;
      Hardness_Putz=0.75;
      Hardness_RB=2.53;
      shear_modulus=82;
      Young_modulus=211;
      bulk_modulus=170;
      Poisson_ratio_PT=0.29;
      Miedema_BVm=12.0;
      Magnetic_Type_PT="Ferromagnetic";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=1043;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=2400;
      HHIR=1400;
      xray_scatt=24.85;
      //xray_scatt=24.6830; Fe JUNKAI CHANGED VALENCE// DX and CO -20170904 radius_covalent[i] uses high spin configuration (most frequent)
    }
    // [AFLOW]STOP=Iron
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Cobalt
    // Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt
    if(ZZ==27) { // Cobalt
      Z=ZZ;
      symbol="Co";
      name="Cobalt";
      Period=4;
      Group=9;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*58.9332;
      MolarVolume=6.62E-6;
      volume=10.3205;
      Miedema_Vm=3.5;
      valence_std=9;
      valence_iupac=5;
      valence_PT=4;
      Density_PT=8.9;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=250.71;lattice_constant[2]=250.71;lattice_constant[3]=406.95;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.125;
      radius_PT=152;
      radius_covalent=1.26;
      radius_covalent_PT=126;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.8575;
      radii_Slatter=1.35;
      radii_Pyykko=1.11;
      electrical_conductivity=1.7E7;
      electronegativity_vec=1.88;
      hardness_Ghosh=3.4556;
      electronegativity_Pearson=4.3;
      electronegativity_Ghosh=3.897;
      electron_affinity_PT=63.7;
      Miedema_phi_star=5.10;
      Miedema_nws=1.75;
      Miedema_gamma_s=2550;
      Pettifor_scale=1.04;
      boiling_point=2927;
      melting_point=1495;
      vaporization_heat_PT=375;
      specific_heat_PT=421;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000013;
      thermal_conductivity=100;
      Brinelll_hardness=700;
      Mohs_hardness=5;
      Vickers_hardness=1043;
      Hardness_Pearson=3.60;
      Hardness_Putz=0.88;
      Hardness_RB=3.53;
      shear_modulus=76;
      Young_modulus=209;
      bulk_modulus=180;
      Poisson_ratio_PT=0.31;
      Miedema_BVm=13.0;
      Magnetic_Type_PT="Ferromagnetic";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=1394;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=3100;
      HHIR=2700;
      xray_scatt=24.59;
      //Co JUNKAI CHANGED VALENCE// DX and CO -20170904 radius_covalent[i] uses low spin configuration (most frequent)
    }
    // [AFLOW]STOP=Cobalt
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Nickel
    // Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel
    if(ZZ==28) { // Nickel
      Z=ZZ;
      symbol="Ni";
      name="Nickel";
      Period=4;
      Group=10;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*58.69;
      MolarVolume=6.5888E-6;
      volume=10.8664;
      Miedema_Vm=3.5;
      valence_std=10;
      valence_iupac=4;
      valence_PT=2;
      Density_PT=8.908;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.000430773;
      lattice_constant[1]=352.4;lattice_constant[2]=352.4;lattice_constant[3]=352.4;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.125;
      radius_PT=149;
      radius_covalent=1.24;
      radius_covalent_PT=124;
      radius_VanDerWaals_PT=163;
      radii_Ghosh08=1.7888;
      radii_Slatter=1.35;
      radii_Pyykko=1.10;
      electrical_conductivity=1.4E7;
      electronegativity_vec=1.91;
      hardness_Ghosh=3.5550;
      electronegativity_Pearson=4.40;
      electronegativity_Ghosh=4.005;
      electron_affinity_PT=112;
      Miedema_phi_star=5.20;
      Miedema_nws=1.75;
      Miedema_gamma_s=2450;
      Pettifor_scale=1.09;
      boiling_point=2913;
      melting_point=1455;
      vaporization_heat_PT=378;
      specific_heat_PT=445;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000134;
      thermal_conductivity=91;
      Brinelll_hardness=700;
      Mohs_hardness=4;
      Vickers_hardness=638;
      Hardness_Pearson=3.25;
      Hardness_Putz=1.02;
      Hardness_RB=4.08;
      shear_modulus=76;
      Young_modulus=200;
      bulk_modulus=180;
      Poisson_ratio_PT=0.31;
      Miedema_BVm=12.0;
      Magnetic_Type_PT="Ferromagnetic";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=631;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=1000;
      HHIR=1500;
      xray_scatt=25.02;
      //Ni
    }
    // [AFLOW]STOP=Nickel
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Copper
    // Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper
    if(ZZ==29) { // Copper
      Z=ZZ;
      symbol="Cu";
      name="Copper";
      Period=4;
      Group=11;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*63.546;
      MolarVolume=7.0922E-6;
      volume=12.0159;
      Miedema_Vm=3.7;
      valence_std=11;
      valence_iupac=4;
      valence_PT=2;
      Density_PT=8.96;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.00021086;
      lattice_constant[1]=361.49;lattice_constant[2]=361.49;lattice_constant[3]=361.49;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.128;
      radius_PT=145;
      radius_covalent=1.32;
      radius_covalent_PT=132;
      radius_VanDerWaals_PT=140;
      radii_Ghosh08=1.725;
      radii_Slatter=1.35;
      radii_Pyykko=1.12;
      electrical_conductivity=5.9E7;
      electronegativity_vec=1.90;
      hardness_Ghosh=3.6544;
      electronegativity_Pearson=4.48;
      electronegativity_Ghosh=4.113;
      electron_affinity_PT=118.4;
      Miedema_phi_star=4.55;
      Miedema_nws=1.47;
      Miedema_gamma_s=1850;
      Pettifor_scale=1.20;
      boiling_point=2562;
      melting_point=1084.62;
      vaporization_heat_PT=300;
      specific_heat_PT=384.4;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000165;
      thermal_conductivity=400;
      Brinelll_hardness=874;
      Mohs_hardness=3;
      Vickers_hardness=369;
      Hardness_Pearson=3.25;
      Hardness_Putz=1.21;
      Hardness_RB=NNN;
      shear_modulus=48;
      Young_modulus=130;
      bulk_modulus=140;
      Poisson_ratio_PT=0.34;
      Miedema_BVm=9.3;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-1.08E-9;
      Volume_Magnetic_Susceptibility=-9.63E-6;
      Molar_Magnetic_Susceptibility=-6.86E-11;
      Curie_point=NNN;
      color_PT="COPPER";
      refractive_index=NNN;
      HHIP=1600;
      HHIR=1500;
      xray_scatt=27.03;
      //Cu JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Copper
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Zinc
    // Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc
    if(ZZ==30) { // Zinc
      Z=ZZ;
      symbol="Zn";
      name="Zinc";
      Period=4;
      Group=12;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*65.38;
      MolarVolume=9.157E-6;
      volume=15.0827;
      Miedema_Vm=4.4;
      valence_std=12;
      valence_iupac=2;
      valence_PT=2;
      Density_PT=7.14;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.000595597;
      lattice_constant[1]=266.49;lattice_constant[2]=266.49;lattice_constant[3]=494.68;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.133;
      radius_PT=142;
      radius_covalent=1.22;
      radius_covalent_PT=122;
      radius_VanDerWaals_PT=139;
      radii_Ghosh08=1.6654;
      radii_Slatter=1.35;
      radii_Pyykko=1.18;
      electrical_conductivity=1.7E7;
      electronegativity_vec=1.65;
      hardness_Ghosh=3.7542;
      electronegativity_Pearson=4.45;
      electronegativity_Ghosh=4.222;
      electron_affinity_PT=0;
      Miedema_phi_star=4.10;
      Miedema_nws=1.32;
      Miedema_gamma_s=1020;
      Pettifor_scale=1.44;
      boiling_point=907;
      melting_point=419.53;
      vaporization_heat_PT=119;
      specific_heat_PT=388;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000302;
      thermal_conductivity=120;
      Brinelll_hardness=412;
      Mohs_hardness=2.5;
      Vickers_hardness=NNN;
      Hardness_Pearson=4.94;
      Hardness_Putz=1.39;
      Hardness_RB=6.01;
      shear_modulus=43;
      Young_modulus=108;
      bulk_modulus=70;
      Poisson_ratio_PT=0.25;
      Miedema_BVm=5.5;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-2.21E-9;
      Volume_Magnetic_Susceptibility=-0.0000158;
      Molar_Magnetic_Susceptibility=-1.45E-10;
      Curie_point=NNN;
      color_PT="SLATEGRAY";
      refractive_index=1.00205;
      HHIP=1600;
      HHIR=1900;
      xray_scatt=28.44;
      //Zn
    }
    // [AFLOW]STOP=Zinc
    // ********************************************************************************************************************************************************

    // p-electron systems 
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Gallium
    // Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium
    if(ZZ==31) { // Gallium
      Z=ZZ;
      symbol="Ga";
      name="Gallium";
      Period=4;
      Group=13;
      Series="PoorMetal";
      Block="p";
      mass=AMU2KILOGRAM*69.737;
      MolarVolume=0.000011809;
      volume=18.9039;
      Miedema_Vm=5.2;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=5.904;
      crystal="orc";
      CrystalStr_PT="Base_Orthorhombic";
      space_group="Cmca";
      space_group_number=64;
      Pearson_coefficient=0.000197588;
      lattice_constant[1]=451.97;lattice_constant[2]=766.33;lattice_constant[3]=452.6;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.135;
      radius_PT=136;
      radius_covalent=1.22;
      radius_covalent_PT=122;
      radius_VanDerWaals_PT=187;
      radii_Ghosh08=1.4489;
      radii_Slatter=1.30;
      radii_Pyykko=1.24;
      electrical_conductivity=7.1E6;
      electronegativity_vec=1.81;
      hardness_Ghosh=4.1855;
      electronegativity_Pearson=3.2;
      electronegativity_Ghosh=4.690;
      electron_affinity_PT=28.9;
      Miedema_phi_star=4.10;
      Miedema_nws=1.31;
      Miedema_gamma_s=830;
      Pettifor_scale=1.68;
      boiling_point=2204;
      melting_point=29.76;
      vaporization_heat_PT=256;
      specific_heat_PT=371;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.00012;
      thermal_conductivity=29;
      Brinelll_hardness=60;
      Mohs_hardness=1.5;
      Vickers_hardness=NNN;
      Hardness_Pearson=2.90;
      Hardness_Putz=1.59;
      Hardness_RB=3.03;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=6.7;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-3E-9;
      Volume_Magnetic_Susceptibility=-0.0000177;
      Molar_Magnetic_Susceptibility=-2.09E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=5500;
      HHIR=1900;
      /*xray_scatt=NNN;*/
      //Ga
    }
    // [AFLOW]STOP=Gallium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Germanium
    // Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium
    if(ZZ==32) { // Germanium
      Z=ZZ;
      symbol="Ge";
      name="Germanium";
      Period=4;
      Group=14;
      Series="Metalloid";
      Block="p";
      mass=AMU2KILOGRAM*72.59;
      MolarVolume=0.000013645;
      volume=19.2948;
      Miedema_Vm=4.6;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=5.323;
      crystal="dia";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.00058782;
      lattice_constant[1]=565.75;lattice_constant[2]=565.75;lattice_constant[3]=565.75;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.122;
      radius_PT=125;
      radius_covalent=1.20;
      radius_covalent_PT=120;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.2823;
      radii_Slatter=1.25;
      radii_Pyykko=1.21;
      electrical_conductivity=2000;
      electronegativity_vec=2.01;
      hardness_Ghosh=4.6166;
      electronegativity_Pearson=4.6;
      electronegativity_Ghosh=5.159;
      electron_affinity_PT=119;
      Miedema_phi_star=4.55;
      Miedema_nws=1.37;
      Miedema_gamma_s=1030;
      Pettifor_scale=1.92;
      boiling_point=2820;
      melting_point=938.3;
      vaporization_heat_PT=334;
      specific_heat_PT=321.4;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6E-6;
      thermal_conductivity=60;
      Brinelll_hardness=7273.402498871;
      Mohs_hardness=6;
      Vickers_hardness=8012.03305;
      Hardness_Pearson=3.40;
      Hardness_Putz=1.94;
      Hardness_RB=3.52;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=10.5;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-1.5E-9;
      Volume_Magnetic_Susceptibility=-7.98E-6;
      Molar_Magnetic_Susceptibility=-1.09E-10;
      Curie_point=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=5300;
      HHIR=1900;
      /*xray_scatt=NNN;*/
      //Ge
    }
    // [AFLOW]STOP=Germanium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Arsenic
    // Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic
    if(ZZ==33) { // Arsenic
      Z=ZZ;
      symbol="As";
      name="Arsenic";
      Period=4;
      Group=15;
      Series="Metalloid";
      Block="p";
      mass=AMU2KILOGRAM*74.9216;
      MolarVolume=0.000013082;
      volume=19.0677;
      Miedema_Vm=5.2;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      Density_PT=5.727;
      crystal="rhl";
      CrystalStr_PT="Simple_Trigonal";
      space_group="R_3m";
      space_group_number=166;
      Pearson_coefficient=0.0;
      lattice_constant[1]=375.98;lattice_constant[2]=375.98;lattice_constant[3]=1054.75;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.125;
      radius_PT=114;
      radius_covalent=1.19;
      radius_covalent_PT=119;
      radius_VanDerWaals_PT=185;
      radii_Ghosh08=1.1450;
      radii_Slatter=1.15;
      radii_Pyykko=1.21;
      electrical_conductivity=3.3E6;
      electronegativity_vec=2.18;
      hardness_Ghosh=5.0662;
      electronegativity_Pearson=5.3;
      electronegativity_Ghosh=5.628;
      electron_affinity_PT=78;
      Miedema_phi_star=4.80;
      Miedema_nws=1.44;
      Miedema_gamma_s=1000;
      Pettifor_scale=2.16;
      boiling_point=614;
      melting_point=817;
      vaporization_heat_PT=32.4;
      specific_heat_PT=328;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=50;
      Brinelll_hardness=1440;
      Mohs_hardness=3.5;
      Vickers_hardness=1510;
      Hardness_Pearson=4.50;
      Hardness_Putz=2.35;
      Hardness_RB=5.04;
      shear_modulus=NNN;
      Young_modulus=8;
      bulk_modulus=22;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=5.1;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-3.9E-9;
      Volume_Magnetic_Susceptibility=-0.0000223;
      Molar_Magnetic_Susceptibility=-2.92E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=1.001552;
      HHIP=3300;
      HHIR=4000;
      /*xray_scatt=NNN;*/
      //As
    }
    // [AFLOW]STOP=Arsenic
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Selenium
    // Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium
    if(ZZ==34) { // Selenium
      Z=ZZ;
      symbol="Se";
      name="Selenium";
      Period=4;
      Group=16;
      Series="Chalcogen";
      Block="p";
      mass=AMU2KILOGRAM*78.96;
      MolarVolume=0.000016387;
      volume=20.3733;
      Miedema_Vm=5.172;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      Density_PT=4.819;
      crystal="hex";
      CrystalStr_PT="Simple_Monoclinic";
      space_group="P12_1/c1";
      space_group_number=14;
      Pearson_coefficient=0.00046279;
      lattice_constant[1]=905.4;lattice_constant[2]=908.3;lattice_constant[3]=1160.1;
      lattice_angle[1]=PI/2;lattice_angle[2]=1.58493;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.116;
      radius_PT=103;
      radius_covalent=1.20;
      radius_covalent_PT=120;
      radius_VanDerWaals_PT=190;
      radii_Ghosh08=1.0424;
      radii_Slatter=1.15;
      radii_Pyykko=1.16;
      electrical_conductivity=NNN;
      electronegativity_vec=2.55;
      hardness_Ghosh=5.4795;
      electronegativity_Pearson=5.89;
      electronegativity_Ghosh=6.096;
      electron_affinity_PT=195;
      Miedema_phi_star=5.17;
      Miedema_nws=1.40;
      Miedema_gamma_s=NNN;
      Pettifor_scale=2.40;
      boiling_point=685;
      melting_point=221;
      vaporization_heat_PT=26;
      specific_heat_PT=321.2;
      critical_Pressure=268.4;
      critical_Temperature_PT=1766;
      thermal_expansion=NNN;
      thermal_conductivity=0.52;
      Brinelll_hardness=736;
      Mohs_hardness=2;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.87;
      Hardness_Putz=2.87;
      Hardness_RB=3.95;
      shear_modulus=3.7;
      Young_modulus=10;
      bulk_modulus=8.3;
      Poisson_ratio_PT=0.33;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-4E-9;
      Volume_Magnetic_Susceptibility=-0.0000193;
      Molar_Magnetic_Susceptibility=-3.16E-10;
      Curie_point=NNN;
      color_PT="GRAY";
      refractive_index=1.000895;
      HHIP=2200;
      HHIR=1900;
      /*xray_scatt=NNN;*/
      //Se Table 27 of JUNKAI
    }
    // [AFLOW]STOP=Selenium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Bromine
    // Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine
    if(ZZ==35) { // Bromine
      Z=ZZ;
      symbol="Br";
      name="Bromine";
      Period=4;
      Group=17;
      Series="Halogen";
      Block="p";
      mass=AMU2KILOGRAM*79.904;
      MolarVolume=0.00002561;
      volume=26.3292;
      Miedema_Vm=7.31;
      valence_std=7;
      valence_iupac=7;
      valence_PT=5;
      Density_PT=3.12;
      crystal="orc";
      CrystalStr_PT="Base_Orthorhombic";
      space_group="Cmca";
      space_group_number=64;
      Pearson_coefficient=0.000156277;
      lattice_constant[1]=672.65;lattice_constant[2]=464.51;lattice_constant[3]=870.23;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Liquid";
      radius=0.119;
      radius_PT=94;
      radius_covalent=1.20;
      radius_covalent_PT=120;
      radius_VanDerWaals_PT=185;
      radii_Ghosh08=0.9532;
      radii_Slatter=1.15;
      radii_Pyykko=1.14;
      electrical_conductivity=1E-10;
      electronegativity_vec=2.96;
      hardness_Ghosh=5.9111;
      electronegativity_Pearson=7.59;
      electronegativity_Ghosh=6.565;
      electron_affinity_PT=324.6;
      Miedema_phi_star=5.20;
      Miedema_nws=1.35;
      Miedema_gamma_s=943;
      Pettifor_scale=2.64;
      boiling_point=59;
      melting_point=-7.3;
      vaporization_heat_PT=14.8;
      specific_heat_PT=947.3;
      critical_Pressure=102;
      critical_Temperature_PT=588;
      thermal_expansion=NNN;
      thermal_conductivity=0.12;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=4.22;
      Hardness_Putz=3.39;
      Hardness_RB=4.4;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=1.9;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=3.4;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-4.9E-9;
      Volume_Magnetic_Susceptibility=-0.0000153;
      Molar_Magnetic_Susceptibility=-7.83E-10;
      Curie_point=NNN;
      color_PT="RED";
      refractive_index=1.001132;
      HHIP=3300;
      HHIR=6900;
      /* xray_scatt=NNN;*/
      //Br interpolation phi_star, nws, Vm, gamma, BVm JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Bromine
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Krypton
    // Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton
    if(ZZ==36) { // Krypton
      Z=ZZ;
      symbol="Kr";
      name="Krypton";
      Period=4;
      Group=18;
      Series="NobleGas";
      Block="p";
      mass=AMU2KILOGRAM*83.8;
      MolarVolume=0.02235;
      volume=-1.0000;
      Miedema_Vm=NNN;
      valence_std=0;
      valence_iupac=2;
      valence_PT=2;
      Density_PT=37.5E-4;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.000248482;
      lattice_constant[1]=570.6;lattice_constant[2]=570.6;lattice_constant[3]=570.6;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Gas";
      radius=0.197;
      radius_PT=87;
      radius_covalent=1.16;
      radius_covalent_PT=116;
      radius_VanDerWaals_PT=202;
      radii_Ghosh08=0.8782;
      radii_Slatter=NNN;
      radii_Pyykko=1.17;
      electrical_conductivity=NNN;
      electronegativity_vec=3;
      hardness_Ghosh=6.3418;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=7.033;
      electron_affinity_PT=0;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=-153.22;
      melting_point=-157.36;
      vaporization_heat_PT=9.02;
      specific_heat_PT=248.05;
      critical_Pressure=54.28;
      critical_Temperature_PT=209.41;
      thermal_expansion=NNN;
      thermal_conductivity=0.00943;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=3.98;
      Hardness_RB=9.45;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-4.4E-9;
      Volume_Magnetic_Susceptibility=-1.65E-8;
      Molar_Magnetic_Susceptibility=-3.69E-10;
      Curie_point=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000427;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Kr
    }
    // [AFLOW]STOP=Krypton
    // ********************************************************************************************************************************************************

    // ROW5
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Rubidium
    // Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium
    if(ZZ==37) { // Rubidium
      Z=ZZ;
      symbol="Rb";
      name="Rubidium";
      Period=5;
      Group=1;
      Series="AlkaliMetal";
      Block="s";
      mass=AMU2KILOGRAM*85.4678;
      MolarVolume=0.000055788;
      volume=91.2738;
      Miedema_Vm=14.6;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      Density_PT=1.532;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=0.000109697;
      lattice_constant[1]=558.5;lattice_constant[2]=558.5;lattice_constant[3]=558.5;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.251;
      radius_PT=265;
      radius_covalent=2.20;
      radius_covalent_PT=220;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=3.8487;
      radii_Slatter=2.35;
      radii_Pyykko=2.10;
      electrical_conductivity=8.3E6;
      electronegativity_vec=0.82;
      hardness_Ghosh=2.1204;
      electronegativity_Pearson=2.34;
      electronegativity_Ghosh=2.849;
      electron_affinity_PT=46.9;
      Miedema_phi_star=2.10;
      Miedema_nws=0.60;
      Miedema_gamma_s=120;
      Pettifor_scale=0.30;
      boiling_point=688;
      melting_point=39.31;
      vaporization_heat_PT=71;
      specific_heat_PT=364;
      critical_Pressure=157.9;
      critical_Temperature_PT=2093;
      thermal_expansion=NNN;
      thermal_conductivity=58;
      Brinelll_hardness=0.216;
      Mohs_hardness=0.3;
      Vickers_hardness=NNN;
      Hardness_Pearson=1.85;
      Hardness_Putz=0.08;
      Hardness_RB=2.21;
      shear_modulus=NNN;
      Young_modulus=2.4;
      bulk_modulus=2.5;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=1.8;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=2.6E-9;
      Volume_Magnetic_Susceptibility=3.98E-6;
      Molar_Magnetic_Susceptibility=2.22E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=6000;
      HHIR=6000;
      /*xray_scatt=NNN;*/
      //Rb
    }
    // [AFLOW]STOP=Rubidium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Strontium
    // Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium
    if(ZZ==38) { // Strontium
      Z=ZZ;
      symbol="Sr";
      name="Strontium";
      Period=5;
      Group=2;
      Series="AlkalineEarthMetal";
      Block="s";
      mass=AMU2KILOGRAM*87.62;
      MolarVolume=0.000033316;
      volume=55.4105;
      Miedema_Vm=10.2;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      Density_PT=2.63;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=6.09969E-05;
      lattice_constant[1]=608.49;lattice_constant[2]=608.49;lattice_constant[3]=608.49;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.215;
      radius_PT=219;
      radius_covalent=1.95;
      radius_covalent_PT=195;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.9709;
      radii_Slatter=2.00;
      radii_Pyykko=1.85;
      electrical_conductivity=7.7E6;
      electronegativity_vec=0.95;
      hardness_Ghosh=2.5374;
      electronegativity_Pearson=2.0;
      electronegativity_Ghosh=3.225;
      electron_affinity_PT=5.03;
      Miedema_phi_star=2.40;
      Miedema_nws=0.84;
      Miedema_gamma_s=430;
      Pettifor_scale=0.55;
      boiling_point=1382;
      melting_point=777;
      vaporization_heat_PT=137;
      specific_heat_PT=300;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000225;
      thermal_conductivity=35;
      Brinelll_hardness=NNN;
      Mohs_hardness=1.5;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.70;
      Hardness_Putz=0.11;
      Hardness_RB=3.08;
      shear_modulus=6.1;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=0.28;
      Miedema_BVm=3.9;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.32E-9;
      Volume_Magnetic_Susceptibility=3.47E-6;
      Molar_Magnetic_Susceptibility=1.16E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=4200;
      HHIR=3000;
      /*xray_scatt=NNN;*/
      //Sr
    }
    // [AFLOW]STOP=Strontium
    // ********************************************************************************************************************************************************

    // d-electron systems: transition metals
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Yttrium
    // Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium
    if(ZZ==39) { // Yttrium
      Z=ZZ;
      symbol="Y";
      name="Yttrium";
      Period=5;
      Group=3;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*88.9059;
      MolarVolume=0.000019881;
      volume=32.4546;
      Miedema_Vm=7.3;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=4.472;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=364.74;lattice_constant[2]=364.74;lattice_constant[3]=573.06;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.181;
      radius_PT=212;
      radius_covalent=1.90;
      radius_covalent_PT=190;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.8224;
      radii_Slatter=1.80;
      radii_Pyykko=1.63;
      electrical_conductivity=1.8E6;
      electronegativity_vec=1.22;
      hardness_Ghosh=2.6335;
      electronegativity_Pearson=3.19;
      electronegativity_Ghosh=3.311;
      electron_affinity_PT=29.6;
      Miedema_phi_star=3.20;
      Miedema_nws=1.21;
      Miedema_gamma_s=1100;
      Pettifor_scale=0.70;
      boiling_point=3345;
      melting_point=1526;
      vaporization_heat_PT=380;
      specific_heat_PT=298;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000106;
      thermal_conductivity=17;
      Brinelll_hardness=588;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.19;
      Hardness_Putz=0.14;
      Hardness_RB=3.67;
      shear_modulus=26;
      Young_modulus=64;
      bulk_modulus=41;
      Poisson_ratio_PT=0.24;
      Miedema_BVm=7.2;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=6.66E-8;
      Volume_Magnetic_Susceptibility=0.0002978;
      Molar_Magnetic_Susceptibility=5.921E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9800;
      HHIR=2600;
      /*xray_scatt=NNN;*/
      //Y
    }
    // [AFLOW]STOP=Yttrium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Zirconium
    // Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium
    if(ZZ==40) { // Zirconium
      Z=ZZ;
      symbol="Zr";
      name="Zirconium";
      Period=5;
      Group=4;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*91.22;
      MolarVolume=0.000014011;
      volume=23.2561;
      Miedema_Vm=5.8;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=6.511;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.000342629;
      lattice_constant[1]=323.2;lattice_constant[2]=323.2;lattice_constant[3]=514.7;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.158;
      radius_PT=206;
      radius_covalent=1.75;
      radius_covalent_PT=175;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.688;
      radii_Slatter=1.55;
      radii_Pyykko=1.54;
      electrical_conductivity=2.4E6;
      electronegativity_vec=1.33;
      hardness_Ghosh=2.7298;
      electronegativity_Pearson=3.64;
      electronegativity_Ghosh=3.398;
      electron_affinity_PT=41.1;
      Miedema_phi_star=3.40;
      Miedema_nws=1.39;
      Miedema_gamma_s=1950;
      Pettifor_scale=0.76;
      boiling_point=4409;
      melting_point=1855;
      vaporization_heat_PT=580;
      specific_heat_PT=278;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=5.7E-6;
      thermal_conductivity=23;
      Brinelll_hardness=650;
      Mohs_hardness=5;
      Vickers_hardness=904;
      Hardness_Pearson=3.21;
      Hardness_Putz=0.17;
      Hardness_RB=2.09;
      shear_modulus=33;
      Young_modulus=67;
      bulk_modulus=NNN;
      Poisson_ratio_PT=0.34;
      Miedema_BVm=12.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.68E-8;
      Volume_Magnetic_Susceptibility=0.000109;
      Molar_Magnetic_Susceptibility=1.53E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3400;
      HHIR=2600;
      /*xray_scatt=NNN;*/
      //Zr
    }
    // [AFLOW]STOP=Zirconium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Niobium
    // Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium
    if(ZZ==41) { // Niobium
      Z=ZZ;
      symbol="Nb";
      name="Niobium";
      Period=5;
      Group=5;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*92.9064;
      MolarVolume=0.000010841;
      volume=18.3132;
      Miedema_Vm=4.9;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      Density_PT=8.57;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=0.0;
      lattice_constant[1]=330.04;lattice_constant[2]=330.04;lattice_constant[3]=330.04;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.143;
      radius_PT=198;
      radius_covalent=1.64;
      radius_covalent_PT=164;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.5658;
      radii_Slatter=1.45;
      radii_Pyykko=1.47;
      electrical_conductivity=6.7E6;
      electronegativity_vec=1.60;
      hardness_Ghosh=2.8260;
      electronegativity_Pearson=4.0;
      electronegativity_Ghosh=3.485;
      electron_affinity_PT=86.1;
      Miedema_phi_star=4.00;
      Miedema_nws=1.62;
      Miedema_gamma_s=2700;
      Pettifor_scale=0.82;
      boiling_point=4744;
      melting_point=2477;
      vaporization_heat_PT=690;
      specific_heat_PT=265;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=7.3E-6;
      thermal_conductivity=54;
      Brinelll_hardness=736;
      Mohs_hardness=6;
      Vickers_hardness=1320;
      Hardness_Pearson=3.00;
      Hardness_Putz=0.21;
      Hardness_RB=3.67;
      shear_modulus=38;
      Young_modulus=105;
      bulk_modulus=170;
      Poisson_ratio_PT=0.4;
      Miedema_BVm=18.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=2.76E-8;
      Volume_Magnetic_Susceptibility=0.000237;
      Molar_Magnetic_Susceptibility=2.56E-9;
      Curie_point=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=8500;
      HHIR=8800;
      /*xray_scatt=NNN;*/
      //Nb
    }
    // [AFLOW]STOP=Niobium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Molybdenum
    // Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum
    if(ZZ==42) { // Molybdenum
      Z=ZZ;
      symbol="Mo";
      name="Molybdenum";
      Period=5;
      Group=6;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*95.94;
      MolarVolume=9.334E-6;
      volume=15.6175;
      Miedema_Vm=4.4;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      Density_PT=10.28;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=0.000598128;
      lattice_constant[1]=314.7;lattice_constant[2]=314.7;lattice_constant[3]=314.7;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.136;
      radius_PT=190;
      radius_covalent=1.54;
      radius_covalent_PT=154;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.4543;
      radii_Slatter=1.45;
      radii_Pyykko=1.38;
      electrical_conductivity=2E7;
      electronegativity_vec=2.16;
      hardness_Ghosh=2.9221;
      electronegativity_Pearson=3.9;
      electronegativity_Ghosh=3.572;
      electron_affinity_PT=71.9;
      Miedema_phi_star=4.65;
      Miedema_nws=1.77;
      Miedema_gamma_s=2950;
      Pettifor_scale=0.88;
      boiling_point=4639;
      melting_point=2623;
      vaporization_heat_PT=600;
      specific_heat_PT=251;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=4.8E-6;
      thermal_conductivity=139;
      Brinelll_hardness=1500;
      Mohs_hardness=5.5;
      Vickers_hardness=1530;
      Hardness_Pearson=3.10;
      Hardness_Putz=0.25;
      Hardness_RB=NNN;
      shear_modulus=20;
      Young_modulus=329;
      bulk_modulus=230;
      Poisson_ratio_PT=0.31;
      Miedema_BVm=26.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.17E-8;
      Volume_Magnetic_Susceptibility=0.0001203;
      Molar_Magnetic_Susceptibility=1.122E-9;
      Curie_point=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=2400;
      HHIR=5300;
      /*xray_scatt=NNN;*/
      //Mo
    }
    // [AFLOW]STOP=Molybdenum
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Technetium
    // Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium
    if(ZZ==43) { // Technetium
      Z=ZZ;
      symbol="Tc";
      name="Technetium";
      Period=5;
      Group=7;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*98.9062;
      MolarVolume=8.434782608696E-6;
      volume=14.4670;
      Miedema_Vm=4.2;
      valence_std=7;
      valence_iupac=7;
      valence_PT=6;
      Density_PT=11.5;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=273.5;lattice_constant[2]=273.5;lattice_constant[3]=438.8;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=NNN;
      radius_PT=183;
      radius_covalent=1.47;
      radius_covalent_PT=147;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.352;
      radii_Slatter=1.35;
      radii_Pyykko=1.28;
      electrical_conductivity=5E6;
      electronegativity_vec=1.90;
      hardness_Ghosh=3.0184;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=3.659;
      electron_affinity_PT=53;
      Miedema_phi_star=5.30;
      Miedema_nws=1.81;
      Miedema_gamma_s=3050;
      Pettifor_scale=0.94;
      boiling_point=4265;
      melting_point=2157;
      vaporization_heat_PT=550;
      specific_heat_PT=63;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=51;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=0.29;
      Hardness_RB=2.05;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=26.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=3.42E-8;
      Volume_Magnetic_Susceptibility=0.0003933;
      Molar_Magnetic_Susceptibility=3.352E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Tc JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Technetium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Ruthenium
    // Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium
    if(ZZ==44) { // Ruthenium
      Z=ZZ;
      symbol="Ru";
      name="Ruthenium";
      Period=5;
      Group=8;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*101.07;
      MolarVolume=8.1706E-6;
      volume=13.8390;
      Miedema_Vm=4.1;
      valence_std=8;
      valence_iupac=8;
      valence_PT=6;
      Density_PT=12.37;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.000406665;
      lattice_constant[1]=270.59;lattice_constant[2]=270.59;lattice_constant[3]=428.15;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.134;
      radius_PT=178;
      radius_covalent=1.46;
      radius_covalent_PT=146;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.2579;
      radii_Slatter=1.30;
      radii_Pyykko=1.25;
      electrical_conductivity=1.4E7;
      electronegativity_vec=2.20;
      hardness_Ghosh=3.1146;
      electronegativity_Pearson=4.5;
      electronegativity_Ghosh=3.745;
      electron_affinity_PT=101.3;
      Miedema_phi_star=5.40;
      Miedema_nws=1.83;
      Miedema_gamma_s=3050;
      Pettifor_scale=1.00;
      boiling_point=4150;
      melting_point=2334;
      vaporization_heat_PT=580;
      specific_heat_PT=238;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.4E-6;
      thermal_conductivity=120;
      Brinelll_hardness=2160;
      Mohs_hardness=6.5;
      Vickers_hardness=2298.138766667;
      Hardness_Pearson=3.00;
      Hardness_Putz=0.35;
      Hardness_RB=NNN;
      shear_modulus=173;
      Young_modulus=447;
      bulk_modulus=220;
      Poisson_ratio_PT=0.3;
      Miedema_BVm=26.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=5.42E-9;
      Volume_Magnetic_Susceptibility=0.000067;
      Molar_Magnetic_Susceptibility=5.48E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3200;
      HHIR=8000;
      /*xray_scatt=NNN;*/
      //Ru JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Ruthenium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Rhodium
    // Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium
    if(ZZ==45) { // Rhodium
      Z=ZZ;
      symbol="Rh";
      name="Rhodium";
      Period=5;
      Group=9;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*102.9055;
      MolarVolume=8.2655E-6;
      volume=14.1731;
      Miedema_Vm=4.1;
      valence_std=9;
      valence_iupac=6;
      valence_PT=6;
      Density_PT=12.45;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=1.90706E-32;
      lattice_constant[1]=380.34;lattice_constant[2]=380.34;lattice_constant[3]=380.34;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.134;
      radius_PT=173;
      radius_covalent=1.42;
      radius_covalent_PT=142;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.1711;
      radii_Slatter=1.35;
      radii_Pyykko=1.25;
      electrical_conductivity=2.3E7;
      electronegativity_vec=2.28;
      hardness_Ghosh=3.2108;
      electronegativity_Pearson=4.3;
      electronegativity_Ghosh=3.832;
      electron_affinity_PT=109.7;
      Miedema_phi_star=5.40;
      Miedema_nws=1.76;
      Miedema_gamma_s=2750;
      Pettifor_scale=1.06;
      boiling_point=3695;
      melting_point=1964;
      vaporization_heat_PT=495;
      specific_heat_PT=240;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=8E-6;
      thermal_conductivity=150;
      Brinelll_hardness=1100;
      Mohs_hardness=6;
      Vickers_hardness=1246;
      Hardness_Pearson=3.16;
      Hardness_Putz=0.41;
      Hardness_RB=NNN;
      shear_modulus=150;
      Young_modulus=275;
      bulk_modulus=380;
      Poisson_ratio_PT=0.26;
      Miedema_BVm=23.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.36E-8;
      Volume_Magnetic_Susceptibility=0.0001693;
      Molar_Magnetic_Susceptibility=1.4E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3200;
      HHIR=8000;
      /*xray_scatt=NNN;*/
      //Rh
    }
    // [AFLOW]STOP=Rhodium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Palladium
    // Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium
    if(ZZ==46) { // Palladium
      Z=ZZ;
      symbol="Pd";
      name="Palladium";
      Period=5;
      Group=10;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*106.4;
      MolarVolume=8.8514E-6;
      volume=15.4596;
      Miedema_Vm=4.3;
      valence_std=10;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=12.023;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.000309478;
      lattice_constant[1]=389.07;lattice_constant[2]=389.07;lattice_constant[3]=389.07;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.137;
      radius_PT=169;
      radius_covalent=1.39;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=163;
      radii_Ghosh08=2.0907;
      radii_Slatter=1.40;
      radii_Pyykko=1.20;
      electrical_conductivity=1E7;
      electronegativity_vec=2.20;
      hardness_Ghosh=3.3069;
      electronegativity_Pearson=4.45;
      electronegativity_Ghosh=3.919;
      electron_affinity_PT=53.7;
      Miedema_phi_star=5.45;
      Miedema_nws=1.67;
      Miedema_gamma_s=2100;
      Pettifor_scale=1.12;
      boiling_point=2963;
      melting_point=1554.9;
      vaporization_heat_PT=380;
      specific_heat_PT=240;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000118;
      thermal_conductivity=71;
      Brinelll_hardness=37.2;
      Mohs_hardness=4.75;
      Vickers_hardness=461;
      Hardness_Pearson=3.89;
      Hardness_Putz=0.47;
      Hardness_RB=6.32;
      shear_modulus=44;
      Young_modulus=121;
      bulk_modulus=180;
      Poisson_ratio_PT=0.39;
      Miedema_BVm=16.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=6.57E-8;
      Volume_Magnetic_Susceptibility=0.0007899;
      Molar_Magnetic_Susceptibility=6.992E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3200;
      HHIR=8000;
      /*xray_scatt=NNN;*/
      //Pd
    }
    // [AFLOW]STOP=Palladium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Silver
    // Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver
    if(ZZ==47) { // Silver
      Z=ZZ;
      symbol="Ag";
      name="Silver";
      Period=5;
      Group=11;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*107.8682;
      MolarVolume=0.000010283;
      volume=18.0678;
      Miedema_Vm=4.7;
      valence_std=11;
      valence_iupac=4;
      valence_PT=1;
      Density_PT=10.49;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=8.57985E-05;
      lattice_constant[1]=408.53;lattice_constant[2]=408.53;lattice_constant[3]=408.53;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.144;
      radius_PT=165;
      radius_covalent=1.45;
      radius_covalent_PT=145;
      radius_VanDerWaals_PT=172;
      radii_Ghosh08=2.016;
      radii_Slatter=1.60;
      radii_Pyykko=1.28;
      electrical_conductivity=6.2E7;
      electronegativity_vec=1.93;
      hardness_Ghosh=3.4032;
      electronegativity_Pearson=4.44;
      electronegativity_Ghosh=4.006;
      electron_affinity_PT=125.6;
      Miedema_phi_star=4.45;
      Miedema_nws=1.39;
      Miedema_gamma_s=1250;
      Pettifor_scale=1.18;
      boiling_point=2162;
      melting_point=961.78;
      vaporization_heat_PT=255;
      specific_heat_PT=235;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000189;
      thermal_conductivity=430;
      Brinelll_hardness=24.5;
      Mohs_hardness=2.5;
      Vickers_hardness=251;
      Hardness_Pearson=3.14;
      Hardness_Putz=0.55;
      Hardness_RB=3.5;
      shear_modulus=30;
      Young_modulus=85;
      bulk_modulus=100;
      Poisson_ratio_PT=0.37;
      Miedema_BVm=10.0;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-2.27E-9;
      Volume_Magnetic_Susceptibility=-0.0000238;
      Molar_Magnetic_Susceptibility=-2.45E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1200;
      HHIR=1400;
      xray_scatt=47.18;
      //Ag JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Silver
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Cadmium
    // Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium
    if(ZZ==48) { // Cadmium
      Z=ZZ;
      symbol="Cd";
      name="Cadmium";
      Period=5;
      Group=12;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*112.41;
      MolarVolume=0.000012996;
      volume=22.0408;
      Miedema_Vm=5.5;
      valence_std=12;
      valence_iupac=2;
      valence_PT=2;
      Density_PT=8.65;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.000271603;
      lattice_constant[1]=297.94;lattice_constant[2]=297.94;lattice_constant[3]=561.86;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.150;
      radius_PT=161;
      radius_covalent=1.44;
      radius_covalent_PT=144;
      radius_VanDerWaals_PT=158;
      radii_Ghosh08=1.9465;
      radii_Slatter=1.55;
      radii_Pyykko=1.36;
      electrical_conductivity=1.4E7;
      electronegativity_vec=1.69;
      hardness_Ghosh=3.4994;
      electronegativity_Pearson=4.33;
      electronegativity_Ghosh=4.093;
      electron_affinity_PT=0;
      Miedema_phi_star=4.05;
      Miedema_nws=1.24;
      Miedema_gamma_s=780;
      Pettifor_scale=1.36;
      boiling_point=767;
      melting_point=321.07;
      vaporization_heat_PT=100;
      specific_heat_PT=230;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000308;
      thermal_conductivity=96;
      Brinelll_hardness=203;
      Mohs_hardness=2;
      Vickers_hardness=NNN;
      Hardness_Pearson=4.66;
      Hardness_Putz=0.63;
      Hardness_RB=5.35;
      shear_modulus=19;
      Young_modulus=50;
      bulk_modulus=42;
      Poisson_ratio_PT=0.3;
      Miedema_BVm=6.10;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-2.3E-9;
      Volume_Magnetic_Susceptibility=-0.0000199;
      Molar_Magnetic_Susceptibility=-2.59E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1700;
      HHIR=1300;
      /*xray_scatt=NNN;*/
      //Cd
    }
    // [AFLOW]STOP=Cadmium
    // ********************************************************************************************************************************************************

    // p-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Indium
    // Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium
    if(ZZ==49) { // Indium
      Z=ZZ;
      symbol="In";
      name="Indium";
      Period=5;
      Group=13;
      Series="PoorMetal";
      Block="p";
      mass=AMU2KILOGRAM*114.82;
      MolarVolume=0.000015707;
      volume=27.5233;
      Miedema_Vm=6.3;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=7.31;
      crystal="fct";
      CrystalStr_PT="Centered_Tetragonal";
      space_group="I4/mmm";
      space_group_number=139;
      Pearson_coefficient=1.24494E-05;
      lattice_constant[1]=325.23;lattice_constant[2]=325.23;lattice_constant[3]=494.61;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.157;
      radius_PT=156;
      radius_covalent=1.42;
      radius_covalent_PT=142;
      radius_VanDerWaals_PT=193;
      radii_Ghosh08=1.6934;
      radii_Slatter=1.55;
      radii_Pyykko=1.42;
      electrical_conductivity=1.2E7;
      electronegativity_vec=1.78;
      hardness_Ghosh=3.9164;
      electronegativity_Pearson=3.1;
      electronegativity_Ghosh=4.469;
      electron_affinity_PT=28.9;
      Miedema_phi_star=3.90;
      Miedema_nws=1.17;
      Miedema_gamma_s=690;
      Pettifor_scale=1.60;
      boiling_point=2072;
      melting_point=156.6;
      vaporization_heat_PT=230;
      specific_heat_PT=233;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000321;
      thermal_conductivity=82;
      Brinelll_hardness=8.83;
      Mohs_hardness=1.2;
      Vickers_hardness=NNN;
      Hardness_Pearson=2.80;
      Hardness_Putz=0.73;
      Hardness_RB=2.77;
      shear_modulus=NNN;
      Young_modulus=11;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=6.4;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-1.4E-9;
      Volume_Magnetic_Susceptibility=-0.0000102;
      Molar_Magnetic_Susceptibility=-1.61E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3300;
      HHIR=2000;
      /*xray_scatt=NNN;*/
      //In
    }
    // [AFLOW]STOP=Indium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Tin
    // Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin
    if(ZZ==50) { // Tin
      Z=ZZ;
      symbol="Sn";
      name="Tin";
      Period=5;
      Group=14;
      Series="PoorMetal";
      Block="p";
      mass=AMU2KILOGRAM*118.69;
      MolarVolume=0.000016239;
      volume=27.5555;
      Miedema_Vm=6.4;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=7.31;
      crystal="bct";
      CrystalStr_PT="Centered_Tetragonal";
      space_group="I4_1/amd";
      space_group_number=141;
      Pearson_coefficient=0.000334085;
      lattice_constant[1]=583.18;lattice_constant[2]=583.18;lattice_constant[3]=318.19;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.158;
      radius_PT=145;
      radius_covalent=1.39;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=217;
      radii_Ghosh08=1.4986;
      radii_Slatter=1.45;
      radii_Pyykko=1.40;
      electrical_conductivity=9.1E6;
      electronegativity_vec=1.96;
      hardness_Ghosh=4.3332;
      electronegativity_Pearson=4.3;
      electronegativity_Ghosh=4.845;
      electron_affinity_PT=107.3;
      Miedema_phi_star=4.15;
      Miedema_nws=1.24;
      Miedema_gamma_s=710;
      Pettifor_scale=1.84;
      boiling_point=2602;
      melting_point=231.93;
      vaporization_heat_PT=290;
      specific_heat_PT=217;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000022;
      thermal_conductivity=67;
      Brinelll_hardness=51;
      Mohs_hardness=1.5;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.05;
      Hardness_Putz=0.88;
      Hardness_RB=3.15;
      shear_modulus=18;
      Young_modulus=50;
      bulk_modulus=58;
      Poisson_ratio_PT=0.36;
      Miedema_BVm=8.8;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-3.1E-9;
      Volume_Magnetic_Susceptibility=-0.0000227;
      Molar_Magnetic_Susceptibility=-3.68E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=2600;
      HHIR=1600;
      /*xray_scatt=NNN;*/
      //Sn
    }
    // [AFLOW]STOP=Tin
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Antimony
    // Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony
    if(ZZ==51) { // Antimony
      Z=ZZ;
      symbol="Sb";
      name="Antimony";
      Period=5;
      Group=15;
      Series="Metalloid";
      Block="p";
      mass=AMU2KILOGRAM*121.75;
      MolarVolume=0.000018181;
      volume=27.1823;
      Miedema_Vm=6.6;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      Density_PT=6.697;
      crystal="rhl";
      CrystalStr_PT="Simple_Trigonal";
      space_group="R_3m";
      space_group_number=166;
      Pearson_coefficient=6.60751E-05;
      lattice_constant[1]=430.7;lattice_constant[2]=430.7;lattice_constant[3]=1127.3;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.161;
      radius_PT=133;
      radius_covalent=1.39;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.344;
      radii_Slatter=1.45;
      radii_Pyykko=1.40;
      electrical_conductivity=2.5E6;
      electronegativity_vec=2.05;
      hardness_Ghosh=4.7501;
      electronegativity_Pearson=4.85;
      electronegativity_Ghosh=5.221;
      electron_affinity_PT=103.2;
      Miedema_phi_star=4.40;
      Miedema_nws=1.26;
      Miedema_gamma_s=680;
      Pettifor_scale=2.08;
      boiling_point=1587;
      melting_point=630.63;
      vaporization_heat_PT=67;
      specific_heat_PT=207;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000011;
      thermal_conductivity=24;
      Brinelll_hardness=294;
      Mohs_hardness=3;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.80;
      Hardness_Putz=1.10;
      Hardness_RB=4.39;
      shear_modulus=20;
      Young_modulus=55;
      bulk_modulus=42;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=7.0;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-1.09E-8;
      Volume_Magnetic_Susceptibility=-0.000073;
      Molar_Magnetic_Susceptibility=-1.327E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=7900;
      HHIR=3400;
      /*xray_scatt=NNN;*/
      //Sb
    }
    // [AFLOW]STOP=Antimony
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Tellurium
    // Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium
    if(ZZ==52) { // Tellurium
      Z=ZZ;
      symbol="Te";
      name="Tellurium";
      Period=5;
      Group=16;
      Series="Chalcogen";
      Block="p";
      mass=AMU2KILOGRAM*127.6;
      MolarVolume=0.000020449;
      volume=28.1993;
      Miedema_Vm=6.439;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      Density_PT=6.24;
      crystal="hex";
      CrystalStr_PT="Simple_Trigonal";
      space_group="P3_121";
      space_group_number=152;
      Pearson_coefficient=0.000283934;
      lattice_constant[1]=445.72;lattice_constant[2]=445.72;lattice_constant[3]=592.9;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.143;
      radius_PT=123;
      radius_covalent=1.38;
      radius_covalent_PT=138;
      radius_VanDerWaals_PT=206;
      radii_Ghosh08=1.2183;
      radii_Slatter=1.40;
      radii_Pyykko=1.36;
      electrical_conductivity=10000;
      electronegativity_vec=2.10;
      hardness_Ghosh=5.1670;
      electronegativity_Pearson=5.49;
      electronegativity_Ghosh=5.597;
      electron_affinity_PT=190.2;
      Miedema_phi_star=4.72;
      Miedema_nws=1.31;
      Miedema_gamma_s=NNN;
      Pettifor_scale=2.32;
      boiling_point=988;
      melting_point=449.51;
      vaporization_heat_PT=48;
      specific_heat_PT=201;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=3;
      Brinelll_hardness=180;
      Mohs_hardness=2.25;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.52;
      Hardness_Putz=1.34;
      Hardness_RB=3.47;
      shear_modulus=16;
      Young_modulus=43;
      bulk_modulus=64;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-3.9E-9;
      Volume_Magnetic_Susceptibility=-0.0000243;
      Molar_Magnetic_Susceptibility=-4.98E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=1.000991;
      HHIP=2900;
      HHIR=4900;
      /*xray_scatt=NNN;*/
      //Te Table 27 of JUNKAI
    }
    // [AFLOW]STOP=Tellurium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Iodine
    // Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine
    if(ZZ==53) { // Iodine
      Z=ZZ;
      symbol="I";
      name="Iodine";
      Period=5;
      Group=17;
      Series="Halogen";
      Block="p";
      mass=AMU2KILOGRAM*126.9045;
      MolarVolume=0.000025689;
      volume=34.9784;
      Miedema_Vm=8.72;
      valence_std=7;
      valence_iupac=7;
      valence_PT=7;
      Density_PT=4.94;
      crystal="orc";
      CrystalStr_PT="Base_Orthorhombic";
      space_group="Cmca";
      space_group_number=64;
      Pearson_coefficient=0.0;
      lattice_constant[1]=718.02;lattice_constant[2]=471.02;lattice_constant[3]=981.03;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.136;
      radius_PT=115;
      radius_covalent=1.39;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=198;
      radii_Ghosh08=1.1141;
      radii_Slatter=1.40;
      radii_Pyykko=1.33;
      electrical_conductivity=1E-7;
      electronegativity_vec=2.66;
      hardness_Ghosh=5.5839;
      electronegativity_Pearson=6.76;
      electronegativity_Ghosh=5.973;
      electron_affinity_PT=295.2;
      Miedema_phi_star=5.33;
      Miedema_nws=0.17;
      Miedema_gamma_s=NNN;
      Pettifor_scale=2.56;
      boiling_point=184.3;
      melting_point=113.7;
      vaporization_heat_PT=20.9;
      specific_heat_PT=429;
      critical_Pressure=115.5;
      critical_Temperature_PT=819;
      thermal_expansion=NNN;
      thermal_conductivity=0.449;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.69;
      Hardness_Putz=1.62;
      Hardness_RB=3.81;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=7.7;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-4.5E-9;
      Volume_Magnetic_Susceptibility=-0.0000222;
      Molar_Magnetic_Susceptibility=-1.14E-9;
      Curie_point=NNN;
      color_PT="SLATEGRAY";
      refractive_index=NNN;
      HHIP=4900;
      HHIR=4800;
      /*xray_scatt=NNN;*/
      //I interpolation phi_star, nws, Vm,
    }
    // [AFLOW]STOP=Iodine
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Xenon
    // Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon
    if(ZZ==54) { // Xenon
      Z=ZZ;
      symbol="Xe";
      name="Xenon";
      Period=5;
      Group=18;
      Series="NobleGas";
      Block="p";
      mass=AMU2KILOGRAM*131.3;
      MolarVolume=0.0223;
      volume=-1.0000;
      Miedema_Vm=NNN;
      valence_std=0;
      valence_iupac=8;
      valence_PT=6;
      Density_PT=59E-4;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.000267781;
      lattice_constant[1]=620.23;lattice_constant[2]=620.23;lattice_constant[3]=620.23;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Gas";
      radius=0.218;
      radius_PT=108;
      radius_covalent=1.40;
      radius_covalent_PT=140;
      radius_VanDerWaals_PT=216;
      radii_Ghosh08=1.0263;
      radii_Slatter=NNN;
      radii_Pyykko=1.31;
      electrical_conductivity=NNN;
      electronegativity_vec=2.60;
      hardness_Ghosh=6.0009;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.349;
      electron_affinity_PT=0;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=-108;
      melting_point=-111.8;
      vaporization_heat_PT=12.64;
      specific_heat_PT=158.32;
      critical_Pressure=57.65;
      critical_Temperature_PT=289.77;
      thermal_expansion=NNN;
      thermal_conductivity=0.00565;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=1.92;
      Hardness_RB=8.23;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-4.3E-9;
      Volume_Magnetic_Susceptibility=-2.54E-8;
      Molar_Magnetic_Susceptibility=-5.65E-10;
      Curie_point=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000702;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Xe JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Xenon
    // ********************************************************************************************************************************************************

    // ROW6
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Cesium
    // Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium
    if(ZZ==55) { // Cesium
      Z=ZZ;
      symbol="Cs";
      name="Cesium";
      Period=6;
      Group=1;
      Series="AlkaliMetal";
      Block="s";
      mass=AMU2KILOGRAM*132.9054;
      MolarVolume=0.000070732;
      volume=117.281;
      Miedema_Vm=16.8;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      Density_PT=1.879;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=0.0;
      lattice_constant[1]=614.1;lattice_constant[2]=614.1;lattice_constant[3]=614.1;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.265;
      radius_PT=298;
      radius_covalent=2.44;
      radius_covalent_PT=244;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=4.2433;
      radii_Slatter=2.60;
      radii_Pyykko=2.32;
      electrical_conductivity=5E6;
      electronegativity_vec=0.79;
      hardness_Ghosh=0.6829;
      electronegativity_Pearson=2.18;
      electronegativity_Ghosh=4.196;
      electron_affinity_PT=45.5;
      Miedema_phi_star=1.95;
      Miedema_nws=0.55;
      Miedema_gamma_s=95;
      Pettifor_scale=0.25;
      boiling_point=671;
      melting_point=28.44;
      vaporization_heat_PT=64;
      specific_heat_PT=242;
      critical_Pressure=92.77;
      critical_Temperature_PT=1938;
      thermal_expansion=NNN;
      thermal_conductivity=36;
      Brinelll_hardness=0.14;
      Mohs_hardness=0.2;
      Vickers_hardness=NNN;
      Hardness_Pearson=1.71;
      Hardness_Putz=NNN;
      Hardness_RB=1.98;
      shear_modulus=NNN;
      Young_modulus=1.7;
      bulk_modulus=1.6;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=1.4;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=2.8E-9;
      Volume_Magnetic_Susceptibility=5.26E-6;
      Molar_Magnetic_Susceptibility=3.72E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=6000;
      HHIR=6000;
      /*xray_scatt=NNN;*/
      //Cs
    }
    // [AFLOW]STOP=Cesium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Barium
    // Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium
    if(ZZ==56) { // Barium
      Z=ZZ;
      symbol="Ba";
      name="Barium";
      Period=6;
      Group=2;
      Series="AlkalineEarthMetal";
      Block="s";
      mass=AMU2KILOGRAM*137.33;
      MolarVolume=0.000039125;
      volume=62.6649;
      Miedema_Vm=11.3;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      Density_PT=3.51;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=6.23705E-05;
      lattice_constant[1]=502.8;lattice_constant[2]=502.8;lattice_constant[3]=502.8;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.217;
      radius_PT=253;
      radius_covalent=2.15;
      radius_covalent_PT=215;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=3.2753;
      radii_Slatter=2.15;
      radii_Pyykko=1.96;
      electrical_conductivity=2.9E6;
      electronegativity_vec=0.89;
      hardness_Ghosh=0.9201;
      electronegativity_Pearson=2.4;
      electronegativity_Ghosh=4.318;
      electron_affinity_PT=13.95;
      Miedema_phi_star=2.32;
      Miedema_nws=0.81;
      Miedema_gamma_s=370;
      Pettifor_scale=0.50;
      boiling_point=1870;
      melting_point=727;
      vaporization_heat_PT=140;
      specific_heat_PT=205;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000206;
      thermal_conductivity=18;
      Brinelll_hardness=NNN;
      Mohs_hardness=1.25;
      Vickers_hardness=NNN;
      Hardness_Pearson=2.90;
      Hardness_Putz=NNN;
      Hardness_RB=2.16;
      shear_modulus=4.9;
      Young_modulus=13;
      bulk_modulus=9.4;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=3.9;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.13E-8;
      Volume_Magnetic_Susceptibility=0.00003966;
      Molar_Magnetic_Susceptibility=1.552E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3000;
      HHIR=2300;
      /*xray_scatt=NNN;*/
      //Ba
    }
    // [AFLOW]STOP=Barium
    // ********************************************************************************************************************************************************

    // d-electron systems: transition metals
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Lanthanium
    // Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium
    if(ZZ==57) { // Lanthanium
      Z=ZZ;
      symbol="La";
      name="Lanthanium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*138.9055;
      MolarVolume=0.000022601;
      volume=36.8495;
      Miedema_Vm=8.0;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=6.146;
      crystal="hex";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=4.65323E-08;
      lattice_constant[1]=377.2;lattice_constant[2]=377.2;lattice_constant[3]=1214.4;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.187;
      radius_PT=NNN;
      radius_covalent=2.07;
      radius_covalent_PT=207;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.6673;
      radii_Slatter=1.95;
      radii_Pyykko=1.80;
      electrical_conductivity=1.6E6;
      electronegativity_vec=1.10;
      hardness_Ghosh=1.1571;
      electronegativity_Pearson=3.1;
      electronegativity_Ghosh=4.439;
      electron_affinity_PT=48;
      Miedema_phi_star=3.05;
      Miedema_nws=1.09;
      Miedema_gamma_s=900;
      Pettifor_scale=0.7480;
      boiling_point=3464;
      melting_point=919;
      vaporization_heat_PT=400;
      specific_heat_PT=195;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000121;
      thermal_conductivity=13;
      Brinelll_hardness=363;
      Mohs_hardness=2.5;
      Vickers_hardness=491;
      Hardness_Pearson=2.60;
      Hardness_Putz=NNN;
      Hardness_RB=2.46;
      shear_modulus=14;
      Young_modulus=37;
      bulk_modulus=28;
      Poisson_ratio_PT=0.28;
      Miedema_BVm=5.5;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.1E-8;
      Volume_Magnetic_Susceptibility=0.00006761;
      Molar_Magnetic_Susceptibility=1.528E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //La
    }
    // [AFLOW]STOP=Lanthanium
    // ********************************************************************************************************************************************************

    // lantanidies
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Cerium
    // Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium
    if(ZZ==58) { // Cerium
      Z=ZZ;
      symbol="Ce";
      name="Cerium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*140.12;
      MolarVolume=0.000020947;
      volume=26.4729;
      Miedema_Vm=7.76;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=6.689;
      crystal="fcc";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=2.24956E-05;
      lattice_constant[1]=362;lattice_constant[2]=362;lattice_constant[3]=599;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.182;
      radius_PT=NNN;
      radius_covalent=2.04;
      radius_covalent_PT=204;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.2494;
      radii_Slatter=1.85;
      radii_Pyykko=1.63;
      electrical_conductivity=1.4E6;
      electronegativity_vec=1.12;
      hardness_Ghosh=1.3943;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.561;
      electron_affinity_PT=50;
      Miedema_phi_star=3.18;
      Miedema_nws=1.19;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7460;
      boiling_point=3360;
      melting_point=798;
      vaporization_heat_PT=350;
      specific_heat_PT=192;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.3E-6;
      thermal_conductivity=11;
      Brinelll_hardness=412;
      Mohs_hardness=2.5;
      Vickers_hardness=270;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=1.8;
      shear_modulus=14;
      Young_modulus=34;
      bulk_modulus=22;
      Poisson_ratio_PT=0.24;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=2.2E-7;
      Volume_Magnetic_Susceptibility=0.0014716;
      Molar_Magnetic_Susceptibility=3.0826E-8;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Ce Pettifor linear interpolation// Miedema from Alonso-March.
    }
    // [AFLOW]STOP=Cerium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Praseodymium
    // Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium
    if(ZZ==59) { // Praseodymium
      Z=ZZ;
      symbol="Pr";
      name="Praseodymium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*140.9077;
      MolarVolume=0.000021221;
      volume=36.4987;
      Miedema_Vm=7.56;
      valence_std=5;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=6.64;
      crystal="hex";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=367.25;lattice_constant[2]=367.25;lattice_constant[3]=1183.54;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.183;
      radius_PT=247;
      radius_covalent=2.03;
      radius_covalent_PT=203;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.9447;
      radii_Slatter=1.85;
      radii_Pyykko=1.76;
      electrical_conductivity=1.4E6;
      electronegativity_vec=1.13;
      hardness_Ghosh=1.6315;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.682;
      electron_affinity_PT=50;
      Miedema_phi_star=3.19;
      Miedema_nws=1.20;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7440;
      boiling_point=3290;
      melting_point=931;
      vaporization_heat_PT=330;
      specific_heat_PT=193;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.7E-6;
      thermal_conductivity=13;
      Brinelll_hardness=481;
      Mohs_hardness=1.41;
      Vickers_hardness=400;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=1.11;
      shear_modulus=15;
      Young_modulus=37;
      bulk_modulus=29;
      Poisson_ratio_PT=0.28;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=4.23E-7;
      Volume_Magnetic_Susceptibility=0.0028087;
      Molar_Magnetic_Susceptibility=5.9604E-8;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Pr Pettifor linear interpolation
    }
    // [AFLOW]STOP=Praseodymium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Neodymium
    // Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium
    if(ZZ==60) { // Neodymium
      Z=ZZ;
      symbol="Nd";
      name="Neodymium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*144.24;
      MolarVolume=0.000020577;
      volume=29.6719;
      Miedema_Vm=7.51;
      valence_std=6;
      valence_iupac=4;
      valence_PT=3;
      Density_PT=7.01;
      crystal="hex";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.000231599;
      lattice_constant[1]=365.8;lattice_constant[2]=365.8;lattice_constant[3]=1179.9;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.182;
      radius_PT=206;
      radius_covalent=2.01;
      radius_covalent_PT=201;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.7129;
      radii_Slatter=1.85;
      radii_Pyykko=1.74;
      electrical_conductivity=1.6E6;
      electronegativity_vec=1.14;
      hardness_Ghosh=1.8684;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.804;
      electron_affinity_PT=50;
      Miedema_phi_star=3.19;
      Miedema_nws=1.20;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7420;
      boiling_point=3100;
      melting_point=1021;
      vaporization_heat_PT=285;
      specific_heat_PT=190;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=9.6E-6;
      thermal_conductivity=17;
      Brinelll_hardness=265;
      Mohs_hardness=1.23;
      Vickers_hardness=343;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=0.7;
      shear_modulus=16;
      Young_modulus=41;
      bulk_modulus=32;
      Poisson_ratio_PT=0.28;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=4.8E-7;
      Volume_Magnetic_Susceptibility=0.0033648;
      Molar_Magnetic_Susceptibility=6.9235E-8;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Nd Pettifor linear interpolation JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Neodymium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Promethium
    // Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium
    if(ZZ==61) { // Promethium
      Z=ZZ;
      symbol="Pm";
      name="Promethium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*146.92;
      MolarVolume=0.00001996145374449;
      volume=34.6133;
      Miedema_Vm=7.43;
      valence_std=7;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=7.264;
      crystal="hex";
      CrystalStr_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      Pearson_coefficient=0.0;
      lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
      lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
      phase="Solid";
      radius=NNN;
      radius_PT=205;
      radius_covalent=1.99;
      radius_covalent_PT=199;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.5303;
      radii_Slatter=1.85;
      radii_Pyykko=1.73;
      electrical_conductivity=1.3E6;
      electronegativity_vec=1.13;
      hardness_Ghosh=2.1056;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.925;
      electron_affinity_PT=50;
      Miedema_phi_star=3.19;
      Miedema_nws=1.21;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7400;
      boiling_point=3000;
      melting_point=1100;
      vaporization_heat_PT=290;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000011;
      thermal_conductivity=15;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=0.33;
      shear_modulus=18;
      Young_modulus=46;
      bulk_modulus=33;
      Poisson_ratio_PT=0.28;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      // Pm Pettifor linear interpolation
    }
    // [AFLOW]STOP=Promethium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Samarium
    // Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium
    if(ZZ==62) { // Samarium
      Z=ZZ;
      symbol="Sm";
      name="Samarium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*150.4;
      MolarVolume=0.000020449;
      volume=33.9484;
      Miedema_Vm=7.37;
      valence_std=8;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=7.353;
      crystal="rhl";
      CrystalStr_PT="Simple_Trigonal";
      space_group="R_3m";
      space_group_number=166;
      Pearson_coefficient=0.000334686;
      lattice_constant[1]=362.1;lattice_constant[2]=362.1;lattice_constant[3]=2625;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.181;
      radius_PT=238;
      radius_covalent=1.98;
      radius_covalent_PT=198;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.3830;
      radii_Slatter=1.85;
      radii_Pyykko=1.72;
      electrical_conductivity=1.1E6;
      electronegativity_vec=1.17;
      hardness_Ghosh=2.3427;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.047;
      electron_affinity_PT=50;
      Miedema_phi_star=3.20;
      Miedema_nws=1.21;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7380;
      boiling_point=1803;
      melting_point=1072;
      vaporization_heat_PT=175;
      specific_heat_PT=196;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000127;
      thermal_conductivity=13;
      Brinelll_hardness=441;
      Mohs_hardness=1.44;
      Vickers_hardness=412;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=0.02;
      shear_modulus=20;
      Young_modulus=50;
      bulk_modulus=38;
      Poisson_ratio_PT=0.27;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.11E-7;
      Volume_Magnetic_Susceptibility=0.00081618;
      Molar_Magnetic_Susceptibility=1.669E-8;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Sm Pettifor linear interpolation
    }
    // [AFLOW]STOP=Samarium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Europium
    // Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium
    if(ZZ==63) { // Europium
      Z=ZZ;
      symbol="Eu";
      name="Europium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*151.96;
      MolarVolume=0.000028979;
      volume=43.1719;
      Miedema_Vm=7.36;
      valence_std=9;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=5.244;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=4.32857E-05;
      lattice_constant[1]=458.1;lattice_constant[2]=458.1;lattice_constant[3]=458.1;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.204;
      radius_PT=231;
      radius_covalent=1.98;
      radius_covalent_PT=198;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.2615;
      radii_Slatter=1.85;
      radii_Pyykko=1.68;
      electrical_conductivity=1.1E6;
      electronegativity_vec=1.20;
      hardness_Ghosh=2.5798;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.168;
      electron_affinity_PT=50;
      Miedema_phi_star=3.20;
      Miedema_nws=1.21;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7360;
      boiling_point=1527;
      melting_point=822;
      vaporization_heat_PT=175;
      specific_heat_PT=182;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000035;
      thermal_conductivity=14;
      Brinelll_hardness=NNN;
      Mohs_hardness=3.07;
      Vickers_hardness=167;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=2.42;
      shear_modulus=7.9;
      Young_modulus=18;
      bulk_modulus=8.3;
      Poisson_ratio_PT=0.15;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=2.76E-7;
      Volume_Magnetic_Susceptibility=0.0014473;
      Molar_Magnetic_Susceptibility=4.1942E-8;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Eu Pettifor linear interpolation
    }
    // [AFLOW]STOP=Europium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Gadolinium
    // Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium
    if(ZZ==64) { // Gadolinium
      Z=ZZ;
      symbol="Gd";
      name="Gadolinium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*157.25;
      MolarVolume=0.000019903;
      volume=32.5777;
      Miedema_Vm=7.34;
      valence_std=10;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=7.901;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.000127674;
      lattice_constant[1]=363.6;lattice_constant[2]=363.6;lattice_constant[3]=578.26;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.180;
      radius_PT=233;
      radius_covalent=1.96;
      radius_covalent_PT=196;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.1596;
      radii_Slatter=1.80;
      radii_Pyykko=1.69;
      electrical_conductivity=770000;
      electronegativity_vec=1.20;
      hardness_Ghosh=2.8170;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.290;
      electron_affinity_PT=50;
      Miedema_phi_star=3.20;
      Miedema_nws=1.21;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7340;
      boiling_point=3250;
      melting_point=1313;
      vaporization_heat_PT=305;
      specific_heat_PT=240;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=9.4E-6;
      thermal_conductivity=11;
      Brinelll_hardness=NNN;
      Mohs_hardness=5.13;
      Vickers_hardness=570;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=-1.02;
      shear_modulus=22;
      Young_modulus=55;
      bulk_modulus=38;
      Poisson_ratio_PT=0.26;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Ferromagnetic";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=292;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      // Gd Pettifor linear interpolation
    }
    // [AFLOW]STOP=Gadolinium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Terbium
    // Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium
    if(ZZ==65) { // Terbium
      Z=ZZ;
      symbol="Tb";
      name="Terbium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*158.9254;
      MolarVolume=0.000019336;
      volume=32.0200;
      Miedema_Vm=7.20;
      valence_std=11;
      valence_iupac=4;
      valence_PT=3;
      Density_PT=8.219;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=360.1;lattice_constant[2]=360.1;lattice_constant[3]=569.36;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.177;
      radius_PT=225;
      radius_covalent=1.94;
      radius_covalent_PT=194;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.0730;
      radii_Slatter=1.75;
      radii_Pyykko=1.68;
      electrical_conductivity=830000;
      electronegativity_vec=1.10;
      hardness_Ghosh=3.0540;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.411;
      electron_affinity_PT=50;
      Miedema_phi_star=3.21;
      Miedema_nws=1.22;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7320;
      boiling_point=3230;
      melting_point=1356;
      vaporization_heat_PT=295;
      specific_heat_PT=182;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000103;
      thermal_conductivity=11;
      Brinelll_hardness=677;
      Mohs_hardness=2.33;
      Vickers_hardness=863;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=1.36;
      shear_modulus=22;
      Young_modulus=56;
      bulk_modulus=38.7;
      Poisson_ratio_PT=0.26;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=0.0000136;
      Volume_Magnetic_Susceptibility=0.1117784;
      Molar_Magnetic_Susceptibility=2.161385E-6;
      Curie_point=222;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      // Tb Pettifor linear interpolation
    }
    // [AFLOW]STOP=Terbium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Dysprosium
    // Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium
    if(ZZ==66) { // Dysprosium
      Z=ZZ;
      symbol="Dy";
      name="Dysprosium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*162.5;
      MolarVolume=0.000019004;
      volume=31.5096;
      Miedema_Vm=7.12;
      valence_std=12;
      valence_iupac=4;
      valence_PT=3;
      Density_PT=8.551;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=5.20771E-05;
      lattice_constant[1]=359.3;lattice_constant[2]=359.3;lattice_constant[3]=565.37;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.177;
      radius_PT=228;
      radius_covalent=1.92;
      radius_covalent_PT=192;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.9984;
      radii_Slatter=1.75;
      radii_Pyykko=1.67;
      electrical_conductivity=1.1E6;
      electronegativity_vec=1.22;
      hardness_Ghosh=3.2912;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.533;
      electron_affinity_PT=50;
      Miedema_phi_star=3.21;
      Miedema_nws=1.22;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7300;
      boiling_point=2567;
      melting_point=1412;
      vaporization_heat_PT=280;
      specific_heat_PT=167;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.00001;
      thermal_conductivity=11;
      Brinelll_hardness=500;
      Mohs_hardness=1.8;
      Vickers_hardness=540;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=1.06;
      shear_modulus=25;
      Young_modulus=61;
      bulk_modulus=41;
      Poisson_ratio_PT=0.25;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=5.45E-6;
      Volume_Magnetic_Susceptibility=0.046603;
      Molar_Magnetic_Susceptibility=8.85625E-7;
      Curie_point=87;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Dy Pettifor linear interpolation JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Dysprosium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Holmium
    // Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium
    if(ZZ==67) { // Holmium
      Z=ZZ;
      symbol="Ho";
      name="Holmium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*164.9304;
      MolarVolume=0.000018753;
      volume=31.0155;
      Miedema_Vm=7.06;
      valence_std=13;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=8.795;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=2.96961E-32;
      lattice_constant[1]=357.73;lattice_constant[2]=357.73;lattice_constant[3]=561.58;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.176;
      radius_PT=226;
      radius_covalent=1.92;
      radius_covalent_PT=192;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.9335;
      radii_Slatter=1.75;
      radii_Pyykko=1.66;
      electrical_conductivity=1.1E6;
      electronegativity_vec=1.23;
      hardness_Ghosh=3.5283;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.654;
      electron_affinity_PT=50;
      Miedema_phi_star=3.22;
      Miedema_nws=1.22;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7280;
      boiling_point=2700;
      melting_point=1474;
      vaporization_heat_PT=265;
      specific_heat_PT=165;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000112;
      thermal_conductivity=16;
      Brinelll_hardness=746;
      Mohs_hardness=1.65;
      Vickers_hardness=481;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=0.78;
      shear_modulus=26;
      Young_modulus=64;
      bulk_modulus=40;
      Poisson_ratio_PT=0.23;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=5.49E-6;
      Volume_Magnetic_Susceptibility=0.0482845;
      Molar_Magnetic_Susceptibility=9.05467E-7;
      Curie_point=20;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Ho Pettifor linear interpolation
    }
    // [AFLOW]STOP=Holmium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Erbium
    // Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium
    if(ZZ==68) { // Erbium
      Z=ZZ;
      symbol="Er";
      name="Erbium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*167.26;
      MolarVolume=0.000018449;
      volume=30.5431;
      Miedema_Vm=6.98;
      valence_std=14;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=9.066;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=7.24618E-05;
      lattice_constant[1]=355.88;lattice_constant[2]=355.88;lattice_constant[3]=558.74;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.175;
      radius_PT=226;
      radius_covalent=1.89;
      radius_covalent_PT=189;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8765;
      radii_Slatter=1.75;
      radii_Pyykko=1.65;
      electrical_conductivity=1.2E6;
      electronegativity_vec=1.24;
      hardness_Ghosh=3.7655;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.776;
      electron_affinity_PT=50;
      Miedema_phi_star=3.22;
      Miedema_nws=1.23;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7260;
      boiling_point=2868;
      melting_point=1497;
      vaporization_heat_PT=285;
      specific_heat_PT=168;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000122;
      thermal_conductivity=15;
      Brinelll_hardness=814;
      Mohs_hardness=1.97;
      Vickers_hardness=588;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=0.54;
      shear_modulus=28;
      Young_modulus=70;
      bulk_modulus=44;
      Poisson_ratio_PT=0.24;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=3.77E-6;
      Volume_Magnetic_Susceptibility=0.0341788;
      Molar_Magnetic_Susceptibility=6.30566E-7;
      Curie_point=32;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Er Pettifor linear interpolation
    }
    // [AFLOW]STOP=Erbium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Thulium
    // Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium
    if(ZZ==69) { // Thulium
      Z=ZZ;
      symbol="Tm";
      name="Thulium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*168.9342;
      MolarVolume=0.000018126;
      volume=30.0016;
      Miedema_Vm=6.90;
      valence_std=15;
      valence_iupac=4;
      valence_PT=3;
      Density_PT=9.32;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=353.75;lattice_constant[2]=353.75;lattice_constant[3]=555.46;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.174;
      radius_PT=222;
      radius_covalent=1.90;
      radius_covalent_PT=190;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8261;
      radii_Slatter=1.75;
      radii_Pyykko=1.64;
      electrical_conductivity=1.4E6;
      electronegativity_vec=1.25;
      hardness_Ghosh=4.0026;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.897;
      electron_affinity_PT=50;
      Miedema_phi_star=3.22;
      Miedema_nws=1.23;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7240;
      boiling_point=1950;
      melting_point=1545;
      vaporization_heat_PT=250;
      specific_heat_PT=160;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000133;
      thermal_conductivity=17;
      Brinelll_hardness=471;
      Mohs_hardness=1.77;
      Vickers_hardness=520;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=0.32;
      shear_modulus=31;
      Young_modulus=74;
      bulk_modulus=45;
      Poisson_ratio_PT=0.21;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.99E-6;
      Volume_Magnetic_Susceptibility=0.0185488;
      Molar_Magnetic_Susceptibility=3.36179E-7;
      Curie_point=25;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Tm Pettifor linear interpolation JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Thulium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Ytterbium
    // Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium
    if(ZZ==70) { // Ytterbium
      Z=ZZ;
      symbol="Yb";
      name="Ytterbium";
      Period=6;
      Group=NNN;
      Series="Lanthanide";
      Block="f";
      mass=AMU2KILOGRAM*173.04;
      MolarVolume=0.000026339;
      volume=39.4395;
      Miedema_Vm=6.86;
      valence_std=16;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=6.57;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=8.54557E-05;
      lattice_constant[1]=548.47;lattice_constant[2]=548.47;lattice_constant[3]=548.47;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.193;
      radius_PT=222;
      radius_covalent=1.87;
      radius_covalent_PT=187;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.7812;
      radii_Slatter=1.75;
      radii_Pyykko=1.70;
      electrical_conductivity=3.6E6;
      electronegativity_vec=1.10;
      hardness_Ghosh=4.2395;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.019;
      electron_affinity_PT=50;
      Miedema_phi_star=3.22;
      Miedema_nws=1.23;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7220;
      boiling_point=1196;
      melting_point=819;
      vaporization_heat_PT=160;
      specific_heat_PT=154;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000263;
      thermal_conductivity=39;
      Brinelll_hardness=343;
      Mohs_hardness=NNN;
      Vickers_hardness=206;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=3.27;
      shear_modulus=10;
      Young_modulus=24;
      bulk_modulus=31;
      Poisson_ratio_PT=0.21;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=5.9E-9;
      Volume_Magnetic_Susceptibility=0.0000388;
      Molar_Magnetic_Susceptibility=1.02E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Yb Pettifor linear interpolation
    }
    // [AFLOW]STOP=Ytterbium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Lutetium
    // Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium
    if(ZZ==71) { // Lutetium
      Z=ZZ;
      symbol="Lu";
      name="Lutetium";
      Period=6;
      Group=3;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*174.967;
      MolarVolume=0.000017779;
      volume=29.3515;
      Miedema_Vm=6.81;
      valence_std=17;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=9.841;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=8.27273E-07;
      lattice_constant[1]=350.31;lattice_constant[2]=350.31;lattice_constant[3]=555.09;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.173;
      radius_PT=217;
      radius_covalent=1.87;
      radius_covalent_PT=187;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.7409;
      radii_Slatter=1.75;
      radii_Pyykko=1.62;
      electrical_conductivity=1.8E6;
      electronegativity_vec=1.27;
      hardness_Ghosh=4.4766;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.140;
      electron_affinity_PT=50;
      Miedema_phi_star=3.22;
      Miedema_nws=1.24;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0.7200;
      boiling_point=3402;
      melting_point=1663;
      vaporization_heat_PT=415;
      specific_heat_PT=154;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.00001;
      thermal_conductivity=16;
      Brinelll_hardness=893;
      Mohs_hardness=2.6;
      Vickers_hardness=1160;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=3.64;
      shear_modulus=27;
      Young_modulus=67;
      bulk_modulus=48;
      Poisson_ratio_PT=0.26;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.2E-9;
      Volume_Magnetic_Susceptibility=0.0000118;
      Molar_Magnetic_Susceptibility=2.1E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Lu
    }
    // [AFLOW]STOP=Lutetium
    // ********************************************************************************************************************************************************

    // d-electron systems: transition metals
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Hafnium
    // Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium
    if(ZZ==72) { // Hafnium
      Z=ZZ;
      symbol="Hf";
      name="Hafnium";
      Period=6;
      Group=4;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*178.49;
      MolarVolume=0.0000134102;
      volume=22.0408;
      Miedema_Vm=5.6;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=13.31;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=5.25384E-05;
      lattice_constant[1]=319.64;lattice_constant[2]=319.64;lattice_constant[3]=505.11;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.159;
      radius_PT=208;
      radius_covalent=1.75;
      radius_covalent_PT=175;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.7056;
      radii_Slatter=1.55;
      radii_Pyykko=1.52;
      electrical_conductivity=3.3E6;
      electronegativity_vec=1.30;
      hardness_Ghosh=4.7065;
      electronegativity_Pearson=3.8;
      electronegativity_Ghosh=6.258;
      electron_affinity_PT=0;
      Miedema_phi_star=3.55;
      Miedema_nws=1.43;
      Miedema_gamma_s=2200;
      Pettifor_scale=0.775;
      boiling_point=4603;
      melting_point=2233;
      vaporization_heat_PT=630;
      specific_heat_PT=144;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=5.9E-6;
      thermal_conductivity=23;
      Brinelll_hardness=1700;
      Mohs_hardness=5.5;
      Vickers_hardness=1760;
      Hardness_Pearson=3.00;
      Hardness_Putz=NNN;
      Hardness_RB=3.94;
      shear_modulus=30;
      Young_modulus=78;
      bulk_modulus=110;
      Poisson_ratio_PT=0.37;
      Miedema_BVm=15.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=5.3E-9;
      Volume_Magnetic_Susceptibility=0.0000705;
      Molar_Magnetic_Susceptibility=9.46E-10;
      Curie_point=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=3400;
      HHIR=2600;
      /*xray_scatt=NNN;*/
      //Hf
    }
    // [AFLOW]STOP=Hafnium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Tantalum
    // Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum
    if(ZZ==73) { // Tantalum
      Z=ZZ;
      symbol="Ta";
      name="Tantalum";
      Period=6;
      Group=5;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*180.9479;
      MolarVolume=0.0000108677;
      volume=18.1100;
      Miedema_Vm=4.9;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      Density_PT=16.65;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=3.66845E-09;
      lattice_constant[1]=330.13;lattice_constant[2]=330.13;lattice_constant[3]=330.13;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.147;
      radius_PT=200;
      radius_covalent=1.70;
      radius_covalent_PT=170;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.6716;
      radii_Slatter=1.45;
      radii_Pyykko=1.46;
      electrical_conductivity=7.7E6;
      electronegativity_vec=1.50;
      hardness_Ghosh=4.9508;
      electronegativity_Pearson=4.11;
      electronegativity_Ghosh=6.383;
      electron_affinity_PT=31;
      Miedema_phi_star=4.05;
      Miedema_nws=1.63;
      Miedema_gamma_s=3050;
      Pettifor_scale=0.83;
      boiling_point=5458;
      melting_point=3017;
      vaporization_heat_PT=736;
      specific_heat_PT=140;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.3E-6;
      thermal_conductivity=57;
      Brinelll_hardness=800;
      Mohs_hardness=6.5;
      Vickers_hardness=873;
      Hardness_Pearson=3.79;
      Hardness_Putz=NNN;
      Hardness_RB=1.75;
      shear_modulus=67;
      Young_modulus=186;
      bulk_modulus=200;
      Poisson_ratio_PT=0.34;
      Miedema_BVm=22.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.07E-8;
      Volume_Magnetic_Susceptibility=0.0001782;
      Molar_Magnetic_Susceptibility=1.936E-9;
      Curie_point=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=2300;
      HHIR=4800;
      /*xray_scatt=NNN;*/
      //Ta
    }
    // [AFLOW]STOP=Tantalum
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Tungsten
    // Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten
    if(ZZ==74) { // Tungsten
      Z=ZZ;
      symbol="W";
      name="Tungsten";
      Period=6;
      Group=6;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*183.85;
      MolarVolume=9.5501E-6;
      volume=15.9387;
      Miedema_Vm=4.5;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      Density_PT=19.25;
      crystal="bcc";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=6.96679E-05;
      lattice_constant[1]=316.52;lattice_constant[2]=316.52;lattice_constant[3]=316.52;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.137;
      radius_PT=193;
      radius_covalent=1.62;
      radius_covalent_PT=162;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.6416;
      radii_Slatter=1.35;
      radii_Pyykko=1.37;
      electrical_conductivity=2E7;
      electronegativity_vec=2.36;
      hardness_Ghosh=5.1879;
      electronegativity_Pearson=4.40;
      electronegativity_Ghosh=6.505;
      electron_affinity_PT=78.6;
      Miedema_phi_star=4.80;
      Miedema_nws=1.81;
      Miedema_gamma_s=3300;
      Pettifor_scale=0.885;
      boiling_point=5555;
      melting_point=3422;
      vaporization_heat_PT=800;
      specific_heat_PT=132;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=4.5E-6;
      thermal_conductivity=170;
      Brinelll_hardness=2570;
      Mohs_hardness=7.5;
      Vickers_hardness=3430;
      Hardness_Pearson=3.58;
      Hardness_Putz=NNN;
      Hardness_RB=1.23;
      shear_modulus=161;
      Young_modulus=411;
      bulk_modulus=310;
      Poisson_ratio_PT=0.28;
      Miedema_BVm=31.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=4.59E-9;
      Volume_Magnetic_Susceptibility=0.0000884;
      Molar_Magnetic_Susceptibility=8.44E-10;
      Curie_point=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=7000;
      HHIR=4300;
      /*xray_scatt=NNN;*/
      //W
    }
    // [AFLOW]STOP=Tungsten
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Rhenium
    // Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium
    if(ZZ==75) { // Rhenium
      Z=ZZ;
      symbol="Re";
      name="Rhenium";
      Period=6;
      Group=7;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*186.2;
      MolarVolume=8.85856E-6;
      volume=14.8941;
      Miedema_Vm=4.3;
      valence_std=7;
      valence_iupac=7;
      valence_PT=7;
      Density_PT=21.02;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=2.70849E-05;
      lattice_constant[1]=276.1;lattice_constant[2]=276.1;lattice_constant[3]=445.6;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.138;
      radius_PT=188;
      radius_covalent=1.51;
      radius_covalent_PT=151;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.6141;
      radii_Slatter=1.35;
      radii_Pyykko=1.31;
      electrical_conductivity=5.6E6;
      electronegativity_vec=1.90;
      hardness_Ghosh=5.4256;
      electronegativity_Pearson=4.02;
      electronegativity_Ghosh=6.626;
      electron_affinity_PT=14.5;
      Miedema_phi_star=5.40;
      Miedema_nws=1.86;
      Miedema_gamma_s=3650;
      Pettifor_scale=0.94;
      boiling_point=5596;
      melting_point=3186;
      vaporization_heat_PT=705;
      specific_heat_PT=137;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.2E-6;
      thermal_conductivity=48;
      Brinelll_hardness=1320;
      Mohs_hardness=7;
      Vickers_hardness=2450;
      Hardness_Pearson=3.87;
      Hardness_Putz=NNN;
      Hardness_RB=2.13;
      shear_modulus=178;
      Young_modulus=463;
      bulk_modulus=370;
      Poisson_ratio_PT=0.3;
      Miedema_BVm=33.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=4.56E-9;
      Volume_Magnetic_Susceptibility=0.0000959;
      Molar_Magnetic_Susceptibility=8.49E-10;
      Curie_point=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=3300;
      HHIR=3300;
      /*xray_scatt=NNN;*/
      //Re
    }
    // [AFLOW]STOP=Rhenium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Osmium
    // Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium
    if(ZZ==76) { // Osmium
      Z=ZZ;
      symbol="Os";
      name="Osmium";
      Period=6;
      Group=8;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*190.2;
      MolarVolume=8.421E-6;
      volume=14.2403;
      Miedema_Vm=4.2;
      valence_std=8;
      valence_iupac=8;
      valence_PT=6;
      Density_PT=22.59;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=7.45234E-05;
      lattice_constant[1]=273.44;lattice_constant[2]=273.44;lattice_constant[3]=431.73;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.135;
      radius_PT=185;
      radius_covalent=1.44;
      radius_covalent_PT=144;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.5890;
      radii_Slatter=1.30;
      radii_Pyykko=1.29;
      electrical_conductivity=1.2E7;
      electronegativity_vec=2.20;
      hardness_Ghosh=5.6619;
      electronegativity_Pearson=4.9;
      electronegativity_Ghosh=6.748;
      electron_affinity_PT=106.1;
      Miedema_phi_star=5.40;
      Miedema_nws=1.85;
      Miedema_gamma_s=3500;
      Pettifor_scale=0.995;
      boiling_point=5012;
      melting_point=3033;
      vaporization_heat_PT=630;
      specific_heat_PT=130;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=5.1E-6;
      thermal_conductivity=87;
      Brinelll_hardness=3920;
      Mohs_hardness=7;
      Vickers_hardness=4137.063913415;
      Hardness_Pearson=3.80;
      Hardness_Putz=NNN;
      Hardness_RB=1.72;
      shear_modulus=222;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=0.25;
      Miedema_BVm=35.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=6E-10;
      Volume_Magnetic_Susceptibility=0.000014;
      Molar_Magnetic_Susceptibility=1.1E-10;
      Curie_point=NNN;
      color_PT="SLATEGRAY";
      refractive_index=NNN;
      HHIP=5500;
      HHIR=9100;
      /*xray_scatt=NNN;*/
      //Os JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Osmium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Iridium
    // Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium
    if(ZZ==77) { // Iridium
      Z=ZZ;
      symbol="Ir";
      name="Iridium";
      Period=6;
      Group=9;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*192.22;
      MolarVolume=8.5203E-6;
      volume=14.5561;
      Miedema_Vm=4.2;
      valence_std=9;
      valence_iupac=8;
      valence_PT=6;
      Density_PT=22.56;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=2.53787E-05;
      lattice_constant[1]=383.9;lattice_constant[2]=383.9;lattice_constant[3]=383.9;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.135;
      radius_PT=180;
      radius_covalent=1.41;
      radius_covalent_PT=141;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.5657;
      radii_Slatter=1.35;
      radii_Pyykko=1.22;
      electrical_conductivity=2.1E7;
      electronegativity_vec=2.20;
      hardness_Ghosh=5.9000;
      electronegativity_Pearson=5.4;
      electronegativity_Ghosh=6.831;
      electron_affinity_PT=151;
      Miedema_phi_star=5.55;
      Miedema_nws=1.83;
      Miedema_gamma_s=3100;
      Pettifor_scale=1.05;
      boiling_point=4428;
      melting_point=2466;
      vaporization_heat_PT=560;
      specific_heat_PT=131;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.4E-6;
      thermal_conductivity=150;
      Brinelll_hardness=1670;
      Mohs_hardness=6.5;
      Vickers_hardness=1760;
      Hardness_Pearson=3.80;
      Hardness_Putz=NNN;
      Hardness_RB=1.27;
      shear_modulus=210;
      Young_modulus=528;
      bulk_modulus=320;
      Poisson_ratio_PT=0.26;
      Miedema_BVm=25.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.67E-9;
      Volume_Magnetic_Susceptibility=0.0000377;
      Molar_Magnetic_Susceptibility=3.21E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=5500;
      HHIR=9100;
      /*xray_scatt=NNN;*/
      //Ir JUNKAI CHANGED VALENCE
    }
    // [AFLOW]STOP=Iridium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Platinum
    // Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum
    if(ZZ==78) { // Platinum
      Z=ZZ;
      symbol="Pt";
      name="Platinum";
      Period=6;
      Group=10;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*195.09;
      MolarVolume=9.0948E-6;
      volume=15.7298;
      Miedema_Vm=4.4;
      valence_std=10;
      valence_iupac=6;
      valence_PT=6;
      Density_PT=21.45;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=3.39206E-05;
      lattice_constant[1]=392.42;lattice_constant[2]=392.42;lattice_constant[3]=392.42;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.138;
      radius_PT=177;
      radius_covalent=1.36;
      radius_covalent_PT=136;
      radius_VanDerWaals_PT=175;
      radii_Ghosh08=0.5443;
      radii_Slatter=1.35;
      radii_Pyykko=1.23;
      electrical_conductivity=9.4E6;
      electronegativity_vec=2.28;
      hardness_Ghosh=6.1367;
      electronegativity_Pearson=5.6;
      electronegativity_Ghosh=6.991;
      electron_affinity_PT=205.3;
      Miedema_phi_star=5.65;
      Miedema_nws=1.78;
      Miedema_gamma_s=2550;
      Pettifor_scale=1.105;
      boiling_point=3825;
      melting_point=1768.3;
      vaporization_heat_PT=490;
      specific_heat_PT=133;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=8.9E-6;
      thermal_conductivity=71;
      Brinelll_hardness=392;
      Mohs_hardness=3.5;
      Vickers_hardness=549;
      Hardness_Pearson=3.50;
      Hardness_Putz=NNN;
      Hardness_RB=3.5;
      shear_modulus=61;
      Young_modulus=168;
      bulk_modulus=230;
      Poisson_ratio_PT=0.38;
      Miedema_BVm=18.0;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=1.22E-8;
      Volume_Magnetic_Susceptibility=0.0002573;
      Molar_Magnetic_Susceptibility=2.38E-9;
      Curie_point=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=5500;
      HHIR=9100;
      /*xray_scatt=NNN;*/
      //Pt
    }
    // [AFLOW]STOP=Platinum
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Gold
    // Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold
    if(ZZ==79) { // Gold
      Z=ZZ;
      symbol="Au";
      name="Gold";
      Period=6;
      Group=11;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*196.9665;
      MolarVolume=0.00001021;
      volume=18.1904;
      Miedema_Vm=4.7;
      valence_std=11;
      valence_iupac=5;
      valence_PT=5;
      Density_PT=19.3;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=2.08217E-32;
      lattice_constant[1]=407.82;lattice_constant[2]=407.82;lattice_constant[3]=407.82;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.144;
      radius_PT=174;
      radius_covalent=1.36;
      radius_covalent_PT=136;
      radius_VanDerWaals_PT=166;
      radii_Ghosh08=0.5244;
      radii_Slatter=1.35;
      radii_Pyykko=1.24;
      electrical_conductivity=4.5E7;
      electronegativity_vec=2.54;
      hardness_Ghosh=6.3741;
      electronegativity_Pearson=5.77;
      electronegativity_Ghosh=7.112;
      electron_affinity_PT=222.8;
      Miedema_phi_star=5.15;
      Miedema_nws=1.57;
      Miedema_gamma_s=1550;
      Pettifor_scale=1.16;
      boiling_point=2856;
      melting_point=1064.18;
      vaporization_heat_PT=330;
      specific_heat_PT=129.1;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000142;
      thermal_conductivity=320;
      Brinelll_hardness=2450;
      Mohs_hardness=2.5;
      Vickers_hardness=216;
      Hardness_Pearson=3.46;
      Hardness_Putz=NNN;
      Hardness_RB=3.44;
      shear_modulus=27;
      Young_modulus=78;
      bulk_modulus=220;
      Poisson_ratio_PT=0.44;
      Miedema_BVm=18.0;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-1.78E-9;
      Volume_Magnetic_Susceptibility=-0.0000344;
      Molar_Magnetic_Susceptibility=-3.51E-10;
      Curie_point=NNN;
      color_PT="GOLD";
      refractive_index=NNN;
      HHIP=1100;
      HHIR=1000;
      xray_scatt=74.99;
      //Au
    }
    // [AFLOW]STOP=Gold
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Mercury
    // Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury
    if(ZZ==80) { // Mercury
      Z=ZZ;
      symbol="Hg";
      name="Mercury";
      Period=6;
      Group=12;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*200.59;
      MolarVolume=0.0000148213;
      volume=29.7156;
      Miedema_Vm=5.8;
      valence_std=12;
      valence_iupac=4;
      valence_PT=2;
      Density_PT=13.534;
      crystal="rhl";
      CrystalStr_PT="Simple_Trigonal";
      space_group="R_3m";
      space_group_number=166;
      Pearson_coefficient=6.52519E-05;
      lattice_constant[1]=300.5;lattice_constant[2]=300.5;lattice_constant[3]=300.5;
      lattice_angle[1]=1.23081;lattice_angle[2]=1.23081;lattice_angle[3]=1.23081;
      phase="Liquid";
      radius=0.150;
      radius_PT=171;
      radius_covalent=1.32;
      radius_covalent_PT=132;
      radius_VanDerWaals_PT=155;
      radii_Ghosh08=0.5060;
      radii_Slatter=1.50;
      radii_Pyykko=1.33;
      electrical_conductivity=1E6;
      electronegativity_vec=2.00;
      hardness_Ghosh=6.6103;
      electronegativity_Pearson=4.91;
      electronegativity_Ghosh=7.233;
      electron_affinity_PT=0;
      Miedema_phi_star=4.20;
      Miedema_nws=1.24;
      Miedema_gamma_s=610;
      Pettifor_scale=1.32;
      boiling_point=356.73;
      melting_point=-38.83;
      vaporization_heat_PT=59.2;
      specific_heat_PT=139.5;
      critical_Pressure=1698;
      critical_Temperature_PT=1750;
      thermal_expansion=0.000181;
      thermal_conductivity=8.3;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=5.54;
      Hardness_Putz=NNN;
      Hardness_RB=5.29;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=25;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=4.0;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-2.1E-9;
      Volume_Magnetic_Susceptibility=-0.0000284;
      Molar_Magnetic_Susceptibility=-4.21E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=1.000933;
      HHIP=5500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Hg
    }
    // [AFLOW]STOP=Mercury
    // ********************************************************************************************************************************************************

    // p-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Thallium
    // Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium
    if(ZZ==81) { // Thallium
      Z=ZZ;
      symbol="Tl";
      name="Thallium";
      Period=6;
      Group=13;
      Series="PoorMetal";
      Block="p";
      mass=AMU2KILOGRAM*204.37;
      MolarVolume=0.0000172473;
      volume=31.0721;
      Miedema_Vm=6.6;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=11.85;
      crystal="hcp";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=1.99659E-05;
      lattice_constant[1]=345.66;lattice_constant[2]=345.66;lattice_constant[3]=552.48;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=0.171;
      radius_PT=156;
      radius_covalent=1.45;
      radius_covalent_PT=145;
      radius_VanDerWaals_PT=196;
      radii_Ghosh08=1.8670;
      radii_Slatter=1.90;
      radii_Pyykko=1.44;
      electrical_conductivity=6.7E6;
      electronegativity_vec=1.62;
      hardness_Ghosh=1.7043;
      electronegativity_Pearson=3.2;
      electronegativity_Ghosh=4.719;
      electron_affinity_PT=19.2;
      Miedema_phi_star=3.90;
      Miedema_nws=1.12;
      Miedema_gamma_s=610;
      Pettifor_scale=1.56;
      boiling_point=1473;
      melting_point=304;
      vaporization_heat_PT=165;
      specific_heat_PT=129;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000299;
      thermal_conductivity=46;
      Brinelll_hardness=26.4;
      Mohs_hardness=1.2;
      Vickers_hardness=NNN;
      Hardness_Pearson=2.90;
      Hardness_Putz=NNN;
      Hardness_RB=2.69;
      shear_modulus=2.8;
      Young_modulus=8;
      bulk_modulus=43;
      Poisson_ratio_PT=0.45;
      Miedema_BVm=6.2;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-3E-9;
      Volume_Magnetic_Susceptibility=-0.0000356;
      Molar_Magnetic_Susceptibility=-6.13E-10;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=6500;
      HHIR=6500;
      /*xray_scatt=NNN;*/
      //Tl
    }
    // [AFLOW]STOP=Thallium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Lead
    // Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead
    if(ZZ==82) { // Lead
      Z=ZZ;
      symbol="Pb";
      name="Lead";
      Period=6;
      Group=14;
      Series="PoorMetal";
      Block="p";
      mass=AMU2KILOGRAM*207.2;
      MolarVolume=0.000018272;
      volume=31.6649;
      Miedema_Vm=6.9;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=11.34;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=1.94378E-05;
      lattice_constant[1]=495.08;lattice_constant[2]=495.08;lattice_constant[3]=495.08;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.175;
      radius_PT=154;
      radius_covalent=1.46;
      radius_covalent_PT=146;
      radius_VanDerWaals_PT=202;
      radii_Ghosh08=1.6523;
      radii_Slatter=NNN;
      radii_Pyykko=1.44;
      electrical_conductivity=4.8E6;
      electronegativity_vec=2.33;
      hardness_Ghosh=1.9414;
      electronegativity_Pearson=3.90;
      electronegativity_Ghosh=4.841;
      electron_affinity_PT=35.1;
      Miedema_phi_star=4.10;
      Miedema_nws=1.15;
      Miedema_gamma_s=610;
      Pettifor_scale=1.80;
      boiling_point=1749;
      melting_point=327.46;
      vaporization_heat_PT=178;
      specific_heat_PT=127;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000289;
      thermal_conductivity=35;
      Brinelll_hardness=38.3;
      Mohs_hardness=1.5;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.53;
      Hardness_Putz=NNN;
      Hardness_RB=3.02;
      shear_modulus=5.6;
      Young_modulus=16;
      bulk_modulus=46;
      Poisson_ratio_PT=0.44;
      Miedema_BVm=7.9;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-1.5E-9;
      Volume_Magnetic_Susceptibility=-0.000017;
      Molar_Magnetic_Susceptibility=-3.11E-10;
      Curie_point=NNN;
      color_PT="SLATEGRAY";
      refractive_index=NNN;
      HHIP=2700;
      HHIR=1800;
      /*xray_scatt=NNN;*/
      //Pb
    }
    // [AFLOW]STOP=Lead
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Bismuth
    // Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth
    if(ZZ==83) { // Bismuth
      Z=ZZ;
      symbol="Bi";
      name="Bismuth";
      Period=6;
      Group=15;
      Series="PoorMetal";
      Block="p";
      mass=AMU2KILOGRAM*208.9804;
      MolarVolume=0.000021368;
      volume=31.5691;
      Miedema_Vm=7.2;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      Density_PT=9.78;
      crystal="rhl";
      CrystalStr_PT="Base-centered_Monoclinic";
      space_group="C12/m1";
      space_group_number=12;
      Pearson_coefficient=0.0;
      lattice_constant[1]=667.4;lattice_constant[2]=611.7;lattice_constant[3]=330.4;
      lattice_angle[1]=PI/2;lattice_angle[2]=1.925622;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.182;
      radius_PT=143;
      radius_covalent=1.48;
      radius_covalent_PT=148;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.4818;
      radii_Slatter=1.60;
      radii_Pyykko=1.51;
      electrical_conductivity=770000;
      electronegativity_vec=2.02;
      hardness_Ghosh=2.1785;
      electronegativity_Pearson=4.69;
      electronegativity_Ghosh=4.962;
      electron_affinity_PT=91.2;
      Miedema_phi_star=4.15;
      Miedema_nws=1.16;
      Miedema_gamma_s=550;
      Pettifor_scale=2.04;
      boiling_point=1564;
      melting_point=271.3;
      vaporization_heat_PT=160;
      specific_heat_PT=122;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000134;
      thermal_conductivity=8;
      Brinelll_hardness=94.2;
      Mohs_hardness=2.25;
      Vickers_hardness=NNN;
      Hardness_Pearson=3.74;
      Hardness_Putz=NNN;
      Hardness_RB=4.14;
      shear_modulus=12;
      Young_modulus=32;
      bulk_modulus=31;
      Poisson_ratio_PT=0.33;
      Miedema_BVm=6.7;
      Magnetic_Type_PT="Diamagnetic";
      Mass_Magnetic_Susceptibility=-1.7E-8;
      Volume_Magnetic_Susceptibility=-0.00017;
      Molar_Magnetic_Susceptibility=-3.6E-9;
      Curie_point=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Bi
    }
    // [AFLOW]STOP=Bismuth
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Polonium
    // Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium
    if(ZZ==84) { // Polonium
      Z=ZZ;
      symbol="Po";
      name="Polonium";
      Period=6;
      Group=16;
      Series="Chalcogen";
      Block="p";
      mass=AMU2KILOGRAM*209.98;
      MolarVolume=0.00002272727272727;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      Density_PT=9.196;
      crystal="sc";
      CrystalStr_PT="Simple_Cubic";
      space_group="Pm-3m";
      space_group_number=221;
      Pearson_coefficient=0.0;
      lattice_constant[1]=335.9;lattice_constant[2]=335.9;lattice_constant[3]=335.9;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.140;
      radius_PT=135;
      radius_covalent=1.40;
      radius_covalent_PT=140;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.3431;
      radii_Slatter=1.90;
      radii_Pyykko=1.45;
      electrical_conductivity=2.3E6;
      electronegativity_vec=2.00;
      hardness_Ghosh=2.4158;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.084;
      electron_affinity_PT=183.3;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=2.28;
      boiling_point=962;
      melting_point=254;
      vaporization_heat_PT=100;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=NNN;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=3.28;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Po
    }
    // [AFLOW]STOP=Polonium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Astatine
    // Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine
    if(ZZ==85) { // Astatine
      Z=ZZ;
      symbol="At";
      name="Astatine";
      Period=6;
      Group=17;
      Series="Halogen";
      Block="p";
      mass=AMU2KILOGRAM*210;
      MolarVolume=NNN;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=7;
      valence_iupac=7;
      valence_PT=7;
      Density_PT=NNN;
      crystal="nnn";
      CrystalStr_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      Pearson_coefficient=0.0;
      lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
      lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
      phase="Solid";
      radius=NNN;
      radius_PT=127;
      radius_covalent=1.50;
      radius_covalent_PT=150;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.2283;
      radii_Slatter=NNN;
      radii_Pyykko=1.47;
      electrical_conductivity=NNN;
      electronegativity_vec=2.20;
      hardness_Ghosh=2.6528;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.206;
      electron_affinity_PT=270.1;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=2.52;
      boiling_point=NNN;
      melting_point=302;
      vaporization_heat_PT=40;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=2;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=3.57;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //At
    }
    // [AFLOW]STOP=Astatine
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Radon
    // Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon
    if(ZZ==86) { // Radon
      Z=ZZ;
      symbol="Rn";
      name="Radon";
      Period=6;
      Group=18;
      Series="NobleGas";
      Block="p";
      mass=AMU2KILOGRAM*222;
      MolarVolume=0.02281603288798;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=0;
      valence_iupac=6;
      valence_PT=6;
      Density_PT=97.3E-4;
      crystal="fcc";
      CrystalStr_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      Pearson_coefficient=0.0;
      lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
      lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
      phase="Gas";
      radius=NNN;
      radius_PT=120;
      radius_covalent=1.50;
      radius_covalent_PT=150;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.1315;
      radii_Slatter=NNN;
      radii_Pyykko=1.42;
      electrical_conductivity=NNN;
      electronegativity_vec=2.2;
      hardness_Ghosh=2.8900;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.327;
      electron_affinity_PT=0;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=-61.7;
      melting_point=-71;
      vaporization_heat_PT=17;
      specific_heat_PT=93.65;
      critical_Pressure=61.98;
      critical_Temperature_PT=377;
      thermal_expansion=NNN;
      thermal_conductivity=0.00361;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=7.69;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="COLORLESS";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Rn
    }
    // [AFLOW]STOP=Radon
    // ********************************************************************************************************************************************************

    // ROW7
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Francium
    // Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium
    if(ZZ==87) { // Francium
      Z=ZZ;
      symbol="Fr";
      name="Francium";
      Period=7;
      Group=1;
      Series="AlkaliMetal";
      Block="s";
      mass=AMU2KILOGRAM*223.02;
      MolarVolume=NNN;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      Density_PT=NNN;
      crystal="bcc";
      CrystalStr_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      Pearson_coefficient=0.0;
      lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
      lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=2.60;
      radius_covalent_PT=260;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=4.4479;
      radii_Slatter=NNN;
      radii_Pyykko=2.23;
      electrical_conductivity=NNN;
      electronegativity_vec=0.70;
      hardness_Ghosh=0.9882;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=2.376;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=NNN;
      melting_point=NNN;
      vaporization_heat_PT=64;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=NNN;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Fr
    }
    // [AFLOW]STOP=Francium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Radium
    // Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium
    if(ZZ==88) { // Radium
      Z=ZZ;
      symbol="Ra";
      name="Radium";
      Period=7;
      Group=2;
      Series="AlkalineEarthMetal";
      Block="s";
      mass=AMU2KILOGRAM*226.0254;
      MolarVolume=0.0000452;
      volume=-1.0000;
      Miedema_Vm=NNN;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      Density_PT=5;
      crystal="bct";
      CrystalStr_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      Pearson_coefficient=0.0;
      lattice_constant[1]=514.8;lattice_constant[2]=514.8;lattice_constant[3]=514.8;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=2.21;
      radius_covalent_PT=221;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=3.4332;
      radii_Slatter=2.15;
      radii_Pyykko=2.01;
      electrical_conductivity=1E6;
      electronegativity_vec=0.89;
      hardness_Ghosh=1.2819;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=2.664;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=1737;
      melting_point=700;
      vaporization_heat_PT=125;
      specific_heat_PT=92;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=19;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Ra
    }
    // [AFLOW]STOP=Radium
    // ********************************************************************************************************************************************************

    // d-electron systems: transition metals
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Actinium
    // Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium
    if(ZZ==89) { // Actinium
      Z=ZZ;
      symbol="Ac";
      name="Actinium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*227.03;
      MolarVolume=0.00002254220456802;
      volume=45.2437;
      Miedema_Vm=NNN;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=10.07;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.0;
      lattice_constant[1]=567;lattice_constant[2]=567;lattice_constant[3]=567;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=2.15;
      radius_covalent_PT=215;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=3.2615;
      radii_Slatter=1.95;
      radii_Pyykko=1.86;
      electrical_conductivity=NNN;
      electronegativity_vec=1.10;
      hardness_Ghosh=1.3497;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=2.730;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=3200;
      melting_point=1050;
      vaporization_heat_PT=400;
      specific_heat_PT=120;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=12;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Ac
    }
    // [AFLOW]STOP=Actinium
    // ********************************************************************************************************************************************************

    // actinidies
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Thorium
    // Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium
    if(ZZ==90) { // Thorium
      Z=ZZ;
      symbol="Th";
      name="Thorium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*232.0381;
      MolarVolume=0.0000197917;
      volume=31.9586;
      Miedema_Vm=7.3;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=11.724;
      crystal="fcc";
      CrystalStr_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      Pearson_coefficient=0.0;
      lattice_constant[1]=508.42;lattice_constant[2]=508.42;lattice_constant[3]=508.42;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.180;
      radius_PT=NNN;
      radius_covalent=2.06;
      radius_covalent_PT=206;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=3.1061;
      radii_Slatter=1.80;
      radii_Pyykko=1.75;
      electrical_conductivity=6.7E6;
      electronegativity_vec=1.30;
      hardness_Ghosh=1.4175;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=2.796;
      electron_affinity_PT=NNN;
      Miedema_phi_star=3.30;
      Miedema_nws=1.28;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=4820;
      melting_point=1750;
      vaporization_heat_PT=530;
      specific_heat_PT=118;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000011;
      thermal_conductivity=54;
      Brinelll_hardness=400;
      Mohs_hardness=3;
      Vickers_hardness=350;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=31;
      Young_modulus=79;
      bulk_modulus=54;
      Poisson_ratio_PT=0.27;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=7.2E-9;
      Volume_Magnetic_Susceptibility=0.000084;
      Molar_Magnetic_Susceptibility=1.7E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      xray_scatt=86.64;
      //Th
    }
    // [AFLOW]STOP=Thorium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Protoactinium
    // Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium
    if(ZZ==91) { // Protoactinium
      Z=ZZ;
      symbol="Pa";
      name="Protoactinium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*231.04;
      MolarVolume=0.0000150316;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      Density_PT=15.37;
      crystal="bct";
      CrystalStr_PT="Centered_Tetragonal";
      space_group="I4/mmm";
      space_group_number=139;
      Pearson_coefficient=0.0;
      lattice_constant[1]=392.5;lattice_constant[2]=392.5;lattice_constant[3]=323.8;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=2.00;
      radius_covalent_PT=200;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.2756;
      radii_Slatter=1.80;
      radii_Pyykko=1.69;
      electrical_conductivity=5.6E6;
      electronegativity_vec=1.50;
      hardness_Ghosh=1.9369;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=3.306;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=4000;
      melting_point=1572;
      vaporization_heat_PT=470;
      specific_heat_PT=99.1;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=47;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=3.25E-8;
      Volume_Magnetic_Susceptibility=0.0004995;
      Molar_Magnetic_Susceptibility=7.509E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Pa
    }
    // [AFLOW]STOP=Protoactinium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Uranium
    // Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium
    if(ZZ==92) { // Uranium
      Z=ZZ;
      symbol="U";
      name="Uranium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*238.03;
      MolarVolume=0.000012495;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      Density_PT=19.05;
      crystal="orc";
      CrystalStr_PT="Base_Orthorhombic";
      space_group="Cmcm";
      space_group_number=63;
      Pearson_coefficient=1.15611E-06;
      lattice_constant[1]=285.37;lattice_constant[2]=586.95;lattice_constant[3]=495.48;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=0.138;
      radius_PT=NNN;
      radius_covalent=1.96;
      radius_covalent_PT=196;
      radius_VanDerWaals_PT=186;
      radii_Ghosh08=1.9767;
      radii_Slatter=1.75;
      radii_Pyykko=1.70;
      electrical_conductivity=3.6E6;
      electronegativity_vec=1.38;
      hardness_Ghosh=2.2306;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=3.594;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=3927;
      melting_point=1135;
      vaporization_heat_PT=420;
      specific_heat_PT=116;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000139;
      thermal_conductivity=27;
      Brinelll_hardness=2400;
      Mohs_hardness=6;
      Vickers_hardness=1960;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=111;
      Young_modulus=208;
      bulk_modulus=100;
      Poisson_ratio_PT=0.23;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=2.16E-8;
      Volume_Magnetic_Susceptibility=0.000411;
      Molar_Magnetic_Susceptibility=5.14E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //U
    }
    // [AFLOW]STOP=Uranium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Neptunium
    // Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium
    if(ZZ==93) { // Neptunium
      Z=ZZ;
      symbol="Np";
      name="Neptunium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*237.05;
      MolarVolume=0.00001158924205379;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=7;
      valence_iupac=7;
      valence_PT=6;
      Density_PT=20.45;
      crystal="nnn";
      CrystalStr_PT="Simple_Orthorhombic";
      space_group="Pnma";
      space_group_number=62;
      Pearson_coefficient=0.0;
      lattice_constant[1]=666.3;lattice_constant[2]=472.3;lattice_constant[3]=488.7;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=1.90;
      radius_covalent_PT=190;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.7473;
      radii_Slatter=1.75;
      radii_Pyykko=1.71;
      electrical_conductivity=830000;
      electronegativity_vec=NNN;
      hardness_Ghosh=2.5241;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=3.882;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=4000;
      melting_point=644;
      vaporization_heat_PT=335;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=6;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Np
    }
    // [AFLOW]STOP=Neptunium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Plutonium
    // Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium
    if(ZZ==94) { // Plutonium
      Z=ZZ;
      symbol="Pu";
      name="Plutonium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*244.06;
      MolarVolume=0.00001231328219621;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=8;
      valence_iupac=7;
      valence_PT=6;
      Density_PT=19.816;
      crystal="nnn";
      CrystalStr_PT="Simple_Monoclinic";
      space_group="P12_1/m1";
      space_group_number=11;
      Pearson_coefficient=0.0;
      lattice_constant[1]=618.3;lattice_constant[2]=482.2;lattice_constant[3]=1096.3;
      lattice_angle[1]=PI/2;lattice_angle[2]=1.776571;lattice_angle[3]=PI/2;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=1.87;
      radius_covalent_PT=187;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.4496;
      radii_Slatter=1.75;
      radii_Pyykko=1.72;
      electrical_conductivity=670000;
      electronegativity_vec=NNN;
      hardness_Ghosh=3.0436;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.391;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=3230;
      melting_point=640;
      vaporization_heat_PT=325;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=6;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=43;
      Young_modulus=96;
      bulk_modulus=NNN;
      Poisson_ratio_PT=0.21;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=3.17E-8;
      Volume_Magnetic_Susceptibility=0.0006282;
      Molar_Magnetic_Susceptibility=7.735E-9;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Pu
    }
    // [AFLOW]STOP=Plutonium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Americium
    // Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium
    if(ZZ==95) { // Americium
      Z=ZZ;
      symbol="Am";
      name="Americium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*243.06;
      MolarVolume=0.00001777615215801;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=9;
      valence_iupac=7;
      valence_PT=4;
      Density_PT=13.67;
      crystal="nnn";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=346.81;lattice_constant[2]=346.81;lattice_constant[3]=1124.1;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=1.80;
      radius_covalent_PT=180;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.2915;
      radii_Slatter=1.75;
      radii_Pyykko=1.66;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=3.4169;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.678;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=2011;
      melting_point=1176;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=10;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="Paramagnetic";
      Mass_Magnetic_Susceptibility=5.15E-8;
      Volume_Magnetic_Susceptibility=0.000704;
      Molar_Magnetic_Susceptibility=1.251E-8;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Am
    }
    // [AFLOW]STOP=Americium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Curium
    // Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium
    if(ZZ==96) { // Curium
      Z=ZZ;
      symbol="Cm";
      name="Curium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*247.07;
      MolarVolume=0.00001828275351591;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=10;
      valence_iupac=8;
      valence_PT=4;
      Density_PT=13.51;
      crystal="nnn";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=349.6;lattice_constant[2]=349.6;lattice_constant[3]=1133.1;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=1.69;
      radius_covalent_PT=169;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.2960;
      radii_Slatter=NNN;
      radii_Pyykko=1.66;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=3.4050;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.745;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=3110;
      melting_point=1345;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=NNN;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Cm
    }
    // [AFLOW]STOP=Curium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Berkelium
    // Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium
    if(ZZ==97) { // Berkelium
      Z=ZZ;
      symbol="Bk";
      name="Berkelium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*247.07;
      MolarVolume=0.00001671177266576;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=11;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=14.78;
      crystal="nnn";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=341.6;lattice_constant[2]=341.6;lattice_constant[3]=1106.9;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.1247;
      radii_Slatter=NNN;
      radii_Pyykko=1.68;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=3.9244;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.256;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=NNN;
      melting_point=1050;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=10;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Bk
    }
    // [AFLOW]STOP=Berkelium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Californium
    // Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium
    if(ZZ==98) { // Californium
      Z=ZZ;
      symbol="Cf";
      name="Californium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*251.08;
      MolarVolume=0.00001662251655629;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=12;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=15.1;
      crystal="nnn";
      CrystalStr_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      Pearson_coefficient=0.0;
      lattice_constant[1]=338;lattice_constant[2]=338;lattice_constant[3]=1102.5;
      lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.0465;
      radii_Slatter=NNN;
      radii_Pyykko=1.68;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=4.2181;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.542;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=NNN;
      melting_point=900;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=NNN;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Cf
    }
    // [AFLOW]STOP=Californium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Einsteinium
    // Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium
    if(ZZ==99) { // Einsteinium
      Z=ZZ;
      symbol="Es";
      name="Einsteinium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*252.08;
      MolarVolume=NNN;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=13;
      valence_iupac=4;
      valence_PT=4;
      Density_PT=NNN;
      crystal="nnn";
      CrystalStr_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      Pearson_coefficient=0.0;
      lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
      lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.9785;
      radii_Slatter=NNN;
      radii_Pyykko=1.65;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=4.5116;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.830;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=NNN;
      melting_point=860;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=NNN;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Es
    }
    // [AFLOW]STOP=Einsteinium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Fermium
    // Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium
    if(ZZ==100) { // Fermium
      Z=ZZ;
      symbol="Fm";
      name="Fermium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*257.1;
      MolarVolume=NNN;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=14;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=NNN;
      crystal="nnn";
      CrystalStr_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      Pearson_coefficient=0.0;
      lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
      lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
      phase="nnn";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.9188;
      radii_Slatter=NNN;
      radii_Pyykko=1.67;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=4.8051;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.118;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=NNN;
      melting_point=1527;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=NNN;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Fm
    }
    // [AFLOW]STOP=Fermium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Mendelevium
    // Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium
    if(ZZ==101) { // Mendelevium
      Z=ZZ;
      symbol="Md";
      name="Mendelevium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*258.1;
      MolarVolume=NNN;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=15;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=NNN;
      crystal="nnn";
      CrystalStr_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      Pearson_coefficient=0.0;
      lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
      lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
      phase="nnn";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8659;
      radii_Slatter=NNN;
      radii_Pyykko=1.73;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=5.0990;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.406;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=NNN;
      melting_point=828;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=NNN;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Md
    }
    // [AFLOW]STOP=Mendelevium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Nobelium
    // Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium
    if(ZZ==102) { // Nobelium
      Z=ZZ;
      symbol="No";
      name="Nobelium";
      Period=7;
      Group=NNN;
      Series="Actinide";
      Block="f";
      mass=AMU2KILOGRAM*259.1;
      MolarVolume=NNN;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=16;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=NNN;
      crystal="nnn";
      CrystalStr_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      Pearson_coefficient=0.0;
      lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
      lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
      phase="nnn";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8188;
      radii_Slatter=NNN;
      radii_Pyykko=1.76;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=5.3926;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.694;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=NNN;
      melting_point=828;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=NNN;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //No
    }
    // [AFLOW]STOP=Nobelium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Lawrencium
    // Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium
    if(ZZ==103) { // Lawrencium
      Z=ZZ;
      symbol="Lr";
      name="Lawrencium";
      Period=7;
      Group=3;
      Series="TransitionMetal";
      Block="d";
      mass=AMU2KILOGRAM*262.11;
      MolarVolume=NNN;
      volume=NNN;
      Miedema_Vm=NNN;
      valence_std=17;
      valence_iupac=3;
      valence_PT=3;
      Density_PT=NNN;
      crystal="nnn";
      CrystalStr_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      Pearson_coefficient=0.0;
      lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
      lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
      phase="nnn";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8086;
      radii_Slatter=NNN;
      radii_Pyykko=1.61;
      electrical_conductivity=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=5.4607;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.760;
      electron_affinity_PT=NNN;
      Miedema_phi_star=NNN;
      Miedema_nws=NNN;
      Miedema_gamma_s=NNN;
      Pettifor_scale=0;
      boiling_point=NNN;
      melting_point=1627;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      thermal_conductivity=NNN;
      Brinelll_hardness=NNN;
      Mohs_hardness=NNN;
      Vickers_hardness=NNN;
      Hardness_Pearson=NNN;
      Hardness_Putz=NNN;
      Hardness_RB=NNN;
      shear_modulus=NNN;
      Young_modulus=NNN;
      bulk_modulus=NNN;
      Poisson_ratio_PT=NNN;
      Miedema_BVm=NNN;
      Magnetic_Type_PT="NNN";
      Mass_Magnetic_Susceptibility=NNN;
      Volume_Magnetic_Susceptibility=NNN;
      Molar_Magnetic_Susceptibility=NNN;
      Curie_point=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Lr
    }
    // [AFLOW]STOP=Lawrencium
    // ********************************************************************************************************************************************************
  }
} // namespace xelement

#endif // _AFLOW_XELEMENT_CPP

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2020              *
// *                                                                        *
// **************************************************************************
