// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2019
#ifndef _AFLOW_XELEMENT_H
#define _AFLOW_XELEMENT_H
//#include "aflow.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// _XELEMENT_PROTOTYPES_

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// XELEMENT
// look into aflow.h for the definitions

class xelement { // simple class.. nothing fancy
public:
  // constructor destructor                              // constructor/destructor
  xelement();                                            // default, just allocate
  xelement(uint);                                            // default, just allocate
  ~xelement();                                           // kill everything
  const xelement& operator=(const xelement &b);          // copy
  void clear();
  // content                                             // content
  bool verbose;
  // [AFLOW]START=DECLARATION
  string symbol;                          // http://periodictable.com      // DU 2019/05/17   // DONE SC(2019/05/24)
  string name;                            // http://periodictable.com      // DU 2019/05/17   // DONE SC(2019/05/24)
  double Period;                          // http://periodictable.com      // DU 2019/05/17
  double Group;                           // http://periodictable.com      // DU 2019/05/17
  string Series;                          // http://periodictable.com For Nh,Fl,Mc,Lv,Ts Value is a guess based on periodic table trend.      // DU 2019/05/17 
  string Block;                           // http://periodictable.com      // DU 2019/05/17
  //                                          
  double mass;                            // (kg)     // DONE SC(2019/05/24)
  double MolarVolume;                     // (m^3/mol) http://periodictable.com      // DU 2019/05/17
  double volume;                          // atomic volume in A^3 from the FCC vasp table and/or successive calculations // DONE SC(2019/05/24)
  double Miedema_Vm;                      // (V_m^{2/3} in (cm^2)) Miedema Rule Table 1a Physica 100B (1980) 1-28
  // for lanthines from J.A. Alonso and N.H. March. Electrons in Metals and Alloys, Academic Press, London (1989) (except La)
  double valence_std;                     // http://en.wikipedia.org/wiki/Valence_(chemistry) standard: number electrons minus closed shell at leff (noble gas)
  double valence_iupac;                   // http://en.wikipedia.org/wiki/Valence_(chemistry) IUPAC Maximum number of univalent atoms that may combine with an atom of the element under consideration, or with a fragment, or for which an atom of this element can be substituted.
  double valence_PT;                      //           http://periodictable.com      // DU 2019/05/17
  double Density_PT;                      // (g/cm^3)  http://periodictable.com      // DU 2019/05/17
  string crystal;                         // Ashcroft-Mermin                                                                                                                   
  string CrystalStr_PT;                   // http://periodictable.com      // DU 2019/05/17
  string space_group;                     // http://periodictable.com      // DU 2019/05/17
  uint space_group_number;                // http://periodictable.com      // DU 2019/05/17
  double Pearson_coefficient;             // Pearson mass deviation coefficient // ME181020
  xvector<double> lattice_constant;       // (pm) http://periodictable.com      // DU 2019/05/17
  xvector<double> lattice_angle;          // (rad) http://periodictable.com      // DU 2019/05/17
  string phase;                           //      http://periodictable.com      // DU 2019/05/17
  double radius;                          // Saxena (nm)
  double radius_PT;                       // (pm)       http://periodictable.com      // DU 2019/05/17
  double radius_covalent_PT;              // (pm)       http://periodictable.com      // DU 2019/05/17
  double radius_covalent;                 // (Angstrom) Dalton Trans. 2836, 2832-2838 (2008) // DX and CO - 9/4/17
  double radius_VanDerWaals_PT;           // (pm)       http://periodictable.com      // DU 2019/05/17
  double radii_Ghosh08;                    // (Angstrom) Journal of Molecular Structure: THEOCHEM 865, 60–67 (2008)      // DU 2019/05/17
  double radii_Slatter;                    // (Angstrom) J. of Chem. Phys. 41, 3199 (1964)      // DU 2019/05/17
  double radii_Pyykko;                     // (pm) single bond covalent radii  Chem. Eur. J. 15, 186-197 (2009)      // DU 2019/05/17
  //                                          
  double electrical_conductivity;          // (S/m)  http://periodictable.com  Value given for graphite. Diamond electrical conductivity is approximately 0.001.      // DU 2019/05/17
  double electronegativity_vec;           // Saxena
  double hardness_Ghosh;                   // (eV) Int. J. Quantum Chem 110, 1206-1213 (2010) Table III       // DU 2019/05/17
  double electronegativityPearson;                  // (eV) Inorg. Chem., 27(4), 734–740 (1988)      // DU 2019/05/17
  double electronegativityGhosh;                    // (eV) Journal of Theoretical and Computational Chemistry, 4, 21-33 (2005)      // DU 2019/05/17
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
  string Magnetic_Type_PT;                 //           http://periodictable.com  // DU 2019/05/17
  double Mass_Magnetic_Susceptibility;      // (m^3/K)   http://periodictable.com // DU 2019/05/17
  double Volume_Magnetic_Susceptibility;    //           http://periodictable.com // DU 2019/05/17
  double Molar_Magnetic_Susceptibility;     // (m^3/mol) http://periodictable.com // DU 2019/05/17
  double Curie_point;                     // (K)       http://periodictable.com   // DU 2019/05/17
  //
  double refractive_index;                 // http://periodictable.com C:diamond      // DU 2019/05/17
  string color_PT;                        // http://periodictable.com      // DU 2019/05/17
  //
  double HHIP;                            // Chem. Mater. 25(15), 2911–2920 (2013) Herfindahl–Hirschman Index (HHI), HHIP: for elemental production, Uncertinities in HHI_P: C,O,F,Cl,Sc,Ga,Rb,Ru,Rh,Cs,Hf,Os,Ir,Tl.      // DU 2019/05/17
  double HHIR;                            // Chem. Mater. 25(15), 2911–2920 (2013) Herfindahl–Hirschman Index (HHI), HHIR: for elemental reserves,   Uncertinities in HHI_R: Be,C,N,O,F,Na,Mg,Al,Si,S,Cl,Ca,Sc,Ga,Ge,As,Rb,Sr,Ru,Rh,Pd,In,Cs,Hf,Os,Ir,Pt,Tl.      // DU 2019/05/17
  double xray_scatt;                      // shift+1 // All data collected from the NIST online tables: http://physics.nist.gov/PhysRefData/FFast/html/form.html//

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
  xelement element_initialize(uint Z);                      // function to clean up the name

private:                                                    //
  void free();                                              // free space
};

void elements_initialize(void);

extern std::vector<xelement> velement;        // store starting from ONE

#endif

/*

subst "atom_symbol_vec" "vatom_symbol" `find . -name "aflow.h"`
subst "atom_name_vec" "vatom_name" `find . -name "aflow.h"`
subst "atom_mass_vec" "vatom_mass" `find . -name "aflow.h"`
subst "atom_volume_vec" "vatom_volume" `find . -name "aflow.h"`
subst "atom_valence_iupac_vec" "vatom_valence_iupac" `find . -name "aflow.h"`
subst "atom_valence_std_vec" "vatom_valence_std" `find . -name "aflow.h"`
subst "atom_radius_vec" "vatom_radius" `find . -name "aflow.h"`
subst "atom_radius_covalent_vec" "vatom_radius_covalent" `find . -name "aflow.h"`
subst "atom_electronegativity_vec" "vatom_electronegativity" `find . -name "aflow.h"`
subst "atom_crystal_vec" "vatom_crystal" `find . -name "aflow.h"`
subst "atom_miedema_phi_star" "vatom_miedema_phi_star" `find . -name "aflow.h"`
subst "atom_miedema_nws" "vatom_miedema_nws" `find . -name "aflow.h"`
subst "atom_miedema_Vm" "vatom_miedema_Vm" `find . -name "aflow.h"`
subst "atom_miedema_gamma_s" "vatom_miedema_gamma_s" `find . -name "aflow.h"`
subst "atom_miedema_BVm" "vatom_miedema_BVm" `find . -name "aflow.h"`
subst "xray_scatt_vec" "vatom_xray_scatt" `find . -name "aflow.h"`
subst "pettifor_scale" "vatom_pettifor_scale" `find . -name "aflow.h"`
subst "pearson_coefficient" "vatom_pearson_coefficient" `find . -name "aflow.h"`


*/
