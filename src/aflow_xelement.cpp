// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2019
#ifndef _AFLOW_XELEMENT_CPP
#define _AFLOW_XELEMENT_CPP
#include "aflow.h"

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
  string symbol_vec;                      // http://periodictable.com      //DU 5/17/19
  string name_vec;                        // http://periodictable.com      //DU 5/17/19
  double Period;                          // http://periodictable.com      //DU 5/17/19
  double Group;                           // http://periodictable.com      //DU 5/17/19
  string Series;                          // http://periodictable.com For Nh,Fl,Mc,Lv,Ts Value is a guess based on periodic table trend.      //DU 5/17/19 
  string Block;                           // http://periodictable.com      //DU 5/17/19
  //                                          
  double mass_vec;                        // (kg)
  double MolarVolume;                     // (m^3/mol) http://periodictable.com      //DU 5/17/19
  double volume_vec;                      // atomic volume in A^3 from the FCC vasp table and/or successive calculations
  double Miedema_Vm;                      // (V_m^{2/3} in (cm^2)) Miedema Rule Table 1a Physica 100B (1980) 1-28
  // for lanthines from J.A. Alonso and N.H. March. Electrons in Metals and Alloys, Academic Press, London (1989) (except La)
  double valence_std_vec;                 // http://en.wikipedia.org/wiki/Valence_(chemistry) standard: number electrons minus closed shell at leff (noble gas)
  double valence_iupac_vec;               // http://en.wikipedia.org/wiki/Valence_(chemistry) IUPAC Maximum number of univalent atoms that may combine with an atom of the element under consideration, or with a fragment, or for which an atom of this element can be substituted.
  double valence_PT;                      //           http://periodictable.com      //DU 5/17/19
  double Density_PT;                      // (g/cm^3)  http://periodictable.com      //DU 5/17/19
  string crystal_vec;                     // Ashcroft-Mermin                                                                                                                   
  string CrystalStr_PT;                   // http://periodictable.com      //DU 5/17/19
  string space_group;                     // http://periodictable.com      //DU 5/17/19
  uint space_group_number;                // http://periodictable.com      //DU 5/17/19
  double Pearson_coefficient;             // Pearson mass deviation coefficient // ME181020
  xvector<double> lattice_constant;       // (pm) http://periodictable.com      //DU 5/17/19
  xvector<double> lattice_angle;          // (rad) http://periodictable.com      //DU 5/17/19
  string phase;                           //      http://periodictable.com      //DU 5/17/19
  double radius_vec;                      // Saxena (nm)
  double radius_PT;                       // (pm)       http://periodictable.com      //DU 5/17/19
  double radius_covalent_PT;              // (pm)       http://periodictable.com      //DU 5/17/19
  double radius_covalent_vec;             // (Angstrom) Dalton Trans. 2836, 2832-2838 (2008) // DX and CO - 9/4/17
  double radius_VanDerWaals_PT;           // (pm)       http://periodictable.com      //DU 5/17/19
  double RadiiGhosh08;                    // (Angstrom) Journal of Molecular Structure: THEOCHEM 865, 60–67 (2008)      //DU 5/17/19
  double RadiiSlatter;                    // (Angstrom) J. of Chem. Phys. 41, 3199 (1964)      //DU 5/17/19
  double RadiiPyykko;                     // (pm) single bond covalent radii  Chem. Eur. J. 15, 186-197 (2009)      //DU 5/17/19
  //                                          
  double ElectricalConductivity;          // (S/m)  http://periodictable.com  Value given for graphite. Diamond electrical conductivity is approximately 0.001.      //DU 5/17/19
  double electronegativity_vec;           // Saxena
  double HardnessGhosh;                   // (eV) Int. J. Quantum Chem 110, 1206-1213 (2010) Table III       //DU 5/17/19
  double ElecNegPearson;                  // (eV) Inorg. Chem., 27(4), 734–740 (1988)      //DU 5/17/19
  double ElecNegGhosh;                    // (eV) Journal of Theoretical and Computational Chemistry, 4, 21-33 (2005)      //DU 5/17/19
  double ElectronAffinity_PT;             // (kJ/mol)  http://periodictable.com       //DU 5/17/19
  double Miedema_phi_star;                // (V)        (phi^\star   Miedema Rule Table 1a Physica 100B 1-28 (1980)
  double Miedema_nws;                     // (d.u.)^1/3 n_{ws}^{1/3} Miedema Rule Table 1a Physica 100B 1-28 (1980)
  double Miedema_gamma_s;                 // (mJ/m^2)   \gamma_s^0   Miedema Rule Table 1a Physica 100B 1-28 (1980)
  double Pettifor_scale;                  // Chemical Scale Pettifor Solid State Communications 51 31-34 (1984)
  //                                          
  double boiling_point;                   // (Celsius), http://periodictable.com C:diamond, P:"YELLOW" Phosphorus, As:sublimates at this T.      //DU 5/17/19
  double melting_point;                   // (Celsius), http://periodictable.com He does not solidify at standard pressure,C: Value given for diamond form, P : Value given for "YELLOW" phosphorus form, S : Value given for monoclinic, beta form, Se: Value given for hexagonal, gray form, Bk: Value given for alpha form.           //DU 5/17/19
  double VaporizationHeat_PT;             // (kJ/mol)   http://periodictable.com      //DU 5/17/19
  double SpecificHeat_PT;                 // (J/(kg.K)) http://periodictable.com Gas_Phase:H(H2),He,N(N2),O(O2),F(F2),Ne,Cl(Cl2),Ar,Kr,Tc,Xe,Rn,Ra,Pa -- Liquid_Phase:Br,Hg -- Solid Phase: B(rhombic),C(graphite),S(rhombic),P(phase of P.4),As(alpha),Se(hexagonal),Cd(gamma),Sn(gray),Li,In,Be,Na,Mg,Al,Si,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,Rb,Sr,Y,Zr,Nb,Mo,Ru,Rh,Pd,Ag,Sb,Te,I,Cs,Ba,La,Ce,Pr,Nd,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Tl,Pb,Bi,Ac,Th,U.      //DU 5/17/19 
  double critical_Pressure;                // (Atm)      http://periodictable.com Li,Na,K,Rb: Value estimated based on extrapolation.      //DU 5/17/19
  double critical_Temperature_PT;          // (K)        http://periodictable.com Li,Na,K,Rb: Value estimated based on extrapolation.      //DU 5/17/19
  double thermal_expansion;               // (K^{-1})   http://periodictable.com C:graphite      //DU 5/17/19
  double thermal_conductivity;            // (W/(mK))   http://periodictable.com      //DU 5/17/19
  //                                         
  double Brinelll_hardness;               // (MPa)  http://periodictable.com For Ge value is converted from Mohs scale      //DU 5/17/19
  double Mohs_hardness;                   //        http://periodictable.com For C, value given for graphite. Diamond value is 10.0; For Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Lu converted from Vickers scale.      //DU 5/17/19
  double Vickers_hardness;                // (MPa)  http://periodictable.com For Si,Ge,As,Ru,Os converted from Brinell scale.      //DU 5/17/19
  double Hardness_Pearson;                // (eV)   Inorg. Chem. 27(4) 734-740 (1988).      //DU 5/17/19
  double Hardness_Putz;                   // (eV/atom) International Journal of Quantum Chemistry, Vol 106, 361–389 (2006), TABLE-V.      //DU 5/17/19
  double Hardness_RB;                     // (eV)   Robles and Bartolotti, J. Am. Chem. Soc. 106, 3723-3727 (1984).      //DU 5/17/19
  double ShearModulus;                    // (GPa)  http://periodictable.com      //DU 5/17/19
  double YoungModulus;                    // (GPa)  http://periodictable.com      //DU 5/17/19
  double BulkModulus;                     // (GPa)  http://periodictable.com      //DU 5/17/19
  double PoissonRatio_PT;                 // (--)   http://periodictable.com      //DU 5/17/19
  double Miedema_BVm;                     // (kJ/mole) BV_m Miedema Rule Table 1a Physica 100B 1-28 (1980) 
  //
  string MagneticType_PT;                 //           http://periodictable.com       //DU 5/17/19
  double MassMagneticSusceptibility;      // (m^3/K)   http://periodictable.com       //DU 5/17/19
  double VolumeMagneticSusceptibility;    //           http://periodictable.com       //DU 5/17/19
  double MolarMagneticSusceptibility;     // (m^3/mol) http://periodictable.com       //DU 5/17/19
  double Curie_point;                     // (K)       http://periodictable.com       //DU 5/17/19
  //
  double RefractiveIndex;                 // http://periodictable.com C:diamond      //DU 5/17/19
  string color_PT;                        // http://periodictable.com      //DU 5/17/19
  //
  double HHIP;                            // Chem. Mater. 25(15), 2911–2920 (2013) Herfindahl–Hirschman Index (HHI), HHIP: for elemental production, Uncertinities in HHI_P: C,O,F,Cl,Sc,Ga,Rb,Ru,Rh,Cs,Hf,Os,Ir,Tl.      //DU 5/17/19
  double HHIR;                            // Chem. Mater. 25(15), 2911–2920 (2013) Herfindahl–Hirschman Index (HHI), HHIR: for elemental reserves,   Uncertinities in HHI_R: Be,C,N,O,F,Na,Mg,Al,Si,S,Cl,Ca,Sc,Ga,Ge,As,Rb,Sr,Ru,Rh,Pd,In,Cs,Hf,Os,Ir,Pt,Tl.      //DU 5/17/19
  double xray_scatt_vec;                  // shift+1 // All data collected from the NIST online tables: http://physics.nist.gov/PhysRefData/FFast/html/form.html//

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

// constructors
  xelement::xelement() {
    // DEFAULT
    verbose=FALSE;
    // [AFLOW]START=CONSTRUCTOR
    symbol_vec="UNDEFINED";
    name_vec="UNDEFINED";
    Period=NNN;
    Group=NNN; 
    Series="UNDEFINED";
    Block="nnn";      
    //                                          
    mass_vec=NNN; //  AMU2KILOGRAM goes inside.
    MolarVolume=NNN;  
    volume_vec=NNN;      
    Miedema_Vm=NNN;      
    //
    valence_std_vec=NNN;  
    valence_iupac_vec=NNN;
    valence_PT=NNN;       
    Density_PT=NNN;       
    crystal_vec="nnn";    
    CrystalStr_PT="UNDEFINED";
    space_group="nnn";     
    space_group_number=NNN;
    Pearson_coefficient=NNN;
    lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
    lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN; 
    phase="nnn";         
    radius_vec=NNN;         
    radius_PT=NNN;          
    radius_covalent_PT=NNN;   
    radius_covalent_vec=NNN;  
    radius_VanDerWaals_PT=NNN;
    RadiiGhosh08=NNN;         
    RadiiSlatter=NNN;         
    RadiiPyykko=NNN;          
    //                                          
    ElectricalConductivity=NNN;
    electronegativity_vec=NNN;    
    HardnessGhosh=NNN;            
    ElecNegPearson=NNN;           
    ElecNegGhosh=NNN;             
    ElectronAffinity_PT=NNN;      
    Miedema_phi_star=NNN;         
    Miedema_nws=NNN;              
    Miedema_gamma_s=NNN;          
    //
    Pettifor_scale=NNN;          
    //
    boiling_point=NNN;         
    melting_point=NNN;         
    VaporizationHeat_PT=NNN;     
    SpecificHeat_PT=NNN;         
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
    ShearModulus=NNN;    
    YoungModulus=NNN;    
    BulkModulus=NNN;     
    PoissonRatio_PT=NNN;    
    Miedema_BVm=NNN;        
    //
    MagneticType_PT="UNDEFINED";     
    MassMagneticSusceptibility=NNN;
    VolumeMagneticSusceptibility=NNN;
    MolarMagneticSusceptibility=NNN; 
    Curie_point=NNN;                  
    //
    RefractiveIndex=NNN;             
    color_PT="UNDEFINED";               
    //
    HHIP=NNN;                           
    HHIR=NNN;                           
    xray_scatt_vec=NNN;   
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

const xelement& xelement::operator=(const xelement& b) {       // operator=
  if(this != &b) {
    free();
    // will populate
    verbose=b.verbose;
    // [AFLOW]START=ASSIGNMENT
    symbol_vec=b.symbol_vec;
    name_vec=b.name_vec;
    Period=b.Period;
    Group=b.Group; 
    Series=b.Series;
    Block=b.Block;      
    //                                          
    mass_vec=b.mass_vec;
    MolarVolume=b.MolarVolume;  
    volume_vec=b.volume_vec;      
    Miedema_Vm=b.Miedema_Vm;      
    //
    valence_std_vec=b.valence_std_vec;  
    valence_iupac_vec=b.valence_iupac_vec;
    valence_PT=b.valence_PT;       
    Density_PT=b.Density_PT;       
    crystal_vec=b.crystal_vec;    
    CrystalStr_PT=b.CrystalStr_PT;
    space_group=b.space_group;
    space_group_number=b.space_group_number;    
    Pearson_coefficient=b.Pearson_coefficient;
    lattice_constant=b.lattice_constant; 
    lattice_angle=b.lattice_angle;   
    phase=b.phase;
    radius_vec=b.radius_vec;         
    radius_PT=b.radius_PT;          
    radius_covalent_PT=b.radius_covalent_PT;   
    radius_covalent_vec=b.radius_covalent_vec;  
    radius_VanDerWaals_PT=b.radius_VanDerWaals_PT;
    RadiiGhosh08=b.RadiiGhosh08;         
    RadiiSlatter=b.RadiiSlatter;         
    RadiiPyykko=b.RadiiPyykko;          
    //                                          
    ElectricalConductivity=b.ElectricalConductivity;
    electronegativity_vec=b.electronegativity_vec;    
    HardnessGhosh=b.HardnessGhosh;            
    ElecNegPearson=b.ElecNegPearson;           
    ElecNegGhosh=b.ElecNegGhosh;             
    ElectronAffinity_PT=b.ElectronAffinity_PT;      
    Miedema_phi_star=b.Miedema_phi_star;         
    Miedema_nws=b.Miedema_nws;              
    Miedema_gamma_s=b.Miedema_gamma_s;          
    //
    Pettifor_scale=b.Pettifor_scale;          
    //
    boiling_point=b.boiling_point;         
    melting_point=b.melting_point;         
    VaporizationHeat_PT=b.VaporizationHeat_PT;     
    SpecificHeat_PT=b.SpecificHeat_PT;         
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
    ShearModulus=b.ShearModulus;    
    YoungModulus=b.YoungModulus;    
    BulkModulus=b.BulkModulus;     
    PoissonRatio_PT=b.PoissonRatio_PT;    
    Miedema_BVm=b.Miedema_BVm;        
    //
    MagneticType_PT=b.MagneticType_PT;
    MassMagneticSusceptibility=b.MassMagneticSusceptibility;
    VolumeMagneticSusceptibility=b.VolumeMagneticSusceptibility;
    MolarMagneticSusceptibility=b.MolarMagneticSusceptibility; 
    Curie_point=b.Curie_point;                  
    //
    RefractiveIndex=b.RefractiveIndex;             
    color_PT=b.color_PT;         
    //
    HHIP=b.HHIP;                           
    HHIR=b.HHIR;                           
    xray_scatt_vec=b.xray_scatt_vec;    
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
 


// constructors
  xelement::xelement(uint Z) {
    // DEFAULT
    verbose=FALSE;

    // ROW 1
  // s-electron systems

// ********************************************************************************************************************************************************
 // [AFLOW]START=Hydrogen
 // Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen
 if(Z==1) {  // Hydrogen
    symbol_vec="H";
    name_vec="Hydrogen";
    Period=1;
    Group=1;
    Series="Nonmetal";
    Block="s";
    mass_vec=AMU2KILOGRAM*1.0079;
    MolarVolume=0.01121;
    volume_vec=0.75110;
    Miedema_Vm=NNN;
    valence_std_vec=1;
    valence_iupac_vec=1;
    valence_PT=1;
    Density_PT=0.899E-4;
    crystal_vec="hex";
    CrystalStr_PT="Simple_Hexagonal";
    space_group="P6_3/mmc";
    space_group_number=194;
    Pearson_coefficient=0.00011460743;
    lattice_constant[1]=470;lattice_constant[2]=470;lattice_constant[3]=340;
    lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
    phase="Gas";
    radius_vec=0.046;
    radius_PT=53;
    radius_covalent_vec=0.31;
    radius_covalent_PT=31;
    radius_VanDerWaals_PT=120;
    RadiiGhosh08=0.5292;
    RadiiSlatter=0.25;
    RadiiPyykko=0.32;
    ElectricalConductivity=NNN;
    electronegativity_vec=2.10;
    HardnessGhosh=6.4299;
    ElecNegPearson=7.18;
    ElecNegGhosh=7.178;
    ElectronAffinity_PT=72.8;
    Miedema_phi_star=5.2;
    Miedema_nws=1.5;
    Miedema_gamma_s=NNN;
    Pettifor_scale=0;
    boiling_point=-252.87;
    melting_point=-259.14;
    VaporizationHeat_PT=0.452;
    SpecificHeat_PT=14300;
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
    ShearModulus=NNN;
    YoungModulus=NNN;
    BulkModulus=NNN;
    PoissonRatio_PT=NNN;
    Miedema_BVm=NNN;
    MagneticType_PT="Diamagnetic";
    MassMagneticSusceptibility=-2.48E-8;
    VolumeMagneticSusceptibility=-2.23E-9;
    MolarMagneticSusceptibility=-4.999E-11;
    Curie_point=NNN;
    color_PT="COLORLESS";
    RefractiveIndex=1.000132;
    HHIP=NNN;
    HHIR=NNN;
    xray_scatt_vec=1.000;
    // H volume wrong *dimer* MIEDEMA =PAUL VAN DER PUT book
 }
 // [AFLOW]STOP=Hydrogen
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Helium
 // Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium
 if(Z==2) {  // Helium
     symbol_vec="He";
     name_vec="Helium";
     Period=1;
     Group=18;
     Series="NobleGas";
     Block="s";
     mass_vec=AMU2KILOGRAM*4.0026;
     MolarVolume=0.022424;
     volume_vec=-1.000;
     Miedema_Vm=NNN;
     valence_std_vec=0;
     valence_iupac_vec=0;
     valence_PT=0;
     Density_PT=1.785E-4;
     crystal_vec="hcp";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=8.32328E-8;
     lattice_constant[1]=424.2;lattice_constant[2]=424.2;lattice_constant[3]=424.2;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Gas";
     radius_vec=NNN;
     radius_PT=31;
     radius_covalent_vec=0.28;
     radius_covalent_PT=28;
     radius_VanDerWaals_PT=140;
     RadiiGhosh08=0.3113;
     RadiiSlatter=NNN;
     RadiiPyykko=0.46;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=12.5449;
     ElecNegPearson=NNN;
     ElecNegGhosh=12.046;
     ElectronAffinity_PT=0;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=-268.93;
     melting_point=NNN;
     VaporizationHeat_PT=0.083;
     SpecificHeat_PT=5193.1;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-5.9E-9;
     VolumeMagneticSusceptibility=-1.05E-9;
     MolarMagneticSusceptibility=-2.36E-11;
     Curie_point=NNN;
     color_PT="COLORLESS";
     RefractiveIndex=1.000035;
     HHIP=3200;
     HHIR=3900;
     xray_scatt_vec=2.000;
     // He
 }
 // [AFLOW]STOP=Helium
// ********************************************************************************************************************************************************

// ROW2
// s-electron systems
// ********************************************************************************************************************************************************
 // [AFLOW]START=Lithium
 // Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium
 if(Z==3) {  // Lithium
     symbol_vec="Li";
     name_vec="Lithium";
     Period=2;
     Group=1;
     Series="AlkaliMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*6.941;
     MolarVolume=0.00001297;
     volume_vec=20.24110;
     Miedema_Vm=5.5;
     valence_std_vec=1;
     valence_iupac_vec=1;
     valence_PT=1;
     Density_PT=0.535;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=0.0014588232;
     lattice_constant[1]=351;lattice_constant[2]=351;lattice_constant[3]=351;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.152;
     radius_PT=167;
     radius_covalent_vec=1.28;
     radius_covalent_PT=128;
     radius_VanDerWaals_PT=182;
     RadiiGhosh08=1.6283;
     RadiiSlatter=1.45;
     RadiiPyykko=1.33;
     ElectricalConductivity=1.1E7;
     electronegativity_vec=0.98;
     HardnessGhosh=2.3746;
     ElecNegPearson=3.01;
     ElecNegGhosh=2.860;
     ElectronAffinity_PT=59.6;
     Miedema_phi_star=2.85;
     Miedema_nws=0.98;
     Miedema_gamma_s=530;
     Pettifor_scale=0.45;
     boiling_point=1342;
     melting_point=180.54;
     VaporizationHeat_PT=147;
     SpecificHeat_PT=3570;
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
     ShearModulus=4.2;
     YoungModulus=4.9;
     BulkModulus=11;
     PoissonRatio_PT=NNN;
     Miedema_BVm=1.5;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=2.56E-8;
     VolumeMagneticSusceptibility=0.0000137;
     MolarMagneticSusceptibility=1.78E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=2900;
     HHIR=4200;
     xray_scatt_vec=3.00145;
     // Li
 }
 // [AFLOW]STOP=Lithium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Beryllium
 // Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium
 if(Z==4) {  // Beryllium
     symbol_vec="Be";
     name_vec="Beryllium";
     Period=2;
     Group=2;
     Series="AlkalineEarthMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*9.0122;
     MolarVolume=4.8767E-6;
     volume_vec=7.83290;
     Miedema_Vm=2.9;
     valence_std_vec=2;
     valence_iupac_vec=2;
     valence_PT=2;
     Density_PT=1.848;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=228.58;lattice_constant[2]=228.58;lattice_constant[3]=358.43;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.114;
     radius_PT=112;
     radius_covalent_vec=0.96;
     radius_covalent_PT=96;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.0855;
     RadiiSlatter=1.05;
     RadiiPyykko=1.02;
     ElectricalConductivity=2.5E7;
     electronegativity_vec=1.57;
     HardnessGhosh=3.4968;
     ElecNegPearson=4.90;
     ElecNegGhosh=3.945;
     ElectronAffinity_PT=0;
     Miedema_phi_star=4.20;
     Miedema_nws=1.60;
     Miedema_gamma_s=1900;
     Pettifor_scale=1.50;
     boiling_point=2470;
     melting_point=1287;
     VaporizationHeat_PT=297;
     SpecificHeat_PT=1820;
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
     ShearModulus=132;
     YoungModulus=287;
     BulkModulus=130;
     PoissonRatio_PT=0.032;
     Miedema_BVm=4.9;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-1.26E-8;
     VolumeMagneticSusceptibility=-0.00002328;
     MolarMagneticSusceptibility=-1.136E-10;
     Curie_point=NNN;
     color_PT="SLATEGRAY";
     RefractiveIndex=NNN;
     HHIP=8000;
     HHIR=4000;
     /*xray_scatt_vec=NNN;*/
     // Be
 }
 // [AFLOW]STOP=Beryllium
// ********************************************************************************************************************************************************

// p-electron systems
// ********************************************************************************************************************************************************
 // [AFLOW]START=Boron
 // Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron
 if(Z==5) {  // Boron
     symbol_vec="B";
     name_vec="Boron";
     Period=2;
     Group=13;
     Series="Metalloid";
     Block="p";
     mass_vec=AMU2KILOGRAM*10.81;
     MolarVolume=4.3943E-6;
     volume_vec=5.88420;
     Miedema_Vm=2.8;
     valence_std_vec=3;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=2.46;
     crystal_vec="tet";
     CrystalStr_PT="Simple_Trigonal";
     space_group="R_3m";
     space_group_number=166;
     Pearson_coefficient=0.00135391428;
     lattice_constant[1]=506;lattice_constant[2]=506;lattice_constant[3]=506;
     lattice_angle[1]=1.01334;lattice_angle[2]=1.01334;lattice_angle[3]=1.01334;
     phase="Solid";
     radius_vec=0.097;
     radius_PT=87;
     radius_covalent_vec=0.84;
     radius_covalent_PT=85;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.8141;
     RadiiSlatter=0.85;
     RadiiPyykko=0.85;
     ElectricalConductivity=0.0001;
     electronegativity_vec=2.04;
     HardnessGhosh=4.6190;
     ElecNegPearson=4.29;
     ElecNegGhosh=5.031;
     ElectronAffinity_PT=26.7;
     Miedema_phi_star=4.75;
     Miedema_nws=1.55;
     Miedema_gamma_s=NNN;
     Pettifor_scale=2.00;
     boiling_point=4000;
     melting_point=2075;
     VaporizationHeat_PT=507;
     SpecificHeat_PT=1030;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=320;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-8.7E-9;
     VolumeMagneticSusceptibility=-0.0000214;
     MolarMagneticSusceptibility=-9.41E-11;
     Curie_point=NNN;
     color_PT="BLACK";
     RefractiveIndex=NNN;
     HHIP=2900;
     HHIR=2000;
     /*xray_scatt_vec=NNN;*/
     // B
 }
 // [AFLOW]STOP=Boron
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Carbon
 // Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon
 if(Z==6) {  // Carbon
     symbol_vec="C";
     name_vec="Carbon";
     Period=2;
     Group=14;
     Series="Nonmetal";
     Block="p";
     mass_vec=AMU2KILOGRAM*12.011;
     MolarVolume=5.3146E-6;
     volume_vec=5.59490;
     Miedema_Vm=1.8;
     valence_std_vec=4;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=2.26;
     crystal_vec="dia";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.00007387218;
     lattice_constant[1]=246.4;lattice_constant[2]=246.4;lattice_constant[3]=671.1;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.077;
     radius_PT=67;
     radius_covalent_vec=0.76;
     radius_covalent_PT=76;
     radius_VanDerWaals_PT=170;
     RadiiGhosh08=0.6513;
     RadiiSlatter=0.70;
     RadiiPyykko=0.75;
     ElectricalConductivity=100000;
     electronegativity_vec=2.55;
     HardnessGhosh=5.7410;
     ElecNegPearson=6.27;
     ElecNegGhosh=6.116;
     ElectronAffinity_PT=153.9;
     Miedema_phi_star=6.20;
     Miedema_nws=1.90;
     Miedema_gamma_s=NNN;
     Pettifor_scale=2.50;
     boiling_point=4027;
     melting_point=3550;
     VaporizationHeat_PT=715;
     SpecificHeat_PT=710;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=33;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-6.2E-9;
     VolumeMagneticSusceptibility=-0.000014;
     MolarMagneticSusceptibility=-7.45E-11;
     Curie_point=NNN;
     color_PT="BLACK";
     RefractiveIndex=2.417;
     HHIP=500;
     HHIR=500;
     xray_scatt_vec=6.019;
     // C //DX and CO -9/4/17 atom_radius_covalent_vec_vec uses sp3 hybridization (most common)
 }
 // [AFLOW]STOP=Carbon
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Nitrogen
 // Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen
 if(Z==7) {  // Nitrogen
     symbol_vec="N";
     name_vec="Nitrogen";
     Period=2;
     Group=15;
     Series="Nonmetal";
     Block="p";
     mass_vec=AMU2KILOGRAM*14.0067;
     MolarVolume=0.011197;
     volume_vec=7.59940;
     Miedema_Vm=2.2;
     valence_std_vec=5;
     valence_iupac_vec=5;
     valence_PT=3;
     Density_PT=12.51E-4;
     crystal_vec="hex";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.00001857771;
     lattice_constant[1]=386.1;lattice_constant[2]=386.1;lattice_constant[3]=626.5;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Gas";
     radius_vec=0.071;
     radius_PT=56;
     radius_covalent_vec=0.71;
     radius_covalent_PT=71;
     radius_VanDerWaals_PT=155;
     RadiiGhosh08=0.5428;
     RadiiSlatter=0.65;
     RadiiPyykko=0.71;
     ElectricalConductivity=NNN;
     electronegativity_vec=3.04;
     HardnessGhosh=6.8625;
     ElecNegPearson=7.30;
     ElecNegGhosh=7.209;
     ElectronAffinity_PT=7;
     Miedema_phi_star=7.00;
     Miedema_nws=1.60;
     Miedema_gamma_s=NNN;
     Pettifor_scale=3.00;
     boiling_point=-195.79;
     melting_point=-210.1;
     VaporizationHeat_PT=2.79;
     SpecificHeat_PT=1040;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-5.4E-9;
     VolumeMagneticSusceptibility=-6.8E-9;
     MolarMagneticSusceptibility=-1.5E-10;
     Curie_point=NNN;
     color_PT="COLORLESS";
     RefractiveIndex=1.000298;
     HHIP=1300;
     HHIR=500;
     /*xray_scatt_vec=NNN;*/
     //N JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Nitrogen
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Oxygen
 // Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen
 if(Z==8) {  // Oxygen
     symbol_vec="O";
     name_vec="Oxygen";
     Period=2;
     Group=16;
     Series="Chalcogen";
     Block="p";
     mass_vec=AMU2KILOGRAM*15.9994;
     MolarVolume=0.011196;
     volume_vec=7.78230;
     Miedema_Vm=2.656;
     valence_std_vec=6;
     valence_iupac_vec=2;
     valence_PT=2;
     Density_PT=14.29E-4;
     crystal_vec="cub";
     CrystalStr_PT="Base-centered_Monoclinic";
     space_group="C12/m1";
     space_group_number=12;
     Pearson_coefficient=0.00003358805;
     lattice_constant[1]=540.3;lattice_constant[2]=342.9;lattice_constant[3]=508.6;
     lattice_angle[1]=PI/2;lattice_angle[2]=2.313085;lattice_angle[3]=PI/2;
     phase="Gas";
     radius_vec=0.060;
     radius_PT=48;
     radius_covalent_vec=0.66;
     radius_covalent_PT=66;
     radius_VanDerWaals_PT=152;
     RadiiGhosh08=0.4652;
     RadiiSlatter=0.60;
     RadiiPyykko=0.63;
     ElectricalConductivity=NNN;
     electronegativity_vec=3.44;
     HardnessGhosh=7.9854;
     ElecNegPearson=7.54;
     ElecNegGhosh=8.287;
     ElectronAffinity_PT=141;
     Miedema_phi_star=6.97;
     Miedema_nws=1.70;
     Miedema_gamma_s=NNN;
     Pettifor_scale=3.50;
     boiling_point=-182.9;
     melting_point=-218.3;
     VaporizationHeat_PT=3.41;
     SpecificHeat_PT=919;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.335E-6;
     VolumeMagneticSusceptibility=1.90772E-6;
     MolarMagneticSusceptibility=4.27184E-8;
     Curie_point=NNN;
     color_PT="COLORLESS";
     RefractiveIndex=1.000271;
     HHIP=500;
     HHIR=500;
     xray_scatt_vec=8.052;
     // O Table 27 of JUNKAI
 }
 // [AFLOW]STOP=Oxygen
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Fluorine
 // Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine
 if(Z==9) {  // Fluorine
     symbol_vec="F";
     name_vec="Fluorine";
     Period=2;
     Group=17;
     Series="Halogen";
     Block="p";
     mass_vec=AMU2KILOGRAM*18.9984;
     MolarVolume=0.011202;
     volume_vec=9.99090;
     Miedema_Vm=NNN;
     valence_std_vec=7;
     valence_iupac_vec=1;
     valence_PT=1;
     Density_PT=16.96E-4;
     crystal_vec="mcl";
     CrystalStr_PT="Base-centered_Monoclinic";
     space_group="C12/c1";
     space_group_number=15;
     Pearson_coefficient=0.0;
     lattice_constant[1]=550;lattice_constant[2]=328;lattice_constant[3]=728;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Gas";
     radius_vec=NNN;
     radius_PT=42;
     radius_covalent_vec=0.57;
     radius_covalent_PT=57;
     radius_VanDerWaals_PT=147;
     RadiiGhosh08=0.4071;
     RadiiSlatter=0.50;
     RadiiPyykko=0.64;
     ElectricalConductivity=NNN;
     electronegativity_vec=3.98;
     HardnessGhosh=9.1065;
     ElecNegPearson=10.41;
     ElecNegGhosh=9.372;
     ElectronAffinity_PT=328;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=4.00;
     boiling_point=-188.12;
     melting_point=-219.6;
     VaporizationHeat_PT=3.27;
     SpecificHeat_PT=824;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="COLORLESS";
     RefractiveIndex=1.000195;
     HHIP=1500;
     HHIR=1500;
     /*xray_scatt_vec=NNN;*/
     //F
 }
 // [AFLOW]STOP=Fluorine
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Neon
 // Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon
 if(Z==10) {  // Neon
     symbol_vec="Ne";
     name_vec="Neon";
     Period=2;
     Group=18;
     Series="NobleGas";
     Block="p";
     mass_vec=AMU2KILOGRAM*20.179;
     MolarVolume=0.02242;
     volume_vec=19.9052;
     Miedema_Vm=NNN;
     valence_std_vec=0;
     valence_iupac_vec=0;
     valence_PT=0;
     Density_PT=9E-4;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.00082783369;
     lattice_constant[1]=442.9;lattice_constant[2]=442.9;lattice_constant[3]=442.9;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Gas";
     radius_vec=0.160;
     radius_PT=38;
     radius_covalent_vec=0.58;
     radius_covalent_PT=58;
     radius_VanDerWaals_PT=154;
     RadiiGhosh08=0.3618;
     RadiiSlatter=NNN;
     RadiiPyykko=0.67;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=10.2303;
     ElecNegPearson=NNN;
     ElecNegGhosh=10.459;
     ElectronAffinity_PT=0;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=-246.08;
     melting_point=-248.59;
     VaporizationHeat_PT=1.75;
     SpecificHeat_PT=1030;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-4.1E-9;
     VolumeMagneticSusceptibility=-3.69E-9;
     MolarMagneticSusceptibility=-8.27E-11;
     Curie_point=NNN;
     color_PT="COLORLESS";
     RefractiveIndex=1.000067;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Ne volume calculated with fcc-pawpbe
 }
 // [AFLOW]STOP=Neon
// ********************************************************************************************************************************************************

// ROW3
// s-electron systems
// ********************************************************************************************************************************************************
 // [AFLOW]START=Sodium
 // Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium
 if(Z==11) {  // Sodium
     symbol_vec="Na";
     name_vec="Sodium";
     Period=3;
     Group=1;
     Series="AlkaliMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*22.9898;
     MolarVolume=0.00002375;
     volume_vec=36.9135;
     Miedema_Vm=8.3;
     valence_std_vec=1;
     valence_iupac_vec=1;
     valence_PT=1;
     Density_PT=0.968;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=0.0;
     lattice_constant[1]=429.06;lattice_constant[2]=429.06;lattice_constant[3]=429.06;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.186;
     radius_PT=190;
     radius_covalent_vec=1.66;
     radius_covalent_PT=166;
     radius_VanDerWaals_PT=227;
     RadiiGhosh08=2.1650;
     RadiiSlatter=1.80;
     RadiiPyykko=1.55;
     ElectricalConductivity=2.1E7;
     electronegativity_vec=0.93;
     HardnessGhosh=2.4441;
     ElecNegPearson=2.85;
     ElecNegGhosh=2.536;
     ElectronAffinity_PT=52.8;
     Miedema_phi_star=2.70;
     Miedema_nws=0.82;
     Miedema_gamma_s=260;
     Pettifor_scale=0.40;
     boiling_point=883;
     melting_point=97.72;
     VaporizationHeat_PT=97.7;
     SpecificHeat_PT=1230;
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
     ShearModulus=3.3;
     YoungModulus=10;
     BulkModulus=6.3;
     PoissonRatio_PT=NNN;
     Miedema_BVm=1.6;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=8.8E-9;
     VolumeMagneticSusceptibility=8.6E-6;
     MolarMagneticSusceptibility=2E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=1100;
     HHIR=500;
     /*xray_scatt_vec=NNN;*/
     // Na
 }
 // [AFLOW]STOP=Sodium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Magnesium
 // Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium
 if(Z==12) {  // Magnesium
     symbol_vec="Mg";
     name_vec="Magnesium";
     Period=3;
     Group=2;
     Series="AlkalineEarthMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*24.305;
     MolarVolume=0.000013984;
     volume_vec=22.8178;
     Miedema_Vm=5.8;
     valence_std_vec=2;
     valence_iupac_vec=2;
     valence_PT=2;
     Density_PT=1.738;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.00073988271;
     lattice_constant[1]=320.94;lattice_constant[2]=320.94;lattice_constant[3]=521.08;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.160;
     radius_PT=145;
     radius_covalent_vec=1.41;
     radius_covalent_PT=141;
     radius_VanDerWaals_PT=173;
     RadiiGhosh08=1.6711;
     RadiiSlatter=1.50;
     RadiiPyykko=1.39;
     ElectricalConductivity=2.3E7;
     electronegativity_vec=1.31;
     HardnessGhosh=3.0146;
     ElecNegPearson=3.75;
     ElecNegGhosh=3.310;
     ElectronAffinity_PT=0;
     Miedema_phi_star=3.45;
     Miedema_nws=1.17;
     Miedema_gamma_s=790;
     Pettifor_scale=1.28;
     boiling_point=1090;
     melting_point=650;
     VaporizationHeat_PT=128;
     SpecificHeat_PT=1020;
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
     ShearModulus=17;
     YoungModulus=45;
     BulkModulus=45;
     PoissonRatio_PT=0.29;
     Miedema_BVm=5.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=6.9E-9;
     VolumeMagneticSusceptibility=0.000012;
     MolarMagneticSusceptibility=1.68E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=5300;
     HHIR=500;
     /*xray_scatt_vec=NNN;*/
     //Mg
 }
 // [AFLOW]STOP=Magnesium
// ********************************************************************************************************************************************************

// p-electron systems
// ********************************************************************************************************************************************************
 // [AFLOW]START=Aluminium
 //Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium
 if(Z==13) {  // Aluminium
     symbol_vec="Al";
     name_vec="Aluminium";
     Period=3;
     Group=13;
     Series="PoorMetal";
     Block="p";
     mass_vec=AMU2KILOGRAM*26.9815;
     MolarVolume=9.99E-6;
     volume_vec=16.4000;
     Miedema_Vm=4.6;
     valence_std_vec=3;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=2.7;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.0;
     lattice_constant[1]=404.95;lattice_constant[2]=404.95;lattice_constant[3]=404.95;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.143;
     radius_PT=118;
     radius_covalent_vec=1.21;
     radius_covalent_PT=121;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.3608;
     RadiiSlatter=1.25;
     RadiiPyykko=1.26;
     ElectricalConductivity=3.8E7;
     electronegativity_vec=1.61;
     HardnessGhosh=3.5849;
     ElecNegPearson=3.23;
     ElecNegGhosh=4.084;
     ElectronAffinity_PT=42.5;
     Miedema_phi_star=4.20;
     Miedema_nws=1.39;
     Miedema_gamma_s=1200;
     Pettifor_scale=1.66;
     boiling_point=2519;
     melting_point=660.32;
     VaporizationHeat_PT=293;
     SpecificHeat_PT=904;
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
     ShearModulus=26;
     YoungModulus=70;
     BulkModulus=76;
     PoissonRatio_PT=0.35;
     Miedema_BVm=7.2;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=7.8E-9;
     VolumeMagneticSusceptibility=0.0000211;
     MolarMagneticSusceptibility=2.1E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=1600;
     HHIR=1000;
     /*xray_scatt_vec=NNN;*/
     //Al
 }
 // [AFLOW]STOP=Aluminium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Silicon
 // Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon
 if(Z==14) {  // Silicon
     symbol_vec="Si";
     name_vec="Silicon";
     Period=3;
     Group=14;
     Series="Metalloid";
     Block="p";
     mass_vec=AMU2KILOGRAM*28.0855;
     MolarVolume=0.000012054;
     volume_vec=14.3536;
     Miedema_Vm=4.2;
     valence_std_vec=4;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=2.33;
     crystal_vec="dia";
     CrystalStr_PT="Tetrahedral_Packing";
     space_group="Fd_3m";
     space_group_number=227;
     Pearson_coefficient=0.00020046752;
     lattice_constant[1]=543.09;lattice_constant[2]=543.09;lattice_constant[3]=543.09;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.117;
     radius_PT=111;
     radius_covalent_vec=1.11;
     radius_covalent_PT=111;
     radius_VanDerWaals_PT=210;
     RadiiGhosh08=1.1477;
     RadiiSlatter=1.10;
     RadiiPyykko=1.16;
     ElectricalConductivity=1000;
     electronegativity_vec=1.90;
     HardnessGhosh=4.1551;
     ElecNegPearson=4.77;
     ElecNegGhosh=4.857;
     ElectronAffinity_PT=133.6;
     Miedema_phi_star=4.70;
     Miedema_nws=1.50;
     Miedema_gamma_s=1290;
     Pettifor_scale=1.92;
     boiling_point=2900;
     melting_point=1414;
     VaporizationHeat_PT=359;
     SpecificHeat_PT=710;
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
     ShearModulus=NNN;
     YoungModulus=47;
     BulkModulus=100;
     PoissonRatio_PT=NNN;
     Miedema_BVm=11.9;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-1.6E-9;
     VolumeMagneticSusceptibility=-3.73E-6;
     MolarMagneticSusceptibility=-4.49E-11;
     Curie_point=NNN;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=4700;
     HHIR=1000;
     xray_scatt_vec=14.43;
     //Si ???
 }
 // [AFLOW]STOP=Silicon
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Phosphorus
 // Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus
 if(Z==15) {  // Phosphorus
     symbol_vec="P";
     name_vec="Phosphorus";
     Period=3;
     Group=15;
     Series="Nonmetal";
     Block="p";
     mass_vec=AMU2KILOGRAM*30.9738;
     MolarVolume=0.000016991;
     volume_vec=14.1995;
     Miedema_Vm=NNN;
     valence_std_vec=5;
     valence_iupac_vec=5;
     valence_PT=5;
     Density_PT=1.823;
     crystal_vec="cub";
     CrystalStr_PT="Simple_Triclinic";
     space_group="P-1";
     space_group_number=2;
     Pearson_coefficient=0.0;
     lattice_constant[1]=1145;lattice_constant[2]=550.3;lattice_constant[3]=1126.1;
     lattice_angle[1]=1.25384;lattice_angle[2]=1.57725;lattice_angle[3]=1.24896;
     phase="Solid";
     radius_vec=0.109;
     radius_PT=98;
     radius_covalent_vec=1.07;
     radius_covalent_PT=107;
     radius_VanDerWaals_PT=180;
     RadiiGhosh08=0.9922;
     RadiiSlatter=1.00;
     RadiiPyykko=1.11;
     ElectricalConductivity=1E7;
     electronegativity_vec=2.19;
     HardnessGhosh=4.7258;
     ElecNegPearson=5.62;
     ElecNegGhosh=5.631;
     ElectronAffinity_PT=71;
     Miedema_phi_star=5.5;
     Miedema_nws=1.65;
     Miedema_gamma_s=NNN;
     Pettifor_scale=2.18;
     boiling_point=280.5;
     melting_point=44.2;
     VaporizationHeat_PT=12.4;
     SpecificHeat_PT=769.7;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=11;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-1.13E-8;
     VolumeMagneticSusceptibility=-0.0000206;
     MolarMagneticSusceptibility=-3.5E-10;
     Curie_point=NNN;
     color_PT="COLORLESS";
     RefractiveIndex=1.001212;
     HHIP=2000;
     HHIR=5100;
     xray_scatt_vec=15.3133;
     //P MIEDEMA =PAUL VAN DER PUT book
 }
 // [AFLOW]STOP=Phosphorus
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Sulphur
 // Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur
 if(Z==16) {  // Sulphur
     symbol_vec="S";
     name_vec="Sulphur";
     Period=3;
     Group=16;
     Series="Chalcogen";
     Block="p";
     mass_vec=AMU2KILOGRAM*32.06;
     MolarVolume=0.000016357;
     volume_vec=15.7301;
     Miedema_Vm=4.376;
     valence_std_vec=6;
     valence_iupac_vec=6;
     valence_PT=6;
     Density_PT=1.96;
     crystal_vec="orc";
     CrystalStr_PT="Face-centered_Orthorhombic";
     space_group="Fddd";
     space_group_number=70;
     Pearson_coefficient=0.00016807795;
     lattice_constant[1]=1043.7;lattice_constant[2]=1284.5;lattice_constant[3]=2436.9;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.106;
     radius_PT=87;
     radius_covalent_vec=1.05;
     radius_covalent_PT=105;
     radius_VanDerWaals_PT=180;
     RadiiGhosh08=0.8739;
     RadiiSlatter=1.00;
     RadiiPyykko=1.03;
     ElectricalConductivity=1E-15;
     electronegativity_vec=2.58;
     HardnessGhosh=5.2960;
     ElecNegPearson=6.22;
     ElecNegGhosh=6.420;
     ElectronAffinity_PT=200;
     Miedema_phi_star=5.6;
     Miedema_nws=1.46;
     Miedema_gamma_s=NNN;
     Pettifor_scale=2.44;
     boiling_point=444.72;
     melting_point=115.21;
     VaporizationHeat_PT=9.8;
     SpecificHeat_PT=705;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=7.7;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-6.2E-9;
     VolumeMagneticSusceptibility=-0.0000122;
     MolarMagneticSusceptibility=-1.99E-10;
     Curie_point=NNN;
     color_PT="YELLOW";
     RefractiveIndex=1.001111;
     HHIP=700;
     HHIR=1000;
     /*xray_scatt_vec=NNN;*/
     //S Table 27 of JUNKAI
 }
 // [AFLOW]STOP=Sulphur
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Chlorine
 // Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine
 if(Z==17) {  // Chlorine
     symbol_vec="Cl";
     name_vec="Chlorine";
     Period=3;
     Group=17;
     Series="Halogen";
     Block="p";
     mass_vec=AMU2KILOGRAM*35.453;
     MolarVolume=0.01103;
     volume_vec=21.2947;
     Miedema_Vm=6.71;
     valence_std_vec=7;
     valence_iupac_vec=7;
     valence_PT=5;
     Density_PT=32.14E-4;
     crystal_vec="orc";
     CrystalStr_PT="Base_Orthorhombic";
     space_group="Cmca";
     space_group_number=64;
     Pearson_coefficient=0.00058238731;
     lattice_constant[1]=622.35;lattice_constant[2]=445.61;lattice_constant[3]=817.85;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Gas";
     radius_vec=0.107;
     radius_PT=79;
     radius_covalent_vec=1.02;
     radius_covalent_PT=102;
     radius_VanDerWaals_PT=175;
     RadiiGhosh08=0.7808;
     RadiiSlatter=1.00;
     RadiiPyykko=0.99;
     ElectricalConductivity=0.01;
     electronegativity_vec=3.16;
     HardnessGhosh=5.8662;
     ElecNegPearson=8.30;
     ElecNegGhosh=7.178;
     ElectronAffinity_PT=349;
     Miedema_phi_star=5.32;
     Miedema_nws=0.34;
     Miedema_gamma_s=1013;
     Pettifor_scale=2.70;
     boiling_point=-34.04;
     melting_point=-101.5;
     VaporizationHeat_PT=10.2;
     SpecificHeat_PT=478.2;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=1.1;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-7.2E-9;
     VolumeMagneticSusceptibility=-2.31E-8;
     MolarMagneticSusceptibility=-5.11E-10;
     Curie_point=NNN;
     color_PT="YELLOW";
     RefractiveIndex=1.000773;
     HHIP=1500;
     HHIR=1500;
     /*xray_scatt_vec=NNN;*/
     //Cl interpolation phi_star, nws, Vm, gamma JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Chlorine
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Argon
 //Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon
 if(Z==18) {  // Argon
     symbol_vec="Ar";
     name_vec="Argon";
     Period=3;
     Group=18;
     Series="NobleGas";
     Block="p";
     mass_vec=AMU2KILOGRAM*39.948;
     MolarVolume=0.022392;
     volume_vec=22.000;
     Miedema_Vm=NNN;
     valence_std_vec=0;
     valence_iupac_vec=2;
     valence_PT=0;
     Density_PT=17.84E-4;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.00003509919;
     lattice_constant[1]=525.6;lattice_constant[2]=525.6;lattice_constant[3]=525.6;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Gas";
     radius_vec=0.192;
     radius_PT=71;
     radius_covalent_vec=1.06;
     radius_covalent_PT=106;
     radius_VanDerWaals_PT=188;
     RadiiGhosh08=0.7056;
     RadiiSlatter=NNN;
     RadiiPyykko=0.96;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=6.4366;
     ElecNegPearson=NNN;
     ElecNegGhosh=7.951;
     ElectronAffinity_PT=0;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=-185.8;
     melting_point=-189.3;
     VaporizationHeat_PT=6.5;
     SpecificHeat_PT=520.33;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-6E-9;
     VolumeMagneticSusceptibility=-1.07E-8;
     MolarMagneticSusceptibility=-2.4E-10;
     Curie_point=NNN;
     color_PT="COLORLESS";
     RefractiveIndex=1.000281;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Ar guessed volume, must double check from results JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Argon
// ********************************************************************************************************************************************************

// ROW4
// s-electron systems
// ********************************************************************************************************************************************************
 // [AFLOW]START=Potassium
 // Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium
 if(Z==19) {  // Potassium
     symbol_vec="K";
     name_vec="Potassium";
     Period=4;
     Group=1;
     Series="AlkaliMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*39.0983;
     MolarVolume=0.00004568;
     volume_vec=73.9091;
     Miedema_Vm=12.8;
     valence_std_vec=1;
     valence_iupac_vec=1;
     valence_PT=1;
     Density_PT=0.856;
     crystal_vec="fcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=0.000164;
     lattice_constant[1]=532.8;lattice_constant[2]=532.8;lattice_constant[3]=532.8;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.231;
     radius_PT=243;
     radius_covalent_vec=2.03;
     radius_covalent_PT=203;
     radius_VanDerWaals_PT=275;
     RadiiGhosh08=3.2930;
     RadiiSlatter=2.20;
     RadiiPyykko=1.96;
     ElectricalConductivity=1.4E7;
     electronegativity_vec=0.82;
     HardnessGhosh=2.3273;
     ElecNegPearson=2.42;
     ElecNegGhosh=2.672;
     ElectronAffinity_PT=48.4;
     Miedema_phi_star=2.25;
     Miedema_nws=0.65;
     Miedema_gamma_s=150;
     Pettifor_scale=0.35;
     boiling_point=759;
     melting_point=63.38;
     VaporizationHeat_PT=76.9;
     SpecificHeat_PT=757;
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
     ShearModulus=1.3;
     YoungModulus=NNN;
     BulkModulus=3.1;
     PoissonRatio_PT=NNN;
     Miedema_BVm=1.5;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=6.7E-9;
     VolumeMagneticSusceptibility=5.74E-6;
     MolarMagneticSusceptibility=2.62E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=1700;
     HHIR=7200;
     /*xray_scatt_vec=NNN;*/
     //K
 }
 // [AFLOW]STOP=Potassium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Calcium
 // Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium
 if(Z==20) {  // Calcium
     symbol_vec="Ca";
     name_vec="Calcium";
     Period=4;
     Group=2;
     Series="AlkalineEarthMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*40.08;
     MolarVolume=0.000025857;
     volume_vec=42.1927;
     Miedema_Vm=8.8;
     valence_std_vec=2;
     valence_iupac_vec=2;
     valence_PT=2;
     Density_PT=1.55;
     crystal_vec="bcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.000297564;
     lattice_constant[1]=558.84;lattice_constant[2]=558.84;lattice_constant[3]=558.84;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.197;
     radius_PT=194;
     radius_covalent_vec=1.76;
     radius_covalent_PT=176;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.5419;
     RadiiSlatter=1.80;
     RadiiPyykko=1.71;
     ElectricalConductivity=2.9E7;
     electronegativity_vec=1.00;
     HardnessGhosh=2.7587;
     ElecNegPearson=2.2;
     ElecNegGhosh=3.140;
     ElectronAffinity_PT=2.37;
     Miedema_phi_star=2.55;
     Miedema_nws=0.91;
     Miedema_gamma_s=490;
     Pettifor_scale=0.60;
     boiling_point=1484;
     melting_point=842;
     VaporizationHeat_PT=155;
     SpecificHeat_PT=631;
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
     ShearModulus=7.4;
     YoungModulus=20;
     BulkModulus=17;
     PoissonRatio_PT=0.31;
     Miedema_BVm=4.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.38E-8;
     VolumeMagneticSusceptibility=0.00002139;
     MolarMagneticSusceptibility=5.531E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=3900;
     HHIR=1500;
     /*xray_scatt_vec=NNN;*/
     //Ca
 }
 // [AFLOW]STOP=Calcium
// ********************************************************************************************************************************************************

// d-electron systems: transition metals
// ********************************************************************************************************************************************************
 // [AFLOW]START=Scandium
 // Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium
 if(Z==21) {  // Scandium
     symbol_vec="Sc";
     name_vec="Scandium";
     Period=4;
     Group=3;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*44.9559;
     MolarVolume=0.000015061;
     volume_vec=24.6739;
     Miedema_Vm=6.1;
     valence_std_vec=3;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=2.985;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=330.9;lattice_constant[2]=330.9;lattice_constant[3]=527.33;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.160;
     radius_PT=184;
     radius_covalent_vec=1.70;
     radius_covalent_PT=170;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.4149;
     RadiiSlatter=1.60;
     RadiiPyykko=1.48;
     ElectricalConductivity=1.8E6;
     electronegativity_vec=1.36;
     HardnessGhosh=2.8582;
     ElecNegPearson=3.34;
     ElecNegGhosh=3.248;
     ElectronAffinity_PT=18.1;
     Miedema_phi_star=3.25;
     Miedema_nws=1.27;
     Miedema_gamma_s=1200;
     Pettifor_scale=0.74;
     boiling_point=2830;
     melting_point=1541;
     VaporizationHeat_PT=318;
     SpecificHeat_PT=567;
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
     ShearModulus=29;
     YoungModulus=74;
     BulkModulus=57;
     PoissonRatio_PT=0.28;
     Miedema_BVm=6.6;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=8.8E-8;
     VolumeMagneticSusceptibility=0.0002627;
     MolarMagneticSusceptibility=3.956E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=5500;
     HHIR=4500;
     xray_scatt_vec=21.34;
     //Sc
 }
 // [AFLOW]STOP=Scandium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Titanium
 // Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium
 if(Z==22) {  // Titanium
     symbol_vec="Ti";
     name_vec="Titanium";
     Period=4;
     Group=4;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*47.9;
     MolarVolume=0.000010621;
     volume_vec=17.1035;
     Miedema_Vm=4.8;
     valence_std_vec=4;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=4.507;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.000286456;
     lattice_constant[1]=295.08;lattice_constant[2]=295.08;lattice_constant[3]=468.55;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.147;
     radius_PT=176;
     radius_covalent_vec=1.60;
     radius_covalent_PT=160;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.2998;
     RadiiSlatter=1.40;
     RadiiPyykko=1.36;
     ElectricalConductivity=2.5E6;
     electronegativity_vec=1.54;
     HardnessGhosh=2.9578;
     ElecNegPearson=3.45;
     ElecNegGhosh=3.357;
     ElectronAffinity_PT=7.6;
     Miedema_phi_star=3.65;
     Miedema_nws=1.47;
     Miedema_gamma_s=2050;
     Pettifor_scale=0.79;
     boiling_point=3287;
     melting_point=1668;
     VaporizationHeat_PT=425;
     SpecificHeat_PT=520;
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
     ShearModulus=44;
     YoungModulus=116;
     BulkModulus=110;
     PoissonRatio_PT=0.32;
     Miedema_BVm=11.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=4.01E-8;
     VolumeMagneticSusceptibility=0.0001807;
     MolarMagneticSusceptibility=1.919E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=1100;
     HHIR=1600;
     xray_scatt_vec=22.24;
     //Ti
 }
 // [AFLOW]STOP=Titanium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Vanadium
 // Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium
 if(Z==23) {  // Vanadium
     symbol_vec="V";
     name_vec="Vanadium";
     Period=4;
     Group=5;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*50.9415;
     MolarVolume=8.3374E-6;
     volume_vec=13.2086;
     Miedema_Vm=4.1;
     valence_std_vec=5;
     valence_iupac_vec=5;
     valence_PT=5;
     Density_PT=6.11;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=9.54831E-07;
     lattice_constant[1]=303;lattice_constant[2]=303;lattice_constant[3]=303;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.132;
     radius_PT=171;
     radius_covalent_vec=1.53;
     radius_covalent_PT=153;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.1953;
     RadiiSlatter=1.35;
     RadiiPyykko=1.34;
     ElectricalConductivity=5E6;
     electronegativity_vec=1.63;
     HardnessGhosh=3.0573;
     ElecNegPearson=3.6;
     ElecNegGhosh=3.465;
     ElectronAffinity_PT=50.6;
     Miedema_phi_star=4.25;
     Miedema_nws=1.64;
     Miedema_gamma_s=2600;
     Pettifor_scale=0.84;
     boiling_point=3407;
     melting_point=1910;
     VaporizationHeat_PT=453;
     SpecificHeat_PT=489;
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
     ShearModulus=47;
     YoungModulus=128;
     BulkModulus=160;
     PoissonRatio_PT=0.37;
     Miedema_BVm=14.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=6.28E-8;
     VolumeMagneticSusceptibility=0.0003837;
     MolarMagneticSusceptibility=3.199E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=3300;
     HHIR=3400;
     /*xray_scatt_vec=NNN;*/
     //V
 }
 // [AFLOW]STOP=Vanadium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Chromium
 // Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium
 if(Z==24) {  // Chromium
     symbol_vec="Cr";
     name_vec="Chromium";
     Period=4;
     Group=6;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*51.996;
     MolarVolume=7.2317E-6;
     volume_vec=11.4136;
     Miedema_Vm=3.7;
     valence_std_vec=6;
     valence_iupac_vec=6;
     valence_PT=6;
     Density_PT=7.19;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=0.00013287;
     lattice_constant[1]=291;lattice_constant[2]=291;lattice_constant[3]=291;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.125;
     radius_PT=166;
     radius_covalent_vec=1.39;
     radius_covalent_PT=139;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.1000;
     RadiiSlatter=1.40;
     RadiiPyykko=1.22;
     ElectricalConductivity=7.9E6;
     electronegativity_vec=1.66;
     HardnessGhosh=3.1567;
     ElecNegPearson=3.72;
     ElecNegGhosh=3.573;
     ElectronAffinity_PT=64.3;
     Miedema_phi_star=4.65;
     Miedema_nws=1.74;
     Miedema_gamma_s=2400;
     Pettifor_scale=0.89;
     boiling_point=2671;
     melting_point=1907;
     VaporizationHeat_PT=339;
     SpecificHeat_PT=448;
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
     ShearModulus=115;
     YoungModulus=279;
     BulkModulus=160;
     PoissonRatio_PT=0.21;
     Miedema_BVm=14.0;
     MagneticType_PT="Antiferromagnetic";
     MassMagneticSusceptibility=4.45E-8;
     VolumeMagneticSusceptibility=0.0003177;
     MolarMagneticSusceptibility=2.314E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=3100;
     HHIR=4100;
     xray_scatt_vec=23.84;
     //Cr
 }
 // [AFLOW]STOP=Chromium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Manganese
 // Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese
 if(Z==25) {  // Manganese
     symbol_vec="Mn";
     name_vec="Manganese";
     Period=4;
     Group=7;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*54.93805;
     MolarVolume=7.3545E-6;
     volume_vec=10.6487;
     Miedema_Vm=3.8;
     valence_std_vec=7;
     valence_iupac_vec=7;
     valence_PT=4;
     Density_PT=7.47;
     crystal_vec="cub";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="I_43m";
     space_group_number=217;
     Pearson_coefficient=1.67276E-32;
     lattice_constant[1]=891.25;lattice_constant[2]=891.25;lattice_constant[3]=891.25;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.112;
     radius_PT=161;
     radius_covalent_vec=1.61;
     radius_covalent_PT=139;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.0124;
     RadiiSlatter=1.40;
     RadiiPyykko=1.19;
     ElectricalConductivity=620000;
     electronegativity_vec=1.55;
     HardnessGhosh=3.2564;
     ElecNegPearson=3.72;
     ElecNegGhosh=3.681;
     ElectronAffinity_PT=0;
     Miedema_phi_star=4.45;
     Miedema_nws=1.61;
     Miedema_gamma_s=1600;
     Pettifor_scale=0.94;
     boiling_point=2061;
     melting_point=1246;
     VaporizationHeat_PT=220;
     SpecificHeat_PT=479;
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
     ShearModulus=NNN;
     YoungModulus=198;
     BulkModulus=120;
     PoissonRatio_PT=NNN;
     Miedema_BVm=4.4;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.21E-7;
     VolumeMagneticSusceptibility=0.00090387;
     MolarMagneticSusceptibility=6.6475E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=1600;
     HHIR=1800;
     xray_scatt_vec=24.46;
     //xray_scatt_vec=24.3589; Mn JUNKAI CHANGED VALENCE // DX and CO- 9/4/17 atom_radius_covalent_vec_vec[i] uses high spin configuration (most frequent)
 }
 // [AFLOW]STOP=Manganese
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Iron
 // Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron
 if(Z==26) {  // Iron
     symbol_vec="Fe";
     name_vec="Iron";
     Period=4;
     Group=8;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*55.847;
     MolarVolume=7.0923E-6;
     volume_vec=10.2315;
     Miedema_Vm=3.7;
     valence_std_vec=8;
     valence_iupac_vec=6;
     valence_PT=3;
     Density_PT=7.874;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=9.17912E-05;
     lattice_constant[1]=286.65;lattice_constant[2]=286.65;lattice_constant[3]=286.65;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.124;
     radius_PT=156;
     radius_covalent_vec=1.52;
     radius_covalent_PT=132;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.9319;
     RadiiSlatter=1.40;
     RadiiPyykko=1.16;
     ElectricalConductivity=1E7;
     electronegativity_vec=1.83;
     HardnessGhosh=3.3559;
     ElecNegPearson=4.06;
     ElecNegGhosh=3.789;
     ElectronAffinity_PT=15.7;
     Miedema_phi_star=4.93;
     Miedema_nws=1.77;
     Miedema_gamma_s=2550;
     Pettifor_scale=0.99;
     boiling_point=2861;
     melting_point=1538;
     VaporizationHeat_PT=347;
     SpecificHeat_PT=449;
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
     ShearModulus=82;
     YoungModulus=211;
     BulkModulus=170;
     PoissonRatio_PT=0.29;
     Miedema_BVm=12.0;
     MagneticType_PT="Ferromagnetic";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=1043;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=2400;
     HHIR=1400;
     xray_scatt_vec=24.85;
     //xray_scatt_vec=24.6830; Fe JUNKAI CHANGED VALENCE // DX and CO - 9/4/17 atom_radius_covalent_vec_vec[i] uses high spin configuration (most frequent)
 }
 // [AFLOW]STOP=Iron
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Cobalt
 // Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt
 if(Z==27) {  // Cobalt
     symbol_vec="Co";
     name_vec="Cobalt";
     Period=4;
     Group=9;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*58.9332;
     MolarVolume=6.62E-6;
     volume_vec=10.3205;
     Miedema_Vm=3.5;
     valence_std_vec=9;
     valence_iupac_vec=5;
     valence_PT=4;
     Density_PT=8.9;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=250.71;lattice_constant[2]=250.71;lattice_constant[3]=406.95;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.125;
     radius_PT=152;
     radius_covalent_vec=1.26;
     radius_covalent_PT=126;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.8575;
     RadiiSlatter=1.35;
     RadiiPyykko=1.11;
     ElectricalConductivity=1.7E7;
     electronegativity_vec=1.88;
     HardnessGhosh=3.4556;
     ElecNegPearson=4.3;
     ElecNegGhosh=3.897;
     ElectronAffinity_PT=63.7;
     Miedema_phi_star=5.10;
     Miedema_nws=1.75;
     Miedema_gamma_s=2550;
     Pettifor_scale=1.04;
     boiling_point=2927;
     melting_point=1495;
     VaporizationHeat_PT=375;
     SpecificHeat_PT=421;
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
     ShearModulus=76;
     YoungModulus=209;
     BulkModulus=180;
     PoissonRatio_PT=0.31;
     Miedema_BVm=13.0;
     MagneticType_PT="Ferromagnetic";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=1394;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=3100;
     HHIR=2700;
     xray_scatt_vec=24.59;
     //Co JUNKAI CHANGED VALENCE // DX and CO - 9/4/17 atom_radius_covalent_vec_vec[i] uses low spin configuration (most frequent)
 }
 // [AFLOW]STOP=Cobalt
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Nickel
 // Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel
 if(Z==28) {  // Nickel
     symbol_vec="Ni";
     name_vec="Nickel";
     Period=4;
     Group=10;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*58.69;
     MolarVolume=6.5888E-6;
     volume_vec=10.8664;
     Miedema_Vm=3.5;
     valence_std_vec=10;
     valence_iupac_vec=4;
     valence_PT=2;
     Density_PT=8.908;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.000430773;
     lattice_constant[1]=352.4;lattice_constant[2]=352.4;lattice_constant[3]=352.4;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.125;
     radius_PT=149;
     radius_covalent_vec=1.24;
     radius_covalent_PT=124;
     radius_VanDerWaals_PT=163;
     RadiiGhosh08=1.7888;
     RadiiSlatter=1.35;
     RadiiPyykko=1.10;
     ElectricalConductivity=1.4E7;
     electronegativity_vec=1.91;
     HardnessGhosh=3.5550;
     ElecNegPearson=4.40;
     ElecNegGhosh=4.005;
     ElectronAffinity_PT=112;
     Miedema_phi_star=5.20;
     Miedema_nws=1.75;
     Miedema_gamma_s=2450;
     Pettifor_scale=1.09;
     boiling_point=2913;
     melting_point=1455;
     VaporizationHeat_PT=378;
     SpecificHeat_PT=445;
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
     ShearModulus=76;
     YoungModulus=200;
     BulkModulus=180;
     PoissonRatio_PT=0.31;
     Miedema_BVm=12.0;
     MagneticType_PT="Ferromagnetic";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=631;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=1000;
     HHIR=1500;
     xray_scatt_vec=25.02;
     //Ni
 }
 // [AFLOW]STOP=Nickel
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Copper
 // Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper
 if(Z==29) {  // Copper
     symbol_vec="Cu";
     name_vec="Copper";
     Period=4;
     Group=11;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*63.546;
     MolarVolume=7.0922E-6;
     volume_vec=12.0159;
     Miedema_Vm=3.7;
     valence_std_vec=11;
     valence_iupac_vec=4;
     valence_PT=2;
     Density_PT=8.96;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.00021086;
     lattice_constant[1]=361.49;lattice_constant[2]=361.49;lattice_constant[3]=361.49;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.128;
     radius_PT=145;
     radius_covalent_vec=1.32;
     radius_covalent_PT=132;
     radius_VanDerWaals_PT=140;
     RadiiGhosh08=1.725;
     RadiiSlatter=1.35;
     RadiiPyykko=1.12;
     ElectricalConductivity=5.9E7;
     electronegativity_vec=1.90;
     HardnessGhosh=3.6544;
     ElecNegPearson=4.48;
     ElecNegGhosh=4.113;
     ElectronAffinity_PT=118.4;
     Miedema_phi_star=4.55;
     Miedema_nws=1.47;
     Miedema_gamma_s=1850;
     Pettifor_scale=1.20;
     boiling_point=2562;
     melting_point=1084.62;
     VaporizationHeat_PT=300;
     SpecificHeat_PT=384.4;
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
     ShearModulus=48;
     YoungModulus=130;
     BulkModulus=140;
     PoissonRatio_PT=0.34;
     Miedema_BVm=9.3;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-1.08E-9;
     VolumeMagneticSusceptibility=-9.63E-6;
     MolarMagneticSusceptibility=-6.86E-11;
     Curie_point=NNN;
     color_PT="COPPER";
     RefractiveIndex=NNN;
     HHIP=1600;
     HHIR=1500;
     xray_scatt_vec=27.03;
     //Cu JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Copper
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Zinc
 // Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc
 if(Z==30) {  // Zinc
     symbol_vec="Zn";
     name_vec="Zinc";
     Period=4;
     Group=12;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*65.38;
     MolarVolume=9.157E-6;
     volume_vec=15.0827;
     Miedema_Vm=4.4;
     valence_std_vec=12;
     valence_iupac_vec=2;
     valence_PT=2;
     Density_PT=7.14;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.000595597;
     lattice_constant[1]=266.49;lattice_constant[2]=266.49;lattice_constant[3]=494.68;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.133;
     radius_PT=142;
     radius_covalent_vec=1.22;
     radius_covalent_PT=122;
     radius_VanDerWaals_PT=139;
     RadiiGhosh08=1.6654;
     RadiiSlatter=1.35;
     RadiiPyykko=1.18;
     ElectricalConductivity=1.7E7;
     electronegativity_vec=1.65;
     HardnessGhosh=3.7542;
     ElecNegPearson=4.45;
     ElecNegGhosh=4.222;
     ElectronAffinity_PT=0;
     Miedema_phi_star=4.10;
     Miedema_nws=1.32;
     Miedema_gamma_s=1020;
     Pettifor_scale=1.44;
     boiling_point=907;
     melting_point=419.53;
     VaporizationHeat_PT=119;
     SpecificHeat_PT=388;
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
     ShearModulus=43;
     YoungModulus=108;
     BulkModulus=70;
     PoissonRatio_PT=0.25;
     Miedema_BVm=5.5;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-2.21E-9;
     VolumeMagneticSusceptibility=-0.0000158;
     MolarMagneticSusceptibility=-1.45E-10;
     Curie_point=NNN;
     color_PT="SLATEGRAY";
     RefractiveIndex=1.00205;
     HHIP=1600;
     HHIR=1900;
     xray_scatt_vec=28.44;
     //Zn
 }
 // [AFLOW]STOP=Zinc
// ********************************************************************************************************************************************************

// p-electron systems 
// ********************************************************************************************************************************************************
 // [AFLOW]START=Gallium
 // Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium
 if(Z==31) {  // Gallium
     symbol_vec="Ga";
     name_vec="Gallium";
     Period=4;
     Group=13;
     Series="PoorMetal";
     Block="p";
     mass_vec=AMU2KILOGRAM*69.737;
     MolarVolume=0.000011809;
     volume_vec=18.9039;
     Miedema_Vm=5.2;
     valence_std_vec=3;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=5.904;
     crystal_vec="orc";
     CrystalStr_PT="Base_Orthorhombic";
     space_group="Cmca";
     space_group_number=64;
     Pearson_coefficient=0.000197588;
     lattice_constant[1]=451.97;lattice_constant[2]=766.33;lattice_constant[3]=452.6;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.135;
     radius_PT=136;
     radius_covalent_vec=1.22;
     radius_covalent_PT=122;
     radius_VanDerWaals_PT=187;
     RadiiGhosh08=1.4489;
     RadiiSlatter=1.30;
     RadiiPyykko=1.24;
     ElectricalConductivity=7.1E6;
     electronegativity_vec=1.81;
     HardnessGhosh=4.1855;
     ElecNegPearson=3.2;
     ElecNegGhosh=4.690;
     ElectronAffinity_PT=28.9;
     Miedema_phi_star=4.10;
     Miedema_nws=1.31;
     Miedema_gamma_s=830;
     Pettifor_scale=1.68;
     boiling_point=2204;
     melting_point=29.76;
     VaporizationHeat_PT=256;
     SpecificHeat_PT=371;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=6.7;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-3E-9;
     VolumeMagneticSusceptibility=-0.0000177;
     MolarMagneticSusceptibility=-2.09E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=5500;
     HHIR=1900;
     /*xray_scatt_vec=NNN;*/
     //Ga
 }
 // [AFLOW]STOP=Gallium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Germanium
 // Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium
 if(Z==32) {  // Germanium
     symbol_vec="Ge";
     name_vec="Germanium";
     Period=4;
     Group=14;
     Series="Metalloid";
     Block="p";
     mass_vec=AMU2KILOGRAM*72.59;
     MolarVolume=0.000013645;
     volume_vec=19.2948;
     Miedema_Vm=4.6;
     valence_std_vec=4;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=5.323;
     crystal_vec="dia";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.00058782;
     lattice_constant[1]=565.75;lattice_constant[2]=565.75;lattice_constant[3]=565.75;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.122;
     radius_PT=125;
     radius_covalent_vec=1.20;
     radius_covalent_PT=120;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.2823;
     RadiiSlatter=1.25;
     RadiiPyykko=1.21;
     ElectricalConductivity=2000;
     electronegativity_vec=2.01;
     HardnessGhosh=4.6166;
     ElecNegPearson=4.6;
     ElecNegGhosh=5.159;
     ElectronAffinity_PT=119;
     Miedema_phi_star=4.55;
     Miedema_nws=1.37;
     Miedema_gamma_s=1030;
     Pettifor_scale=1.92;
     boiling_point=2820;
     melting_point=938.3;
     VaporizationHeat_PT=334;
     SpecificHeat_PT=321.4;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=10.5;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-1.5E-9;
     VolumeMagneticSusceptibility=-7.98E-6;
     MolarMagneticSusceptibility=-1.09E-10;
     Curie_point=NNN;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=5300;
     HHIR=1900;
     /*xray_scatt_vec=NNN;*/
     //Ge
 }
 // [AFLOW]STOP=Germanium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Arsenic
 // Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic
 if(Z==33) {  // Arsenic
     symbol_vec="As";
     name_vec="Arsenic";
     Period=4;
     Group=15;
     Series="Metalloid";
     Block="p";
     mass_vec=AMU2KILOGRAM*74.9216;
     MolarVolume=0.000013082;
     volume_vec=19.0677;
     Miedema_Vm=5.2;
     valence_std_vec=5;
     valence_iupac_vec=5;
     valence_PT=5;
     Density_PT=5.727;
     crystal_vec="rhl";
     CrystalStr_PT="Simple_Trigonal";
     space_group="R_3m";
     space_group_number=166;
     Pearson_coefficient=0.0;
     lattice_constant[1]=375.98;lattice_constant[2]=375.98;lattice_constant[3]=1054.75;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.125;
     radius_PT=114;
     radius_covalent_vec=1.19;
     radius_covalent_PT=119;
     radius_VanDerWaals_PT=185;
     RadiiGhosh08=1.1450;
     RadiiSlatter=1.15;
     RadiiPyykko=1.21;
     ElectricalConductivity=3.3E6;
     electronegativity_vec=2.18;
     HardnessGhosh=5.0662;
     ElecNegPearson=5.3;
     ElecNegGhosh=5.628;
     ElectronAffinity_PT=78;
     Miedema_phi_star=4.80;
     Miedema_nws=1.44;
     Miedema_gamma_s=1000;
     Pettifor_scale=2.16;
     boiling_point=614;
     melting_point=817;
     VaporizationHeat_PT=32.4;
     SpecificHeat_PT=328;
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
     ShearModulus=NNN;
     YoungModulus=8;
     BulkModulus=22;
     PoissonRatio_PT=NNN;
     Miedema_BVm=5.1;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-3.9E-9;
     VolumeMagneticSusceptibility=-0.0000223;
     MolarMagneticSusceptibility=-2.92E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=1.001552;
     HHIP=3300;
     HHIR=4000;
     /*xray_scatt_vec=NNN;*/
     //As
 }
 // [AFLOW]STOP=Arsenic
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Selenium
 // Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium
 if(Z==34) {  // Selenium
     symbol_vec="Se";
     name_vec="Selenium";
     Period=4;
     Group=16;
     Series="Chalcogen";
     Block="p";
     mass_vec=AMU2KILOGRAM*78.96;
     MolarVolume=0.000016387;
     volume_vec=20.3733;
     Miedema_Vm=5.172;
     valence_std_vec=6;
     valence_iupac_vec=6;
     valence_PT=6;
     Density_PT=4.819;
     crystal_vec="hex";
     CrystalStr_PT="Simple_Monoclinic";
     space_group="P12_1/c1";
     space_group_number=14;
     Pearson_coefficient=0.00046279;
     lattice_constant[1]=905.4;lattice_constant[2]=908.3;lattice_constant[3]=1160.1;
     lattice_angle[1]=PI/2;lattice_angle[2]=1.58493;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.116;
     radius_PT=103;
     radius_covalent_vec=1.20;
     radius_covalent_PT=120;
     radius_VanDerWaals_PT=190;
     RadiiGhosh08=1.0424;
     RadiiSlatter=1.15;
     RadiiPyykko=1.16;
     ElectricalConductivity=NNN;
     electronegativity_vec=2.55;
     HardnessGhosh=5.4795;
     ElecNegPearson=5.89;
     ElecNegGhosh=6.096;
     ElectronAffinity_PT=195;
     Miedema_phi_star=5.17;
     Miedema_nws=1.40;
     Miedema_gamma_s=NNN;
     Pettifor_scale=2.40;
     boiling_point=685;
     melting_point=221;
     VaporizationHeat_PT=26;
     SpecificHeat_PT=321.2;
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
     ShearModulus=3.7;
     YoungModulus=10;
     BulkModulus=8.3;
     PoissonRatio_PT=0.33;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-4E-9;
     VolumeMagneticSusceptibility=-0.0000193;
     MolarMagneticSusceptibility=-3.16E-10;
     Curie_point=NNN;
     color_PT="GRAY";
     RefractiveIndex=1.000895;
     HHIP=2200;
     HHIR=1900;
     /*xray_scatt_vec=NNN;*/
     //Se Table 27 of JUNKAI
 }
 // [AFLOW]STOP=Selenium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Bromine
 // Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine
 if(Z==35) {  // Bromine
     symbol_vec="Br";
     name_vec="Bromine";
     Period=4;
     Group=17;
     Series="Halogen";
     Block="p";
     mass_vec=AMU2KILOGRAM*79.904;
     MolarVolume=0.00002561;
     volume_vec=26.3292;
     Miedema_Vm=7.31;
     valence_std_vec=7;
     valence_iupac_vec=7;
     valence_PT=5;
     Density_PT=3.12;
     crystal_vec="orc";
     CrystalStr_PT="Base_Orthorhombic";
     space_group="Cmca";
     space_group_number=64;
     Pearson_coefficient=0.000156277;
     lattice_constant[1]=672.65;lattice_constant[2]=464.51;lattice_constant[3]=870.23;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Liquid";
     radius_vec=0.119;
     radius_PT=94;
     radius_covalent_vec=1.20;
     radius_covalent_PT=120;
     radius_VanDerWaals_PT=185;
     RadiiGhosh08=0.9532;
     RadiiSlatter=1.15;
     RadiiPyykko=1.14;
     ElectricalConductivity=1E-10;
     electronegativity_vec=2.96;
     HardnessGhosh=5.9111;
     ElecNegPearson=7.59;
     ElecNegGhosh=6.565;
     ElectronAffinity_PT=324.6;
     Miedema_phi_star=5.20;
     Miedema_nws=1.35;
     Miedema_gamma_s=943;
     Pettifor_scale=2.64;
     boiling_point=59;
     melting_point=-7.3;
     VaporizationHeat_PT=14.8;
     SpecificHeat_PT=947.3;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=1.9;
     PoissonRatio_PT=NNN;
     Miedema_BVm=3.4;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-4.9E-9;
     VolumeMagneticSusceptibility=-0.0000153;
     MolarMagneticSusceptibility=-7.83E-10;
     Curie_point=NNN;
     color_PT="RED";
     RefractiveIndex=1.001132;
     HHIP=3300;
     HHIR=6900;
     /* xray_scatt_vec=NNN;*/
     //Br interpolation phi_star, nws, Vm, gamma, BVm JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Bromine
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Krypton
 // Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton
 if(Z==36) {  // Krypton
     symbol_vec="Kr";
     name_vec="Krypton";
     Period=4;
     Group=18;
     Series="NobleGas";
     Block="p";
     mass_vec=AMU2KILOGRAM*83.8;
     MolarVolume=0.02235;
     volume_vec=-1.0000;
     Miedema_Vm=NNN;
     valence_std_vec=0;
     valence_iupac_vec=2;
     valence_PT=2;
     Density_PT=37.5E-4;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.000248482;
     lattice_constant[1]=570.6;lattice_constant[2]=570.6;lattice_constant[3]=570.6;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Gas";
     radius_vec=0.197;
     radius_PT=87;
     radius_covalent_vec=1.16;
     radius_covalent_PT=116;
     radius_VanDerWaals_PT=202;
     RadiiGhosh08=0.8782;
     RadiiSlatter=NNN;
     RadiiPyykko=1.17;
     ElectricalConductivity=NNN;
     electronegativity_vec=3;
     HardnessGhosh=6.3418;
     ElecNegPearson=NNN;
     ElecNegGhosh=7.033;
     ElectronAffinity_PT=0;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=-153.22;
     melting_point=-157.36;
     VaporizationHeat_PT=9.02;
     SpecificHeat_PT=248.05;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-4.4E-9;
     VolumeMagneticSusceptibility=-1.65E-8;
     MolarMagneticSusceptibility=-3.69E-10;
     Curie_point=NNN;
     color_PT="COLORLESS";
     RefractiveIndex=1.000427;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Kr
 }
 // [AFLOW]STOP=Krypton
// ********************************************************************************************************************************************************

// ROW5
// s-electron systems
// ********************************************************************************************************************************************************
 // [AFLOW]START=Rubidium
 // Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium
 if(Z==37) {  // Rubidium
     symbol_vec="Rb";
     name_vec="Rubidium";
     Period=5;
     Group=1;
     Series="AlkaliMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*85.4678;
     MolarVolume=0.000055788;
     volume_vec=91.2738;
     Miedema_Vm=14.6;
     valence_std_vec=1;
     valence_iupac_vec=1;
     valence_PT=1;
     Density_PT=1.532;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=0.000109697;
     lattice_constant[1]=558.5;lattice_constant[2]=558.5;lattice_constant[3]=558.5;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.251;
     radius_PT=265;
     radius_covalent_vec=2.20;
     radius_covalent_PT=220;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=3.8487;
     RadiiSlatter=2.35;
     RadiiPyykko=2.10;
     ElectricalConductivity=8.3E6;
     electronegativity_vec=0.82;
     HardnessGhosh=2.1204;
     ElecNegPearson=2.34;
     ElecNegGhosh=2.849;
     ElectronAffinity_PT=46.9;
     Miedema_phi_star=2.10;
     Miedema_nws=0.60;
     Miedema_gamma_s=120;
     Pettifor_scale=0.30;
     boiling_point=688;
     melting_point=39.31;
     VaporizationHeat_PT=71;
     SpecificHeat_PT=364;
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
     ShearModulus=NNN;
     YoungModulus=2.4;
     BulkModulus=2.5;
     PoissonRatio_PT=NNN;
     Miedema_BVm=1.8;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=2.6E-9;
     VolumeMagneticSusceptibility=3.98E-6;
     MolarMagneticSusceptibility=2.22E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=6000;
     HHIR=6000;
     /*xray_scatt_vec=NNN;*/
     //Rb
 }
 // [AFLOW]STOP=Rubidium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Strontium
 // Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium
 if(Z==38) {  // Strontium
     symbol_vec="Sr";
     name_vec="Strontium";
     Period=5;
     Group=2;
     Series="AlkalineEarthMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*87.62;
     MolarVolume=0.000033316;
     volume_vec=55.4105;
     Miedema_Vm=10.2;
     valence_std_vec=2;
     valence_iupac_vec=2;
     valence_PT=2;
     Density_PT=2.63;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=6.09969E-05;
     lattice_constant[1]=608.49;lattice_constant[2]=608.49;lattice_constant[3]=608.49;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.215;
     radius_PT=219;
     radius_covalent_vec=1.95;
     radius_covalent_PT=195;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.9709;
     RadiiSlatter=2.00;
     RadiiPyykko=1.85;
     ElectricalConductivity=7.7E6;
     electronegativity_vec=0.95;
     HardnessGhosh=2.5374;
     ElecNegPearson=2.0;
     ElecNegGhosh=3.225;
     ElectronAffinity_PT=5.03;
     Miedema_phi_star=2.40;
     Miedema_nws=0.84;
     Miedema_gamma_s=430;
     Pettifor_scale=0.55;
     boiling_point=1382;
     melting_point=777;
     VaporizationHeat_PT=137;
     SpecificHeat_PT=300;
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
     ShearModulus=6.1;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=0.28;
     Miedema_BVm=3.9;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.32E-9;
     VolumeMagneticSusceptibility=3.47E-6;
     MolarMagneticSusceptibility=1.16E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=4200;
     HHIR=3000;
     /*xray_scatt_vec=NNN;*/
     //Sr
 }
 // [AFLOW]STOP=Strontium
// ********************************************************************************************************************************************************

// d-electron systems: transition metals
// ********************************************************************************************************************************************************
 // [AFLOW]START=Yttrium
 // Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium
 if(Z==39) {  // Yttrium
     symbol_vec="Y";
     name_vec="Yttrium";
     Period=5;
     Group=3;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*88.9059;
     MolarVolume=0.000019881;
     volume_vec=32.4546;
     Miedema_Vm=7.3;
     valence_std_vec=3;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=4.472;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=364.74;lattice_constant[2]=364.74;lattice_constant[3]=573.06;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.181;
     radius_PT=212;
     radius_covalent_vec=1.90;
     radius_covalent_PT=190;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.8224;
     RadiiSlatter=1.80;
     RadiiPyykko=1.63;
     ElectricalConductivity=1.8E6;
     electronegativity_vec=1.22;
     HardnessGhosh=2.6335;
     ElecNegPearson=3.19;
     ElecNegGhosh=3.311;
     ElectronAffinity_PT=29.6;
     Miedema_phi_star=3.20;
     Miedema_nws=1.21;
     Miedema_gamma_s=1100;
     Pettifor_scale=0.70;
     boiling_point=3345;
     melting_point=1526;
     VaporizationHeat_PT=380;
     SpecificHeat_PT=298;
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
     ShearModulus=26;
     YoungModulus=64;
     BulkModulus=41;
     PoissonRatio_PT=0.24;
     Miedema_BVm=7.2;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=6.66E-8;
     VolumeMagneticSusceptibility=0.0002978;
     MolarMagneticSusceptibility=5.921E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9800;
     HHIR=2600;
     /*xray_scatt_vec=NNN;*/
     //Y
 }
 // [AFLOW]STOP=Yttrium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Zirconium
 // Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium
 if(Z==40) {  // Zirconium
     symbol_vec="Zr";
     name_vec="Zirconium";
     Period=5;
     Group=4;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*91.22;
     MolarVolume=0.000014011;
     volume_vec=23.2561;
     Miedema_Vm=5.8;
     valence_std_vec=4;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=6.511;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.000342629;
     lattice_constant[1]=323.2;lattice_constant[2]=323.2;lattice_constant[3]=514.7;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.158;
     radius_PT=206;
     radius_covalent_vec=1.75;
     radius_covalent_PT=175;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.688;
     RadiiSlatter=1.55;
     RadiiPyykko=1.54;
     ElectricalConductivity=2.4E6;
     electronegativity_vec=1.33;
     HardnessGhosh=2.7298;
     ElecNegPearson=3.64;
     ElecNegGhosh=3.398;
     ElectronAffinity_PT=41.1;
     Miedema_phi_star=3.40;
     Miedema_nws=1.39;
     Miedema_gamma_s=1950;
     Pettifor_scale=0.76;
     boiling_point=4409;
     melting_point=1855;
     VaporizationHeat_PT=580;
     SpecificHeat_PT=278;
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
     ShearModulus=33;
     YoungModulus=67;
     BulkModulus=NNN;
     PoissonRatio_PT=0.34;
     Miedema_BVm=12.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.68E-8;
     VolumeMagneticSusceptibility=0.000109;
     MolarMagneticSusceptibility=1.53E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=3400;
     HHIR=2600;
     /*xray_scatt_vec=NNN;*/
     //Zr
 }
 // [AFLOW]STOP=Zirconium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Niobium
 // Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium
 if(Z==41) {  // Niobium
     symbol_vec="Nb";
     name_vec="Niobium";
     Period=5;
     Group=5;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*92.9064;
     MolarVolume=0.000010841;
     volume_vec=18.3132;
     Miedema_Vm=4.9;
     valence_std_vec=5;
     valence_iupac_vec=5;
     valence_PT=5;
     Density_PT=8.57;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=0.0;
     lattice_constant[1]=330.04;lattice_constant[2]=330.04;lattice_constant[3]=330.04;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.143;
     radius_PT=198;
     radius_covalent_vec=1.64;
     radius_covalent_PT=164;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.5658;
     RadiiSlatter=1.45;
     RadiiPyykko=1.47;
     ElectricalConductivity=6.7E6;
     electronegativity_vec=1.60;
     HardnessGhosh=2.8260;
     ElecNegPearson=4.0;
     ElecNegGhosh=3.485;
     ElectronAffinity_PT=86.1;
     Miedema_phi_star=4.00;
     Miedema_nws=1.62;
     Miedema_gamma_s=2700;
     Pettifor_scale=0.82;
     boiling_point=4744;
     melting_point=2477;
     VaporizationHeat_PT=690;
     SpecificHeat_PT=265;
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
     ShearModulus=38;
     YoungModulus=105;
     BulkModulus=170;
     PoissonRatio_PT=0.4;
     Miedema_BVm=18.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=2.76E-8;
     VolumeMagneticSusceptibility=0.000237;
     MolarMagneticSusceptibility=2.56E-9;
     Curie_point=NNN;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=8500;
     HHIR=8800;
     /*xray_scatt_vec=NNN;*/
     //Nb
 }
 // [AFLOW]STOP=Niobium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Molybdenum
 // Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum
 if(Z==42) {  // Molybdenum
     symbol_vec="Mo";
     name_vec="Molybdenum";
     Period=5;
     Group=6;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*95.94;
     MolarVolume=9.334E-6;
     volume_vec=15.6175;
     Miedema_Vm=4.4;
     valence_std_vec=6;
     valence_iupac_vec=6;
     valence_PT=6;
     Density_PT=10.28;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=0.000598128;
     lattice_constant[1]=314.7;lattice_constant[2]=314.7;lattice_constant[3]=314.7;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.136;
     radius_PT=190;
     radius_covalent_vec=1.54;
     radius_covalent_PT=154;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.4543;
     RadiiSlatter=1.45;
     RadiiPyykko=1.38;
     ElectricalConductivity=2E7;
     electronegativity_vec=2.16;
     HardnessGhosh=2.9221;
     ElecNegPearson=3.9;
     ElecNegGhosh=3.572;
     ElectronAffinity_PT=71.9;
     Miedema_phi_star=4.65;
     Miedema_nws=1.77;
     Miedema_gamma_s=2950;
     Pettifor_scale=0.88;
     boiling_point=4639;
     melting_point=2623;
     VaporizationHeat_PT=600;
     SpecificHeat_PT=251;
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
     ShearModulus=20;
     YoungModulus=329;
     BulkModulus=230;
     PoissonRatio_PT=0.31;
     Miedema_BVm=26.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.17E-8;
     VolumeMagneticSusceptibility=0.0001203;
     MolarMagneticSusceptibility=1.122E-9;
     Curie_point=NNN;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=2400;
     HHIR=5300;
     /*xray_scatt_vec=NNN;*/
     //Mo
 }
 // [AFLOW]STOP=Molybdenum
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Technetium
 // Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium
 if(Z==43) {  // Technetium
     symbol_vec="Tc";
     name_vec="Technetium";
     Period=5;
     Group=7;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*98.9062;
     MolarVolume=8.434782608696E-6;
     volume_vec=14.4670;
     Miedema_Vm=4.2;
     valence_std_vec=7;
     valence_iupac_vec=7;
     valence_PT=6;
     Density_PT=11.5;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=273.5;lattice_constant[2]=273.5;lattice_constant[3]=438.8;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=183;
     radius_covalent_vec=1.47;
     radius_covalent_PT=147;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.352;
     RadiiSlatter=1.35;
     RadiiPyykko=1.28;
     ElectricalConductivity=5E6;
     electronegativity_vec=1.90;
     HardnessGhosh=3.0184;
     ElecNegPearson=NNN;
     ElecNegGhosh=3.659;
     ElectronAffinity_PT=53;
     Miedema_phi_star=5.30;
     Miedema_nws=1.81;
     Miedema_gamma_s=3050;
     Pettifor_scale=0.94;
     boiling_point=4265;
     melting_point=2157;
     VaporizationHeat_PT=550;
     SpecificHeat_PT=63;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=26.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=3.42E-8;
     VolumeMagneticSusceptibility=0.0003933;
     MolarMagneticSusceptibility=3.352E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Tc JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Technetium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Ruthenium
 // Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium
 if(Z==44) {  // Ruthenium
     symbol_vec="Ru";
     name_vec="Ruthenium";
     Period=5;
     Group=8;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*101.07;
     MolarVolume=8.1706E-6;
     volume_vec=13.8390;
     Miedema_Vm=4.1;
     valence_std_vec=8;
     valence_iupac_vec=8;
     valence_PT=6;
     Density_PT=12.37;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.000406665;
     lattice_constant[1]=270.59;lattice_constant[2]=270.59;lattice_constant[3]=428.15;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.134;
     radius_PT=178;
     radius_covalent_vec=1.46;
     radius_covalent_PT=146;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.2579;
     RadiiSlatter=1.30;
     RadiiPyykko=1.25;
     ElectricalConductivity=1.4E7;
     electronegativity_vec=2.20;
     HardnessGhosh=3.1146;
     ElecNegPearson=4.5;
     ElecNegGhosh=3.745;
     ElectronAffinity_PT=101.3;
     Miedema_phi_star=5.40;
     Miedema_nws=1.83;
     Miedema_gamma_s=3050;
     Pettifor_scale=1.00;
     boiling_point=4150;
     melting_point=2334;
     VaporizationHeat_PT=580;
     SpecificHeat_PT=238;
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
     ShearModulus=173;
     YoungModulus=447;
     BulkModulus=220;
     PoissonRatio_PT=0.3;
     Miedema_BVm=26.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=5.42E-9;
     VolumeMagneticSusceptibility=0.000067;
     MolarMagneticSusceptibility=5.48E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=3200;
     HHIR=8000;
     /*xray_scatt_vec=NNN;*/
     //Ru JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Ruthenium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Rhodium
 // Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium
 if(Z==45) {  // Rhodium
     symbol_vec="Rh";
     name_vec="Rhodium";
     Period=5;
     Group=9;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*102.9055;
     MolarVolume=8.2655E-6;
     volume_vec=14.1731;
     Miedema_Vm=4.1;
     valence_std_vec=9;
     valence_iupac_vec=6;
     valence_PT=6;
     Density_PT=12.45;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=1.90706E-32;
     lattice_constant[1]=380.34;lattice_constant[2]=380.34;lattice_constant[3]=380.34;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.134;
     radius_PT=173;
     radius_covalent_vec=1.42;
     radius_covalent_PT=142;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.1711;
     RadiiSlatter=1.35;
     RadiiPyykko=1.25;
     ElectricalConductivity=2.3E7;
     electronegativity_vec=2.28;
     HardnessGhosh=3.2108;
     ElecNegPearson=4.3;
     ElecNegGhosh=3.832;
     ElectronAffinity_PT=109.7;
     Miedema_phi_star=5.40;
     Miedema_nws=1.76;
     Miedema_gamma_s=2750;
     Pettifor_scale=1.06;
     boiling_point=3695;
     melting_point=1964;
     VaporizationHeat_PT=495;
     SpecificHeat_PT=240;
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
     ShearModulus=150;
     YoungModulus=275;
     BulkModulus=380;
     PoissonRatio_PT=0.26;
     Miedema_BVm=23.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.36E-8;
     VolumeMagneticSusceptibility=0.0001693;
     MolarMagneticSusceptibility=1.4E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=3200;
     HHIR=8000;
     /*xray_scatt_vec=NNN;*/
     //Rh
 }
 // [AFLOW]STOP=Rhodium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Palladium
 // Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium
 if(Z==46) {  // Palladium
     symbol_vec="Pd";
     name_vec="Palladium";
     Period=5;
     Group=10;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*106.4;
     MolarVolume=8.8514E-6;
     volume_vec=15.4596;
     Miedema_Vm=4.3;
     valence_std_vec=10;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=12.023;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.000309478;
     lattice_constant[1]=389.07;lattice_constant[2]=389.07;lattice_constant[3]=389.07;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.137;
     radius_PT=169;
     radius_covalent_vec=1.39;
     radius_covalent_PT=139;
     radius_VanDerWaals_PT=163;
     RadiiGhosh08=2.0907;
     RadiiSlatter=1.40;
     RadiiPyykko=1.20;
     ElectricalConductivity=1E7;
     electronegativity_vec=2.20;
     HardnessGhosh=3.3069;
     ElecNegPearson=4.45;
     ElecNegGhosh=3.919;
     ElectronAffinity_PT=53.7;
     Miedema_phi_star=5.45;
     Miedema_nws=1.67;
     Miedema_gamma_s=2100;
     Pettifor_scale=1.12;
     boiling_point=2963;
     melting_point=1554.9;
     VaporizationHeat_PT=380;
     SpecificHeat_PT=240;
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
     ShearModulus=44;
     YoungModulus=121;
     BulkModulus=180;
     PoissonRatio_PT=0.39;
     Miedema_BVm=16.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=6.57E-8;
     VolumeMagneticSusceptibility=0.0007899;
     MolarMagneticSusceptibility=6.992E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=3200;
     HHIR=8000;
     /*xray_scatt_vec=NNN;*/
     //Pd
 }
 // [AFLOW]STOP=Palladium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Silver
 // Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver
 if(Z==47) {  // Silver
     symbol_vec="Ag";
     name_vec="Silver";
     Period=5;
     Group=11;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*107.8682;
     MolarVolume=0.000010283;
     volume_vec=18.0678;
     Miedema_Vm=4.7;
     valence_std_vec=11;
     valence_iupac_vec=4;
     valence_PT=1;
     Density_PT=10.49;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=8.57985E-05;
     lattice_constant[1]=408.53;lattice_constant[2]=408.53;lattice_constant[3]=408.53;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.144;
     radius_PT=165;
     radius_covalent_vec=1.45;
     radius_covalent_PT=145;
     radius_VanDerWaals_PT=172;
     RadiiGhosh08=2.016;
     RadiiSlatter=1.60;
     RadiiPyykko=1.28;
     ElectricalConductivity=6.2E7;
     electronegativity_vec=1.93;
     HardnessGhosh=3.4032;
     ElecNegPearson=4.44;
     ElecNegGhosh=4.006;
     ElectronAffinity_PT=125.6;
     Miedema_phi_star=4.45;
     Miedema_nws=1.39;
     Miedema_gamma_s=1250;
     Pettifor_scale=1.18;
     boiling_point=2162;
     melting_point=961.78;
     VaporizationHeat_PT=255;
     SpecificHeat_PT=235;
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
     ShearModulus=30;
     YoungModulus=85;
     BulkModulus=100;
     PoissonRatio_PT=0.37;
     Miedema_BVm=10.0;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-2.27E-9;
     VolumeMagneticSusceptibility=-0.0000238;
     MolarMagneticSusceptibility=-2.45E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=1200;
     HHIR=1400;
     xray_scatt_vec=47.18;
     //Ag JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Silver
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Cadmium
 // Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium
 if(Z==48) {  // Cadmium
     symbol_vec="Cd";
     name_vec="Cadmium";
     Period=5;
     Group=12;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*112.41;
     MolarVolume=0.000012996;
     volume_vec=22.0408;
     Miedema_Vm=5.5;
     valence_std_vec=12;
     valence_iupac_vec=2;
     valence_PT=2;
     Density_PT=8.65;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.000271603;
     lattice_constant[1]=297.94;lattice_constant[2]=297.94;lattice_constant[3]=561.86;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.150;
     radius_PT=161;
     radius_covalent_vec=1.44;
     radius_covalent_PT=144;
     radius_VanDerWaals_PT=158;
     RadiiGhosh08=1.9465;
     RadiiSlatter=1.55;
     RadiiPyykko=1.36;
     ElectricalConductivity=1.4E7;
     electronegativity_vec=1.69;
     HardnessGhosh=3.4994;
     ElecNegPearson=4.33;
     ElecNegGhosh=4.093;
     ElectronAffinity_PT=0;
     Miedema_phi_star=4.05;
     Miedema_nws=1.24;
     Miedema_gamma_s=780;
     Pettifor_scale=1.36;
     boiling_point=767;
     melting_point=321.07;
     VaporizationHeat_PT=100;
     SpecificHeat_PT=230;
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
     ShearModulus=19;
     YoungModulus=50;
     BulkModulus=42;
     PoissonRatio_PT=0.3;
     Miedema_BVm=6.10;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-2.3E-9;
     VolumeMagneticSusceptibility=-0.0000199;
     MolarMagneticSusceptibility=-2.59E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=1700;
     HHIR=1300;
     /*xray_scatt_vec=NNN;*/
     //Cd
 }
 // [AFLOW]STOP=Cadmium
// ********************************************************************************************************************************************************

// p-electron systems
// ********************************************************************************************************************************************************
 // [AFLOW]START=Indium
 // Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium
 if(Z==49) {  // Indium
     symbol_vec="In";
     name_vec="Indium";
     Period=5;
     Group=13;
     Series="PoorMetal";
     Block="p";
     mass_vec=AMU2KILOGRAM*114.82;
     MolarVolume=0.000015707;
     volume_vec=27.5233;
     Miedema_Vm=6.3;
     valence_std_vec=3;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=7.31;
     crystal_vec="fct";
     CrystalStr_PT="Centered_Tetragonal";
     space_group="I4/mmm";
     space_group_number=139;
     Pearson_coefficient=1.24494E-05;
     lattice_constant[1]=325.23;lattice_constant[2]=325.23;lattice_constant[3]=494.61;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.157;
     radius_PT=156;
     radius_covalent_vec=1.42;
     radius_covalent_PT=142;
     radius_VanDerWaals_PT=193;
     RadiiGhosh08=1.6934;
     RadiiSlatter=1.55;
     RadiiPyykko=1.42;
     ElectricalConductivity=1.2E7;
     electronegativity_vec=1.78;
     HardnessGhosh=3.9164;
     ElecNegPearson=3.1;
     ElecNegGhosh=4.469;
     ElectronAffinity_PT=28.9;
     Miedema_phi_star=3.90;
     Miedema_nws=1.17;
     Miedema_gamma_s=690;
     Pettifor_scale=1.60;
     boiling_point=2072;
     melting_point=156.6;
     VaporizationHeat_PT=230;
     SpecificHeat_PT=233;
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
     ShearModulus=NNN;
     YoungModulus=11;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=6.4;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-1.4E-9;
     VolumeMagneticSusceptibility=-0.0000102;
     MolarMagneticSusceptibility=-1.61E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=3300;
     HHIR=2000;
     /*xray_scatt_vec=NNN;*/
     //In
 }
 // [AFLOW]STOP=Indium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Tin
 // Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin
 if(Z==50) {  // Tin
     symbol_vec="Sn";
     name_vec="Tin";
     Period=5;
     Group=14;
     Series="PoorMetal";
     Block="p";
     mass_vec=AMU2KILOGRAM*118.69;
     MolarVolume=0.000016239;
     volume_vec=27.5555;
     Miedema_Vm=6.4;
     valence_std_vec=4;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=7.31;
     crystal_vec="bct";
     CrystalStr_PT="Centered_Tetragonal";
     space_group="I4_1/amd";
     space_group_number=141;
     Pearson_coefficient=0.000334085;
     lattice_constant[1]=583.18;lattice_constant[2]=583.18;lattice_constant[3]=318.19;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.158;
     radius_PT=145;
     radius_covalent_vec=1.39;
     radius_covalent_PT=139;
     radius_VanDerWaals_PT=217;
     RadiiGhosh08=1.4986;
     RadiiSlatter=1.45;
     RadiiPyykko=1.40;
     ElectricalConductivity=9.1E6;
     electronegativity_vec=1.96;
     HardnessGhosh=4.3332;
     ElecNegPearson=4.3;
     ElecNegGhosh=4.845;
     ElectronAffinity_PT=107.3;
     Miedema_phi_star=4.15;
     Miedema_nws=1.24;
     Miedema_gamma_s=710;
     Pettifor_scale=1.84;
     boiling_point=2602;
     melting_point=231.93;
     VaporizationHeat_PT=290;
     SpecificHeat_PT=217;
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
     ShearModulus=18;
     YoungModulus=50;
     BulkModulus=58;
     PoissonRatio_PT=0.36;
     Miedema_BVm=8.8;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-3.1E-9;
     VolumeMagneticSusceptibility=-0.0000227;
     MolarMagneticSusceptibility=-3.68E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=2600;
     HHIR=1600;
     /*xray_scatt_vec=NNN;*/
     //Sn
 }
 // [AFLOW]STOP=Tin
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Antimony
 // Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony
 if(Z==51) {  // Antimony
     symbol_vec="Sb";
     name_vec="Antimony";
     Period=5;
     Group=15;
     Series="Metalloid";
     Block="p";
     mass_vec=AMU2KILOGRAM*121.75;
     MolarVolume=0.000018181;
     volume_vec=27.1823;
     Miedema_Vm=6.6;
     valence_std_vec=5;
     valence_iupac_vec=5;
     valence_PT=5;
     Density_PT=6.697;
     crystal_vec="rhl";
     CrystalStr_PT="Simple_Trigonal";
     space_group="R_3m";
     space_group_number=166;
     Pearson_coefficient=6.60751E-05;
     lattice_constant[1]=430.7;lattice_constant[2]=430.7;lattice_constant[3]=1127.3;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.161;
     radius_PT=133;
     radius_covalent_vec=1.39;
     radius_covalent_PT=139;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.344;
     RadiiSlatter=1.45;
     RadiiPyykko=1.40;
     ElectricalConductivity=2.5E6;
     electronegativity_vec=2.05;
     HardnessGhosh=4.7501;
     ElecNegPearson=4.85;
     ElecNegGhosh=5.221;
     ElectronAffinity_PT=103.2;
     Miedema_phi_star=4.40;
     Miedema_nws=1.26;
     Miedema_gamma_s=680;
     Pettifor_scale=2.08;
     boiling_point=1587;
     melting_point=630.63;
     VaporizationHeat_PT=67;
     SpecificHeat_PT=207;
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
     ShearModulus=20;
     YoungModulus=55;
     BulkModulus=42;
     PoissonRatio_PT=NNN;
     Miedema_BVm=7.0;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-1.09E-8;
     VolumeMagneticSusceptibility=-0.000073;
     MolarMagneticSusceptibility=-1.327E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=7900;
     HHIR=3400;
     /*xray_scatt_vec=NNN;*/
     //Sb
 }
 // [AFLOW]STOP=Antimony
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Tellurium
 // Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium
 if(Z==52) {  // Tellurium
     symbol_vec="Te";
     name_vec="Tellurium";
     Period=5;
     Group=16;
     Series="Chalcogen";
     Block="p";
     mass_vec=AMU2KILOGRAM*127.6;
     MolarVolume=0.000020449;
     volume_vec=28.1993;
     Miedema_Vm=6.439;
     valence_std_vec=6;
     valence_iupac_vec=6;
     valence_PT=6;
     Density_PT=6.24;
     crystal_vec="hex";
     CrystalStr_PT="Simple_Trigonal";
     space_group="P3_121";
     space_group_number=152;
     Pearson_coefficient=0.000283934;
     lattice_constant[1]=445.72;lattice_constant[2]=445.72;lattice_constant[3]=592.9;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.143;
     radius_PT=123;
     radius_covalent_vec=1.38;
     radius_covalent_PT=138;
     radius_VanDerWaals_PT=206;
     RadiiGhosh08=1.2183;
     RadiiSlatter=1.40;
     RadiiPyykko=1.36;
     ElectricalConductivity=10000;
     electronegativity_vec=2.10;
     HardnessGhosh=5.1670;
     ElecNegPearson=5.49;
     ElecNegGhosh=5.597;
     ElectronAffinity_PT=190.2;
     Miedema_phi_star=4.72;
     Miedema_nws=1.31;
     Miedema_gamma_s=NNN;
     Pettifor_scale=2.32;
     boiling_point=988;
     melting_point=449.51;
     VaporizationHeat_PT=48;
     SpecificHeat_PT=201;
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
     ShearModulus=16;
     YoungModulus=43;
     BulkModulus=64;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-3.9E-9;
     VolumeMagneticSusceptibility=-0.0000243;
     MolarMagneticSusceptibility=-4.98E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=1.000991;
     HHIP=2900;
     HHIR=4900;
     /*xray_scatt_vec=NNN;*/
     //Te Table 27 of JUNKAI
 }
 // [AFLOW]STOP=Tellurium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Iodine
 // Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine
 if(Z==53) {  // Iodine
     symbol_vec="I";
     name_vec="Iodine";
     Period=5;
     Group=17;
     Series="Halogen";
     Block="p";
     mass_vec=AMU2KILOGRAM*126.9045;
     MolarVolume=0.000025689;
     volume_vec=34.9784;
     Miedema_Vm=8.72;
     valence_std_vec=7;
     valence_iupac_vec=7;
     valence_PT=7;
     Density_PT=4.94;
     crystal_vec="orc";
     CrystalStr_PT="Base_Orthorhombic";
     space_group="Cmca";
     space_group_number=64;
     Pearson_coefficient=0.0;
     lattice_constant[1]=718.02;lattice_constant[2]=471.02;lattice_constant[3]=981.03;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.136;
     radius_PT=115;
     radius_covalent_vec=1.39;
     radius_covalent_PT=139;
     radius_VanDerWaals_PT=198;
     RadiiGhosh08=1.1141;
     RadiiSlatter=1.40;
     RadiiPyykko=1.33;
     ElectricalConductivity=1E-7;
     electronegativity_vec=2.66;
     HardnessGhosh=5.5839;
     ElecNegPearson=6.76;
     ElecNegGhosh=5.973;
     ElectronAffinity_PT=295.2;
     Miedema_phi_star=5.33;
     Miedema_nws=0.17;
     Miedema_gamma_s=NNN;
     Pettifor_scale=2.56;
     boiling_point=184.3;
     melting_point=113.7;
     VaporizationHeat_PT=20.9;
     SpecificHeat_PT=429;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=7.7;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-4.5E-9;
     VolumeMagneticSusceptibility=-0.0000222;
     MolarMagneticSusceptibility=-1.14E-9;
     Curie_point=NNN;
     color_PT="SLATEGRAY";
     RefractiveIndex=NNN;
     HHIP=4900;
     HHIR=4800;
     /*xray_scatt_vec=NNN;*/
     //I interpolation phi_star, nws, Vm,
 }
 // [AFLOW]STOP=Iodine
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Xenon
 // Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon
 if(Z==54) {  // Xenon
     symbol_vec="Xe";
     name_vec="Xenon";
     Period=5;
     Group=18;
     Series="NobleGas";
     Block="p";
     mass_vec=AMU2KILOGRAM*131.3;
     MolarVolume=0.0223;
     volume_vec=-1.0000;
     Miedema_Vm=NNN;
     valence_std_vec=0;
     valence_iupac_vec=8;
     valence_PT=6;
     Density_PT=59E-4;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.000267781;
     lattice_constant[1]=620.23;lattice_constant[2]=620.23;lattice_constant[3]=620.23;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Gas";
     radius_vec=0.218;
     radius_PT=108;
     radius_covalent_vec=1.40;
     radius_covalent_PT=140;
     radius_VanDerWaals_PT=216;
     RadiiGhosh08=1.0263;
     RadiiSlatter=NNN;
     RadiiPyykko=1.31;
     ElectricalConductivity=NNN;
     electronegativity_vec=2.60;
     HardnessGhosh=6.0009;
     ElecNegPearson=NNN;
     ElecNegGhosh=6.349;
     ElectronAffinity_PT=0;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=-108;
     melting_point=-111.8;
     VaporizationHeat_PT=12.64;
     SpecificHeat_PT=158.32;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-4.3E-9;
     VolumeMagneticSusceptibility=-2.54E-8;
     MolarMagneticSusceptibility=-5.65E-10;
     Curie_point=NNN;
     color_PT="COLORLESS";
     RefractiveIndex=1.000702;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Xe JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Xenon
// ********************************************************************************************************************************************************

// ROW6
// s-electron systems
// ********************************************************************************************************************************************************
 // [AFLOW]START=Cesium
 // Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium
 if(Z==55) {  // Cesium
     symbol_vec="Cs";
     name_vec="Cesium";
     Period=6;
     Group=1;
     Series="AlkaliMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*132.9054;
     MolarVolume=0.000070732;
     volume_vec=117.281;
     Miedema_Vm=16.8;
     valence_std_vec=1;
     valence_iupac_vec=1;
     valence_PT=1;
     Density_PT=1.879;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=0.0;
     lattice_constant[1]=614.1;lattice_constant[2]=614.1;lattice_constant[3]=614.1;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.265;
     radius_PT=298;
     radius_covalent_vec=2.44;
     radius_covalent_PT=244;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=4.2433;
     RadiiSlatter=2.60;
     RadiiPyykko=2.32;
     ElectricalConductivity=5E6;
     electronegativity_vec=0.79;
     HardnessGhosh=0.6829;
     ElecNegPearson=2.18;
     ElecNegGhosh=4.196;
     ElectronAffinity_PT=45.5;
     Miedema_phi_star=1.95;
     Miedema_nws=0.55;
     Miedema_gamma_s=95;
     Pettifor_scale=0.25;
     boiling_point=671;
     melting_point=28.44;
     VaporizationHeat_PT=64;
     SpecificHeat_PT=242;
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
     ShearModulus=NNN;
     YoungModulus=1.7;
     BulkModulus=1.6;
     PoissonRatio_PT=NNN;
     Miedema_BVm=1.4;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=2.8E-9;
     VolumeMagneticSusceptibility=5.26E-6;
     MolarMagneticSusceptibility=3.72E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=6000;
     HHIR=6000;
     /*xray_scatt_vec=NNN;*/
     //Cs
 }
 // [AFLOW]STOP=Cesium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Barium
 // Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium
 if(Z==56) {  // Barium
     symbol_vec="Ba";
     name_vec="Barium";
     Period=6;
     Group=2;
     Series="AlkalineEarthMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*137.33;
     MolarVolume=0.000039125;
     volume_vec=62.6649;
     Miedema_Vm=11.3;
     valence_std_vec=2;
     valence_iupac_vec=2;
     valence_PT=2;
     Density_PT=3.51;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=6.23705E-05;
     lattice_constant[1]=502.8;lattice_constant[2]=502.8;lattice_constant[3]=502.8;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.217;
     radius_PT=253;
     radius_covalent_vec=2.15;
     radius_covalent_PT=215;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=3.2753;
     RadiiSlatter=2.15;
     RadiiPyykko=1.96;
     ElectricalConductivity=2.9E6;
     electronegativity_vec=0.89;
     HardnessGhosh=0.9201;
     ElecNegPearson=2.4;
     ElecNegGhosh=4.318;
     ElectronAffinity_PT=13.95;
     Miedema_phi_star=2.32;
     Miedema_nws=0.81;
     Miedema_gamma_s=370;
     Pettifor_scale=0.50;
     boiling_point=1870;
     melting_point=727;
     VaporizationHeat_PT=140;
     SpecificHeat_PT=205;
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
     ShearModulus=4.9;
     YoungModulus=13;
     BulkModulus=9.4;
     PoissonRatio_PT=NNN;
     Miedema_BVm=3.9;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.13E-8;
     VolumeMagneticSusceptibility=0.00003966;
     MolarMagneticSusceptibility=1.552E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=3000;
     HHIR=2300;
     /*xray_scatt_vec=NNN;*/
     //Ba
 }
 // [AFLOW]STOP=Barium
// ********************************************************************************************************************************************************

// d-electron systems: transition metals
// ********************************************************************************************************************************************************
 // [AFLOW]START=Lanthanium
 // Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium
 if(Z==57) {  // Lanthanium
     symbol_vec="La";
     name_vec="Lanthanium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*138.9055;
     MolarVolume=0.000022601;
     volume_vec=36.8495;
     Miedema_Vm=8.0;
     valence_std_vec=3;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=6.146;
     crystal_vec="hex";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=4.65323E-08;
     lattice_constant[1]=377.2;lattice_constant[2]=377.2;lattice_constant[3]=1214.4;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.187;
     radius_PT=NNN;
     radius_covalent_vec=2.07;
     radius_covalent_PT=207;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.6673;
     RadiiSlatter=1.95;
     RadiiPyykko=1.80;
     ElectricalConductivity=1.6E6;
     electronegativity_vec=1.10;
     HardnessGhosh=1.1571;
     ElecNegPearson=3.1;
     ElecNegGhosh=4.439;
     ElectronAffinity_PT=48;
     Miedema_phi_star=3.05;
     Miedema_nws=1.09;
     Miedema_gamma_s=900;
     Pettifor_scale=0.7480;
     boiling_point=3464;
     melting_point=919;
     VaporizationHeat_PT=400;
     SpecificHeat_PT=195;
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
     ShearModulus=14;
     YoungModulus=37;
     BulkModulus=28;
     PoissonRatio_PT=0.28;
     Miedema_BVm=5.5;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.1E-8;
     VolumeMagneticSusceptibility=0.00006761;
     MolarMagneticSusceptibility=1.528E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //La
 }
 // [AFLOW]STOP=Lanthanium
// ********************************************************************************************************************************************************

// lantanidies
// ********************************************************************************************************************************************************
 // [AFLOW]START=Cerium
 // Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium
 if(Z==58) {  // Cerium
     symbol_vec="Ce";
     name_vec="Cerium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*140.12;
     MolarVolume=0.000020947;
     volume_vec=26.4729;
     Miedema_Vm=7.76;
     valence_std_vec=4;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=6.689;
     crystal_vec="fcc";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=2.24956E-05;
     lattice_constant[1]=362;lattice_constant[2]=362;lattice_constant[3]=599;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.182;
     radius_PT=NNN;
     radius_covalent_vec=2.04;
     radius_covalent_PT=204;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.2494;
     RadiiSlatter=1.85;
     RadiiPyykko=1.63;
     ElectricalConductivity=1.4E6;
     electronegativity_vec=1.12;
     HardnessGhosh=1.3943;
     ElecNegPearson=NNN;
     ElecNegGhosh=4.561;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.18;
     Miedema_nws=1.19;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7460;
     boiling_point=3360;
     melting_point=798;
     VaporizationHeat_PT=350;
     SpecificHeat_PT=192;
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
     ShearModulus=14;
     YoungModulus=34;
     BulkModulus=22;
     PoissonRatio_PT=0.24;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=2.2E-7;
     VolumeMagneticSusceptibility=0.0014716;
     MolarMagneticSusceptibility=3.0826E-8;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Ce Pettifor linear interpolation // Miedema from Alonso-March.
 }
 // [AFLOW]STOP=Cerium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Praseodymium
 // Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium
 if(Z==59) {  // Praseodymium
     symbol_vec="Pr";
     name_vec="Praseodymium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*140.9077;
     MolarVolume=0.000021221;
     volume_vec=36.4987;
     Miedema_Vm=7.56;
     valence_std_vec=5;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=6.64;
     crystal_vec="hex";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=367.25;lattice_constant[2]=367.25;lattice_constant[3]=1183.54;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.183;
     radius_PT=247;
     radius_covalent_vec=2.03;
     radius_covalent_PT=203;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.9447;
     RadiiSlatter=1.85;
     RadiiPyykko=1.76;
     ElectricalConductivity=1.4E6;
     electronegativity_vec=1.13;
     HardnessGhosh=1.6315;
     ElecNegPearson=NNN;
     ElecNegGhosh=4.682;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.19;
     Miedema_nws=1.20;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7440;
     boiling_point=3290;
     melting_point=931;
     VaporizationHeat_PT=330;
     SpecificHeat_PT=193;
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
     ShearModulus=15;
     YoungModulus=37;
     BulkModulus=29;
     PoissonRatio_PT=0.28;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=4.23E-7;
     VolumeMagneticSusceptibility=0.0028087;
     MolarMagneticSusceptibility=5.9604E-8;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Pr Pettifor linear interpolation
 }
 // [AFLOW]STOP=Praseodymium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Neodymium
 // Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium
 if(Z==60) {  // Neodymium
     symbol_vec="Nd";
     name_vec="Neodymium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*144.24;
     MolarVolume=0.000020577;
     volume_vec=29.6719;
     Miedema_Vm=7.51;
     valence_std_vec=6;
     valence_iupac_vec=4;
     valence_PT=3;
     Density_PT=7.01;
     crystal_vec="hex";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.000231599;
     lattice_constant[1]=365.8;lattice_constant[2]=365.8;lattice_constant[3]=1179.9;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.182;
     radius_PT=206;
     radius_covalent_vec=2.01;
     radius_covalent_PT=201;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.7129;
     RadiiSlatter=1.85;
     RadiiPyykko=1.74;
     ElectricalConductivity=1.6E6;
     electronegativity_vec=1.14;
     HardnessGhosh=1.8684;
     ElecNegPearson=NNN;
     ElecNegGhosh=4.804;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.19;
     Miedema_nws=1.20;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7420;
     boiling_point=3100;
     melting_point=1021;
     VaporizationHeat_PT=285;
     SpecificHeat_PT=190;
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
     ShearModulus=16;
     YoungModulus=41;
     BulkModulus=32;
     PoissonRatio_PT=0.28;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=4.8E-7;
     VolumeMagneticSusceptibility=0.0033648;
     MolarMagneticSusceptibility=6.9235E-8;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Nd Pettifor linear interpolation JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Neodymium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Promethium
 // Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium
 if(Z==61) {  // Promethium
     symbol_vec="Pm";
     name_vec="Promethium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*146.92;
     MolarVolume=0.00001996145374449;
     volume_vec=34.6133;
     Miedema_Vm=7.43;
     valence_std_vec=7;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=7.264;
     crystal_vec="hex";
     CrystalStr_PT="NNN";
     space_group="NNN";
     space_group_number=NNN;
     Pearson_coefficient=0.0;
     lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
     lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=205;
     radius_covalent_vec=1.99;
     radius_covalent_PT=199;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.5303;
     RadiiSlatter=1.85;
     RadiiPyykko=1.73;
     ElectricalConductivity=1.3E6;
     electronegativity_vec=1.13;
     HardnessGhosh=2.1056;
     ElecNegPearson=NNN;
     ElecNegGhosh=4.925;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.19;
     Miedema_nws=1.21;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7400;
     boiling_point=3000;
     melting_point=1100;
     VaporizationHeat_PT=290;
     SpecificHeat_PT=NNN;
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
     ShearModulus=18;
     YoungModulus=46;
     BulkModulus=33;
     PoissonRatio_PT=0.28;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     // Pm Pettifor linear interpolation
 }
 // [AFLOW]STOP=Promethium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Samarium
 // Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium
 if(Z==62) {  // Samarium
     symbol_vec="Sm";
     name_vec="Samarium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*150.4;
     MolarVolume=0.000020449;
     volume_vec=33.9484;
     Miedema_Vm=7.37;
     valence_std_vec=8;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=7.353;
     crystal_vec="rhl";
     CrystalStr_PT="Simple_Trigonal";
     space_group="R_3m";
     space_group_number=166;
     Pearson_coefficient=0.000334686;
     lattice_constant[1]=362.1;lattice_constant[2]=362.1;lattice_constant[3]=2625;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.181;
     radius_PT=238;
     radius_covalent_vec=1.98;
     radius_covalent_PT=198;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.3830;
     RadiiSlatter=1.85;
     RadiiPyykko=1.72;
     ElectricalConductivity=1.1E6;
     electronegativity_vec=1.17;
     HardnessGhosh=2.3427;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.047;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.20;
     Miedema_nws=1.21;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7380;
     boiling_point=1803;
     melting_point=1072;
     VaporizationHeat_PT=175;
     SpecificHeat_PT=196;
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
     ShearModulus=20;
     YoungModulus=50;
     BulkModulus=38;
     PoissonRatio_PT=0.27;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.11E-7;
     VolumeMagneticSusceptibility=0.00081618;
     MolarMagneticSusceptibility=1.669E-8;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Sm Pettifor linear interpolation
 }
 // [AFLOW]STOP=Samarium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Europium
 // Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium
 if(Z==63) {  // Europium
     symbol_vec="Eu";
     name_vec="Europium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*151.96;
     MolarVolume=0.000028979;
     volume_vec=43.1719;
     Miedema_Vm=7.36;
     valence_std_vec=9;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=5.244;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=4.32857E-05;
     lattice_constant[1]=458.1;lattice_constant[2]=458.1;lattice_constant[3]=458.1;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.204;
     radius_PT=231;
     radius_covalent_vec=1.98;
     radius_covalent_PT=198;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.2615;
     RadiiSlatter=1.85;
     RadiiPyykko=1.68;
     ElectricalConductivity=1.1E6;
     electronegativity_vec=1.20;
     HardnessGhosh=2.5798;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.168;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.20;
     Miedema_nws=1.21;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7360;
     boiling_point=1527;
     melting_point=822;
     VaporizationHeat_PT=175;
     SpecificHeat_PT=182;
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
     ShearModulus=7.9;
     YoungModulus=18;
     BulkModulus=8.3;
     PoissonRatio_PT=0.15;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=2.76E-7;
     VolumeMagneticSusceptibility=0.0014473;
     MolarMagneticSusceptibility=4.1942E-8;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Eu Pettifor linear interpolation
 }
 // [AFLOW]STOP=Europium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Gadolinium
 // Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium
 if(Z==64) {  // Gadolinium
     symbol_vec="Gd";
     name_vec="Gadolinium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*157.25;
     MolarVolume=0.000019903;
     volume_vec=32.5777;
     Miedema_Vm=7.34;
     valence_std_vec=10;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=7.901;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.000127674;
     lattice_constant[1]=363.6;lattice_constant[2]=363.6;lattice_constant[3]=578.26;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.180;
     radius_PT=233;
     radius_covalent_vec=1.96;
     radius_covalent_PT=196;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.1596;
     RadiiSlatter=1.80;
     RadiiPyykko=1.69;
     ElectricalConductivity=770000;
     electronegativity_vec=1.20;
     HardnessGhosh=2.8170;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.290;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.20;
     Miedema_nws=1.21;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7340;
     boiling_point=3250;
     melting_point=1313;
     VaporizationHeat_PT=305;
     SpecificHeat_PT=240;
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
     ShearModulus=22;
     YoungModulus=55;
     BulkModulus=38;
     PoissonRatio_PT=0.26;
     Miedema_BVm=NNN;
     MagneticType_PT="Ferromagnetic";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=292;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     // Gd Pettifor linear interpolation
 }
 // [AFLOW]STOP=Gadolinium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Terbium
 // Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium
 if(Z==65) {  // Terbium
     symbol_vec="Tb";
     name_vec="Terbium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*158.9254;
     MolarVolume=0.000019336;
     volume_vec=32.0200;
     Miedema_Vm=7.20;
     valence_std_vec=11;
     valence_iupac_vec=4;
     valence_PT=3;
     Density_PT=8.219;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=360.1;lattice_constant[2]=360.1;lattice_constant[3]=569.36;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.177;
     radius_PT=225;
     radius_covalent_vec=1.94;
     radius_covalent_PT=194;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.0730;
     RadiiSlatter=1.75;
     RadiiPyykko=1.68;
     ElectricalConductivity=830000;
     electronegativity_vec=1.10;
     HardnessGhosh=3.0540;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.411;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.21;
     Miedema_nws=1.22;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7320;
     boiling_point=3230;
     melting_point=1356;
     VaporizationHeat_PT=295;
     SpecificHeat_PT=182;
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
     ShearModulus=22;
     YoungModulus=56;
     BulkModulus=38.7;
     PoissonRatio_PT=0.26;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=0.0000136;
     VolumeMagneticSusceptibility=0.1117784;
     MolarMagneticSusceptibility=2.161385E-6;
     Curie_point=222;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     // Tb Pettifor linear interpolation
 }
 // [AFLOW]STOP=Terbium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Dysprosium
 // Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium
 if(Z==66) {  // Dysprosium
     symbol_vec="Dy";
     name_vec="Dysprosium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*162.5;
     MolarVolume=0.000019004;
     volume_vec=31.5096;
     Miedema_Vm=7.12;
     valence_std_vec=12;
     valence_iupac_vec=4;
     valence_PT=3;
     Density_PT=8.551;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=5.20771E-05;
     lattice_constant[1]=359.3;lattice_constant[2]=359.3;lattice_constant[3]=565.37;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.177;
     radius_PT=228;
     radius_covalent_vec=1.92;
     radius_covalent_PT=192;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.9984;
     RadiiSlatter=1.75;
     RadiiPyykko=1.67;
     ElectricalConductivity=1.1E6;
     electronegativity_vec=1.22;
     HardnessGhosh=3.2912;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.533;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.21;
     Miedema_nws=1.22;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7300;
     boiling_point=2567;
     melting_point=1412;
     VaporizationHeat_PT=280;
     SpecificHeat_PT=167;
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
     ShearModulus=25;
     YoungModulus=61;
     BulkModulus=41;
     PoissonRatio_PT=0.25;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=5.45E-6;
     VolumeMagneticSusceptibility=0.046603;
     MolarMagneticSusceptibility=8.85625E-7;
     Curie_point=87;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Dy Pettifor linear interpolation JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Dysprosium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Holmium
 // Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium
 if(Z==67) {  // Holmium
     symbol_vec="Ho";
     name_vec="Holmium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*164.9304;
     MolarVolume=0.000018753;
     volume_vec=31.0155;
     Miedema_Vm=7.06;
     valence_std_vec=13;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=8.795;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=2.96961E-32;
     lattice_constant[1]=357.73;lattice_constant[2]=357.73;lattice_constant[3]=561.58;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.176;
     radius_PT=226;
     radius_covalent_vec=1.92;
     radius_covalent_PT=192;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.9335;
     RadiiSlatter=1.75;
     RadiiPyykko=1.66;
     ElectricalConductivity=1.1E6;
     electronegativity_vec=1.23;
     HardnessGhosh=3.5283;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.654;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.22;
     Miedema_nws=1.22;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7280;
     boiling_point=2700;
     melting_point=1474;
     VaporizationHeat_PT=265;
     SpecificHeat_PT=165;
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
     ShearModulus=26;
     YoungModulus=64;
     BulkModulus=40;
     PoissonRatio_PT=0.23;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=5.49E-6;
     VolumeMagneticSusceptibility=0.0482845;
     MolarMagneticSusceptibility=9.05467E-7;
     Curie_point=20;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Ho Pettifor linear interpolation
 }
 // [AFLOW]STOP=Holmium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Erbium
 // Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium
 if(Z==68) {  // Erbium
     symbol_vec="Er";
     name_vec="Erbium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*167.26;
     MolarVolume=0.000018449;
     volume_vec=30.5431;
     Miedema_Vm=6.98;
     valence_std_vec=14;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=9.066;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=7.24618E-05;
     lattice_constant[1]=355.88;lattice_constant[2]=355.88;lattice_constant[3]=558.74;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.175;
     radius_PT=226;
     radius_covalent_vec=1.89;
     radius_covalent_PT=189;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.8765;
     RadiiSlatter=1.75;
     RadiiPyykko=1.65;
     ElectricalConductivity=1.2E6;
     electronegativity_vec=1.24;
     HardnessGhosh=3.7655;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.776;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.22;
     Miedema_nws=1.23;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7260;
     boiling_point=2868;
     melting_point=1497;
     VaporizationHeat_PT=285;
     SpecificHeat_PT=168;
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
     ShearModulus=28;
     YoungModulus=70;
     BulkModulus=44;
     PoissonRatio_PT=0.24;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=3.77E-6;
     VolumeMagneticSusceptibility=0.0341788;
     MolarMagneticSusceptibility=6.30566E-7;
     Curie_point=32;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Er Pettifor linear interpolation
 }
 // [AFLOW]STOP=Erbium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Thulium
 // Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium
 if(Z==69) {  // Thulium
     symbol_vec="Tm";
     name_vec="Thulium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*168.9342;
     MolarVolume=0.000018126;
     volume_vec=30.0016;
     Miedema_Vm=6.90;
     valence_std_vec=15;
     valence_iupac_vec=4;
     valence_PT=3;
     Density_PT=9.32;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=353.75;lattice_constant[2]=353.75;lattice_constant[3]=555.46;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.174;
     radius_PT=222;
     radius_covalent_vec=1.90;
     radius_covalent_PT=190;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.8261;
     RadiiSlatter=1.75;
     RadiiPyykko=1.64;
     ElectricalConductivity=1.4E6;
     electronegativity_vec=1.25;
     HardnessGhosh=4.0026;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.897;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.22;
     Miedema_nws=1.23;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7240;
     boiling_point=1950;
     melting_point=1545;
     VaporizationHeat_PT=250;
     SpecificHeat_PT=160;
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
     ShearModulus=31;
     YoungModulus=74;
     BulkModulus=45;
     PoissonRatio_PT=0.21;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.99E-6;
     VolumeMagneticSusceptibility=0.0185488;
     MolarMagneticSusceptibility=3.36179E-7;
     Curie_point=25;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Tm Pettifor linear interpolation JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Thulium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Ytterbium
 // Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium
 if(Z==70) {  // Ytterbium
     symbol_vec="Yb";
     name_vec="Ytterbium";
     Period=6;
     Group=NNN;
     Series="Lanthanide";
     Block="f";
     mass_vec=AMU2KILOGRAM*173.04;
     MolarVolume=0.000026339;
     volume_vec=39.4395;
     Miedema_Vm=6.86;
     valence_std_vec=16;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=6.57;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=8.54557E-05;
     lattice_constant[1]=548.47;lattice_constant[2]=548.47;lattice_constant[3]=548.47;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.193;
     radius_PT=222;
     radius_covalent_vec=1.87;
     radius_covalent_PT=187;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.7812;
     RadiiSlatter=1.75;
     RadiiPyykko=1.70;
     ElectricalConductivity=3.6E6;
     electronegativity_vec=1.10;
     HardnessGhosh=4.2395;
     ElecNegPearson=NNN;
     ElecNegGhosh=6.019;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.22;
     Miedema_nws=1.23;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7220;
     boiling_point=1196;
     melting_point=819;
     VaporizationHeat_PT=160;
     SpecificHeat_PT=154;
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
     ShearModulus=10;
     YoungModulus=24;
     BulkModulus=31;
     PoissonRatio_PT=0.21;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=5.9E-9;
     VolumeMagneticSusceptibility=0.0000388;
     MolarMagneticSusceptibility=1.02E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Yb Pettifor linear interpolation
 }
 // [AFLOW]STOP=Ytterbium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Lutetium
 // Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium
 if(Z==71) {  // Lutetium
     symbol_vec="Lu";
     name_vec="Lutetium";
     Period=6;
     Group=3;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*174.967;
     MolarVolume=0.000017779;
     volume_vec=29.3515;
     Miedema_Vm=6.81;
     valence_std_vec=17;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=9.841;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=8.27273E-07;
     lattice_constant[1]=350.31;lattice_constant[2]=350.31;lattice_constant[3]=555.09;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.173;
     radius_PT=217;
     radius_covalent_vec=1.87;
     radius_covalent_PT=187;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.7409;
     RadiiSlatter=1.75;
     RadiiPyykko=1.62;
     ElectricalConductivity=1.8E6;
     electronegativity_vec=1.27;
     HardnessGhosh=4.4766;
     ElecNegPearson=NNN;
     ElecNegGhosh=6.140;
     ElectronAffinity_PT=50;
     Miedema_phi_star=3.22;
     Miedema_nws=1.24;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0.7200;
     boiling_point=3402;
     melting_point=1663;
     VaporizationHeat_PT=415;
     SpecificHeat_PT=154;
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
     ShearModulus=27;
     YoungModulus=67;
     BulkModulus=48;
     PoissonRatio_PT=0.26;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.2E-9;
     VolumeMagneticSusceptibility=0.0000118;
     MolarMagneticSusceptibility=2.1E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=9500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Lu
 }
 // [AFLOW]STOP=Lutetium
// ********************************************************************************************************************************************************

// d-electron systems: transition metals
// ********************************************************************************************************************************************************
 // [AFLOW]START=Hafnium
 // Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium
 if(Z==72) {  // Hafnium
     symbol_vec="Hf";
     name_vec="Hafnium";
     Period=6;
     Group=4;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*178.49;
     MolarVolume=0.0000134102;
     volume_vec=22.0408;
     Miedema_Vm=5.6;
     valence_std_vec=4;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=13.31;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=5.25384E-05;
     lattice_constant[1]=319.64;lattice_constant[2]=319.64;lattice_constant[3]=505.11;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.159;
     radius_PT=208;
     radius_covalent_vec=1.75;
     radius_covalent_PT=175;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.7056;
     RadiiSlatter=1.55;
     RadiiPyykko=1.52;
     ElectricalConductivity=3.3E6;
     electronegativity_vec=1.30;
     HardnessGhosh=4.7065;
     ElecNegPearson=3.8;
     ElecNegGhosh=6.258;
     ElectronAffinity_PT=0;
     Miedema_phi_star=3.55;
     Miedema_nws=1.43;
     Miedema_gamma_s=2200;
     Pettifor_scale=0.775;
     boiling_point=4603;
     melting_point=2233;
     VaporizationHeat_PT=630;
     SpecificHeat_PT=144;
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
     ShearModulus=30;
     YoungModulus=78;
     BulkModulus=110;
     PoissonRatio_PT=0.37;
     Miedema_BVm=15.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=5.3E-9;
     VolumeMagneticSusceptibility=0.0000705;
     MolarMagneticSusceptibility=9.46E-10;
     Curie_point=NNN;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=3400;
     HHIR=2600;
     /*xray_scatt_vec=NNN;*/
     //Hf
 }
 // [AFLOW]STOP=Hafnium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Tantalum
 // Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum
 if(Z==73) {  // Tantalum
     symbol_vec="Ta";
     name_vec="Tantalum";
     Period=6;
     Group=5;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*180.9479;
     MolarVolume=0.0000108677;
     volume_vec=18.1100;
     Miedema_Vm=4.9;
     valence_std_vec=5;
     valence_iupac_vec=5;
     valence_PT=5;
     Density_PT=16.65;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=3.66845E-09;
     lattice_constant[1]=330.13;lattice_constant[2]=330.13;lattice_constant[3]=330.13;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.147;
     radius_PT=200;
     radius_covalent_vec=1.70;
     radius_covalent_PT=170;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.6716;
     RadiiSlatter=1.45;
     RadiiPyykko=1.46;
     ElectricalConductivity=7.7E6;
     electronegativity_vec=1.50;
     HardnessGhosh=4.9508;
     ElecNegPearson=4.11;
     ElecNegGhosh=6.383;
     ElectronAffinity_PT=31;
     Miedema_phi_star=4.05;
     Miedema_nws=1.63;
     Miedema_gamma_s=3050;
     Pettifor_scale=0.83;
     boiling_point=5458;
     melting_point=3017;
     VaporizationHeat_PT=736;
     SpecificHeat_PT=140;
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
     ShearModulus=67;
     YoungModulus=186;
     BulkModulus=200;
     PoissonRatio_PT=0.34;
     Miedema_BVm=22.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.07E-8;
     VolumeMagneticSusceptibility=0.0001782;
     MolarMagneticSusceptibility=1.936E-9;
     Curie_point=NNN;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=2300;
     HHIR=4800;
     /*xray_scatt_vec=NNN;*/
     //Ta
 }
 // [AFLOW]STOP=Tantalum
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Tungsten
 // Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten
 if(Z==74) {  // Tungsten
     symbol_vec="W";
     name_vec="Tungsten";
     Period=6;
     Group=6;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*183.85;
     MolarVolume=9.5501E-6;
     volume_vec=15.9387;
     Miedema_Vm=4.5;
     valence_std_vec=6;
     valence_iupac_vec=6;
     valence_PT=6;
     Density_PT=19.25;
     crystal_vec="bcc";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=6.96679E-05;
     lattice_constant[1]=316.52;lattice_constant[2]=316.52;lattice_constant[3]=316.52;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.137;
     radius_PT=193;
     radius_covalent_vec=1.62;
     radius_covalent_PT=162;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.6416;
     RadiiSlatter=1.35;
     RadiiPyykko=1.37;
     ElectricalConductivity=2E7;
     electronegativity_vec=2.36;
     HardnessGhosh=5.1879;
     ElecNegPearson=4.40;
     ElecNegGhosh=6.505;
     ElectronAffinity_PT=78.6;
     Miedema_phi_star=4.80;
     Miedema_nws=1.81;
     Miedema_gamma_s=3300;
     Pettifor_scale=0.885;
     boiling_point=5555;
     melting_point=3422;
     VaporizationHeat_PT=800;
     SpecificHeat_PT=132;
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
     ShearModulus=161;
     YoungModulus=411;
     BulkModulus=310;
     PoissonRatio_PT=0.28;
     Miedema_BVm=31.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=4.59E-9;
     VolumeMagneticSusceptibility=0.0000884;
     MolarMagneticSusceptibility=8.44E-10;
     Curie_point=NNN;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=7000;
     HHIR=4300;
     /*xray_scatt_vec=NNN;*/
     //W
 }
 // [AFLOW]STOP=Tungsten
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Rhenium
 // Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium
 if(Z==75) {  // Rhenium
     symbol_vec="Re";
     name_vec="Rhenium";
     Period=6;
     Group=7;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*186.2;
     MolarVolume=8.85856E-6;
     volume_vec=14.8941;
     Miedema_Vm=4.3;
     valence_std_vec=7;
     valence_iupac_vec=7;
     valence_PT=7;
     Density_PT=21.02;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=2.70849E-05;
     lattice_constant[1]=276.1;lattice_constant[2]=276.1;lattice_constant[3]=445.6;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.138;
     radius_PT=188;
     radius_covalent_vec=1.51;
     radius_covalent_PT=151;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.6141;
     RadiiSlatter=1.35;
     RadiiPyykko=1.31;
     ElectricalConductivity=5.6E6;
     electronegativity_vec=1.90;
     HardnessGhosh=5.4256;
     ElecNegPearson=4.02;
     ElecNegGhosh=6.626;
     ElectronAffinity_PT=14.5;
     Miedema_phi_star=5.40;
     Miedema_nws=1.86;
     Miedema_gamma_s=3650;
     Pettifor_scale=0.94;
     boiling_point=5596;
     melting_point=3186;
     VaporizationHeat_PT=705;
     SpecificHeat_PT=137;
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
     ShearModulus=178;
     YoungModulus=463;
     BulkModulus=370;
     PoissonRatio_PT=0.3;
     Miedema_BVm=33.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=4.56E-9;
     VolumeMagneticSusceptibility=0.0000959;
     MolarMagneticSusceptibility=8.49E-10;
     Curie_point=NNN;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=3300;
     HHIR=3300;
     /*xray_scatt_vec=NNN;*/
     //Re
 }
 // [AFLOW]STOP=Rhenium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Osmium
 // Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium
 if(Z==76) {  // Osmium
     symbol_vec="Os";
     name_vec="Osmium";
     Period=6;
     Group=8;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*190.2;
     MolarVolume=8.421E-6;
     volume_vec=14.2403;
     Miedema_Vm=4.2;
     valence_std_vec=8;
     valence_iupac_vec=8;
     valence_PT=6;
     Density_PT=22.59;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=7.45234E-05;
     lattice_constant[1]=273.44;lattice_constant[2]=273.44;lattice_constant[3]=431.73;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.135;
     radius_PT=185;
     radius_covalent_vec=1.44;
     radius_covalent_PT=144;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.5890;
     RadiiSlatter=1.30;
     RadiiPyykko=1.29;
     ElectricalConductivity=1.2E7;
     electronegativity_vec=2.20;
     HardnessGhosh=5.6619;
     ElecNegPearson=4.9;
     ElecNegGhosh=6.748;
     ElectronAffinity_PT=106.1;
     Miedema_phi_star=5.40;
     Miedema_nws=1.85;
     Miedema_gamma_s=3500;
     Pettifor_scale=0.995;
     boiling_point=5012;
     melting_point=3033;
     VaporizationHeat_PT=630;
     SpecificHeat_PT=130;
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
     ShearModulus=222;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=0.25;
     Miedema_BVm=35.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=6E-10;
     VolumeMagneticSusceptibility=0.000014;
     MolarMagneticSusceptibility=1.1E-10;
     Curie_point=NNN;
     color_PT="SLATEGRAY";
     RefractiveIndex=NNN;
     HHIP=5500;
     HHIR=9100;
     /*xray_scatt_vec=NNN;*/
     //Os JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Osmium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Iridium
 // Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium
 if(Z==77) {  // Iridium
     symbol_vec="Ir";
     name_vec="Iridium";
     Period=6;
     Group=9;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*192.22;
     MolarVolume=8.5203E-6;
     volume_vec=14.5561;
     Miedema_Vm=4.2;
     valence_std_vec=9;
     valence_iupac_vec=8;
     valence_PT=6;
     Density_PT=22.56;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=2.53787E-05;
     lattice_constant[1]=383.9;lattice_constant[2]=383.9;lattice_constant[3]=383.9;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.135;
     radius_PT=180;
     radius_covalent_vec=1.41;
     radius_covalent_PT=141;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.5657;
     RadiiSlatter=1.35;
     RadiiPyykko=1.22;
     ElectricalConductivity=2.1E7;
     electronegativity_vec=2.20;
     HardnessGhosh=5.9000;
     ElecNegPearson=5.4;
     ElecNegGhosh=6.831;
     ElectronAffinity_PT=151;
     Miedema_phi_star=5.55;
     Miedema_nws=1.83;
     Miedema_gamma_s=3100;
     Pettifor_scale=1.05;
     boiling_point=4428;
     melting_point=2466;
     VaporizationHeat_PT=560;
     SpecificHeat_PT=131;
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
     ShearModulus=210;
     YoungModulus=528;
     BulkModulus=320;
     PoissonRatio_PT=0.26;
     Miedema_BVm=25.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.67E-9;
     VolumeMagneticSusceptibility=0.0000377;
     MolarMagneticSusceptibility=3.21E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=5500;
     HHIR=9100;
     /*xray_scatt_vec=NNN;*/
     //Ir JUNKAI CHANGED VALENCE
 }
 // [AFLOW]STOP=Iridium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Platinum
 // Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum
 if(Z==78) {  // Platinum
     symbol_vec="Pt";
     name_vec="Platinum";
     Period=6;
     Group=10;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*195.09;
     MolarVolume=9.0948E-6;
     volume_vec=15.7298;
     Miedema_Vm=4.4;
     valence_std_vec=10;
     valence_iupac_vec=6;
     valence_PT=6;
     Density_PT=21.45;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=3.39206E-05;
     lattice_constant[1]=392.42;lattice_constant[2]=392.42;lattice_constant[3]=392.42;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.138;
     radius_PT=177;
     radius_covalent_vec=1.36;
     radius_covalent_PT=136;
     radius_VanDerWaals_PT=175;
     RadiiGhosh08=0.5443;
     RadiiSlatter=1.35;
     RadiiPyykko=1.23;
     ElectricalConductivity=9.4E6;
     electronegativity_vec=2.28;
     HardnessGhosh=6.1367;
     ElecNegPearson=5.6;
     ElecNegGhosh=6.991;
     ElectronAffinity_PT=205.3;
     Miedema_phi_star=5.65;
     Miedema_nws=1.78;
     Miedema_gamma_s=2550;
     Pettifor_scale=1.105;
     boiling_point=3825;
     melting_point=1768.3;
     VaporizationHeat_PT=490;
     SpecificHeat_PT=133;
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
     ShearModulus=61;
     YoungModulus=168;
     BulkModulus=230;
     PoissonRatio_PT=0.38;
     Miedema_BVm=18.0;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=1.22E-8;
     VolumeMagneticSusceptibility=0.0002573;
     MolarMagneticSusceptibility=2.38E-9;
     Curie_point=NNN;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=5500;
     HHIR=9100;
     /*xray_scatt_vec=NNN;*/
     //Pt
 }
 // [AFLOW]STOP=Platinum
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Gold
 // Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold
 if(Z==79) {  // Gold
     symbol_vec="Au";
     name_vec="Gold";
     Period=6;
     Group=11;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*196.9665;
     MolarVolume=0.00001021;
     volume_vec=18.1904;
     Miedema_Vm=4.7;
     valence_std_vec=11;
     valence_iupac_vec=5;
     valence_PT=5;
     Density_PT=19.3;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=2.08217E-32;
     lattice_constant[1]=407.82;lattice_constant[2]=407.82;lattice_constant[3]=407.82;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.144;
     radius_PT=174;
     radius_covalent_vec=1.36;
     radius_covalent_PT=136;
     radius_VanDerWaals_PT=166;
     RadiiGhosh08=0.5244;
     RadiiSlatter=1.35;
     RadiiPyykko=1.24;
     ElectricalConductivity=4.5E7;
     electronegativity_vec=2.54;
     HardnessGhosh=6.3741;
     ElecNegPearson=5.77;
     ElecNegGhosh=7.112;
     ElectronAffinity_PT=222.8;
     Miedema_phi_star=5.15;
     Miedema_nws=1.57;
     Miedema_gamma_s=1550;
     Pettifor_scale=1.16;
     boiling_point=2856;
     melting_point=1064.18;
     VaporizationHeat_PT=330;
     SpecificHeat_PT=129.1;
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
     ShearModulus=27;
     YoungModulus=78;
     BulkModulus=220;
     PoissonRatio_PT=0.44;
     Miedema_BVm=18.0;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-1.78E-9;
     VolumeMagneticSusceptibility=-0.0000344;
     MolarMagneticSusceptibility=-3.51E-10;
     Curie_point=NNN;
     color_PT="GOLD";
     RefractiveIndex=NNN;
     HHIP=1100;
     HHIR=1000;
     xray_scatt_vec=74.99;
     //Au
 }
 // [AFLOW]STOP=Gold
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Mercury
 // Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury
 if(Z==80) {  // Mercury
     symbol_vec="Hg";
     name_vec="Mercury";
     Period=6;
     Group=12;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*200.59;
     MolarVolume=0.0000148213;
     volume_vec=29.7156;
     Miedema_Vm=5.8;
     valence_std_vec=12;
     valence_iupac_vec=4;
     valence_PT=2;
     Density_PT=13.534;
     crystal_vec="rhl";
     CrystalStr_PT="Simple_Trigonal";
     space_group="R_3m";
     space_group_number=166;
     Pearson_coefficient=6.52519E-05;
     lattice_constant[1]=300.5;lattice_constant[2]=300.5;lattice_constant[3]=300.5;
     lattice_angle[1]=1.23081;lattice_angle[2]=1.23081;lattice_angle[3]=1.23081;
     phase="Liquid";
     radius_vec=0.150;
     radius_PT=171;
     radius_covalent_vec=1.32;
     radius_covalent_PT=132;
     radius_VanDerWaals_PT=155;
     RadiiGhosh08=0.5060;
     RadiiSlatter=1.50;
     RadiiPyykko=1.33;
     ElectricalConductivity=1E6;
     electronegativity_vec=2.00;
     HardnessGhosh=6.6103;
     ElecNegPearson=4.91;
     ElecNegGhosh=7.233;
     ElectronAffinity_PT=0;
     Miedema_phi_star=4.20;
     Miedema_nws=1.24;
     Miedema_gamma_s=610;
     Pettifor_scale=1.32;
     boiling_point=356.73;
     melting_point=-38.83;
     VaporizationHeat_PT=59.2;
     SpecificHeat_PT=139.5;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=25;
     PoissonRatio_PT=NNN;
     Miedema_BVm=4.0;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-2.1E-9;
     VolumeMagneticSusceptibility=-0.0000284;
     MolarMagneticSusceptibility=-4.21E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=1.000933;
     HHIP=5500;
     HHIR=3100;
     /*xray_scatt_vec=NNN;*/
     //Hg
 }
 // [AFLOW]STOP=Mercury
// ********************************************************************************************************************************************************

// p-electron systems
// ********************************************************************************************************************************************************
 // [AFLOW]START=Thallium
 // Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium
 if(Z==81) {  // Thallium
     symbol_vec="Tl";
     name_vec="Thallium";
     Period=6;
     Group=13;
     Series="PoorMetal";
     Block="p";
     mass_vec=AMU2KILOGRAM*204.37;
     MolarVolume=0.0000172473;
     volume_vec=31.0721;
     Miedema_Vm=6.6;
     valence_std_vec=3;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=11.85;
     crystal_vec="hcp";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=1.99659E-05;
     lattice_constant[1]=345.66;lattice_constant[2]=345.66;lattice_constant[3]=552.48;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=0.171;
     radius_PT=156;
     radius_covalent_vec=1.45;
     radius_covalent_PT=145;
     radius_VanDerWaals_PT=196;
     RadiiGhosh08=1.8670;
     RadiiSlatter=1.90;
     RadiiPyykko=1.44;
     ElectricalConductivity=6.7E6;
     electronegativity_vec=1.62;
     HardnessGhosh=1.7043;
     ElecNegPearson=3.2;
     ElecNegGhosh=4.719;
     ElectronAffinity_PT=19.2;
     Miedema_phi_star=3.90;
     Miedema_nws=1.12;
     Miedema_gamma_s=610;
     Pettifor_scale=1.56;
     boiling_point=1473;
     melting_point=304;
     VaporizationHeat_PT=165;
     SpecificHeat_PT=129;
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
     ShearModulus=2.8;
     YoungModulus=8;
     BulkModulus=43;
     PoissonRatio_PT=0.45;
     Miedema_BVm=6.2;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-3E-9;
     VolumeMagneticSusceptibility=-0.0000356;
     MolarMagneticSusceptibility=-6.13E-10;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=6500;
     HHIR=6500;
     /*xray_scatt_vec=NNN;*/
     //Tl
 }
 // [AFLOW]STOP=Thallium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Lead
 // Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead
 if(Z==82) {  // Lead
     symbol_vec="Pb";
     name_vec="Lead";
     Period=6;
     Group=14;
     Series="PoorMetal";
     Block="p";
     mass_vec=AMU2KILOGRAM*207.2;
     MolarVolume=0.000018272;
     volume_vec=31.6649;
     Miedema_Vm=6.9;
     valence_std_vec=4;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=11.34;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=1.94378E-05;
     lattice_constant[1]=495.08;lattice_constant[2]=495.08;lattice_constant[3]=495.08;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.175;
     radius_PT=154;
     radius_covalent_vec=1.46;
     radius_covalent_PT=146;
     radius_VanDerWaals_PT=202;
     RadiiGhosh08=1.6523;
     RadiiSlatter=NNN;
     RadiiPyykko=1.44;
     ElectricalConductivity=4.8E6;
     electronegativity_vec=2.33;
     HardnessGhosh=1.9414;
     ElecNegPearson=3.90;
     ElecNegGhosh=4.841;
     ElectronAffinity_PT=35.1;
     Miedema_phi_star=4.10;
     Miedema_nws=1.15;
     Miedema_gamma_s=610;
     Pettifor_scale=1.80;
     boiling_point=1749;
     melting_point=327.46;
     VaporizationHeat_PT=178;
     SpecificHeat_PT=127;
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
     ShearModulus=5.6;
     YoungModulus=16;
     BulkModulus=46;
     PoissonRatio_PT=0.44;
     Miedema_BVm=7.9;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-1.5E-9;
     VolumeMagneticSusceptibility=-0.000017;
     MolarMagneticSusceptibility=-3.11E-10;
     Curie_point=NNN;
     color_PT="SLATEGRAY";
     RefractiveIndex=NNN;
     HHIP=2700;
     HHIR=1800;
     /*xray_scatt_vec=NNN;*/
     //Pb
 }
 // [AFLOW]STOP=Lead
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Bismuth
 // Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth
 if(Z==83) {  // Bismuth
     symbol_vec="Bi";
     name_vec="Bismuth";
     Period=6;
     Group=15;
     Series="PoorMetal";
     Block="p";
     mass_vec=AMU2KILOGRAM*208.9804;
     MolarVolume=0.000021368;
     volume_vec=31.5691;
     Miedema_Vm=7.2;
     valence_std_vec=5;
     valence_iupac_vec=5;
     valence_PT=5;
     Density_PT=9.78;
     crystal_vec="rhl";
     CrystalStr_PT="Base-centered_Monoclinic";
     space_group="C12/m1";
     space_group_number=12;
     Pearson_coefficient=0.0;
     lattice_constant[1]=667.4;lattice_constant[2]=611.7;lattice_constant[3]=330.4;
     lattice_angle[1]=PI/2;lattice_angle[2]=1.925622;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.182;
     radius_PT=143;
     radius_covalent_vec=1.48;
     radius_covalent_PT=148;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.4818;
     RadiiSlatter=1.60;
     RadiiPyykko=1.51;
     ElectricalConductivity=770000;
     electronegativity_vec=2.02;
     HardnessGhosh=2.1785;
     ElecNegPearson=4.69;
     ElecNegGhosh=4.962;
     ElectronAffinity_PT=91.2;
     Miedema_phi_star=4.15;
     Miedema_nws=1.16;
     Miedema_gamma_s=550;
     Pettifor_scale=2.04;
     boiling_point=1564;
     melting_point=271.3;
     VaporizationHeat_PT=160;
     SpecificHeat_PT=122;
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
     ShearModulus=12;
     YoungModulus=32;
     BulkModulus=31;
     PoissonRatio_PT=0.33;
     Miedema_BVm=6.7;
     MagneticType_PT="Diamagnetic";
     MassMagneticSusceptibility=-1.7E-8;
     VolumeMagneticSusceptibility=-0.00017;
     MolarMagneticSusceptibility=-3.6E-9;
     Curie_point=NNN;
     color_PT="GRAY";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Bi
 }
 // [AFLOW]STOP=Bismuth
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Polonium
 // Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium
 if(Z==84) {  // Polonium
     symbol_vec="Po";
     name_vec="Polonium";
     Period=6;
     Group=16;
     Series="Chalcogen";
     Block="p";
     mass_vec=AMU2KILOGRAM*209.98;
     MolarVolume=0.00002272727272727;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=6;
     valence_iupac_vec=6;
     valence_PT=6;
     Density_PT=9.196;
     crystal_vec="sc";
     CrystalStr_PT="Simple_Cubic";
     space_group="Pm-3m";
     space_group_number=221;
     Pearson_coefficient=0.0;
     lattice_constant[1]=335.9;lattice_constant[2]=335.9;lattice_constant[3]=335.9;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.140;
     radius_PT=135;
     radius_covalent_vec=1.40;
     radius_covalent_PT=140;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.3431;
     RadiiSlatter=1.90;
     RadiiPyykko=1.45;
     ElectricalConductivity=2.3E6;
     electronegativity_vec=2.00;
     HardnessGhosh=2.4158;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.084;
     ElectronAffinity_PT=183.3;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=2.28;
     boiling_point=962;
     melting_point=254;
     VaporizationHeat_PT=100;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Po
 }
 // [AFLOW]STOP=Polonium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Astatine
 // Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine
 if(Z==85) {  // Astatine
     symbol_vec="At";
     name_vec="Astatine";
     Period=6;
     Group=17;
     Series="Halogen";
     Block="p";
     mass_vec=AMU2KILOGRAM*210;
     MolarVolume=NNN;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=7;
     valence_iupac_vec=7;
     valence_PT=7;
     Density_PT=NNN;
     crystal_vec="nnn";
     CrystalStr_PT="NNN";
     space_group="NNN";
     space_group_number=NNN;
     Pearson_coefficient=0.0;
     lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
     lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=127;
     radius_covalent_vec=1.50;
     radius_covalent_PT=150;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.2283;
     RadiiSlatter=NNN;
     RadiiPyykko=1.47;
     ElectricalConductivity=NNN;
     electronegativity_vec=2.20;
     HardnessGhosh=2.6528;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.206;
     ElectronAffinity_PT=270.1;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=2.52;
     boiling_point=NNN;
     melting_point=302;
     VaporizationHeat_PT=40;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //At
 }
 // [AFLOW]STOP=Astatine
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Radon
 // Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon
 if(Z==86) {  // Radon
     symbol_vec="Rn";
     name_vec="Radon";
     Period=6;
     Group=18;
     Series="NobleGas";
     Block="p";
     mass_vec=AMU2KILOGRAM*222;
     MolarVolume=0.02281603288798;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=0;
     valence_iupac_vec=6;
     valence_PT=6;
     Density_PT=97.3E-4;
     crystal_vec="fcc";
     CrystalStr_PT="NNN";
     space_group="NNN";
     space_group_number=NNN;
     Pearson_coefficient=0.0;
     lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
     lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
     phase="Gas";
     radius_vec=NNN;
     radius_PT=120;
     radius_covalent_vec=1.50;
     radius_covalent_PT=150;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.1315;
     RadiiSlatter=NNN;
     RadiiPyykko=1.42;
     ElectricalConductivity=NNN;
     electronegativity_vec=2.2;
     HardnessGhosh=2.8900;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.327;
     ElectronAffinity_PT=0;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=-61.7;
     melting_point=-71;
     VaporizationHeat_PT=17;
     SpecificHeat_PT=93.65;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="COLORLESS";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Rn
 }
 // [AFLOW]STOP=Radon
// ********************************************************************************************************************************************************

// ROW7
// s-electron systems
// ********************************************************************************************************************************************************
 // [AFLOW]START=Francium
 // Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium
 if(Z==87) {  // Francium
     symbol_vec="Fr";
     name_vec="Francium";
     Period=7;
     Group=1;
     Series="AlkaliMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*223.02;
     MolarVolume=NNN;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=1;
     valence_iupac_vec=1;
     valence_PT=1;
     Density_PT=NNN;
     crystal_vec="bcc";
     CrystalStr_PT="NNN";
     space_group="NNN";
     space_group_number=NNN;
     Pearson_coefficient=0.0;
     lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
     lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=2.60;
     radius_covalent_PT=260;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=4.4479;
     RadiiSlatter=NNN;
     RadiiPyykko=2.23;
     ElectricalConductivity=NNN;
     electronegativity_vec=0.70;
     HardnessGhosh=0.9882;
     ElecNegPearson=NNN;
     ElecNegGhosh=2.376;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=NNN;
     melting_point=NNN;
     VaporizationHeat_PT=64;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Fr
 }
 // [AFLOW]STOP=Francium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Radium
 // Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium
 if(Z==88) {  // Radium
     symbol_vec="Ra";
     name_vec="Radium";
     Period=7;
     Group=2;
     Series="AlkalineEarthMetal";
     Block="s";
     mass_vec=AMU2KILOGRAM*226.0254;
     MolarVolume=0.0000452;
     volume_vec=-1.0000;
     Miedema_Vm=NNN;
     valence_std_vec=2;
     valence_iupac_vec=2;
     valence_PT=2;
     Density_PT=5;
     crystal_vec="bct";
     CrystalStr_PT="Body-centered_Cubic";
     space_group="Im_3m";
     space_group_number=229;
     Pearson_coefficient=0.0;
     lattice_constant[1]=514.8;lattice_constant[2]=514.8;lattice_constant[3]=514.8;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=2.21;
     radius_covalent_PT=221;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=3.4332;
     RadiiSlatter=2.15;
     RadiiPyykko=2.01;
     ElectricalConductivity=1E6;
     electronegativity_vec=0.89;
     HardnessGhosh=1.2819;
     ElecNegPearson=NNN;
     ElecNegGhosh=2.664;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=1737;
     melting_point=700;
     VaporizationHeat_PT=125;
     SpecificHeat_PT=92;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Ra
 }
 // [AFLOW]STOP=Radium
// ********************************************************************************************************************************************************

// d-electron systems: transition metals
// ********************************************************************************************************************************************************
 // [AFLOW]START=Actinium
 // Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium
 if(Z==89) {  // Actinium
     symbol_vec="Ac";
     name_vec="Actinium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*227.03;
     MolarVolume=0.00002254220456802;
     volume_vec=45.2437;
     Miedema_Vm=NNN;
     valence_std_vec=3;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=10.07;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.0;
     lattice_constant[1]=567;lattice_constant[2]=567;lattice_constant[3]=567;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=2.15;
     radius_covalent_PT=215;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=3.2615;
     RadiiSlatter=1.95;
     RadiiPyykko=1.86;
     ElectricalConductivity=NNN;
     electronegativity_vec=1.10;
     HardnessGhosh=1.3497;
     ElecNegPearson=NNN;
     ElecNegGhosh=2.730;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=3200;
     melting_point=1050;
     VaporizationHeat_PT=400;
     SpecificHeat_PT=120;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Ac
 }
 // [AFLOW]STOP=Actinium
// ********************************************************************************************************************************************************

// actinidies
// ********************************************************************************************************************************************************
 // [AFLOW]START=Thorium
 // Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium
 if(Z==90) {  // Thorium
     symbol_vec="Th";
     name_vec="Thorium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*232.0381;
     MolarVolume=0.0000197917;
     volume_vec=31.9586;
     Miedema_Vm=7.3;
     valence_std_vec=4;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=11.724;
     crystal_vec="fcc";
     CrystalStr_PT="Face-centered_Cubic";
     space_group="Fm_3m";
     space_group_number=225;
     Pearson_coefficient=0.0;
     lattice_constant[1]=508.42;lattice_constant[2]=508.42;lattice_constant[3]=508.42;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.180;
     radius_PT=NNN;
     radius_covalent_vec=2.06;
     radius_covalent_PT=206;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=3.1061;
     RadiiSlatter=1.80;
     RadiiPyykko=1.75;
     ElectricalConductivity=6.7E6;
     electronegativity_vec=1.30;
     HardnessGhosh=1.4175;
     ElecNegPearson=NNN;
     ElecNegGhosh=2.796;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=3.30;
     Miedema_nws=1.28;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=4820;
     melting_point=1750;
     VaporizationHeat_PT=530;
     SpecificHeat_PT=118;
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
     ShearModulus=31;
     YoungModulus=79;
     BulkModulus=54;
     PoissonRatio_PT=0.27;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=7.2E-9;
     VolumeMagneticSusceptibility=0.000084;
     MolarMagneticSusceptibility=1.7E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     xray_scatt_vec=86.64;
     //Th
 }
 // [AFLOW]STOP=Thorium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Protoactinium
 // Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium
 if(Z==91) {  // Protoactinium
     symbol_vec="Pa";
     name_vec="Protoactinium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*231.04;
     MolarVolume=0.0000150316;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=5;
     valence_iupac_vec=5;
     valence_PT=5;
     Density_PT=15.37;
     crystal_vec="bct";
     CrystalStr_PT="Centered_Tetragonal";
     space_group="I4/mmm";
     space_group_number=139;
     Pearson_coefficient=0.0;
     lattice_constant[1]=392.5;lattice_constant[2]=392.5;lattice_constant[3]=323.8;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=2.00;
     radius_covalent_PT=200;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=2.2756;
     RadiiSlatter=1.80;
     RadiiPyykko=1.69;
     ElectricalConductivity=5.6E6;
     electronegativity_vec=1.50;
     HardnessGhosh=1.9369;
     ElecNegPearson=NNN;
     ElecNegGhosh=3.306;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=4000;
     melting_point=1572;
     VaporizationHeat_PT=470;
     SpecificHeat_PT=99.1;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=3.25E-8;
     VolumeMagneticSusceptibility=0.0004995;
     MolarMagneticSusceptibility=7.509E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Pa
 }
 // [AFLOW]STOP=Protoactinium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Uranium
 // Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium
 if(Z==92) {  // Uranium
     symbol_vec="U";
     name_vec="Uranium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*238.03;
     MolarVolume=0.000012495;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=6;
     valence_iupac_vec=6;
     valence_PT=6;
     Density_PT=19.05;
     crystal_vec="orc";
     CrystalStr_PT="Base_Orthorhombic";
     space_group="Cmcm";
     space_group_number=63;
     Pearson_coefficient=1.15611E-06;
     lattice_constant[1]=285.37;lattice_constant[2]=586.95;lattice_constant[3]=495.48;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=0.138;
     radius_PT=NNN;
     radius_covalent_vec=1.96;
     radius_covalent_PT=196;
     radius_VanDerWaals_PT=186;
     RadiiGhosh08=1.9767;
     RadiiSlatter=1.75;
     RadiiPyykko=1.70;
     ElectricalConductivity=3.6E6;
     electronegativity_vec=1.38;
     HardnessGhosh=2.2306;
     ElecNegPearson=NNN;
     ElecNegGhosh=3.594;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=3927;
     melting_point=1135;
     VaporizationHeat_PT=420;
     SpecificHeat_PT=116;
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
     ShearModulus=111;
     YoungModulus=208;
     BulkModulus=100;
     PoissonRatio_PT=0.23;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=2.16E-8;
     VolumeMagneticSusceptibility=0.000411;
     MolarMagneticSusceptibility=5.14E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //U
 }
 // [AFLOW]STOP=Uranium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Neptunium
 // Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium
 if(Z==93) {  // Neptunium
     symbol_vec="Np";
     name_vec="Neptunium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*237.05;
     MolarVolume=0.00001158924205379;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=7;
     valence_iupac_vec=7;
     valence_PT=6;
     Density_PT=20.45;
     crystal_vec="nnn";
     CrystalStr_PT="Simple_Orthorhombic";
     space_group="Pnma";
     space_group_number=62;
     Pearson_coefficient=0.0;
     lattice_constant[1]=666.3;lattice_constant[2]=472.3;lattice_constant[3]=488.7;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=1.90;
     radius_covalent_PT=190;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.7473;
     RadiiSlatter=1.75;
     RadiiPyykko=1.71;
     ElectricalConductivity=830000;
     electronegativity_vec=NNN;
     HardnessGhosh=2.5241;
     ElecNegPearson=NNN;
     ElecNegGhosh=3.882;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=4000;
     melting_point=644;
     VaporizationHeat_PT=335;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Np
 }
 // [AFLOW]STOP=Neptunium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Plutonium
 // Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium
 if(Z==94) {  // Plutonium
     symbol_vec="Pu";
     name_vec="Plutonium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*244.06;
     MolarVolume=0.00001231328219621;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=8;
     valence_iupac_vec=7;
     valence_PT=6;
     Density_PT=19.816;
     crystal_vec="nnn";
     CrystalStr_PT="Simple_Monoclinic";
     space_group="P12_1/m1";
     space_group_number=11;
     Pearson_coefficient=0.0;
     lattice_constant[1]=618.3;lattice_constant[2]=482.2;lattice_constant[3]=1096.3;
     lattice_angle[1]=PI/2;lattice_angle[2]=1.776571;lattice_angle[3]=PI/2;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=1.87;
     radius_covalent_PT=187;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.4496;
     RadiiSlatter=1.75;
     RadiiPyykko=1.72;
     ElectricalConductivity=670000;
     electronegativity_vec=NNN;
     HardnessGhosh=3.0436;
     ElecNegPearson=NNN;
     ElecNegGhosh=4.391;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=3230;
     melting_point=640;
     VaporizationHeat_PT=325;
     SpecificHeat_PT=NNN;
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
     ShearModulus=43;
     YoungModulus=96;
     BulkModulus=NNN;
     PoissonRatio_PT=0.21;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=3.17E-8;
     VolumeMagneticSusceptibility=0.0006282;
     MolarMagneticSusceptibility=7.735E-9;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Pu
 }
 // [AFLOW]STOP=Plutonium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Americium
 // Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium
 if(Z==95) {  // Americium
     symbol_vec="Am";
     name_vec="Americium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*243.06;
     MolarVolume=0.00001777615215801;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=9;
     valence_iupac_vec=7;
     valence_PT=4;
     Density_PT=13.67;
     crystal_vec="nnn";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=346.81;lattice_constant[2]=346.81;lattice_constant[3]=1124.1;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=1.80;
     radius_covalent_PT=180;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.2915;
     RadiiSlatter=1.75;
     RadiiPyykko=1.66;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=3.4169;
     ElecNegPearson=NNN;
     ElecNegGhosh=4.678;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=2011;
     melting_point=1176;
     VaporizationHeat_PT=NNN;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="Paramagnetic";
     MassMagneticSusceptibility=5.15E-8;
     VolumeMagneticSusceptibility=0.000704;
     MolarMagneticSusceptibility=1.251E-8;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Am
 }
 // [AFLOW]STOP=Americium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Curium
 // Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium
 if(Z==96) {  // Curium
     symbol_vec="Cm";
     name_vec="Curium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*247.07;
     MolarVolume=0.00001828275351591;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=10;
     valence_iupac_vec=8;
     valence_PT=4;
     Density_PT=13.51;
     crystal_vec="nnn";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=349.6;lattice_constant[2]=349.6;lattice_constant[3]=1133.1;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=1.69;
     radius_covalent_PT=169;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.2960;
     RadiiSlatter=NNN;
     RadiiPyykko=1.66;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=3.4050;
     ElecNegPearson=NNN;
     ElecNegGhosh=4.745;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=3110;
     melting_point=1345;
     VaporizationHeat_PT=NNN;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="SILVER";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Cm
 }
 // [AFLOW]STOP=Curium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Berkelium
 // Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium
 if(Z==97) {  // Berkelium
     symbol_vec="Bk";
     name_vec="Berkelium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*247.07;
     MolarVolume=0.00001671177266576;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=11;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=14.78;
     crystal_vec="nnn";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=341.6;lattice_constant[2]=341.6;lattice_constant[3]=1106.9;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=NNN;
     radius_covalent_PT=NNN;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.1247;
     RadiiSlatter=NNN;
     RadiiPyykko=1.68;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=3.9244;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.256;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=NNN;
     melting_point=1050;
     VaporizationHeat_PT=NNN;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="NNN";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Bk
 }
 // [AFLOW]STOP=Berkelium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Californium
 // Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium
 if(Z==98) {  // Californium
     symbol_vec="Cf";
     name_vec="Californium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*251.08;
     MolarVolume=0.00001662251655629;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=12;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=15.1;
     crystal_vec="nnn";
     CrystalStr_PT="Simple_Hexagonal";
     space_group="P6_3/mmc";
     space_group_number=194;
     Pearson_coefficient=0.0;
     lattice_constant[1]=338;lattice_constant[2]=338;lattice_constant[3]=1102.5;
     lattice_angle[1]=PI/2;lattice_angle[2]=PI/2;lattice_angle[3]=2*PI/3;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=NNN;
     radius_covalent_PT=NNN;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=1.0465;
     RadiiSlatter=NNN;
     RadiiPyykko=1.68;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=4.2181;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.542;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=NNN;
     melting_point=900;
     VaporizationHeat_PT=NNN;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="NNN";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Cf
 }
 // [AFLOW]STOP=Californium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Einsteinium
 // Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium
 if(Z==99) {  // Einsteinium
     symbol_vec="Es";
     name_vec="Einsteinium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*252.08;
     MolarVolume=NNN;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=13;
     valence_iupac_vec=4;
     valence_PT=4;
     Density_PT=NNN;
     crystal_vec="nnn";
     CrystalStr_PT="NNN";
     space_group="NNN";
     space_group_number=NNN;
     Pearson_coefficient=0.0;
     lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
     lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
     phase="Solid";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=NNN;
     radius_covalent_PT=NNN;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.9785;
     RadiiSlatter=NNN;
     RadiiPyykko=1.65;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=4.5116;
     ElecNegPearson=NNN;
     ElecNegGhosh=5.830;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=NNN;
     melting_point=860;
     VaporizationHeat_PT=NNN;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="NNN";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Es
 }
 // [AFLOW]STOP=Einsteinium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Fermium
 // Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium
 if(Z==100) {  // Fermium
     symbol_vec="Fm";
     name_vec="Fermium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*257.1;
     MolarVolume=NNN;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=14;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=NNN;
     crystal_vec="nnn";
     CrystalStr_PT="NNN";
     space_group="NNN";
     space_group_number=NNN;
     Pearson_coefficient=0.0;
     lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
     lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
     phase="nnn";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=NNN;
     radius_covalent_PT=NNN;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.9188;
     RadiiSlatter=NNN;
     RadiiPyykko=1.67;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=4.8051;
     ElecNegPearson=NNN;
     ElecNegGhosh=6.118;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=NNN;
     melting_point=1527;
     VaporizationHeat_PT=NNN;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="NNN";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Fm
 }
 // [AFLOW]STOP=Fermium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Mendelevium
 // Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium
 if(Z==101) {  // Mendelevium
     symbol_vec="Md";
     name_vec="Mendelevium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*258.1;
     MolarVolume=NNN;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=15;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=NNN;
     crystal_vec="nnn";
     CrystalStr_PT="NNN";
     space_group="NNN";
     space_group_number=NNN;
     Pearson_coefficient=0.0;
     lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
     lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
     phase="nnn";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=NNN;
     radius_covalent_PT=NNN;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.8659;
     RadiiSlatter=NNN;
     RadiiPyykko=1.73;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=5.0990;
     ElecNegPearson=NNN;
     ElecNegGhosh=6.406;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=NNN;
     melting_point=828;
     VaporizationHeat_PT=NNN;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="NNN";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Md
 }
 // [AFLOW]STOP=Mendelevium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Nobelium
 // Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium
 if(Z==102) {  // Nobelium
     symbol_vec="No";
     name_vec="Nobelium";
     Period=7;
     Group=NNN;
     Series="Actinide";
     Block="f";
     mass_vec=AMU2KILOGRAM*259.1;
     MolarVolume=NNN;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=16;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=NNN;
     crystal_vec="nnn";
     CrystalStr_PT="NNN";
     space_group="NNN";
     space_group_number=NNN;
     Pearson_coefficient=0.0;
     lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
     lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
     phase="nnn";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=NNN;
     radius_covalent_PT=NNN;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.8188;
     RadiiSlatter=NNN;
     RadiiPyykko=1.76;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=5.3926;
     ElecNegPearson=NNN;
     ElecNegGhosh=6.694;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=NNN;
     melting_point=828;
     VaporizationHeat_PT=NNN;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="NNN";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //No
 }
 // [AFLOW]STOP=Nobelium
// ********************************************************************************************************************************************************

// ********************************************************************************************************************************************************
 // [AFLOW]START=Lawrencium
 // Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium
 if(Z==103) {  // Lawrencium
     symbol_vec="Lr";
     name_vec="Lawrencium";
     Period=7;
     Group=3;
     Series="TransitionMetal";
     Block="d";
     mass_vec=AMU2KILOGRAM*262.11;
     MolarVolume=NNN;
     volume_vec=NNN;
     Miedema_Vm=NNN;
     valence_std_vec=17;
     valence_iupac_vec=3;
     valence_PT=3;
     Density_PT=NNN;
     crystal_vec="nnn";
     CrystalStr_PT="NNN";
     space_group="NNN";
     space_group_number=NNN;
     Pearson_coefficient=0.0;
     lattice_constant[1]=NNN;lattice_constant[2]=NNN;lattice_constant[3]=NNN;
     lattice_angle[1]=NNN;lattice_angle[2]=NNN;lattice_angle[3]=NNN;
     phase="nnn";
     radius_vec=NNN;
     radius_PT=NNN;
     radius_covalent_vec=NNN;
     radius_covalent_PT=NNN;
     radius_VanDerWaals_PT=NNN;
     RadiiGhosh08=0.8086;
     RadiiSlatter=NNN;
     RadiiPyykko=1.61;
     ElectricalConductivity=NNN;
     electronegativity_vec=NNN;
     HardnessGhosh=5.4607;
     ElecNegPearson=NNN;
     ElecNegGhosh=6.760;
     ElectronAffinity_PT=NNN;
     Miedema_phi_star=NNN;
     Miedema_nws=NNN;
     Miedema_gamma_s=NNN;
     Pettifor_scale=0;
     boiling_point=NNN;
     melting_point=1627;
     VaporizationHeat_PT=NNN;
     SpecificHeat_PT=NNN;
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
     ShearModulus=NNN;
     YoungModulus=NNN;
     BulkModulus=NNN;
     PoissonRatio_PT=NNN;
     Miedema_BVm=NNN;
     MagneticType_PT="NNN";
     MassMagneticSusceptibility=NNN;
     VolumeMagneticSusceptibility=NNN;
     MolarMagneticSusceptibility=NNN;
     Curie_point=NNN;
     color_PT="NNN";
     RefractiveIndex=NNN;
     HHIP=NNN;
     HHIR=NNN;
     /*xray_scatt_vec=NNN;*/
     //Lr
 }
 // [AFLOW]STOP=Lawrencium
 // ********************************************************************************************************************************************************
}

#endif  // _AFLOW_X(*THIS)_CPP

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2019              *
// *                                                                        *
// **************************************************************************
