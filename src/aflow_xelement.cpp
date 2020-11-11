// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2019
#ifndef _AFLOW_XELEMENT_CPP_
#define _AFLOW_XELEMENT_CPP_
#include "aflow.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// _XELEMENT_PROTOTYPES_

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// XELEMENT
// look into aflow.h for the definitions

#define _UTYPE2STRING_PREC_ 12

std::vector<xelement::xelement> velement(NUM_ELEMENTS);        // store starting from ONE

namespace pflow {
  void XelementPrint(const string& options,ostream& oss) {
    bool LDEBUG=0;//(FALSE || XHOST.DEBUG);
    string soliloquy=XPID+"pflow::XelementPrint():";
    if(LDEBUG) cerr << soliloquy << " [BEGIN]" << endl;
    if(LDEBUG) cerr << "options=" << options << endl;
    if(LDEBUG) cerr << "velement.size()=" << velement.size() << endl;

    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(LDEBUG) cerr << "tokens.size()=" << tokens.size() << endl;
    if(tokens.size()==0) {
      init::ErrorOption(options,soliloquy,"aflow --element=Z|name|symbol[,property[,property]....]");
    } 
    // move on
    string species=tokens.at(0);
    uint Z=0; // some defaults
    // try with number
    if(tokens.size()>=1) if(aurostd::string2utype<uint>(species)>0) Z=aurostd::string2utype<uint>(species);
    if(Z>103) {
      init::ErrorOption(options,soliloquy,aurostd::liststring2string("aflow --element=Z|name|symbol[,property[,property]....]","Z outside [1,103] or name or symbol unrecognized"));
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

      //      cerr <<  xelement::xelement(3).name << endl;
      //      cerr <<  xelement::xelement("Li").name << endl;
      //      cerr <<  xelement::xelement("LiThIuM").name << endl;

      // now look at properties
      if(tokens.size()>=2) {
        string c="";
        vector<string> vs;
        uint len=52;
        for(uint i=1;i<tokens.size();i++) {
          c=aurostd::toupper(tokens.at(i));
          vs.clear();
          if(c=="ALL" || c==aurostd::toupper("name")) vs.push_back(aurostd::PaddedPOST("name="+xelement::xelement(Z).name,len));
          if(c=="ALL" || c==aurostd::toupper("symbol")) vs.push_back(aurostd::PaddedPOST("symbol="+xelement::xelement(Z).symbol,len));
          if(c=="ALL" || c==aurostd::toupper("Z")) vs.push_back(aurostd::PaddedPOST("Z="+aurostd::utype2string(xelement::xelement(Z).Z),len));
          if(c=="ALL" || c==aurostd::toupper("period")) vs.push_back(aurostd::PaddedPOST("period="+aurostd::utype2string(xelement::xelement(Z).period),len));
          if(c=="ALL" || c==aurostd::toupper("group")) vs.push_back(aurostd::PaddedPOST("group="+aurostd::utype2string(xelement::xelement(Z).group),len));
          if(c=="ALL" || c==aurostd::toupper("series")) vs.push_back(aurostd::PaddedPOST("series="+xelement::xelement(Z).series,len));
          if(c=="ALL" || c==aurostd::toupper("block")) vs.push_back(aurostd::PaddedPOST("block="+xelement::xelement(Z).block,len));
          //
          if(c=="ALL" || c==aurostd::toupper("mass")) vs.push_back(aurostd::PaddedPOST("mass="+aurostd::utype2string(xelement::xelement(Z).mass,_UTYPE2STRING_PREC_),len)+"// (kg)");
          if(c=="ALL" || c==aurostd::toupper("molar_volume")) vs.push_back(aurostd::PaddedPOST("molar_volume="+aurostd::utype2string(xelement::xelement(Z).molar_volume,_UTYPE2STRING_PREC_),len)+"// (m^3/mol)");
          if(c=="ALL" || c==aurostd::toupper("volume")) vs.push_back(aurostd::PaddedPOST("volume="+aurostd::utype2string(xelement::xelement(Z).volume,_UTYPE2STRING_PREC_),len)+"// A^3");
          if(c=="ALL" || c==aurostd::toupper("Vm_Miedema")) vs.push_back(aurostd::PaddedPOST("Vm_Miedema="+aurostd::utype2string(xelement::xelement(Z).Vm_Miedema,_UTYPE2STRING_PREC_),len)+"// (V_m^{2/3} in (cm^2))");
          //
          if(c=="ALL" || c==aurostd::toupper("valence_std")) vs.push_back(aurostd::PaddedPOST("valence_std="+aurostd::utype2string(xelement::xelement(Z).valence_std),len));
          if(c=="ALL" || c==aurostd::toupper("valence_iupac")) vs.push_back(aurostd::PaddedPOST("valence_iupac="+aurostd::utype2string(xelement::xelement(Z).valence_iupac),len));
          if(c=="ALL" || c==aurostd::toupper("valence_PT")) vs.push_back(aurostd::PaddedPOST("valence_PT="+aurostd::utype2string(xelement::xelement(Z).valence_PT),len));
          if(c=="ALL" || c==aurostd::toupper("density_PT")) vs.push_back(aurostd::PaddedPOST("density_PT="+aurostd::utype2string(xelement::xelement(Z).density_PT,_UTYPE2STRING_PREC_),len)+"// (g/cm^3)");
          if(c=="ALL" || c==aurostd::toupper("crystal")) vs.push_back(aurostd::PaddedPOST("crystal="+xelement::xelement(Z).crystal,len));
          if(c=="ALL" || c==aurostd::toupper("crystal_structure_PT")) vs.push_back(aurostd::PaddedPOST("crystal_structure_PT="+xelement::xelement(Z).crystal_structure_PT,len));
          if(c=="ALL" || c==aurostd::toupper("space_group")) vs.push_back(aurostd::PaddedPOST("space_group="+xelement::xelement(Z).space_group,len));
          if(c=="ALL" || c==aurostd::toupper("space_group_number")) vs.push_back(aurostd::PaddedPOST("space_group_number="+aurostd::utype2string(xelement::xelement(Z).space_group_number),len));
          if(c=="ALL" || c==aurostd::toupper("variance_parameter_mass")) vs.push_back(aurostd::PaddedPOST("variance_parameter_mass="+aurostd::utype2string(xelement::xelement(Z).variance_parameter_mass,_UTYPE2STRING_PREC_),len));
          if(c=="ALL" || c==aurostd::toupper("lattice_constants")) vs.push_back(aurostd::PaddedPOST("lattice_constants="+aurostd::utype2string(xelement::xelement(Z).lattice_constants[1],_UTYPE2STRING_PREC_)+","+aurostd::utype2string(xelement::xelement(Z).lattice_constants[2],_UTYPE2STRING_PREC_)+","+aurostd::utype2string(xelement::xelement(Z).lattice_constants[3],_UTYPE2STRING_PREC_),len)+"// (pm)");
          if(c=="ALL" || c==aurostd::toupper("lattice_angles")) vs.push_back(aurostd::PaddedPOST("lattice_angles="+aurostd::utype2string(xelement::xelement(Z).lattice_angles[1],_UTYPE2STRING_PREC_)+","+aurostd::utype2string(xelement::xelement(Z).lattice_angles[2],_UTYPE2STRING_PREC_)+","+aurostd::utype2string(xelement::xelement(Z).lattice_angles[3],_UTYPE2STRING_PREC_),len)+"// (rad)");
          if(c=="ALL" || c==aurostd::toupper("phase")) vs.push_back(aurostd::PaddedPOST("phase="+xelement::xelement(Z).phase,len));
          if(c=="ALL" || c==aurostd::toupper("radius")) vs.push_back(aurostd::PaddedPOST("radius="+aurostd::utype2string(xelement::xelement(Z).radius,_UTYPE2STRING_PREC_),len)+"// (nm)");
          if(c=="ALL" || c==aurostd::toupper("radius_PT")) vs.push_back(aurostd::PaddedPOST("radius_PT="+aurostd::utype2string(xelement::xelement(Z).radius_PT,_UTYPE2STRING_PREC_),len)+"// (pm)");
          if(c=="ALL" || c==aurostd::toupper("radius_covalent_PT")) vs.push_back(aurostd::PaddedPOST("radius_covalent_PT="+aurostd::utype2string(xelement::xelement(Z).radius_covalent_PT,_UTYPE2STRING_PREC_),len)+"// (pm)");
          if(c=="ALL" || c==aurostd::toupper("radius_covalent")) vs.push_back(aurostd::PaddedPOST("radius_covalent="+aurostd::utype2string(xelement::xelement(Z).radius_covalent,_UTYPE2STRING_PREC_),len)+"// (Angstrom)");
          if(c=="ALL" || c==aurostd::toupper("radius_VanDerWaals_PT")) vs.push_back(aurostd::PaddedPOST("radius_VanDerWaals_PT="+aurostd::utype2string(xelement::xelement(Z).radius_VanDerWaals_PT,_UTYPE2STRING_PREC_),len)+"// (pm)");
          if(c=="ALL" || c==aurostd::toupper("radii_Ghosh08")) vs.push_back(aurostd::PaddedPOST("radii_Ghosh08="+aurostd::utype2string(xelement::xelement(Z).radii_Ghosh08,_UTYPE2STRING_PREC_),len)+"// (Angstrom)");
          if(c=="ALL" || c==aurostd::toupper("radii_Slatter")) vs.push_back(aurostd::PaddedPOST("radii_Slatter="+aurostd::utype2string(xelement::xelement(Z).radii_Slatter,_UTYPE2STRING_PREC_),len)+"// (Angstrom)");
          if(c=="ALL" || c==aurostd::toupper("radii_Pyykko")) vs.push_back(aurostd::PaddedPOST("radii_Pyykko="+aurostd::utype2string(xelement::xelement(Z).radii_Pyykko,_UTYPE2STRING_PREC_),len)+"// (pm)");
          //
          if(c=="ALL" || c==aurostd::toupper("conductivity_electrical")) vs.push_back(aurostd::PaddedPOST("conductivity_electrical="+aurostd::utype2string(xelement::xelement(Z).conductivity_electrical,_UTYPE2STRING_PREC_),len)+"// (S/m)");
          if(c=="ALL" || c==aurostd::toupper("electronegativity_vec")) vs.push_back(aurostd::PaddedPOST("electronegativity_vec="+aurostd::utype2string(xelement::xelement(Z).electronegativity_vec,_UTYPE2STRING_PREC_),len));
          if(c=="ALL" || c==aurostd::toupper("hardness_Ghosh")) vs.push_back(aurostd::PaddedPOST("hardness_Ghosh="+aurostd::utype2string(xelement::xelement(Z).hardness_Ghosh,_UTYPE2STRING_PREC_),len)+"// (eV)");
          if(c=="ALL" || c==aurostd::toupper("electronegativity_Pearson")) vs.push_back(aurostd::PaddedPOST("electronegativity_Pearson="+aurostd::utype2string(xelement::xelement(Z).electronegativity_Pearson,_UTYPE2STRING_PREC_),len)); //+"// (eV)
          if(c=="ALL" || c==aurostd::toupper("electronegativity_Ghosh")) vs.push_back(aurostd::PaddedPOST("electronegativity_Ghosh="+aurostd::utype2string(xelement::xelement(Z).electronegativity_Ghosh,_UTYPE2STRING_PREC_),len)); //+"// (eV)
          if(c=="ALL" || c==aurostd::toupper("electronegativity_Allen")) vs.push_back(aurostd::PaddedPOST("electronegativity_Allen="+aurostd::utype2string(xelement::xelement(Z).electronegativity_Allen,_UTYPE2STRING_PREC_),len)); //CO20200731
          if(c=="ALL" || c==aurostd::toupper("oxidation_states")) vs.push_back(aurostd::PaddedPOST("oxidation_states="+aurostd::joinWDelimiter(aurostd::vecDouble2vecString(xelement::xelement(Z).oxidation_states,_UTYPE2STRING_PREC_),","),len)); //CO20200731
          if(c=="ALL" || c==aurostd::toupper("oxidation_states_preferred")) vs.push_back(aurostd::PaddedPOST("oxidation_states_preferred="+aurostd::joinWDelimiter(aurostd::vecDouble2vecString(xelement::xelement(Z).oxidation_states_preferred,_UTYPE2STRING_PREC_),","),len)); //CO20200731
          if(c=="ALL" || c==aurostd::toupper("electron_affinity_PT")) vs.push_back(aurostd::PaddedPOST("electron_affinity_PT="+aurostd::utype2string(xelement::xelement(Z).electron_affinity_PT,_UTYPE2STRING_PREC_),len)+"// (kJ/mol)");
          if(c=="ALL" || c==aurostd::toupper("phi_star_Miedema")) vs.push_back(aurostd::PaddedPOST("phi_star_Miedema="+aurostd::utype2string(xelement::xelement(Z).phi_star_Miedema,_UTYPE2STRING_PREC_),len)+"// (V) (phi^star)");
          if(c=="ALL" || c==aurostd::toupper("nws_Miedema")) vs.push_back(aurostd::PaddedPOST("nws_Miedema="+aurostd::utype2string(xelement::xelement(Z).nws_Miedema,_UTYPE2STRING_PREC_),len)+"// (d.u.)^1/3 n_{ws}^{1/3}");
          if(c=="ALL" || c==aurostd::toupper("gamma_s_Miedema")) vs.push_back(aurostd::PaddedPOST("gamma_s_Miedema="+aurostd::utype2string(xelement::xelement(Z).gamma_s_Miedema,_UTYPE2STRING_PREC_),len)+"// (mJ/m^2)");
          //
          if(c=="ALL" || c==aurostd::toupper("scale_Pettifor")) vs.push_back(aurostd::PaddedPOST("scale_Pettifor="+aurostd::utype2string(xelement::xelement(Z).scale_Pettifor,_UTYPE2STRING_PREC_),len)); 
          //
          if(c=="ALL" || c==aurostd::toupper("temperature_boiling")) vs.push_back(aurostd::PaddedPOST("temperature_boiling="+aurostd::utype2string(xelement::xelement(Z).temperature_boiling,_UTYPE2STRING_PREC_),len)+"// (Celsius)");
          if(c=="ALL" || c==aurostd::toupper("temperature_melting")) vs.push_back(aurostd::PaddedPOST("temperature_melting="+aurostd::utype2string(xelement::xelement(Z).temperature_melting,_UTYPE2STRING_PREC_),len)+"// (Celsius)");
          if(c=="ALL" || c==aurostd::toupper("vaporization_heat_PT")) vs.push_back(aurostd::PaddedPOST("vaporization_heat_PT="+aurostd::utype2string(xelement::xelement(Z).vaporization_heat_PT,_UTYPE2STRING_PREC_),len)+"// (kJ/mol)");
          if(c=="ALL" || c==aurostd::toupper("specific_heat_PT")) vs.push_back(aurostd::PaddedPOST("specific_heat_PT="+aurostd::utype2string(xelement::xelement(Z).specific_heat_PT,_UTYPE2STRING_PREC_),len)+"// (J/(kg.K))");
          if(c=="ALL" || c==aurostd::toupper("critical_Pressure")) vs.push_back(aurostd::PaddedPOST("critical_Pressure="+aurostd::utype2string(xelement::xelement(Z).critical_Pressure,_UTYPE2STRING_PREC_),len)+"// (Atm) "); 
          if(c=="ALL" || c==aurostd::toupper("critical_Temperature_PT")) vs.push_back(aurostd::PaddedPOST("critical_Temperature_PT="+aurostd::utype2string(xelement::xelement(Z).critical_Temperature_PT,_UTYPE2STRING_PREC_),len)+"// (K)"); 
          if(c=="ALL" || c==aurostd::toupper("thermal_expansion")) vs.push_back(aurostd::PaddedPOST("thermal_expansion="+aurostd::utype2string(xelement::xelement(Z).thermal_expansion,_UTYPE2STRING_PREC_),len)+"// (K^{-1})");
          if(c=="ALL" || c==aurostd::toupper("conductivity_thermal")) vs.push_back(aurostd::PaddedPOST("conductivity_thermal="+aurostd::utype2string(xelement::xelement(Z).conductivity_thermal,_UTYPE2STRING_PREC_),len)+"// (W/(mK))");
          //
          if(c=="ALL" || c==aurostd::toupper("hardness_Brinell")) vs.push_back(aurostd::PaddedPOST("hardness_Brinell="+aurostd::utype2string(xelement::xelement(Z).hardness_Brinell,_UTYPE2STRING_PREC_),len)+"// (MPa)");
          if(c=="ALL" || c==aurostd::toupper("hardness_Mohs")) vs.push_back(aurostd::PaddedPOST("hardness_Mohs="+aurostd::utype2string(xelement::xelement(Z).hardness_Mohs,_UTYPE2STRING_PREC_),len));
          if(c=="ALL" || c==aurostd::toupper("hardness_Vickers")) vs.push_back(aurostd::PaddedPOST("hardness_Vickers="+aurostd::utype2string(xelement::xelement(Z).hardness_Vickers,_UTYPE2STRING_PREC_),len)+"// (MPa)");
          if(c=="ALL" || c==aurostd::toupper("hardness_Pearson")) vs.push_back(aurostd::PaddedPOST("hardness_Pearson="+aurostd::utype2string(xelement::xelement(Z).hardness_Pearson,_UTYPE2STRING_PREC_),len)+"// (eV)");
          if(c=="ALL" || c==aurostd::toupper("hardness_Putz")) vs.push_back(aurostd::PaddedPOST("hardness_Putz="+aurostd::utype2string(xelement::xelement(Z).hardness_Putz,_UTYPE2STRING_PREC_),len)+"// (eV/atom)");
          if(c=="ALL" || c==aurostd::toupper("hardness_RB")) vs.push_back(aurostd::PaddedPOST("hardness_RB="+aurostd::utype2string(xelement::xelement(Z).hardness_RB,_UTYPE2STRING_PREC_),len)+"// (eV)");
          if(c=="ALL" || c==aurostd::toupper("modulus_shear")) vs.push_back(aurostd::PaddedPOST("modulus_shear="+aurostd::utype2string(xelement::xelement(Z).modulus_shear,_UTYPE2STRING_PREC_),len)+"// (GPa)");
          if(c=="ALL" || c==aurostd::toupper("modulus_Young")) vs.push_back(aurostd::PaddedPOST("modulus_Young="+aurostd::utype2string(xelement::xelement(Z).modulus_Young,_UTYPE2STRING_PREC_),len)+"// (GPa)");
          if(c=="ALL" || c==aurostd::toupper("modulus_bulk")) vs.push_back(aurostd::PaddedPOST("modulus_bulk="+aurostd::utype2string(xelement::xelement(Z).modulus_bulk,_UTYPE2STRING_PREC_),len)+"// (GPa)");
          if(c=="ALL" || c==aurostd::toupper("Poisson_ratio_PT")) vs.push_back(aurostd::PaddedPOST("Poisson_ratio_PT="+aurostd::utype2string(xelement::xelement(Z).Poisson_ratio_PT,_UTYPE2STRING_PREC_),len));
          if(c=="ALL" || c==aurostd::toupper("BVm_Miedema")) vs.push_back(aurostd::PaddedPOST("BVm_Miedema="+aurostd::utype2string(xelement::xelement(Z).BVm_Miedema,_UTYPE2STRING_PREC_),len)+"// (kJ/mole)");
          //
          if(c=="ALL" || c==aurostd::toupper("magnetic_type_PT")) vs.push_back(aurostd::PaddedPOST("magnetic_type_PT="+xelement::xelement(Z).magnetic_type_PT,len));
          if(c=="ALL" || c==aurostd::toupper("susceptibility_magnetic_mass")) vs.push_back(aurostd::PaddedPOST("susceptibility_magnetic_mass="+aurostd::utype2string(xelement::xelement(Z).susceptibility_magnetic_mass,_UTYPE2STRING_PREC_),len)+"// (m^3/K)");
          if(c=="ALL" || c==aurostd::toupper("susceptibility_magnetic_volume")) vs.push_back(aurostd::PaddedPOST("susceptibility_magnetic_volume="+aurostd::utype2string(xelement::xelement(Z).susceptibility_magnetic_volume,_UTYPE2STRING_PREC_),len));
          if(c=="ALL" || c==aurostd::toupper("susceptibility_magnetic_molar")) vs.push_back(aurostd::PaddedPOST("susceptibility_magnetic_molar="+aurostd::utype2string(xelement::xelement(Z).susceptibility_magnetic_molar,_UTYPE2STRING_PREC_),len)+"// (m^3/mol)");
          if(c=="ALL" || c==aurostd::toupper("temperature_Curie")) vs.push_back(aurostd::PaddedPOST("temperature_Curie="+aurostd::utype2string(xelement::xelement(Z).temperature_Curie,_UTYPE2STRING_PREC_),len)+"// (K)");
          //
          if(c=="ALL" || c==aurostd::toupper("refractive_index")) vs.push_back(aurostd::PaddedPOST("refractive_index="+aurostd::utype2string(xelement::xelement(Z).refractive_index,_UTYPE2STRING_PREC_),len));
          if(c=="ALL" || c==aurostd::toupper("color_PT")) vs.push_back(aurostd::PaddedPOST("color_PT="+xelement::xelement(Z).color_PT,len));
          //
          if(c=="ALL" || c==aurostd::toupper("HHIP")) vs.push_back(aurostd::PaddedPOST("HHIP="+aurostd::utype2string(xelement::xelement(Z).HHIP,_UTYPE2STRING_PREC_),len));
          if(c=="ALL" || c==aurostd::toupper("HHIR")) vs.push_back(aurostd::PaddedPOST("HHIR="+aurostd::utype2string(xelement::xelement(Z).HHIR,_UTYPE2STRING_PREC_),len));
          if(c=="ALL" || c==aurostd::toupper("xray_scatt")) vs.push_back(aurostd::PaddedPOST("xray_scatt="+aurostd::utype2string(xelement::xelement(Z).xray_scatt,_UTYPE2STRING_PREC_),len)+"// shift+1");

          if(vs.size())
            for(uint j=0;j<vs.size();j++)
              oss << vs.at(j) << endl;
        }
      }
    }

    if(LDEBUG) cerr << soliloquy << " [END]" << endl;
  }
}

//   std::vector<string> vatom_symbol(NUM_ELEMENTS);   // store starting from ONE // DONE
//   std::vector<string> vatom_name(NUM_ELEMENTS);   // store starting from ONE // DONE
//   std::vector<double> vatom_mass(NUM_ELEMENTS);     // store starting from ONE // DONE
//   std::vector<double> vatom_volume(NUM_ELEMENTS);       // store starting from ONE // DONE
//   std::vector<int> vatom_valence_iupac(NUM_ELEMENTS);   // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry) // DONE
//   std::vector<int> vatom_valence_std(NUM_ELEMENTS);     // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry) // DONE
//   std::vector<double> vatom_miedema_phi_star(NUM_ELEMENTS); // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28  
//   std::vector<double> vatom_miedema_nws(NUM_ELEMENTS);      // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
//   std::vector<double> vatom_miedema_Vm(NUM_ELEMENTS);       // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
//   std::vector<double> vatom_miedema_gamma_s(NUM_ELEMENTS);  // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
//   std::vector<double> vatom_miedema_BVm(NUM_ELEMENTS);      // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
//   // for lanthines from J.A. Alonso and N.H. March. Electrons in Metals and Alloys, Academic Press, London (1989) (except La)
//   std::vector<double> vatom_radius(NUM_ELEMENTS);       // store starting from ONE  // DONE
//   std::vector<double> vatom_radius_covalent(NUM_ELEMENTS);// store starting from ONE//DX+CO20170904 
//   std::vector<double> vatom_electronegativity(NUM_ELEMENTS);       // store starting from ONE
//   std::vector<string> vatom_crystal(NUM_ELEMENTS);       // store starting from ONE  // DONE
//   std::vector<double> vatom_xray_scatt(NUM_ELEMENTS);        // store starting from ONE
//   std::vector<double> vatom_pettifor_scale(NUM_ELEMENTS);        // store starting from ONE Chemical Scale Pettifor Solid State Communications 51 31-34 1984
//   std::vector<double> vatom_pearson_coefficient(NUM_ELEMENTS);   //ME20181020 Pearson mass deviation coefficient

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
  xelement::xelement() {free();}  //CO20200520

  // destructor
  xelement::~xelement() {free();} //CO20200520

  void xelement::free() { //CO20200520
    // will populate
    // [AFLOW]START=FREE
    // [AFLOW]STOP=FREE
    // DEFAULT
    verbose=FALSE;
    // [AFLOW]START=CONSTRUCTOR
    Z=0;
    symbol="XX";//"UNDEFINED";
    name="UNDEFINED";
    period=NNN;
    group=NNN; 
    series="UNDEFINED";
    block="nnn";      
    //                                          
    mass=NNN;//  AMU2KILOGRAM goes inside.
    molar_volume=NNN;  
    volume=NNN;      
    Vm_Miedema=NNN;      
    //
    valence_std=NNN;  
    valence_iupac=NNN;
    valence_PT=NNN;       
    density_PT=NNN;       
    crystal="nnn";    
    crystal_structure_PT="UNDEFINED";
    space_group="nnn";     
    space_group_number=NNN;
    variance_parameter_mass=NNN;
    lattice_constants[1]=NNN;lattice_constants[2]=NNN;lattice_constants[3]=NNN;
    lattice_angles[1]=NNN;lattice_angles[2]=NNN;lattice_angles[3]=NNN; 
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
    conductivity_electrical=NNN;
    electronegativity_vec=NNN;    
    hardness_Ghosh=NNN;            
    electronegativity_Pearson=NNN;           
    electronegativity_Ghosh=NNN;             
    electronegativity_Allen=NNN;  //RF+SK20200410
    oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
    oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
    electron_affinity_PT=NNN;      
    phi_star_Miedema=NNN;         
    nws_Miedema=NNN;              
    gamma_s_Miedema=NNN;          
    //
    scale_Pettifor=NNN;          
    //
    temperature_boiling=NNN;         
    temperature_melting=NNN;         
    vaporization_heat_PT=NNN;     
    specific_heat_PT=NNN;         
    critical_Pressure=NNN;     
    critical_Temperature_PT=NNN;  
    thermal_expansion=NNN;     
    conductivity_thermal=NNN;  
    //                                         
    hardness_Brinell=NNN;
    hardness_Mohs=NNN;    
    hardness_Vickers=NNN; 
    hardness_Pearson=NNN;   
    hardness_Putz=NNN;      
    hardness_RB=NNN;        
    modulus_shear=NNN;    
    modulus_Young=NNN;    
    modulus_bulk=NNN;     
    Poisson_ratio_PT=NNN;    
    BVm_Miedema=NNN;        
    //
    magnetic_type_PT="UNDEFINED";     
    susceptibility_magnetic_mass=NNN;
    susceptibility_magnetic_volume=NNN;
    susceptibility_magnetic_molar=NNN; 
    temperature_Curie=NNN;                  
    //
    refractive_index=NNN;             
    color_PT="UNDEFINED";               
    //
    HHIP=NNN;                           
    HHIR=NNN;                           
    xray_scatt=NNN;   
    // [AFLOW]STOP=CONSTRUCTOR
  }

  const xelement& xelement::operator=(const xelement& b) {      // operator=
    if(this!=&b) {free();copy(b);}  //CO20200520
    return *this;
  }

  void xelement::copy(const xelement& b) {  //copy PRIVATE  //CO20200520
    // will populate
    verbose=b.verbose;
    // [AFLOW]START=ASSIGNMENT
    Z=b.Z;
    symbol=b.symbol;
    name=b.name;
    period=b.period;
    group=b.group; 
    series=b.series;
    block=b.block;      
    //                                          
    mass=b.mass;
    molar_volume=b.molar_volume;  
    volume=b.volume;      
    Vm_Miedema=b.Vm_Miedema;      
    //
    valence_std=b.valence_std;  
    valence_iupac=b.valence_iupac;
    valence_PT=b.valence_PT;       
    density_PT=b.density_PT;       
    crystal=b.crystal;    
    crystal_structure_PT=b.crystal_structure_PT;
    space_group=b.space_group;
    space_group_number=b.space_group_number;    
    variance_parameter_mass=b.variance_parameter_mass;
    lattice_constants=b.lattice_constants; 
    lattice_angles=b.lattice_angles;   
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
    conductivity_electrical=b.conductivity_electrical;
    electronegativity_vec=b.electronegativity_vec;    
    hardness_Ghosh=b.hardness_Ghosh;            
    electronegativity_Pearson=b.electronegativity_Pearson;           
    electronegativity_Ghosh=b.electronegativity_Ghosh;             
    electronegativity_Allen=b.electronegativity_Allen;  //RF+SK20200410
    oxidation_states=b.oxidation_states;  //RF+SK20200410
    oxidation_states_preferred=b.oxidation_states_preferred;  //RF+SK20200410
    electron_affinity_PT=b.electron_affinity_PT;      
    phi_star_Miedema=b.phi_star_Miedema;         
    nws_Miedema=b.nws_Miedema;              
    gamma_s_Miedema=b.gamma_s_Miedema;          
    //
    scale_Pettifor=b.scale_Pettifor;          
    //
    temperature_boiling=b.temperature_boiling;         
    temperature_melting=b.temperature_melting;         
    vaporization_heat_PT=b.vaporization_heat_PT;     
    specific_heat_PT=b.specific_heat_PT;         
    critical_Pressure=b.critical_Pressure;     
    critical_Temperature_PT=b.critical_Temperature_PT;  
    thermal_expansion=b.thermal_expansion;     
    conductivity_thermal=b.conductivity_thermal;  
    //                                         
    hardness_Brinell=b.hardness_Brinell;
    hardness_Mohs=b.hardness_Mohs;    
    hardness_Vickers=b.hardness_Vickers; 
    hardness_Pearson=b.hardness_Pearson;   
    hardness_Putz=b.hardness_Putz;      
    hardness_RB=b.hardness_RB;        
    modulus_shear=b.modulus_shear;    
    modulus_Young=b.modulus_Young;    
    modulus_bulk=b.modulus_bulk;     
    Poisson_ratio_PT=b.Poisson_ratio_PT;    
    BVm_Miedema=b.BVm_Miedema;        
    //
    magnetic_type_PT=b.magnetic_type_PT;
    susceptibility_magnetic_mass=b.susceptibility_magnetic_mass;
    susceptibility_magnetic_volume=b.susceptibility_magnetic_volume;
    susceptibility_magnetic_molar=b.susceptibility_magnetic_molar; 
    temperature_Curie=b.temperature_Curie;                  
    //
    refractive_index=b.refractive_index;             
    color_PT=b.color_PT;         
    //
    HHIP=b.HHIP;                           
    HHIR=b.HHIR;                           
    xray_scatt=b.xray_scatt;    
    // [AFLOW]STOP=ASSIGNMENT
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

  xelement::xelement(const string& element) {free();populate(element);}  //CO20200520
  xelement::xelement(uint ZZ) {free();populate(ZZ);} //CO20200520

  string xelement::getProperty(const string& property,const string& delim) const { //CO20201111
    if(property==aurostd::toupper("name")) return name;
    if(property==aurostd::toupper("symbol")) return symbol;
    if(property==aurostd::toupper("Z")) return aurostd::utype2string(Z);
    if(property==aurostd::toupper("period")) return aurostd::utype2string(period);
    if(property==aurostd::toupper("group")) return aurostd::utype2string(group);
    if(property==aurostd::toupper("series")) return series;
    if(property==aurostd::toupper("block")) return block;
    //
    if(property==aurostd::toupper("mass")) return aurostd::utype2string(mass,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("molar_volume")) return aurostd::utype2string(molar_volume,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("volume")) return aurostd::utype2string(volume,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("Vm_Miedema")) return aurostd::utype2string(Vm_Miedema,_UTYPE2STRING_PREC_);
    //
    if(property==aurostd::toupper("valence_std")) return aurostd::utype2string(valence_std);
    if(property==aurostd::toupper("valence_iupac")) return aurostd::utype2string(valence_iupac);
    if(property==aurostd::toupper("valence_PT")) return aurostd::utype2string(valence_PT);
    if(property==aurostd::toupper("density_PT")) return aurostd::utype2string(density_PT,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("crystal")) return crystal;
    if(property==aurostd::toupper("crystal_structure_PT")) return crystal_structure_PT;
    if(property==aurostd::toupper("space_group")) return space_group;
    if(property==aurostd::toupper("space_group_number")) return aurostd::utype2string(space_group_number);
    if(property==aurostd::toupper("variance_parameter_mass")) return aurostd::utype2string(variance_parameter_mass,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("lattice_constants")) return aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(lattice_constants,_UTYPE2STRING_PREC_),delim);
    if(property==aurostd::toupper("lattice_angles")) return aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(lattice_angles,_UTYPE2STRING_PREC_),delim);
    if(property==aurostd::toupper("phase")) return phase;
    if(property==aurostd::toupper("radius")) return aurostd::utype2string(radius,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("radius_PT")) return aurostd::utype2string(radius_PT,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("radius_covalent_PT")) return aurostd::utype2string(radius_covalent_PT,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("radius_covalent")) return aurostd::utype2string(radius_covalent,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("radius_VanDerWaals_PT")) return aurostd::utype2string(radius_VanDerWaals_PT,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("radii_Ghosh08")) return aurostd::utype2string(radii_Ghosh08,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("radii_Slatter")) return aurostd::utype2string(radii_Slatter,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("radii_Pyykko")) return aurostd::utype2string(radii_Pyykko,_UTYPE2STRING_PREC_);
    //
    if(property==aurostd::toupper("conductivity_electrical")) return aurostd::utype2string(conductivity_electrical,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("electronegativity_vec")) return aurostd::utype2string(electronegativity_vec,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("hardness_Ghosh")) return aurostd::utype2string(hardness_Ghosh,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("electronegativity_Pearson")) return aurostd::utype2string(electronegativity_Pearson,_UTYPE2STRING_PREC_); //+"// (eV)
    if(property==aurostd::toupper("electronegativity_Ghosh")) return aurostd::utype2string(electronegativity_Ghosh,_UTYPE2STRING_PREC_); //+"// (eV)
    if(property==aurostd::toupper("electronegativity_Allen")) return aurostd::utype2string(electronegativity_Allen,_UTYPE2STRING_PREC_); //CO20200731
    if(property==aurostd::toupper("oxidation_states")) return aurostd::joinWDelimiter(aurostd::vecDouble2vecString(oxidation_states,_UTYPE2STRING_PREC_),delim); //CO20200731
    if(property==aurostd::toupper("oxidation_states_preferred")) return aurostd::joinWDelimiter(aurostd::vecDouble2vecString(oxidation_states_preferred,_UTYPE2STRING_PREC_),delim); //CO20200731
    if(property==aurostd::toupper("electron_affinity_PT")) return aurostd::utype2string(electron_affinity_PT,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("phi_star_Miedema")) return aurostd::utype2string(phi_star_Miedema,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("nws_Miedema")) return aurostd::utype2string(nws_Miedema,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("gamma_s_Miedema")) return aurostd::utype2string(gamma_s_Miedema,_UTYPE2STRING_PREC_);
    //
    if(property==aurostd::toupper("scale_Pettifor")) return aurostd::utype2string(scale_Pettifor,_UTYPE2STRING_PREC_); 
    //
    if(property==aurostd::toupper("temperature_boiling")) return aurostd::utype2string(temperature_boiling,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("temperature_melting")) return aurostd::utype2string(temperature_melting,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("vaporization_heat_PT")) return aurostd::utype2string(vaporization_heat_PT,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("specific_heat_PT")) return aurostd::utype2string(specific_heat_PT,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("critical_Pressure")) return aurostd::utype2string(critical_Pressure,_UTYPE2STRING_PREC_); 
    if(property==aurostd::toupper("critical_Temperature_PT")) return aurostd::utype2string(critical_Temperature_PT,_UTYPE2STRING_PREC_); 
    if(property==aurostd::toupper("thermal_expansion")) return aurostd::utype2string(thermal_expansion,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("conductivity_thermal")) return aurostd::utype2string(conductivity_thermal,_UTYPE2STRING_PREC_);
    //
    if(property==aurostd::toupper("hardness_Brinell")) return aurostd::utype2string(hardness_Brinell,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("hardness_Mohs")) return aurostd::utype2string(hardness_Mohs,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("hardness_Vickers")) return aurostd::utype2string(hardness_Vickers,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("hardness_Pearson")) return aurostd::utype2string(hardness_Pearson,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("hardness_Putz")) return aurostd::utype2string(hardness_Putz,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("hardness_RB")) return aurostd::utype2string(hardness_RB,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("modulus_shear")) return aurostd::utype2string(modulus_shear,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("modulus_Young")) return aurostd::utype2string(modulus_Young,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("modulus_bulk")) return aurostd::utype2string(modulus_bulk,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("Poisson_ratio_PT")) return aurostd::utype2string(Poisson_ratio_PT,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("BVm_Miedema")) return aurostd::utype2string(BVm_Miedema,_UTYPE2STRING_PREC_);
    //
    if(property==aurostd::toupper("magnetic_type_PT")) return magnetic_type_PT;
    if(property==aurostd::toupper("susceptibility_magnetic_mass")) return aurostd::utype2string(susceptibility_magnetic_mass,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("susceptibility_magnetic_volume")) return aurostd::utype2string(susceptibility_magnetic_volume,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("susceptibility_magnetic_molar")) return aurostd::utype2string(susceptibility_magnetic_molar,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("temperature_Curie")) return aurostd::utype2string(temperature_Curie,_UTYPE2STRING_PREC_);
    //
    if(property==aurostd::toupper("refractive_index")) return aurostd::utype2string(refractive_index,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("color_PT")) return color_PT;
    //
    if(property==aurostd::toupper("HHIP")) return aurostd::utype2string(HHIP,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("HHIR")) return aurostd::utype2string(HHIR,_UTYPE2STRING_PREC_);
    if(property==aurostd::toupper("xray_scatt")) return aurostd::utype2string(xray_scatt,_UTYPE2STRING_PREC_);
    return "";
  }

  // ********************************************************************************************************************************************************
  // populate by name or symbol
  void xelement::populate(const string& element) {  //CO20200520
    free();
    // DEFAULT
    verbose=FALSE;
    uint Z=0;

    // try with symbol
    if(Z==0) {
      for(uint i=1;i<=103&&Z==0;i++)  //CO20200520
        if(aurostd::toupper(element)==aurostd::toupper(xelement(i).symbol)) Z=i;
    }
    // try with name
    if(Z==0) {
      for(uint i=1;i<=103&&Z==0;i++)  //CO20200520
        if(aurostd::toupper(element)==aurostd::toupper(xelement(i).name)) Z=i;
    }
    if(Z!=0) {(*this)=xelement(Z);return;}  //CO20200520

    throw aurostd::xerror(_AFLOW_FILE_NAME_,"xelement::xelement():","Element symbol/name does not exist: "+element,_VALUE_ILLEGAL_); //CO20200520
  }

  // ********************************************************************************************************************************************************
  // populate by Z
  void xelement::populate(uint ZZ) {  //CO20200520
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
      return; //CO20200520
    }
    // ROW 1
    // s-electron systems

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Hydrogen
    // Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen Hydrogen
    else if(ZZ==1) { // Hydrogen
      Z=ZZ;
      symbol="H";
      name="Hydrogen";
      period=1;
      group=1;
      series="Nonmetal";
      block="s";
      mass=AMU2KILOGRAM*1.0079;
      molar_volume=0.01121;
      volume=0.75110;
      Vm_Miedema=NNN;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      density_PT=0.899E-4;
      crystal="hex";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.00011460743;
      lattice_constants[1]=470;lattice_constants[2]=470;lattice_constants[3]=340;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Gas";
      radius=0.046;
      radius_PT=53;
      radius_covalent=0.31;
      radius_covalent_PT=31;
      radius_VanDerWaals_PT=120;
      radii_Ghosh08=0.5292;
      radii_Slatter=0.25;
      radii_Pyykko=0.32;
      conductivity_electrical=NNN;
      electronegativity_vec=2.10;
      hardness_Ghosh=6.4299;
      electronegativity_Pearson=7.18;
      electronegativity_Ghosh=7.178;
      electronegativity_Allen=2.300;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(1); oxidation_states.push_back(-1);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(1);  //RF+SK20200410
      electron_affinity_PT=72.8;
      phi_star_Miedema=5.2;
      nws_Miedema=1.5;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=-252.87;
      temperature_melting=-259.14;
      vaporization_heat_PT=0.452;
      specific_heat_PT=14300;
      critical_Pressure=12.76;
      critical_Temperature_PT=32.97;
      thermal_expansion=NNN;
      conductivity_thermal=0.1805;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=6.43;
      hardness_Putz=6.45;
      hardness_RB=6.83;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-2.48E-8;
      susceptibility_magnetic_volume=-2.23E-9;
      susceptibility_magnetic_molar=-4.999E-11;
      temperature_Curie=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000132;
      HHIP=NNN;
      HHIR=NNN;
      xray_scatt=1.000;
      // H volume wrong *dimer* MIEDEMA =PAUL VAN DER PUT book
      return; //CO20200520
    }
    // [AFLOW]STOP=Hydrogen
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Helium
    // Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium Helium
    else if(ZZ==2) { // Helium
      Z=ZZ;
      symbol="He";
      name="Helium";
      period=1;
      group=18;
      series="NobleGas";
      block="s";
      mass=AMU2KILOGRAM*4.0026;
      molar_volume=0.022424;
      volume=-1.000;
      Vm_Miedema=NNN;
      valence_std=0;
      valence_iupac=0;
      valence_PT=0;
      density_PT=1.785E-4;
      crystal="hcp";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=8.32328E-8;
      lattice_constants[1]=424.2;lattice_constants[2]=424.2;lattice_constants[3]=424.2;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Gas";
      radius=NNN;
      radius_PT=31;
      radius_covalent=0.28;
      radius_covalent_PT=28;
      radius_VanDerWaals_PT=140;
      radii_Ghosh08=0.3113;
      radii_Slatter=NNN;
      radii_Pyykko=0.46;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=12.5449;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=12.046;
      electronegativity_Allen=4.160;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=-268.93;
      temperature_melting=NNN;
      vaporization_heat_PT=0.083;
      specific_heat_PT=5193.1;
      critical_Pressure=2.24;
      critical_Temperature_PT=5.19;
      thermal_expansion=NNN;
      conductivity_thermal=0.1513;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=25.79;
      hardness_RB=16.88;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-5.9E-9;
      susceptibility_magnetic_volume=-1.05E-9;
      susceptibility_magnetic_molar=-2.36E-11;
      temperature_Curie=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000035;
      HHIP=3200;
      HHIR=3900;
      xray_scatt=2.000;
      // He
      return; //CO20200520
    }
    // [AFLOW]STOP=Helium
    // ********************************************************************************************************************************************************

    // ROW2
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Lithium
    // Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium Lithium
    else if(ZZ==3) { // Lithium
      Z=ZZ;
      symbol="Li";
      name="Lithium";
      period=2;
      group=1;
      series="AlkaliMetal";
      block="s";
      mass=AMU2KILOGRAM*6.941;
      molar_volume=0.00001297;
      volume=20.24110;
      Vm_Miedema=5.5;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      density_PT=0.535;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=0.0014588232;
      lattice_constants[1]=351;lattice_constants[2]=351;lattice_constants[3]=351;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.152;
      radius_PT=167;
      radius_covalent=1.28;
      radius_covalent_PT=128;
      radius_VanDerWaals_PT=182;
      radii_Ghosh08=1.6283;
      radii_Slatter=1.45;
      radii_Pyykko=1.33;
      conductivity_electrical=1.1E7;
      electronegativity_vec=0.98;
      hardness_Ghosh=2.3746;
      electronegativity_Pearson=3.01;
      electronegativity_Ghosh=2.860;
      electronegativity_Allen=0.912;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(1);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(1);  //RF+SK20200410
      electron_affinity_PT=59.6;
      phi_star_Miedema=2.85;
      nws_Miedema=0.98;
      gamma_s_Miedema=530;
      scale_Pettifor=0.45;
      temperature_boiling=1342;
      temperature_melting=180.54;
      vaporization_heat_PT=147;
      specific_heat_PT=3570;
      critical_Pressure=661.2;
      critical_Temperature_PT=3223;
      thermal_expansion=0.000046;
      conductivity_thermal=85;
      hardness_Brinell=NNN;
      hardness_Mohs=0.6;
      hardness_Vickers=NNN;
      hardness_Pearson=2.39;
      hardness_Putz=0.65;
      hardness_RB=3.06;
      modulus_shear=4.2;
      modulus_Young=4.9;
      modulus_bulk=11;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=1.5;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=2.56E-8;
      susceptibility_magnetic_volume=0.0000137;
      susceptibility_magnetic_molar=1.78E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=2900;
      HHIR=4200;
      xray_scatt=3.00145;
      // Li
      return; //CO20200520
    }
    // [AFLOW]STOP=Lithium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Beryllium
    // Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium Beryllium
    else if(ZZ==4) { // Beryllium
      Z=ZZ;
      symbol="Be";
      name="Beryllium";
      period=2;
      group=2;
      series="AlkalineEarthMetal";
      block="s";
      mass=AMU2KILOGRAM*9.0122;
      molar_volume=4.8767E-6;
      volume=7.83290;
      Vm_Miedema=2.9;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      density_PT=1.848;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=228.58;lattice_constants[2]=228.58;lattice_constants[3]=358.43;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.114;
      radius_PT=112;
      radius_covalent=0.96;
      radius_covalent_PT=96;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.0855;
      radii_Slatter=1.05;
      radii_Pyykko=1.02;
      conductivity_electrical=2.5E7;
      electronegativity_vec=1.57;
      hardness_Ghosh=3.4968;
      electronegativity_Pearson=4.90;
      electronegativity_Ghosh=3.945;
      electronegativity_Allen=1.576;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=4.20;
      nws_Miedema=1.60;
      gamma_s_Miedema=1900;
      scale_Pettifor=1.50;
      temperature_boiling=2470;
      temperature_melting=1287;
      vaporization_heat_PT=297;
      specific_heat_PT=1820;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000113;
      conductivity_thermal=190;
      hardness_Brinell=600;
      hardness_Mohs=5.5;
      hardness_Vickers=1670;
      hardness_Pearson=4.50;
      hardness_Putz=1.69;
      hardness_RB=5.16;
      modulus_shear=132;
      modulus_Young=287;
      modulus_bulk=130;
      Poisson_ratio_PT=0.032;
      BVm_Miedema=4.9;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-1.26E-8;
      susceptibility_magnetic_volume=-0.00002328;
      susceptibility_magnetic_molar=-1.136E-10;
      temperature_Curie=NNN;
      color_PT="SLATEGRAY";
      refractive_index=NNN;
      HHIP=8000;
      HHIR=4000;
      /*xray_scatt=NNN;*/
      // Be
      return; //CO20200520
    }
    // [AFLOW]STOP=Beryllium
    // ********************************************************************************************************************************************************

    // p-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Boron
    // Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron Boron
    else if(ZZ==5) { // Boron
      Z=ZZ;
      symbol="B";
      name="Boron";
      period=2;
      group=13;
      series="Metalloid";
      block="p";
      mass=AMU2KILOGRAM*10.81;
      molar_volume=4.3943E-6;
      volume=5.88420;
      Vm_Miedema=2.8;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      density_PT=2.46;
      crystal="tet";
      crystal_structure_PT="Simple_Trigonal";
      space_group="R_3m";
      space_group_number=166;
      variance_parameter_mass=0.00135391428;
      lattice_constants[1]=506;lattice_constants[2]=506;lattice_constants[3]=506;
      lattice_angles[1]=1.01334;lattice_angles[2]=1.01334;lattice_angles[3]=1.01334;
      phase="Solid";
      radius=0.097;
      radius_PT=87;
      radius_covalent=0.84;
      radius_covalent_PT=85;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8141;
      radii_Slatter=0.85;
      radii_Pyykko=0.85;
      conductivity_electrical=0.0001;
      electronegativity_vec=2.04;
      hardness_Ghosh=4.6190;
      electronegativity_Pearson=4.29;
      electronegativity_Ghosh=5.031;
      electronegativity_Allen=2.051;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(3);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);  //RF+SK20200410
      electron_affinity_PT=26.7;
      phi_star_Miedema=4.75;
      nws_Miedema=1.55;
      gamma_s_Miedema=NNN;
      scale_Pettifor=2.00;
      temperature_boiling=4000;
      temperature_melting=2075;
      vaporization_heat_PT=507;
      specific_heat_PT=1030;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6E-6;
      conductivity_thermal=27;
      hardness_Brinell=NNN;
      hardness_Mohs=9.3;
      hardness_Vickers=49000;
      hardness_Pearson=4.01;
      hardness_Putz=3.46;
      hardness_RB=4.39;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=320;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-8.7E-9;
      susceptibility_magnetic_volume=-0.0000214;
      susceptibility_magnetic_molar=-9.41E-11;
      temperature_Curie=NNN;
      color_PT="BLACK";
      refractive_index=NNN;
      HHIP=2900;
      HHIR=2000;
      /*xray_scatt=NNN;*/
      // B
      return; //CO20200520
    }
    // [AFLOW]STOP=Boron
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Carbon
    // Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon Carbon
    else if(ZZ==6) { // Carbon
      Z=ZZ;
      symbol="C";
      name="Carbon";
      period=2;
      group=14;
      series="Nonmetal";
      block="p";
      mass=AMU2KILOGRAM*12.011;
      molar_volume=5.3146E-6;
      volume=5.59490;
      Vm_Miedema=1.8;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      density_PT=2.26;
      crystal="dia";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.00007387218;
      lattice_constants[1]=246.4;lattice_constants[2]=246.4;lattice_constants[3]=671.1;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.077;
      radius_PT=67;
      radius_covalent=0.76;
      radius_covalent_PT=76;
      radius_VanDerWaals_PT=170;
      radii_Ghosh08=0.6513;
      radii_Slatter=0.70;
      radii_Pyykko=0.75;
      conductivity_electrical=100000;
      electronegativity_vec=2.55;
      hardness_Ghosh=5.7410;
      electronegativity_Pearson=6.27;
      electronegativity_Ghosh=6.116;
      electronegativity_Allen=2.544;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(2); oxidation_states.push_back(-4); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4); oxidation_states_preferred.push_back(-4);  //RF+SK20200410
      electron_affinity_PT=153.9;
      phi_star_Miedema=6.20;
      nws_Miedema=1.90;
      gamma_s_Miedema=NNN;
      scale_Pettifor=2.50;
      temperature_boiling=4027;
      temperature_melting=3550;
      vaporization_heat_PT=715;
      specific_heat_PT=710;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=7.1E-6;
      conductivity_thermal=140;
      hardness_Brinell=NNN;
      hardness_Mohs=0.5;
      hardness_Vickers=NNN;
      hardness_Pearson=5.00;
      hardness_Putz=6.21;
      hardness_RB=5.49;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=33;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-6.2E-9;
      susceptibility_magnetic_volume=-0.000014;
      susceptibility_magnetic_molar=-7.45E-11;
      temperature_Curie=NNN;
      color_PT="BLACK";
      refractive_index=2.417;
      HHIP=500;
      HHIR=500;
      xray_scatt=6.019;
      // C//DX+CO20170904 radius_covalent uses sp3 hybridization (most common)
      return; //CO20200520
    }
    // [AFLOW]STOP=Carbon
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Nitrogen
    // Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen Nitrogen
    else if(ZZ==7) { // Nitrogen
      Z=ZZ;
      symbol="N";
      name="Nitrogen";
      period=2;
      group=15;
      series="Nonmetal";
      block="p";
      mass=AMU2KILOGRAM*14.0067;
      molar_volume=0.011197;
      volume=7.59940;
      Vm_Miedema=2.2;
      valence_std=5;
      valence_iupac=5;
      valence_PT=3;
      density_PT=12.51E-4;
      crystal="hex";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.00001857771;
      lattice_constants[1]=386.1;lattice_constants[2]=386.1;lattice_constants[3]=626.5;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Gas";
      radius=0.071;
      radius_PT=56;
      radius_covalent=0.71;
      radius_covalent_PT=71;
      radius_VanDerWaals_PT=155;
      radii_Ghosh08=0.5428;
      radii_Slatter=0.65;
      radii_Pyykko=0.71;
      conductivity_electrical=NNN;
      electronegativity_vec=3.04;
      hardness_Ghosh=6.8625;
      electronegativity_Pearson=7.30;
      electronegativity_Ghosh=7.209;
      electronegativity_Allen=3.066;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(5); oxidation_states.push_back(4); oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(-3); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(-3); //RF+SK20200410
      electron_affinity_PT=7;
      phi_star_Miedema=7.00;
      nws_Miedema=1.60;
      gamma_s_Miedema=NNN;
      scale_Pettifor=3.00;
      temperature_boiling=-195.79;
      temperature_melting=-210.1;
      vaporization_heat_PT=2.79;
      specific_heat_PT=1040;
      critical_Pressure=33.46;
      critical_Temperature_PT=126.21;
      thermal_expansion=NNN;
      conductivity_thermal=0.02583;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=7.23;
      hardness_Putz=9.59;
      hardness_RB=8.59;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-5.4E-9;
      susceptibility_magnetic_volume=-6.8E-9;
      susceptibility_magnetic_molar=-1.5E-10;
      temperature_Curie=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000298;
      HHIP=1300;
      HHIR=500;
      /*xray_scatt=NNN;*/
      //N JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Nitrogen
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Oxygen
    // Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen Oxygen
    else if(ZZ==8) { // Oxygen
      Z=ZZ;
      symbol="O";
      name="Oxygen";
      period=2;
      group=16;
      series="Chalcogen";
      block="p";
      mass=AMU2KILOGRAM*15.9994;
      molar_volume=0.011196;
      volume=7.78230;
      Vm_Miedema=2.656;
      valence_std=6;
      valence_iupac=2;
      valence_PT=2;
      density_PT=14.29E-4;
      crystal="cub";
      crystal_structure_PT="Base-centered_Monoclinic";
      space_group="C12/m1";
      space_group_number=12;
      variance_parameter_mass=0.00003358805;
      lattice_constants[1]=540.3;lattice_constants[2]=342.9;lattice_constants[3]=508.6;
      lattice_angles[1]=PI/2;lattice_angles[2]=2.313085;lattice_angles[3]=PI/2;
      phase="Gas";
      radius=0.060;
      radius_PT=48;
      radius_covalent=0.66;
      radius_covalent_PT=66;
      radius_VanDerWaals_PT=152;
      radii_Ghosh08=0.4652;
      radii_Slatter=0.60;
      radii_Pyykko=0.63;
      conductivity_electrical=NNN;
      electronegativity_vec=3.44;
      hardness_Ghosh=7.9854;
      electronegativity_Pearson=7.54;
      electronegativity_Ghosh=8.287;
      electronegativity_Allen=3.610;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(-0.5); oxidation_states.push_back(-1); oxidation_states.push_back(-2); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(-2); //RF+SK20200410
      electron_affinity_PT=141;
      phi_star_Miedema=6.97;
      nws_Miedema=1.70;
      gamma_s_Miedema=NNN;
      scale_Pettifor=3.50;
      temperature_boiling=-182.9;
      temperature_melting=-218.3;
      vaporization_heat_PT=3.41;
      specific_heat_PT=919;
      critical_Pressure=49.77;
      critical_Temperature_PT=154.59;
      thermal_expansion=NNN;
      conductivity_thermal=0.02658;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=6.08;
      hardness_Putz=13.27;
      hardness_RB=6.42;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.335E-6;
      susceptibility_magnetic_volume=1.90772E-6;
      susceptibility_magnetic_molar=4.27184E-8;
      temperature_Curie=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000271;
      HHIP=500;
      HHIR=500;
      xray_scatt=8.052;
      // O Table 27 of JX
      return; //CO20200520
    }
    // [AFLOW]STOP=Oxygen
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Fluorine
    // Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine Fluorine
    else if(ZZ==9) { // Fluorine
      Z=ZZ;
      symbol="F";
      name="Fluorine";
      period=2;
      group=17;
      series="Halogen";
      block="p";
      mass=AMU2KILOGRAM*18.9984;
      molar_volume=0.011202;
      volume=9.99090;
      Vm_Miedema=NNN;
      valence_std=7;
      valence_iupac=1;
      valence_PT=1;
      density_PT=16.96E-4;
      crystal="mcl";
      crystal_structure_PT="Base-centered_Monoclinic";
      space_group="C12/c1";
      space_group_number=15;
      variance_parameter_mass=0.0;
      lattice_constants[1]=550;lattice_constants[2]=328;lattice_constants[3]=728;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Gas";
      radius=NNN;
      radius_PT=42;
      radius_covalent=0.57;
      radius_covalent_PT=57;
      radius_VanDerWaals_PT=147;
      radii_Ghosh08=0.4071;
      radii_Slatter=0.50;
      radii_Pyykko=0.64;
      conductivity_electrical=NNN;
      electronegativity_vec=3.98;
      hardness_Ghosh=9.1065;
      electronegativity_Pearson=10.41;
      electronegativity_Ghosh=9.372;
      electronegativity_Allen=4.193;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(-1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(-1); //RF+SK20200410
      electron_affinity_PT=328;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=4.00;
      temperature_boiling=-188.12;
      temperature_melting=-219.6;
      vaporization_heat_PT=3.27;
      specific_heat_PT=824;
      critical_Pressure=51.04;
      critical_Temperature_PT=144.13;
      thermal_expansion=NNN;
      conductivity_thermal=0.0277;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=7.01;
      hardness_Putz=16.16;
      hardness_RB=7.52;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000195;
      HHIP=1500;
      HHIR=1500;
      /*xray_scatt=NNN;*/
      //F
      return; //CO20200520
    }
    // [AFLOW]STOP=Fluorine
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Neon
    // Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon Neon
    else if(ZZ==10) { // Neon
      Z=ZZ;
      symbol="Ne";
      name="Neon";
      period=2;
      group=18;
      series="NobleGas";
      block="p";
      mass=AMU2KILOGRAM*20.179;
      molar_volume=0.02242;
      volume=19.9052;
      Vm_Miedema=NNN;
      valence_std=0;
      valence_iupac=0;
      valence_PT=0;
      density_PT=9E-4;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.00082783369;
      lattice_constants[1]=442.9;lattice_constants[2]=442.9;lattice_constants[3]=442.9;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Gas";
      radius=0.160;
      radius_PT=38;
      radius_covalent=0.58;
      radius_covalent_PT=58;
      radius_VanDerWaals_PT=154;
      radii_Ghosh08=0.3618;
      radii_Slatter=NNN;
      radii_Pyykko=0.67;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=10.2303;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=10.459;
      electronegativity_Allen=4.787;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=-246.08;
      temperature_melting=-248.59;
      vaporization_heat_PT=1.75;
      specific_heat_PT=1030;
      critical_Pressure=27.24;
      critical_Temperature_PT=44.4;
      thermal_expansion=NNN;
      conductivity_thermal=0.0491;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=17.87;
      hardness_RB=15.45;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-4.1E-9;
      susceptibility_magnetic_volume=-3.69E-9;
      susceptibility_magnetic_molar=-8.27E-11;
      temperature_Curie=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000067;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Ne volume calculated with fcc-pawpbe
      return; //CO20200520
    }
    // [AFLOW]STOP=Neon
    // ********************************************************************************************************************************************************

    // ROW3
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Sodium
    // Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium Sodium
    else if(ZZ==11) { // Sodium
      Z=ZZ;
      symbol="Na";
      name="Sodium";
      period=3;
      group=1;
      series="AlkaliMetal";
      block="s";
      mass=AMU2KILOGRAM*22.9898;
      molar_volume=0.00002375;
      volume=36.9135;
      Vm_Miedema=8.3;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      density_PT=0.968;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=0.0;
      lattice_constants[1]=429.06;lattice_constants[2]=429.06;lattice_constants[3]=429.06;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.186;
      radius_PT=190;
      radius_covalent=1.66;
      radius_covalent_PT=166;
      radius_VanDerWaals_PT=227;
      radii_Ghosh08=2.1650;
      radii_Slatter=1.80;
      radii_Pyykko=1.55;
      conductivity_electrical=2.1E7;
      electronegativity_vec=0.93;
      hardness_Ghosh=2.4441;
      electronegativity_Pearson=2.85;
      electronegativity_Ghosh=2.536;
      electronegativity_Allen=0.869;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(1);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(1);  //RF+SK20200410
      electron_affinity_PT=52.8;
      phi_star_Miedema=2.70;
      nws_Miedema=0.82;
      gamma_s_Miedema=260;
      scale_Pettifor=0.40;
      temperature_boiling=883;
      temperature_melting=97.72;
      vaporization_heat_PT=97.7;
      specific_heat_PT=1230;
      critical_Pressure=345.4;
      critical_Temperature_PT=2573;
      thermal_expansion=0.00007;
      conductivity_thermal=140;
      hardness_Brinell=0.69;
      hardness_Mohs=0.5;
      hardness_Vickers=NNN;
      hardness_Pearson=2.30;
      hardness_Putz=0.66;
      hardness_RB=2.91;
      modulus_shear=3.3;
      modulus_Young=10;
      modulus_bulk=6.3;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=1.6;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=8.8E-9;
      susceptibility_magnetic_volume=8.6E-6;
      susceptibility_magnetic_molar=2E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1100;
      HHIR=500;
      /*xray_scatt=NNN;*/
      // Na
      return; //CO20200520
    }
    // [AFLOW]STOP=Sodium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Magnesium
    // Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium Magnesium
    else if(ZZ==12) { // Magnesium
      Z=ZZ;
      symbol="Mg";
      name="Magnesium";
      period=3;
      group=2;
      series="AlkalineEarthMetal";
      block="s";
      mass=AMU2KILOGRAM*24.305;
      molar_volume=0.000013984;
      volume=22.8178;
      Vm_Miedema=5.8;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      density_PT=1.738;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.00073988271;
      lattice_constants[1]=320.94;lattice_constants[2]=320.94;lattice_constants[3]=521.08;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.160;
      radius_PT=145;
      radius_covalent=1.41;
      radius_covalent_PT=141;
      radius_VanDerWaals_PT=173;
      radii_Ghosh08=1.6711;
      radii_Slatter=1.50;
      radii_Pyykko=1.39;
      conductivity_electrical=2.3E7;
      electronegativity_vec=1.31;
      hardness_Ghosh=3.0146;
      electronegativity_Pearson=3.75;
      electronegativity_Ghosh=3.310;
      electronegativity_Allen=1.293;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=3.45;
      nws_Miedema=1.17;
      gamma_s_Miedema=790;
      scale_Pettifor=1.28;
      temperature_boiling=1090;
      temperature_melting=650;
      vaporization_heat_PT=128;
      specific_heat_PT=1020;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000248;
      conductivity_thermal=160;
      hardness_Brinell=260;
      hardness_Mohs=2.5;
      hardness_Vickers=NNN;
      hardness_Pearson=3.90;
      hardness_Putz=0.93;
      hardness_RB=4.63;
      modulus_shear=17;
      modulus_Young=45;
      modulus_bulk=45;
      Poisson_ratio_PT=0.29;
      BVm_Miedema=5.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=6.9E-9;
      susceptibility_magnetic_volume=0.000012;
      susceptibility_magnetic_molar=1.68E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=5300;
      HHIR=500;
      /*xray_scatt=NNN;*/
      //Mg
      return; //CO20200520
    }
    // [AFLOW]STOP=Magnesium
    // ********************************************************************************************************************************************************

    // p-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Aluminium
    //Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium Aluminium
    else if(ZZ==13) { // Aluminium
      Z=ZZ;
      symbol="Al";
      name="Aluminium";
      period=3;
      group=13;
      series="PoorMetal";
      block="p";
      mass=AMU2KILOGRAM*26.9815;
      molar_volume=9.99E-6;
      volume=16.4000;
      Vm_Miedema=4.6;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      density_PT=2.7;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.0;
      lattice_constants[1]=404.95;lattice_constants[2]=404.95;lattice_constants[3]=404.95;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.143;
      radius_PT=118;
      radius_covalent=1.21;
      radius_covalent_PT=121;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.3608;
      radii_Slatter=1.25;
      radii_Pyykko=1.26;
      conductivity_electrical=3.8E7;
      electronegativity_vec=1.61;
      hardness_Ghosh=3.5849;
      electronegativity_Pearson=3.23;
      electronegativity_Ghosh=4.084;
      electronegativity_Allen=1.613;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(3);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);  //RF+SK20200410
      electron_affinity_PT=42.5;
      phi_star_Miedema=4.20;
      nws_Miedema=1.39;
      gamma_s_Miedema=1200;
      scale_Pettifor=1.66;
      temperature_boiling=2519;
      temperature_melting=660.32;
      vaporization_heat_PT=293;
      specific_heat_PT=904;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000231;
      conductivity_thermal=235;
      hardness_Brinell=245;
      hardness_Mohs=2.75;
      hardness_Vickers=167;
      hardness_Pearson=2.77;
      hardness_Putz=1.42;
      hardness_RB=2.94;
      modulus_shear=26;
      modulus_Young=70;
      modulus_bulk=76;
      Poisson_ratio_PT=0.35;
      BVm_Miedema=7.2;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=7.8E-9;
      susceptibility_magnetic_volume=0.0000211;
      susceptibility_magnetic_molar=2.1E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1600;
      HHIR=1000;
      /*xray_scatt=NNN;*/
      //Al
      return; //CO20200520
    }
    // [AFLOW]STOP=Aluminium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Silicon
    // Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon Silicon
    else if(ZZ==14) { // Silicon
      Z=ZZ;
      symbol="Si";
      name="Silicon";
      period=3;
      group=14;
      series="Metalloid";
      block="p";
      mass=AMU2KILOGRAM*28.0855;
      molar_volume=0.000012054;
      volume=14.3536;
      Vm_Miedema=4.2;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      density_PT=2.33;
      crystal="dia";
      crystal_structure_PT="Tetrahedral_Packing";
      space_group="Fd_3m";
      space_group_number=227;
      variance_parameter_mass=0.00020046752;
      lattice_constants[1]=543.09;lattice_constants[2]=543.09;lattice_constants[3]=543.09;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.117;
      radius_PT=111;
      radius_covalent=1.11;
      radius_covalent_PT=111;
      radius_VanDerWaals_PT=210;
      radii_Ghosh08=1.1477;
      radii_Slatter=1.10;
      radii_Pyykko=1.16;
      conductivity_electrical=1000;
      electronegativity_vec=1.90;
      hardness_Ghosh=4.1551;
      electronegativity_Pearson=4.77;
      electronegativity_Ghosh=4.857;
      electronegativity_Allen=1.916;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(-4);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4);  //RF+SK20200410
      electron_affinity_PT=133.6;
      phi_star_Miedema=4.70;
      nws_Miedema=1.50;
      gamma_s_Miedema=1290;
      scale_Pettifor=1.92;
      temperature_boiling=2900;
      temperature_melting=1414;
      vaporization_heat_PT=359;
      specific_heat_PT=710;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=2.6E-6;
      conductivity_thermal=150;
      hardness_Brinell=NNN;
      hardness_Mohs=6.5;
      hardness_Vickers=9630.1303;
      hardness_Pearson=3.38;
      hardness_Putz=2.10;
      hardness_RB=3.61;
      modulus_shear=NNN;
      modulus_Young=47;
      modulus_bulk=100;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=11.9;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-1.6E-9;
      susceptibility_magnetic_volume=-3.73E-6;
      susceptibility_magnetic_molar=-4.49E-11;
      temperature_Curie=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=4700;
      HHIR=1000;
      xray_scatt=14.43;
      //Si ???
      return; //CO20200520
    }
    // [AFLOW]STOP=Silicon
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Phosphorus
    // Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus Phosphorus
    else if(ZZ==15) { // Phosphorus
      Z=ZZ;
      symbol="P";
      name="Phosphorus";
      period=3;
      group=15;
      series="Nonmetal";
      block="p";
      mass=AMU2KILOGRAM*30.9738;
      molar_volume=0.000016991;
      volume=14.1995;
      Vm_Miedema=NNN;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      density_PT=1.823;
      crystal="cub";
      crystal_structure_PT="Simple_Triclinic";
      space_group="P-1";
      space_group_number=2;
      variance_parameter_mass=0.0;
      lattice_constants[1]=1145;lattice_constants[2]=550.3;lattice_constants[3]=1126.1;
      lattice_angles[1]=1.25384;lattice_angles[2]=1.57725;lattice_angles[3]=1.24896;
      phase="Solid";
      radius=0.109;
      radius_PT=98;
      radius_covalent=1.07;
      radius_covalent_PT=107;
      radius_VanDerWaals_PT=180;
      radii_Ghosh08=0.9922;
      radii_Slatter=1.00;
      radii_Pyykko=1.11;
      conductivity_electrical=1E7;
      electronegativity_vec=2.19;
      hardness_Ghosh=4.7258;
      electronegativity_Pearson=5.62;
      electronegativity_Ghosh=5.631;
      electronegativity_Allen=2.253;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(5); oxidation_states.push_back(3); oxidation_states.push_back(-3); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(5);  //RF+SK20200410
      electron_affinity_PT=71;
      phi_star_Miedema=5.5;
      nws_Miedema=1.65;
      gamma_s_Miedema=NNN;
      scale_Pettifor=2.18;
      temperature_boiling=280.5;
      temperature_melting=44.2;
      vaporization_heat_PT=12.4;
      specific_heat_PT=769.7;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=0.236;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=4.88;
      hardness_Putz=2.92;
      hardness_RB=5.42;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=11;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-1.13E-8;
      susceptibility_magnetic_volume=-0.0000206;
      susceptibility_magnetic_molar=-3.5E-10;
      temperature_Curie=NNN;
      color_PT="COLORLESS";
      refractive_index=1.001212;
      HHIP=2000;
      HHIR=5100;
      xray_scatt=15.3133;
      //P MIEDEMA =PAUL VAN DER PUT book
      return; //CO20200520
    }
    // [AFLOW]STOP=Phosphorus
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Sulphur
    // Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur Sulphur
    else if(ZZ==16) { // Sulphur
      Z=ZZ;
      symbol="S";
      name="Sulphur";
      period=3;
      group=16;
      series="Chalcogen";
      block="p";
      mass=AMU2KILOGRAM*32.06;
      molar_volume=0.000016357;
      volume=15.7301;
      Vm_Miedema=4.376;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      density_PT=1.96;
      crystal="orc";
      crystal_structure_PT="Face-centered_Orthorhombic";
      space_group="Fddd";
      space_group_number=70;
      variance_parameter_mass=0.00016807795;
      lattice_constants[1]=1043.7;lattice_constants[2]=1284.5;lattice_constants[3]=2436.9;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.106;
      radius_PT=87;
      radius_covalent=1.05;
      radius_covalent_PT=105;
      radius_VanDerWaals_PT=180;
      radii_Ghosh08=0.8739;
      radii_Slatter=1.00;
      radii_Pyykko=1.03;
      conductivity_electrical=1E-15;
      electronegativity_vec=2.58;
      hardness_Ghosh=5.2960;
      electronegativity_Pearson=6.22;
      electronegativity_Ghosh=6.420;
      electronegativity_Allen=2.589;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(6); oxidation_states.push_back(4); oxidation_states.push_back(2); oxidation_states.push_back(-2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(6);  //RF+SK20200410
      electron_affinity_PT=200;
      phi_star_Miedema=5.6;
      nws_Miedema=1.46;
      gamma_s_Miedema=NNN;
      scale_Pettifor=2.44;
      temperature_boiling=444.72;
      temperature_melting=115.21;
      vaporization_heat_PT=9.8;
      specific_heat_PT=705;
      critical_Pressure=204.3;
      critical_Temperature_PT=1314;
      thermal_expansion=NNN;
      conductivity_thermal=0.205;
      hardness_Brinell=NNN;
      hardness_Mohs=2;
      hardness_Vickers=NNN;
      hardness_Pearson=4.14;
      hardness_Putz=3.82;
      hardness_RB=4.28;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=7.7;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-6.2E-9;
      susceptibility_magnetic_volume=-0.0000122;
      susceptibility_magnetic_molar=-1.99E-10;
      temperature_Curie=NNN;
      color_PT="YELLOW";
      refractive_index=1.001111;
      HHIP=700;
      HHIR=1000;
      /*xray_scatt=NNN;*/
      //S Table 27 of JX
      return; //CO20200520
    }
    // [AFLOW]STOP=Sulphur
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Chlorine
    // Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine Chlorine
    else if(ZZ==17) { // Chlorine
      Z=ZZ;
      symbol="Cl";
      name="Chlorine";
      period=3;
      group=17;
      series="Halogen";
      block="p";
      mass=AMU2KILOGRAM*35.453;
      molar_volume=0.01103;
      volume=21.2947;
      Vm_Miedema=6.71;
      valence_std=7;
      valence_iupac=7;
      valence_PT=5;
      density_PT=32.14E-4;
      crystal="orc";
      crystal_structure_PT="Base_Orthorhombic";
      space_group="Cmca";
      space_group_number=64;
      variance_parameter_mass=0.00058238731;
      lattice_constants[1]=622.35;lattice_constants[2]=445.61;lattice_constants[3]=817.85;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Gas";
      radius=0.107;
      radius_PT=79;
      radius_covalent=1.02;
      radius_covalent_PT=102;
      radius_VanDerWaals_PT=175;
      radii_Ghosh08=0.7808;
      radii_Slatter=1.00;
      radii_Pyykko=0.99;
      conductivity_electrical=0.01;
      electronegativity_vec=3.16;
      hardness_Ghosh=5.8662;
      electronegativity_Pearson=8.30;
      electronegativity_Ghosh=7.178;
      electronegativity_Allen=2.869;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(7); oxidation_states.push_back(5); oxidation_states.push_back(3); oxidation_states.push_back(1); oxidation_states.push_back(-1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(-1); //RF+SK20200410
      electron_affinity_PT=349;
      phi_star_Miedema=5.32;
      nws_Miedema=0.34;
      gamma_s_Miedema=1013;
      scale_Pettifor=2.70;
      temperature_boiling=-34.04;
      temperature_melting=-101.5;
      vaporization_heat_PT=10.2;
      specific_heat_PT=478.2;
      critical_Pressure=78.87;
      critical_Temperature_PT=416.9;
      thermal_expansion=NNN;
      conductivity_thermal=0.0089;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=4.68;
      hardness_Putz=5.01;
      hardness_RB=4.91;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=1.1;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-7.2E-9;
      susceptibility_magnetic_volume=-2.31E-8;
      susceptibility_magnetic_molar=-5.11E-10;
      temperature_Curie=NNN;
      color_PT="YELLOW";
      refractive_index=1.000773;
      HHIP=1500;
      HHIR=1500;
      /*xray_scatt=NNN;*/
      //Cl interpolation phi_star, nws, Vm, gamma JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Chlorine
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Argon
    //Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon Argon
    else if(ZZ==18) { // Argon
      Z=ZZ;
      symbol="Ar";
      name="Argon";
      period=3;
      group=18;
      series="NobleGas";
      block="p";
      mass=AMU2KILOGRAM*39.948;
      molar_volume=0.022392;
      volume=22.000;
      Vm_Miedema=NNN;
      valence_std=0;
      valence_iupac=2;
      valence_PT=0;
      density_PT=17.84E-4;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.00003509919;
      lattice_constants[1]=525.6;lattice_constants[2]=525.6;lattice_constants[3]=525.6;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Gas";
      radius=0.192;
      radius_PT=71;
      radius_covalent=1.06;
      radius_covalent_PT=106;
      radius_VanDerWaals_PT=188;
      radii_Ghosh08=0.7056;
      radii_Slatter=NNN;
      radii_Pyykko=0.96;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=6.4366;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=7.951;
      electronegativity_Allen=3.242;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=-185.8;
      temperature_melting=-189.3;
      vaporization_heat_PT=6.5;
      specific_heat_PT=520.33;
      critical_Pressure=48.34;
      critical_Temperature_PT=150.87;
      thermal_expansion=NNN;
      conductivity_thermal=0.01772;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=6.16;
      hardness_RB=10.69;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-6E-9;
      susceptibility_magnetic_volume=-1.07E-8;
      susceptibility_magnetic_molar=-2.4E-10;
      temperature_Curie=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000281;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Ar guessed volume, must double check from results JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Argon
    // ********************************************************************************************************************************************************

    // ROW4
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Potassium
    // Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium Potassium
    else if(ZZ==19) { // Potassium
      Z=ZZ;
      symbol="K";
      name="Potassium";
      period=4;
      group=1;
      series="AlkaliMetal";
      block="s";
      mass=AMU2KILOGRAM*39.0983;
      molar_volume=0.00004568;
      volume=73.9091;
      Vm_Miedema=12.8;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      density_PT=0.856;
      crystal="fcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=0.000164;
      lattice_constants[1]=532.8;lattice_constants[2]=532.8;lattice_constants[3]=532.8;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.231;
      radius_PT=243;
      radius_covalent=2.03;
      radius_covalent_PT=203;
      radius_VanDerWaals_PT=275;
      radii_Ghosh08=3.2930;
      radii_Slatter=2.20;
      radii_Pyykko=1.96;
      conductivity_electrical=1.4E7;
      electronegativity_vec=0.82;
      hardness_Ghosh=2.3273;
      electronegativity_Pearson=2.42;
      electronegativity_Ghosh=2.672;
      electronegativity_Allen=0.734;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(1);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(1);  //RF+SK20200410
      electron_affinity_PT=48.4;
      phi_star_Miedema=2.25;
      nws_Miedema=0.65;
      gamma_s_Miedema=150;
      scale_Pettifor=0.35;
      temperature_boiling=759;
      temperature_melting=63.38;
      vaporization_heat_PT=76.9;
      specific_heat_PT=757;
      critical_Pressure=157.9;
      critical_Temperature_PT=2223;
      thermal_expansion=NNN;
      conductivity_thermal=100;
      hardness_Brinell=0.363;
      hardness_Mohs=0.4;
      hardness_Vickers=NNN;
      hardness_Pearson=1.92;
      hardness_Putz=0.18;
      hardness_RB=2.35;
      modulus_shear=1.3;
      modulus_Young=NNN;
      modulus_bulk=3.1;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=1.5;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=6.7E-9;
      susceptibility_magnetic_volume=5.74E-6;
      susceptibility_magnetic_molar=2.62E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1700;
      HHIR=7200;
      /*xray_scatt=NNN;*/
      //K
      return; //CO20200520
    }
    // [AFLOW]STOP=Potassium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Calcium
    // Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium Calcium
    else if(ZZ==20) { // Calcium
      Z=ZZ;
      symbol="Ca";
      name="Calcium";
      period=4;
      group=2;
      series="AlkalineEarthMetal";
      block="s";
      mass=AMU2KILOGRAM*40.08;
      molar_volume=0.000025857;
      volume=42.1927;
      Vm_Miedema=8.8;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      density_PT=1.55;
      crystal="bcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.000297564;
      lattice_constants[1]=558.84;lattice_constants[2]=558.84;lattice_constants[3]=558.84;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.197;
      radius_PT=194;
      radius_covalent=1.76;
      radius_covalent_PT=176;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.5419;
      radii_Slatter=1.80;
      radii_Pyykko=1.71;
      conductivity_electrical=2.9E7;
      electronegativity_vec=1.00;
      hardness_Ghosh=2.7587;
      electronegativity_Pearson=2.2;
      electronegativity_Ghosh=3.140;
      electronegativity_Allen=1.034;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=2.37;
      phi_star_Miedema=2.55;
      nws_Miedema=0.91;
      gamma_s_Miedema=490;
      scale_Pettifor=0.60;
      temperature_boiling=1484;
      temperature_melting=842;
      vaporization_heat_PT=155;
      specific_heat_PT=631;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000223;
      conductivity_thermal=200;
      hardness_Brinell=167;
      hardness_Mohs=1.75;
      hardness_Vickers=NNN;
      hardness_Pearson=4.00;
      hardness_Putz=0.25;
      hardness_RB=3.07;
      modulus_shear=7.4;
      modulus_Young=20;
      modulus_bulk=17;
      Poisson_ratio_PT=0.31;
      BVm_Miedema=4.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.38E-8;
      susceptibility_magnetic_volume=0.00002139;
      susceptibility_magnetic_molar=5.531E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3900;
      HHIR=1500;
      /*xray_scatt=NNN;*/
      //Ca
      return; //CO20200520
    }
    // [AFLOW]STOP=Calcium
    // ********************************************************************************************************************************************************

    // d-electron systems: transition metals
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Scandium
    // Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium Scandium
    else if(ZZ==21) { // Scandium
      Z=ZZ;
      symbol="Sc";
      name="Scandium";
      period=4;
      group=3;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*44.9559;
      molar_volume=0.000015061;
      volume=24.6739;
      Vm_Miedema=6.1;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      density_PT=2.985;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=330.9;lattice_constants[2]=330.9;lattice_constants[3]=527.33;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.160;
      radius_PT=184;
      radius_covalent=1.70;
      radius_covalent_PT=170;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.4149;
      radii_Slatter=1.60;
      radii_Pyykko=1.48;
      conductivity_electrical=1.8E6;
      electronegativity_vec=1.36;
      hardness_Ghosh=2.8582;
      electronegativity_Pearson=3.34;
      electronegativity_Ghosh=3.248;
      electronegativity_Allen=1.190;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(3);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);  //RF+SK20200410
      electron_affinity_PT=18.1;
      phi_star_Miedema=3.25;
      nws_Miedema=1.27;
      gamma_s_Miedema=1200;
      scale_Pettifor=0.74;
      temperature_boiling=2830;
      temperature_melting=1541;
      vaporization_heat_PT=318;
      specific_heat_PT=567;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000102;
      conductivity_thermal=16;
      hardness_Brinell=750;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=3.20;
      hardness_Putz=0.31;
      hardness_RB=2.52;
      modulus_shear=29;
      modulus_Young=74;
      modulus_bulk=57;
      Poisson_ratio_PT=0.28;
      BVm_Miedema=6.6;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=8.8E-8;
      susceptibility_magnetic_volume=0.0002627;
      susceptibility_magnetic_molar=3.956E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=5500;
      HHIR=4500;
      xray_scatt=21.34;
      //Sc
      return; //CO20200520
    }
    // [AFLOW]STOP=Scandium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Titanium
    // Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium Titanium
    else if(ZZ==22) { // Titanium
      Z=ZZ;
      symbol="Ti";
      name="Titanium";
      period=4;
      group=4;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*47.9;
      molar_volume=0.000010621;
      volume=17.1035;
      Vm_Miedema=4.8;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      density_PT=4.507;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.000286456;
      lattice_constants[1]=295.08;lattice_constants[2]=295.08;lattice_constants[3]=468.55;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.147;
      radius_PT=176;
      radius_covalent=1.60;
      radius_covalent_PT=160;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.2998;
      radii_Slatter=1.40;
      radii_Pyykko=1.36;
      conductivity_electrical=2.5E6;
      electronegativity_vec=1.54;
      hardness_Ghosh=2.9578;
      electronegativity_Pearson=3.45;
      electronegativity_Ghosh=3.357;
      electronegativity_Allen=1.38; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(3); oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4);  //RF+SK20200410
      electron_affinity_PT=7.6;
      phi_star_Miedema=3.65;
      nws_Miedema=1.47;
      gamma_s_Miedema=2050;
      scale_Pettifor=0.79;
      temperature_boiling=3287;
      temperature_melting=1668;
      vaporization_heat_PT=425;
      specific_heat_PT=520;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=8.6E-6;
      conductivity_thermal=22;
      hardness_Brinell=715;
      hardness_Mohs=6;
      hardness_Vickers=970;
      hardness_Pearson=3.37;
      hardness_Putz=0.38;
      hardness_RB=2.03;
      modulus_shear=44;
      modulus_Young=116;
      modulus_bulk=110;
      Poisson_ratio_PT=0.32;
      BVm_Miedema=11.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=4.01E-8;
      susceptibility_magnetic_volume=0.0001807;
      susceptibility_magnetic_molar=1.919E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1100;
      HHIR=1600;
      xray_scatt=22.24;
      //Ti
      return; //CO20200520
    }
    // [AFLOW]STOP=Titanium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Vanadium
    // Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium Vanadium
    else if(ZZ==23) { // Vanadium
      Z=ZZ;
      symbol="V";
      name="Vanadium";
      period=4;
      group=5;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*50.9415;
      molar_volume=8.3374E-6;
      volume=13.2086;
      Vm_Miedema=4.1;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      density_PT=6.11;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=9.54831E-07;
      lattice_constants[1]=303;lattice_constants[2]=303;lattice_constants[3]=303;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.132;
      radius_PT=171;
      radius_covalent=1.53;
      radius_covalent_PT=153;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.1953;
      radii_Slatter=1.35;
      radii_Pyykko=1.34;
      conductivity_electrical=5E6;
      electronegativity_vec=1.63;
      hardness_Ghosh=3.0573;
      electronegativity_Pearson=3.6;
      electronegativity_Ghosh=3.465;
      electronegativity_Allen=1.53; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(5); oxidation_states.push_back(4); oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(0);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(5);  //RF+SK20200410
      electron_affinity_PT=50.6;
      phi_star_Miedema=4.25;
      nws_Miedema=1.64;
      gamma_s_Miedema=2600;
      scale_Pettifor=0.84;
      temperature_boiling=3407;
      temperature_melting=1910;
      vaporization_heat_PT=453;
      specific_heat_PT=489;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=8.4E-6;
      conductivity_thermal=31;
      hardness_Brinell=628;
      hardness_Mohs=7;
      hardness_Vickers=628;
      hardness_Pearson=3.10;
      hardness_Putz=0.45;
      hardness_RB=NNN;
      modulus_shear=47;
      modulus_Young=128;
      modulus_bulk=160;
      Poisson_ratio_PT=0.37;
      BVm_Miedema=14.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=6.28E-8;
      susceptibility_magnetic_volume=0.0003837;
      susceptibility_magnetic_molar=3.199E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3300;
      HHIR=3400;
      /*xray_scatt=NNN;*/
      //V
      return; //CO20200520
    }
    // [AFLOW]STOP=Vanadium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Chromium
    // Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium Chromium
    else if(ZZ==24) { // Chromium
      Z=ZZ;
      symbol="Cr";
      name="Chromium";
      period=4;
      group=6;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*51.996;
      molar_volume=7.2317E-6;
      volume=11.4136;
      Vm_Miedema=3.7;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      density_PT=7.19;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=0.00013287;
      lattice_constants[1]=291;lattice_constants[2]=291;lattice_constants[3]=291;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.125;
      radius_PT=166;
      radius_covalent=1.39;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.1000;
      radii_Slatter=1.40;
      radii_Pyykko=1.22;
      conductivity_electrical=7.9E6;
      electronegativity_vec=1.66;
      hardness_Ghosh=3.1567;
      electronegativity_Pearson=3.72;
      electronegativity_Ghosh=3.573;
      electronegativity_Allen=1.650;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(6); oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(0); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3); oxidation_states_preferred.push_back(6); //RF+SK20200410 //Cr+3 most preferred oxidation number
      electron_affinity_PT=64.3;
      phi_star_Miedema=4.65;
      nws_Miedema=1.74;
      gamma_s_Miedema=2400;
      scale_Pettifor=0.89;
      temperature_boiling=2671;
      temperature_melting=1907;
      vaporization_heat_PT=339;
      specific_heat_PT=448;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=4.9E-6;
      conductivity_thermal=94;
      hardness_Brinell=1120;
      hardness_Mohs=8.5;
      hardness_Vickers=1060;
      hardness_Pearson=3.06;
      hardness_Putz=0.54;
      hardness_RB=4.06;
      modulus_shear=115;
      modulus_Young=279;
      modulus_bulk=160;
      Poisson_ratio_PT=0.21;
      BVm_Miedema=14.0;
      magnetic_type_PT="Antiferromagnetic";
      susceptibility_magnetic_mass=4.45E-8;
      susceptibility_magnetic_volume=0.0003177;
      susceptibility_magnetic_molar=2.314E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3100;
      HHIR=4100;
      xray_scatt=23.84;
      //Cr
      return; //CO20200520
    }
    // [AFLOW]STOP=Chromium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Manganese
    // Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese Manganese
    else if(ZZ==25) { // Manganese
      Z=ZZ;
      symbol="Mn";
      name="Manganese";
      period=4;
      group=7;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*54.93805;
      molar_volume=7.3545E-6;
      volume=10.6487;
      Vm_Miedema=3.8;
      valence_std=7;
      valence_iupac=7;
      valence_PT=4;
      density_PT=7.47;
      crystal="cub";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="I_43m";
      space_group_number=217;
      variance_parameter_mass=1.67276E-32;
      lattice_constants[1]=891.25;lattice_constants[2]=891.25;lattice_constants[3]=891.25;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.112;
      radius_PT=161;
      radius_covalent=1.61;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.0124;
      radii_Slatter=1.40;
      radii_Pyykko=1.19;
      conductivity_electrical=620000;
      electronegativity_vec=1.55;
      hardness_Ghosh=3.2564;
      electronegativity_Pearson=3.72;
      electronegativity_Ghosh=3.681;
      electronegativity_Allen=1.75; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(7); oxidation_states.push_back(6); oxidation_states.push_back(4); oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(0); oxidation_states.push_back(-1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=4.45;
      nws_Miedema=1.61;
      gamma_s_Miedema=1600;
      scale_Pettifor=0.94;
      temperature_boiling=2061;
      temperature_melting=1246;
      vaporization_heat_PT=220;
      specific_heat_PT=479;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000217;
      conductivity_thermal=7.7;
      hardness_Brinell=196;
      hardness_Mohs=6;
      hardness_Vickers=NNN;
      hardness_Pearson=3.72;
      hardness_Putz=0.64;
      hardness_RB=2.88;
      modulus_shear=NNN;
      modulus_Young=198;
      modulus_bulk=120;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=4.4;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.21E-7;
      susceptibility_magnetic_volume=0.00090387;
      susceptibility_magnetic_molar=6.6475E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1600;
      HHIR=1800;
      xray_scatt=24.46;
      //xray_scatt=24.3589; Mn JX CHANGED VALENCE //DX+CO20170904 radius_covalent[i] uses high spin configuration (most frequent)
      return; //CO20200520
    }
    // [AFLOW]STOP=Manganese
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Iron
    // Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron Iron
    else if(ZZ==26) { // Iron
      Z=ZZ;
      symbol="Fe";
      name="Iron";
      period=4;
      group=8;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*55.847;
      molar_volume=7.0923E-6;
      volume=10.2315;
      Vm_Miedema=3.7;
      valence_std=8;
      valence_iupac=6;
      valence_PT=3;
      density_PT=7.874;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=9.17912E-05;
      lattice_constants[1]=286.65;lattice_constants[2]=286.65;lattice_constants[3]=286.65;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.124;
      radius_PT=156;
      radius_covalent=1.52;
      radius_covalent_PT=132;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.9319;
      radii_Slatter=1.40;
      radii_Pyykko=1.16;
      conductivity_electrical=1E7;
      electronegativity_vec=1.83;
      hardness_Ghosh=3.3559;
      electronegativity_Pearson=4.06;
      electronegativity_Ghosh=3.789;
      electronegativity_Allen=1.80; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(6); oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(0); oxidation_states.push_back(-2); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3); oxidation_states_preferred.push_back(2); //RF+SK20200410 //Fe+3 most preferred oxidation number
      electron_affinity_PT=15.7;
      phi_star_Miedema=4.93;
      nws_Miedema=1.77;
      gamma_s_Miedema=2550;
      scale_Pettifor=0.99;
      temperature_boiling=2861;
      temperature_melting=1538;
      vaporization_heat_PT=347;
      specific_heat_PT=449;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000118;
      conductivity_thermal=79;
      hardness_Brinell=490;
      hardness_Mohs=4;
      hardness_Vickers=608;
      hardness_Pearson=3.81;
      hardness_Putz=0.75;
      hardness_RB=2.53;
      modulus_shear=82;
      modulus_Young=211;
      modulus_bulk=170;
      Poisson_ratio_PT=0.29;
      BVm_Miedema=12.0;
      magnetic_type_PT="Ferromagnetic";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=1043;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=2400;
      HHIR=1400;
      xray_scatt=24.85;
      //xray_scatt=24.6830; Fe JX CHANGED VALENCE //DX+CO20170904 radius_covalent[i] uses high spin configuration (most frequent)
      return; //CO20200520
    }
    // [AFLOW]STOP=Iron
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Cobalt
    // Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt Cobalt
    else if(ZZ==27) { // Cobalt
      Z=ZZ;
      symbol="Co";
      name="Cobalt";
      period=4;
      group=9;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*58.9332;
      molar_volume=6.62E-6;
      volume=10.3205;
      Vm_Miedema=3.5;
      valence_std=9;
      valence_iupac=5;
      valence_PT=4;
      density_PT=8.9;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=250.71;lattice_constants[2]=250.71;lattice_constants[3]=406.95;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.125;
      radius_PT=152;
      radius_covalent=1.26;
      radius_covalent_PT=126;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.8575;
      radii_Slatter=1.35;
      radii_Pyykko=1.11;
      conductivity_electrical=1.7E7;
      electronegativity_vec=1.88;
      hardness_Ghosh=3.4556;
      electronegativity_Pearson=4.3;
      electronegativity_Ghosh=3.897;
      electronegativity_Allen=1.84; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(0); oxidation_states.push_back(-1);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=63.7;
      phi_star_Miedema=5.10;
      nws_Miedema=1.75;
      gamma_s_Miedema=2550;
      scale_Pettifor=1.04;
      temperature_boiling=2927;
      temperature_melting=1495;
      vaporization_heat_PT=375;
      specific_heat_PT=421;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000013;
      conductivity_thermal=100;
      hardness_Brinell=700;
      hardness_Mohs=5;
      hardness_Vickers=1043;
      hardness_Pearson=3.60;
      hardness_Putz=0.88;
      hardness_RB=3.53;
      modulus_shear=76;
      modulus_Young=209;
      modulus_bulk=180;
      Poisson_ratio_PT=0.31;
      BVm_Miedema=13.0;
      magnetic_type_PT="Ferromagnetic";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=1394;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=3100;
      HHIR=2700;
      xray_scatt=24.59;
      //Co JX CHANGED VALENCE //DX+CO20170904 radius_covalent[i] uses low spin configuration (most frequent)
      return; //CO20200520
    }
    // [AFLOW]STOP=Cobalt
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Nickel
    // Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel Nickel
    else if(ZZ==28) { // Nickel
      Z=ZZ;
      symbol="Ni";
      name="Nickel";
      period=4;
      group=10;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*58.69;
      molar_volume=6.5888E-6;
      volume=10.8664;
      Vm_Miedema=3.5;
      valence_std=10;
      valence_iupac=4;
      valence_PT=2;
      density_PT=8.908;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.000430773;
      lattice_constants[1]=352.4;lattice_constants[2]=352.4;lattice_constants[3]=352.4;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.125;
      radius_PT=149;
      radius_covalent=1.24;
      radius_covalent_PT=124;
      radius_VanDerWaals_PT=163;
      radii_Ghosh08=1.7888;
      radii_Slatter=1.35;
      radii_Pyykko=1.10;
      conductivity_electrical=1.4E7;
      electronegativity_vec=1.91;
      hardness_Ghosh=3.5550;
      electronegativity_Pearson=4.40;
      electronegativity_Ghosh=4.005;
      electronegativity_Allen=1.88; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(0);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=112;
      phi_star_Miedema=5.20;
      nws_Miedema=1.75;
      gamma_s_Miedema=2450;
      scale_Pettifor=1.09;
      temperature_boiling=2913;
      temperature_melting=1455;
      vaporization_heat_PT=378;
      specific_heat_PT=445;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000134;
      conductivity_thermal=91;
      hardness_Brinell=700;
      hardness_Mohs=4;
      hardness_Vickers=638;
      hardness_Pearson=3.25;
      hardness_Putz=1.02;
      hardness_RB=4.08;
      modulus_shear=76;
      modulus_Young=200;
      modulus_bulk=180;
      Poisson_ratio_PT=0.31;
      BVm_Miedema=12.0;
      magnetic_type_PT="Ferromagnetic";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=631;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=1000;
      HHIR=1500;
      xray_scatt=25.02;
      //Ni
      return; //CO20200520
    }
    // [AFLOW]STOP=Nickel
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Copper
    // Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper Copper
    else if(ZZ==29) { // Copper
      Z=ZZ;
      symbol="Cu";
      name="Copper";
      period=4;
      group=11;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*63.546;
      molar_volume=7.0922E-6;
      volume=12.0159;
      Vm_Miedema=3.7;
      valence_std=11;
      valence_iupac=4;
      valence_PT=2;
      density_PT=8.96;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.00021086;
      lattice_constants[1]=361.49;lattice_constants[2]=361.49;lattice_constants[3]=361.49;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.128;
      radius_PT=145;
      radius_covalent=1.32;
      radius_covalent_PT=132;
      radius_VanDerWaals_PT=140;
      radii_Ghosh08=1.725;
      radii_Slatter=1.35;
      radii_Pyykko=1.12;
      conductivity_electrical=5.9E7;
      electronegativity_vec=1.90;
      hardness_Ghosh=3.6544;
      electronegativity_Pearson=4.48;
      electronegativity_Ghosh=4.113;
      electronegativity_Allen=1.85; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2); oxidation_states.push_back(1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2); oxidation_states_preferred.push_back(1); //RF+SK20200410 //Cu+2 most preferred oxidation number
      electron_affinity_PT=118.4;
      phi_star_Miedema=4.55;
      nws_Miedema=1.47;
      gamma_s_Miedema=1850;
      scale_Pettifor=1.20;
      temperature_boiling=2562;
      temperature_melting=1084.62;
      vaporization_heat_PT=300;
      specific_heat_PT=384.4;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000165;
      conductivity_thermal=400;
      hardness_Brinell=874;
      hardness_Mohs=3;
      hardness_Vickers=369;
      hardness_Pearson=3.25;
      hardness_Putz=1.21;
      hardness_RB=NNN;
      modulus_shear=48;
      modulus_Young=130;
      modulus_bulk=140;
      Poisson_ratio_PT=0.34;
      BVm_Miedema=9.3;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-1.08E-9;
      susceptibility_magnetic_volume=-9.63E-6;
      susceptibility_magnetic_molar=-6.86E-11;
      temperature_Curie=NNN;
      color_PT="COPPER";
      refractive_index=NNN;
      HHIP=1600;
      HHIR=1500;
      xray_scatt=27.03;
      //Cu JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Copper
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Zinc
    // Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc Zinc
    else if(ZZ==30) { // Zinc
      Z=ZZ;
      symbol="Zn";
      name="Zinc";
      period=4;
      group=12;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*65.38;
      molar_volume=9.157E-6;
      volume=15.0827;
      Vm_Miedema=4.4;
      valence_std=12;
      valence_iupac=2;
      valence_PT=2;
      density_PT=7.14;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.000595597;
      lattice_constants[1]=266.49;lattice_constants[2]=266.49;lattice_constants[3]=494.68;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.133;
      radius_PT=142;
      radius_covalent=1.22;
      radius_covalent_PT=122;
      radius_VanDerWaals_PT=139;
      radii_Ghosh08=1.6654;
      radii_Slatter=1.35;
      radii_Pyykko=1.18;
      conductivity_electrical=1.7E7;
      electronegativity_vec=1.65;
      hardness_Ghosh=3.7542;
      electronegativity_Pearson=4.45;
      electronegativity_Ghosh=4.222;
      electronegativity_Allen=1.59; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=4.10;
      nws_Miedema=1.32;
      gamma_s_Miedema=1020;
      scale_Pettifor=1.44;
      temperature_boiling=907;
      temperature_melting=419.53;
      vaporization_heat_PT=119;
      specific_heat_PT=388;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000302;
      conductivity_thermal=120;
      hardness_Brinell=412;
      hardness_Mohs=2.5;
      hardness_Vickers=NNN;
      hardness_Pearson=4.94;
      hardness_Putz=1.39;
      hardness_RB=6.01;
      modulus_shear=43;
      modulus_Young=108;
      modulus_bulk=70;
      Poisson_ratio_PT=0.25;
      BVm_Miedema=5.5;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-2.21E-9;
      susceptibility_magnetic_volume=-0.0000158;
      susceptibility_magnetic_molar=-1.45E-10;
      temperature_Curie=NNN;
      color_PT="SLATEGRAY";
      refractive_index=1.00205;
      HHIP=1600;
      HHIR=1900;
      xray_scatt=28.44;
      //Zn
      return; //CO20200520
    }
    // [AFLOW]STOP=Zinc
    // ********************************************************************************************************************************************************

    // p-electron systems 
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Gallium
    // Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium Gallium
    else if(ZZ==31) { // Gallium
      Z=ZZ;
      symbol="Ga";
      name="Gallium";
      period=4;
      group=13;
      series="PoorMetal";
      block="p";
      mass=AMU2KILOGRAM*69.737;
      molar_volume=0.000011809;
      volume=18.9039;
      Vm_Miedema=5.2;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      density_PT=5.904;
      crystal="orc";
      crystal_structure_PT="Base_Orthorhombic";
      space_group="Cmca";
      space_group_number=64;
      variance_parameter_mass=0.000197588;
      lattice_constants[1]=451.97;lattice_constants[2]=766.33;lattice_constants[3]=452.6;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.135;
      radius_PT=136;
      radius_covalent=1.22;
      radius_covalent_PT=122;
      radius_VanDerWaals_PT=187;
      radii_Ghosh08=1.4489;
      radii_Slatter=1.30;
      radii_Pyykko=1.24;
      conductivity_electrical=7.1E6;
      electronegativity_vec=1.81;
      hardness_Ghosh=4.1855;
      electronegativity_Pearson=3.2;
      electronegativity_Ghosh=4.690;
      electronegativity_Allen=1.756;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(3);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);  //RF+SK20200410
      electron_affinity_PT=28.9;
      phi_star_Miedema=4.10;
      nws_Miedema=1.31;
      gamma_s_Miedema=830;
      scale_Pettifor=1.68;
      temperature_boiling=2204;
      temperature_melting=29.76;
      vaporization_heat_PT=256;
      specific_heat_PT=371;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.00012;
      conductivity_thermal=29;
      hardness_Brinell=60;
      hardness_Mohs=1.5;
      hardness_Vickers=NNN;
      hardness_Pearson=2.90;
      hardness_Putz=1.59;
      hardness_RB=3.03;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=6.7;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-3E-9;
      susceptibility_magnetic_volume=-0.0000177;
      susceptibility_magnetic_molar=-2.09E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=5500;
      HHIR=1900;
      /*xray_scatt=NNN;*/
      //Ga
      return; //CO20200520
    }
    // [AFLOW]STOP=Gallium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Germanium
    // Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium Germanium
    else if(ZZ==32) { // Germanium
      Z=ZZ;
      symbol="Ge";
      name="Germanium";
      period=4;
      group=14;
      series="Metalloid";
      block="p";
      mass=AMU2KILOGRAM*72.59;
      molar_volume=0.000013645;
      volume=19.2948;
      Vm_Miedema=4.6;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      density_PT=5.323;
      crystal="dia";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.00058782;
      lattice_constants[1]=565.75;lattice_constants[2]=565.75;lattice_constants[3]=565.75;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.122;
      radius_PT=125;
      radius_covalent=1.20;
      radius_covalent_PT=120;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.2823;
      radii_Slatter=1.25;
      radii_Pyykko=1.21;
      conductivity_electrical=2000;
      electronegativity_vec=2.01;
      hardness_Ghosh=4.6166;
      electronegativity_Pearson=4.6;
      electronegativity_Ghosh=5.159;
      electronegativity_Allen=1.994;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(4);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4);  //RF+SK20200410
      electron_affinity_PT=119;
      phi_star_Miedema=4.55;
      nws_Miedema=1.37;
      gamma_s_Miedema=1030;
      scale_Pettifor=1.92;
      temperature_boiling=2820;
      temperature_melting=938.3;
      vaporization_heat_PT=334;
      specific_heat_PT=321.4;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6E-6;
      conductivity_thermal=60;
      hardness_Brinell=7273.402498871;
      hardness_Mohs=6;
      hardness_Vickers=8012.03305;
      hardness_Pearson=3.40;
      hardness_Putz=1.94;
      hardness_RB=3.52;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=10.5;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-1.5E-9;
      susceptibility_magnetic_volume=-7.98E-6;
      susceptibility_magnetic_molar=-1.09E-10;
      temperature_Curie=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=5300;
      HHIR=1900;
      /*xray_scatt=NNN;*/
      //Ge
      return; //CO20200520
    }
    // [AFLOW]STOP=Germanium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Arsenic
    // Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic Arsenic
    else if(ZZ==33) { // Arsenic
      Z=ZZ;
      symbol="As";
      name="Arsenic";
      period=4;
      group=15;
      series="Metalloid";
      block="p";
      mass=AMU2KILOGRAM*74.9216;
      molar_volume=0.000013082;
      volume=19.0677;
      Vm_Miedema=5.2;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      density_PT=5.727;
      crystal="rhl";
      crystal_structure_PT="Simple_Trigonal";
      space_group="R_3m";
      space_group_number=166;
      variance_parameter_mass=0.0;
      lattice_constants[1]=375.98;lattice_constants[2]=375.98;lattice_constants[3]=1054.75;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.125;
      radius_PT=114;
      radius_covalent=1.19;
      radius_covalent_PT=119;
      radius_VanDerWaals_PT=185;
      radii_Ghosh08=1.1450;
      radii_Slatter=1.15;
      radii_Pyykko=1.21;
      conductivity_electrical=3.3E6;
      electronegativity_vec=2.18;
      hardness_Ghosh=5.0662;
      electronegativity_Pearson=5.3;
      electronegativity_Ghosh=5.628;
      electronegativity_Allen=2.211;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(5); oxidation_states.push_back(3); oxidation_states.push_back(-3); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);  //RF+SK20200410
      electron_affinity_PT=78;
      phi_star_Miedema=4.80;
      nws_Miedema=1.44;
      gamma_s_Miedema=1000;
      scale_Pettifor=2.16;
      temperature_boiling=614;
      temperature_melting=817;
      vaporization_heat_PT=32.4;
      specific_heat_PT=328;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=50;
      hardness_Brinell=1440;
      hardness_Mohs=3.5;
      hardness_Vickers=1510;
      hardness_Pearson=4.50;
      hardness_Putz=2.35;
      hardness_RB=5.04;
      modulus_shear=NNN;
      modulus_Young=8;
      modulus_bulk=22;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=5.1;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-3.9E-9;
      susceptibility_magnetic_volume=-0.0000223;
      susceptibility_magnetic_molar=-2.92E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=1.001552;
      HHIP=3300;
      HHIR=4000;
      /*xray_scatt=NNN;*/
      //As
      return; //CO20200520
    }
    // [AFLOW]STOP=Arsenic
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Selenium
    // Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium Selenium
    else if(ZZ==34) { // Selenium
      Z=ZZ;
      symbol="Se";
      name="Selenium";
      period=4;
      group=16;
      series="Chalcogen";
      block="p";
      mass=AMU2KILOGRAM*78.96;
      molar_volume=0.000016387;
      volume=20.3733;
      Vm_Miedema=5.172;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      density_PT=4.819;
      crystal="hex";
      crystal_structure_PT="Simple_Monoclinic";
      space_group="P12_1/c1";
      space_group_number=14;
      variance_parameter_mass=0.00046279;
      lattice_constants[1]=905.4;lattice_constants[2]=908.3;lattice_constants[3]=1160.1;
      lattice_angles[1]=PI/2;lattice_angles[2]=1.58493;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.116;
      radius_PT=103;
      radius_covalent=1.20;
      radius_covalent_PT=120;
      radius_VanDerWaals_PT=190;
      radii_Ghosh08=1.0424;
      radii_Slatter=1.15;
      radii_Pyykko=1.16;
      conductivity_electrical=NNN;
      electronegativity_vec=2.55;
      hardness_Ghosh=5.4795;
      electronegativity_Pearson=5.89;
      electronegativity_Ghosh=6.096;
      electronegativity_Allen=2.424;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(6); oxidation_states.push_back(4); oxidation_states.push_back(-2); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4);  //RF+SK20200410
      electron_affinity_PT=195;
      phi_star_Miedema=5.17;
      nws_Miedema=1.40;
      gamma_s_Miedema=NNN;
      scale_Pettifor=2.40;
      temperature_boiling=685;
      temperature_melting=221;
      vaporization_heat_PT=26;
      specific_heat_PT=321.2;
      critical_Pressure=268.4;
      critical_Temperature_PT=1766;
      thermal_expansion=NNN;
      conductivity_thermal=0.52;
      hardness_Brinell=736;
      hardness_Mohs=2;
      hardness_Vickers=NNN;
      hardness_Pearson=3.87;
      hardness_Putz=2.87;
      hardness_RB=3.95;
      modulus_shear=3.7;
      modulus_Young=10;
      modulus_bulk=8.3;
      Poisson_ratio_PT=0.33;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-4E-9;
      susceptibility_magnetic_volume=-0.0000193;
      susceptibility_magnetic_molar=-3.16E-10;
      temperature_Curie=NNN;
      color_PT="GRAY";
      refractive_index=1.000895;
      HHIP=2200;
      HHIR=1900;
      /*xray_scatt=NNN;*/
      //Se Table 27 of JX
      return; //CO20200520
    }
    // [AFLOW]STOP=Selenium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Bromine
    // Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine Bromine
    else if(ZZ==35) { // Bromine
      Z=ZZ;
      symbol="Br";
      name="Bromine";
      period=4;
      group=17;
      series="Halogen";
      block="p";
      mass=AMU2KILOGRAM*79.904;
      molar_volume=0.00002561;
      volume=26.3292;
      Vm_Miedema=7.31;
      valence_std=7;
      valence_iupac=7;
      valence_PT=5;
      density_PT=3.12;
      crystal="orc";
      crystal_structure_PT="Base_Orthorhombic";
      space_group="Cmca";
      space_group_number=64;
      variance_parameter_mass=0.000156277;
      lattice_constants[1]=672.65;lattice_constants[2]=464.51;lattice_constants[3]=870.23;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Liquid";
      radius=0.119;
      radius_PT=94;
      radius_covalent=1.20;
      radius_covalent_PT=120;
      radius_VanDerWaals_PT=185;
      radii_Ghosh08=0.9532;
      radii_Slatter=1.15;
      radii_Pyykko=1.14;
      conductivity_electrical=1E-10;
      electronegativity_vec=2.96;
      hardness_Ghosh=5.9111;
      electronegativity_Pearson=7.59;
      electronegativity_Ghosh=6.565;
      electronegativity_Allen=2.685;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(7); oxidation_states.push_back(5); oxidation_states.push_back(3); oxidation_states.push_back(1); oxidation_states.push_back(-1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(-1); //RF+SK20200410
      electron_affinity_PT=324.6;
      phi_star_Miedema=5.20;
      nws_Miedema=1.35;
      gamma_s_Miedema=943;
      scale_Pettifor=2.64;
      temperature_boiling=59;
      temperature_melting=-7.3;
      vaporization_heat_PT=14.8;
      specific_heat_PT=947.3;
      critical_Pressure=102;
      critical_Temperature_PT=588;
      thermal_expansion=NNN;
      conductivity_thermal=0.12;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=4.22;
      hardness_Putz=3.39;
      hardness_RB=4.4;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=1.9;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=3.4;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-4.9E-9;
      susceptibility_magnetic_volume=-0.0000153;
      susceptibility_magnetic_molar=-7.83E-10;
      temperature_Curie=NNN;
      color_PT="RED";
      refractive_index=1.001132;
      HHIP=3300;
      HHIR=6900;
      /* xray_scatt=NNN;*/
      //Br interpolation phi_star, nws, Vm, gamma, BVm JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Bromine
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Krypton
    // Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton Krypton
    else if(ZZ==36) { // Krypton
      Z=ZZ;
      symbol="Kr";
      name="Krypton";
      period=4;
      group=18;
      series="NobleGas";
      block="p";
      mass=AMU2KILOGRAM*83.8;
      molar_volume=0.02235;
      volume=-1.0000;
      Vm_Miedema=NNN;
      valence_std=0;
      valence_iupac=2;
      valence_PT=2;
      density_PT=37.5E-4;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.000248482;
      lattice_constants[1]=570.6;lattice_constants[2]=570.6;lattice_constants[3]=570.6;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Gas";
      radius=0.197;
      radius_PT=87;
      radius_covalent=1.16;
      radius_covalent_PT=116;
      radius_VanDerWaals_PT=202;
      radii_Ghosh08=0.8782;
      radii_Slatter=NNN;
      radii_Pyykko=1.17;
      conductivity_electrical=NNN;
      electronegativity_vec=3;
      hardness_Ghosh=6.3418;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=7.033;
      electronegativity_Allen=2.966;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=-153.22;
      temperature_melting=-157.36;
      vaporization_heat_PT=9.02;
      specific_heat_PT=248.05;
      critical_Pressure=54.28;
      critical_Temperature_PT=209.41;
      thermal_expansion=NNN;
      conductivity_thermal=0.00943;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=3.98;
      hardness_RB=9.45;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-4.4E-9;
      susceptibility_magnetic_volume=-1.65E-8;
      susceptibility_magnetic_molar=-3.69E-10;
      temperature_Curie=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000427;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Kr
      return; //CO20200520
    }
    // [AFLOW]STOP=Krypton
    // ********************************************************************************************************************************************************

    // ROW5
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Rubidium
    // Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium Rubidium
    else if(ZZ==37) { // Rubidium
      Z=ZZ;
      symbol="Rb";
      name="Rubidium";
      period=5;
      group=1;
      series="AlkaliMetal";
      block="s";
      mass=AMU2KILOGRAM*85.4678;
      molar_volume=0.000055788;
      volume=91.2738;
      Vm_Miedema=14.6;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      density_PT=1.532;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=0.000109697;
      lattice_constants[1]=558.5;lattice_constants[2]=558.5;lattice_constants[3]=558.5;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.251;
      radius_PT=265;
      radius_covalent=2.20;
      radius_covalent_PT=220;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=3.8487;
      radii_Slatter=2.35;
      radii_Pyykko=2.10;
      conductivity_electrical=8.3E6;
      electronegativity_vec=0.82;
      hardness_Ghosh=2.1204;
      electronegativity_Pearson=2.34;
      electronegativity_Ghosh=2.849;
      electronegativity_Allen=0.706;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(1);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(1);  //RF+SK20200410
      electron_affinity_PT=46.9;
      phi_star_Miedema=2.10;
      nws_Miedema=0.60;
      gamma_s_Miedema=120;
      scale_Pettifor=0.30;
      temperature_boiling=688;
      temperature_melting=39.31;
      vaporization_heat_PT=71;
      specific_heat_PT=364;
      critical_Pressure=157.9;
      critical_Temperature_PT=2093;
      thermal_expansion=NNN;
      conductivity_thermal=58;
      hardness_Brinell=0.216;
      hardness_Mohs=0.3;
      hardness_Vickers=NNN;
      hardness_Pearson=1.85;
      hardness_Putz=0.08;
      hardness_RB=2.21;
      modulus_shear=NNN;
      modulus_Young=2.4;
      modulus_bulk=2.5;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=1.8;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=2.6E-9;
      susceptibility_magnetic_volume=3.98E-6;
      susceptibility_magnetic_molar=2.22E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=6000;
      HHIR=6000;
      /*xray_scatt=NNN;*/
      //Rb
      return; //CO20200520
    }
    // [AFLOW]STOP=Rubidium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Strontium
    // Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium Strontium
    else if(ZZ==38) { // Strontium
      Z=ZZ;
      symbol="Sr";
      name="Strontium";
      period=5;
      group=2;
      series="AlkalineEarthMetal";
      block="s";
      mass=AMU2KILOGRAM*87.62;
      molar_volume=0.000033316;
      volume=55.4105;
      Vm_Miedema=10.2;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      density_PT=2.63;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=6.09969E-05;
      lattice_constants[1]=608.49;lattice_constants[2]=608.49;lattice_constants[3]=608.49;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.215;
      radius_PT=219;
      radius_covalent=1.95;
      radius_covalent_PT=195;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.9709;
      radii_Slatter=2.00;
      radii_Pyykko=1.85;
      conductivity_electrical=7.7E6;
      electronegativity_vec=0.95;
      hardness_Ghosh=2.5374;
      electronegativity_Pearson=2.0;
      electronegativity_Ghosh=3.225;
      electronegativity_Allen=0.963;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=5.03;
      phi_star_Miedema=2.40;
      nws_Miedema=0.84;
      gamma_s_Miedema=430;
      scale_Pettifor=0.55;
      temperature_boiling=1382;
      temperature_melting=777;
      vaporization_heat_PT=137;
      specific_heat_PT=300;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000225;
      conductivity_thermal=35;
      hardness_Brinell=NNN;
      hardness_Mohs=1.5;
      hardness_Vickers=NNN;
      hardness_Pearson=3.70;
      hardness_Putz=0.11;
      hardness_RB=3.08;
      modulus_shear=6.1;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=0.28;
      BVm_Miedema=3.9;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.32E-9;
      susceptibility_magnetic_volume=3.47E-6;
      susceptibility_magnetic_molar=1.16E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=4200;
      HHIR=3000;
      /*xray_scatt=NNN;*/
      //Sr
      return; //CO20200520
    }
    // [AFLOW]STOP=Strontium
    // ********************************************************************************************************************************************************

    // d-electron systems: transition metals
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Yttrium
    // Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium Yttrium
    else if(ZZ==39) { // Yttrium
      Z=ZZ;
      symbol="Y";
      name="Yttrium";
      period=5;
      group=3;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*88.9059;
      molar_volume=0.000019881;
      volume=32.4546;
      Vm_Miedema=7.3;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      density_PT=4.472;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=364.74;lattice_constants[2]=364.74;lattice_constants[3]=573.06;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.181;
      radius_PT=212;
      radius_covalent=1.90;
      radius_covalent_PT=190;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.8224;
      radii_Slatter=1.80;
      radii_Pyykko=1.63;
      conductivity_electrical=1.8E6;
      electronegativity_vec=1.22;
      hardness_Ghosh=2.6335;
      electronegativity_Pearson=3.19;
      electronegativity_Ghosh=3.311;
      electronegativity_Allen=1.12; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(3);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);  //RF+SK20200410
      electron_affinity_PT=29.6;
      phi_star_Miedema=3.20;
      nws_Miedema=1.21;
      gamma_s_Miedema=1100;
      scale_Pettifor=0.70;
      temperature_boiling=3345;
      temperature_melting=1526;
      vaporization_heat_PT=380;
      specific_heat_PT=298;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000106;
      conductivity_thermal=17;
      hardness_Brinell=588;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=3.19;
      hardness_Putz=0.14;
      hardness_RB=3.67;
      modulus_shear=26;
      modulus_Young=64;
      modulus_bulk=41;
      Poisson_ratio_PT=0.24;
      BVm_Miedema=7.2;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=6.66E-8;
      susceptibility_magnetic_volume=0.0002978;
      susceptibility_magnetic_molar=5.921E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9800;
      HHIR=2600;
      /*xray_scatt=NNN;*/
      //Y
      return; //CO20200520
    }
    // [AFLOW]STOP=Yttrium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Zirconium
    // Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium Zirconium
    else if(ZZ==40) { // Zirconium
      Z=ZZ;
      symbol="Zr";
      name="Zirconium";
      period=5;
      group=4;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*91.22;
      molar_volume=0.000014011;
      volume=23.2561;
      Vm_Miedema=5.8;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      density_PT=6.511;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.000342629;
      lattice_constants[1]=323.2;lattice_constants[2]=323.2;lattice_constants[3]=514.7;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.158;
      radius_PT=206;
      radius_covalent=1.75;
      radius_covalent_PT=175;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.688;
      radii_Slatter=1.55;
      radii_Pyykko=1.54;
      conductivity_electrical=2.4E6;
      electronegativity_vec=1.33;
      hardness_Ghosh=2.7298;
      electronegativity_Pearson=3.64;
      electronegativity_Ghosh=3.398;
      electronegativity_Allen=1.32; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(3); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4);  //RF+SK20200410
      electron_affinity_PT=41.1;
      phi_star_Miedema=3.40;
      nws_Miedema=1.39;
      gamma_s_Miedema=1950;
      scale_Pettifor=0.76;
      temperature_boiling=4409;
      temperature_melting=1855;
      vaporization_heat_PT=580;
      specific_heat_PT=278;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=5.7E-6;
      conductivity_thermal=23;
      hardness_Brinell=650;
      hardness_Mohs=5;
      hardness_Vickers=904;
      hardness_Pearson=3.21;
      hardness_Putz=0.17;
      hardness_RB=2.09;
      modulus_shear=33;
      modulus_Young=67;
      modulus_bulk=NNN;
      Poisson_ratio_PT=0.34;
      BVm_Miedema=12.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.68E-8;
      susceptibility_magnetic_volume=0.000109;
      susceptibility_magnetic_molar=1.53E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3400;
      HHIR=2600;
      /*xray_scatt=NNN;*/
      //Zr
      return; //CO20200520
    }
    // [AFLOW]STOP=Zirconium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Niobium
    // Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium Niobium
    else if(ZZ==41) { // Niobium
      Z=ZZ;
      symbol="Nb";
      name="Niobium";
      period=5;
      group=5;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*92.9064;
      molar_volume=0.000010841;
      volume=18.3132;
      Vm_Miedema=4.9;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      density_PT=8.57;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=0.0;
      lattice_constants[1]=330.04;lattice_constants[2]=330.04;lattice_constants[3]=330.04;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.143;
      radius_PT=198;
      radius_covalent=1.64;
      radius_covalent_PT=164;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.5658;
      radii_Slatter=1.45;
      radii_Pyykko=1.47;
      conductivity_electrical=6.7E6;
      electronegativity_vec=1.60;
      hardness_Ghosh=2.8260;
      electronegativity_Pearson=4.0;
      electronegativity_Ghosh=3.485;
      electronegativity_Allen=1.41; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(5); oxidation_states.push_back(3); oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(5);  //RF+SK20200410
      electron_affinity_PT=86.1;
      phi_star_Miedema=4.00;
      nws_Miedema=1.62;
      gamma_s_Miedema=2700;
      scale_Pettifor=0.82;
      temperature_boiling=4744;
      temperature_melting=2477;
      vaporization_heat_PT=690;
      specific_heat_PT=265;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=7.3E-6;
      conductivity_thermal=54;
      hardness_Brinell=736;
      hardness_Mohs=6;
      hardness_Vickers=1320;
      hardness_Pearson=3.00;
      hardness_Putz=0.21;
      hardness_RB=3.67;
      modulus_shear=38;
      modulus_Young=105;
      modulus_bulk=170;
      Poisson_ratio_PT=0.4;
      BVm_Miedema=18.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=2.76E-8;
      susceptibility_magnetic_volume=0.000237;
      susceptibility_magnetic_molar=2.56E-9;
      temperature_Curie=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=8500;
      HHIR=8800;
      /*xray_scatt=NNN;*/
      //Nb
      return; //CO20200520
    }
    // [AFLOW]STOP=Niobium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Molybdenum
    // Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum Molybdenum
    else if(ZZ==42) { // Molybdenum
      Z=ZZ;
      symbol="Mo";
      name="Molybdenum";
      period=5;
      group=6;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*95.94;
      molar_volume=9.334E-6;
      volume=15.6175;
      Vm_Miedema=4.4;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      density_PT=10.28;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=0.000598128;
      lattice_constants[1]=314.7;lattice_constants[2]=314.7;lattice_constants[3]=314.7;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.136;
      radius_PT=190;
      radius_covalent=1.54;
      radius_covalent_PT=154;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.4543;
      radii_Slatter=1.45;
      radii_Pyykko=1.38;
      conductivity_electrical=2E7;
      electronegativity_vec=2.16;
      hardness_Ghosh=2.9221;
      electronegativity_Pearson=3.9;
      electronegativity_Ghosh=3.572;
      electronegativity_Allen=1.47; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(6); oxidation_states.push_back(5); oxidation_states.push_back(4); oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(0); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(6);  //RF+SK20200410
      electron_affinity_PT=71.9;
      phi_star_Miedema=4.65;
      nws_Miedema=1.77;
      gamma_s_Miedema=2950;
      scale_Pettifor=0.88;
      temperature_boiling=4639;
      temperature_melting=2623;
      vaporization_heat_PT=600;
      specific_heat_PT=251;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=4.8E-6;
      conductivity_thermal=139;
      hardness_Brinell=1500;
      hardness_Mohs=5.5;
      hardness_Vickers=1530;
      hardness_Pearson=3.10;
      hardness_Putz=0.25;
      hardness_RB=NNN;
      modulus_shear=20;
      modulus_Young=329;
      modulus_bulk=230;
      Poisson_ratio_PT=0.31;
      BVm_Miedema=26.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.17E-8;
      susceptibility_magnetic_volume=0.0001203;
      susceptibility_magnetic_molar=1.122E-9;
      temperature_Curie=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=2400;
      HHIR=5300;
      /*xray_scatt=NNN;*/
      //Mo
      return; //CO20200520
    }
    // [AFLOW]STOP=Molybdenum
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Technetium
    // Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium Technetium
    else if(ZZ==43) { // Technetium
      Z=ZZ;
      symbol="Tc";
      name="Technetium";
      period=5;
      group=7;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*98.9062;
      molar_volume=8.434782608696E-6;
      volume=14.4670;
      Vm_Miedema=4.2;
      valence_std=7;
      valence_iupac=7;
      valence_PT=6;
      density_PT=11.5;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=273.5;lattice_constants[2]=273.5;lattice_constants[3]=438.8;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=NNN;
      radius_PT=183;
      radius_covalent=1.47;
      radius_covalent_PT=147;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.352;
      radii_Slatter=1.35;
      radii_Pyykko=1.28;
      conductivity_electrical=5E6;
      electronegativity_vec=1.90;
      hardness_Ghosh=3.0184;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=3.659;
      electronegativity_Allen=1.51; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(7);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(7);  //RF+SK20200410
      electron_affinity_PT=53;
      phi_star_Miedema=5.30;
      nws_Miedema=1.81;
      gamma_s_Miedema=3050;
      scale_Pettifor=0.94;
      temperature_boiling=4265;
      temperature_melting=2157;
      vaporization_heat_PT=550;
      specific_heat_PT=63;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=51;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=0.29;
      hardness_RB=2.05;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=26.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=3.42E-8;
      susceptibility_magnetic_volume=0.0003933;
      susceptibility_magnetic_molar=3.352E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Tc JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Technetium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Ruthenium
    // Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium Ruthenium
    else if(ZZ==44) { // Ruthenium
      Z=ZZ;
      symbol="Ru";
      name="Ruthenium";
      period=5;
      group=8;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*101.07;
      molar_volume=8.1706E-6;
      volume=13.8390;
      Vm_Miedema=4.1;
      valence_std=8;
      valence_iupac=8;
      valence_PT=6;
      density_PT=12.37;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.000406665;
      lattice_constants[1]=270.59;lattice_constants[2]=270.59;lattice_constants[3]=428.15;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.134;
      radius_PT=178;
      radius_covalent=1.46;
      radius_covalent_PT=146;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.2579;
      radii_Slatter=1.30;
      radii_Pyykko=1.25;
      conductivity_electrical=1.4E7;
      electronegativity_vec=2.20;
      hardness_Ghosh=3.1146;
      electronegativity_Pearson=4.5;
      electronegativity_Ghosh=3.745;
      electronegativity_Allen=1.54; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(8); oxidation_states.push_back(6); oxidation_states.push_back(4); oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(0); oxidation_states.push_back(-2); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4); oxidation_states_preferred.push_back(3); //RF+SK20200410
      electron_affinity_PT=101.3;
      phi_star_Miedema=5.40;
      nws_Miedema=1.83;
      gamma_s_Miedema=3050;
      scale_Pettifor=1.00;
      temperature_boiling=4150;
      temperature_melting=2334;
      vaporization_heat_PT=580;
      specific_heat_PT=238;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.4E-6;
      conductivity_thermal=120;
      hardness_Brinell=2160;
      hardness_Mohs=6.5;
      hardness_Vickers=2298.138766667;
      hardness_Pearson=3.00;
      hardness_Putz=0.35;
      hardness_RB=NNN;
      modulus_shear=173;
      modulus_Young=447;
      modulus_bulk=220;
      Poisson_ratio_PT=0.3;
      BVm_Miedema=26.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=5.42E-9;
      susceptibility_magnetic_volume=0.000067;
      susceptibility_magnetic_molar=5.48E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3200;
      HHIR=8000;
      /*xray_scatt=NNN;*/
      //Ru JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Ruthenium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Rhodium
    // Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium Rhodium
    else if(ZZ==45) { // Rhodium
      Z=ZZ;
      symbol="Rh";
      name="Rhodium";
      period=5;
      group=9;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*102.9055;
      molar_volume=8.2655E-6;
      volume=14.1731;
      Vm_Miedema=4.1;
      valence_std=9;
      valence_iupac=6;
      valence_PT=6;
      density_PT=12.45;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=1.90706E-32;
      lattice_constants[1]=380.34;lattice_constants[2]=380.34;lattice_constants[3]=380.34;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.134;
      radius_PT=173;
      radius_covalent=1.42;
      radius_covalent_PT=142;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.1711;
      radii_Slatter=1.35;
      radii_Pyykko=1.25;
      conductivity_electrical=2.3E7;
      electronegativity_vec=2.28;
      hardness_Ghosh=3.2108;
      electronegativity_Pearson=4.3;
      electronegativity_Ghosh=3.832;
      electronegativity_Allen=1.56; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(5); oxidation_states.push_back(4); oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(1); oxidation_states.push_back(0); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3); oxidation_states_preferred.push_back(1); //RF+SK20200410
      electron_affinity_PT=109.7;
      phi_star_Miedema=5.40;
      nws_Miedema=1.76;
      gamma_s_Miedema=2750;
      scale_Pettifor=1.06;
      temperature_boiling=3695;
      temperature_melting=1964;
      vaporization_heat_PT=495;
      specific_heat_PT=240;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=8E-6;
      conductivity_thermal=150;
      hardness_Brinell=1100;
      hardness_Mohs=6;
      hardness_Vickers=1246;
      hardness_Pearson=3.16;
      hardness_Putz=0.41;
      hardness_RB=NNN;
      modulus_shear=150;
      modulus_Young=275;
      modulus_bulk=380;
      Poisson_ratio_PT=0.26;
      BVm_Miedema=23.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.36E-8;
      susceptibility_magnetic_volume=0.0001693;
      susceptibility_magnetic_molar=1.4E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3200;
      HHIR=8000;
      /*xray_scatt=NNN;*/
      //Rh
      return; //CO20200520
    }
    // [AFLOW]STOP=Rhodium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Palladium
    // Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium Palladium
    else if(ZZ==46) { // Palladium
      Z=ZZ;
      symbol="Pd";
      name="Palladium";
      period=5;
      group=10;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*106.4;
      molar_volume=8.8514E-6;
      volume=15.4596;
      Vm_Miedema=4.3;
      valence_std=10;
      valence_iupac=4;
      valence_PT=4;
      density_PT=12.023;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.000309478;
      lattice_constants[1]=389.07;lattice_constants[2]=389.07;lattice_constants[3]=389.07;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.137;
      radius_PT=169;
      radius_covalent=1.39;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=163;
      radii_Ghosh08=2.0907;
      radii_Slatter=1.40;
      radii_Pyykko=1.20;
      conductivity_electrical=1E7;
      electronegativity_vec=2.20;
      hardness_Ghosh=3.3069;
      electronegativity_Pearson=4.45;
      electronegativity_Ghosh=3.919;
      electronegativity_Allen=1.58; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(2); oxidation_states.push_back(0);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=53.7;
      phi_star_Miedema=5.45;
      nws_Miedema=1.67;
      gamma_s_Miedema=2100;
      scale_Pettifor=1.12;
      temperature_boiling=2963;
      temperature_melting=1554.9;
      vaporization_heat_PT=380;
      specific_heat_PT=240;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000118;
      conductivity_thermal=71;
      hardness_Brinell=37.2;
      hardness_Mohs=4.75;
      hardness_Vickers=461;
      hardness_Pearson=3.89;
      hardness_Putz=0.47;
      hardness_RB=6.32;
      modulus_shear=44;
      modulus_Young=121;
      modulus_bulk=180;
      Poisson_ratio_PT=0.39;
      BVm_Miedema=16.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=6.57E-8;
      susceptibility_magnetic_volume=0.0007899;
      susceptibility_magnetic_molar=6.992E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3200;
      HHIR=8000;
      /*xray_scatt=NNN;*/
      //Pd
      return; //CO20200520
    }
    // [AFLOW]STOP=Palladium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Silver
    // Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver Silver
    else if(ZZ==47) { // Silver
      Z=ZZ;
      symbol="Ag";
      name="Silver";
      period=5;
      group=11;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*107.8682;
      molar_volume=0.000010283;
      volume=18.0678;
      Vm_Miedema=4.7;
      valence_std=11;
      valence_iupac=4;
      valence_PT=1;
      density_PT=10.49;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=8.57985E-05;
      lattice_constants[1]=408.53;lattice_constants[2]=408.53;lattice_constants[3]=408.53;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.144;
      radius_PT=165;
      radius_covalent=1.45;
      radius_covalent_PT=145;
      radius_VanDerWaals_PT=172;
      radii_Ghosh08=2.016;
      radii_Slatter=1.60;
      radii_Pyykko=1.28;
      conductivity_electrical=6.2E7;
      electronegativity_vec=1.93;
      hardness_Ghosh=3.4032;
      electronegativity_Pearson=4.44;
      electronegativity_Ghosh=4.006;
      electronegativity_Allen=1.87; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2); oxidation_states.push_back(1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(1);  //RF+SK20200410
      electron_affinity_PT=125.6;
      phi_star_Miedema=4.45;
      nws_Miedema=1.39;
      gamma_s_Miedema=1250;
      scale_Pettifor=1.18;
      temperature_boiling=2162;
      temperature_melting=961.78;
      vaporization_heat_PT=255;
      specific_heat_PT=235;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000189;
      conductivity_thermal=430;
      hardness_Brinell=24.5;
      hardness_Mohs=2.5;
      hardness_Vickers=251;
      hardness_Pearson=3.14;
      hardness_Putz=0.55;
      hardness_RB=3.5;
      modulus_shear=30;
      modulus_Young=85;
      modulus_bulk=100;
      Poisson_ratio_PT=0.37;
      BVm_Miedema=10.0;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-2.27E-9;
      susceptibility_magnetic_volume=-0.0000238;
      susceptibility_magnetic_molar=-2.45E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1200;
      HHIR=1400;
      xray_scatt=47.18;
      //Ag JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Silver
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Cadmium
    // Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium Cadmium
    else if(ZZ==48) { // Cadmium
      Z=ZZ;
      symbol="Cd";
      name="Cadmium";
      period=5;
      group=12;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*112.41;
      molar_volume=0.000012996;
      volume=22.0408;
      Vm_Miedema=5.5;
      valence_std=12;
      valence_iupac=2;
      valence_PT=2;
      density_PT=8.65;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.000271603;
      lattice_constants[1]=297.94;lattice_constants[2]=297.94;lattice_constants[3]=561.86;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.150;
      radius_PT=161;
      radius_covalent=1.44;
      radius_covalent_PT=144;
      radius_VanDerWaals_PT=158;
      radii_Ghosh08=1.9465;
      radii_Slatter=1.55;
      radii_Pyykko=1.36;
      conductivity_electrical=1.4E7;
      electronegativity_vec=1.69;
      hardness_Ghosh=3.4994;
      electronegativity_Pearson=4.33;
      electronegativity_Ghosh=4.093;
      electronegativity_Allen=1.52; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=4.05;
      nws_Miedema=1.24;
      gamma_s_Miedema=780;
      scale_Pettifor=1.36;
      temperature_boiling=767;
      temperature_melting=321.07;
      vaporization_heat_PT=100;
      specific_heat_PT=230;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000308;
      conductivity_thermal=96;
      hardness_Brinell=203;
      hardness_Mohs=2;
      hardness_Vickers=NNN;
      hardness_Pearson=4.66;
      hardness_Putz=0.63;
      hardness_RB=5.35;
      modulus_shear=19;
      modulus_Young=50;
      modulus_bulk=42;
      Poisson_ratio_PT=0.3;
      BVm_Miedema=6.10;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-2.3E-9;
      susceptibility_magnetic_volume=-0.0000199;
      susceptibility_magnetic_molar=-2.59E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=1700;
      HHIR=1300;
      /*xray_scatt=NNN;*/
      //Cd
      return; //CO20200520
    }
    // [AFLOW]STOP=Cadmium
    // ********************************************************************************************************************************************************

    // p-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Indium
    // Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium Indium
    else if(ZZ==49) { // Indium
      Z=ZZ;
      symbol="In";
      name="Indium";
      period=5;
      group=13;
      series="PoorMetal";
      block="p";
      mass=AMU2KILOGRAM*114.82;
      molar_volume=0.000015707;
      volume=27.5233;
      Vm_Miedema=6.3;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      density_PT=7.31;
      crystal="fct";
      crystal_structure_PT="Centered_Tetragonal";
      space_group="I4/mmm";
      space_group_number=139;
      variance_parameter_mass=1.24494E-05;
      lattice_constants[1]=325.23;lattice_constants[2]=325.23;lattice_constants[3]=494.61;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.157;
      radius_PT=156;
      radius_covalent=1.42;
      radius_covalent_PT=142;
      radius_VanDerWaals_PT=193;
      radii_Ghosh08=1.6934;
      radii_Slatter=1.55;
      radii_Pyykko=1.42;
      conductivity_electrical=1.2E7;
      electronegativity_vec=1.78;
      hardness_Ghosh=3.9164;
      electronegativity_Pearson=3.1;
      electronegativity_Ghosh=4.469;
      electronegativity_Allen=1.656;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(3);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);  //RF+SK20200410
      electron_affinity_PT=28.9;
      phi_star_Miedema=3.90;
      nws_Miedema=1.17;
      gamma_s_Miedema=690;
      scale_Pettifor=1.60;
      temperature_boiling=2072;
      temperature_melting=156.6;
      vaporization_heat_PT=230;
      specific_heat_PT=233;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000321;
      conductivity_thermal=82;
      hardness_Brinell=8.83;
      hardness_Mohs=1.2;
      hardness_Vickers=NNN;
      hardness_Pearson=2.80;
      hardness_Putz=0.73;
      hardness_RB=2.77;
      modulus_shear=NNN;
      modulus_Young=11;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=6.4;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-1.4E-9;
      susceptibility_magnetic_volume=-0.0000102;
      susceptibility_magnetic_molar=-1.61E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3300;
      HHIR=2000;
      /*xray_scatt=NNN;*/
      //In
      return; //CO20200520
    }
    // [AFLOW]STOP=Indium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Tin
    // Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin Tin
    else if(ZZ==50) { // Tin
      Z=ZZ;
      symbol="Sn";
      name="Tin";
      period=5;
      group=14;
      series="PoorMetal";
      block="p";
      mass=AMU2KILOGRAM*118.69;
      molar_volume=0.000016239;
      volume=27.5555;
      Vm_Miedema=6.4;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      density_PT=7.31;
      crystal="bct";
      crystal_structure_PT="Centered_Tetragonal";
      space_group="I4_1/amd";
      space_group_number=141;
      variance_parameter_mass=0.000334085;
      lattice_constants[1]=583.18;lattice_constants[2]=583.18;lattice_constants[3]=318.19;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.158;
      radius_PT=145;
      radius_covalent=1.39;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=217;
      radii_Ghosh08=1.4986;
      radii_Slatter=1.45;
      radii_Pyykko=1.40;
      conductivity_electrical=9.1E6;
      electronegativity_vec=1.96;
      hardness_Ghosh=4.3332;
      electronegativity_Pearson=4.3;
      electronegativity_Ghosh=4.845;
      electronegativity_Allen=1.824;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(2); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4); oxidation_states_preferred.push_back(2); //RF+SK20200410
      electron_affinity_PT=107.3;
      phi_star_Miedema=4.15;
      nws_Miedema=1.24;
      gamma_s_Miedema=710;
      scale_Pettifor=1.84;
      temperature_boiling=2602;
      temperature_melting=231.93;
      vaporization_heat_PT=290;
      specific_heat_PT=217;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000022;
      conductivity_thermal=67;
      hardness_Brinell=51;
      hardness_Mohs=1.5;
      hardness_Vickers=NNN;
      hardness_Pearson=3.05;
      hardness_Putz=0.88;
      hardness_RB=3.15;
      modulus_shear=18;
      modulus_Young=50;
      modulus_bulk=58;
      Poisson_ratio_PT=0.36;
      BVm_Miedema=8.8;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-3.1E-9;
      susceptibility_magnetic_volume=-0.0000227;
      susceptibility_magnetic_molar=-3.68E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=2600;
      HHIR=1600;
      /*xray_scatt=NNN;*/
      //Sn
      return; //CO20200520
    }
    // [AFLOW]STOP=Tin
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Antimony
    // Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony Antimony
    else if(ZZ==51) { // Antimony
      Z=ZZ;
      symbol="Sb";
      name="Antimony";
      period=5;
      group=15;
      series="Metalloid";
      block="p";
      mass=AMU2KILOGRAM*121.75;
      molar_volume=0.000018181;
      volume=27.1823;
      Vm_Miedema=6.6;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      density_PT=6.697;
      crystal="rhl";
      crystal_structure_PT="Simple_Trigonal";
      space_group="R_3m";
      space_group_number=166;
      variance_parameter_mass=6.60751E-05;
      lattice_constants[1]=430.7;lattice_constants[2]=430.7;lattice_constants[3]=1127.3;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.161;
      radius_PT=133;
      radius_covalent=1.39;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.344;
      radii_Slatter=1.45;
      radii_Pyykko=1.40;
      conductivity_electrical=2.5E6;
      electronegativity_vec=2.05;
      hardness_Ghosh=4.7501;
      electronegativity_Pearson=4.85;
      electronegativity_Ghosh=5.221;
      electronegativity_Allen=1.984;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(5); oxidation_states.push_back(3); oxidation_states.push_back(-3); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);  //RF+SK20200410
      electron_affinity_PT=103.2;
      phi_star_Miedema=4.40;
      nws_Miedema=1.26;
      gamma_s_Miedema=680;
      scale_Pettifor=2.08;
      temperature_boiling=1587;
      temperature_melting=630.63;
      vaporization_heat_PT=67;
      specific_heat_PT=207;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000011;
      conductivity_thermal=24;
      hardness_Brinell=294;
      hardness_Mohs=3;
      hardness_Vickers=NNN;
      hardness_Pearson=3.80;
      hardness_Putz=1.10;
      hardness_RB=4.39;
      modulus_shear=20;
      modulus_Young=55;
      modulus_bulk=42;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=7.0;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-1.09E-8;
      susceptibility_magnetic_volume=-0.000073;
      susceptibility_magnetic_molar=-1.327E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=7900;
      HHIR=3400;
      /*xray_scatt=NNN;*/
      //Sb
      return; //CO20200520
    }
    // [AFLOW]STOP=Antimony
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Tellurium
    // Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium Tellurium
    else if(ZZ==52) { // Tellurium
      Z=ZZ;
      symbol="Te";
      name="Tellurium";
      period=5;
      group=16;
      series="Chalcogen";
      block="p";
      mass=AMU2KILOGRAM*127.6;
      molar_volume=0.000020449;
      volume=28.1993;
      Vm_Miedema=6.439;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      density_PT=6.24;
      crystal="hex";
      crystal_structure_PT="Simple_Trigonal";
      space_group="P3_121";
      space_group_number=152;
      variance_parameter_mass=0.000283934;
      lattice_constants[1]=445.72;lattice_constants[2]=445.72;lattice_constants[3]=592.9;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.143;
      radius_PT=123;
      radius_covalent=1.38;
      radius_covalent_PT=138;
      radius_VanDerWaals_PT=206;
      radii_Ghosh08=1.2183;
      radii_Slatter=1.40;
      radii_Pyykko=1.36;
      conductivity_electrical=10000;
      electronegativity_vec=2.10;
      hardness_Ghosh=5.1670;
      electronegativity_Pearson=5.49;
      electronegativity_Ghosh=5.597;
      electronegativity_Allen=2.158;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(6); oxidation_states.push_back(4); oxidation_states.push_back(-2); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4);  //RF+SK20200410
      electron_affinity_PT=190.2;
      phi_star_Miedema=4.72;
      nws_Miedema=1.31;
      gamma_s_Miedema=NNN;
      scale_Pettifor=2.32;
      temperature_boiling=988;
      temperature_melting=449.51;
      vaporization_heat_PT=48;
      specific_heat_PT=201;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=3;
      hardness_Brinell=180;
      hardness_Mohs=2.25;
      hardness_Vickers=NNN;
      hardness_Pearson=3.52;
      hardness_Putz=1.34;
      hardness_RB=3.47;
      modulus_shear=16;
      modulus_Young=43;
      modulus_bulk=64;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-3.9E-9;
      susceptibility_magnetic_volume=-0.0000243;
      susceptibility_magnetic_molar=-4.98E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=1.000991;
      HHIP=2900;
      HHIR=4900;
      /*xray_scatt=NNN;*/
      //Te Table 27 of JX
      return; //CO20200520
    }
    // [AFLOW]STOP=Tellurium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Iodine
    // Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine Iodine
    else if(ZZ==53) { // Iodine
      Z=ZZ;
      symbol="I";
      name="Iodine";
      period=5;
      group=17;
      series="Halogen";
      block="p";
      mass=AMU2KILOGRAM*126.9045;
      molar_volume=0.000025689;
      volume=34.9784;
      Vm_Miedema=8.72;
      valence_std=7;
      valence_iupac=7;
      valence_PT=7;
      density_PT=4.94;
      crystal="orc";
      crystal_structure_PT="Base_Orthorhombic";
      space_group="Cmca";
      space_group_number=64;
      variance_parameter_mass=0.0;
      lattice_constants[1]=718.02;lattice_constants[2]=471.02;lattice_constants[3]=981.03;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.136;
      radius_PT=115;
      radius_covalent=1.39;
      radius_covalent_PT=139;
      radius_VanDerWaals_PT=198;
      radii_Ghosh08=1.1141;
      radii_Slatter=1.40;
      radii_Pyykko=1.33;
      conductivity_electrical=1E-7;
      electronegativity_vec=2.66;
      hardness_Ghosh=5.5839;
      electronegativity_Pearson=6.76;
      electronegativity_Ghosh=5.973;
      electronegativity_Allen=2.359;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(7); oxidation_states.push_back(5); oxidation_states.push_back(1); oxidation_states.push_back(-1);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(-1); //RF+SK20200410
      electron_affinity_PT=295.2;
      phi_star_Miedema=5.33;
      nws_Miedema=0.17;
      gamma_s_Miedema=NNN;
      scale_Pettifor=2.56;
      temperature_boiling=184.3;
      temperature_melting=113.7;
      vaporization_heat_PT=20.9;
      specific_heat_PT=429;
      critical_Pressure=115.5;
      critical_Temperature_PT=819;
      thermal_expansion=NNN;
      conductivity_thermal=0.449;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=3.69;
      hardness_Putz=1.62;
      hardness_RB=3.81;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=7.7;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-4.5E-9;
      susceptibility_magnetic_volume=-0.0000222;
      susceptibility_magnetic_molar=-1.14E-9;
      temperature_Curie=NNN;
      color_PT="SLATEGRAY";
      refractive_index=NNN;
      HHIP=4900;
      HHIR=4800;
      /*xray_scatt=NNN;*/
      //I interpolation phi_star, nws, Vm,
      return; //CO20200520
    }
    // [AFLOW]STOP=Iodine
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Xenon
    // Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon Xenon
    else if(ZZ==54) { // Xenon
      Z=ZZ;
      symbol="Xe";
      name="Xenon";
      period=5;
      group=18;
      series="NobleGas";
      block="p";
      mass=AMU2KILOGRAM*131.3;
      molar_volume=0.0223;
      volume=-1.0000;
      Vm_Miedema=NNN;
      valence_std=0;
      valence_iupac=8;
      valence_PT=6;
      density_PT=59E-4;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.000267781;
      lattice_constants[1]=620.23;lattice_constants[2]=620.23;lattice_constants[3]=620.23;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Gas";
      radius=0.218;
      radius_PT=108;
      radius_covalent=1.40;
      radius_covalent_PT=140;
      radius_VanDerWaals_PT=216;
      radii_Ghosh08=1.0263;
      radii_Slatter=NNN;
      radii_Pyykko=1.31;
      conductivity_electrical=NNN;
      electronegativity_vec=2.60;
      hardness_Ghosh=6.0009;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.349;
      electronegativity_Allen=2.582;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(8); oxidation_states.push_back(6); oxidation_states.push_back(4); oxidation_states.push_back(2); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=-108;
      temperature_melting=-111.8;
      vaporization_heat_PT=12.64;
      specific_heat_PT=158.32;
      critical_Pressure=57.65;
      critical_Temperature_PT=289.77;
      thermal_expansion=NNN;
      conductivity_thermal=0.00565;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=1.92;
      hardness_RB=8.23;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-4.3E-9;
      susceptibility_magnetic_volume=-2.54E-8;
      susceptibility_magnetic_molar=-5.65E-10;
      temperature_Curie=NNN;
      color_PT="COLORLESS";
      refractive_index=1.000702;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Xe JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Xenon
    // ********************************************************************************************************************************************************

    // ROW6
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Cesium
    // Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium Cesium
    else if(ZZ==55) { // Cesium
      Z=ZZ;
      symbol="Cs";
      name="Cesium";
      period=6;
      group=1;
      series="AlkaliMetal";
      block="s";
      mass=AMU2KILOGRAM*132.9054;
      molar_volume=0.000070732;
      volume=117.281;
      Vm_Miedema=16.8;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      density_PT=1.879;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=0.0;
      lattice_constants[1]=614.1;lattice_constants[2]=614.1;lattice_constants[3]=614.1;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.265;
      radius_PT=298;
      radius_covalent=2.44;
      radius_covalent_PT=244;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=4.2433;
      radii_Slatter=2.60;
      radii_Pyykko=2.32;
      conductivity_electrical=5E6;
      electronegativity_vec=0.79;
      hardness_Ghosh=0.6829;
      electronegativity_Pearson=2.18;
      electronegativity_Ghosh=4.196;
      electronegativity_Allen=0.659;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(1);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(1);  //RF+SK20200410
      electron_affinity_PT=45.5;
      phi_star_Miedema=1.95;
      nws_Miedema=0.55;
      gamma_s_Miedema=95;
      scale_Pettifor=0.25;
      temperature_boiling=671;
      temperature_melting=28.44;
      vaporization_heat_PT=64;
      specific_heat_PT=242;
      critical_Pressure=92.77;
      critical_Temperature_PT=1938;
      thermal_expansion=NNN;
      conductivity_thermal=36;
      hardness_Brinell=0.14;
      hardness_Mohs=0.2;
      hardness_Vickers=NNN;
      hardness_Pearson=1.71;
      hardness_Putz=NNN;
      hardness_RB=1.98;
      modulus_shear=NNN;
      modulus_Young=1.7;
      modulus_bulk=1.6;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=1.4;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=2.8E-9;
      susceptibility_magnetic_volume=5.26E-6;
      susceptibility_magnetic_molar=3.72E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=6000;
      HHIR=6000;
      /*xray_scatt=NNN;*/
      //Cs
      return; //CO20200520
    }
    // [AFLOW]STOP=Cesium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Barium
    // Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium Barium
    else if(ZZ==56) { // Barium
      Z=ZZ;
      symbol="Ba";
      name="Barium";
      period=6;
      group=2;
      series="AlkalineEarthMetal";
      block="s";
      mass=AMU2KILOGRAM*137.33;
      molar_volume=0.000039125;
      volume=62.6649;
      Vm_Miedema=11.3;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      density_PT=3.51;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=6.23705E-05;
      lattice_constants[1]=502.8;lattice_constants[2]=502.8;lattice_constants[3]=502.8;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.217;
      radius_PT=253;
      radius_covalent=2.15;
      radius_covalent_PT=215;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=3.2753;
      radii_Slatter=2.15;
      radii_Pyykko=1.96;
      conductivity_electrical=2.9E6;
      electronegativity_vec=0.89;
      hardness_Ghosh=0.9201;
      electronegativity_Pearson=2.4;
      electronegativity_Ghosh=4.318;
      electronegativity_Allen=0.881;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=13.95;
      phi_star_Miedema=2.32;
      nws_Miedema=0.81;
      gamma_s_Miedema=370;
      scale_Pettifor=0.50;
      temperature_boiling=1870;
      temperature_melting=727;
      vaporization_heat_PT=140;
      specific_heat_PT=205;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000206;
      conductivity_thermal=18;
      hardness_Brinell=NNN;
      hardness_Mohs=1.25;
      hardness_Vickers=NNN;
      hardness_Pearson=2.90;
      hardness_Putz=NNN;
      hardness_RB=2.16;
      modulus_shear=4.9;
      modulus_Young=13;
      modulus_bulk=9.4;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=3.9;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.13E-8;
      susceptibility_magnetic_volume=0.00003966;
      susceptibility_magnetic_molar=1.552E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=3000;
      HHIR=2300;
      /*xray_scatt=NNN;*/
      //Ba
      return; //CO20200520
    }
    // [AFLOW]STOP=Barium
    // ********************************************************************************************************************************************************

    // d-electron systems: transition metals
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Lanthanium
    // Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium Lanthanium
    else if(ZZ==57) { // Lanthanium
      Z=ZZ;
      symbol="La";
      name="Lanthanium";
      period=6;
      group=3;  //[CO20200930]NNN
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*138.9055;
      molar_volume=0.000022601;
      volume=36.8495;
      Vm_Miedema=8.0;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      density_PT=6.146;
      crystal="hex";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=4.65323E-08;
      lattice_constants[1]=377.2;lattice_constants[2]=377.2;lattice_constants[3]=1214.4;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.187;
      radius_PT=NNN;
      radius_covalent=2.07;
      radius_covalent_PT=207;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.6673;
      radii_Slatter=1.95;
      radii_Pyykko=1.80;
      conductivity_electrical=1.6E6;
      electronegativity_vec=1.10;
      hardness_Ghosh=1.1571;
      electronegativity_Pearson=3.1;
      electronegativity_Ghosh=4.439;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(3);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=48;
      phi_star_Miedema=3.05;
      nws_Miedema=1.09;
      gamma_s_Miedema=900;
      scale_Pettifor=0.7480;
      temperature_boiling=3464;
      temperature_melting=919;
      vaporization_heat_PT=400;
      specific_heat_PT=195;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000121;
      conductivity_thermal=13;
      hardness_Brinell=363;
      hardness_Mohs=2.5;
      hardness_Vickers=491;
      hardness_Pearson=2.60;
      hardness_Putz=NNN;
      hardness_RB=2.46;
      modulus_shear=14;
      modulus_Young=37;
      modulus_bulk=28;
      Poisson_ratio_PT=0.28;
      BVm_Miedema=5.5;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.1E-8;
      susceptibility_magnetic_volume=0.00006761;
      susceptibility_magnetic_molar=1.528E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //La
      return; //CO20200520
    }
    // [AFLOW]STOP=Lanthanium
    // ********************************************************************************************************************************************************

    // lantanidies
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Cerium
    // Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium Cerium
    else if(ZZ==58) { // Cerium
      Z=ZZ;
      symbol="Ce";
      name="Cerium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*140.12;
      molar_volume=0.000020947;
      volume=26.4729;
      Vm_Miedema=7.76;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      density_PT=6.689;
      crystal="fcc";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=2.24956E-05;
      lattice_constants[1]=362;lattice_constants[2]=362;lattice_constants[3]=599;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.182;
      radius_PT=NNN;
      radius_covalent=2.04;
      radius_covalent_PT=204;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.2494;
      radii_Slatter=1.85;
      radii_Pyykko=1.63;
      conductivity_electrical=1.4E6;
      electronegativity_vec=1.12;
      hardness_Ghosh=1.3943;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.561;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(3);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.18;
      nws_Miedema=1.19;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7460;
      temperature_boiling=3360;
      temperature_melting=798;
      vaporization_heat_PT=350;
      specific_heat_PT=192;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.3E-6;
      conductivity_thermal=11;
      hardness_Brinell=412;
      hardness_Mohs=2.5;
      hardness_Vickers=270;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=1.8;
      modulus_shear=14;
      modulus_Young=34;
      modulus_bulk=22;
      Poisson_ratio_PT=0.24;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=2.2E-7;
      susceptibility_magnetic_volume=0.0014716;
      susceptibility_magnetic_molar=3.0826E-8;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Ce Pettifor linear interpolation// Miedema from Alonso-March.
      return; //CO20200520
    }
    // [AFLOW]STOP=Cerium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Praseodymium
    // Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium Praseodymium
    else if(ZZ==59) { // Praseodymium
      Z=ZZ;
      symbol="Pr";
      name="Praseodymium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*140.9077;
      molar_volume=0.000021221;
      volume=36.4987;
      Vm_Miedema=7.56;
      valence_std=5;
      valence_iupac=4;
      valence_PT=4;
      density_PT=6.64;
      crystal="hex";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=367.25;lattice_constants[2]=367.25;lattice_constants[3]=1183.54;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.183;
      radius_PT=247;
      radius_covalent=2.03;
      radius_covalent_PT=203;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.9447;
      radii_Slatter=1.85;
      radii_Pyykko=1.76;
      conductivity_electrical=1.4E6;
      electronegativity_vec=1.13;
      hardness_Ghosh=1.6315;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.682;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(3);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.19;
      nws_Miedema=1.20;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7440;
      temperature_boiling=3290;
      temperature_melting=931;
      vaporization_heat_PT=330;
      specific_heat_PT=193;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.7E-6;
      conductivity_thermal=13;
      hardness_Brinell=481;
      hardness_Mohs=1.41;
      hardness_Vickers=400;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=1.11;
      modulus_shear=15;
      modulus_Young=37;
      modulus_bulk=29;
      Poisson_ratio_PT=0.28;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=4.23E-7;
      susceptibility_magnetic_volume=0.0028087;
      susceptibility_magnetic_molar=5.9604E-8;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Pr Pettifor linear interpolation
      return; //CO20200520
    }
    // [AFLOW]STOP=Praseodymium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Neodymium
    // Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium Neodymium
    else if(ZZ==60) { // Neodymium
      Z=ZZ;
      symbol="Nd";
      name="Neodymium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*144.24;
      molar_volume=0.000020577;
      volume=29.6719;
      Vm_Miedema=7.51;
      valence_std=6;
      valence_iupac=4;
      valence_PT=3;
      density_PT=7.01;
      crystal="hex";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.000231599;
      lattice_constants[1]=365.8;lattice_constants[2]=365.8;lattice_constants[3]=1179.9;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.182;
      radius_PT=206;
      radius_covalent=2.01;
      radius_covalent_PT=201;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.7129;
      radii_Slatter=1.85;
      radii_Pyykko=1.74;
      conductivity_electrical=1.6E6;
      electronegativity_vec=1.14;
      hardness_Ghosh=1.8684;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.804;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(3);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.19;
      nws_Miedema=1.20;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7420;
      temperature_boiling=3100;
      temperature_melting=1021;
      vaporization_heat_PT=285;
      specific_heat_PT=190;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=9.6E-6;
      conductivity_thermal=17;
      hardness_Brinell=265;
      hardness_Mohs=1.23;
      hardness_Vickers=343;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=0.7;
      modulus_shear=16;
      modulus_Young=41;
      modulus_bulk=32;
      Poisson_ratio_PT=0.28;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=4.8E-7;
      susceptibility_magnetic_volume=0.0033648;
      susceptibility_magnetic_molar=6.9235E-8;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Nd Pettifor linear interpolation JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Neodymium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Promethium
    // Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium Promethium
    else if(ZZ==61) { // Promethium
      Z=ZZ;
      symbol="Pm";
      name="Promethium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*146.92;
      molar_volume=0.00001996145374449;
      volume=34.6133;
      Vm_Miedema=7.43;
      valence_std=7;
      valence_iupac=3;
      valence_PT=3;
      density_PT=7.264;
      crystal="hex";
      crystal_structure_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      variance_parameter_mass=0.0;
      lattice_constants[1]=NNN;lattice_constants[2]=NNN;lattice_constants[3]=NNN;
      lattice_angles[1]=NNN;lattice_angles[2]=NNN;lattice_angles[3]=NNN;
      phase="Solid";
      radius=NNN;
      radius_PT=205;
      radius_covalent=1.99;
      radius_covalent_PT=199;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.5303;
      radii_Slatter=1.85;
      radii_Pyykko=1.73;
      conductivity_electrical=1.3E6;
      electronegativity_vec=1.13;
      hardness_Ghosh=2.1056;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.925;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(3);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.19;
      nws_Miedema=1.21;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7400;
      temperature_boiling=3000;
      temperature_melting=1100;
      vaporization_heat_PT=290;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000011;
      conductivity_thermal=15;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=0.33;
      modulus_shear=18;
      modulus_Young=46;
      modulus_bulk=33;
      Poisson_ratio_PT=0.28;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      // Pm Pettifor linear interpolation
      return; //CO20200520
    }
    // [AFLOW]STOP=Promethium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Samarium
    // Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium Samarium
    else if(ZZ==62) { // Samarium
      Z=ZZ;
      symbol="Sm";
      name="Samarium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*150.4;
      molar_volume=0.000020449;
      volume=33.9484;
      Vm_Miedema=7.37;
      valence_std=8;
      valence_iupac=3;
      valence_PT=3;
      density_PT=7.353;
      crystal="rhl";
      crystal_structure_PT="Simple_Trigonal";
      space_group="R_3m";
      space_group_number=166;
      variance_parameter_mass=0.000334686;
      lattice_constants[1]=362.1;lattice_constants[2]=362.1;lattice_constants[3]=2625;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.181;
      radius_PT=238;
      radius_covalent=1.98;
      radius_covalent_PT=198;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.3830;
      radii_Slatter=1.85;
      radii_Pyykko=1.72;
      conductivity_electrical=1.1E6;
      electronegativity_vec=1.17;
      hardness_Ghosh=2.3427;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.047;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(3); oxidation_states.push_back(2);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.20;
      nws_Miedema=1.21;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7380;
      temperature_boiling=1803;
      temperature_melting=1072;
      vaporization_heat_PT=175;
      specific_heat_PT=196;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000127;
      conductivity_thermal=13;
      hardness_Brinell=441;
      hardness_Mohs=1.44;
      hardness_Vickers=412;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=0.02;
      modulus_shear=20;
      modulus_Young=50;
      modulus_bulk=38;
      Poisson_ratio_PT=0.27;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.11E-7;
      susceptibility_magnetic_volume=0.00081618;
      susceptibility_magnetic_molar=1.669E-8;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Sm Pettifor linear interpolation
      return; //CO20200520
    }
    // [AFLOW]STOP=Samarium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Europium
    // Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium Europium
    else if(ZZ==63) { // Europium
      Z=ZZ;
      symbol="Eu";
      name="Europium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*151.96;
      molar_volume=0.000028979;
      volume=43.1719;
      Vm_Miedema=7.36;
      valence_std=9;
      valence_iupac=3;
      valence_PT=3;
      density_PT=5.244;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=4.32857E-05;
      lattice_constants[1]=458.1;lattice_constants[2]=458.1;lattice_constants[3]=458.1;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.204;
      radius_PT=231;
      radius_covalent=1.98;
      radius_covalent_PT=198;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.2615;
      radii_Slatter=1.85;
      radii_Pyykko=1.68;
      conductivity_electrical=1.1E6;
      electronegativity_vec=1.20;
      hardness_Ghosh=2.5798;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.168;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(3); oxidation_states.push_back(2);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.20;
      nws_Miedema=1.21;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7360;
      temperature_boiling=1527;
      temperature_melting=822;
      vaporization_heat_PT=175;
      specific_heat_PT=182;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000035;
      conductivity_thermal=14;
      hardness_Brinell=NNN;
      hardness_Mohs=3.07;
      hardness_Vickers=167;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=2.42;
      modulus_shear=7.9;
      modulus_Young=18;
      modulus_bulk=8.3;
      Poisson_ratio_PT=0.15;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=2.76E-7;
      susceptibility_magnetic_volume=0.0014473;
      susceptibility_magnetic_molar=4.1942E-8;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Eu Pettifor linear interpolation
      return; //CO20200520
    }
    // [AFLOW]STOP=Europium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Gadolinium
    // Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium Gadolinium
    else if(ZZ==64) { // Gadolinium
      Z=ZZ;
      symbol="Gd";
      name="Gadolinium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*157.25;
      molar_volume=0.000019903;
      volume=32.5777;
      Vm_Miedema=7.34;
      valence_std=10;
      valence_iupac=3;
      valence_PT=3;
      density_PT=7.901;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.000127674;
      lattice_constants[1]=363.6;lattice_constants[2]=363.6;lattice_constants[3]=578.26;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.180;
      radius_PT=233;
      radius_covalent=1.96;
      radius_covalent_PT=196;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.1596;
      radii_Slatter=1.80;
      radii_Pyykko=1.69;
      conductivity_electrical=770000;
      electronegativity_vec=1.20;
      hardness_Ghosh=2.8170;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.290;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(3);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.20;
      nws_Miedema=1.21;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7340;
      temperature_boiling=3250;
      temperature_melting=1313;
      vaporization_heat_PT=305;
      specific_heat_PT=240;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=9.4E-6;
      conductivity_thermal=11;
      hardness_Brinell=NNN;
      hardness_Mohs=5.13;
      hardness_Vickers=570;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=-1.02;
      modulus_shear=22;
      modulus_Young=55;
      modulus_bulk=38;
      Poisson_ratio_PT=0.26;
      BVm_Miedema=NNN;
      magnetic_type_PT="Ferromagnetic";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=292;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      // Gd Pettifor linear interpolation
      return; //CO20200520
    }
    // [AFLOW]STOP=Gadolinium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Terbium
    // Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium Terbium
    else if(ZZ==65) { // Terbium
      Z=ZZ;
      symbol="Tb";
      name="Terbium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*158.9254;
      molar_volume=0.000019336;
      volume=32.0200;
      Vm_Miedema=7.20;
      valence_std=11;
      valence_iupac=4;
      valence_PT=3;
      density_PT=8.219;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=360.1;lattice_constants[2]=360.1;lattice_constants[3]=569.36;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.177;
      radius_PT=225;
      radius_covalent=1.94;
      radius_covalent_PT=194;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.0730;
      radii_Slatter=1.75;
      radii_Pyykko=1.68;
      conductivity_electrical=830000;
      electronegativity_vec=1.10;
      hardness_Ghosh=3.0540;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.411;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(3);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.21;
      nws_Miedema=1.22;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7320;
      temperature_boiling=3230;
      temperature_melting=1356;
      vaporization_heat_PT=295;
      specific_heat_PT=182;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000103;
      conductivity_thermal=11;
      hardness_Brinell=677;
      hardness_Mohs=2.33;
      hardness_Vickers=863;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=1.36;
      modulus_shear=22;
      modulus_Young=56;
      modulus_bulk=38.7;
      Poisson_ratio_PT=0.26;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=0.0000136;
      susceptibility_magnetic_volume=0.1117784;
      susceptibility_magnetic_molar=2.161385E-6;
      temperature_Curie=222;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      // Tb Pettifor linear interpolation
      return; //CO20200520
    }
    // [AFLOW]STOP=Terbium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Dysprosium
    // Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium Dysprosium
    else if(ZZ==66) { // Dysprosium
      Z=ZZ;
      symbol="Dy";
      name="Dysprosium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*162.5;
      molar_volume=0.000019004;
      volume=31.5096;
      Vm_Miedema=7.12;
      valence_std=12;
      valence_iupac=4;
      valence_PT=3;
      density_PT=8.551;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=5.20771E-05;
      lattice_constants[1]=359.3;lattice_constants[2]=359.3;lattice_constants[3]=565.37;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.177;
      radius_PT=228;
      radius_covalent=1.92;
      radius_covalent_PT=192;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.9984;
      radii_Slatter=1.75;
      radii_Pyykko=1.67;
      conductivity_electrical=1.1E6;
      electronegativity_vec=1.22;
      hardness_Ghosh=3.2912;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.533;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(3);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.21;
      nws_Miedema=1.22;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7300;
      temperature_boiling=2567;
      temperature_melting=1412;
      vaporization_heat_PT=280;
      specific_heat_PT=167;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.00001;
      conductivity_thermal=11;
      hardness_Brinell=500;
      hardness_Mohs=1.8;
      hardness_Vickers=540;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=1.06;
      modulus_shear=25;
      modulus_Young=61;
      modulus_bulk=41;
      Poisson_ratio_PT=0.25;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=5.45E-6;
      susceptibility_magnetic_volume=0.046603;
      susceptibility_magnetic_molar=8.85625E-7;
      temperature_Curie=87;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Dy Pettifor linear interpolation JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Dysprosium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Holmium
    // Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium Holmium
    else if(ZZ==67) { // Holmium
      Z=ZZ;
      symbol="Ho";
      name="Holmium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*164.9304;
      molar_volume=0.000018753;
      volume=31.0155;
      Vm_Miedema=7.06;
      valence_std=13;
      valence_iupac=3;
      valence_PT=3;
      density_PT=8.795;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=2.96961E-32;
      lattice_constants[1]=357.73;lattice_constants[2]=357.73;lattice_constants[3]=561.58;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.176;
      radius_PT=226;
      radius_covalent=1.92;
      radius_covalent_PT=192;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.9335;
      radii_Slatter=1.75;
      radii_Pyykko=1.66;
      conductivity_electrical=1.1E6;
      electronegativity_vec=1.23;
      hardness_Ghosh=3.5283;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.654;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(3);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.22;
      nws_Miedema=1.22;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7280;
      temperature_boiling=2700;
      temperature_melting=1474;
      vaporization_heat_PT=265;
      specific_heat_PT=165;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000112;
      conductivity_thermal=16;
      hardness_Brinell=746;
      hardness_Mohs=1.65;
      hardness_Vickers=481;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=0.78;
      modulus_shear=26;
      modulus_Young=64;
      modulus_bulk=40;
      Poisson_ratio_PT=0.23;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=5.49E-6;
      susceptibility_magnetic_volume=0.0482845;
      susceptibility_magnetic_molar=9.05467E-7;
      temperature_Curie=20;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Ho Pettifor linear interpolation
      return; //CO20200520
    }
    // [AFLOW]STOP=Holmium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Erbium
    // Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium Erbium
    else if(ZZ==68) { // Erbium
      Z=ZZ;
      symbol="Er";
      name="Erbium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*167.26;
      molar_volume=0.000018449;
      volume=30.5431;
      Vm_Miedema=6.98;
      valence_std=14;
      valence_iupac=3;
      valence_PT=3;
      density_PT=9.066;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=7.24618E-05;
      lattice_constants[1]=355.88;lattice_constants[2]=355.88;lattice_constants[3]=558.74;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.175;
      radius_PT=226;
      radius_covalent=1.89;
      radius_covalent_PT=189;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8765;
      radii_Slatter=1.75;
      radii_Pyykko=1.65;
      conductivity_electrical=1.2E6;
      electronegativity_vec=1.24;
      hardness_Ghosh=3.7655;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.776;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(3);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.22;
      nws_Miedema=1.23;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7260;
      temperature_boiling=2868;
      temperature_melting=1497;
      vaporization_heat_PT=285;
      specific_heat_PT=168;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000122;
      conductivity_thermal=15;
      hardness_Brinell=814;
      hardness_Mohs=1.97;
      hardness_Vickers=588;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=0.54;
      modulus_shear=28;
      modulus_Young=70;
      modulus_bulk=44;
      Poisson_ratio_PT=0.24;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=3.77E-6;
      susceptibility_magnetic_volume=0.0341788;
      susceptibility_magnetic_molar=6.30566E-7;
      temperature_Curie=32;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Er Pettifor linear interpolation
      return; //CO20200520
    }
    // [AFLOW]STOP=Erbium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Thulium
    // Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium Thulium
    else if(ZZ==69) { // Thulium
      Z=ZZ;
      symbol="Tm";
      name="Thulium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*168.9342;
      molar_volume=0.000018126;
      volume=30.0016;
      Vm_Miedema=6.90;
      valence_std=15;
      valence_iupac=4;
      valence_PT=3;
      density_PT=9.32;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=353.75;lattice_constants[2]=353.75;lattice_constants[3]=555.46;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.174;
      radius_PT=222;
      radius_covalent=1.90;
      radius_covalent_PT=190;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8261;
      radii_Slatter=1.75;
      radii_Pyykko=1.64;
      conductivity_electrical=1.4E6;
      electronegativity_vec=1.25;
      hardness_Ghosh=4.0026;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.897;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(3); oxidation_states.push_back(2);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.22;
      nws_Miedema=1.23;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7240;
      temperature_boiling=1950;
      temperature_melting=1545;
      vaporization_heat_PT=250;
      specific_heat_PT=160;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000133;
      conductivity_thermal=17;
      hardness_Brinell=471;
      hardness_Mohs=1.77;
      hardness_Vickers=520;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=0.32;
      modulus_shear=31;
      modulus_Young=74;
      modulus_bulk=45;
      Poisson_ratio_PT=0.21;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.99E-6;
      susceptibility_magnetic_volume=0.0185488;
      susceptibility_magnetic_molar=3.36179E-7;
      temperature_Curie=25;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Tm Pettifor linear interpolation JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Thulium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Ytterbium
    // Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium Ytterbium
    else if(ZZ==70) { // Ytterbium
      Z=ZZ;
      symbol="Yb";
      name="Ytterbium";
      period=6;
      group=NNN;
      series="Lanthanide";
      block="f";
      mass=AMU2KILOGRAM*173.04;
      molar_volume=0.000026339;
      volume=39.4395;
      Vm_Miedema=6.86;
      valence_std=16;
      valence_iupac=3;
      valence_PT=3;
      density_PT=6.57;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=8.54557E-05;
      lattice_constants[1]=548.47;lattice_constants[2]=548.47;lattice_constants[3]=548.47;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.193;
      radius_PT=222;
      radius_covalent=1.87;
      radius_covalent_PT=187;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.7812;
      radii_Slatter=1.75;
      radii_Pyykko=1.70;
      conductivity_electrical=3.6E6;
      electronegativity_vec=1.10;
      hardness_Ghosh=4.2395;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.019;
      electronegativity_Allen=1.09; // RF+SK20200410; use value for Lu since no values are available for other lanthanides and they are all usually very similar chemically
      oxidation_states.clear();oxidation_states.push_back(3); oxidation_states.push_back(2);
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);
      electron_affinity_PT=50;
      phi_star_Miedema=3.22;
      nws_Miedema=1.23;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7220;
      temperature_boiling=1196;
      temperature_melting=819;
      vaporization_heat_PT=160;
      specific_heat_PT=154;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000263;
      conductivity_thermal=39;
      hardness_Brinell=343;
      hardness_Mohs=NNN;
      hardness_Vickers=206;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=3.27;
      modulus_shear=10;
      modulus_Young=24;
      modulus_bulk=31;
      Poisson_ratio_PT=0.21;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=5.9E-9;
      susceptibility_magnetic_volume=0.0000388;
      susceptibility_magnetic_molar=1.02E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Yb Pettifor linear interpolation
      return; //CO20200520
    }
    // [AFLOW]STOP=Ytterbium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Lutetium
    // Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium Lutetium
    else if(ZZ==71) { // Lutetium
      Z=ZZ;
      symbol="Lu";
      name="Lutetium";
      period=6;
      group=3;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*174.967;
      molar_volume=0.000017779;
      volume=29.3515;
      Vm_Miedema=6.81;
      valence_std=17;
      valence_iupac=3;
      valence_PT=3;
      density_PT=9.841;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=8.27273E-07;
      lattice_constants[1]=350.31;lattice_constants[2]=350.31;lattice_constants[3]=555.09;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.173;
      radius_PT=217;
      radius_covalent=1.87;
      radius_covalent_PT=187;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.7409;
      radii_Slatter=1.75;
      radii_Pyykko=1.62;
      conductivity_electrical=1.8E6;
      electronegativity_vec=1.27;
      hardness_Ghosh=4.4766;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.140;
      electronegativity_Allen=1.09; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(3);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);  //RF+SK20200410
      electron_affinity_PT=50;
      phi_star_Miedema=3.22;
      nws_Miedema=1.24;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0.7200;
      temperature_boiling=3402;
      temperature_melting=1663;
      vaporization_heat_PT=415;
      specific_heat_PT=154;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.00001;
      conductivity_thermal=16;
      hardness_Brinell=893;
      hardness_Mohs=2.6;
      hardness_Vickers=1160;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=3.64;
      modulus_shear=27;
      modulus_Young=67;
      modulus_bulk=48;
      Poisson_ratio_PT=0.26;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.2E-9;
      susceptibility_magnetic_volume=0.0000118;
      susceptibility_magnetic_molar=2.1E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=9500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Lu
      return; //CO20200520
    }
    // [AFLOW]STOP=Lutetium
    // ********************************************************************************************************************************************************

    // d-electron systems: transition metals
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Hafnium
    // Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium Hafnium
    else if(ZZ==72) { // Hafnium
      Z=ZZ;
      symbol="Hf";
      name="Hafnium";
      period=6;
      group=4;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*178.49;
      molar_volume=0.0000134102;
      volume=22.0408;
      Vm_Miedema=5.6;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      density_PT=13.31;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=5.25384E-05;
      lattice_constants[1]=319.64;lattice_constants[2]=319.64;lattice_constants[3]=505.11;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.159;
      radius_PT=208;
      radius_covalent=1.75;
      radius_covalent_PT=175;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.7056;
      radii_Slatter=1.55;
      radii_Pyykko=1.52;
      conductivity_electrical=3.3E6;
      electronegativity_vec=1.30;
      hardness_Ghosh=4.7065;
      electronegativity_Pearson=3.8;
      electronegativity_Ghosh=6.258;
      electronegativity_Allen=1.16; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(3); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=3.55;
      nws_Miedema=1.43;
      gamma_s_Miedema=2200;
      scale_Pettifor=0.775;
      temperature_boiling=4603;
      temperature_melting=2233;
      vaporization_heat_PT=630;
      specific_heat_PT=144;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=5.9E-6;
      conductivity_thermal=23;
      hardness_Brinell=1700;
      hardness_Mohs=5.5;
      hardness_Vickers=1760;
      hardness_Pearson=3.00;
      hardness_Putz=NNN;
      hardness_RB=3.94;
      modulus_shear=30;
      modulus_Young=78;
      modulus_bulk=110;
      Poisson_ratio_PT=0.37;
      BVm_Miedema=15.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=5.3E-9;
      susceptibility_magnetic_volume=0.0000705;
      susceptibility_magnetic_molar=9.46E-10;
      temperature_Curie=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=3400;
      HHIR=2600;
      /*xray_scatt=NNN;*/
      //Hf
      return; //CO20200520
    }
    // [AFLOW]STOP=Hafnium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Tantalum
    // Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum Tantalum
    else if(ZZ==73) { // Tantalum
      Z=ZZ;
      symbol="Ta";
      name="Tantalum";
      period=6;
      group=5;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*180.9479;
      molar_volume=0.0000108677;
      volume=18.1100;
      Vm_Miedema=4.9;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      density_PT=16.65;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=3.66845E-09;
      lattice_constants[1]=330.13;lattice_constants[2]=330.13;lattice_constants[3]=330.13;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.147;
      radius_PT=200;
      radius_covalent=1.70;
      radius_covalent_PT=170;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.6716;
      radii_Slatter=1.45;
      radii_Pyykko=1.46;
      conductivity_electrical=7.7E6;
      electronegativity_vec=1.50;
      hardness_Ghosh=4.9508;
      electronegativity_Pearson=4.11;
      electronegativity_Ghosh=6.383;
      electronegativity_Allen=1.34; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(5); oxidation_states.push_back(3); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(5);  //RF+SK20200410
      electron_affinity_PT=31;
      phi_star_Miedema=4.05;
      nws_Miedema=1.63;
      gamma_s_Miedema=3050;
      scale_Pettifor=0.83;
      temperature_boiling=5458;
      temperature_melting=3017;
      vaporization_heat_PT=736;
      specific_heat_PT=140;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.3E-6;
      conductivity_thermal=57;
      hardness_Brinell=800;
      hardness_Mohs=6.5;
      hardness_Vickers=873;
      hardness_Pearson=3.79;
      hardness_Putz=NNN;
      hardness_RB=1.75;
      modulus_shear=67;
      modulus_Young=186;
      modulus_bulk=200;
      Poisson_ratio_PT=0.34;
      BVm_Miedema=22.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.07E-8;
      susceptibility_magnetic_volume=0.0001782;
      susceptibility_magnetic_molar=1.936E-9;
      temperature_Curie=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=2300;
      HHIR=4800;
      /*xray_scatt=NNN;*/
      //Ta
      return; //CO20200520
    }
    // [AFLOW]STOP=Tantalum
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Tungsten
    // Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten Tungsten
    else if(ZZ==74) { // Tungsten
      Z=ZZ;
      symbol="W";
      name="Tungsten";
      period=6;
      group=6;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*183.85;
      molar_volume=9.5501E-6;
      volume=15.9387;
      Vm_Miedema=4.5;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      density_PT=19.25;
      crystal="bcc";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=6.96679E-05;
      lattice_constants[1]=316.52;lattice_constants[2]=316.52;lattice_constants[3]=316.52;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.137;
      radius_PT=193;
      radius_covalent=1.62;
      radius_covalent_PT=162;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.6416;
      radii_Slatter=1.35;
      radii_Pyykko=1.37;
      conductivity_electrical=2E7;
      electronegativity_vec=2.36;
      hardness_Ghosh=5.1879;
      electronegativity_Pearson=4.40;
      electronegativity_Ghosh=6.505;
      electronegativity_Allen=1.47; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(6); oxidation_states.push_back(5); oxidation_states.push_back(4); oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(0); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(6);  //RF+SK20200410
      electron_affinity_PT=78.6;
      phi_star_Miedema=4.80;
      nws_Miedema=1.81;
      gamma_s_Miedema=3300;
      scale_Pettifor=0.885;
      temperature_boiling=5555;
      temperature_melting=3422;
      vaporization_heat_PT=800;
      specific_heat_PT=132;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=4.5E-6;
      conductivity_thermal=170;
      hardness_Brinell=2570;
      hardness_Mohs=7.5;
      hardness_Vickers=3430;
      hardness_Pearson=3.58;
      hardness_Putz=NNN;
      hardness_RB=1.23;
      modulus_shear=161;
      modulus_Young=411;
      modulus_bulk=310;
      Poisson_ratio_PT=0.28;
      BVm_Miedema=31.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=4.59E-9;
      susceptibility_magnetic_volume=0.0000884;
      susceptibility_magnetic_molar=8.44E-10;
      temperature_Curie=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=7000;
      HHIR=4300;
      /*xray_scatt=NNN;*/
      //W
      return; //CO20200520
    }
    // [AFLOW]STOP=Tungsten
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Rhenium
    // Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium Rhenium
    else if(ZZ==75) { // Rhenium
      Z=ZZ;
      symbol="Re";
      name="Rhenium";
      period=6;
      group=7;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*186.2;
      molar_volume=8.85856E-6;
      volume=14.8941;
      Vm_Miedema=4.3;
      valence_std=7;
      valence_iupac=7;
      valence_PT=7;
      density_PT=21.02;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=2.70849E-05;
      lattice_constants[1]=276.1;lattice_constants[2]=276.1;lattice_constants[3]=445.6;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.138;
      radius_PT=188;
      radius_covalent=1.51;
      radius_covalent_PT=151;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.6141;
      radii_Slatter=1.35;
      radii_Pyykko=1.31;
      conductivity_electrical=5.6E6;
      electronegativity_vec=1.90;
      hardness_Ghosh=5.4256;
      electronegativity_Pearson=4.02;
      electronegativity_Ghosh=6.626;
      electronegativity_Allen=1.60; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(7); oxidation_states.push_back(6); oxidation_states.push_back(4); oxidation_states.push_back(2); oxidation_states.push_back(-1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(7);  //RF+SK20200410
      electron_affinity_PT=14.5;
      phi_star_Miedema=5.40;
      nws_Miedema=1.86;
      gamma_s_Miedema=3650;
      scale_Pettifor=0.94;
      temperature_boiling=5596;
      temperature_melting=3186;
      vaporization_heat_PT=705;
      specific_heat_PT=137;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.2E-6;
      conductivity_thermal=48;
      hardness_Brinell=1320;
      hardness_Mohs=7;
      hardness_Vickers=2450;
      hardness_Pearson=3.87;
      hardness_Putz=NNN;
      hardness_RB=2.13;
      modulus_shear=178;
      modulus_Young=463;
      modulus_bulk=370;
      Poisson_ratio_PT=0.3;
      BVm_Miedema=33.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=4.56E-9;
      susceptibility_magnetic_volume=0.0000959;
      susceptibility_magnetic_molar=8.49E-10;
      temperature_Curie=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=3300;
      HHIR=3300;
      /*xray_scatt=NNN;*/
      //Re
      return; //CO20200520
    }
    // [AFLOW]STOP=Rhenium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Osmium
    // Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium Osmium
    else if(ZZ==76) { // Osmium
      Z=ZZ;
      symbol="Os";
      name="Osmium";
      period=6;
      group=8;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*190.2;
      molar_volume=8.421E-6;
      volume=14.2403;
      Vm_Miedema=4.2;
      valence_std=8;
      valence_iupac=8;
      valence_PT=6;
      density_PT=22.59;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=7.45234E-05;
      lattice_constants[1]=273.44;lattice_constants[2]=273.44;lattice_constants[3]=431.73;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.135;
      radius_PT=185;
      radius_covalent=1.44;
      radius_covalent_PT=144;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.5890;
      radii_Slatter=1.30;
      radii_Pyykko=1.29;
      conductivity_electrical=1.2E7;
      electronegativity_vec=2.20;
      hardness_Ghosh=5.6619;
      electronegativity_Pearson=4.9;
      electronegativity_Ghosh=6.748;
      electronegativity_Allen=1.65; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(8); oxidation_states.push_back(6); oxidation_states.push_back(4); oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(0); oxidation_states.push_back(-2); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4);  //RF+SK20200410
      electron_affinity_PT=106.1;
      phi_star_Miedema=5.40;
      nws_Miedema=1.85;
      gamma_s_Miedema=3500;
      scale_Pettifor=0.995;
      temperature_boiling=5012;
      temperature_melting=3033;
      vaporization_heat_PT=630;
      specific_heat_PT=130;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=5.1E-6;
      conductivity_thermal=87;
      hardness_Brinell=3920;
      hardness_Mohs=7;
      hardness_Vickers=4137.063913415;
      hardness_Pearson=3.80;
      hardness_Putz=NNN;
      hardness_RB=1.72;
      modulus_shear=222;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=0.25;
      BVm_Miedema=35.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=6E-10;
      susceptibility_magnetic_volume=0.000014;
      susceptibility_magnetic_molar=1.1E-10;
      temperature_Curie=NNN;
      color_PT="SLATEGRAY";
      refractive_index=NNN;
      HHIP=5500;
      HHIR=9100;
      /*xray_scatt=NNN;*/
      //Os JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Osmium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Iridium
    // Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium Iridium
    else if(ZZ==77) { // Iridium
      Z=ZZ;
      symbol="Ir";
      name="Iridium";
      period=6;
      group=9;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*192.22;
      molar_volume=8.5203E-6;
      volume=14.5561;
      Vm_Miedema=4.2;
      valence_std=9;
      valence_iupac=8;
      valence_PT=6;
      density_PT=22.56;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=2.53787E-05;
      lattice_constants[1]=383.9;lattice_constants[2]=383.9;lattice_constants[3]=383.9;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.135;
      radius_PT=180;
      radius_covalent=1.41;
      radius_covalent_PT=141;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.5657;
      radii_Slatter=1.35;
      radii_Pyykko=1.22;
      conductivity_electrical=2.1E7;
      electronegativity_vec=2.20;
      hardness_Ghosh=5.9000;
      electronegativity_Pearson=5.4;
      electronegativity_Ghosh=6.831;
      electronegativity_Allen=1.68; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(6); oxidation_states.push_back(4); oxidation_states.push_back(3); oxidation_states.push_back(2); oxidation_states.push_back(1); oxidation_states.push_back(0); oxidation_states.push_back(-1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4); oxidation_states_preferred.push_back(1); //RF+SK20200410
      electron_affinity_PT=151;
      phi_star_Miedema=5.55;
      nws_Miedema=1.83;
      gamma_s_Miedema=3100;
      scale_Pettifor=1.05;
      temperature_boiling=4428;
      temperature_melting=2466;
      vaporization_heat_PT=560;
      specific_heat_PT=131;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=6.4E-6;
      conductivity_thermal=150;
      hardness_Brinell=1670;
      hardness_Mohs=6.5;
      hardness_Vickers=1760;
      hardness_Pearson=3.80;
      hardness_Putz=NNN;
      hardness_RB=1.27;
      modulus_shear=210;
      modulus_Young=528;
      modulus_bulk=320;
      Poisson_ratio_PT=0.26;
      BVm_Miedema=25.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.67E-9;
      susceptibility_magnetic_volume=0.0000377;
      susceptibility_magnetic_molar=3.21E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=5500;
      HHIR=9100;
      /*xray_scatt=NNN;*/
      //Ir JX CHANGED VALENCE
      return; //CO20200520
    }
    // [AFLOW]STOP=Iridium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Platinum
    // Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum Platinum
    else if(ZZ==78) { // Platinum
      Z=ZZ;
      symbol="Pt";
      name="Platinum";
      period=6;
      group=10;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*195.09;
      molar_volume=9.0948E-6;
      volume=15.7298;
      Vm_Miedema=4.4;
      valence_std=10;
      valence_iupac=6;
      valence_PT=6;
      density_PT=21.45;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=3.39206E-05;
      lattice_constants[1]=392.42;lattice_constants[2]=392.42;lattice_constants[3]=392.42;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.138;
      radius_PT=177;
      radius_covalent=1.36;
      radius_covalent_PT=136;
      radius_VanDerWaals_PT=175;
      radii_Ghosh08=0.5443;
      radii_Slatter=1.35;
      radii_Pyykko=1.23;
      conductivity_electrical=9.4E6;
      electronegativity_vec=2.28;
      hardness_Ghosh=6.1367;
      electronegativity_Pearson=5.6;
      electronegativity_Ghosh=6.991;
      electronegativity_Allen=1.72; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(2); oxidation_states.push_back(0);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4); oxidation_states_preferred.push_back(2); //RF+SK20200410
      electron_affinity_PT=205.3;
      phi_star_Miedema=5.65;
      nws_Miedema=1.78;
      gamma_s_Miedema=2550;
      scale_Pettifor=1.105;
      temperature_boiling=3825;
      temperature_melting=1768.3;
      vaporization_heat_PT=490;
      specific_heat_PT=133;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=8.9E-6;
      conductivity_thermal=71;
      hardness_Brinell=392;
      hardness_Mohs=3.5;
      hardness_Vickers=549;
      hardness_Pearson=3.50;
      hardness_Putz=NNN;
      hardness_RB=3.5;
      modulus_shear=61;
      modulus_Young=168;
      modulus_bulk=230;
      Poisson_ratio_PT=0.38;
      BVm_Miedema=18.0;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=1.22E-8;
      susceptibility_magnetic_volume=0.0002573;
      susceptibility_magnetic_molar=2.38E-9;
      temperature_Curie=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=5500;
      HHIR=9100;
      /*xray_scatt=NNN;*/
      //Pt
      return; //CO20200520
    }
    // [AFLOW]STOP=Platinum
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Gold
    // Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold Gold
    else if(ZZ==79) { // Gold
      Z=ZZ;
      symbol="Au";
      name="Gold";
      period=6;
      group=11;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*196.9665;
      molar_volume=0.00001021;
      volume=18.1904;
      Vm_Miedema=4.7;
      valence_std=11;
      valence_iupac=5;
      valence_PT=5;
      density_PT=19.3;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=2.08217E-32;
      lattice_constants[1]=407.82;lattice_constants[2]=407.82;lattice_constants[3]=407.82;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.144;
      radius_PT=174;
      radius_covalent=1.36;
      radius_covalent_PT=136;
      radius_VanDerWaals_PT=166;
      radii_Ghosh08=0.5244;
      radii_Slatter=1.35;
      radii_Pyykko=1.24;
      conductivity_electrical=4.5E7;
      electronegativity_vec=2.54;
      hardness_Ghosh=6.3741;
      electronegativity_Pearson=5.77;
      electronegativity_Ghosh=7.112;
      electronegativity_Allen=1.92; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(3); oxidation_states.push_back(1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);  //RF+SK20200410
      electron_affinity_PT=222.8;
      phi_star_Miedema=5.15;
      nws_Miedema=1.57;
      gamma_s_Miedema=1550;
      scale_Pettifor=1.16;
      temperature_boiling=2856;
      temperature_melting=1064.18;
      vaporization_heat_PT=330;
      specific_heat_PT=129.1;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000142;
      conductivity_thermal=320;
      hardness_Brinell=2450;
      hardness_Mohs=2.5;
      hardness_Vickers=216;
      hardness_Pearson=3.46;
      hardness_Putz=NNN;
      hardness_RB=3.44;
      modulus_shear=27;
      modulus_Young=78;
      modulus_bulk=220;
      Poisson_ratio_PT=0.44;
      BVm_Miedema=18.0;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-1.78E-9;
      susceptibility_magnetic_volume=-0.0000344;
      susceptibility_magnetic_molar=-3.51E-10;
      temperature_Curie=NNN;
      color_PT="GOLD";
      refractive_index=NNN;
      HHIP=1100;
      HHIR=1000;
      xray_scatt=74.99;
      //Au
      return; //CO20200520
    }
    // [AFLOW]STOP=Gold
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Mercury
    // Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury Mercury
    else if(ZZ==80) { // Mercury
      Z=ZZ;
      symbol="Hg";
      name="Mercury";
      period=6;
      group=12;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*200.59;
      molar_volume=0.0000148213;
      volume=29.7156;
      Vm_Miedema=5.8;
      valence_std=12;
      valence_iupac=4;
      valence_PT=2;
      density_PT=13.534;
      crystal="rhl";
      crystal_structure_PT="Simple_Trigonal";
      space_group="R_3m";
      space_group_number=166;
      variance_parameter_mass=6.52519E-05;
      lattice_constants[1]=300.5;lattice_constants[2]=300.5;lattice_constants[3]=300.5;
      lattice_angles[1]=1.23081;lattice_angles[2]=1.23081;lattice_angles[3]=1.23081;
      phase="Liquid";
      radius=0.150;
      radius_PT=171;
      radius_covalent=1.32;
      radius_covalent_PT=132;
      radius_VanDerWaals_PT=155;
      radii_Ghosh08=0.5060;
      radii_Slatter=1.50;
      radii_Pyykko=1.33;
      conductivity_electrical=1E6;
      electronegativity_vec=2.00;
      hardness_Ghosh=6.6103;
      electronegativity_Pearson=4.91;
      electronegativity_Ghosh=7.233;
      electronegativity_Allen=1.76; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2); oxidation_states.push_back(1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=4.20;
      nws_Miedema=1.24;
      gamma_s_Miedema=610;
      scale_Pettifor=1.32;
      temperature_boiling=356.73;
      temperature_melting=-38.83;
      vaporization_heat_PT=59.2;
      specific_heat_PT=139.5;
      critical_Pressure=1698;
      critical_Temperature_PT=1750;
      thermal_expansion=0.000181;
      conductivity_thermal=8.3;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=5.54;
      hardness_Putz=NNN;
      hardness_RB=5.29;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=25;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=4.0;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-2.1E-9;
      susceptibility_magnetic_volume=-0.0000284;
      susceptibility_magnetic_molar=-4.21E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=1.000933;
      HHIP=5500;
      HHIR=3100;
      /*xray_scatt=NNN;*/
      //Hg
      return; //CO20200520
    }
    // [AFLOW]STOP=Mercury
    // ********************************************************************************************************************************************************

    // p-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Thallium
    // Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium Thallium
    else if(ZZ==81) { // Thallium
      Z=ZZ;
      symbol="Tl";
      name="Thallium";
      period=6;
      group=13;
      series="PoorMetal";
      block="p";
      mass=AMU2KILOGRAM*204.37;
      molar_volume=0.0000172473;
      volume=31.0721;
      Vm_Miedema=6.6;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      density_PT=11.85;
      crystal="hcp";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=1.99659E-05;
      lattice_constants[1]=345.66;lattice_constants[2]=345.66;lattice_constants[3]=552.48;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=0.171;
      radius_PT=156;
      radius_covalent=1.45;
      radius_covalent_PT=145;
      radius_VanDerWaals_PT=196;
      radii_Ghosh08=1.8670;
      radii_Slatter=1.90;
      radii_Pyykko=1.44;
      conductivity_electrical=6.7E6;
      electronegativity_vec=1.62;
      hardness_Ghosh=1.7043;
      electronegativity_Pearson=3.2;
      electronegativity_Ghosh=4.719;
      electronegativity_Allen=1.789;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(3); oxidation_states.push_back(1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(1);  //RF+SK20200410
      electron_affinity_PT=19.2;
      phi_star_Miedema=3.90;
      nws_Miedema=1.12;
      gamma_s_Miedema=610;
      scale_Pettifor=1.56;
      temperature_boiling=1473;
      temperature_melting=304;
      vaporization_heat_PT=165;
      specific_heat_PT=129;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000299;
      conductivity_thermal=46;
      hardness_Brinell=26.4;
      hardness_Mohs=1.2;
      hardness_Vickers=NNN;
      hardness_Pearson=2.90;
      hardness_Putz=NNN;
      hardness_RB=2.69;
      modulus_shear=2.8;
      modulus_Young=8;
      modulus_bulk=43;
      Poisson_ratio_PT=0.45;
      BVm_Miedema=6.2;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-3E-9;
      susceptibility_magnetic_volume=-0.0000356;
      susceptibility_magnetic_molar=-6.13E-10;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=6500;
      HHIR=6500;
      /*xray_scatt=NNN;*/
      //Tl
      return; //CO20200520
    }
    // [AFLOW]STOP=Thallium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Lead
    // Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead Lead
    else if(ZZ==82) { // Lead
      Z=ZZ;
      symbol="Pb";
      name="Lead";
      period=6;
      group=14;
      series="PoorMetal";
      block="p";
      mass=AMU2KILOGRAM*207.2;
      molar_volume=0.000018272;
      volume=31.6649;
      Vm_Miedema=6.9;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      density_PT=11.34;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=1.94378E-05;
      lattice_constants[1]=495.08;lattice_constants[2]=495.08;lattice_constants[3]=495.08;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.175;
      radius_PT=154;
      radius_covalent=1.46;
      radius_covalent_PT=146;
      radius_VanDerWaals_PT=202;
      radii_Ghosh08=1.6523;
      radii_Slatter=NNN;
      radii_Pyykko=1.44;
      conductivity_electrical=4.8E6;
      electronegativity_vec=2.33;
      hardness_Ghosh=1.9414;
      electronegativity_Pearson=3.90;
      electronegativity_Ghosh=4.841;
      electronegativity_Allen=1.854;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(4); oxidation_states.push_back(2); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=35.1;
      phi_star_Miedema=4.10;
      nws_Miedema=1.15;
      gamma_s_Miedema=610;
      scale_Pettifor=1.80;
      temperature_boiling=1749;
      temperature_melting=327.46;
      vaporization_heat_PT=178;
      specific_heat_PT=127;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000289;
      conductivity_thermal=35;
      hardness_Brinell=38.3;
      hardness_Mohs=1.5;
      hardness_Vickers=NNN;
      hardness_Pearson=3.53;
      hardness_Putz=NNN;
      hardness_RB=3.02;
      modulus_shear=5.6;
      modulus_Young=16;
      modulus_bulk=46;
      Poisson_ratio_PT=0.44;
      BVm_Miedema=7.9;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-1.5E-9;
      susceptibility_magnetic_volume=-0.000017;
      susceptibility_magnetic_molar=-3.11E-10;
      temperature_Curie=NNN;
      color_PT="SLATEGRAY";
      refractive_index=NNN;
      HHIP=2700;
      HHIR=1800;
      /*xray_scatt=NNN;*/
      //Pb
      return; //CO20200520
    }
    // [AFLOW]STOP=Lead
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Bismuth
    // Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth Bismuth
    else if(ZZ==83) { // Bismuth
      Z=ZZ;
      symbol="Bi";
      name="Bismuth";
      period=6;
      group=15;
      series="PoorMetal";
      block="p";
      mass=AMU2KILOGRAM*208.9804;
      molar_volume=0.000021368;
      volume=31.5691;
      Vm_Miedema=7.2;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      density_PT=9.78;
      crystal="rhl";
      crystal_structure_PT="Base-centered_Monoclinic";
      space_group="C12/m1";
      space_group_number=12;
      variance_parameter_mass=0.0;
      lattice_constants[1]=667.4;lattice_constants[2]=611.7;lattice_constants[3]=330.4;
      lattice_angles[1]=PI/2;lattice_angles[2]=1.925622;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.182;
      radius_PT=143;
      radius_covalent=1.48;
      radius_covalent_PT=148;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.4818;
      radii_Slatter=1.60;
      radii_Pyykko=1.51;
      conductivity_electrical=770000;
      electronegativity_vec=2.02;
      hardness_Ghosh=2.1785;
      electronegativity_Pearson=4.69;
      electronegativity_Ghosh=4.962;
      electronegativity_Allen=2.01; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(5); oxidation_states.push_back(3); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(3);  //RF+SK20200410
      electron_affinity_PT=91.2;
      phi_star_Miedema=4.15;
      nws_Miedema=1.16;
      gamma_s_Miedema=550;
      scale_Pettifor=2.04;
      temperature_boiling=1564;
      temperature_melting=271.3;
      vaporization_heat_PT=160;
      specific_heat_PT=122;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000134;
      conductivity_thermal=8;
      hardness_Brinell=94.2;
      hardness_Mohs=2.25;
      hardness_Vickers=NNN;
      hardness_Pearson=3.74;
      hardness_Putz=NNN;
      hardness_RB=4.14;
      modulus_shear=12;
      modulus_Young=32;
      modulus_bulk=31;
      Poisson_ratio_PT=0.33;
      BVm_Miedema=6.7;
      magnetic_type_PT="Diamagnetic";
      susceptibility_magnetic_mass=-1.7E-8;
      susceptibility_magnetic_volume=-0.00017;
      susceptibility_magnetic_molar=-3.6E-9;
      temperature_Curie=NNN;
      color_PT="GRAY";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Bi
      return; //CO20200520
    }
    // [AFLOW]STOP=Bismuth
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Polonium
    // Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium Polonium
    else if(ZZ==84) { // Polonium
      Z=ZZ;
      symbol="Po";
      name="Polonium";
      period=6;
      group=16;
      series="Chalcogen";
      block="p";
      mass=AMU2KILOGRAM*209.98;
      molar_volume=0.00002272727272727;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      density_PT=9.196;
      crystal="sc";
      crystal_structure_PT="Simple_Cubic";
      space_group="Pm-3m";
      space_group_number=221;
      variance_parameter_mass=0.0;
      lattice_constants[1]=335.9;lattice_constants[2]=335.9;lattice_constants[3]=335.9;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.140;
      radius_PT=135;
      radius_covalent=1.40;
      radius_covalent_PT=140;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.3431;
      radii_Slatter=1.90;
      radii_Pyykko=1.45;
      conductivity_electrical=2.3E6;
      electronegativity_vec=2.00;
      hardness_Ghosh=2.4158;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.084;
      electronegativity_Allen=2.19; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(6); oxidation_states.push_back(4); oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(4);  //RF+SK20200410
      electron_affinity_PT=183.3;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=2.28;
      temperature_boiling=962;
      temperature_melting=254;
      vaporization_heat_PT=100;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=NNN;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=3.28;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Po
      return; //CO20200520
    }
    // [AFLOW]STOP=Polonium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Astatine
    // Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine Astatine
    else if(ZZ==85) { // Astatine
      Z=ZZ;
      symbol="At";
      name="Astatine";
      period=6;
      group=17;
      series="Halogen";
      block="p";
      mass=AMU2KILOGRAM*210;
      molar_volume=NNN;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=7;
      valence_iupac=7;
      valence_PT=7;
      density_PT=NNN;
      crystal="nnn";
      crystal_structure_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      variance_parameter_mass=0.0;
      lattice_constants[1]=NNN;lattice_constants[2]=NNN;lattice_constants[3]=NNN;
      lattice_angles[1]=NNN;lattice_angles[2]=NNN;lattice_angles[3]=NNN;
      phase="Solid";
      radius=NNN;
      radius_PT=127;
      radius_covalent=1.50;
      radius_covalent_PT=150;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.2283;
      radii_Slatter=NNN;
      radii_Pyykko=1.47;
      conductivity_electrical=NNN;
      electronegativity_vec=2.20;
      hardness_Ghosh=2.6528;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.206;
      electronegativity_Allen=2.39; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(7); oxidation_states.push_back(5); oxidation_states.push_back(3); oxidation_states.push_back(1); oxidation_states.push_back(-1); //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(-1); //RF+SK20200410
      electron_affinity_PT=270.1;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=2.52;
      temperature_boiling=NNN;
      temperature_melting=302;
      vaporization_heat_PT=40;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=2;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=3.57;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //At
      return; //CO20200520
    }
    // [AFLOW]STOP=Astatine
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Radon
    // Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon Radon
    else if(ZZ==86) { // Radon
      Z=ZZ;
      symbol="Rn";
      name="Radon";
      period=6;
      group=18;
      series="NobleGas";
      block="p";
      mass=AMU2KILOGRAM*222;
      molar_volume=0.02281603288798;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=0;
      valence_iupac=6;
      valence_PT=6;
      density_PT=97.3E-4;
      crystal="fcc";
      crystal_structure_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      variance_parameter_mass=0.0;
      lattice_constants[1]=NNN;lattice_constants[2]=NNN;lattice_constants[3]=NNN;
      lattice_angles[1]=NNN;lattice_angles[2]=NNN;lattice_angles[3]=NNN;
      phase="Gas";
      radius=NNN;
      radius_PT=120;
      radius_covalent=1.50;
      radius_covalent_PT=150;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.1315;
      radii_Slatter=NNN;
      radii_Pyykko=1.42;
      conductivity_electrical=NNN;
      electronegativity_vec=2.2;
      hardness_Ghosh=2.8900;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.327;
      electronegativity_Allen=2.60; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=0;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=-61.7;
      temperature_melting=-71;
      vaporization_heat_PT=17;
      specific_heat_PT=93.65;
      critical_Pressure=61.98;
      critical_Temperature_PT=377;
      thermal_expansion=NNN;
      conductivity_thermal=0.00361;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=7.69;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="COLORLESS";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Rn
      return; //CO20200520
    }
    // [AFLOW]STOP=Radon
    // ********************************************************************************************************************************************************

    // ROW7
    // s-electron systems
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Francium
    // Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium Francium
    else if(ZZ==87) { // Francium
      Z=ZZ;
      symbol="Fr";
      name="Francium";
      period=7;
      group=1;
      series="AlkaliMetal";
      block="s";
      mass=AMU2KILOGRAM*223.02;
      molar_volume=NNN;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=1;
      valence_iupac=1;
      valence_PT=1;
      density_PT=NNN;
      crystal="bcc";
      crystal_structure_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      variance_parameter_mass=0.0;
      lattice_constants[1]=NNN;lattice_constants[2]=NNN;lattice_constants[3]=NNN;
      lattice_angles[1]=NNN;lattice_angles[2]=NNN;lattice_angles[3]=NNN;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=2.60;
      radius_covalent_PT=260;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=4.4479;
      radii_Slatter=NNN;
      radii_Pyykko=2.23;
      conductivity_electrical=NNN;
      electronegativity_vec=0.70;
      hardness_Ghosh=0.9882;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=2.376;
      electronegativity_Allen=0.67; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(1);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(1);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=NNN;
      temperature_melting=NNN;
      vaporization_heat_PT=64;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=NNN;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Fr
      return; //CO20200520
    }
    // [AFLOW]STOP=Francium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Radium
    // Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium Radium
    else if(ZZ==88) { // Radium
      Z=ZZ;
      symbol="Ra";
      name="Radium";
      period=7;
      group=2;
      series="AlkalineEarthMetal";
      block="s";
      mass=AMU2KILOGRAM*226.0254;
      molar_volume=0.0000452;
      volume=-1.0000;
      Vm_Miedema=NNN;
      valence_std=2;
      valence_iupac=2;
      valence_PT=2;
      density_PT=5;
      crystal="bct";
      crystal_structure_PT="Body-centered_Cubic";
      space_group="Im_3m";
      space_group_number=229;
      variance_parameter_mass=0.0;
      lattice_constants[1]=514.8;lattice_constants[2]=514.8;lattice_constants[3]=514.8;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=2.21;
      radius_covalent_PT=221;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=3.4332;
      radii_Slatter=2.15;
      radii_Pyykko=2.01;
      conductivity_electrical=1E6;
      electronegativity_vec=0.89;
      hardness_Ghosh=1.2819;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=2.664;
      electronegativity_Allen=0.89; //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(2);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(2);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=1737;
      temperature_melting=700;
      vaporization_heat_PT=125;
      specific_heat_PT=92;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=19;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Ra
      return; //CO20200520
    }
    // [AFLOW]STOP=Radium
    // ********************************************************************************************************************************************************

    // d-electron systems: transition metals
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Actinium
    // Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium Actinium
    else if(ZZ==89) { // Actinium
      Z=ZZ;
      symbol="Ac";
      name="Actinium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*227.03;
      molar_volume=0.00002254220456802;
      volume=45.2437;
      Vm_Miedema=NNN;
      valence_std=3;
      valence_iupac=3;
      valence_PT=3;
      density_PT=10.07;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.0;
      lattice_constants[1]=567;lattice_constants[2]=567;lattice_constants[3]=567;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=2.15;
      radius_covalent_PT=215;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=3.2615;
      radii_Slatter=1.95;
      radii_Pyykko=1.86;
      conductivity_electrical=NNN;
      electronegativity_vec=1.10;
      hardness_Ghosh=1.3497;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=2.730;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=3200;
      temperature_melting=1050;
      vaporization_heat_PT=400;
      specific_heat_PT=120;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=12;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Ac
      return; //CO20200520
    }
    // [AFLOW]STOP=Actinium
    // ********************************************************************************************************************************************************

    // actinidies
    // ********************************************************************************************************************************************************
    // [AFLOW]START=Thorium
    // Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium Thorium
    else if(ZZ==90) { // Thorium
      Z=ZZ;
      symbol="Th";
      name="Thorium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*232.0381;
      molar_volume=0.0000197917;
      volume=31.9586;
      Vm_Miedema=7.3;
      valence_std=4;
      valence_iupac=4;
      valence_PT=4;
      density_PT=11.724;
      crystal="fcc";
      crystal_structure_PT="Face-centered_Cubic";
      space_group="Fm_3m";
      space_group_number=225;
      variance_parameter_mass=0.0;
      lattice_constants[1]=508.42;lattice_constants[2]=508.42;lattice_constants[3]=508.42;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.180;
      radius_PT=NNN;
      radius_covalent=2.06;
      radius_covalent_PT=206;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=3.1061;
      radii_Slatter=1.80;
      radii_Pyykko=1.75;
      conductivity_electrical=6.7E6;
      electronegativity_vec=1.30;
      hardness_Ghosh=1.4175;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=2.796;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=3.30;
      nws_Miedema=1.28;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=4820;
      temperature_melting=1750;
      vaporization_heat_PT=530;
      specific_heat_PT=118;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.000011;
      conductivity_thermal=54;
      hardness_Brinell=400;
      hardness_Mohs=3;
      hardness_Vickers=350;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=31;
      modulus_Young=79;
      modulus_bulk=54;
      Poisson_ratio_PT=0.27;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=7.2E-9;
      susceptibility_magnetic_volume=0.000084;
      susceptibility_magnetic_molar=1.7E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      xray_scatt=86.64;
      //Th
      return; //CO20200520
    }
    // [AFLOW]STOP=Thorium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Protoactinium
    // Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium Protoactinium
    else if(ZZ==91) { // Protoactinium
      Z=ZZ;
      symbol="Pa";
      name="Protoactinium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*231.04;
      molar_volume=0.0000150316;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=5;
      valence_iupac=5;
      valence_PT=5;
      density_PT=15.37;
      crystal="bct";
      crystal_structure_PT="Centered_Tetragonal";
      space_group="I4/mmm";
      space_group_number=139;
      variance_parameter_mass=0.0;
      lattice_constants[1]=392.5;lattice_constants[2]=392.5;lattice_constants[3]=323.8;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=2.00;
      radius_covalent_PT=200;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=2.2756;
      radii_Slatter=1.80;
      radii_Pyykko=1.69;
      conductivity_electrical=5.6E6;
      electronegativity_vec=1.50;
      hardness_Ghosh=1.9369;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=3.306;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=4000;
      temperature_melting=1572;
      vaporization_heat_PT=470;
      specific_heat_PT=99.1;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=47;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=3.25E-8;
      susceptibility_magnetic_volume=0.0004995;
      susceptibility_magnetic_molar=7.509E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Pa
      return; //CO20200520
    }
    // [AFLOW]STOP=Protoactinium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Uranium
    // Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium Uranium
    else if(ZZ==92) { // Uranium
      Z=ZZ;
      symbol="U";
      name="Uranium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*238.03;
      molar_volume=0.000012495;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=6;
      valence_iupac=6;
      valence_PT=6;
      density_PT=19.05;
      crystal="orc";
      crystal_structure_PT="Base_Orthorhombic";
      space_group="Cmcm";
      space_group_number=63;
      variance_parameter_mass=1.15611E-06;
      lattice_constants[1]=285.37;lattice_constants[2]=586.95;lattice_constants[3]=495.48;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=0.138;
      radius_PT=NNN;
      radius_covalent=1.96;
      radius_covalent_PT=196;
      radius_VanDerWaals_PT=186;
      radii_Ghosh08=1.9767;
      radii_Slatter=1.75;
      radii_Pyykko=1.70;
      conductivity_electrical=3.6E6;
      electronegativity_vec=1.38;
      hardness_Ghosh=2.2306;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=3.594;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=3927;
      temperature_melting=1135;
      vaporization_heat_PT=420;
      specific_heat_PT=116;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=0.0000139;
      conductivity_thermal=27;
      hardness_Brinell=2400;
      hardness_Mohs=6;
      hardness_Vickers=1960;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=111;
      modulus_Young=208;
      modulus_bulk=100;
      Poisson_ratio_PT=0.23;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=2.16E-8;
      susceptibility_magnetic_volume=0.000411;
      susceptibility_magnetic_molar=5.14E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //U
      return; //CO20200520
    }
    // [AFLOW]STOP=Uranium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Neptunium
    // Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium Neptunium
    else if(ZZ==93) { // Neptunium
      Z=ZZ;
      symbol="Np";
      name="Neptunium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*237.05;
      molar_volume=0.00001158924205379;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=7;
      valence_iupac=7;
      valence_PT=6;
      density_PT=20.45;
      crystal="nnn";
      crystal_structure_PT="Simple_Orthorhombic";
      space_group="Pnma";
      space_group_number=62;
      variance_parameter_mass=0.0;
      lattice_constants[1]=666.3;lattice_constants[2]=472.3;lattice_constants[3]=488.7;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=1.90;
      radius_covalent_PT=190;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.7473;
      radii_Slatter=1.75;
      radii_Pyykko=1.71;
      conductivity_electrical=830000;
      electronegativity_vec=NNN;
      hardness_Ghosh=2.5241;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=3.882;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=4000;
      temperature_melting=644;
      vaporization_heat_PT=335;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=6;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Np
      return; //CO20200520
    }
    // [AFLOW]STOP=Neptunium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Plutonium
    // Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium Plutonium
    else if(ZZ==94) { // Plutonium
      Z=ZZ;
      symbol="Pu";
      name="Plutonium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*244.06;
      molar_volume=0.00001231328219621;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=8;
      valence_iupac=7;
      valence_PT=6;
      density_PT=19.816;
      crystal="nnn";
      crystal_structure_PT="Simple_Monoclinic";
      space_group="P12_1/m1";
      space_group_number=11;
      variance_parameter_mass=0.0;
      lattice_constants[1]=618.3;lattice_constants[2]=482.2;lattice_constants[3]=1096.3;
      lattice_angles[1]=PI/2;lattice_angles[2]=1.776571;lattice_angles[3]=PI/2;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=1.87;
      radius_covalent_PT=187;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.4496;
      radii_Slatter=1.75;
      radii_Pyykko=1.72;
      conductivity_electrical=670000;
      electronegativity_vec=NNN;
      hardness_Ghosh=3.0436;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.391;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=3230;
      temperature_melting=640;
      vaporization_heat_PT=325;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=6;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=43;
      modulus_Young=96;
      modulus_bulk=NNN;
      Poisson_ratio_PT=0.21;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=3.17E-8;
      susceptibility_magnetic_volume=0.0006282;
      susceptibility_magnetic_molar=7.735E-9;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Pu
      return; //CO20200520
    }
    // [AFLOW]STOP=Plutonium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Americium
    // Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium Americium
    else if(ZZ==95) { // Americium
      Z=ZZ;
      symbol="Am";
      name="Americium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*243.06;
      molar_volume=0.00001777615215801;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=9;
      valence_iupac=7;
      valence_PT=4;
      density_PT=13.67;
      crystal="nnn";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=346.81;lattice_constants[2]=346.81;lattice_constants[3]=1124.1;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=1.80;
      radius_covalent_PT=180;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.2915;
      radii_Slatter=1.75;
      radii_Pyykko=1.66;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=3.4169;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.678;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=2011;
      temperature_melting=1176;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=10;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="Paramagnetic";
      susceptibility_magnetic_mass=5.15E-8;
      susceptibility_magnetic_volume=0.000704;
      susceptibility_magnetic_molar=1.251E-8;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Am
      return; //CO20200520
    }
    // [AFLOW]STOP=Americium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Curium
    // Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium Curium
    else if(ZZ==96) { // Curium
      Z=ZZ;
      symbol="Cm";
      name="Curium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*247.07;
      molar_volume=0.00001828275351591;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=10;
      valence_iupac=8;
      valence_PT=4;
      density_PT=13.51;
      crystal="nnn";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=349.6;lattice_constants[2]=349.6;lattice_constants[3]=1133.1;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=1.69;
      radius_covalent_PT=169;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.2960;
      radii_Slatter=NNN;
      radii_Pyykko=1.66;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=3.4050;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=4.745;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=3110;
      temperature_melting=1345;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=NNN;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="SILVER";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Cm
      return; //CO20200520
    }
    // [AFLOW]STOP=Curium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Berkelium
    // Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium Berkelium
    else if(ZZ==97) { // Berkelium
      Z=ZZ;
      symbol="Bk";
      name="Berkelium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*247.07;
      molar_volume=0.00001671177266576;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=11;
      valence_iupac=4;
      valence_PT=4;
      density_PT=14.78;
      crystal="nnn";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=341.6;lattice_constants[2]=341.6;lattice_constants[3]=1106.9;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.1247;
      radii_Slatter=NNN;
      radii_Pyykko=1.68;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=3.9244;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.256;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=NNN;
      temperature_melting=1050;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=10;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Bk
      return; //CO20200520
    }
    // [AFLOW]STOP=Berkelium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Californium
    // Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium Californium
    else if(ZZ==98) { // Californium
      Z=ZZ;
      symbol="Cf";
      name="Californium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*251.08;
      molar_volume=0.00001662251655629;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=12;
      valence_iupac=4;
      valence_PT=4;
      density_PT=15.1;
      crystal="nnn";
      crystal_structure_PT="Simple_Hexagonal";
      space_group="P6_3/mmc";
      space_group_number=194;
      variance_parameter_mass=0.0;
      lattice_constants[1]=338;lattice_constants[2]=338;lattice_constants[3]=1102.5;
      lattice_angles[1]=PI/2;lattice_angles[2]=PI/2;lattice_angles[3]=2*PI/3;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=1.0465;
      radii_Slatter=NNN;
      radii_Pyykko=1.68;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=4.2181;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.542;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=NNN;
      temperature_melting=900;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=NNN;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Cf
      return; //CO20200520
    }
    // [AFLOW]STOP=Californium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Einsteinium
    // Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium Einsteinium
    else if(ZZ==99) { // Einsteinium
      Z=ZZ;
      symbol="Es";
      name="Einsteinium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*252.08;
      molar_volume=NNN;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=13;
      valence_iupac=4;
      valence_PT=4;
      density_PT=NNN;
      crystal="nnn";
      crystal_structure_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      variance_parameter_mass=0.0;
      lattice_constants[1]=NNN;lattice_constants[2]=NNN;lattice_constants[3]=NNN;
      lattice_angles[1]=NNN;lattice_angles[2]=NNN;lattice_angles[3]=NNN;
      phase="Solid";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.9785;
      radii_Slatter=NNN;
      radii_Pyykko=1.65;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=4.5116;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=5.830;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=NNN;
      temperature_melting=860;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=NNN;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Es
      return; //CO20200520
    }
    // [AFLOW]STOP=Einsteinium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Fermium
    // Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium Fermium
    else if(ZZ==100) { // Fermium
      Z=ZZ;
      symbol="Fm";
      name="Fermium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*257.1;
      molar_volume=NNN;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=14;
      valence_iupac=3;
      valence_PT=3;
      density_PT=NNN;
      crystal="nnn";
      crystal_structure_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      variance_parameter_mass=0.0;
      lattice_constants[1]=NNN;lattice_constants[2]=NNN;lattice_constants[3]=NNN;
      lattice_angles[1]=NNN;lattice_angles[2]=NNN;lattice_angles[3]=NNN;
      phase="nnn";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.9188;
      radii_Slatter=NNN;
      radii_Pyykko=1.67;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=4.8051;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.118;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=NNN;
      temperature_melting=1527;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=NNN;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Fm
      return; //CO20200520
    }
    // [AFLOW]STOP=Fermium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Mendelevium
    // Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium Mendelevium
    else if(ZZ==101) { // Mendelevium
      Z=ZZ;
      symbol="Md";
      name="Mendelevium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*258.1;
      molar_volume=NNN;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=15;
      valence_iupac=3;
      valence_PT=3;
      density_PT=NNN;
      crystal="nnn";
      crystal_structure_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      variance_parameter_mass=0.0;
      lattice_constants[1]=NNN;lattice_constants[2]=NNN;lattice_constants[3]=NNN;
      lattice_angles[1]=NNN;lattice_angles[2]=NNN;lattice_angles[3]=NNN;
      phase="nnn";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8659;
      radii_Slatter=NNN;
      radii_Pyykko=1.73;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=5.0990;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.406;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=NNN;
      temperature_melting=828;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=NNN;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Md
      return; //CO20200520
    }
    // [AFLOW]STOP=Mendelevium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Nobelium
    // Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium Nobelium
    else if(ZZ==102) { // Nobelium
      Z=ZZ;
      symbol="No";
      name="Nobelium";
      period=7;
      group=NNN;
      series="Actinide";
      block="f";
      mass=AMU2KILOGRAM*259.1;
      molar_volume=NNN;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=16;
      valence_iupac=3;
      valence_PT=3;
      density_PT=NNN;
      crystal="nnn";
      crystal_structure_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      variance_parameter_mass=0.0;
      lattice_constants[1]=NNN;lattice_constants[2]=NNN;lattice_constants[3]=NNN;
      lattice_angles[1]=NNN;lattice_angles[2]=NNN;lattice_angles[3]=NNN;
      phase="nnn";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8188;
      radii_Slatter=NNN;
      radii_Pyykko=1.76;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=5.3926;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.694;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=NNN;
      temperature_melting=828;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=NNN;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //No
      return; //CO20200520
    }
    // [AFLOW]STOP=Nobelium
    // ********************************************************************************************************************************************************

    // ********************************************************************************************************************************************************
    // [AFLOW]START=Lawrencium
    // Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium Lawrencium
    else if(ZZ==103) { // Lawrencium
      Z=ZZ;
      symbol="Lr";
      name="Lawrencium";
      period=7;
      group=3;
      series="TransitionMetal";
      block="d";
      mass=AMU2KILOGRAM*262.11;
      molar_volume=NNN;
      volume=NNN;
      Vm_Miedema=NNN;
      valence_std=17;
      valence_iupac=3;
      valence_PT=3;
      density_PT=NNN;
      crystal="nnn";
      crystal_structure_PT="NNN";
      space_group="NNN";
      space_group_number=NNN;
      variance_parameter_mass=0.0;
      lattice_constants[1]=NNN;lattice_constants[2]=NNN;lattice_constants[3]=NNN;
      lattice_angles[1]=NNN;lattice_angles[2]=NNN;lattice_angles[3]=NNN;
      phase="nnn";
      radius=NNN;
      radius_PT=NNN;
      radius_covalent=NNN;
      radius_covalent_PT=NNN;
      radius_VanDerWaals_PT=NNN;
      radii_Ghosh08=0.8086;
      radii_Slatter=NNN;
      radii_Pyykko=1.61;
      conductivity_electrical=NNN;
      electronegativity_vec=NNN;
      hardness_Ghosh=5.4607;
      electronegativity_Pearson=NNN;
      electronegativity_Ghosh=6.760;
      electronegativity_Allen=NNN;  //RF+SK20200410
      oxidation_states.clear();oxidation_states.push_back(NNN);  //RF+SK20200410
      oxidation_states_preferred.clear();oxidation_states_preferred.push_back(NNN);  //RF+SK20200410
      electron_affinity_PT=NNN;
      phi_star_Miedema=NNN;
      nws_Miedema=NNN;
      gamma_s_Miedema=NNN;
      scale_Pettifor=0;
      temperature_boiling=NNN;
      temperature_melting=1627;
      vaporization_heat_PT=NNN;
      specific_heat_PT=NNN;
      critical_Pressure=NNN;
      critical_Temperature_PT=NNN;
      thermal_expansion=NNN;
      conductivity_thermal=NNN;
      hardness_Brinell=NNN;
      hardness_Mohs=NNN;
      hardness_Vickers=NNN;
      hardness_Pearson=NNN;
      hardness_Putz=NNN;
      hardness_RB=NNN;
      modulus_shear=NNN;
      modulus_Young=NNN;
      modulus_bulk=NNN;
      Poisson_ratio_PT=NNN;
      BVm_Miedema=NNN;
      magnetic_type_PT="NNN";
      susceptibility_magnetic_mass=NNN;
      susceptibility_magnetic_volume=NNN;
      susceptibility_magnetic_molar=NNN;
      temperature_Curie=NNN;
      color_PT="NNN";
      refractive_index=NNN;
      HHIP=NNN;
      HHIR=NNN;
      /*xray_scatt=NNN;*/
      //Lr
      return; //CO20200520
    }
    // [AFLOW]STOP=Lawrencium
    // ********************************************************************************************************************************************************

    throw aurostd::xerror(_AFLOW_FILE_NAME_,"xelement::xelement():","Element number does not exist: "+aurostd::utype2string(ZZ),_FILE_NOT_FOUND_);  //CO20200520
  }
} // namespace xelement

#endif // _AFLOW_XELEMENT_CPP

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2020              *
// *                                                                        *
// **************************************************************************
