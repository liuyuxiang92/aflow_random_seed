// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu

#ifndef _AFLOW_ML_CPP_
#define _AFLOW_ML_CPP_

#include "aflow.h"
#include "aflow_cce.h"

#define _DEBUG_ML_ false

#define _ENERGIES_IONIZATION_MAX_ 5

namespace aflowML {
  double getLIB0EnergyPAWPBE(const string& species_pp){
    if(species_pp=="Ac") return -0.330411;
    if(species_pp=="Ac_s") return -0.985025;
    if(species_pp=="Ag") return -0.337763;
    if(species_pp=="Al") return -0.313186;
    if(species_pp=="Al_h") return -0.645972;
    if(species_pp=="Ar") return -0.051018;
    if(species_pp=="As") return -1.6833;
    if(species_pp=="Au") return -0.286055;
    if(species_pp=="B_h") return -0.444479;
    if(species_pp=="B_s") return -0.485521;
    if(species_pp=="Ba_sv") return -0.034562;
    if(species_pp=="Be") return -0.025499;
    if(species_pp=="Be_sv") return -0.017215;
    if(species_pp=="Bi") return -1.32564;
    if(species_pp=="Bi_d") return -1.51376;
    if(species_pp=="Br") return -0.261686;
    if(species_pp=="C") return -1.37179;
    if(species_pp=="C_h") return -1.3699;
    if(species_pp=="C_s") return -1.35209;
    if(species_pp=="Ca_pv") return -0.061212;
    if(species_pp=="Ca_sv") return -0.093014;
    if(species_pp=="Cd") return -0.166094;
    if(species_pp=="Cl") return -0.373905;
    if(species_pp=="Cl_h") return -0.35806;
    if(species_pp=="Co") return -1.99048;
    if(species_pp=="Cr") return -5.48369;
    if(species_pp=="Cr_pv") return -5.59721;
    if(species_pp=="Cs_sv") return -0.137636;
    if(species_pp=="Cu") return -0.241525;
    if(species_pp=="Cu_pv") return -0.599875;
    if(species_pp=="Dy_3") return -0.387234;
    if(species_pp=="Er_3") return -0.390044;
    if(species_pp=="Eu") return -8.36334;
    if(species_pp=="Eu_2") return -0.045949;
    if(species_pp=="F") return -0.704599;
    if(species_pp=="F_h") return -0.820093;
    if(species_pp=="F_s") return -0.663234;
    if(species_pp=="Fe") return -3.45752;
    if(species_pp=="Fe_pv") return -3.46819;
    if(species_pp=="Ga") return -0.277948;
    if(species_pp=="Ga_d") return -0.397361;
    if(species_pp=="Ga_h") return -0.282091;
    if(species_pp=="Ge") return -0.771177;
    if(species_pp=="Ge_d") return -0.885122;
    if(species_pp=="Ge_h") return -0.78229;
    if(species_pp=="H") return -1.11416;
    if(species_pp=="H_h") return -1.11179;
    if(species_pp=="He") return -0.000538;
    if(species_pp=="Hf") return -3.47402;
    if(species_pp=="Hf_pv") return -3.53246;
    if(species_pp=="Hg") return -0.122356;
    if(species_pp=="Ho_3") return -0.381068;
    if(species_pp=="I") return -0.212174;
    if(species_pp=="In") return -0.22805;
    if(species_pp=="In_d") return -0.407416;
    if(species_pp=="Ir") return -1.62733;
    if(species_pp=="K_pv") return -0.158564;
    if(species_pp=="K_sv") return -0.227694;
    if(species_pp=="La") return -0.5631;
    if(species_pp=="La_s") return -0.516331;
    if(species_pp=="Li") return -0.293945;
    if(species_pp=="Li_sv") return -0.297272;
    if(species_pp=="Lu") return -0.479887;
    if(species_pp=="Lu_3") return -0.382234;
    if(species_pp=="Mg") return -0.038135;
    if(species_pp=="Mg_pv") return -0.094786;
    if(species_pp=="Mn") return -5.15867;
    if(species_pp=="Mn_pv") return -5.32776;
    if(species_pp=="Mo") return -4.59003;
    if(species_pp=="Mo_pv") return -4.60432;
    if(species_pp=="N") return -3.12444;
    if(species_pp=="N_h") return -3.11812;
    if(species_pp=="N_s") return -3.13766;
    if(species_pp=="Na") return -0.21998;
    if(species_pp=="Na_pv") return -0.228379;
    if(species_pp=="Na_sv") return -0.235575;
    if(species_pp=="Nb_sv") return -3.21181;
    if(species_pp=="Nd") return -4.15081;
    if(species_pp=="Ne") return -0.012381;
    if(species_pp=="O") return -1.90885;
    if(species_pp=="O_h") return -1.96399;
    if(species_pp=="O_s") return -1.89585;
    if(species_pp=="Os") return -2.90915;
    if(species_pp=="Os_pv") return -2.92055;
    if(species_pp=="P") return -1.88731;
    if(species_pp=="P_h") return -1.88281;
    if(species_pp=="Pb") return -0.58418;
    if(species_pp=="Pb_d") return -0.765497;
    if(species_pp=="Pd") return -1.47315;
    if(species_pp=="Pd_pv") return -1.74771;
    if(species_pp=="Pm_3") return -0.499376;
    if(species_pp=="Pr") return -2.48756;
    if(species_pp=="Pr_3") return -0.558962;
    if(species_pp=="Pt") return -0.604572;
    if(species_pp=="Rb_pv") return -0.156737;
    if(species_pp=="Rb_sv") return -0.188051;
    if(species_pp=="Re") return -4.60422;
    if(species_pp=="Re_pv") return -4.69957;
    if(species_pp=="Rh") return -1.64118;
    if(species_pp=="Rh_pv") return -1.7322;
    if(species_pp=="Ru") return -2.52896;
    if(species_pp=="Ru_pv") return -2.57847;
    if(species_pp=="S") return -1.07971;
    if(species_pp=="S_h") return -1.06433;
    if(species_pp=="Sb") return -1.41051;
    if(species_pp=="Sc_sv") return -2.18958;
    if(species_pp=="Se") return -0.878847;
    if(species_pp=="Si") return -0.871529;
    if(species_pp=="Si_h") return -0.878813;
    if(species_pp=="Sm_3") return -0.478431;
    if(species_pp=="Sn") return -0.643434;
    if(species_pp=="Sn_d") return -0.832202;
    if(species_pp=="Sr_sv") return -0.068251;
    if(species_pp=="Ta") return -3.66606;
    if(species_pp=="Ta_pv") return -3.74451;
    if(species_pp=="Tc") return -3.44648;
    if(species_pp=="Tc_pv") return -3.53766;
    if(species_pp=="Te") return -0.730289;
    if(species_pp=="Th_s") return -0.948771;
    if(species_pp=="Ti") return -2.23219;
    if(species_pp=="Ti_pv") return -2.61834;
    if(species_pp=="Ti_sv") return -2.64228;
    if(species_pp=="Tl") return -0.196234;
    if(species_pp=="Tl_d") return -0.360034;
    if(species_pp=="Tm_3") return -0.287414;
    if(species_pp=="W") return -4.53997;
    if(species_pp=="W_pv") return -4.65963;
    if(species_pp=="Xe") return -0.008491;
    if(species_pp=="Y_sv") return -2.2892;
    if(species_pp=="Yb") return -0.026493;
    if(species_pp=="Yb_2") return -0.062556;
    if(species_pp=="Zn") return -0.163774;
    if(species_pp=="Zr") return -2.11254;
    if(species_pp=="Zr_sv") return -2.30965;
    return NNN;
  }

  double getLIB1EnergyPAWPBE(const string& species_pp){
    if(species_pp=="Ac") return -4.09393;
    if(species_pp=="Ac_s") return -4.03196;
    if(species_pp=="Ag") return -2.82746;
    if(species_pp=="Al") return -3.74356;
    if(species_pp=="Al_h") return -3.80021;
    if(species_pp=="Ar") return -0.065575;
    if(species_pp=="As") return -4.65243;
    if(species_pp=="Au") return -3.27212;
    if(species_pp=="Ba_sv") return -1.924;
    if(species_pp=="Be") return -3.75366;
    if(species_pp=="Be_sv") return -3.74102;
    if(species_pp=="Bi") return -3.87274;
    if(species_pp=="Bi_d") return -4.03716;
    if(species_pp=="C") return -9.22034;
    if(species_pp=="C_h") return -9.19568;
    if(species_pp=="C_s") return -9.19508;
    if(species_pp=="Ca_pv") return -1.97638;
    if(species_pp=="Ca_sv") return -2.00107;
    if(species_pp=="Cd") return -0.906233;
    if(species_pp=="Cl") return -1.78721;
    if(species_pp=="Cl_h") return -1.77694;
    if(species_pp=="Co") return -7.10876;
    if(species_pp=="Cr") return -9.51285;
    if(species_pp=="Cr_pv") return -9.6294;
    if(species_pp=="Cs_sv") return -0.852451;
    if(species_pp=="Cu") return -3.71935;
    if(species_pp=="Cu_pv") return -4.09701;
    if(species_pp=="Dy_3") return -4.58813;
    if(species_pp=="Eu") return -10.2377;
    if(species_pp=="F") return -1.85953;
    if(species_pp=="F_h") return -1.90899;
    if(species_pp=="F_s") return -1.79087;
    if(species_pp=="Fe") return -8.31138;
    if(species_pp=="Fe_pv") return -8.45502;
    if(species_pp=="Ga_h") return -2.90153;
    if(species_pp=="Ge") return -4.49261;
    if(species_pp=="Ge_d") return -4.62213;
    if(species_pp=="Ge_h") return -4.50372;
    if(species_pp=="H") return -3.38635;
    if(species_pp=="H_h") return -3.36923;
    if(species_pp=="He") return 0.236864;
    if(species_pp=="Hf") return -9.95711;
    if(species_pp=="Hf_pv") return -9.95294;
    if(species_pp=="Ho_3") return -4.56866;
    if(species_pp=="I") return -1.51735;
    if(species_pp=="In_d") return -2.72115;
    if(species_pp=="Ir") return -8.85711;
    if(species_pp=="K_pv") return -1.02692;
    if(species_pp=="K_sv") return -1.09647;
    if(species_pp=="La") return -4.91782;
    if(species_pp=="La_s") return -4.86223;
    if(species_pp=="Li") return -1.89739;
    if(species_pp=="Mg") return -1.54148;
    if(species_pp=="Mg_pv") return -1.59348;
    if(species_pp=="Mn") return -9.028;
    if(species_pp=="Mo") return -10.9465;
    if(species_pp=="Mo_pv") return -10.8439;
    if(species_pp=="N") return -8.3187;
    if(species_pp=="N_h") return -8.3242;
    if(species_pp=="Na") return -1.3064;
    if(species_pp=="Na_sv") return -1.3139;
    if(species_pp=="Nb_pv") return -10.083;
    if(species_pp=="Nb_sv") return -10.2253;
    if(species_pp=="Ne") return -0.032434;
    if(species_pp=="Ni") return -5.57108;
    if(species_pp=="Ni_pv") return -5.77783;
    if(species_pp=="O") return -4.93114;
    if(species_pp=="O_h") return -5.01797;
    if(species_pp=="O_s") return -4.70103;
    if(species_pp=="Os") return -11.244;
    if(species_pp=="Os_pv") return -11.2193;
    if(species_pp=="P") return -5.32404;
    if(species_pp=="Pb") return -3.57056;
    if(species_pp=="Pb_d") return -3.7047;
    if(species_pp=="Pd") return -5.17847;
    if(species_pp=="Pd_pv") return -5.38197;
    if(species_pp=="Pt") return -6.05454;
    if(species_pp=="Re") return -12.4113;
    if(species_pp=="Re_pv") return -12.4325;
    if(species_pp=="Rh") return -7.27044;
    if(species_pp=="Rh_pv") return -7.34058;
    if(species_pp=="Ru") return -9.20414;
    if(species_pp=="Ru_pv") return -9.27135;
    if(species_pp=="S") return -4.12636;
    if(species_pp=="Sc_sv") return -6.33212;
    if(species_pp=="Se") return -3.48266;
    if(species_pp=="Si") return -5.42373;
    if(species_pp=="Si_h") return -5.44184;
    if(species_pp=="Sm_3") return -4.71227;
    if(species_pp=="Sn") return -3.79514;
    if(species_pp=="Sn_d") return -3.96372;
    if(species_pp=="Sr_sv") return -1.68354;
    if(species_pp=="Ta") return -11.86;
    if(species_pp=="Ta_pv") return -11.8489;
    if(species_pp=="Tc") return -10.3047;
    if(species_pp=="Tc_pv") return -10.3597;
    if(species_pp=="Te") return -3.14141;
    if(species_pp=="Ti") return -7.76396;
    if(species_pp=="Ti_pv") return -7.89056;
    if(species_pp=="Ti_sv") return -7.93863;
    if(species_pp=="Tl") return -2.2435;
    if(species_pp=="Tl_d") return -2.36274;
    if(species_pp=="V") return -8.94402;
    if(species_pp=="V_pv") return -9.07822;
    if(species_pp=="V_sv") return -9.11496;
    if(species_pp=="W") return -13.0124;
    if(species_pp=="W_pv") return -12.9546;
    if(species_pp=="Y_sv") return -6.46317;
    if(species_pp=="Yb") return -1.67034;
    if(species_pp=="Yb_2") return -1.51978;
    if(species_pp=="Zn") return -1.26581;
    if(species_pp=="Zr") return -8.47756;
    if(species_pp=="Zr_sv") return -8.54365;
    return NNN;
  }

  void insertElementalPropertiesCCE(const vector<string>& vproperties,const xelement::xelement& xel,vector<string>& vitems) {
    uint i=0,j=0;
    int index=0;
    for(i=0;i<vproperties.size();i++){
      if(vproperties[i]=="lattice_constants"||vproperties[i]=="lattice_angles"){
        const xvector<double>& xvec=xel.getPropertyXVectorDouble(vproperties[i]);
        for(index=xvec.lrows;index<=xvec.urows;index++){
          vitems.push_back(aurostd::utype2string(xvec[index],_DOUBLE_WRITE_PRECISION_));
        }
      }
      else if(vproperties[i]=="energies_ionization"){
        //vitems.push_back(xel.getPropertyString(vproperties[i],",",2)); //only go to second ionization
        const vector<double> vec=xel.energies_ionization;
        for(j=0;j<_ENERGIES_IONIZATION_MAX_;j++){
          if(j>vec.size()-1){vitems.push_back(aurostd::utype2string(NNN,_DOUBLE_WRITE_PRECISION_));continue;}
          vitems.push_back(aurostd::utype2string(vec[j],_DOUBLE_WRITE_PRECISION_));
        }
      }
      else{vitems.push_back(xel.getPropertyString(vproperties[i],","));}
    }
  }

  void insertCrystalPropertiesCCE(const string& structure_path,const string& anion,const vector<string>& vheaders,vector<string>& vitems) {
    string soliloquy=XPID+"aflowML::insertCrystalPropertiesCCE():";
    uint i=0,index_cation=0,index_anion=0;
    aflowlib::_aflowlib_entry entry;
    string entry_path="/common/LIB2/RAW/"+structure_path+"/aflowlib.out";
    if(!aurostd::FileExist(entry_path)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,entry_path+" not found",_FILE_NOT_FOUND_);}
    entry.file2aflowlib(entry_path);
    //
    //stoich features
    vector<string> vheaders_stoich;
    vector<double> vfeatures_stoich;
    entry.getStoichFeatures(vheaders_stoich,vfeatures_stoich);
    for(i=0;i<vfeatures_stoich.size();i++){vitems.push_back(aurostd::utype2string(vfeatures_stoich[i],_DOUBLE_WRITE_PRECISION_));}
    //
    for(i=0;i<entry.vgeometry.size();i++){vitems.push_back(aurostd::utype2string(entry.vgeometry[i],_DOUBLE_WRITE_PRECISION_));}
    vitems.push_back(aurostd::utype2string(entry.vgeometry[0]/entry.vgeometry[1],_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.vgeometry[1]/entry.vgeometry[2],_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.vgeometry[0]/entry.vgeometry[2],_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.natoms));
    index_cation=index_anion=AUROSTD_MAX_UINT;
    if(entry.vspecies.size()!=2){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"non-binary found",_FILE_CORRUPT_);}
    if(!(entry.vspecies[0]==anion||entry.vspecies[1]==anion)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"no N found",_FILE_CORRUPT_);}
    if(entry.vspecies[0]==anion){index_cation=1;index_anion=0;}else{index_cation=0;index_anion=1;}
    vitems.push_back(aurostd::utype2string(entry.vcomposition[index_cation]));
    vitems.push_back(aurostd::utype2string(entry.vcomposition[index_anion]));
    vitems.push_back(aurostd::utype2string(entry.vstoichiometry[index_cation]));
    vitems.push_back(aurostd::utype2string(entry.vstoichiometry[index_anion]));
    vitems.push_back(aurostd::utype2string(entry.density,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.vspecies_pp_ZVAL[index_cation]));
    vitems.push_back(aurostd::utype2string(entry.vspecies_pp_ZVAL[index_anion]));
    vitems.push_back(aurostd::utype2string(entry.valence_cell_iupac));
    vitems.push_back(aurostd::utype2string(entry.valence_cell_std));
    vitems.push_back(aurostd::utype2string(entry.volume_cell,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.volume_atom,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.energy_cell,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.energy_atom,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.energy_cutoff,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.enthalpy_formation_cell,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.enthalpy_formation_atom,_DOUBLE_WRITE_PRECISION_));
    if(entry.vnbondxx.size()!=3){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"non-binary found (nbondxx)",_FILE_CORRUPT_);}
    vitems.push_back(aurostd::utype2string(entry.vnbondxx[ (index_cation==0?0:2) ],_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.vnbondxx[ 1 ],_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.vnbondxx[ (index_cation==0?2:0) ],_DOUBLE_WRITE_PRECISION_));
    //for(i=0;i<entry.vnbondxx.size();i++){vitems.push_back(aurostd::utype2string(entry.vnbondxx[i],_DOUBLE_WRITE_PRECISION_));}
    vitems.push_back(aurostd::utype2string(entry.spacegroup_relax));
    vitems.push_back(aurostd::utype2string(entry.point_group_order));
    vitems.push_back(entry.Bravais_lattice_relax);
    //
  }
  void writeCCECSV() {
    bool LDEBUG=(TRUE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"aflowML::writeCCECSV():";

    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    uint ielement=0,ioxidation=0,i=0,j=0,k=0,ipp=0;
    string correction_line="",bader_line="",input_pre="",input="";
    xelement::xelement xel_cation,xel_N,xel_O;xel_N.populate("N");xel_O.populate("O");
    //cce::get_corrections_line_O("Al_+3_N");
    vector<vector<string> > vlines;
    vector<string> vitems_pre,vitems,vtokens;
    vlines.resize(2); //N,O
    
    vector<string> vheaders;
    vector<string> vions;aurostd::string2tokens("cation,anion",vions,",");
    vector<string> vproperties_elements,vproperties_elements_full;

    //aurostd::string2tokens("symbol,Z,period,group,mass,volume_molar,volume,valence_std,valence_iupac,valence_PT,valence_s,valence_p,valence_d,valence_f,density_PT,spacegroup_number,variance_parameter_mass,lattice_constants,lattice_angles,radius_Saxena,radius_PT,radius_covalent_PT,radius_covalent,radius_VanDerWaals_PT,radii_Ghosh08,radii_Slatter,radii_Pyykko,conductivity_electrical,electronegativity_Pauling,hardness_chemical_Ghosh,electronegativity_Pearson,electronegativity_Ghosh,electronegativity_Allen,electron_affinity_PT,energies_ionization,scale_Pettifor,temperature_boiling,temperature_melting,enthalpy_fusion,enthalpy_vaporization,enthalpy_atomization_WE,energy_cohesive,specific_heat_PT,critical_pressure,critical_temperature_PT,thermal_expansion,conductivity_thermal,hardness_mechanical_Brinell,hardness_mechanical_Mohs,hardness_mechanical_Vickers,hardness_chemical_Pearson,hardness_chemical_Putz,hardness_chemical_RB,modulus_shear,modulus_Young,modulus_bulk,Poisson_ratio_PT,refractive_index",vproperties_elements,",");
    //string properties first
    aurostd::string2tokens("symbol,period,group",vproperties_elements,",");
    //number properties next
    aurostd::string2tokens(_AFLOW_XELEMENT_PROPERTIES_ALL_,vproperties_elements_full,",");
    for(j=0;j<vproperties_elements_full.size();j++){
      if(xel_N.getType(vproperties_elements_full[j])=="number"||xel_N.getType(vproperties_elements_full[j])=="numbers"){
        if(vproperties_elements_full[j]!="oxidation_states" && vproperties_elements_full[j]!="oxidation_states_preferred"){ //exclude these entirely
          vproperties_elements.push_back(vproperties_elements_full[j]);
        }
      }
    }
    for(i=0;i<vions.size();i++){
      for(j=0;j<vproperties_elements.size();j++){
        if(vproperties_elements[j]=="lattice_constants"){
          vheaders.push_back(vproperties_elements[j]+"_a_"+vions[i]);
          vheaders.push_back(vproperties_elements[j]+"_b_"+vions[i]);
          vheaders.push_back(vproperties_elements[j]+"_c_"+vions[i]);
        }
        else if(vproperties_elements[j]=="lattice_angles"){
          vheaders.push_back(vproperties_elements[j]+"_alpha_"+vions[i]);
          vheaders.push_back(vproperties_elements[j]+"_beta_"+vions[i]);
          vheaders.push_back(vproperties_elements[j]+"_gamma_"+vions[i]);
        }
        else if(vproperties_elements[j]=="energies_ionization"){
          for(k=0;k<_ENERGIES_IONIZATION_MAX_;k++){vheaders.push_back(vproperties_elements[j]+"_"+aurostd::utype2string(k+1)+"_"+vions[i]);}
        }
        else{vheaders.push_back(vproperties_elements[j]+"_"+vions[i]);}
      }
      vheaders.push_back("oxidation_"+vions[i]);
      vheaders.push_back("energy_groundstate_PBE_"+vions[i]);
      vheaders.push_back("EATOM_PBE_"+vions[i]);
    }
    //
    vheaders.push_back("charge_bader_cation_PBE");
    //
    vheaders.push_back("PBE_298.15K");
    vheaders.push_back("PBE_0K");
    vheaders.push_back("LDA_298.15K");
    vheaders.push_back("LDA_0K");
    vheaders.push_back("SCAN_298.15K");
    vheaders.push_back("SCAN_0K");
    vheaders.push_back("PBE+U_298.15K");
    vheaders.push_back("PBE+U_0K");
    vheaders.push_back("exp_298.15K");
    vheaders.push_back("M-X_bonds");
    vheaders.push_back("natoms_per_f.u._cation");
    vheaders.push_back("natoms_per_f.u._anion");
    vheaders.push_back("enthalpy_formation_atom_exp");
    //
    vector<string> vheaders_stoich;
    aflowlib::_aflowlib_entry entry;
    entry.getStoichFeatures(vheaders_stoich);
    vheaders.insert(vheaders.end(),vheaders_stoich.begin(),vheaders_stoich.end());
    //
    vheaders.push_back("geometry_a_crystal");
    vheaders.push_back("geometry_b_crystal");
    vheaders.push_back("geometry_c_crystal");
    vheaders.push_back("geometry_alpha_crystal");
    vheaders.push_back("geometry_beta_crystal");
    vheaders.push_back("geometry_gamma_crystal");
    vheaders.push_back("geometry_a_/_b_crystal");
    vheaders.push_back("geometry_b_/_c_crystal");
    vheaders.push_back("geometry_a_/_c_crystal");
    vheaders.push_back("natoms_crystal");
    vheaders.push_back("natoms_cation_crystal");
    vheaders.push_back("natoms_anion_crystal");
    vheaders.push_back("stoich_cation_crystal");
    vheaders.push_back("stoich_anion_crystal");
    vheaders.push_back("density_crystal");
    vheaders.push_back("ZVAL_cation_crystal");  //not Z
    vheaders.push_back("ZVAL_anion_crystal"); //not Z
    vheaders.push_back("ZVAL_iupac_crystal");  //IUPAC
    vheaders.push_back("ZVAL_std_crystal");
    vheaders.push_back("volume_cell_crystal");
    vheaders.push_back("volume_atom_crystal");
    vheaders.push_back("energy_cell_crystal");
    vheaders.push_back("energy_atom_crystal");
    vheaders.push_back("energy_cutoff_crystal");
    vheaders.push_back("enthalpy_formation_cell_crystal");
    vheaders.push_back("enthalpy_formation_atom_crystal");
    vheaders.push_back("distance_cation_cation_crystal");
    vheaders.push_back("distance_cation_anion_crystal");
    vheaders.push_back("distance_anion_anion_crystal");
    vheaders.push_back("spacegroup_crystal");
    vheaders.push_back("point_group_order_crystal");
    vheaders.push_back("Bravais_lattice_crystal");

    for(i=0;i<vlines.size();i++){vlines[i].push_back(aurostd::joinWDelimiter(vheaders,","));}
    
    string species_pp="";
    bool found_pp=false;
    string structure_path="";
    for(ielement=0;ielement<100;ielement++){
      for(ioxidation=0;ioxidation<10;ioxidation++){
        xel_cation.populate(ielement);
        input_pre=xel_cation.symbol;
        input_pre+="_+"+aurostd::utype2string(ioxidation);
        if(LDEBUG){cerr << soliloquy << " input=" << input_pre << endl;}
        vitems_pre.clear();
        //cation
        insertElementalPropertiesCCE(vproperties_elements,xel_cation,vitems_pre);
        //
        vitems_pre.push_back(aurostd::utype2string(ioxidation,_DOUBLE_WRITE_PRECISION_)); //ioxidation
        try{species_pp=AVASP_Get_PseudoPotential_PAW_PBE(xel_cation.symbol);}
        catch(aurostd::xerror& excpt){continue;}
        found_pp=false;
        for(ipp=0;ipp<vxpseudopotential.size()&&found_pp==false;ipp++) {
          if(vxpseudopotential[ipp].species_pp_type[0]=="PAW_PBE" && vxpseudopotential[ipp].species_pp[0]==species_pp){
            found_pp=true;
            //vitems_pre.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_Z[0],_DOUBLE_WRITE_PRECISION_)); //Z
            vitems_pre.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_pp_groundstate_energy[0],_DOUBLE_WRITE_PRECISION_)); //groundstate_energy
            vitems_pre.push_back(aurostd::utype2string(vxpseudopotential[ipp].vEATOM[0],_DOUBLE_WRITE_PRECISION_)); //EATOM
          }
        }

        //N
        input=input_pre+"_N";
        correction_line=cce::get_corrections_line_N(input);
        if(!correction_line.empty()){
          vitems.clear();for(i=0;i<vitems_pre.size();i++){vitems.push_back(vitems_pre[i]);}
          //anion
          insertElementalPropertiesCCE(vproperties_elements,xel_N,vitems);
          //
          vitems.push_back(aurostd::utype2string(-3,_DOUBLE_WRITE_PRECISION_)); //ioxidation
          try{species_pp=AVASP_Get_PseudoPotential_PAW_PBE(xel_N.symbol);}
          catch(aurostd::xerror& excpt){continue;}
          found_pp=false;
          for(ipp=0;ipp<vxpseudopotential.size()&&found_pp==false;ipp++) {
            if(vxpseudopotential[ipp].species_pp_type[0]=="PAW_PBE" && vxpseudopotential[ipp].species_pp[0]==species_pp){
              found_pp=true;
              if(vxpseudopotential[ipp].species_pp_groundstate_structure[0]!="diatom"){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"pp_groundstate_structure for "+xel_N.symbol+" is NOT diatom",_INPUT_ILLEGAL_);}
              //vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_Z[0],_DOUBLE_WRITE_PRECISION_)); //Z
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_pp_groundstate_energy[0],_DOUBLE_WRITE_PRECISION_)); //groundstate_energy
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].vEATOM[0],_DOUBLE_WRITE_PRECISION_)); //EATOM
            }
          }
          //
          vitems.push_back("0.0");  //no bader yet

          //
          if(LDEBUG){cerr << soliloquy << " correction_line=\"" << correction_line << "\"" << endl;}
          aurostd::string2tokens(correction_line,vtokens," ");
          if(vtokens.size()<16){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vtokens.size()<16",_FILE_CORRUPT_);}
          for(i=2;i<=11;i++){vitems.push_back(vtokens[i]);}
          //
          structure_path=vtokens[12];
          if(LDEBUG){cerr << soliloquy << " structure_path=" << structure_path << endl;}
          //
          vitems.push_back(vtokens[13]);  //ncations_per_f.u.
          vitems.push_back(vtokens[14]);  //nanions_per_f.u.
          vitems.push_back(vtokens[15]);  //H_f^exp_298.15K
          //
          insertCrystalPropertiesCCE(structure_path,"N",vheaders,vitems);
          //
          if(vitems.size()!=vheaders.size()){
            for(uint ii=0;ii<vheaders.size()&&ii<vitems.size();ii++){
              cerr << vheaders[ii] << "=" << vitems[ii] << endl;
            }
            throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vitems.size()["+aurostd::utype2string(vitems.size())+"]!=vheaders.size()["+aurostd::utype2string(vheaders.size())+"]",_FILE_CORRUPT_);
          }
          vlines[0].push_back(aurostd::joinWDelimiter(vitems,","));
        }

        //O
        input=input_pre+"_O";
        correction_line=cce::get_corrections_line_O(input);
        if(!correction_line.empty()){
          vitems.clear();for(i=0;i<vitems_pre.size();i++){vitems.push_back(vitems_pre[i]);}
          //anion
          insertElementalPropertiesCCE(vproperties_elements,xel_O,vitems);
          //
          vitems.push_back(aurostd::utype2string(-2,_DOUBLE_WRITE_PRECISION_)); //ioxidation
          try{species_pp=AVASP_Get_PseudoPotential_PAW_PBE(xel_O.symbol);}
          catch(aurostd::xerror& excpt){continue;}
          found_pp=false;
          for(ipp=0;ipp<vxpseudopotential.size()&&found_pp==false;ipp++) {
            if(vxpseudopotential[ipp].species_pp_type[0]=="PAW_PBE" && vxpseudopotential[ipp].species_pp[0]==species_pp){
              found_pp=true;
              if(vxpseudopotential[ipp].species_pp_groundstate_structure[0]!="diatom"){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"pp_groundstate_structure for "+xel_O.symbol+" is NOT diatom",_INPUT_ILLEGAL_);}
              //vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_Z[0],_DOUBLE_WRITE_PRECISION_)); //Z
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_pp_groundstate_energy[0],_DOUBLE_WRITE_PRECISION_)); //groundstate_energy
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].vEATOM[0],_DOUBLE_WRITE_PRECISION_)); //EATOM
            }
          }
          //
          bader_line=cce::get_Bader_templates(xel_cation.symbol);
          if(!bader_line.empty()){
            if(LDEBUG){cerr << soliloquy << " bader_line=\"" << bader_line << "\"" << endl;}
            aurostd::string2tokens(bader_line,vtokens," ");
            if(vtokens.size()<6){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vtokens.size()<13",_FILE_CORRUPT_);}
            vitems.push_back(vtokens[2]);
          }

          //
          if(LDEBUG){cerr << soliloquy << " correction_line=\"" << correction_line << "\"" << endl;}
          aurostd::string2tokens(correction_line,vtokens," ");
          if(vtokens.size()<16){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vtokens.size()<16",_FILE_CORRUPT_);}
          for(i=2;i<=11;i++){vitems.push_back(vtokens[i]);}
          //
          structure_path=vtokens[12];
          if(LDEBUG){cerr << soliloquy << " structure_path=" << structure_path << endl;}
          //
          vitems.push_back(vtokens[13]);  //ncations_per_f.u.
          vitems.push_back(vtokens[14]);  //nanions_per_f.u.
          vitems.push_back(vtokens[15]);  //H_f^exp_298.15K
          //
          insertCrystalPropertiesCCE(structure_path,"O",vheaders,vitems);
          //
          if(vitems.size()!=vheaders.size()){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vitems.size()!=vheaders.size()",_FILE_CORRUPT_);}
          vlines[1].push_back(aurostd::joinWDelimiter(vitems,","));
        }
      }
    }

    stringstream file;
    //N
    aurostd::StringstreamClean(file);
    for(i=0;i<vlines[0].size();i++){file << vlines[0][i] << endl;}
    aurostd::stringstream2file(file,"cce_data_N.csv");
    //O
    aurostd::StringstreamClean(file);
    for(i=0;i<vlines[1].size();i++){file << vlines[1][i] << endl;}
    aurostd::stringstream2file(file,"cce_data_O.csv");
    //total
    aurostd::StringstreamClean(file);
    for(i=0;i<vlines[0].size();i++){file << vlines[0][i] << endl;}
    for(i=1;i<vlines[1].size();i++){file << vlines[1][i] << endl;} //skip header
    aurostd::stringstream2file(file,"cce_data_NO.csv");

  }
} // namespace aflowML

namespace aflowlib {
  void insertStoichStats(const vector<string> vstats,const xvector<double>& nspecies_xv,const xvector<double>& stoich_xv,vector<double>& vfeatures){
    bool LDEBUG=(FALSE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"aflowlib::insertStoichStats():";

    uint k=0,l=0;
    int index=0;
    double d_tmp;
    vector<uint> vi_tmp;

    //check for NNN or AUROSTD_NAN
    bool has_NaN=false;
    for(index=nspecies_xv.lrows;index<=nspecies_xv.urows&&!has_NaN;index++){
      if(nspecies_xv[index]==NNN||nspecies_xv[index]==AUROSTD_NAN||nspecies_xv[index]==AUROSTD_MAX_DOUBLE){has_NaN=true;}
    }

    if(LDEBUG){cerr << soliloquy << " nspecies_xv=" << nspecies_xv << endl;}
    for(k=0;k<vstats.size();k++){
      if(vstats[k]=="min"){if(has_NaN){vfeatures.push_back(NNN);continue;} vfeatures.push_back( aurostd::min(nspecies_xv) );}
      else if(vstats[k]=="max"){if(has_NaN){vfeatures.push_back(NNN);continue;} vfeatures.push_back( aurostd::max(nspecies_xv) );}
      else if(vstats[k]=="range"){if(has_NaN){vfeatures.push_back(NNN);continue;} vfeatures.push_back( aurostd::max(nspecies_xv) - aurostd::min(nspecies_xv) );}
      else if(vstats[k]=="mean"){if(has_NaN){vfeatures.push_back(NNN);continue;} vfeatures.push_back( aurostd::scalar_product(stoich_xv,nspecies_xv) );}
      else if(vstats[k]=="dev"){
        if(has_NaN){vfeatures.push_back(NNN);continue;}
        d_tmp=aurostd::scalar_product(stoich_xv,nspecies_xv);  //mean
        vfeatures.push_back( aurostd::scalar_product(stoich_xv,aurostd::abs(nspecies_xv-d_tmp)) );
      }
      else if(vstats[k]=="mode"){ //property of most promiment species
        if(has_NaN){vfeatures.push_back(NNN);continue;}
        d_tmp=aurostd::max(stoich_xv);  //stoich_max
        vi_tmp.clear();
        for(index=stoich_xv.lrows;index<=stoich_xv.urows;index++){
          if(aurostd::isequal(stoich_xv[index],d_tmp)){vi_tmp.push_back(index);}
        }
        if(vi_tmp.size()==1){vfeatures.push_back( nspecies_xv[vi_tmp[0]] );}  //easy case
        else{
          //take average
          d_tmp=0;
          for(l=0;l<vi_tmp.size();l++){d_tmp+=nspecies_xv[vi_tmp[l]];}
          d_tmp/=(double)vi_tmp.size();
          vfeatures.push_back( d_tmp );
        }
      }
      else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown statistic type: "+vstats[k],_RUNTIME_ERROR_);}
    }
  }

  void _aflowlib_entry::getStoichFeatures(vector<string>& vheaders){
    vector<double> vfeatures; //dummy
    return getStoichFeatures(vheaders,vfeatures,true);
  }
  void _aflowlib_entry::getStoichFeatures(vector<string>& vheaders,vector<double>& vfeatures,bool vheaders_only){
    //follows supplementary of 10.1038/npjcompumats.2016.28
    bool LDEBUG=(FALSE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"_aflowlib_entry::getStoichFeatures():";
    vheaders.clear();vfeatures.clear();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //headers
    stringstream tmp_ss;
    xelement::xelement xel;
    vector<xelement::xelement> vxel;
    uint i=0,j=0,k=0;

    //L^p norms
    vector<uint> vp;
    aurostd::string2tokens("0,2,3,4,5,6,7,8,9,10",vp,",");
    for(i=0;i<vp.size();i++){
      aurostd::StringstreamClean(tmp_ss);
      tmp_ss << "stoich_norm_p_" << std::setfill('0') << std::setw(2) << vp[i];
      vheaders.push_back(tmp_ss.str());
    }

    //element-property-based
    //get which properties to average
    xel.populate(1);  //dummy to get properties
    vector<string> vproperties_full,vproperties;
    vector<string> vstats;
    aurostd::string2tokens(_AFLOW_XELEMENT_PROPERTIES_ALL_,vproperties_full,",");
    aurostd::string2tokens("min,max,range,mean,dev,mode",vstats,",");
    //load up vheaders
    for(i=0;i<vproperties_full.size();i++){
      if(vproperties_full[i]=="oxidation_states"){continue;} //skip this
      if(xel.getType(vproperties_full[i])=="number"||xel.getType(vproperties_full[i])=="numbers"){
        vproperties.push_back(vproperties_full[i]);
        
        if(xel.getType(vproperties.back())=="number"){
          if(LDEBUG){cerr << soliloquy << " " << vproperties.back() << " is a number" << endl;}
          for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_stoich_"+vstats[j]);}
        }
        else if(xel.getType(vproperties.back())=="numbers"){
          if(LDEBUG){cerr << soliloquy << " " << vproperties.back() << " are numbers" << endl;}
          if(vproperties.back()=="lattice_constants"){
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_a_stoich_"+vstats[j]);}
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_b_stoich_"+vstats[j]);}
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_c_stoich_"+vstats[j]);}
          }
          else if(vproperties.back()=="lattice_angles"){
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_alpha_stoich_"+vstats[j]);}
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_beta_stoich_"+vstats[j]);}
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_gamma_stoich_"+vstats[j]);}
          }
          else if(vproperties.back()=="oxidation_states_preferred"){
            for(j=0;j<vstats.size();j++){
              vheaders.push_back(vproperties.back()+"_stoich_"+vstats[j]);  //only use 0th oxidation_state_preferred
            }
          }
          else if(vproperties.back()=="energies_ionization"){
            for(k=0;k<_ENERGIES_IONIZATION_MAX_;k++){
              for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_"+aurostd::utype2string(k+1)+"_stoich_"+vstats[j]);}
            }
            //for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_2_stoich_"+vstats[j]);}
            //for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_3_stoich_"+vstats[j]);}
            //for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_4_stoich_"+vstats[j]);}
            //for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_5_stoich_"+vstats[j]);}
          }
          else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown numbers type: "+vproperties.back(),_RUNTIME_ERROR_);}
        }
      }
    }

    //valence (un)occupation
    vector<string> vorbitals;
    aurostd::string2tokens("s,p,d,f",vorbitals,",");
    for(i=0;i<vorbitals.size();i++){
      vheaders.push_back("valence_fraction_occupied_"+vorbitals[i]);
      vheaders.push_back("valence_fraction_unoccupied_"+vorbitals[i]);
    }

    //ionic character
    vheaders.push_back("formability_ionic");
    vector<string> vEN;
    for(i=0;i<vproperties.size();i++){
      if(vproperties[i].find("electronegativity")!=string::npos && xel.getUnits(vproperties[i]).empty()){ //must have NO units (goes in exp)
        vEN.push_back(vproperties[i]);
      }
    }
    for(i=0;i<vEN.size();i++){
      vheaders.push_back("character_ionic_"+vEN[i]+"_max");
      vheaders.push_back("character_ionic_"+vEN[i]+"_mean");
    }

    if(LDEBUG){
      for(i=0;i<vheaders.size();i++){cerr << soliloquy << " vheaders[i=" << i << "]=\"" << vheaders[i] << "\"" << endl;}
    }

    if(vheaders_only) return;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //features
    
    uint nspecies=vspecies.size();
    if(nspecies==0){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"nspecies==0",_RUNTIME_ERROR_);}
    xvector<double> nspecies_xv(nspecies);
    
    //L^p norms
    xvector<double> stoich_xv(nspecies);
    for(j=0;j<nspecies;j++){stoich_xv[stoich_xv.lrows+j]=vcomposition[j]/natoms;}
    if(LDEBUG){cerr << soliloquy << " stoich_xv=" << stoich_xv << endl;}
    for(i=0;i<vp.size();i++){
      for(j=0;j<nspecies;j++){
        nspecies_xv[nspecies_xv.lrows+j]=std::pow(stoich_xv[stoich_xv.lrows+j],(double)vp[i]);
      }
      if(LDEBUG){cerr << soliloquy << " nspecies_xv[\"stoich_norm_p_"+aurostd::utype2string(vp[i])+"\"]=" << nspecies_xv << endl;}
      vfeatures.push_back( std::pow(sum(nspecies_xv),(vp[i]==0?1.0:1.0/vp[i])) );
    }
    
    //element-property-based
    //load up vxel
    int index=0,index_min=0,index_max=0;
    for(j=0;j<nspecies;j++){
      vxel.push_back(xelement::xelement(vspecies[j]));
    }
    for(i=0;i<vproperties.size();i++){
      if(xel.getType(vproperties[i])=="number"){
        for(j=0;j<nspecies;j++){
          nspecies_xv[nspecies_xv.lrows+j]=vxel[j].getPropertyDouble(vproperties[i]);
        }
        if(LDEBUG){cerr << soliloquy << " nspecies_xv[\""+vproperties[i]+"\"]=" << nspecies_xv << endl;}
        insertStoichStats(vstats,nspecies_xv,stoich_xv,vfeatures);
      }
      if(xel.getType(vproperties[i])=="numbers"){
        if(vproperties[i]=="lattice_constants"||vproperties[i]=="lattice_angles"){
          index_min=1;index_max=3;
          for(index=index_min;index<=index_max;index++){
            for(j=0;j<nspecies;j++){
              const xvector<double>& xvec=vxel[j].getPropertyXVectorDouble(vproperties[i]);
              nspecies_xv[nspecies_xv.lrows+j]=xvec[index];
            }
            if(LDEBUG){cerr << soliloquy << " nspecies_xv[\""+vproperties[i]+"_index_"+aurostd::utype2string(index)+"\"]=" << nspecies_xv << endl;}
            insertStoichStats(vstats,nspecies_xv,stoich_xv,vfeatures);
          }
        }
        else if(vproperties[i]=="oxidation_states_preferred"||vproperties[i]=="energies_ionization"){
          if(vproperties[i]=="oxidation_states_preferred"){index_min=0;index_max=0;}
          else if(vproperties[i]=="energies_ionization"){index_min=0;index_max=_ENERGIES_IONIZATION_MAX_-1;}
          else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown numbers property (vector): "+vproperties[i],_RUNTIME_ERROR_);}
          for(index=index_min;index<=index_max;index++){
            for(j=0;j<nspecies;j++){
              const vector<double>& vec=vxel[j].getPropertyVectorDouble(vproperties[i]);
              nspecies_xv[nspecies_xv.lrows+j]=(index<(int)vec.size()?vec[index]:NNN);
            }
            if(LDEBUG){cerr << soliloquy << " nspecies_xv[\""+vproperties[i]+"_index="+aurostd::utype2string(index)+"\"]=" << nspecies_xv << endl;}
            insertStoichStats(vstats,nspecies_xv,stoich_xv,vfeatures);
          }
        }
        else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown numbers property: "+vproperties[i],_RUNTIME_ERROR_);}
      }
    }
    
    //valence (un)occupation
    for(j=0;j<nspecies;j++){nspecies_xv[nspecies_xv.lrows+j]=vxel[j].getPropertyDouble("valence_std");}
    double denom=aurostd::scalar_product(stoich_xv,nspecies_xv); //same for all quantities
    vector<double> vval_total;
    aurostd::string2tokens("2,6,10,14",vval_total,",");
    for(i=0;i<vorbitals.size();i++){
      for(j=0;j<nspecies;j++){nspecies_xv[nspecies_xv.lrows+j]=vxel[j].getPropertyDouble("valence_"+vorbitals[i]);}  //populate with orbital occupation
      vfeatures.push_back( aurostd::scalar_product(stoich_xv,nspecies_xv)/denom );  //occupied
      for(j=0;j<nspecies;j++){nspecies_xv[nspecies_xv.lrows+j]=vval_total[i]-nspecies_xv[nspecies_xv.lrows+j];} //has occupied inside already
      vfeatures.push_back( aurostd::scalar_product(stoich_xv,nspecies_xv)/denom );  //unoccupied
    }
    
    //ionic character
    //ionic formability
    bool formability_ionic=false;
    bool has_NaN=false;
    vector<int> nspecies_v;
    aurostd::xcombos xc;
    for(j=0;j<nspecies&&!has_NaN;j++){
      const vector<double>& oxidation_states=vxel[j].getPropertyVectorDouble("oxidation_states");
      if(oxidation_states.size()==0){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"No oxidation states found for element: "+vxel[j].symbol,_RUNTIME_ERROR_);}  //should have NNN
      if(oxidation_states.size()==1 && oxidation_states[0]==NNN){has_NaN=true;}
      nspecies_v.push_back((int)oxidation_states.size());
    }
    if(has_NaN){cerr << soliloquy << " has NaN" << endl;}
    if(!has_NaN){
      xvector<double> natoms_xv(natoms);
      if(LDEBUG){cerr << soliloquy << " oxidation_states_count=" << aurostd::joinWDelimiter(nspecies_v,",") << endl;}
      xc.reset(nspecies_v,'E');
      while(xc.increment()&&!formability_ionic){
        const vector<int>& indices=xc.getCombo();
        if(LDEBUG){cerr << soliloquy << " indices=" << aurostd::joinWDelimiter(indices,",") << endl;}
        i=0;
        for(j=0;j<nspecies&&!has_NaN;j++){
          const vector<double>& oxidation_states=vxel[j].getPropertyVectorDouble("oxidation_states");
          for(k=0;k<vcomposition[j];k++){
            if(oxidation_states[indices[j]]==NNN){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Found NaN among populated oxidation_states",_RUNTIME_ERROR_);}
            natoms_xv[natoms_xv.lrows+(i++)]=oxidation_states[indices[j]];
          }
        }
        if(LDEBUG){cerr << soliloquy << " natoms_xv[\"oxidation_states\"]=" << natoms_xv << endl;}
        if(aurostd::isequal(aurostd::sum(natoms_xv),0.0)){formability_ionic=true;}
      }
    }
    vfeatures.push_back( (!has_NaN&&formability_ionic?1:0) );
    //character_ionic_max and _mean
    xvector<double> pairs_xv(aurostd::nCk((int)nspecies,2));  //electronegativities
    xvector<double> pairs2_xv(aurostd::nCk((int)nspecies,2)); //stoich
    for(i=0;i<vEN.size();i++){
      if(LDEBUG){cerr << soliloquy << " EN=" << vEN[i] << endl;}
      xc.reset(nspecies,2);
      k=0;
      while(xc.increment()){
        const vector<int>& indices=xc.getIndices();
        if(LDEBUG){cerr << soliloquy << " indices=" << aurostd::joinWDelimiter(indices,",") << endl;}
        pairs_xv[pairs_xv.lrows+k]=1.0-std::exp(-0.25*(std::pow(vxel[indices[0]].getPropertyDouble(vEN[i])-vxel[indices[1]].getPropertyDouble(vEN[i]),2.0)));
        pairs2_xv[pairs2_xv.lrows+k]=stoich_xv[stoich_xv.lrows+indices[0]]*stoich_xv[stoich_xv.lrows+indices[1]];
        k++;
      }
      if(LDEBUG){cerr << soliloquy << " pairs_xv[\""+vEN[i]+"\"]=" << pairs_xv << endl;}
      vfeatures.push_back( aurostd::max(pairs_xv) );
      vfeatures.push_back( aurostd::scalar_product(pairs_xv,pairs2_xv) );
    }

    if(vheaders.size()!=vfeatures.size()){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vheaders.size()!=vfeatures.size()",_RUNTIME_ERROR_);}

    if(LDEBUG){
      for(i=0;i<vheaders.size();i++){cerr << soliloquy << " " << vheaders[i] << "=" << vfeatures[i] << endl;}
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  }
} // namespace aflowlib

#endif  // _AFLOW_ML_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
