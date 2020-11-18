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

namespace aflowML {
  void insertCrystalPropertiesCCE(const string& structure_path,const string& anion,vector<string>& vitems) {
    string soliloquy=XPID+"aflowML::insertCrystalPropertiesCCE():";
    uint i=0,index_cation=0,index_anion=0;
    aflowlib::_aflowlib_entry entry;
    string entry_path="/common/LIB2/RAW/"+structure_path+"/aflowlib.out";
    if(!aurostd::FileExist(entry_path)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,entry_path+" not found",_FILE_NOT_FOUND_);}
    entry.file2aflowlib(entry_path);
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
  }
  void writeCCECSV() {
    bool LDEBUG=(TRUE || XHOST.DEBUG);
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
    vector<string> vproperties_elements;
    aurostd::string2tokens("symbol,Z,period,group,mass,volume_molar,volume,valence_std,valence_iupac,valence_PT,valence_s,valence_p,valence_d,valence_f,density_PT,spacegroup_number,variance_parameter_mass,lattice_constants,lattice_angles,radius_Saxena,radius_PT,radius_covalent_PT,radius_covalent,radius_VanDerWaals_PT,radii_Ghosh08,radii_Slatter,radii_Pyykko,conductivity_electrical,electronegativity_Pauling,hardness_chemical_Ghosh,electronegativity_Pearson,electronegativity_Ghosh,electronegativity_Allen,electron_affinity_PT,energies_ionization,scale_Pettifor,temperature_boiling,temperature_melting,enthalpy_fusion,enthalpy_vaporization,enthalpy_atomization_WE,energy_cohesive,specific_heat_PT,critical_pressure,critical_temperature_PT,thermal_expansion,conductivity_thermal,hardness_mechanical_Brinell,hardness_mechanical_Mohs,hardness_mechanical_Vickers,hardness_chemical_Pearson,hardness_chemical_Putz,hardness_chemical_RB,modulus_shear,modulus_Young,modulus_bulk,Poisson_ratio_PT,refractive_index",vproperties_elements,",");
    uint energy_ionization_max=5;
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
          for(k=0;k<energy_ionization_max;k++){vheaders.push_back(vproperties_elements[j]+"_"+aurostd::utype2string(k+1)+"_"+vions[i]);}
        }
        else{
          vheaders.push_back(vproperties_elements[j]+"_"+vions[i]);
        }
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
        for(i=0;i<vproperties_elements.size();i++){
          if(vproperties_elements[i]=="energies_ionization"){
            //vitems_pre.push_back(xel_cation.getProperty(vproperties_elements[i],",",2)); //only go to second ionization
            for(j=0;j<energy_ionization_max;j++){
              if(j>xel_cation.energies_ionization.size()-1){vitems_pre.push_back(aurostd::utype2string(NNN,_DOUBLE_WRITE_PRECISION_));continue;}
              vitems_pre.push_back(aurostd::utype2string(xel_cation.energies_ionization[j],_DOUBLE_WRITE_PRECISION_));
            }
          }
          else{vitems_pre.push_back(xel_cation.getProperty(vproperties_elements[i],","));}
        }
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
          for(i=0;i<vproperties_elements.size();i++){
            if(vproperties_elements[i]=="energies_ionization"){
              //vitems.push_back(xel_N.getProperty(vproperties_elements[i],",",2)); //only go to second ionization
              for(j=0;j<energy_ionization_max;j++){
                if(j>xel_N.energies_ionization.size()-1){vitems.push_back(aurostd::utype2string(NNN,_DOUBLE_WRITE_PRECISION_));continue;}
                vitems.push_back(aurostd::utype2string(xel_N.energies_ionization[j],_DOUBLE_WRITE_PRECISION_));
              }
            }
            else{vitems.push_back(xel_N.getProperty(vproperties_elements[i],","));}
          }
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
          insertCrystalPropertiesCCE(structure_path,"N",vitems);
          //
          vlines[0].push_back(aurostd::joinWDelimiter(vitems,","));
        }

        //O
        input=input_pre+"_O";
        correction_line=cce::get_corrections_line_O(input);
        if(!correction_line.empty()){
          vitems.clear();for(i=0;i<vitems_pre.size();i++){vitems.push_back(vitems_pre[i]);}
          //anion
          for(i=0;i<vproperties_elements.size();i++){
            if(vproperties_elements[i]=="energies_ionization"){
              //vitems.push_back(xel_O.getProperty(vproperties_elements[i],",",2)); //only go to second ionization
              for(j=0;j<energy_ionization_max;j++){
                if(j>xel_O.energies_ionization.size()-1){vitems.push_back(aurostd::utype2string(NNN,_DOUBLE_WRITE_PRECISION_));continue;}
                vitems.push_back(aurostd::utype2string(xel_O.energies_ionization[j],_DOUBLE_WRITE_PRECISION_));
              }
            }
            else{vitems.push_back(xel_O.getProperty(vproperties_elements[i],","));}
          }
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
          insertCrystalPropertiesCCE(structure_path,"O",vitems);
          //
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
}

#endif  // _AFLOW_ML_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
