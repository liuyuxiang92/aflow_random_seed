// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2018
// FILE "ANRL/aflow_anrl_ABC_hP36_175_jk_jk_jk.cpp"

#ifndef _AFLOW_ANRL_ABC_hP36_175_jk_jk_jk_CPP
#define _AFLOW_ANRL_ABC_hP36_175_jk_jk_jk_CPP
#include "../aflow.h"

namespace anrl {
  uint WebANRL_ABC_hP36_175_jk_jk_jk(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_ABC_hP36_175_jk_jk_jk(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system ABC_hP36_175_jk_jk_jk

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_ABC_hP36_175_jk_jk_jk(web,LDEBUG); // PLUG WEB STUFF
      #ifdef _ANRL_NOWEB_
      cout << "no web" << endl;
      #else
      cout << web.str() << endl;
      #endif
      exit(0);
    }

    vector<double> vparameters;
    aurostd::string2tokens(parameters,vparameters,",");

    uint nspecies,natoms,spacegroup,nunderscores,nparameters;
    string label,Pearson_symbol,params,Strukturbericht,prototype,dialect;

    anrl::vproto2tokens(proto_line,label,nspecies,natoms,spacegroup,nunderscores,nparameters,Pearson_symbol,params,Strukturbericht,prototype,dialect);

    if(!anrl::PrototypeANRL_Consistency(oss,vparameters.size(),nparameters,prototype,label,
                 Strukturbericht,Pearson_symbol,spacegroup, params, print_mode) && print_mode!=1) { exit(0);}    

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: vparameters.size()=" << vparameters.size() << endl;}

    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);

    _atom atom;

    if(print_mode==1 && vparameters.size()==0){
      for(uint n=0;n<nparameters;n++){
        vparameters.push_back(0);
      }
    }

    uint i=0;
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: a=" << a << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: c=" << c << " (c/a=" << covera << ")" << endl;}
    
    double x1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: x1=" << x1 << endl;}
    double y1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: y1=" << y1 << endl;}
    double x2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: x2=" << x2 << endl;}
    double y2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: y2=" << y2 << endl;}
    double x3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: x3=" << x3 << endl;}
    double y3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: y3=" << y3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: x4=" << x4 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: y4=" << y4 << endl;}
    double x5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: x5=" << x5 << endl;}
    double y5=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: y5=" << y5 << endl;}
    double x6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: x6=" << x6 << endl;}
    double y6=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_ABC_hP36_175_jk_jk_jk: y6=" << y6 << endl;}
        
    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(spacegroup)+DOI_ANRL; //CO190520
    str.scale=1.0;

    a1=(1.0/2.0)*a*xn-(sqrt(3.0)/2.0)*a*yn;
    a2=(1.0/2.0)*a*xn+(sqrt(3.0)/2.0)*a*yn;
    a3=c*zn;
    
    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    // symbolic representation of lattice vectors
    vector<string> a1_equation, a2_equation, a3_equation;
    a1_equation.push_back("(1.0/2.0)*a");a1_equation.push_back("-(sqrt(3.0)/2.0)*a");a1_equation.push_back("0");
    a2_equation.push_back("(1.0/2.0)*a");a2_equation.push_back("(sqrt(3.0)/2.0)*a");a2_equation.push_back("0");
    a3_equation.push_back("0");a3_equation.push_back("0");a3_equation.push_back("c");
    str.symbolic_math_lattice.push_back(a1_equation);
    str.symbolic_math_lattice.push_back(a2_equation);
    str.symbolic_math_lattice.push_back(a3_equation);
    
    str.num_lattice_parameters = 2;
    
    str.num_parameters = vparameters.size();
    vector<string> parameter_list; aurostd::string2tokens(params,parameter_list,",");
    str.prototype_parameter_list = parameter_list;
    str.prototype_parameter_values = vparameters;

    if(print_mode!=1){
      str.FixLattices(); // Reciprocal/f2c/c2f
    }
    
    atom.name="A"; atom.type=0;                                       // atom B1
    atom.fpos(1)=x1;atom.fpos(2)=y1;atom.fpos(3)=0.0;                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x1");atom.fpos_equation.push_back("y1");atom.fpos_equation.push_back("0.0");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1
    
    atom.name="A"; atom.type=0;                                       // atom B2
    atom.fpos(1)=-y1;atom.fpos(2)=(x1-y1);atom.fpos(3)=0.0;                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y1");atom.fpos_equation.push_back("(x1-y1)");atom.fpos_equation.push_back("0.0");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2
    
    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=(-x1+y1);atom.fpos(2)=-x1;atom.fpos(3)=0.0;                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x1+y1)");atom.fpos_equation.push_back("-x1");atom.fpos_equation.push_back("0.0");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3
    
    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=-x1;atom.fpos(2)=-y1;atom.fpos(3)=0.0;                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x1");atom.fpos_equation.push_back("-y1");atom.fpos_equation.push_back("0.0");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4
    
    atom.name="A"; atom.type=0;                                       // atom B5
    atom.fpos(1)=y1;atom.fpos(2)=(-x1+y1);atom.fpos(3)=0.0;                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y1");atom.fpos_equation.push_back("(-x1+y1)");atom.fpos_equation.push_back("0.0");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5
    
    atom.name="A"; atom.type=0;                                       // atom B6
    atom.fpos(1)=(x1-y1);atom.fpos(2)=x1;atom.fpos(3)=0.0;                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x1-y1)");atom.fpos_equation.push_back("x1");atom.fpos_equation.push_back("0.0");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6
    
    atom.name="A"; atom.type=0;                                       // atom B19
    atom.fpos(1)=x4;atom.fpos(2)=y4;atom.fpos(3)=(1.0/2.0);                     // atom B19
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("y4");atom.fpos_equation.push_back("(1.0/2.0)");// atom B19 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B19 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B19
    
    atom.name="A"; atom.type=0;                                       // atom B20
    atom.fpos(1)=-y4;atom.fpos(2)=(x4-y4);atom.fpos(3)=(1.0/2.0);                     // atom B20
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y4");atom.fpos_equation.push_back("(x4-y4)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B20 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B20 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B20
    
    atom.name="A"; atom.type=0;                                       // atom B21
    atom.fpos(1)=(-x4+y4);atom.fpos(2)=-x4;atom.fpos(3)=(1.0/2.0);                     // atom B21
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x4+y4)");atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("(1.0/2.0)");// atom B21 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B21 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B21
    
    atom.name="A"; atom.type=0;                                       // atom B22
    atom.fpos(1)=-x4;atom.fpos(2)=-y4;atom.fpos(3)=(1.0/2.0);                     // atom B22
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x4");atom.fpos_equation.push_back("-y4");atom.fpos_equation.push_back("(1.0/2.0)");// atom B22 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B22 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B22
    
    atom.name="A"; atom.type=0;                                       // atom B23
    atom.fpos(1)=y4;atom.fpos(2)=(-x4+y4);atom.fpos(3)=(1.0/2.0);                     // atom B23
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y4");atom.fpos_equation.push_back("(-x4+y4)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B23 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B23 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B23
    
    atom.name="A"; atom.type=0;                                       // atom B24
    atom.fpos(1)=(x4-y4);atom.fpos(2)=x4;atom.fpos(3)=(1.0/2.0);                     // atom B24
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x4-y4)");atom.fpos_equation.push_back("x4");atom.fpos_equation.push_back("(1.0/2.0)");// atom B24 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B24 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B24
    
    atom.name="B"; atom.type=1;                                       // atom B7
    atom.fpos(1)=x2;atom.fpos(2)=y2;atom.fpos(3)=0.0;                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x2");atom.fpos_equation.push_back("y2");atom.fpos_equation.push_back("0.0");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7
    
    atom.name="B"; atom.type=1;                                       // atom B8
    atom.fpos(1)=-y2;atom.fpos(2)=(x2-y2);atom.fpos(3)=0.0;                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y2");atom.fpos_equation.push_back("(x2-y2)");atom.fpos_equation.push_back("0.0");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8
    
    atom.name="B"; atom.type=1;                                       // atom B9
    atom.fpos(1)=(-x2+y2);atom.fpos(2)=-x2;atom.fpos(3)=0.0;                     // atom B9
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x2+y2)");atom.fpos_equation.push_back("-x2");atom.fpos_equation.push_back("0.0");// atom B9 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B9 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B9
    
    atom.name="B"; atom.type=1;                                       // atom B10
    atom.fpos(1)=-x2;atom.fpos(2)=-y2;atom.fpos(3)=0.0;                     // atom B10
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x2");atom.fpos_equation.push_back("-y2");atom.fpos_equation.push_back("0.0");// atom B10 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B10 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B10
    
    atom.name="B"; atom.type=1;                                       // atom B11
    atom.fpos(1)=y2;atom.fpos(2)=(-x2+y2);atom.fpos(3)=0.0;                     // atom B11
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y2");atom.fpos_equation.push_back("(-x2+y2)");atom.fpos_equation.push_back("0.0");// atom B11 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B11 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B11
    
    atom.name="B"; atom.type=1;                                       // atom B12
    atom.fpos(1)=(x2-y2);atom.fpos(2)=x2;atom.fpos(3)=0.0;                     // atom B12
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x2-y2)");atom.fpos_equation.push_back("x2");atom.fpos_equation.push_back("0.0");// atom B12 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B12 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B12
    
    atom.name="B"; atom.type=1;                                       // atom B25
    atom.fpos(1)=x5;atom.fpos(2)=y5;atom.fpos(3)=(1.0/2.0);                     // atom B25
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("(1.0/2.0)");// atom B25 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B25 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B25
    
    atom.name="B"; atom.type=1;                                       // atom B26
    atom.fpos(1)=-y5;atom.fpos(2)=(x5-y5);atom.fpos(3)=(1.0/2.0);                     // atom B26
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("(x5-y5)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B26 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B26 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B26
    
    atom.name="B"; atom.type=1;                                       // atom B27
    atom.fpos(1)=(-x5+y5);atom.fpos(2)=-x5;atom.fpos(3)=(1.0/2.0);                     // atom B27
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x5+y5)");atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("(1.0/2.0)");// atom B27 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B27 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B27
    
    atom.name="B"; atom.type=1;                                       // atom B28
    atom.fpos(1)=-x5;atom.fpos(2)=-y5;atom.fpos(3)=(1.0/2.0);                     // atom B28
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x5");atom.fpos_equation.push_back("-y5");atom.fpos_equation.push_back("(1.0/2.0)");// atom B28 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B28 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B28
    
    atom.name="B"; atom.type=1;                                       // atom B29
    atom.fpos(1)=y5;atom.fpos(2)=(-x5+y5);atom.fpos(3)=(1.0/2.0);                     // atom B29
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y5");atom.fpos_equation.push_back("(-x5+y5)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B29 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B29 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B29
    
    atom.name="B"; atom.type=1;                                       // atom B30
    atom.fpos(1)=(x5-y5);atom.fpos(2)=x5;atom.fpos(3)=(1.0/2.0);                     // atom B30
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x5-y5)");atom.fpos_equation.push_back("x5");atom.fpos_equation.push_back("(1.0/2.0)");// atom B30 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B30 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B30
    
    atom.name="C"; atom.type=2;                                       // atom B13
    atom.fpos(1)=x3;atom.fpos(2)=y3;atom.fpos(3)=0.0;                     // atom B13
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("0.0");// atom B13 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B13 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B13
    
    atom.name="C"; atom.type=2;                                       // atom B14
    atom.fpos(1)=-y3;atom.fpos(2)=(x3-y3);atom.fpos(3)=0.0;                     // atom B14
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("(x3-y3)");atom.fpos_equation.push_back("0.0");// atom B14 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B14 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B14
    
    atom.name="C"; atom.type=2;                                       // atom B15
    atom.fpos(1)=(-x3+y3);atom.fpos(2)=-x3;atom.fpos(3)=0.0;                     // atom B15
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x3+y3)");atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("0.0");// atom B15 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B15 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B15
    
    atom.name="C"; atom.type=2;                                       // atom B16
    atom.fpos(1)=-x3;atom.fpos(2)=-y3;atom.fpos(3)=0.0;                     // atom B16
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x3");atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("0.0");// atom B16 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B16 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B16
    
    atom.name="C"; atom.type=2;                                       // atom B17
    atom.fpos(1)=y3;atom.fpos(2)=(-x3+y3);atom.fpos(3)=0.0;                     // atom B17
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("(-x3+y3)");atom.fpos_equation.push_back("0.0");// atom B17 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B17 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B17
    
    atom.name="C"; atom.type=2;                                       // atom B18
    atom.fpos(1)=(x3-y3);atom.fpos(2)=x3;atom.fpos(3)=0.0;                     // atom B18
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x3-y3)");atom.fpos_equation.push_back("x3");atom.fpos_equation.push_back("0.0");// atom B18 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B18 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B18
    
    atom.name="C"; atom.type=2;                                       // atom B31
    atom.fpos(1)=x6;atom.fpos(2)=y6;atom.fpos(3)=(1.0/2.0);                     // atom B31
    atom.fpos_equation.clear();atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("(1.0/2.0)");// atom B31 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B31 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B31
    
    atom.name="C"; atom.type=2;                                       // atom B32
    atom.fpos(1)=-y6;atom.fpos(2)=(x6-y6);atom.fpos(3)=(1.0/2.0);                     // atom B32
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("(x6-y6)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B32 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B32 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B32
    
    atom.name="C"; atom.type=2;                                       // atom B33
    atom.fpos(1)=(-x6+y6);atom.fpos(2)=-x6;atom.fpos(3)=(1.0/2.0);                     // atom B33
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x6+y6)");atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("(1.0/2.0)");// atom B33 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B33 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B33
    
    atom.name="C"; atom.type=2;                                       // atom B34
    atom.fpos(1)=-x6;atom.fpos(2)=-y6;atom.fpos(3)=(1.0/2.0);                     // atom B34
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-x6");atom.fpos_equation.push_back("-y6");atom.fpos_equation.push_back("(1.0/2.0)");// atom B34 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B34 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B34
    
    atom.name="C"; atom.type=2;                                       // atom B35
    atom.fpos(1)=y6;atom.fpos(2)=(-x6+y6);atom.fpos(3)=(1.0/2.0);                     // atom B35
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y6");atom.fpos_equation.push_back("(-x6+y6)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B35 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B35 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B35
    
    atom.name="C"; atom.type=2;                                       // atom B36
    atom.fpos(1)=(x6-y6);atom.fpos(2)=x6;atom.fpos(3)=(1.0/2.0);                     // atom B36
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x6-y6)");atom.fpos_equation.push_back("x6");atom.fpos_equation.push_back("(1.0/2.0)");// atom B36 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B36 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B36
    

    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_ABC_hP36_175_jk_jk_jk(stringstream& web,bool LDEBUG) {
    #ifndef _ANRL_NOWEB_
    #endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_ABC_hP36_175_jk_jk_jk: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************