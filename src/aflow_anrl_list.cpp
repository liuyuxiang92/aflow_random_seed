// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo, David Hicks - 2016

#ifndef _AFLOW_ANRL_LIST_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_LIST_CPP // AFLOW_REMOVE_GREP
#include "aflow.h" // AFLOW_REMOVE_GREP  //CO20200521

// ***************************************************************************
namespace anrl { 
  uint PrototypeANRL_LoadList(vector<string>& vproto,
      vector<string>& vproto_label,
      vector<uint>& vproto_nspecies,
      vector<uint>& vproto_natoms,
      vector<uint>& vproto_spacegroup,
      vector<uint>& vproto_nunderscores,
      vector<uint>& vproto_nparameters,
      vector<string>& vproto_Pearson_symbol,
      vector<string>& vproto_params,
      vector<string>& vproto_Strukturbericht,
      vector<string>& vproto_prototype,
      vector<string>& vproto_dialect) {

    vproto.clear();
    vproto_label.clear();
    vproto_nspecies.clear();
    vproto_natoms.clear();
    vproto_spacegroup.clear();
    vproto_nunderscores.clear();
    vproto_nparameters.clear();
    vproto_Pearson_symbol.clear();
    vproto_params.clear();
    vproto_Strukturbericht.clear();
    vproto_prototype.clear();
    vproto_dialect.clear();


    //Label     # of Species    # atoms/primcell    #space_group_number   # Underscores     # Parameters    pearson_symbol    params    Strukturbericht     Prototype   Dialect        
    // -------------------------------------------------------------------------
    // Part 1
    // -------------------------------------------------------------------------
    vproto.push_back("AB2_aP12_1_4a_8a;2;12;1;4;42;aP12;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;FeS2;FeS2");
    vproto.push_back("ABC2_aP16_1_4a_4a_8a;3;16;1;5;54;aP16;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;AsKSe2;AsKSe2");
    vproto.push_back("A2B_aP6_2_2i_i;2;6;2;4;15;aP6;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;P2I4;P2I4");
    vproto.push_back("A_aP4_2_aci;1;4;2;3;9;aP4;a,b/a,c/a,alpha,beta,gamma,x3,y3,z3;-;Cf;Cf");
    vproto.push_back("A2B_mP12_3_bc3e_2e;2;12;3;4;21;mP12;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;SiO2;SiO2");
    vproto.push_back("A_mP4_4_2a;1;4;4;3;10;mP4;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2;-;Te;High-Pressure Te");
    vproto.push_back("A_mC12_5_3c;1;6;5;3;13;mC12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;A19;Po;Po");
    vproto.push_back("A3BC_mC10_8_ab_a_a;3;5;8;5;13;mC10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,y4,z4;-;Pb(Zr0.52Ti0.48)O3;Monoclinic PZT [PbO3]");
    vproto.push_back("A2B_mC144_9_24a_12a;2;72;9;4;112;mC144;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30,x31,y31,z31,x32,y32,z32,x33,y33,z33,x34,y34,z34,x35,y35,z35,x36,y36,z36;-;SiO2;Monoclinic Low Tridymite");
    vproto.push_back("AB_mP4_11_e_e;2;4;11;4;8;mP4;a,b/a,c/a,beta,x1,z1,x2,z2;-;NiTi;NiTi2");
    vproto.push_back("ABC3_mP10_11_e_e_ef;3;10;11;5;13;mP10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,y4,z4;G0_6;KClO3;KClO3");
    vproto.push_back("A_mP16_11_8e;1;16;11;3;20;mP16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;alphaPu;alphaPu");
    vproto.push_back("AB2_mC6_12_a_i;2;3;12;4;6;mC6;a,b/a,c/a,beta,x2,z2;C34;AuTe2;Calaverite");
    vproto.push_back("A_mC34_12_ah3i2j;1;17;12;3;17;mC34;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7;-;betaPu;betaPu");
    vproto.push_back("AB3_mC16_12_g_ij;2;8;12;4;10;mC16;a,b/a,c/a,beta,y1,x2,z2,x3,y3,z3;D0_15;AlCl3;AlCl3");
    vproto.push_back("A5B2_mC14_12_a2i_i;2;7;12;4;10;mC14;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4;-;Au5Mn2;Au5Mn2");
    vproto.push_back("A_mC4_12_i;1;2;12;3;6;mC4;a,b/a,c/a,beta,x1,z1;-;O2;alphaO2");
    vproto.push_back("ABC4_mP12_13_e_a_2g;3;12;13;5;11;mP12;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4;E1_b;AgAuTe4;Sylvanite");
    vproto.push_back("A_mP84_13_21g;1;84;13;3;67;mP84;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21;-;P;Monoclinic Phosphorus");
    vproto.push_back("A2B_mP12_14_2e_e;2;12;14;4;13;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;C43/-/-;ZrO2/O2W/N2Os1;Baddeleyite/O2W (ICSD #8217)/N2Os1 (ICSD #240759)"); //DX20200703 - added R. Friedrich binary oxide //DX20201019 - added metal nitride
    vproto.push_back("A_mP32_14_8e;1;32;14;3;28;mP32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;A_l;Se;betaSe");
    vproto.push_back("A_mP64_14_16e;1;64;14;3;52;mP64;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;A_k;Se;Se");
    vproto.push_back("A2B5_mC28_15_f_e2f;2;14;15;4;14;mC28;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-/-;B2Pd5/Sb2O5;B2Pd5/Sb2O5 (ICSD #1422)"); //DX20200703 - added R. Friedrich binary oxide
    vproto.push_back("AB_mC8_15_c_e;2;4;15;4;5;mC8;a,b/a,c/a,beta,y2;B26;CuO;Tenorite");
    vproto.push_back("A2B_mC48_15_ae3f_2f;2;24;15;4;20;mC48;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;SiO2;Coesite");
    vproto.push_back("ABC6D2_mC40_15_e_e_3f_f;4;20;15;6;18;mC40;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;CaFeO6Si2;Esseneite");
    vproto.push_back("ABC4_oP12_16_ag_cd_2u;3;12;16;5;9;oP12;a,b/a,c/a,x5,y5,z5,x6,y6,z6;-;AlPS4;AlPS4");
    vproto.push_back("AB3_oP16_18_ab_3c;2;16;18;4;14;oP16;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;D0_{17};BaS3;BaS3"); //DX 20180925 - moved Strukturberict designation from AB3_tP8_113_a_ce per M. Mehl's suggestion (historical reasons)
    vproto.push_back("A2B_oP12_19_2a_a;2;12;19;4;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;Ag2Se;Naumannite");
    vproto.push_back("A2B_oC24_20_abc_c;2;12;20;4;11;oC24;a,b/a,c/a,x1,y2,x3,y3,z3,x4,y4,z4;-;SiO2;Orthorhombic Tridymite");
    vproto.push_back("AB_oP2_25_b_a;2;2;25;4;5;oP2;a,b/a,c/a,z1,z2;-;CdTe;High-Pressure CdTe");
    vproto.push_back("AB2_oP24_28_acd_2c3d;2;24;28;4;22;oP24;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;C46;AuTe2;Krennerite");
    vproto.push_back("AB3C4_oP16_31_a_ab_2ab;3;16;31;5;17;oP16;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5,x6,y6,z6;H2_5;AsCu3S4;Enargite");
    vproto.push_back("AB_oP8_33_a_a;2;8;33;4;9;oP8;a,b/a,c/a,x1,y1,z1,x2,y2,z2;-;CoAs;Modderite");
    vproto.push_back("AB3C4_oP32_33_a_3a_4a;3;32;33;5;27;oP32;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;AsK3S4;AsK3S4");
    vproto.push_back("A2B_oC12_36_2a_a;2;6;36;4;9;oC12;a,b/a,c/a,y1,z1,y2,z2,y3,z3;C24;HgBr2;HgBr2");
    vproto.push_back("A2BC_oC8_38_e_a_b;3;4;38;5;7;oC8;a,b/a,c/a,z1,z2,y3,z3;-;C2CeNi;C2CeNi");
    vproto.push_back("A2B_oC12_38_de_ab;2;6;38;4;9;oC12;a,b/a,c/a,z1,z2,y3,z3,y4,z4;-;Au2V;Au2V");
    vproto.push_back("AB4_oC20_41_a_2b;2;10;41;4;10;oC20;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;D1_c;PtSn4;PtSn4");
    vproto.push_back("AB2_oC24_41_2a_2b;2;12;41;4;11;oC24;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4;C_e;PdSn2;PdSn2");
    vproto.push_back("AB2_oF72_43_ab_3b;2;18;43;4;16;oF72;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;C44;GeS2;GeS2");
    vproto.push_back("AB_oI4_44_a_b;2;2;44;4;5;oI4;a,b/a,c/a,z1,z2;-;AsGa;High-pressure GaAs");
    vproto.push_back("A2B3C7D_oP13_47_t_aq_eqrs_h;4;13;47;6;8;oP13;a,b/a,c/a,z4,z5,z6,z7,z8;-;YBa2Cu3O7-x;1212C [YBa2Cu3O7-x]");
    vproto.push_back("AB_oP4_51_e_f;2;4;51;4;5;oP4;a,b/a,c/a,z1,z2;B19;AuCd;beta'-AuCd");
    vproto.push_back("A3B2_oP20_56_ce_e;2;20;56;4;10;oP20;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;D5_11;Sb2O3;Sb2O3");
    vproto.push_back("ABCD_oP16_57_d_c_d_d;4;16;57;6;10;oP16;a,b/a,c/a,x1,x2,y2,x3,y3,x4,y4;F5_9;KCNS;KCNS");
    vproto.push_back("AB_oP8_57_d_d;2;8;57;4;7;oP8;a,b/a,c/a,x1,y1,x2,y2;-;TlF;TlF-II, 57");
    vproto.push_back("AB2_oP6_58_a_g;2;6;58;4;5;oP6;a,b/a,c/a,x2,y2;C35/-/C18/-/-;CaCl2/etaFe2C/FeS2/Ir1N2/Au1N2;Hydrophilite/eta-Fe2C/Marcasite/Ir1N2 (ICSD #160620)/Au1N2 (ICSD #166465)"); //DX20201019 - added two metal nitrides
    vproto.push_back("AB_oP4_59_a_b;2;4;59;4;5;oP4;a,b/a,c/a,z1,z2;-;CuTe;Vulcanite");
    vproto.push_back("ABC_oP6_59_a_a_a;3;6;59;5;6;oP6;a,b/a,c/a,z1,z2,z3;-;CNCl;CNCl");
    vproto.push_back("A3B_oP8_59_bf_a;2;8;59;4;7;oP8;a,b/a,c/a,z1,z2,x3,z3;D0_a;TiCu3;betaTiCu3");
    vproto.push_back("AB_oP16_61_c_c;2;16;61;4;9;oP16;a,b/a,c/a,x1,y1,z1,x2,y2,z2;B_e;CdSb;CdSb");
    vproto.push_back("A2B_oP24_61_2c_c;2;24;61;4;12;oP24;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;C21;TiO2;Brookite");
    vproto.push_back("A3B2_oP20_62_3c_2c;2;20;62;4;13;oP20;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;D5_8;Sb2S3;Stibnite");
    vproto.push_back("AB3C_oP20_62_c_cd_a;3;20;62;5;10;oP20;a,b/a,c/a,x2,z2,x3,z3,x4,y4,z4;-/-;CaTiO3/CaO3Zr;CaTiO3 Pnma Perovskite/CaO3Zr (ICSD #97463)"); //DX20200703 - added R. Friedrich ternary oxide
    vproto.push_back("A4B_oP20_62_2cd_c;2;20;62;4;12;oP20;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,y4,z4;-;MgB4;MgB4");
    vproto.push_back("AB2C_oP16_62_c_2c_c;3;16;62;5;11;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4;F5_6;CuSbS2;Chalcostibite");
    vproto.push_back("A2B_oP12_62_2c_c;2;12;62;4;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;C37/C25/C23/C29/-;Co2Si/HgCl2/PbCl2/SrH2/N2Re1;Co2Si/HgCl2/Cotunnite/SrH2/N2Re1 (ICSD #187446)"); //added info from part 2 A2B_oP12_62_2c_c;2;12;62;4;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;C29;SrH2;SrH2 //DX20201019 - added metal-nitride
    vproto.push_back("AB_oP8_62_c_c;2;8;62;4;7;oP8;a,b/a,c/a,x1,z1,x2,z2;B16/B31/B27/B29/B14;GeS/MnP/FeB/SnS/FeAs;GeS/MnP/FeB/SnS/Westerveldite"); //added info from part 2 AB_oP8_62_c_c;2;8;62;4;7;oP8;a,b/a,c/a,x1,z1,x2,z2;B14;FeAs;Westerveldite
    vproto.push_back("AB3_oP16_62_c_cd;2;16;62;4;10;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3;D0_11/-;Fe3C/N1Na3;Cementite/N1Na3 (ICSD #165990)"); //DX20201019 - added metal-nitride
    vproto.push_back("A3B7_oP40_62_cd_3c2d;2;40;62;4;20;oP40;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;D10_1;C3Cr7;C3Cr7");
    vproto.push_back("A_oP8_62_2c;1;8;62;3;7;oP8;a,b/a,c/a,x1,z1,x2,z2;A_c;alphaNp;alphaNp");
    vproto.push_back("AB2C_oC16_63_c_2c_c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,y4;-/-/-/-;SrCuO2/Au1Ca2N1/Cu1Sn2Sr1/Co1Sn2Tb1;SrCuO2/Au1Ca2N1 (ICSD #85528)/Cu1Sn2Sr1 (ICSD #103110)/Co1Sn2Tb1 (ICSD #240096)"); //DX20201016 - added metal-nitride //DX20201028 - added two metals
    vproto.push_back("A2B_oC12_63_2c_c;2;6;63;4;6;oC12;a,b/a,c/a,y1,y2,y3;C49;ZrSi2;ZrSi2");
    vproto.push_back("AB_oC8_63_c_c;2;4;63;4;5;oC8;a,b/a,c/a,y1,y2;B33/-/-;CrB/-/Al1N1;CrB/-/Al1N1 (ICSD #163951)"); //DX20201019 - added metal-nitride
    vproto.push_back("A_oC4_63_c;1;2;63;3;4;oC4;a,b/a,c/a,y1;A20;alpha-U;alpha-Uranium");
    vproto.push_back("A_oC8_64_f;1;4;64;3;5;oC8;a,b/a,c/a,y1,z1;A11/A17/A14;alphaGa/P/I2;alphaGallium/Black Phosphorus/Molecular Iodine");
    vproto.push_back("A2B2C_oC80_64_efg_efg_df;3;40;64;5;18;oC80;a,b/a,c/a,x1,y2,y3,y4,z4,y5,z5,y6,z6,x7,y7,z7,x8,y8,z8;-;MgB2C2;MgB2C2");
    vproto.push_back("AB_oC8_65_j_g;2;4;65;4;5;oC8;a,b/a,c/a,x1,y2;-;alphaIrV;alphaIrV");
    vproto.push_back("A3B5_oC16_65_ah_bej;2;8;65;4;5;oC16;a,b/a,c/a,x4,y5;-;Ga3Pt5;Ga3Pt5");
    vproto.push_back("AB3_oC8_65_a_bf;2;4;65;4;3;oC8;a,b/a,c/a;L1_3;CdPt3;Predicted CdPt3");
    vproto.push_back("AB_oF8_69_a_b;2;2;69;4;3;oF8;a,b/a,c/a;B24;TlF;TlF");
    vproto.push_back("A_oF8_70_a;1;2;70;3;3;oF8;a,b/a,c/a;-;gammaPu;gamma-Pu");
    vproto.push_back("A2B_oF24_70_e_a;2;6;70;4;4;oF24;a,b/a,c/a,x2;C54;TiSi2;TiSi2");
    vproto.push_back("A_oF128_70_4h;1;32;70;3;15;oF128;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;A16;S;alpha-Sulfur");
    vproto.push_back("AB2_oI6_71_a_i;2;3;71;4;4;oI6;a,b/a,c/a,z2;-;ReSi2;ReSi2");
    vproto.push_back("AB2_oI6_71_a_g;2;3;71;4;4;oI6;a,b/a,c/a,y2;-;MoPt2;MoPt2");
    vproto.push_back("A2B_oI12_72_j_a;2;6;72;4;5;oI12;a,b/a,c/a,x2,y2;C42;SiS2;SiS2");
    vproto.push_back("AB4C_tI12_82_c_g_a;3;6;82;5;5;tI12;a,c/a,x3,y3,z3;H0_7;BPO4;BPO4");
    vproto.push_back("A2BC4_tI14_82_bc_a_g;3;7;82;5;5;tI14;a,c/a,x4,y4,z4;E3;CdAl2S4;CdAl2S4");
    vproto.push_back("AB_tP16_84_cej_k;2;16;84;4;7;tP16;a,c/a,x3,y3,x4,y4,z4;B34;PdS;PdS");
    vproto.push_back("A4B5_tI18_87_h_ah;2;9;87;4;6;tI18;a,c/a,x2,y2,x3,y3;-;Ti5Te4;Ti5Te4");
    vproto.push_back("AB4_tI10_87_a_h;2;5;87;4;4;tI10;a,c/a,x2,y2;D1_a;Ni4Mo;Ni4Mo");
    vproto.push_back("A2B_tP12_92_b_a;2;12;92;4;6;tP12;a,c/a,x1,x2,y2,z2;-;SiO2;alpha Cristobalite");
    vproto.push_back("A2B_tP36_96_3b_ab;2;36;96;4;15;tP36;a,c/a,x1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;SiO2;Keatite");
    vproto.push_back("A_tP12_96_ab;1;12;96;3;6;tP12;a,c/a,x1,x2,y2,z2;-;Si;\"ST12'' of Si");
    vproto.push_back("A3BC_tP5_99_bc_a_b;3;5;99;5;6;tP5;a,c/a,z1,z2,z3,z4;-;Pb(Zr0.52Ti0.48)O3;Tetragonal PZT [PbO3]");
    vproto.push_back("AB3_tP8_113_a_ce;2;8;113;4;5;tP8;a,c/a,z2,x3,z3;-;BaS3;BaS3");//DX 20180925 - moved Strukturberict designation to AB3_oP16_18_ab_3c per M. Mehl's suggestion (historical reasons)
    vproto.push_back("A2BC4D_tI16_121_d_a_i_b;4;8;121;6;4;tI16;a,c/a,x4,z4;H2_6;Cu2FeS4Sn;Stannite");
    vproto.push_back("ABC2_tI16_122_a_b_d;3;8;122;5;3;tI16;a,c/a,x3;E1_1;CuFeS2;Chalcopyrite");
    vproto.push_back("AB5C_tP7_123_b_ci_a;3;7;123;5;3;tP7;a,c/a,z4;-;HoCoGa5;HoCoGa5");
    vproto.push_back("AB3_tP4_123_a_ce;2;4;123;4;2;tP4;a,c/a;L6_0;CuTi3;CuTi3");
    vproto.push_back("AB_tP2_123_a_d;2;2;123;4;2;tP2;a,c/a;L1_0;CuAu;CuAu");
    vproto.push_back("ABC2_tP4_123_d_a_f;3;4;123;5;2;tP4;a,c/a;-;CaCuO2;CaCuO2");
    vproto.push_back("A2B3_tP10_127_g_ah;2;10;127;4;4;tP10;a,c/a,x2,x3;D5_a;Si2U3;Si2U3");
    vproto.push_back("ABCD_tP8_129_c_b_a_c;4;8;129;6;4;tP8;a,c/a,z3,z4;-;AsCuSiZr;AsCuSiZr");
    vproto.push_back("A_tP4_129_ac;1;4;129;3;3;tP4;a,c/a,z2;A_d;betaNp;betaNp");
    vproto.push_back("ABC_tP6_129_c_a_c;3;6;129;5;4;tP6;a,c/a,z2,z3;E0_1;PbFCl;Matlockite");
    vproto.push_back("A2B_tP6_129_ac_c;2;6;129;4;4;tP6;a,c/a,z2,z3;C38;Cu2Sb;Cu2Sb");
    vproto.push_back("AB_tP4_129_a_c;2;4;129;4;3;tP4;a,c/a,z2;B10;PbO;PbO");
    vproto.push_back("AB_tP4_129_c_c;2;4;129;4;4;tP4;a,c/a,z1,z2;B11/-;gammaCuTi/-;gammaCuTi/-");
    vproto.push_back("AB_tP4_131_c_e;2;4;131;4;2;tP4;a,c/a;B17/-;PtS/PdO;PtS/PdO (ICSD #26598)"); //DX20200703 - added R. Friedrich binary oxide
    vproto.push_back("A_tP50_134_b2m2n;1;50;134;3;12;tP50;a,c/a,x2,z2,x3,z3,x4,y4,z4,x5,y5,z5;A_g;B;T-50 Boron");
    vproto.push_back("A_tP30_136_bf2ij;1;30;136;3;9;tP30;a,c/a,x2,x3,y3,x4,y4,x5,z5;A_b;betaU;betaU");
    vproto.push_back("AB_tP8_136_g_f;2;8;136;4;4;tP8;a,c/a,x1,x2;-;betaBeO;betaBeO");
    vproto.push_back("A2B_tP6_136_f_a;2;6;136;4;3;tP6;a,c/a,x2;C4;TiO2;CrSi2");
    //vproto.push_back("sigma_tP30_136_bf2ij;5;30;136;3;9;tP30;a,c/a,x2,x3,y3,x4,y4,x5,z5;D8_b;CrFe;sigmaCrFe");
    vproto.push_back("sigma_tP30_136_bf2ij;1;30;136;3;9;tP30;a,c/a,x2,x3,y3,x4,y4,x5,z5;D8_b;CrFe;sigmaCrFe");
    vproto.push_back("A_tP4_136_f;1;4;136;3;3;tP4;a,c/a,x1;-;N2;gammaN2");
    vproto.push_back("A_tP16_138_j;1;16;138;3;5;tP16;a,c/a,x1,y1,z1;A18;Cl2;Cl2");
    vproto.push_back("A3B_tI16_139_cde_e;2;8;139;4;4;tI16;a,c/a,z3,z4;D0_23;Al3Zr;Al3Zr");
    vproto.push_back("A_tI4_139_e;1;2;139;3;3;tI4;a,c/a,z1;-;Si;Hypothetical BCT5 Si");
    vproto.push_back("AB2C4_tI14_139_a_e_ce;3;7;139;5;4;tI14;a,c/a,z3,z4;-/-;(La,Ba)2CuO4/O4Sr2Ti;0201 [(La,Ba)2CuO4]/O4Sr2Ti (ICSD #157402)"); //DX20200703 - added R. Friedrich ternary oxide
    vproto.push_back("A12B_tI26_139_fij_a;2;13;139;4;4;tI26;a,c/a,x3,x4;D2_b;Mn12Th;Mn12Th");
    vproto.push_back("A_tI2_139_a;1;1;139;3;2;tI2;a,c/a;A6/A_a;In/Pa;In/alphaPa");
    vproto.push_back("A_tI8_139_h;1;4;139;3;3;tI8;a,c/a,x1;-;C;Hypothetical Tetrahedrally Bonded Carbon with 4-Member Rings");
    vproto.push_back("A3B_tI8_139_bd_a;2;4;139;4;2;tI8;a,c/a;D0_22;Al3Ti;Al3Ti");
    vproto.push_back("AB2_tI6_139_a_e;2;3;139;4;3;tI6;a,c/a,z2;C11_b/-/-;MoSi2/-/KO2;MoSi2/-/KO2 (ICSD #38245)"); //DX20200703 - added R. Friedrich binary oxide
    vproto.push_back("A4B5_tI18_139_i_ah;2;9;139;4;4;tI18;a,c/a,x2,x3;-;V4Zn5;V4Zn5");
    vproto.push_back("A4B_tI10_139_de_a;2;5;139;4;3;tI10;a,c/a,z3;D1_3;Al4Ba;Al4Ba");
    vproto.push_back("A8B_tI18_139_hi_a;2;9;139;4;4;tI18;a,c/a,x2,x3;-;Pt8Ti;Pt8Ti");
    vproto.push_back("A2B_tI6_139_d_a;2;3;139;4;2;tI6;a,c/a;L\'2;ThH2;ThH2");
    vproto.push_back("A2B_tI12_140_h_a;2;6;140;4;3;tI12;a,c/a,x2;C16;Al2Cu;Khatyrkite");
    vproto.push_back("AB3_tI16_140_b_ah;2;8;140;4;3;tI16;a,c/a,x3;D0_c;SiU3;SiU3");
    vproto.push_back("AB_tI16_140_ab_h;2;8;140;4;3;tI16;a,c/a,x3;B37;SeTl;SeTl");
    vproto.push_back("A4BC_tI24_141_h_b_a;3;12;141;5;4;tI24;a,c/a,y3,z3;-;ZrSiO4;Zircon");
    vproto.push_back("A_tI4_141_a;1;2;141;3;2;tI4;a,c/a;A5;Sn;betaSn");
    vproto.push_back("A3B4_tI28_141_ad_h;2;14;141;4;4;tI28;a,c/a,y3,z3;-;Mn3O4;Hausmannite");
    vproto.push_back("A2B_tI12_141_e_a;2;6;141;4;3;tI12;a,c/a,z2;C5/C_{c};TiO2/alpha-ThSi2;Anatase/alpha-ThSi2"); // added info from part 2 - A2B_tI12_141_e_a;2;6;141;4;3;tI12;a,c/a,z2;C_{c};alpha-ThSi2;alpha-ThSi2
    vproto.push_back("AB_tI16_141_e_e;2;8;141;4;4;tI16;a,c/a,z1,z2;B_g;MoB;MoB");
    vproto.push_back("A2B_tI24_141_2e_e;2;12;141;4;5;tI28;a,c/a,z1,z2,z3;-;Ga2Hf;Ga2Hf");
    vproto.push_back("AB_tI8_141_a_b;2;4;141;4;2;tI8;a,c/a;\"40\"/-;NbP/-;NbP/-");
    vproto.push_back("A2B3_tI80_141_ceh_3h;2;40;141;4;11;tI80;a,c/a,z2,y3,z3,y4,z4,y5,z5,y6,z6;-;In2S3;betaIn2S3");
    vproto.push_back("ABC4_tI96_142_e_ab_2g;3;48;142;5;9;tI96;a,c/a,x3,x4,y4,z4,x5,y5,z5;-;PPrS4;PPrS4");
    vproto.push_back("A2B_hP9_147_g_ad;2;9;147;4;6;hP9;a,c/a,z2,x3,y3,z3;B_b;AgZn;zetaAgZn");
    vproto.push_back("AB_hR16_148_cf_cf;2;16;148;4;10;hR16;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4;-;C8H8;Solid Cubane");
    vproto.push_back("AB3_hR8_148_c_f;2;8;148;4;6;hR8;a,c/a,x1,x2,y2,z2;D0_5;BiI3;BiI3");
    vproto.push_back("AB_hR26_148_b2f_a2f;2;26;148;4;14;hR26;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;PdAl;PdAl");
    vproto.push_back("AB3C_hR10_148_c_f_c;3;10;148;5;7;hR10;a,c/a,x1,x2,x3,y3,z3;-;FeTiO3;Ilmenite");
    vproto.push_back("A2B_hP9_150_ef_bd;2;9;150;4;5;hP9;a,c/a,z2,x3,x4;C22;Fe2P;Original Fe2P");
    vproto.push_back("A3B_hP24_151_3c_2a;2;24;151;4;13;hP24;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5;D0_4;CrCl3;CrCl3");
    vproto.push_back("A2B_hP9_152_c_a;2;9;152;4;6;hP9;a,c/a,x1,x2,y2,z2;-;SiO2;alphaQuartz");
    vproto.push_back("A_hP3_152_a;1;3;152;3;3;hP3;a,c/a,x1;A8;Se;gammaSe");
    vproto.push_back("AB_hP6_154_a_b;2;6;154;4;4;hP6;a,c/a,x1,x2;B9;HgS;Cinnabar");
    vproto.push_back("AB3_hR8_155_c_de;2;8;155;4;5;hR8;a,c/a,x1,y2,y3;D0_14;AlF3;AlF3");
    vproto.push_back("A3B2_hR5_155_e_c;2;5;155;4;4;hR5;a,c/a,x1,y2;D5_e;Ni3S2;Hazelwoodite");
    vproto.push_back("AB_hR6_160_b_b;2;6;160;4;6;hR6;a,c/a,x1,z1,x2,z2;B13;NiS;Millerite");
    vproto.push_back("AB_hR6_160_3a_3a;2;6;160;4;8;hR6;a,c/a,x1,x2,x3,x4,x5,x6;-;CSi;Moissanite 9R structure");
    vproto.push_back("ABC3_hR10_161_a_a_b;3;10;161;5;7;hR10;a,c/a,x1,x2,x3,y3,z3;-;LiNbO3;Ferroelectric LiNbO3");
    vproto.push_back("AB2_hP9_162_ad_k;2;9;162;4;4;hP9;a,c/a,x3,z3;-;NV2;betaV2N");
    vproto.push_back("AB2CD2_hP36_163_h_i_bf_i;4;36;163;6;10;hP36;a,c/a,z2,x3,x4,y4,z4,x5,y5,z5;F5_10;KAg(CN)2;KAg2");
    vproto.push_back("A3B2_hP5_164_ad_d;2;5;164;4;4;hP5;a,c/a,z2,z3;D5_13/-/-/-;Al3Ni2/Ca3N2/Ca3N2/Ca3N2;Al3Ni2/Ca3N2 (ICSD #162795)/Ca3N2 (ICSD #162796)/Ca3N2 (ICSD #169727)"); //DX20201019 - added three metal-nitrides
    vproto.push_back("AB2_hP3_164_a_d;2;3;164;4;3;hP3;a,c/a,z2;C6;CdI2;omega Phase");
    vproto.push_back("A3B_hP24_165_adg_f;2;24;165;4;7;hP24;a,c/a,z2,x3,x4,y4,z4;-;H3Ho;H3Ho");
    vproto.push_back("AB_hR2_166_a_b;2;2;166;4;2;hR2;a,c/a;L1_1;CuPt;CuPt");
    vproto.push_back("A_hR2_166_c;1;2;166;3;3;hR2;a,c/a,x1;A7/-/-;As/C/O2;alphaAs/Rhombohedral Graphite/betaO2");
    vproto.push_back("A_hR1_166_a;1;1;166;3;2;hR1;a,c/a;A_i/A10;Po/Hg;betaPo/alphaHg");
    vproto.push_back("A7B6_hR13_166_ah_3c;2;13;166;4;7;hR13;a,c/a,x2,x3,x4,x5,z5;D8_5;Fe7W6;Fe7W6 mu-phase");
    vproto.push_back("A_hR3_166_ac;1;3;166;3;3;hR3;a,c/a,x2;C19;Sm;alphaSm");
    vproto.push_back("A2B3_hR5_166_c_ac;2;5;166;4;4;hR5;a,c/a,x2,x3;C33;Bi2Te3;Bi2Te3");
    vproto.push_back("A5B2_hR7_166_a2c_c;2;7;166;4;5;hR7;a,c/a,x2,x3,x4;D8_i;B5Mo2;Mo2B5");
    vproto.push_back("A_hR12_166_2h;1;12;166;3;6;hR12;a,c/a,x1,z1,x2,z2;-;B;alphaBoron");
    vproto.push_back("ABC2_hR4_166_a_b_c;3;4;166;5;3;hR4;a,c/a,x3;F5_1/-/-/-;CrNiS2/-/AlLiO2/CuFeO2;Caswellsilverite/-/AlLiO2 (ICSD #28288)/CuFeO2 (ICSD #246912)"); //DX20200703 - added R. Friedrich ternary oxides (2)
    vproto.push_back("A_hR105_166_bc9h4i;1;105;166;3;33;hR105;a,c/a,x2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;B;betaBoron");
    vproto.push_back("A6B_hR7_166_g_a;2;7;166;4;3;hR7;a,c/a,x2;-;CaC6;CaC6");
    vproto.push_back("ABC3_hR10_167_a_b_e;3;10;167;5;3;hR10;a,c/a,x3;-/G0_1;LiNbO3/CaCO3;Paraelectric LiNbO3/Calcite");
    vproto.push_back("A2B3_hR10_167_c_e;2;10;167;4;4;hR10;a,c/a,x1,x2;D5_1/-;Al2O3/O3V2;Corundum/O3V2 (ICSD #94768)"); //DX20200703 - added R. Friedrich binary oxide
    vproto.push_back("A2B_hP18_180_fi_bd;2;18;180;4;4;hP18;a,c/a,z3,x4;C_a;Mg2Ni;Mg2Ni");
    vproto.push_back("AB2_hP9_180_d_j;2;9;180;4;3;hP9;a,c/a,x2;C40;CrSi2;CrSi2");
    vproto.push_back("A2B_hP9_180_j_c;2;9;180;4;3;hP9;a,c/a,x2;C8;SiO2;beta-Quartz");
    vproto.push_back("AB3_hP8_182_c_g;2;8;182;4;3;hP8;a,c/a,x2;-;Fe3C;Bainite");
    vproto.push_back("A_hP4_186_ab;1;4;186;3;4;hP4;a,c/a,z1,z2;-;C;Buckled Graphite");
    vproto.push_back("AB_hP8_186_ab_ab;2;8;186;4;6;hP8;a,c/a,z1,z2,z3,z4;B5;SiC;Moissanite-4H SiC");
    vproto.push_back("AB_hP4_186_b_b;2;4;186;4;4;hP4;a,c/a,z1,z2;B4/-/-/-;ZnS/Ga1N1/La1N1/Mo1N1;Wurtzite/Ga1N1 (ICSD #153888)/La1N1 (ICSD #162195)/Mo1N1 (ICSD #168369)"); //DX20201019 - added three metal-nitride
    vproto.push_back("AB_hP12_186_a2b_a2b;2;12;186;4;8;hP12;a,c/a,z1,z2,z3,z4,z5,z6;B6;SiC;Moissanite-6H SiC");
    vproto.push_back("A5B3C_hP18_186_2a3b_2ab_b;3;18;186;5;11;hP18;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8,z9;E9_4;Al5C3N;Al5C3N");
    vproto.push_back("AB_hP4_186_b_a;2;4;186;4;4;hP4;a,c/a,z1,z2;B12/-;BN/Ga1N1;Original BN/Ga1N1 (ICSD #159250)"); //DX20201019 - added metal-nitride
    vproto.push_back("ABC_hP3_187_a_d_f;3;3;187;5;2;hP3;a,c/a;-;BaPtSb;BaPtSb");
    vproto.push_back("AB_hP2_187_d_a;2;2;187;4;2;hP2;a,c/a;B_h;WC;Tungsten Carbide");
    vproto.push_back("A2B_hP9_189_fg_bc;2;9;189;4;4;hP9;a,c/a,x3,x4;C22;Fe2P;Revised Fe2P");
    vproto.push_back("AB4C_hP6_191_a_h_b;3;6;191;5;3;hP6;a,c/a,z3;-;AlB4Mg;AlB4Mg");
    vproto.push_back("AB5_hP6_191_a_cg;2;6;191;4;2;hP6;a,c/a;D2_d;CaCu5;CaCu5");
    vproto.push_back("A_hP1_191_a;1;1;191;3;2;hP1;a,c/a;A_f;gammaHgSn6-10;Simple Hexagonal Lattice");
    vproto.push_back("A3B_hP4_191_bc_a;2;4;191;4;2;hP4;a,c/a;-;Li3N;Li3 Ni");
    vproto.push_back("AB2_hP3_191_a_d;2;3;191;4;2;hP3;a,c/a;C32;AlB2;Hexagonal omega");
    vproto.push_back("A2B_hP6_191_h_e;2;6;191;4;4;hP6;a,c/a,z1,z2;C_h;Cu2Te;Cu2Te");
    vproto.push_back("AB_hP6_191_f_ad;2;6;191;4;2;hP6;a,c/a;B35/-;CoSn/N1Ta1;CoSn/N1Ta1 (ICSD #25659)"); //DX20201019 - added metal-nitride
    vproto.push_back("AB_hP8_194_ad_f;2;8;194;4;3;hP8;a,c/a,z3;B_i;AsTi;AsTi");
    vproto.push_back("A_hP6_194_h;1;6;194;3;3;hP6;a,c/a,x1;-;C;Hypothetical Tetrahedrally Bonded Carbon with 3-Member Rings");
    vproto.push_back("AB_hP12_194_af_bf;2;12;194;4;4;hP12;a,c/a,z3,z4;-;CMo;CMo structure");
    vproto.push_back("A_hP4_194_ac;1;4;194;3;2;hP4;a,c/a;A3\';La;alphaLa");
    vproto.push_back("AB3_hP8_194_c_bf;2;8;194;4;3;hP8;a,c/a,z3;D0_18;Na3As;Na3As");
    vproto.push_back("AB2_hP6_194_b_f;2;6;194;4;3;hP6;a,c/a,z2;-;CaIn2;CaIn2");
    vproto.push_back("AB_hP4_194_c_d;2;4;194;4;2;hP4;a,c/a;B_k/-;BN/Al1N1;BN/Al1N1 (ICSD #163950)"); //DX20201019 - added metal-nitride
    vproto.push_back("ABC2_hP8_194_d_a_f;3;8;194;5;3;hP8;a,c/a,z3;-;AlCCr2;AlCCr2");
    vproto.push_back("A3B_hP8_194_h_c;2;8;194;4;3;hP8;a,c/a,x2;D0_19;Ni3Sn;Ni3Sn");
    vproto.push_back("A_hP4_194_bc;1;4;194;3;2;hP4;a,c/a;A9;C;Hexagonal Graphite structure");
    vproto.push_back("AB2_hP6_194_c_f;2;6;194;4;3;hP6;a,c/a,z2;C7;MoS2;Molybdenite");
    vproto.push_back("A5B2_hP14_194_abdf_f;2;14;194;4;4;hP14;a,c/a,z4,z5;D8_h;W2B5;W2B5");
    vproto.push_back("AB2_hP12_194_f_ah;2;12;194;4;4;hP12;a,c/a,z2,x3;C14;MgZn2;MgZn2 Hexagonal Laves");
    vproto.push_back("ABC_hP6_194_c_d_a;3;6;194;5;2;hP6;a,c/a;-/-;LiBC/Hg1Pb1Yb1;LiBC/Hg1Pb1Yb1 (ICSD #108545)"); //DX20201028 - added metal ternary
    vproto.push_back("A_hP4_194_f;1;4;194;3;3;hP4;a,c/a,z1;-;C;Lonsdaleite");
    vproto.push_back("AB2_hP6_194_c_ad;2;6;194;4;2;hP6;a,c/a;B8_2;Ni2In;Ni2In");
    vproto.push_back("AB3C4_hP16_194_c_af_ef;3;16;194;5;5;hP16;a,c/a,z3,z4,z5;-;AlN3Ti4;AlN3Ti4");
    vproto.push_back("A_hP2_194_c;1;2;194;3;2;hP2;a,c/a;A3;Mg;Hexagonal Close Packed");
    vproto.push_back("AB2_hP24_194_ef_fgh;2;24;194;4;6;hP24;a,c/a,z1,z2,z3,x5;C36;MgNi2;MgNi2 Hexagonal Laves");
    vproto.push_back("AB_hP12_194_df_ce;2;12;194;4;4;hP12;a,c/a,z3,z4;B18;CuS;Covellite");
    vproto.push_back("AB_hP4_194_c_a;2;4;194;4;2;hP4;a,c/a;B8_1;NiAs;NiAs");
    vproto.push_back("A2B_hP12_194_cg_f;2;12;194;4;3;hP12;a,c/a,z2;C10;SiO2;beta-Tridymite");
    vproto.push_back("A4B_cI40_197_cde_c;2;20;197;4;5;cI40;a,x1,x2,x3,x4;-;Ga4Ni;Ga4Ni3");
    vproto.push_back("ABC_cP12_198_a_a_a;3;12;198;5;4;cP12;a,x1,x2,x3;F0_1;NiSSb;Ullmanite");
    vproto.push_back("A3B_cP16_198_b_a;2;16;198;4;5;cP16;a,x1,x2,y2,z2;D1;NH3;Ammonia");
    vproto.push_back("A_cP8_198_2a;1;8;198;3;3;cP8;a,x1,x2;-;N2;alpha-N2");
    vproto.push_back("AB_cP8_198_a_a;2;8;198;4;3;cP8;a,x1,x2;B21/B20;CO/FeSi;alphaCO/FeSi");
    vproto.push_back("AB_cI16_199_a_a;2;8;199;4;3;cI16;a,x1,x2;B_a;CoU;CoU");
    vproto.push_back("AB32C48_cI162_204_a_2efg_2gh;3;81;204;5;13;cI162;a,x2,x3,x4,y5,z5,y6,z6,y7,z7,x8,y8,z8;-;Mg32(Al,Zn)49;Bergman [Mg32(Al,Zn)49]");
    vproto.push_back("A3B_cI32_204_g_c;2;16;204;4;3;cI32;a,y2,z2;D0_2;CoAs3;Skutterudite");
    vproto.push_back("A12B_cI26_204_g_a;2;13;204;4;3;cI26;a,y2,z2;-;Al12W;Al12W");
    vproto.push_back("A_cP8_205_c;1;8;205;3;2;cP8;a,x1;-;N2;alpha-N2");
    vproto.push_back("AB_cP16_205_c_c;2;16;205;4;3;cP16;a,x1,x2;-;CuCl;SC16");
    vproto.push_back("AB2_cP12_205_a_c;2;12;205;4;2;cP12;a,x2;C2/-;FeS2/NaO2;Pyrite/NaO2 (ICSD #87178"); //DX20200703 - added binary oxide for R. Friedrich 
    vproto.push_back("AB3C6_cI80_206_b_d_e;3;40;206;5;5;cI80;a,x2,x3,y3,z3;D5_3;(Mn,Fe)2O3;Bixbyite");
    vproto.push_back("A_cI16_206_c;1;8;206;3;2;cI16;a,x1;-;Si;BC8");
    vproto.push_back("A_cP20_213_cd;1;20;213;3;3;cP20;a,x1,y2;A13;Mn;betaMn");
    vproto.push_back("A3B4C_cP8_215_d_e_a;3;8;215;5;2;cP8;a,x3;H2_4;Cu3S4V;Sulvanite");
    vproto.push_back("AB4_cP5_215_a_e;2;5;215;4;2;cP5;a,x2;-;Fe4C;Fe4C");
    vproto.push_back("AB3C4_cP8_215_a_c_e;3;8;215;5;2;cP8;a,x3;-;Cu3AsS4;Cubic Lazarevicite");
    vproto.push_back("AB5_cF24_216_a_ce;2;6;216;4;2;cF24;a,x3;C15_b;AuBe5;AuBe5");
    vproto.push_back("ABC_cF12_216_b_c_a;3;3;216;5;1;cF12;a;C1_b;AgAsMg;Half-Heusler");
    vproto.push_back("AB_cF8_216_c_a;2;2;216;4;1;cF8;a;B3;ZnS;NaTl");
    vproto.push_back("A4B_cI10_217_c_a;2;5;217;4;2;cI10;a,x2;-;SiF4;SiF4");
    vproto.push_back("A_cI58_217_ac2g;1;29;217;3;6;cI58;a,x2,x3,z3,x4,z4;A12;Mn;alphaMn");
    vproto.push_back("A5B8_cI52_217_ce_cg;2;26;217;4;6;cI52;a,x1,x2,x3,x4,z4;-;Cu5Zn8;gammaBrass");
    vproto.push_back("A_cI16_220_c;1;8;220;3;2;cI16;a,x1;-;Li;High Pressure cI16 Li");
    vproto.push_back("A3B2_cI40_220_d_c;2;20;220;4;3;cI40;a,x1,x2;D5_c;Pu2C3;Pu2C3");
    vproto.push_back("AB_cP2_221_b_a;2;2;221;4;1;cP2;a;B2;CsCl;CsCl");
    vproto.push_back("AB_cP6_221_c_d;2;6;221;4;1;cP6;a;-;NbO;NbO");
    vproto.push_back("AB3C_cP5_221_a_c_b;3;5;221;5;1;cP5;a;E2_1;CaTiO3;Cubic Perovskite");
    vproto.push_back("AB27CD3_cP32_221_a_dij_b_c;4;32;221;6;3;cP32;a,y5,y6;-;CrFe2525Ni6;Model of Austenite");
    vproto.push_back("AB3_cP4_221_a_c;2;4;221;4;1;cP4;a;L1_2;Cu3Au;Cu3Au");
    vproto.push_back("A_cP1_221_a;1;1;221;3;1;cP1;a;A_h;Po;alphaPo");
    vproto.push_back("AB11_cP36_221_c_agij;2;36;221;4;4;cP36;a,x3,y4,y5;D2_e;BaHg11;BaHg11");
    vproto.push_back("AB11CD3_cP16_221_a_dg_b_c;4;16;221;6;2;cP16;a,x5;-;CrFe11MoNi3;Model of Ferrite");
    vproto.push_back("A3B_cP4_221_d_a;2;4;221;4;1;cP4;a;D0_9;ReO3;alphaReO3");
    vproto.push_back("A6B_cP7_221_f_a;2;7;221;4;2;cP7;a,x2;D2_1;CaB6;CaB6");
    vproto.push_back("A3B_cP8_223_c_a;2;8;223;4;1;cP8;a;A15;Cr3Si;Cr3Si");
    vproto.push_back("A_cP46_223_dik;1;46;223;3;4;cP46;a,x2,y3,z3;-;Si;Si46 Clathrate");
    vproto.push_back("A2B_cP6_224_b_a;2;6;224;4;1;cP6;a;C3;Cu2O;Cuprite");
    vproto.push_back("A7B_cF32_225_bd_a;2;8;225;4;1;cF32;a;-;Ca7Ge;Ca7Ge");
    vproto.push_back("AB3_cF16_225_a_bc;2;4;225;4;1;cF16;a;D0_3;BiF3;BiF3");
    vproto.push_back("A9B16C7_cF128_225_acd_2f_be;3;32;225;5;4;cF128;a,x5,x6,x7;-;Cr9Fe16Ni7;Model of Ferrite");
    vproto.push_back("A12B_cF52_225_i_a;2;13;225;4;2;cF52;a,y2;D2_f;UB12;UB12");
    vproto.push_back("AB2_cF12_225_a_c;2;3;225;4;1;cF12;a;C1;CaF2;Cu2Mg Cubic Laves");
    vproto.push_back("A6B23_cF116_225_e_acfh;2;29;225;4;4;cF116;a,x3,x4,y5;D8_4;Cr23C6;Cr23C6");
    vproto.push_back("AB2C_cF16_225_a_c_b;3;4;225;5;1;cF16;a;L2_1;AlCu2Mn;Heusler");
    vproto.push_back("A_cF4_225_a;1;1;225;3;1;cF4;a;A1;Cu;Face-Centered Cubic");
    vproto.push_back("AB18C8_cF108_225_a_eh_f;3;27;225;5;4;cF108;a,x2,x3,y4;-;CrFe18Ni8;Model of Austenite");
    vproto.push_back("AB_cF8_225_a_b;2;2;225;4;1;cF8;a;B1;NaCl;Rock Salt");
    vproto.push_back("A2B_cF24_227_c_a;2;6;227;4;1;cF24;a;C9;SiO2;Ideal beta-Cristobalite");
    vproto.push_back("AB2_cF96_227_e_cf;2;24;227;4;3;cF96;a,x2,x3;-;NiTi2;NiTi2");
    vproto.push_back("AB_cF16_227_a_b;2;4;227;4;1;cF16;a;B32;NaTl;NaTl");
    vproto.push_back("A_cF136_227_aeg;1;34;227;3;4;cF136;a,x2,x3,z3;-;Si;Si34 Clathrate");
    vproto.push_back("A2B_cF24_227_d_a;2;6;227;4;1;cF24;a;C15;Cu2Mg;Cu2Mg Cubic Laves");
    vproto.push_back("A_cF8_227_a;1;2;227;3;1;cF8;a;A4;C;Diamond");
    vproto.push_back("A2BC4_cF56_227_d_a_e;3;14;227;5;2;cF56;a,x3;H1_1;Al2MgO4;Spinel");
    vproto.push_back("AB2_cF48_227_c_e;2;12;227;4;2;cF48;a,x2;-;CTi2;CTi2");
    vproto.push_back("AB3C3_cF112_227_c_de_f;3;28;227;5;3;cF112;a,x3,x4;E9_3;Fe3W3C;Fe3W3C");
    vproto.push_back("A_cI2_229_a;1;1;229;3;1;cI2;a;A2;W;Body-Centered Cubic");
    vproto.push_back("A3B_cI8_229_b_a;2;4;229;4;1;cI8;a;-;H3S;High Presssure H3S");
    vproto.push_back("A4B3_cI14_229_c_b;2;7;229;4;1;cI14;a;-;Pt3O4;Pt3O4");
    vproto.push_back("A2B7_cI54_229_e_afh;2;27;229;4;4;cI54;a,x2,x3,y4;L2_2;Sb2Tl7;L22");
    vproto.push_back("AB12C3_cI32_229_a_h_b;3;16;229;5;2;cI32;a,y3;-;CrFe12Ni3;Model of Austenite");
    vproto.push_back("AB4C3_cI16_229_a_c_b;3;8;229;5;1;cI16;a;-;CrFe4Ni3;Model of Ferrite");
    vproto.push_back("A4B3_cI112_230_af_g;2;56;230;4;3;cI112;a,x2,y3;-;Ga4Ni3;Ga4Ni3");
    // -------------------------------------------------------------------------
    // Part 2
    // -------------------------------------------------------------------------
    vproto.push_back("A2B_aP6_2_aei_i;2;6;2;4;12;aP6;a,b/a,c/a,alpha,beta,gamma,x3,y3,z3,x4,y4,z4;-;H2S;H2S");
    vproto.push_back("A8B5_mP13_6_a7b_3a2b;2;13;6;4;30;mP13;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,z13;-;Mo8P5;Mo8P5");
    vproto.push_back("AB_mP4_6_2b_2a;2;4;6;4;12;mP4;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;FeNi;FeNi");
    vproto.push_back("A2B_mP12_7_4a_2a;2;12;7;4;22;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;H2S;H2S IV");
    vproto.push_back("A2B_mP18_7_6a_3a;2;18;7;4;31;mP18;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;As2Ba;As2Ba");
    vproto.push_back("A3B_mP16_7_6a_2a;2;16;7;4;28;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;epsilon-WO3;epsilon-WO3");
    vproto.push_back("A9B2_mP22_7_9a_2a;2;22;7;4;37;mP22;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-;Rh2Ga9;Rh2Ga9");
    vproto.push_back("A5B3_mC32_9_5a_3a;2;16;9;4;28;mC32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;alpha-P3N5;alpha-P3N5");
    vproto.push_back("AB3_mC16_9_a_3a;2;8;9;4;16;mC16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;H3Cl;H3Cl");
    vproto.push_back("A2B_mP6_10_mn_bg;2;6;10;4;8;mP6;a,b/a,c/a,beta,x3,z3,x4,z4;-;delta-PdCl2;delta-PdCl2");
    vproto.push_back("AB3_mP16_10_mn_3m3n;2;16;10;4;20;mP16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;H3Cl;H3Cl");
    vproto.push_back("ABC2_mP8_10_ac_eh_mn;3;8;10;5;8;mP8;a,b/a,c/a,beta,x5,z5,x6,z6;-;AuAgTe2;Muthmannite");
    vproto.push_back("AB_mP6_10_en_am;2;6;10;4;8;mP6;a,b/a,c/a,beta,x3,z3,x4,z4;-;LiSn;LiSn");
    vproto.push_back("A_mP8_10_2m2n;1;8;10;3;12;mP8;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;C;S-carbon");
    vproto.push_back("A7B2C2_mC22_12_aij_h_i;3;11;12;5;12;mC22;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,y5,z5;S2_{1};[Sc,Y]2Si2O7;Thortveitite");
    vproto.push_back("A_mC16_12_4i;1;8;12;3;12;mC16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;C;M-carbon");
    vproto.push_back("A2B_mP12_13_2g_ef;2;12;13;4;12;mP12;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4;-;H2S;H2S");
    vproto.push_back("A2B_mP6_14_e_a;2;6;14;4;7;mP6;a,b/a,c/a,beta,x2,y2,z2;-;gamma-PdCl2;gamma-PdCl2");
    vproto.push_back("A7B8_mP120_14_14e_16e;2;120;14;4;94;mP120;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30;-;alpha-C7H8;alpha-Toluene");
    vproto.push_back("AB3_mC16_15_e_cf;2;8;15;4;8;mC16;a,b/a,c/a,beta,y2,x3,y3,z3;-;H3Cl;H3Cl");
    vproto.push_back("A_mC24_15_2e2f;1;12;15;3;12;mC24;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4;-;H;H-III");
    vproto.push_back("A2B_oP12_17_abe_e;2;12;17;4;11;oP12;a,b/a,c/a,x1,x2,x3,y3,z3,x4,y4,z4;-;alpha-Ag2Se;alpha-Naumannite");
    vproto.push_back("AB3_oP16_19_a_3a;2;16;19;4;15;oP16;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;H3Cl;H3Cl");
    vproto.push_back("AB2_oC6_21_a_k;2;3;21;4;4;oC6;a,b/a,c/a,z2;-;Ta2H;Ta2H");
    vproto.push_back("A2BC2_oF40_22_fi_ad_gh;3;10;22;5;7;oF40;a,b/a,c/a,y3,z4,z5,y6;-;CeRu2B2;CeRu2B2");
    vproto.push_back("AB_oF8_22_a_c;2;2;22;4;3;oF8;a,b/a,c/a;-;FeS;FeS");
    vproto.push_back("A3B_oI32_23_ij2k_k;2;16;23;4;14;oI32;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;H3S;H3S");
    vproto.push_back("A8B2C12D2E_oI50_23_bcfk_i_3k_j_a;5;25;23;7;18;oI50;a,b/a,c/a,x4,z5,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Cu8(Fe,Zn)3Sn2S12;Stannoidite");
    vproto.push_back("ABC2_oI16_23_ab_i_k;3;8;23;5;7;oI16;a,b/a,c/a,z3,x4,y4,z4;-;NaFeS2;NaFeS2");
    vproto.push_back("ABC4_oI12_23_a_b_k;3;6;23;5;6;oI12;a,b/a,c/a,x3,y3,z3;-;BPS4;BPS4");
    vproto.push_back("AB7CD2_oI44_24_a_b3d_c_ac;4;22;24;6;17;oI44;a,b/a,c/a,x1,x2,y3,z4,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Na2MgAlF7;Weberite");
    vproto.push_back("A2B_oP12_26_abc_ab;2;12;26;4;14;oP12;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5;-;H2S;H2S");
    vproto.push_back("A2B_oP12_26_abc_ab;2;12;26;4;14;oP12;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5;-;beta-SeO2;beta-SeO2");
    vproto.push_back("A5B_oP24_26_3a3b2c_ab;2;24;26;4;25;oP24;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5,y6,z6,y7,z7,y8,z8,x9,y9,z9,x10,y10,z10;-;TlP5;TlP5");
    vproto.push_back("A6B4C16D_oP108_27_abcd4e_4e_16e_e;4;108;27;6;82;oP108;a,b/a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29;-;Ca4Al6O16S;Ca4Al6O16S");
    vproto.push_back("A2B_oP12_29_2a_a;2;12;29;4;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;ZrO2;ZrO2");
    vproto.push_back("AB2_oP12_29_a_2a;2;12;29;4;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;FeS2;Pyrite");
    vproto.push_back("ABC_oP12_29_a_a_a;3;12;29;5;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;CoAsS;Cobaltite");
    vproto.push_back("A5B3C15_oP46_30_a2c_bc_a7c;3;46;30;5;36;oP46;a,b/a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;-;Bi5Nb3O15;Bi5Nb3O15");
    vproto.push_back("ABC3_oP20_30_2a_c_3c;3;20;30;5;17;oP20;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;CuBrSe3;CuBrSe3");
    vproto.push_back("A13B2C2_oP34_32_a6c_c_c;3;34;32;5;28;oP34;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Re2O5[SO4]2;Re2O52");
    vproto.push_back("A2B3_oP40_33_4a_6a;2;40;33;4;33;oP40;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;kappa-Al2O3;kappa-alumina");
    vproto.push_back("A2B8C_oP22_34_c_4c_a;3;22;34;5;19;oP22;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;TiAl2Br8;TiAl2Br8");
    vproto.push_back("AB2_oP6_34_a_c;2;6;34;4;7;oP6;a,b/a,c/a,z1,x2,y2,z2;-;FeSb2;FeSb2");
    vproto.push_back("AB8C2_oC22_35_a_ab3e_e;3;11;35;5;14;oC22;a,b/a,c/a,z1,z2,z3,y4,z4,y5,z5,y6,z6,y7,z7;-;V2MoO8;V2MoO8");
    vproto.push_back("AB_oC8_36_a_a;2;4;36;4;7;oC8;a,b/a,c/a,y1,z1,y2,z2;-;HCl;HCl");
    vproto.push_back("A2B5C2_oC36_37_d_c2d_d;3;18;37;5;16;oC36;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Li2Si2O5;Li2Si2O5");
    vproto.push_back("A2B3_oC40_39_2d_2c2d;2;20;39;4;19;oC40;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Ta3S2;Ta3S2");
    vproto.push_back("A9BC_oC44_39_3c3d_a_c;3;22;39;5;21;oC44;a,b/a,c/a,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;VPCl9;VPCl9");
    vproto.push_back("AB2C_oC16_40_a_2b_b;3;8;40;5;10;oC16;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4;-;K2CdPb;K2CdPb");
    vproto.push_back("AB3_oC16_40_b_3b;2;8;40;4;11;oC16;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4;-;CeTe3;CeTe3");
    vproto.push_back("A10B3_oF52_42_2abce_ab;2;13;42;4;13;oF52;a,b/a,c/a,z1,z2,z3,z4,z5,y6,z6,x7,y7,z7;-;W3O10;W3O10");
    vproto.push_back("AB_oF8_42_a_a;2;2;42;4;5;oF8;a,b/a,c/a,z1,z2;-;BN;BN");
    vproto.push_back("A2BC2_oI20_45_c_b_c;3;10;45;5;10;oI20;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;-;MnGa2Sb2;MnGa2Sb2");
    vproto.push_back("ABC_oI36_46_ac_bc_3b;3;18;46;5;18;oI36;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4,y5,z5,x6,y6,z6,x7,y7,z7;-;TiFeSi;TiFeSi");
    vproto.push_back("A2B8CD_oP24_48_k_2m_d_b;4;24;48;6;10;oP24;a,b/a,c/a,z3,x4,y4,z4,x5,y5,z5;-;alpha-RbPr[MoO4]2;alpha-RbPr2");
    vproto.push_back("A5B2_oP14_49_dehq_ab;2;14;49;4;5;oP14;a,b/a,c/a,x6,y6;-;beta-Ta2O5;beta-Ta2O5");
    vproto.push_back("AB2C8D_oP24_49_g_q_2qr_e;4;24;49;6;12;oP24;a,b/a,c/a,x3,y3,x4,y4,x5,y5,x6,y6,z6;-;CsPr[MoO4]2;CsPr2");
    vproto.push_back("A2BC4_oP28_50_ij_ac_ijm;3;28;50;5;10;oP28;a,b/a,c/a,y3,y4,y5,y6,x7,y7,z7;-;La2NiO4;La2NiO4");
    vproto.push_back("A3BC2_oP48_50_3m_m_2m;3;48;50;5;21;oP48;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;alpha-Tl2TeO3;alpha-Tl2TeO3");
    vproto.push_back("A2B_oP24_52_2e_cd;2;24;52;4;11;oP24;a,b/a,c/a,z1,x2,x3,y3,z3,x4,y4,z4;-;GaCl2;GaCl2");
    vproto.push_back("A3B2_oP20_52_de_cd;2;20;52;4;9;oP20;a,b/a,c/a,z1,x2,x3,x4,y4,z4;-;Sr2Bi3;Sr2Bi3");
    vproto.push_back("ABC2_oP16_53_h_e_gh;3;16;53;5;9;oP16;a,b/a,c/a,x1,y2,y3,z3,y4,z4;-;TaNiTe2;TaNiTe2");
    vproto.push_back("ABC3_oP20_53_e_g_hi;3;20;53;5;10;oP20;a,b/a,c/a,x1,y2,y3,z3,x4,y4,z4;-;CuBrSe3;CuBrSe3");
    vproto.push_back("ABC3_oP20_54_e_d_cf;3;20;54;5;9;oP20;a,b/a,c/a,y1,z2,z3,x4,y4,z4;-;BiGaO3;BiGaO3");
    vproto.push_back("A2B_oP24_55_2g2h_gh;2;24;55;4;15;oP24;a,b/a,c/a,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6;-;GeAs2;GeAs2");
    vproto.push_back("A3B5_oP16_55_ch_agh;2;16;55;4;9;oP16;a,b/a,c/a,x3,y3,x4,y4,x5,y5;-;Rh5Ge3;Rh5Ge3");
    vproto.push_back("A_oP16_55_2g2h;1;16;55;3;11;oP16;a,b/a,c/a,x1,y1,x2,y2,x3,y3,x4,y4;-;C;R-carbon");
    vproto.push_back("A2B_oP6_58_g_a;2;6;58;4;5;oP6;a,b/a,c/a,x2,y2;C50/-/-;alpha-PdCl2/N2Os1/N2Pt1;alpha-PdCl2/N2Os1 (ICSD #157283)/N2Pt1 (ICSD #166463)"); //DX20201019 - added two metal nitrides
    vproto.push_back("ABC_oP6_59_a_b_a;3;6;59;5;6;oP6;a,b/a,c/a,z1,z2,z3;E0_{5};FeOCl;FeOCl"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A2B3_oP20_60_d_cd;2;20;60;4;10;oP20;a,b/a,c/a,y1,x2,y2,z2,x3,y3,z3;-;Rh2S3;Rh2S3");
    vproto.push_back("A3B_oP32_60_3d_d;2;32;60;4;15;oP32;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;WO3;WO3");
    vproto.push_back("A7B8_oP120_60_7d_8d;2;120;60;4;48;oP120;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;beta-C7H8;beta-Toluene");
    vproto.push_back("AB_oP48_61_3c_3c;2;48;61;4;21;oP48;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Benzene;Benzene");
    vproto.push_back("A2B3_oP20_62_2c_3c;2;20;62;4;13;oP20;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;D5_{10};Cr3C2;Tongbaite");
    vproto.push_back("A2B4C_oP28_62_ac_2cd_c;3;28;62;5;14;oP28;a,b/a,c/a,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6;S1_{2}/-;Mg2SiO4/Ag2CrO4;Forsterite/Ag2CrO4 (ICSD #16298)"); //DX20200703 - added R. Friedrich ternary oxide
    //DX 20180619 - added info to label in part 1 vproto.push_back("A2B_oP12_62_2c_c;2;12;62;4;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;C29;SrH2;SrH2");
    vproto.push_back("A3B_oP16_62_cd_c;2;16;62;4;10;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3;D0_{20};epsilon-NiAl3;epsilon-NiAl3");
    vproto.push_back("AB2C3_oP24_62_c_d_cd;3;24;62;5;13;oP24;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3,x4,y4,z4;E9_{e};CuFe2S3;Cubanite");
    vproto.push_back("AB3_oP16_62_c_3c;2;16;62;4;11;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4;D0_{8};MoO3;Molybdite");
    vproto.push_back("AB4C_oP24_62_c_2cd_c;3;24;62;5;14;oP24;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5;H0_{2};BaSO4;Barite");
    //DX 20180619 - added info to label in part 1 vproto.push_back("AB_oP8_62_c_c;2;8;62;4;7;oP8;a,b/a,c/a,x1,z1,x2,z2;B14;FeAs;Westerveldite");
    vproto.push_back("A2BC3_oC24_63_e_c_cg;3;12;63;5;8;oC24;a,b/a,c/a,y1,y2,x3,x4,y4;-;KFe2S3;Rasvumite");
    vproto.push_back("A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h;3;130;63;5;59;oC260;a,b/a,c/a,y1,y2,y3,x4,y5,z5,y6,z6,y7,z7,y8,z8,y9,z9,y10,z10,y11,z11,y12,z12,y13,z13,y14,z14,y15,z15,y16,z16,x17,y17,x18,y18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26;-;La43Ni17Mg5;La43Ni17Mg5");
    vproto.push_back("A6B_oC28_63_efg_c;2;14;63;4;9;oC28;a,b/a,c/a,y1,x2,y3,z3,x4,y4;D2_{h};MnAl6;MnAl6");
    vproto.push_back("AB3C_oC20_63_a_cf_c;3;10;63;5;7;oC20;a,b/a,c/a,y2,y3,y4,z4;-;MgSiO3;Post-perovskite");
    vproto.push_back("AB4C_oC24_63_a_fg_c;3;12;63;5;8;oC24;a,b/a,c/a,y2,y3,z3,x4,y4;-;MgO4S;MgSO4");
    vproto.push_back("AB4C_oC24_63_c_fg_c;3;12;63;5;9;oC24;a,b/a,c/a,y1,y2,y3,z3,x4,y4;H0_{1};CaSO4;Anhydrite");
    vproto.push_back("A2B_oC24_64_2f_f;2;12;64;4;9;oC24;a,b/a,c/a,y1,z1,y2,z2,y3,z3;-;H2S;H2S");
    vproto.push_back("A2B4C_oC28_66_l_kl_a;3;14;66;5;8;oC28;a,b/a,c/a,z2,x3,y3,x4,y4;-;SrAl2Se4;SrAl2Se4");
    vproto.push_back("A3B_oC64_66_gi2lm_2l;2;32;66;4;16;oC64;a,b/a,c/a,x1,z2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,z7;-;H3S;H3S");
    vproto.push_back("A3B_oC64_66_kl2m_bdl;2;32;66;4;14;oC64;a,b/a,c/a,z3,x4,y4,x5,y5,x6,y6,z6,x7,y7,z7;-;beta-ThI3;beta-ThI3");
    vproto.push_back("A2BC_oC16_67_ag_b_g;3;8;67;5;5;oC16;a,b/a,c/a,z3,z4;-;Al2CuIr;Al2CuIr");
    vproto.push_back("ABC2_oC16_67_b_g_ag;3;8;67;5;5;oC16;a,b/a,c/a,z3,z4;-;HoCuP2;HoCuP2");
    vproto.push_back("AB_oC8_67_a_g;2;4;67;4;4;oC8;a,b/a,c/a,z2;-;alpha-FeSe;alpha-FeSe");
    vproto.push_back("AB_oC8_67_a_g;2;4;67;4;4;oC8;a,b/a,c/a,z2;-;alpha-PbO;alpha-PbO");
    vproto.push_back("AB4_oC20_68_a_i;2;10;68;4;6;oC20;a,b/a,c/a,x2,y2,z2;-;PdSn4;PdSn4");
    vproto.push_back("AB2_oF48_70_f_fg;2;12;70;4;6;oF48;a,b/a,c/a,y1,y2,z3;D1_{f};Mn2B;Mn2B");
    vproto.push_back("A4B3_oI14_71_gh_cg;2;7;71;4;6;oI14;a,b/a,c/a,y2,y3,y4;D7_{b};Ta3B4;Ta3B4");
    vproto.push_back("ABC_oI12_71_h_j_g;3;6;71;5;6;oI12;a,b/a,c/a,y1,y2,z3;-;NbPS;NbPS");
    vproto.push_back("ABCD3_oI48_73_d_e_e_ef;4;24;73;6;10;oI48;a,b/a,c/a,y1,z2,z3,z4,x5,y5,z5;-;KAg[CO3];KAg");
    vproto.push_back("A2B_oI12_74_h_e;2;6;74;4;6;oI12;a,b/a,c/a,z1,y2,z2;-;KHg2;KHg2");
    vproto.push_back("A4B_oI20_74_beh_e;2;10;74;4;7;oI20;a,b/a,c/a,z2,z3,y4,z4;D1_{b};Al4U;Al4U");
    vproto.push_back("AB2C12D4_tP76_75_2a2b_2d_12d_4d;4;76;75;6;60;tP76;a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22;-;BaCr2Ru4O12;BaCr2Ru4O12");
    vproto.push_back("A2BC_tP16_76_2a_a_a;3;16;76;5;14;tP16;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;LaRhC2;LaRhC2");
    vproto.push_back("A3B7_tP40_76_3a_7a;2;40;76;4;32;tP40;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Cs3P7;Cs3P7");
    vproto.push_back("A2B6CD7_tP64_77_2d_6d_d_ab6d;4;64;77;6;49;tP64;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17;-;MgB2O(OH)6;Pinnoite");
    vproto.push_back("A2B_tP48_77_8d_4d;2;48;77;4;38;tP48;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;H2S;H2S III");
    vproto.push_back("A2B7C2_tP88_78_4a_14a_4a;3;88;78;5;68;tP88;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22;-;Sr2As2O7;Sr2As2O7");
    vproto.push_back("A2BC2_tI20_79_c_2a_c;3;10;79;5;10;tI20;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4;-;TlZn2Sb2;TlZn2Sb2");
    vproto.push_back("AB2_tI48_80_2b_4b;2;24;80;4;20;tI48;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;beta-NbO2;beta-NbO2");
    vproto.push_back("AB2_tP12_81_adg_2h;2;12;81;4;9;tP12;a,c/a,z3,x4,y4,z4,x5,y5,z5;-;GeSe2;GeSe2");
    vproto.push_back("A3B_tI32_82_3g_g;2;16;82;4;14;tI32;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;D0_{e};Ni3P;Ni3P");
    vproto.push_back("A3B2_tP10_83_adk_j;2;10;83;4;6;tP10;a,c/a,x3,y3,x4,y4;-;Ti2Ge3;Ti2Ge3");
    vproto.push_back("A2B_tP30_85_ab2g_cg;2;30;85;4;12;tP30;a,c/a,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;SrBr2;SrBr2");
    vproto.push_back("AB3_tP32_86_g_3g;2;32;86;4;14;tP32;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Ti3P;Ti3P");
    vproto.push_back("A4B_tI20_88_f_a;2;10;88;4;5;tI20;a,c/a,x2,y2,z2;-;ThCl4;ThCl4");
    vproto.push_back("AB2_tI96_88_2f_4f;2;48;88;4;20;tI96;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;alpha-NbO2;alpha-NbO2");
    vproto.push_back("A17BC4D_tP184_89_17p_p_4p_io;4;184;89;6;70;tP184;a,c/a,z1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24;-;C17FeO4Pt;C17FeO4Pt");
    vproto.push_back("A4B2C13D_tP40_90_g_d_cef2g_c;4;40;90;6;16;tP40;a,c/a,z1,z2,z3,x4,x5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Na4Ti2Si8O22[H2O]4;Na4Ti2Si8O224");
    vproto.push_back("AB4C17D4E_tP54_90_a_g_c4g_g_c;5;54;90;7;22;tP54;a,c/a,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;BaCu4[VO][PO4]4;BaCu44");
    vproto.push_back("ABC_tP24_91_d_d_d;3;24;91;5;11;tP24;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;ThBC;ThBC");
    vproto.push_back("AB32CD4E8_tP184_93_i_16p_af_2p_4p;5;184;93;7;69;tP184;a,c/a,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25;-;AsPh4CeS8P4Me8;AsPh4CeS8P4Me8");
    vproto.push_back("A14B3C5_tP44_94_c3g_ad_bg;3;44;94;5;16;tP44;a,c/a,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Na5Fe3F14;Na5Fe3F14");
    vproto.push_back("A6B2C_tP18_94_eg_c_a;3;18;94;5;7;tP18;a,c/a,z2,x3,x4,y4,z4;-;Li2MoF6;Li2MoF6");
    vproto.push_back("ABC_tP24_95_d_d_d;3;24;95;5;11;tP24;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;ThBC;ThBC");
    vproto.push_back("A2B8CD_tI24_97_d_k_a_b;4;12;97;6;5;tI24;a,c/a,x4,y4,z4;-;NaGdCu2F8;NaGdCu2F8");
    vproto.push_back("AB8C2_tI44_97_e_2k_cd;3;22;97;5;9;tI44;a,c/a,z3,x4,y4,z4,x5,y5,z5;-;Ta2Se8I;Ta2Se8I");
    vproto.push_back("A2B_tI12_98_f_a;2;6;98;4;3;tI12;a,c/a,x2;-;CdAs2;CdAs2");
    vproto.push_back("A2B8C2D_tP26_100_c_abcd_c_a;4;26;100;6;14;tP26;a,c/a,z1,z2,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;-;Ba2TiSi2O8;Fresnoite");
    vproto.push_back("A3B11C6_tP40_100_ac_bc2d_cd;3;40;100;5;19;tP40;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Ce3Si6N11;Ce3Si6N11");
    vproto.push_back("A7B7C2_tP32_101_bde_ade_d;3;32;101;5;16;tP32;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7;-;gamma-MgNiSn;gamma-MgNiSn");
    vproto.push_back("A2B3_tP20_102_2c_b2c;2;20;102;4;11;tP20;a,c/a,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;Gd3Al2;Gd3Al2");
    vproto.push_back("AB4_tP10_103_a_d;2;10;103;4;6;tP10;a,c/a,z1,x2,y2,z2;-;NbTe4;NbTe4");
    vproto.push_back("A5B5C4_tP28_104_ac_ac_c;3;28;104;5;13;tP28;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Ba5In4Bi5;Ba5In4Bi5");
    vproto.push_back("AB6C4_tP22_104_a_2ac_c;3;22;104;5;11;tP22;a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5;-;Tl4HgI6;Tl4HgI6");
    vproto.push_back("A2BC2_tP20_105_f_ac_2e;3;20;105;5;11;tP20;a,c/a,z1,z2,x3,z3,x4,z4,x5,y5,z5;-;BaGe2As2;BaGe2As2");
    vproto.push_back("A3BC3D_tP64_106_3c_c_3c_c;4;64;106;6;26;tP64;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;NaZn[OH]3;NaZn3");
    vproto.push_back("A5B7_tI24_107_ac_abd;2;12;107;4;9;tI24;a,c/a,z1,z2,z3,x4,z4,x5,z5;-;Co5Ge7;Co5Ge7");
    vproto.push_back("AB_tI4_107_a_a;2;2;107;4;4;tI4;a,c/a,z1,z2;-;GeP;GeP");
    vproto.push_back("A3B5_tI32_108_ac_a2c;2;16;108;4;10;tI32;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5;-;Sr5Si3;Sr5Si3");
    vproto.push_back("ABC_tI12_109_a_a_a;3;6;109;5;5;tI12;a,c/a,z1,z2,z3;-;LaPtSi;LaPtSi");
    vproto.push_back("AB_tI8_109_a_a;2;4;109;4;4;tI8;a,c/a,z1,z2;-;NbAs;NbAs");
    vproto.push_back("A2BC8_tI176_110_2b_b_8b;3;88;110;5;35;tI176;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-;Be[BH4]2;Be2");
    vproto.push_back("A2B_tP12_111_2n_adf;2;12;111;4;6;tP12;a,c/a,x4,z4,x5,z5;-;MnF2;MnF2");
    vproto.push_back("AB_tP8_111_n_n;2;8;111;4;6;tP8;a,c/a,x1,z1,x2,z2;-;VN;VN"); //DX 20180925 - prototype name should be VN not NV
    vproto.push_back("AB4C_tP12_112_b_n_e;3;12;112;5;5;tP12;a,c/a,x3,y3,z3;-;alpha-CuAlCl4;alpha-CuAlCl4");
    vproto.push_back("A2BC7D2_tP24_113_e_a_cef_e;4;24;113;6;12;tP24;a,c/a,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6;S5_{3};Ca2MgSi2O7;Akermanite");
    vproto.push_back("A3B_tP32_114_3e_e;2;32;114;4;14;tP32;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;SeO3;SeO3");
    vproto.push_back("A4B_tP10_114_e_a;2;10;114;4;5;tP10;a,c/a,x2,y2,z2;-;Pd4Se;Pd4Se");
    vproto.push_back("A2B3_tP5_115_g_ag;2;5;115;4;4;tP5;a,c/a,z2,z3;-;Rh3P2;Rh3P2");
    vproto.push_back("AB2_tP12_115_j_egi;2;12;115;4;7;tP12;a,c/a,z1,z2,x3,x4,z4;-;HgI2;HgI2");
    vproto.push_back("A2B3_tP20_116_bci_fj;2;20;116;4;7;tP20;a,c/a,x3,z4,x5,y5,z5;-;Ru2Sn3;Ru2Sn3");
    vproto.push_back("A2B3_tP20_117_i_adgh;2;20;117;4;7;tP20;a,c/a,x3,x4,x5,y5,z5;D5_{12};beta-Bi2O3;beta-Bi2O3"); //DX 20180925 - added Strukturbericht designation
    vproto.push_back("A3B_tP16_118_ei_f;2;16;118;4;7;tP16;a,c/a,z1,x2,x3,y3,z3;-;RuIn3;RuIn3");
    vproto.push_back("A5B3_tP32_118_g2i_aceh;2;32;118;4;11;tP32;a,c/a,z3,x4,z5,x6,y6,z6,x7,y7,z7;-;Ir3Ga5;Ir3Ga5");
    vproto.push_back("A3B_tI24_119_b2i_af;2;12;119;4;7;tI24;a,c/a,z3,x4,z4,x5,z5;-;RbGa3;RbGa3");
    vproto.push_back("AB_tI4_119_c_a;2;2;119;4;2;tI4;a,c/a;-;GaSb;GaSb");
    vproto.push_back("A4BC2_tI28_120_i_d_e;3;14;120;5;6;tI28;a,c/a,x2,x3,y3,z3;-;KAu4Sn2;KAu4Sn2");
    vproto.push_back("A4BC4D_tP10_123_gh_a_i_d;4;10;123;6;5;tP10;a,c/a,z3,z4,z5;-;CaRbFe4As4;CaRbFe4As4"); //DX 20180925 - prototype name should be Ca not Cs
    vproto.push_back("AB4C_tP12_124_a_m_c;3;12;124;5;4;tP12;a,c/a,x3,y3;-;Nb4CoSi;Nb4CoSi");
    vproto.push_back("AB4_tP10_124_a_m;2;10;124;4;4;tP10;a,c/a,x2,y2;-;NbTe4;NbTe4");
    vproto.push_back("A4B_tP10_125_m_a;2;10;125;4;4;tP10;a,c/a,x2,z2;-;PtPb4;PtPb4");
    vproto.push_back("ABC4_tP12_125_a_b_m;3;12;125;5;4;tP12;a,c/a,x3,z3;-;KCeSe4;KCeSe4");
    vproto.push_back("A2BC4_tP28_126_cd_e_k;3;28;126;5;6;tP28;a,c/a,z3,x4,y4,z4;-;BiAl2S4;BiAl2S4");
    vproto.push_back("A4B_tP20_127_ehj_g;2;20;127;4;7;tP20;a,c/a,z1,x2,x3,x4,y4;D1_{e};ThB4;ThB4");
    vproto.push_back("A6B2C_tP18_128_eh_d_a;3;18;128;5;5;tP18;a,c/a,z3,x4,y4;-;K2SnCl6;K2SnCl6");
    vproto.push_back("A7B2C_tP40_128_egi_h_e;3;40;128;5;10;tP40;a,c/a,z1,z2,x3,x4,y4,x5,y5,z5;E9_{a};FeCu2Al7;FeCu2Al7");
    vproto.push_back("A2BC4_tP28_130_f_c_g;3;28;130;5;7;tP28;a,c/a,z1,x2,x3,y3,z3;-;CuBi2O4;CuBi2O4");
    vproto.push_back("A5B3_tP32_130_cg_cf;2;32;130;4;8;tP32;a,c/a,z1,z2,x3,x4,y4,z4;-;Ba5Si3;Ba5Si3");
    vproto.push_back("A2B2C4D_tP18_132_e_i_o_d;4;18;132;6;5;tP18;a,c/a,x3,x4,z4;-;Rb2TiCu2Se4;Rb2TiCu2S4");
    vproto.push_back("AB6C_tP16_132_d_io_a;3;16;132;5;5;tP16;a,c/a,x3,x4,z4;-;AgUF6;AgUF6");
    vproto.push_back("AB3_tP32_133_h_i2j;2;32;133;4;6;tP32;a,c/a,x1,x2,x3,x4;-;beta-V3S;beta-V3S");
    vproto.push_back("A2B_tP24_135_gh_h;2;24;135;4;7;tP24;a,c/a,x1,x2,y2,x3,y3;C47;SeO2;Downeyite");
    vproto.push_back("A4B2C_tP28_135_gh_h_d;3;28;135;5;7;tP28;a,c/a,x2,x3,y3,x4,y4;-;ZnSb2O4;ZnSb2O4");
    vproto.push_back("A2B3_tP40_137_cdf_3g;2;40;137;4;11;tP40;a,c/a,z1,z2,x3,y4,z4,y5,z5,y6,z6;D5_{9};Zn3P2;Zn3P2");
    vproto.push_back("A2B_tP6_137_d_a;2;6;137;4;3;tP6;a,c/a,z2;-;ZrO2;ZrO2");
    vproto.push_back("A4BC4_tP18_137_g_b_g;3;18;137;5;6;tP18;a,c/a,y2,z2,y3,z3;-;CeCo4B4;CeCo4B4");
    vproto.push_back("AB2_tP6_137_a_d;2;6;137;4;3;tP6;a,c/a,z2;C13;HgI2;HgI2");
    vproto.push_back("A_tP12_138_bi;1;12;138;3;4;tP12;a,c/a,x2,z2;-;C;C");
    vproto.push_back("AB_tI8_139_e_e;2;4;139;4;4;tI8;a,c/a,z1,z2;D3_{1}/-;Hg2Cl2/-;Calomel/-");
    vproto.push_back("A3B5_tI32_140_ah_bk;2;16;140;4;5;tI32;a,c/a,x3,x4,y4;D8_{m};W5Si3;W5Si3");
    vproto.push_back("A3B5_tI32_140_ah_cl;2;16;140;4;5;tI32;a,c/a,x3,x4,z4;D8_{l};Cr5B3;Cr5B3");
    //DX 20180619 - added info to label in part 1 vproto.push_back("A2B_tI12_141_e_a;2;6;141;4;3;tI12;a,c/a,z2;C_{c};alpha-ThSi2;alpha-ThSi2");
    vproto.push_back("A_tI16_142_f;1;8;142;3;3;tI16;a,c/a,x1;-;S;S-III");
    vproto.push_back("A4B14C3_hP21_143_bd_ac4d_d;3;21;143;5;23;hP21;a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Ta3Al4O13[OH];Simpsonite");
    vproto.push_back("A4B6C_hP11_143_bd_2d_a;3;11;143;5;13;hP11;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;ScRh6P4;ScRh6P4");
    vproto.push_back("AB2_hP12_143_cd_ab2d;2;12;143;4;14;hP12;a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;MoS2;MoS2");
    vproto.push_back("A4B_hP15_144_4a_a;2;15;144;4;17;hP15;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;IrGe4;IrGe4");
    vproto.push_back("AB_hP6_144_a_a;2;6;144;4;8;hP6;a,c/a,x1,y1,z1,x2,y2,z2;-;ZnTe;ZnTe"); //DX 20180925 - prototype name should be ZnTe not TeZn
    vproto.push_back("A2B3C3DE7_hP48_145_2a_3a_3a_a_7a;5;48;145;7;50;hP48;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;NaCa3[CO3]2F3[H2O];Sheldrickite");
    vproto.push_back("A3BC_hR5_146_b_a_a;3;5;146;5;7;hR5;a,c/a,x1,x2,x3,y3,z3;-;gamma-Ag3SI;gamma-Ag3SI");
    vproto.push_back("ABC3_hR10_146_2a_2a_2b;3;10;146;5;12;hR10;a,c/a,x1,x2,x3,x4,x5,y5,z5,x6,y6,z6;-;FePSe3;FePSe3");
    vproto.push_back("A2B4C_hR42_148_2f_4f_f;3;42;148;5;23;hR42;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;S1_{3};Be2SiO4;Phenakite");
    vproto.push_back("A2B_hR18_148_2f_f;2;18;148;4;11;hR18;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;beta-PdCl2;beta-PdCl2");
    vproto.push_back("AB3_hP24_149_acgi_3l;2;24;149;4;13;hP24;a,c/a,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Ti3O;Ti3O");
    vproto.push_back("A3B_hP24_153_3c_2b;2;24;153;4;13;hP24;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;CrCl3;CrCl3");
    vproto.push_back("A_hP9_154_bc;1;9;154;3;6;hP9;a,c/a,x1,x2,y2,z2;-;S;S-II");
    vproto.push_back("AB2_hP9_156_b2c_3a2bc;2;9;156;4;11;hP9;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8,z9;-;CdI2;CdI2");
    vproto.push_back("AB_hP12_156_2ab3c_2ab3c;2;12;156;4;14;hP12;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12;-;CuI;CuI");
    vproto.push_back("AB_hP4_156_ab_ab;2;4;156;4;6;hP4;a,c/a,z1,z2,z3,z4;-;beta-CuI;beta-CuI");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("AB_hP4_156_ac_ac;2;4;156;4;6;hP4;a,c/a,z1,z2,z3,z4;-;beta-CuI;beta-CuI");
    vproto.push_back("A5B6C2_hP13_157_2ac_2c_b;3;13;157;5;11;hP13;a,c/a,z1,z2,z3,x4,z4,x5,z5,x6,z6;-;Ag5Pb2O6;Ag5Pb2O6");
    vproto.push_back("A3B_hP8_158_d_a;2;8;158;4;6;hP8;a,c/a,z1,x2,y2,z2;-;beta-RuCl3;beta-RuCl3");
    vproto.push_back("A2B3_hP20_159_bc_2c;2;20;159;4;12;hP20;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Bi2O3;Bi2O3");
    vproto.push_back("A4B3_hP28_159_ab2c_2c;2;28;159;4;16;hP28;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;alpha-Si3N4;Nierite");
    vproto.push_back("AB4C7D_hP26_159_b_ac_a2c_b;4;26;159;6;15;hP26;a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;YbBaCo4O7;YbBaCo4O7");
    vproto.push_back("A3B_hR4_160_b_a;2;4;160;4;5;hR4;a,c/a,x1,x2,z2;-;H3S;H3S");
    vproto.push_back("A8B5_hR26_160_a3bc_a3b;2;26;160;4;19;hR26;a,c/a,x1,x2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,y9,z9;D8_{10};Cr5Al8;Cr5Al8"); //DX 20180925 - prototype name should be Cr5Al8 not Al8Cr5
    vproto.push_back("ABC_hR3_160_a_a_a;3;3;160;5;5;hR3;a,c/a,x1,x2,x3;F0_{2};COS;Carbonyl Sulphide");
    vproto.push_back("AB_hR10_160_5a_5a;2;10;160;4;12;hR10;a,c/a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10;B7;SiC;Moissanite-15R");
    vproto.push_back("A2B3_hP5_164_d_ad;2;5;164;4;4;hP5;a,c/a,z2,z3;D5_{2};La2O3;La2O3");
    vproto.push_back("AB2_hP9_164_bd_c2d;2;9;164;4;6;hP9;a,c/a,z2,z3,z4,z5;-;deltaH^II-NW2;deltaH^II-NW2");
    vproto.push_back("ABC2_hP4_164_a_b_d;3;4;164;5;3;hP4;a,c/a,z3;-;CuNiSb2;CuNiSb2");
    vproto.push_back("A3B_hP24_165_bdg_f;2;24;165;4;7;hP24;a,c/a,z2,x3,x4,y4,z4;D0_{21};Cu3P;Cu3P");
    vproto.push_back("A4B3_hR7_166_2c_ac;2;7;166;4;5;hR7;a,c/a,x2,x3,x4;D7_{1};Al4C3;Al4C3");
    vproto.push_back("ABC_hR6_166_c_c_c;3;6;166;5;5;hR6;a,c/a,x1,x2,x3;-;SmSI;SmSI");
    vproto.push_back("AB3C_hR10_167_b_e_a;3;10;167;5;3;hR10;a,c/a,x3;-;PrNiO3;PrNiO3");
    vproto.push_back("ABC2_hR24_167_e_e_2e;3;24;167;5;6;hR24;a,c/a,x1,x2,x3,x4;F5_{13};KBO2;KBO2");
    vproto.push_back("A2B13C4_hP57_168_d_c6d_2d;3;57;168;5;30;hP57;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;K2Ta4O9F4;K2Ta4O9F4");
    vproto.push_back("AB4C_hP72_168_2d_8d_2d;3;72;168;5;38;hP72;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;Al[PO4];Al");
    vproto.push_back("A2B3_hP30_169_2a_3a;2;30;169;4;17;hP30;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;alpha-Al2S3;alpha-Al2S3");
    vproto.push_back("A2B3_hP30_170_2a_3a;2;30;170;4;17;hP30;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Al2S3;Al2S3");
    vproto.push_back("A10B2C_hP39_171_5c_c_a;3;39;171;5;21;hP39;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Sr[S2O6][H2O]4;Sr4");
    vproto.push_back("A10B2C_hP39_172_5c_c_a;3;39;172;5;21;hP39;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Sr[S2O6][H2O]4;Sr4");
    vproto.push_back("A3B_hP8_173_c_b;2;8;173;4;6;hP8;a,c/a,z1,x2,y2,z2;-;PI3;PI3");
    vproto.push_back("A4B3_hP14_173_bc_c;2;14;173;4;9;hP14;a,c/a,z1,x2,y2,z2,x3,y3,z3;-;beta-Si3N4;beta-Si3N4");
    vproto.push_back("A12B7C2_hP21_174_2j2k_ajk_cf;3;21;174;5;14;hP21;a,c/a,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9;-;Fe12Zr2P7;Fe12Zr2P7");
    vproto.push_back("ABC_hP12_174_cj_fk_aj;3;12;174;5;8;hP12;a,c/a,x4,y4,x5,y5,x6,y6;-;GdSI;GdSI");
    vproto.push_back("A8B7C6_hP21_175_ck_aj_k;3;21;175;5;8;hP21;a,c/a,x3,y3,x4,y4,x5,y5;-;Nb7Ru6B8;Nb7Ru6B8");
    vproto.push_back("ABC_hP36_175_jk_jk_jk;3;36;175;5;14;hP36;a,c/a,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6;-;Mg[NH];Mg");
    vproto.push_back("A3B2_hP10_176_h_bc;2;10;176;4;4;hP10;a,c/a,x3,y3;-;Er3Ru2;Er3Ru2");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3B2_hP10_176_h_bd;2;10;176;4;4;hP10;a,c/a,x3,y3;-;Er3Ru2;Er3Ru2");
    vproto.push_back("A3B3C_hP14_176_h_h_c;3;14;176;5;6;hP14;a,c/a,x2,y2,x3,y3;-;Fe3Te3Tl;Fe3Te3Tl");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3B3C_hP14_176_h_h_d;3;14;176;5;6;hP14;a,c/a,x2,y2,x3,y3;-;Fe3Te3Tl;Fe3Te3Tl");
    vproto.push_back("A3B_hP8_176_h_c;2;8;176;4;4;hP8;a,c/a,x2,y2;-;UCl3;UCl3");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3B_hP8_176_h_d;2;8;176;4;4;hP8;a,c/a,x2,y2;-;UCl3;UCl3");
    vproto.push_back("A2B_hP36_177_j2lm_n;2;36;177;4;9;hP36;a,c/a,x1,x2,x3,x4,x5,y5,z5;-;SiO2;SiO2");
    vproto.push_back("AB3_hP24_178_b_ac;2;24;178;4;7;hP24;a,c/a,x1,x2,x3,y3,z3;-;AuF3;AuF3");
    vproto.push_back("A_hP6_178_a;1;6;178;3;3;hP6;a,c/a,x1;-;Sc;Sc-V");
    vproto.push_back("AB3_hP24_179_b_ac;2;24;179;4;7;hP24;a,c/a,x1,x2,x3,y3,z3;-;AuF3;AuF3");
    vproto.push_back("A2B_hP9_181_j_c;2;9;181;4;3;hP9;a,c/a,x2;-;beta-SiO2;beta-SiO2");
    vproto.push_back("ABC_hP3_183_a_a_a;3;3;183;5;5;hP3;a,c/a,z1,z2,z3;-;AuCN;AuCN");
    vproto.push_back("AB_hP6_183_c_ab;2;6;183;4;5;hP6;a,c/a,z1,z2,z3;-;CrFe3NiSn5;CrFe3NiSn5");
    vproto.push_back("AB4C_hP72_184_d_4d_d;3;72;184;5;20;hP72;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Al[PO4];Al");
    vproto.push_back("A3BC_hP30_185_cd_c_ab;3;30;185;5;11;hP30;a,c/a,z1,z2,x3,z3,x4,z4,x5,y5,z5;-;KNiCl3;KNiCl3");
    vproto.push_back("A3B_hP24_185_ab2c_c;2;24;185;4;10;hP24;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5;-;Cu3P;Cu3P");
    vproto.push_back("A3B_hP8_185_c_a;2;8;185;4;5;hP8;a,c/a,z1,x2,z2;-;beta-RuCl3;beta-RuCl3");
    vproto.push_back("AB3_hP24_185_c_ab2c;2;24;185;4;10;hP24;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5;-;Na3As;Na3As");
    vproto.push_back("A3B7_hP20_186_c_b2c;2;20;186;4;9;hP20;a,c/a,z1,x2,z2,x3,z3,x4,z4;D10_{2};Fe3Th7;Fe3Th7");
    vproto.push_back("AB3_hP4_187_e_fh;2;4;187;4;3;hP4;a,c/a,z3;-;Re3N;Re3N");
    vproto.push_back("A3BC_hP10_188_k_c_a;3;10;188;5;4;hP10;a,c/a,x3,y3;-;LiScI3;LiScI3");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3BC_hP10_188_k_a_e;3;10;188;5;4;hP10;a,c/a,x3,y3;-;LiScI3;LiScI3");
    vproto.push_back("AB9C4_hP28_188_e_kl_ak;3;28;188;5;9;hP28;a,c/a,x3,y3,x4,y4,x5,y5,z5;S3_{2};BaSi4O9;BaSi4O9"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A8BC3D6_hP18_189_bfh_a_g_i;4;18;189;6;7;hP18;a,c/a,x3,x4,z5,x6,z6;E9_{b};pi-FeMg3Al8Si6;pi-FeMg3Al8Si6");
    vproto.push_back("A9BC3D5_hP18_189_fi_a_g_bh;4;18;189;6;7;hP18;a,c/a,x3,x4,z5,x6,z6;-;pi-FeMg3Al9Si5;pi-FeMg3Al9Si5");
    vproto.push_back("A2B_hP18_190_gh_bf;2;18;190;4;6;hP18;a,c/a,z2,x3,x4,y4;-;Li2Sb;Li2Sb");
    vproto.push_back("A5B3_hP16_190_bdh_g;2;16;190;4;5;hP16;a,c/a,x3,x4,y4;-;alpha-Sm3Ge5;alpha-Sm3Ge5");
    vproto.push_back("AB_hP24_190_i_afh;2;24;190;4;8;hP24;a,c/a,z2,x3,y3,x4,y4,z4;-;FeS;Troilite");
    vproto.push_back("A2B3C18D6_hP58_192_c_f_lm_l;4;58;192;6;9;hP58;a,c/a,x3,y3,x4,y4,x5,y5,z5;G3_{1};Be3Al2Si6O18;Beryl");
    vproto.push_back("AB2_hP72_192_m_j2kl;2;72;192;4;10;hP72;a,c/a,x1,x2,x3,x4,y4,x5,y5,z5;-;AlPO4;AlPO4");
    vproto.push_back("A5B3_hP16_193_dg_g;2;16;193;4;4;hP16;a,c/a,x2,x3;D8_{8};Mn5Si3;Mavlyanovite"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A3B_hP16_194_gh_ac;2;16;194;4;3;hP16;a,c/a,x4;D0_{24};Ni3Ti;Ni3Ti");
    vproto.push_back("A5B2_hP28_194_ahk_ch;2;28;194;4;6;hP28;a,c/a,x3,x4,x5,z5;D8_{11};Co2Al5;Co2Al5");
    vproto.push_back("A9B3C_hP26_194_hk_h_a;3;26;194;5;6;hP26;a,c/a,x2,x3,x4,z4;E9_{c};Al9Mn3Si;Al9Mn3Si");
    vproto.push_back("A12BC4_cP34_195_2j_ab_2e;3;34;195;5;9;cP34;a,x3,x4,x5,y5,z5,x6,y6,z6;-;PrRu4P12;PrRu4P12");
    vproto.push_back("A12B2C_cF60_196_h_bc_a;3;15;196;5;4;cF60;a,x4,y4,z4;-;Cu2Fe[CN]6;Cu2Fe6");
    //DX 20180925 - wrong space group, should be 210: vproto.push_back("A12B36CD12_cF488_196_2h_6h_ac_fgh;4;122;196;6;30;cF488;a,x3,x4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;-;MgB12H12[H2O]12;MgB12H1212");
    vproto.push_back("ABC3_cP20_198_a_a_b;3;20;198;5;6;cP20;a,x1,x2,x3,y3,z3;G3;NaClO3;Sodium Chlorate");
    vproto.push_back("A2B11_cP39_200_f_aghij;2;39;200;4;7;cP39;a,x2,x3,x4,x5,y6,z6;D8_{c};Mg2Zn11;Mg2Zn11");
    vproto.push_back("AB3C_cP60_201_be_fh_g;3;60;201;5;7;cP60;a,x2,x3,x4,x5,y5,z5;-;KSbO3;KSbO3");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("AB3C_cP60_201_ce_fh_g;3;60;201;5;7;cP60;a,x2,x3,x4,x5,y5,z5;-;KSbO3;KSbO3");
    vproto.push_back("A6B6C_cF104_202_h_h_c;3;26;202;5;5;cF104;a,y2,z2,y3,z3;-;KB6H6;KB6H6");
    vproto.push_back("A_cF240_202_h2i;1;60;202;3;9;cF240;a,y1,z1,x2,y2,z2,x3,y3,z3;-;C;FCC C60 Buckminsterfullerine");
    vproto.push_back("A2BCD3E6_cF208_203_e_c_d_f_g;5;52;203;7;6;cF208;a,x3,x4,x5,y5,z5;-;Na3Co(CO3)2Cl;Pyrochlore");
    vproto.push_back("A4B2C6D16E_cF232_203_e_d_f_eg_a;5;58;203;7;7;cF232;a,x3,x4,x5,x6,y6,z6;-;Na6Mg2(SO4)(CO3)4;Tychite");
    vproto.push_back("AB3C16_cF160_203_a_bc_eg;3;40;203;5;5;cF160;a,x4,x5,y5,z5;-;Rb3AsSe16;Rb3AsSe16");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("AB3C16_cF160_203_b_ad_eg;3;40;203;5;5;cF160;a,x4,x5,y5,z5;-;Rb3AsSe16;Rb3AsSe16");
    vproto.push_back("A2B3C6_cP264_205_2d_ab2c2d_6d;3;264;205;5;33;cP264;a,x3,x4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;Ca3Al2O6;Ca3Al2O6");
    vproto.push_back("A_cP240_205_10d;1;240;205;3;31;cP240;a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;C;Simple Cubic C60 Buckminsterfullerine");
    vproto.push_back("AB3C2_cI96_206_c_e_ad;3;48;206;5;6;cI96;a,x2,x3,x4,y4,z4;E9_{d};AlLi3N2;AlLi3N2");
    vproto.push_back("A17B15_cP64_207_acfk_eij;2;64;207;4;8;cP64;a,x3,x4,y5,y6,x7,y7,z7;-;Pd17Se15;Pd17Se15");
    vproto.push_back("A3B_cP16_208_j_b;2;16;208;4;2;cP16;a,x2;-;PH3;PH3");
    vproto.push_back("A6B2CD6E_cP64_208_m_ad_b_m_c;5;64;208;7;7;cP64;a,x5,y5,z5,x6,y6,z6;-;Cs2ZnFe[CN]6;Cs2ZnFe6");
    vproto.push_back("A24BC_cF104_209_j_a_b;3;26;209;5;4;cF104;a,x3,y3,z3;-;F6KP;F6KP");
    vproto.push_back("A12B36CD12_cF488_210_h_3h_a_fg;4;122;210;6;15;cF488;a,x2,y3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;MgB12H12[H2O]12;MgB12H1212"); //DX 20180925 - moved this structure from SG196 to SG210
    vproto.push_back("A12B6C_cF608_210_4h_2h_e;3;152;210;5;20;cF608;a,x1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Te[OH]6;Te6");
    vproto.push_back("A2B_cI72_211_hi_i;2;36;211;4;4;cI72;a,y1,y2,y3;-;SiO2;SiO2");
    vproto.push_back("A2B_cP12_212_c_a;2;12;212;4;2;cP12;a,x2;-;SrSi2;SrSi2");
    vproto.push_back("A3B3C_cI56_214_g_h_a;3;28;214;5;3;cI56;a,y2,y3;-;Ca3PI3;Ca3PI3");
    vproto.push_back("A3BC2_cI48_214_f_a_e;3;24;214;5;3;cI48;a,x2,x3;-;Ag3AuTe2;Petzite");
    vproto.push_back("A4B9_cP52_215_ei_3efgi;2;52;215;4;11;cP52;a,x1,x2,x3,x4,x5,x6,x7,z7,x8,z8;D8_{3};gamma-Cu9Al4;gamma-brass");
    vproto.push_back("ABCD_cF16_216_c_d_b_a;4;4;216;6;1;cF16;a;-;LiMgAuSn;Quaternary Heusler"); //DX 20180925 - fixed typo in "Quaternary"
    vproto.push_back("A3B4C_cP16_218_c_e_a;3;16;218;5;2;cP16;a,x3;H2_{1};Ag3[PO4];Ag3"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A7BC3D13_cF192_219_de_b_c_ah;4;48;219;6;5;cF192;a,x5,x6,y6,z6;-;Mg3B7ClO13;Boracite");
    vproto.push_back("A15B4_cI76_220_ae_c;2;38;220;4;5;cI76;a,x2,x3,y3,z3;D8_{6};Cu15Si4;Cu15Si4");
    vproto.push_back("A4B3_cI28_220_c_a;2;14;220;4;2;cI28;a,x2;D7_{3};Th3P4;Th3P4");
    vproto.push_back("A2B3C6_cP33_221_cd_ag_fh;3;33;221;5;4;cP33;a,x4,x5,x6;E9_{1};Ca3Al2O6;Ca3Al2O6");
    vproto.push_back("A5B3C16_cP96_222_ce_d_fi;3;96;222;5;6;cP96;a,x3,x4,x5,y5,z5;-;Ce5Mo3O16;Ce5Mo3O16");
    vproto.push_back("A23B6_cF116_225_bd2f_e;2;29;225;4;4;cF116;a,x3,x4,x5;D8_{a};Th6Mn23;Th6Mn23");
    vproto.push_back("A6B2C_cF36_225_e_c_a;3;9;225;5;2;cF36;a,x3;J1_{1};K2PtCl6;K2PtCl6");
    vproto.push_back("AB13_cF112_226_a_bi;2;28;226;4;3;cF112;a,y3,z3;D2_{3};NaZn13;NaZn13");
    vproto.push_back("A2B2C7_cF88_227_c_d_af;3;22;227;5;2;cF88;a,x4;-;Eu2Ir2O7;Pyrochlore Iridate");
    vproto.push_back("A3B4_cF56_227_ad_e;2;14;227;4;2;cF56;a,x3;D7_{2};Co3O4;Spinel");
    vproto.push_back("A5BCD6_cF416_228_eg_c_b_h;4;104;228;6;6;cF416;a,x3,y4,x5,y5,z5;-;CuCrCl5[NH3]6;CuCrCl56");
    vproto.push_back("A6B_cF224_228_h_c;2;56;228;4;4;cF224;a,x2,y2,z2;-;TeO6H6;TeO6H6");
    vproto.push_back("A3B10_cI52_229_e_fh;2;26;229;4;4;cI52;a,x1,x2,y3;D8_{1};gamma-Fe3Zn10;gamma-brass");
    vproto.push_back("A4B_cI10_229_c_a;2;5;229;4;1;cI10;a;-;beta-Hg4Pt;beta-Hg4Pt");
    vproto.push_back("A7B3_cI40_229_df_e;2;20;229;4;3;cI40;a,x2,x3;D8_{f};Ir3Ge7;Ir3Ge7");
    vproto.push_back("A2B3C12D3_cI160_230_a_c_h_d;4;80;230;6;4;cI160;a,x4,y4,z4;S1_{4};Co3Al2Si3O12;Garnet");
    //DX 20181130 - add Ohad's SQS structures - START
    // -------------------------------------------------------------------------
    // SQS (from O. Levy)
    // -------------------------------------------------------------------------
    vproto.push_back("AB_aP16_2_4i_4i;2;16;2;4;30;aP16;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;TaTi;TaTi");
    vproto.push_back("A5B11_mP16_6_2abc_2a3b3c;2;16;6;4;32;mP16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;Ta5Ti11;Ta5Ti11");
    vproto.push_back("AB3_mC32_8_4a_12a;2;16;8;4;36;mC32;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,z13,x14,z14,x15,z15,x16,z16;-;TaTi3;TaTi3");
    vproto.push_back("AB3_mC32_8_4a_4a4b;2;16;8;4;32;mC32;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;TaTi3;TaTi3");
    vproto.push_back("A3B13_oC32_38_ac_a2bcdef;2;16;38;4;18;oC32;a,b/a,c/a,z1,z2,z3,z4,x5,z5,x6,z6,y7,z7,y8,z8,x9,y9,z9;-;Ta3Ti13;Ta3Ti13");
    vproto.push_back("A3B5_oC32_38_abce_abcdf;2;16;38;4;18;oC32;a,b/a,c/a,z1,z2,z3,z4,x5,z5,x6,z6,y7,z7,y8,z8,x9,y9,z9;-;Ta3Ti5;Ta3Ti5");
    vproto.push_back("AB7_hR16_166_c_c2h;2;16;166;4;8;hR16;a,c/a,x1,x2,x3,z3,x4,z4;-;TaTi7;TaTi7");
    //DX 20181130 - add Ohad's SQS structures - END
    //DX 20181211 - add Corey's kesterite structure - START
    vproto.push_back("A2BCD4_tI16_82_ac_b_d_g;4;8;82;6;5;tI16;a,c/a,x5,y5,z5;-;Cu2(Zn,Fe)SnS4;Kesterite");
    //DX 20181211 - add Corey's kesterite structure - END
    // -------------------------------------------------------------------------
    // misc prototypes (from Y. Lederer) //DX 20200203
    // -------------------------------------------------------------------------
    vproto.push_back("AB3_mC8_12_a_di;2;4;12;4;6;mC8;a,b/a,c/a,beta,x3,z3;-/-/-;-/Li1N3/Li1N3;-/Li1N3 (ICSD #34675)/Li1N3 (ICSD #181559)"); //DX20201019 - added two metal nitrides
    vproto.push_back("AB_mC8_12_i_i;2;4;12;4;8;mC8;a,b/a,c/a,beta,x1,z1,x2,z2;-;-;-");
    vproto.push_back("AB3C4_mC16_12_a_di_2i;3;8;12;5;10;mC16;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5;-;-;-");
    vproto.push_back("ABC2_mC16_12_i_i_adi;3;8;12;5;10;mC16;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5;-;-;-");
    vproto.push_back("AB3_oP4_47_a_ct;2;4;47;4;4;oP4;a,b/a,c/a,z3;-;-;-");
    vproto.push_back("A3B_oP4_47_cr_a;2;4;47;4;4;oP4;a,b/a,c/a,z3;-;-;-");
    vproto.push_back("AB3C4_oP8_47_a_ct_egs;3;8;47;5;5;oP8;a,b/a,c/a,z5,z6;-;-;-");
    vproto.push_back("A3BC4_oP8_47_eq_g_bdt;3;8;47;5;5;oP8;a,b/a,c/a,z5,z6;-;-;-");
    vproto.push_back("AB_oP4_51_e_e;2;4;51;4;5;oP4;a,b/a,c/a,z1,z2;-;-;-");
    vproto.push_back("ABC2_oP8_51_e_e_2f;3;8;51;5;7;oP8;a,b/a,c/a,z1,z2,z3,z4;-;-;-");
    vproto.push_back("AB_oP4_59_a_a;2;4;59;4;5;oP4;a,b/a,c/a,z1,z2;-;-;-");
    vproto.push_back("ABC2_oP8_59_a_a_2b;3;8;59;5;7;oP8;a,b/a,c/a,z1,z2,z3,z4;-;-;-");
    vproto.push_back("ABC2_oC16_63_c_c_g;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,x3,y3;-;-;-");
    vproto.push_back("AB2C3_oC12_65_a_i_cj;3;6;65;5;5;oC12;a,b/a,c/a,y3,y4;-;-;-");
    vproto.push_back("A3BC4_oC16_65_ai_b_q;3;8;65;5;6;oC16;a,b/a,c/a,y3,x4,y4;-;-;-");
    vproto.push_back("A3BC4_oC16_65_bj_a_eh;3;8;65;5;5;oC16;a,b/a,c/a,x4,y5;-;-;-");
    vproto.push_back("AB3C4_oC16_65_a_bf_hi;3;8;65;5;5;oC16;a,b/a,c/a,x4,y5;-;-;-");
    vproto.push_back("AB2_oC6_65_a_i;2;3;65;4;4;oC6;a,b/a,c/a,y2;-;-2;-2");
    vproto.push_back("A3B_oC8_65_ai_b;2;4;65;4;4;oC8;a,b/a,c/a,y3;-;-;-");
    vproto.push_back("A3B_oC8_65_bj_a;2;4;65;4;4;oC8;a,b/a,c/a,y3;-;-;-");
    vproto.push_back("AB_oC8_65_i_i;2;4;65;4;5;oC8;a,b/a,c/a,y1,y2;-;-;-");
    vproto.push_back("ABC2_oC16_65_i_i_fh;3;8;65;5;6;oC16;a,b/a,c/a,x2,y3,y4;-;-;-");
    vproto.push_back("AB2C3_oI12_71_a_e_df;3;6;71;5;5;oI12;a,b/a,c/a,x3,x4;-;-;-");
    vproto.push_back("AB_tP2_123_a_b;2;2;123;4;2;tP2;a,c/a;-;-;-");
    vproto.push_back("AB_tP2_123_a_c;2;2;123;4;2;tP2;a,c/a;-;-;-");
    vproto.push_back("A2B_tP3_123_g_a;2;3;123;4;3;tP3;a,c/a,z2;-;-;-");
    vproto.push_back("A3B_tP4_123_abc_d;2;4;123;4;2;tP4;a,c/a;-;-;-");
    vproto.push_back("A3B_tP4_123_ag_b;2;4;123;4;3;tP4;a,c/a,z3;-;-;-");
    vproto.push_back("A3B_tP4_123_cf_a;2;4;123;4;2;tP4;a,c/a;-;-;-");
    vproto.push_back("AB3_tP4_123_a_bh;2;4;123;4;3;tP4;a,c/a,z3;-;-;-");
    vproto.push_back("ABC2_tP4_123_a_b_h;3;4;123;5;3;tP4;a,c/a,z3;-;-;-");
    vproto.push_back("ABC2_tP4_123_a_c_e;3;4;123;5;2;tP4;a,c/a;-;-;-");
    vproto.push_back("ABC2_tP4_123_a_d_bc;3;4;123;5;2;tP4;a,c/a;-;-;-");
    vproto.push_back("AB_tP4_123_g_g;2;4;123;4;4;tP4;a,c/a,z1,z2;-;-;-");
    vproto.push_back("A2BC3_tP6_123_g_b_ch;3;6;123;5;4;tP6;a,c/a,z3,z4;-;-;-");
    vproto.push_back("A3BC4_tP8_123_abc_d_i;3;8;123;5;3;tP8;a,c/a,z5;-;-;-");
    vproto.push_back("A3BC4_tP8_123_ag_b_2h;3;8;123;5;5;tP8;a,c/a,z3,z4,z5;-;-;-");
    vproto.push_back("A3BC4_tP8_123_cf_a_k;3;8;123;5;3;tP8;a,c/a,x4;-;-;-");
    vproto.push_back("AB3C4_tP8_123_a_bh_cdg;3;8;123;5;4;tP8;a,c/a,z5,z6;-;-;-");
    vproto.push_back("ABC2_tP8_123_h_h_abg;3;8;123;5;5;tP8;a,c/a,z3,z4,z5;-;-;-");
    vproto.push_back("ABC2_tP8_129_c_c_2c;3;8;129;5;6;tP8;a,c/a,z1,z2,z3,z4;-;-;-");
    vproto.push_back("A3B_tI8_139_ae_b;2;4;139;4;3;tI8;a,c/a,z3;-;-;-");
    vproto.push_back("AB3_tI8_139_a_bd;2;4;139;4;2;tI8;a,c/a;-;-;-");
    vproto.push_back("AB2C3_tI12_139_a_e_be;3;6;139;5;4;tI12;a,c/a,z3,z4;-;-;-");
    vproto.push_back("A3BC4_tI16_139_ae_b_g;3;8;139;5;4;tI16;a,c/a,z3,z4;-;-;-");
    vproto.push_back("AB3C4_tI16_139_a_bd_ce;3;8;139;5;3;tI16;a,c/a,z5;-;-;-");
    vproto.push_back("ABC2_tI16_139_e_e_cd;3;8;139;5;4;tI16;a,c/a,z3,z4;-;-;-");
    vproto.push_back("ABC2_tI16_141_a_b_e;3;8;141;5;3;tI16;a,c/a,z3;-;-;-");
    vproto.push_back("A2B_hP3_164_d_a;2;3;164;4;3;hP3;a,c/a,z2;-;-;-");
    vproto.push_back("A2BC3_hP6_164_d_a_bd;3;6;164;5;4;hP6;a,c/a,z3,z4;-/-;-/Bi2La1Li3;-/Bi2La1Li3 (ICSD #616769)"); //DX20201028 - added metal prototype
    vproto.push_back("A2BC3_hP6_164_d_b_ad;3;6;164;5;4;hP6;a,c/a,z3,z4;-;-;-");
    vproto.push_back("A3B_hR4_166_bc_a;2;4;166;4;3;hR4;a,c/a,x3;-/-;-/N3Na1;-/N3Na1 (ICSD #1144)"); //DX20201019 - added metal-nitride
    vproto.push_back("AB3_hR4_166_a_bc;2;4;166;4;3;hR4;a,c/a,x3;-;-;-");
    vproto.push_back("AB_hR4_166_c_c;2;4;166;4;4;hR4;a,c/a,x1,x2;-;-;-");
    vproto.push_back("A3BC4_hR8_166_bc_a_2c;3;8;166;5;5;hR8;a,c/a,x3,x4,x5;-;-;-");
    vproto.push_back("AB3C4_hR8_166_a_bc_2c;3;8;166;5;5;hR8;a,c/a,x3,x4,x5;-;-;-");
    vproto.push_back("ABC2_hR8_166_c_c_abc;3;8;166;5;5;hR8;a,c/a,x3,x4,x5;-;-;-");
    vproto.push_back("AB3C4_cP8_221_a_c_bd;3;8;221;5;1;cP8;a;-;-;-");
    // -------------------------------------------------------------------------
    // oxide prototypes (from R. Friedrich) //DX20200624
    // -------------------------------------------------------------------------
    // binaries
    vproto.push_back("AB_mC20_12_b2i_c2i;2;10;12;4;12;mC20;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,z6;-;OTi;OTi (ICSD #56694)");
    vproto.push_back("A2B3_mC20_12_2i_3i;2;10;12;4;14;mC20;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;Ga2O3;Ga2O3 (ICSD #83645)");
    vproto.push_back("A5B3_mC32_12_5i_3i;2;16;12;4;20;mC32;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;O5Ti3;O5Ti3 (ICSD #75193)");
    vproto.push_back("A3B_mP16_14_3e_e;2;16;14;4;16;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;O3W;O3W (ICSD #84848)");
    vproto.push_back("A2B3_mP20_14_2e_3e;2;20;14;4;19;mP20;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Bi2O3;Bi2O3 (ICSD #168806)");
    vproto.push_back("AB_mC8_15_a_e;2;4;15;4;5;mC8;a,b/a,c/a,beta,y2;-;CuO;CuO (ICSD #92368)");
    vproto.push_back("A4B_mC20_15_2f_e;2;10;15;4;11;mC20;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3;-;O4Os;O4Os (ICSD #23803)");
    vproto.push_back("A2B5_oP28_19_2a_5a;2;28;19;4;24;oP28;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;As2O5;As2O5 (ICSD #654040)");
    vproto.push_back("A2B_oP24_33_4a_2a;2;24;33;4;21;oP24;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;O2Sb;O2Sb (ICSD #4109)");
    vproto.push_back("AB3_oC16_40_b_a2b;2;8;40;4;10;oC16;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4;-;CrO3;CrO3 (ICSD #16031)"); //DX20200805 - created with bad geometry file; fixed
    vproto.push_back("A5B2_oP14_59_a2e_e;2;14;59;4;10;oP14;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4;-;O5V2;O5V2 (ICSD #60767)");
    vproto.push_back("AB_oC16_64_e_f;2;8;64;4;6;oC16;a,b/a,c/a,y1,y2,z2;-;KO;KO (ICSD #180559)");
    vproto.push_back("A4B3_tP28_135_gh_dh;2;28;135;4;7;tP28;a,c/a,x2,x3,y3,x4,y4;-;O4Pb3;O4Pb3:vs (ICSD #29094)");
    vproto.push_back("A2B3_hP15_144_2a_3a;2;15;144;4;17;hP15;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;B2O3;B2O3 (ICSD #24649)");
    vproto.push_back("A2B_hR3_166_c_a;2;3;166;4;3;hR3;a,c/a,x2;-/-/-;Cs2O/Ca2N1/Ba2N1;Cs2O (ICSD #27919)/Ca2N1 (ICSD #22231)/Ba2N1 (ICSD #409851)"); //DX20201019 - added two metal-nitrides
    vproto.push_back("AB2_hR6_166_c_2c;2;6;166;4;5;hR6;a,c/a,x1,x2,x3;-;OTl2;OTl2 (ICSD #16220)");
    vproto.push_back("AB_hP12_189_fg_eh;2;12;189;4;6;hP12;a,c/a,z1,x2,x3,z4;-;NaO;NaO (ICSD #25526)");
    vproto.push_back("AB_hP8_194_ac_f;2;8;194;4;3;hP8;a,c/a,z3;-;LiO;LiO (ICSD #152183)");
    vproto.push_back("A2B3_cI80_199_a2b_2c;2;40;199;4;10;cI80;a,x1,x2,x3,x4,y4,z4,x5,y5,z5;-;In2O3;In2O3 (ICSD #33649)");
    vproto.push_back("A3B2_cF80_227_f_e;2;20;227;4;3;cF80;a,x1,x2;-;O3Sb2;O3Sb2 (ICSD #31102)");
    // ternaries
    vproto.push_back("A4B4C_aP18_2_4i_4i_i;3;18;2;5;33;aP18;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Na4O4Si;Na4O4Si (ICSD #15500)");
    vproto.push_back("A2B7C2_aP22_2_2i_7i_2i;3;22;2;5;39;aP22;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-/-;Ca2O7V2/Mg2O7V2;Ca2O7V2 (ICSD #421266)/Mg2O7V2 (ICSD #2321)");
    vproto.push_back("AB3C_aP30_2_3i_9i_3i;3;30;2;5;51;aP30;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;CaO3Si;CaO3Si (ICSD #201537)");
    vproto.push_back("A2B5C_aP32_2_4i_10i_2i;3;32;2;5;54;aP32;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;Al2O5Si;Al2O5Si (ICSD #85742)");
    vproto.push_back("A2B4C_mP28_4_4a_8a_2a;3;28;4;5;46;mP28;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;Al2O4Sr;Al2O4Sr (ICSD #160297)");
    vproto.push_back("A3BC_mC60_5_ab8c_ab2c_3c;3;30;5;5;47;mC60;a,b/a,c/a,beta,y1,y2,y3,y4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17;-;O3SiSr;O3SiSr (ICSD #32542)");
    vproto.push_back("A3B5C_mC54_8_3a3b_9a3b_3a;3;27;8;5;52;mC54;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,z13,x14,z14,x15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21;-;Ca3O5Si;Ca3O5Si (ICSD #81100)");
    vproto.push_back("A2B7C3_mP24_11_2e_7e_3e;3;24;11;5;28;mP24;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12;-;Na2O7Ti3;Na2O7Ti3 (ICSD #15463)");
    vproto.push_back("AB6C2_mC18_12_a_3i_i;3;9;12;5;12;mC18;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5;-/-;CaO6V2/MgO6V2;CaO6V2 (ICSD #166516)/MgO6V2 (ICSD #10391)");
    vproto.push_back("ABC4_mC48_12_gi_hi_2i3j;3;24;12;5;23;mC48;a,b/a,c/a,beta,y1,y2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-/-;MgMoO4/FeMoO4;MgMoO4 (ICSD #20418)/FeMoO4 (ICSD #43012)");
    vproto.push_back("AB4C_mP12_13_f_2g_e;3;12;13;5;12;mP12;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4;-/-;MgO4W/NiO4W;MgO4W (ICSD #67903)/NiO4W (ICSD #16685)");
    vproto.push_back("AB3C_mP40_14_2e_6e_2e;3;40;14;5;34;mP40;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;MgO3Si;MgO3Si (ICSD #30895)");
    vproto.push_back("A2B3C_mC24_15_2e_af_e;3;12;15;5;10;mC24;a,b/a,c/a,beta,y2,y3,y4,x5,y5,z5;-;Li2O3Zr;Li2O3Zr (ICSD #94893)");
    vproto.push_back("AB3C_mC40_15_2e_3f_f;3;20;18;5;18;mC40;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;NaO3V;NaO3V (ICSD #2103)");
    vproto.push_back("A2B3C_mC48_15_aef_3f_2e;3;24;15;5;19;mC48;a,b/a,c/a,beta,y2,y3,y4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Li2O3Ti;Li2O3Ti (ICSD #162215)");
    vproto.push_back("A4BC7_mC48_15_2f_e_e3f;3;24;15;5;21;mC48;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Al4CaO7;Al4CaO7 (ICSD #14270)");
    vproto.push_back("AB3C_mC60_15_cf_e4f_ef;3;30;15;5;24;mC60;a,b/a,c/a,beta,y2,y3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;CaO3Si;CaO3Si (ICSD #87694)");
    vproto.push_back("ABC2_oP16_33_a_a_2a;3;16;33;5;15;oP16;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;FeNaO2;FeNaO2 (ICSD #186309)");
    vproto.push_back("A2B3C_oC24_36_b_ab_a;3;12;36;5;13;oC24;a,b/a,c/a,y1,z1,y2,z2,x3,y3,z3,x4,y4,z4;-/-;Li2O3Si/Na2O3Si;Li2O3Si (ICSD #28192)/Na2O3Si (ICSD #24664)");
    vproto.push_back("A2B5C_oP32_58_eg_3gh_g;3;32;58;5;17;oP32;a,b/a,c/a,z1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,z7;-;Al2O5Si;Al2O5Si (ICSD #30679)");
    vproto.push_back("AB2C4_oC28_63_c_ac_fg;3;14;63;5;9;oC28;a,b/a,c/a,y2,y3,y4,z4,x5,y5;-;CrNa2O4;CrNa2O4 (ICSD #76001)");
    vproto.push_back("AB5C2_oC32_63_c_c2f_f;3;16;63;5;11;oC32;a,b/a,c/a,y1,y2,y3,z3,y4,z4,y5,z5;-;MgO5Ti2;MgO5Ti2 (ICSD #37232)");
    vproto.push_back("ABC4_tI24_88_a_b_f;3;12;88;5;5;tI24;a,c/a,x3,y3,z3;-;BaMoO4;BaMoO4 (ICSD #50821)");
    vproto.push_back("A4BC2_tP28_91_2d_b_ac;3;28;91;5;11;tP28;a,c/a,y1,y2,x3,x4,y4,z4,x5,y5,z5;-;O4TiZn2;O4TiZn2 (ICSD #109093)");
    vproto.push_back("A2BC4_hP56_173_2b2c_ac_b5c;3;56;173;5;30;hP56;a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;Al2BaO4;Al2BaO4 (ICSD #246027)");
    vproto.push_back("AB2C4_cF56_227_b_c_e;3;14;227;5;2;cF56;a,x3;-;MoNa2O4;MoNa2O4 (ICSD #151970)");
    // -------------------------------------------------------------------------
    // nitrides prototypes (from R. Friedrich) //DX20200624
    // -------------------------------------------------------------------------
    // binaries
    vproto.push_back("A2B3_cI80_206_ad_e;2;40;206;4;5;cI80;a,x2,x3,y3,z3;-;N2Zn3;N2Zn3 (ICSD #84918)");
    // ternaries
    vproto.push_back("ABC_oP12_62_c_c_c;3;12;62;5;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;-/-/-/-/-/-/-/-/-/-/-/-/-/-;CaLiN/Li1Mg1N1/Ca1Mg1Sn1/Ni1Sn1Tb1/Co1Dy1Sn1/Nd1Pd1Sn1/Al1Pd1Tb1/Bi1Ca1Li1/Ni1Sc1Sn1/Co1Dy1Sn1/La1Ni1Sn1/La1Ni1Sn1/Ni1Pr1Sn1/Er1Ga1Rh1;CaLiN (ICSD #107304)/Li1Mg1N1 (ICSD #93259)/Ca1Mg1Sn1 (ICSD #42757)/Ni1Sn1Tb1 (ICSD #54301)/Co1Dy1Sn1 (ICSD #54415)/Nd1Pd1Sn1 (ICSD #54460)/Al1Pd1Tb1 (ICSD #54938)/Bi1Ca1Li1 (ICSD #58762)/Ni1Sc1Sn1 (ICSD #105338)/Co1Dy1Sn1 (ICSD #106459)/La1Ni1Sn1 (ICSD #108571)/La1Ni1Sn1 (ICSD #157921)/Ni1Pr1Sn1 (ICSD #157923)/Er1Ga1Rh1 (ICSD #630573)"); //DX20201016 - added another metal-nitride system //DX20201028 - added 12 metal systems
    vproto.push_back("A2B2C_tI10_139_e_e_a;3;5;139;5;4;tI10;a,c/a,z2,z3;-/-/-;Ca2N2Zn/N2Sr2Zn1/Ba2N2Zn1;Ca2N2Zn (ICSD #69049)/N2Sr2Zn1 (ICSD #80376)/Ba2N2Zn1 (ICSD #80377)"); //DX20201016 - added two metal-nitride prototypes
    vproto.push_back("A3BC2_hP6_191_f_a_d;3;6;191;5;2;hP6;a,c/a;-;Na3TaTi2;Na3TaTi2 (ICSD #186418)");
    // -------------------------------------------------------------------------
    // metal-nitride prototypes (from DX) //DX2021016
    // -------------------------------------------------------------------------
    // binaries
    vproto.push_back("A2B_mP6_6_3ab_ab;2;6;6;4;16;mP6;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;N2Re1;N2Re1 (ICSD #187448)");
    vproto.push_back("A2B_mC12_8_4a_2a;2;6;8;4;16;mC12;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;N2Re1;N2Re1 (ICSD #187447)");
    vproto.push_back("A3B_mC16_8_6a_b;2;8;8;4;19;mC16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;-;N3Rb1;N3Rb1 (ICSD #155169)");
    vproto.push_back("A2B_mP6_11_2e_e;2;6;11;4;10;mP6;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3;-;N2Re1;N2Re1 (ICSD #187444)");
    vproto.push_back("AB6_mP14_11_e_6e;2;14;11;4;18;mP14;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;Ba1N6;Ba1N6 (ICSD #14244)");
    vproto.push_back("A3B_mC8_12_di_a;2;4;12;4;6;mC8;a,b/a,c/a,beta,x3,z3;-/-;N3Na1/N3Na1;N3Na1 (ICSD #29370)/N3Na1 (ICSD #29376)");
    vproto.push_back("A2B_mC12_12_2i_i;2;6;12;4;10;mC12;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3;-;N2Re1;N2Re1 (ICSD #187441)");
    vproto.push_back("AB_mC16_12_2i_2i;2;8;12;4;12;mC16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;N1Sr1;N1Sr1 (ICSD #411555)");
    vproto.push_back("A3B2_mC30_12_a4i_3i;2;15;12;4;18;mC30;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-/-;Ca3N2/Ca3N2;Ca3N2 (ICSD #162794)/Ca3N2 (ICSD #169726)");
    vproto.push_back("A2B_mP6_13_g_e;2;6;13;4;8;mP6;a,b/a,c/a,beta,y1,x2,y2,z2;-;N2Re1;N2Re1 (ICSD #187451)");
    vproto.push_back("AB2_mP12_14_e_2e;2;12;14;4;13;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;Ir1N2;Ir1N2 (ICSD #160623)");
    vproto.push_back("AB3_mP16_14_e_3e;2;16;14;4;16;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Hg2N6;Hg2N6 (ICSD #98661)");
    vproto.push_back("AB2_mC12_15_e_f;2;6;15;4;8;mC12;a,b/a,c/a,beta,y1,x2,y2,z2;-;Ba1N2;Ba1N2 (ICSD #280681)");
    vproto.push_back("AB_oC16_36_2a_2a;2;8;36;4;11;oC16;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4;-;N1Os1;N1Os1 (ICSD #167514)");
    vproto.push_back("AB_oP8_53_h_h;2;8;53;4;7;oP8;a,b/a,c/a,y1,z1,y2,z2;-;N1Os1;N1Os1 (ICSD #167513)");
    vproto.push_back("A2B_oP12_57_e_c;2;12;57;4;7;oP12;a,b/a,c/a,x1,x2,y2,z2;-;N2Re1;N2Re1 (ICSD #187450)");
    vproto.push_back("A5B3_oC32_63_c2f_cf;2;16;63;4;11;oC32;a,b/a,c/a,y1,y2,y3,z3,y4,z4,y5,z5;-/-;N5Ta3/N5Ta3;N5Ta3 (ICSD #16253)/N5Ta3 (ICSD #66533)");
    vproto.push_back("A2B_oC12_65_gj_e;2;6;65;4;5;oC12;a,b/a,c/a,x2,y3;-;N2Re1;N2Re1 (ICSD #187452)");
    vproto.push_back("AB_oI4_71_b_c;2;4;71;4;3;oI4;a,b/a,c/a;-;Cr1N1;Cr1N1 (ICSD #53146)");
    vproto.push_back("A2B_oI12_71_2e_f;2;6;71;4;6;oI12;a,b/a,c/a,x1,x2,x3;-;N2Re1;N2Re1 (ICSD #187449)");
    vproto.push_back("AB_oI16_71_abe_n;2;8;71;4;6;oI16;a,b/a,c/a,x3,x4,y4;-;Li2N2;Li2N2 (ICSD #423831)");
    vproto.push_back("AB3_oI16_72_b_cj;2;8;72;4;5;oI16;a,b/a,c/a,x3,y3;-;Ag1N3;Ag1N3 (ICSD #27135)");
    vproto.push_back("A5B4_tI18_87_ah_h;2;9;87;4;6;tI18;a,c/a,x2,y2,x3,y3;-;N5Nb4;N5Nb4 (ICSD #26251)");
    vproto.push_back("AB3_tI32_88_d_cf;2;16;88;4;5;tI32;a,c/a,x3,y3,z3;-;Cu1N3;Cu1N3 (ICSD #30633)");
    vproto.push_back("AB_tP2_123_d_a;2;2;123;4;2;tP2;a,c/a;-;N1Pr1;N1Pr1 (ICSD #168645)");
    vproto.push_back("A3B_tP4_123_ag_d;2;4;123;4;3;tP4;a,c/a,z3;-;N3Rb1;N3Rb1 (ICSD #16963)");
    vproto.push_back("A3B_tP8_123_egh_ab;2;8;123;4;4;tP8;a,c/a,z4,z5;-;Cu3N1;Cu3N1 (ICSD #180237)");
    vproto.push_back("A2B_tP6_127_g_b;2;6;127;4;3;tP6;a,c/a,x2;-;N2Re1;N2Re1 (ICSD #187442)");
    vproto.push_back("AB_tP4_131_e_c;2;4;131;4;2;tP4;a,c/a;-;N1Tc1;N1Tc1 (ICSD #187709)");
    vproto.push_back("AB_tI4_139_a_b;2;2;139;4;2;tI4;a,c/a;-;Mn1N1;Mn1N1 (ICSD #106932)");
    vproto.push_back("A3B2_tI10_139_ae_e;2;5;139;4;4;tI10;a,c/a,z2,z3;-/-;Mn3N2/Mn3N2;Mn3N2 (ICSD #71638)/Mn3N2 (ICSD #84202)");
    vproto.push_back("A3B4_tI14_139_ad_ce;2;7;139;4;3;tI14;a,c/a,z4;-;N3Nb4;N3Nb4 (ICSD #76389)");
    vproto.push_back("A8B_tI18_139_deh_a;2;9;139;4;4;tI18;a,c/a,z3,x4;-;Fe16N2;Fe16N2 (ICSD #41953)");
    vproto.push_back("AB3_tI16_140_a_dh;2;8;140;4;3;tI16;a,c/a,x3;-/-/-;K1N3/Cs1N3/Ag1N3;K1N3 (ICSD #24007)/Cs1N3 (ICSD #25008)/Ag1N3 (ICSD #183201)");
    vproto.push_back("AB3_tI16_140_a_ch;2;8;140;4;3;tI16;a,c/a,x3;-;Cu1N3;Cu1N3 (ICSD #187018)");
    vproto.push_back("AB3_tI16_140_c_ah;2;8;140;4;3;tI16;a,c/a,x3;-;K1N3;K1N3 (ICSD #1145)");
    vproto.push_back("AB2_tI12_141_a_e;2;6;141;4;3;tI12;a,c/a,z2;-;N1Ti2;N1Ti2 (ICSD #23403)");
    vproto.push_back("A3B_hR4_160_3a_a;2;4;160;4;6;hR4;a,c/a,x1,x2,x3,x4;-;N3Na1;N3Na1 (ICSD #644523)");
    vproto.push_back("A3B2_hR5_160_3a_2a;2;5;160;4;7;hR5;a,c/a,x1,x2,x3,x4,x5;-;Be3N2;Be3N2 (ICSD #185490)");
    vproto.push_back("AB_hP16_162_ek_ci;2;16;162;4;6;hP16;a,c/a,z2,x3,x4,z4;-;Mo1N1;Mo1N1 (ICSD #43559)");
    vproto.push_back("A3B2_hP5_164_ac_d;2;5;164;4;4;hP5;a,c/a,z2,z3;-;N3V2;N3V2 (ICSD #182700)");
    vproto.push_back("AB_hP16_164_ci_di;2;16;164;4;8;hP16;a,c/a,z1,z2,x3,z3,x4,z4;-;Mo1N1;Mo1N1 (ICSD #60168)");
    vproto.push_back("AB_hP16_186_ac_bc;2;16;186;4;8;hP16;a,c/a,z1,z2,x3,z3,x4,z4;-;Mo1N1;Mo1N1 (ICSD #76280)");
    vproto.push_back("AB2_hP3_187_a_h;2;3;187;4;3;hP3;a,c/a,z2;-;Hf1N2;Hf1N2 (ICSD #290427)");
    vproto.push_back("A2B_hP3_191_e_a;2;3;191;4;3;hP3;a,c/a,z2;-;N2Os1;N2Os1 (ICSD #260545)");
    vproto.push_back("A2B_hP6_191_af_d;2;6;191;4;2;hP6;a,c/a;-;N4Ta2;N4Ta2 (ICSD #182351)");
    vproto.push_back("AB_hP6_191_ad_f;2;6;191;4;2;hP6;a,c/a;-;N1Ta1;N1Ta1 (ICSD #76457)");
    vproto.push_back("AB_hP6_191_bc_ad;2;6;191;4;2;hP6;a,c/a;-;N1Y1;N1Y1 (ICSD #161079)");
    vproto.push_back("A3B_hP8_193_g_b;2;8;193;4;3;hP8;a,c/a,x2;-;Ba3N1;Ba3N1 (ICSD #77730)");
    vproto.push_back("AB_hP4_194_a_c;2;4;194;4;2;hP4;a,c/a;-;Ir1N1;Ir1N1 (ICSD #183157)");
    vproto.push_back("AB2_hP6_194_c_e;2;6;194;4;3;hP6;a,c/a,z2;-;Hf1N2;Hf1N2 (ICSD #290428)");
    vproto.push_back("A3B2_hP10_194_bf_ac;2;10;194;4;3;hP10;a,c/a,z4;-/-;Be3N2/Ca3N2;Be3N2 (ICSD #25656)/Ca3N2 (ICSD #162797)");
    vproto.push_back("A3B2_hP10_194_cf_f;2;10;194;4;4;hP10;a,c/a,z2,z3;-;N3W2;N3W2 (ICSD #186207)");
    vproto.push_back("AB_hP16_194_bh_ag;2;16;194;4;3;hP16;a,c/a,x4;-;Mo1N1;Mo1N1 (ICSD #106926)");
    vproto.push_back("AB3_cI32_204_c_g;2;16;204;4;3;cI32;a,y2,z2;-/-;Co1N3/Co1N3;Co1N3 (ICSD #162105)/Co1N3 (ICSD #162106)");
    vproto.push_back("AB4_cP5_221_a_bc;2;5;221;4;1;cP5;a;-;N1Ni4;N1Ni4 (ICSD #76403)");

    // -------------------------------------------------------------------------
    // ternaries
    vproto.push_back("AB2C5_aP8_1_a_2a_5a;3;8;1;5;30;aP8;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Li1N2Na5;Li1N2Na5 (ICSD #92316)");
    vproto.push_back("AB2C5_mC16_5_a_c_a2bc;3;8;5;5;14;mC16;a,b/a,c/a,beta,y1,y2,y3,y4,x5,y5,z5,x6,y6,z6;-;Li1N2Na5;Li1N2Na5 (ICSD #92315)");
    vproto.push_back("A3B2C3_mP8_6_2ab_ab_a2b;3;8;6;5;20;mP8;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;Li3N2Na3;Li3N2Na3 (ICSD #92312)");
    vproto.push_back("A5B2C_mP8_6_3a2b_ab_b;3;8;6;5;20;mP8;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;Li5N2Na1;Li5N2Na1 (ICSD #92314)");
    vproto.push_back("A2BC2_mC20_12_2i_i_2i;3;10;12;5;14;mC20;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;Ca2Fe1N2;Ca2Fe1N2 (ICSD #72389)");
    vproto.push_back("AB2C12_mC30_12_a_i_6i;3;15;12;5;18;mC30;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-/-;Cd1K2N12/Cd1K2N12;Cd1K2N12 (ICSD #31297)/Cd1K2N12 (ICSD #659621)");
    vproto.push_back("ABC_mP12_14_e_e_e;3;12;14;5;13;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;Be1Li1N1;Be1Li1N1 (ICSD #402341)");
    vproto.push_back("A2BC_oC16_63_2c_c_c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,y4;-;Ca2In1N1;Ca2In1N1 (ICSD #96228)");
    vproto.push_back("ABC3_oC20_63_c_a_cf;3;10;63;5;7;oC20;a,b/a,c/a,y2,y3,y4,z4;-;Al1N1Zr3;Al1N1Zr3 (ICSD #29521)");
    vproto.push_back("A3BC3_oC28_63_cg_c_cg;3;14;63;5;10;oC28;a,b/a,c/a,y1,y2,y3,x4,y4,x5,y5;-;Ca3Cr1N3;Ca3Cr1N3 (ICSD #40205)");
    vproto.push_back("A2BC_oC8_65_i_c_b;3;4;65;5;4;oC8;a,b/a,c/a,y3;-;Li4N2Na2;Li4N2Na2 (ICSD #92307)");
    vproto.push_back("A7B2C3_oC24_65_a2ij_i_cj;3;12;65;5;8;oC24;a,b/a,c/a,y3,y4,y5,y6,y7;-;Ca7N2Tl3;Ca7N2Tl3 (ICSD #98178)");
    vproto.push_back("A2BC3_oI12_71_e_a_de;3;6;71;5;5;oI12;a,b/a,c/a,x3,x4;-;Ce2Mn1N3;Ce2Mn1N3 (ICSD #50579)");
    vproto.push_back("ABC_tP6_99_ab_ab_ab;3;6;99;5;8;tP6;a,c/a,z1,z2,z3,z4,z5,z6;-;Cr1N1Nb1;Cr1N1Nb1 (ICSD #23779)");
    vproto.push_back("A3B2C3_tP8_115_abc_g_de;3;8;115;5;4;tP8;a,c/a,z5,z6;-;Li3N2Na3;Li3N2Na3 (ICSD #92311)");
    vproto.push_back("ABC_tP3_123_c_a_b;3;3;123;5;2;tP3;a,c/a;-;Fe1N1Ni1;Fe1N1Ni1 (ICSD #53505)");
    vproto.push_back("A5B2C_tP8_123_ai_bc_d;3;8;123;5;3;tP8;a,c/a,z5;-;Li5N2Na1;Li5N2Na1 (ICSD #92313)");
    vproto.push_back("ABC_tP6_129_c_c_c;3;6;129;5;5;tP6;a,c/a,z1,z2,z3;-/-;Ca1Ga1N1/Mo1N1Ta1;Ca1Ga1N1 (ICSD #2027)/Mo1N1Ta1 (ICSD #100437)");
    vproto.push_back("ABC2_tP8_129_c_c_bc;3;8;129;5;5;tP8;a,c/a,z2,z3,z4;-;Ba1Hf1N2;Ba1Hf1N2 (ICSD #50994)");
    vproto.push_back("ABC2_tP8_129_a_c_bc;3;8;129;5;4;tP8;a,c/a,z3,z4;-;Li2N2Na4;Li2N2Na4 (ICSD #92309)");
    vproto.push_back("ABC_tP6_131_e_c_b;3;6;131;5;2;tP6;a,c/a;-;Ca1N1Ni1;Ca1N1Ni1 (ICSD #69044)");
    vproto.push_back("ABC_tP6_131_b_c_e;3;6;131;5;2;tP6;a,c/a;-;Li1N1Sr1;Li1N1Sr1 (ICSD #87414)");
    vproto.push_back("A2BC_tI8_139_d_b_a;3;4;139;5;2;tI8;a,c/a;-;Li4N2Na2;Li4N2Na2 (ICSD #92305)");
    vproto.push_back("A3B6C2_tI22_139_ae_eg_e;3;11;139;5;6;tI22;a,c/a,z2,z3,z4,z5;-/-;La3N6V2/La3N6Nb2;La3N6V2 (ICSD #98477)/La3N6Nb2 (ICSD #411473)");
    vproto.push_back("A2B2C_tI20_140_h_h_a;3;10;140;5;4;tI20;a,c/a,x2,x3;-;Be4N4Sr2;Be4N4Sr2 (ICSD #413356)");
    vproto.push_back("A2BC2_tI20_140_h_a_h;3;10;140;5;4;tI20;a,c/a,x2,x3;-;Be4Ca2N4;Be4Ca2N4 (ICSD #413357)");
    vproto.push_back("AB2C2_tI20_140_a_h_h;3;10;140;5;4;tI20;a,c/a,x2,x3;-;Ba1Be2N2;Ba1Be2N2 (ICSD #415304)");
    vproto.push_back("AB2C_hR4_160_a_2a_a;3;4;160;5;6;hR4;a,c/a,x1,x2,x3,x4;-;Cr1N2W1;Cr1N2W1 (ICSD #84639)");
    vproto.push_back("A2B2C_hP5_164_d_d_a;3;5;164;5;4;hP5;a,c/a,z2,z3;-/-;Li2N2Zr1/Mg2N2Sr1;Li2N2Zr1 (ICSD #16231)/Mg2N2Sr1 (ICSD #410826)");
    vproto.push_back("AB2C2_hP5_164_a_d_d;3;5;164;5;4;hP5;a,c/a,z2,z3;-;Ce1Li2N2;Ce1Li2N2 (ICSD #34003)");
    vproto.push_back("AB2C_hR4_166_a_c_b;3;4;166;5;3;hR4;a,c/a,x3;-/-;Cu1N2Ta1/Ba1Cu2Ga1;Cu1N2Ta1 (ICSD #71136)/Ba1Cu2Ga1 (ICSD #615828)"); //DX20201028 - added metal prototype
    vproto.push_back("A3BC3_hP14_176_h_c_h;3;14;176;5;6;hP14;a,c/a,x2,y2,x3,y3;-;Ba3Fe1N3;Ba3Fe1N3 (ICSD #36502)");
    vproto.push_back("ABC_hP6_186_b_b_a;3;6;186;5;5;hP6;a,c/a,z1,z2,z3;-/-/-;N1Na1Sn1/Cu1Sn1Ti1/Bi1Cu1Yb1;N1Na1Sn1 (ICSD #172471)/Cu1Sn1Ti1 (ICSD #54657)/Bi1Cu1Yb1 (ICSD #416143)"); //DX20201028 - added two metal prototypes
    vproto.push_back("AB2C_hP8_186_a_ab_b;3;8;186;5;6;hP8;a,c/a,z1,z2,z3,z4;-;Mn1N2W1;Mn1N2W1 (ICSD #80029)");
    vproto.push_back("ABC_hP3_187_a_c_d;3;3;187;5;2;hP3;a,c/a;-;Li1N1Ni1;Li1N1Ni1 (ICSD #247028)");
    vproto.push_back("A5B3C3_hP11_189_dg_g_f;3;11;189;5;5;hP11;a,c/a,x2,x3,x4;-;Li5N3Ni3;Li5N3Ni3 (ICSD #411152)");
    vproto.push_back("A2BC_hP4_191_c_a_b;3;4;191;5;2;hP4;a,c/a;-;Li2N1Na1;Li2N1Na1 (ICSD #92308)");
    vproto.push_back("ABC2_hP4_191_b_a_c;3;4;191;5;2;hP4;a,c/a;-;Li1N1Na2;Li1N1Na2 (ICSD #92310)");
    vproto.push_back("ABC2_hP8_194_c_a_f;3;8;194;5;3;hP8;a,c/a,z3;-;Ba1Ce1N2;Ba1Ce1N2 (ICSD #74791)");
    vproto.push_back("AB2C_hP8_194_a_f_b;3;8;194;5;3;hP8;a,c/a,z3;-;Fe1N2W1;Fe1N2W1 (ICSD #75971)");
    vproto.push_back("ABC2_hP8_194_a_c_f;3;8;194;5;3;hP8;a,c/a,z3;-;Mg1Mo1N2;Mg1Mo1N2 (ICSD #185913)");
    vproto.push_back("A3BC_hP10_194_h_a_c;3;10;194;5;3;hP10;a,c/a,x3;-;Ba3N1Na1;Ba3N1Na1 (ICSD #67497)");
    vproto.push_back("A3BC_hP10_194_h_c_a;3;10;194;5;3;hP10;a,c/a,x3;-;Ba3Bi1N1;Ba3Bi1N1 (ICSD #152055)");
    vproto.push_back("AB3C4_hP16_194_c_ae_2f;3;16;194;5;5;hP16;a,c/a,z3,z4,z5;-;Al1N3V4;Al1N3V4 (ICSD #181351)");
    vproto.push_back("AB2C_cF32_227_b_c_a;3;8;227;5;1;cF32;a;-;Cs1N2Nb1;Cs1N2Nb1 (ICSD #72546)");

    // -------------------------------------------------------------------------
    // metal prototypes (from DX) //DX2021028
    // -------------------------------------------------------------------------
    // ternaries
    vproto.push_back("ABC2_mC16_8_b_2a_2ab;3;8;8;5;18;mC16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5,x6,y6,z6;-;Bi1Pb1Pd2;Bi1Pb1Pd2 (ICSD #58830)");
    vproto.push_back("A2B4C3_mP9_10_m_2n_am;3;9;10;5;12;mP9;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5;-;Co2Gd4Mg3;Co2Gd4Mg3 (ICSD #417035)");
    vproto.push_back("A2B5C4_mP11_10_m_a2m_2n;3;11;10;5;14;mP11;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;Co2In5Zr4;Co2In5Zr4 (ICSD #55578)");
    vproto.push_back("AB2C2_mP10_11_e_2e_2e;3;10;11;5;14;mP10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;Ba1Bi2Pd2;Ba1Bi2Pd2 (ICSD #416299)");
    vproto.push_back("A2BC2_mP10_11_2e_e_2e;3;10;11;5;14;mP10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;Sn2Sr1Zn2;Sn2Sr1Zn2 (ICSD #424108)");
    vproto.push_back("A2B4C_mP14_11_2e_4e_e;3;14;11;5;18;mP14;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;Au2In4Yb1;Au2In4Yb1 (ICSD #261042)");
    vproto.push_back("A2BC2_mC10_12_i_a_i;3;5;12;5;8;mC10;a,b/a,c/a,beta,x2,z2,x3,z3;-;Cu2Eu1Sn2;Cu2Eu1Sn2 (ICSD #182050)");
    vproto.push_back("ABC_mC12_12_i_i_i;3;6;12;5;10;mC12;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3;-;Al1Ce1Co1;Al1Ce1Co1 (ICSD #20632)");
    vproto.push_back("A3B2C2_mC14_12_ai_i_i;3;7;12;5;10;mC14;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4;-;In3Rh2Sr2;In3Rh2Sr2 (ICSD #410985)");
    vproto.push_back("A6B2C_mC18_12_3i_i_a;3;9;12;5;12;mC18;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5;-;Na6Sn2Zn1;Na6Sn2Zn1 (ICSD #260159)");
    vproto.push_back("A2B5C4_mC22_12_i_a2i_2i;3;11;12;5;14;mC22;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;Mn2Sn5Yb4;Mn2Sn5Yb4 (ICSD #419134)");
    vproto.push_back("A3B4C8_mC30_12_ai_hi_2ij;3;15;12;5;16;mC30;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;-;Ba3Li4Sn8;Ba3Li4Sn8 (ICSD #240016)");
    vproto.push_back("A6B2C7_mC30_12_3i_i_a3i;3;15;12;5;18;mC30;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;Ca6Cu2Sn7;Ca6Cu2Sn7 (ICSD #171243)");
    vproto.push_back("AB2C2_mC20_15_e_f_f;3;10;15;5;11;mC20;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3;-;In1Pd2Sr2;In1Pd2Sr2 (ICSD #391432)");
    vproto.push_back("A2B4C_mC28_15_f_2f_e;3;14;15;5;14;mC28;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Na2Sn4Sr1;Na2Sn4Sr1 (ICSD #261117)");
    vproto.push_back("ABC2_oC16_20_b_a_c;3;8;20;5;8;oC16;a,b/a,c/a,x1,y2,x3,y3,z3;-;Ag1Bi1K2;Ag1Bi1K2 (ICSD #1156)");
    vproto.push_back("ABC_oP12_33_a_a_a;3;12;33;5;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;Ce1Pd1Sn1;Ce1Pd1Sn1 (ICSD #106418)");
    vproto.push_back("ABC_oC12_36_a_a_a;3;6;36;5;9;oC12;a,b/a,c/a,y1,z1,y2,z2,y3,z3;-;Ag1Ce1Sn1;Ag1Ce1Sn1 (ICSD #55819)");
    vproto.push_back("ABC2_oC32_36_b_2a_2ab;3;16;36;5;17;oC32;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5,x6,y6,z6;-;Bi1Pb1Pd2;Bi1Pb1Pd2 (ICSD #56278)");
    vproto.push_back("ABC_oC18_38_ad_be_bd;3;9;38;5;12;oC18;a,b/a,c/a,z1,z2,z3,y4,z4,y5,z5,y6,z6;-;In1Nd1Rh1;In1Nd1Rh1 (ICSD #59430)");
    vproto.push_back("ABC_oC18_38_bd_be_ad;3;9;38;5;12;oC18;a,b/a,c/a,z1,z2,z3,y4,z4,y5,z5,y6,z6;-;Cu1Dy1In1;Cu1Dy1In1 (ICSD #260992)");
    vproto.push_back("ABC_oC18_38_be_ad_bd;3;9;38;5;12;oC18;a,b/a,c/a,z1,z2,z3,y4,z4,y5,z5,y6,z6;-;Gd1In1Pd1;Gd1In1Pd1 (ICSD #415506)");
    vproto.push_back("A10BC3_oP14_47_2q2rs_h_bt;3;14;47;5;9;oP14;a,b/a,c/a,z3,z4,z5,z6,z7,z8;-;Ga10Ni1Pr3;Ga10Ni1Pr3 (ICSD #20663)");
    vproto.push_back("ABC2_oP8_51_f_e_be;3;8;51;5;6;oP8;a,b/a,c/a,z2,z3,z4;-;In1La1Ni2;In1La1Ni2 (ICSD #106800)");
    vproto.push_back("A2BC_oP8_51_ae_f_f;3;8;51;5;6;oP8;a,b/a,c/a,z2,z3,z4;-;Co2Ga1La1;Co2Ga1La1 (ICSD #623100)");
    vproto.push_back("A2BC_oP8_51_af_f_e;3;8;51;5;6;oP8;a,b/a,c/a,z2,z3,z4;-;Co2Dy1In1;Co2Dy1In1 (ICSD #658675)");
    vproto.push_back("A4BC_oP12_51_afj_e_e;3;12;51;5;8;oP12;a,b/a,c/a,z2,z3,z4,x5,z5;-/-/-;Al4Co1La1/Al4Ce1Co1/In4Rh1Yb1;Al4Co1La1 (ICSD #9986)/Al4Ce1Co1 (ICSD #55598)/In4Rh1Yb1 (ICSD #411393)");
    vproto.push_back("A3B2C_oP12_51_ak_ef_f;3;12;51;5;8;oP12;a,b/a,c/a,z2,z3,z4,y5,z5;-;Au3Rb2Tl1;Au3Rb2Tl1 (ICSD #249924)");
    vproto.push_back("AB4C_oP12_51_e_afj_e;3;12;51;5;8;oP12;a,b/a,c/a,z2,z3,z4,x5,z5;-;Ca1In4Rh1;Ca1In4Rh1 (ICSD #410891)");
    vproto.push_back("ABC_oP12_57_c_d_d;3;12;57;5;8;oP12;a,b/a,c/a,x1,x2,y2,x3,y3;-;Al1Ca1Pd1;Al1Ca1Pd1 (ICSD #370036)");
    vproto.push_back("AB2C_oP16_57_d_2d_c;3;16;57;5;10;oP16;a,b/a,c/a,x1,x2,y2,x3,y3,x4,y4;-;Bi1K2Sn1;Bi1K2Sn1 (ICSD #107616)");
    vproto.push_back("AB2C_oP8_59_a_e_b;3;8;59;5;7;oP8;a,b/a,c/a,z1,z2,y3,z3;-;Cu1Ni2Ti1;Cu1Ni2Ti1 (ICSD #628580)");
    vproto.push_back("A3B3C_oP14_59_ae_ae_b;3;14;59;5;10;oP14;a,b/a,c/a,z1,z2,z3,y4,z4,y5,z5;-;Au3In3Sr1;Au3In3Sr1 (ICSD #245679)");
    vproto.push_back("AB3C4_oP16_59_b_ae_ef;3;16;59;5;11;oP16;a,b/a,c/a,z1,z2,y3,z3,y4,z4,x5,z5;-/-;Au1K3Sn4/Au1Cs3Pb4;Au1K3Sn4 (ICSD #107444)/Au1Cs3Pb4 (ICSD #107448)");
    vproto.push_back("AB2C_oP16_62_c_d_c;3;16;62;5;10;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3;-/-;Ca1Cd2Pd1/Ga1Pd2Y1;Ca1Cd2Pd1 (ICSD #425509)/Ga1Pd2Y1 (ICSD #602790)");
    vproto.push_back("ABC2_oP16_62_c_c_d;3;16;62;5;10;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3;-;Dy1Ga1Pd2;Dy1Ga1Pd2 (ICSD #656966)");
    vproto.push_back("ABC_oC12_63_c_c_c;3;6;63;5;6;oC12;a,b/a,c/a,y1,y2,y3;-;Mg1Ni1Tb1;Mg1Ni1Tb1 (ICSD #166874)");
    vproto.push_back("ABC_oC12_63_c_c_a;3;6;63;5;5;oC12;a,b/a,c/a,y2,y3;-;Pt1Sc1Zn1;Pt1Sc1Zn1 (ICSD #424600)");
    vproto.push_back("A2BC_oC16_63_f_c_c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,z3;-/-/-/-/-/-/-;Al2Co1Y1/Al2Cu1Mg1/Al2Dy1Ni1/Ga2Ni1Sc1/In2Rh1Yb1/Mg2Pd1Sr1/Al2Ce1Pt1;Al2Co1Y1 (ICSD #57645)/Al2Cu1Mg1 (ICSD #57693)/Al2Dy1Ni1 (ICSD #57759)/Ga2Ni1Sc1 (ICSD #103875)/In2Rh1Yb1 (ICSD #152048)/Mg2Pd1Sr1 (ICSD #425493)/Al2Ce1Pt1 (ICSD #658140)");
    vproto.push_back("ABC2_oC16_63_a_c_g;3;8;63;5;6;oC16;a,b/a,c/a,y2,x3,y3;-/-/-/-;Au1Bi1Li2/Au1Bi1Li2/Au1Bi1Na2/Au1Bi1K2;Au1Bi1Li2 (ICSD #261785)/Au1Bi1Li2 (ICSD #261786)/Au1Bi1Na2 (ICSD #261788)/Au1Bi1K2 (ICSD #380341)");
    vproto.push_back("ABC2_oC16_63_c_c_2c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,y4;-/-/-/-/-;Ba1Cu1Sn2/Co1Ho1Sn2/Ce1Ni1Sn2/Cu1Pr1Sn2/La1Rh1Sn2;Ba1Cu1Sn2 (ICSD #58647)/Co1Ho1Sn2 (ICSD #240098)/Ce1Ni1Sn2 (ICSD #621687)/Cu1Pr1Sn2 (ICSD #628742)/La1Rh1Sn2 (ICSD #641754)");
    vproto.push_back("AB2C_oC16_63_c_f_c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,z3;-/-/-/-;La1Mg2Ni1/Ce1In2Pd1/Ba1In2Ir1/Dy1Ga2Pd1;La1Mg2Ni1 (ICSD #96152)/Ce1In2Pd1 (ICSD #150132)/Ba1In2Ir1 (ICSD #414434)/Dy1Ga2Pd1 (ICSD #600150)");
    vproto.push_back("A2BC_oC16_63_g_c_c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,x3,y3;-;Cd2Cu1Er1;Cd2Cu1Er1 (ICSD #99139)");
    vproto.push_back("ABC2_oC16_63_c_c_f;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,z3;-;Eu1Pd1Tl2;Eu1Pd1Tl2 (ICSD #172230)");
    vproto.push_back("A2B2C_oC20_63_g_2c_a;3;10;63;5;7;oC20;a,b/a,c/a,y2,y3,x4,y4;-;Bi2Cs2Pt1;Bi2Cs2Pt1 (ICSD #658701)");
    vproto.push_back("A4BC_oC24_63_acf_c_c;3;12;63;5;8;oC24;a,b/a,c/a,y2,y3,y4,y5,z5;-/-;Al4Dy1Ni1/Al4Ni1Tb1;Al4Dy1Ni1 (ICSD #57760)/Al4Ni1Tb1 (ICSD #58060)");
    vproto.push_back("AB3C2_oC24_63_c_cg_e;3;12;63;5;8;oC24;a,b/a,c/a,y1,y2,x3,x4,y4;-;Ca1Ga3Ni2;Ca1Ga3Ni2 (ICSD #58898)");
    vproto.push_back("A3BC3_oC28_63_cf_a_cf;3;14;63;5;9;oC28;a,b/a,c/a,y2,y3,y4,z4,y5,z5;-;Co3Ga1Y3;Co3Ga1Y3 (ICSD #10044)");
    vproto.push_back("ABC2_oC32_63_ac_f_2cf;3;16;63;5;10;oC32;a,b/a,c/a,y2,y3,y4,y5,z5,y6,z6;-;La1Rh1Sn2;La1Rh1Sn2 (ICSD #410732)");
    vproto.push_back("A2B2C_oC10_65_j_i_a;3;5;65;5;5;oC10;a,b/a,c/a,y2,y3;-;Ho2Ni2Pb1;Ho2Ni2Pb1 (ICSD #54612)");
    vproto.push_back("A2B3C_oC12_65_i_cf_a;3;6;65;5;4;oC12;a,b/a,c/a,y4;-;Al2Ni3Pr1;Al2Ni3Pr1 (ICSD #107864)");
    vproto.push_back("A3B2C2_oC14_65_aj_j_i;3;7;65;5;6;oC14;a,b/a,c/a,y2,y3,y4;-;Mg3Ni2Tb2;Mg3Ni2Tb2 (ICSD #240761)");
    vproto.push_back("A2BC_oC16_65_aci_j_i;3;8;65;5;6;oC16;a,b/a,c/a,y3,y4,y5;-;Ga2Nd1Ni1;Ga2Nd1Ni1 (ICSD #103850)");
    vproto.push_back("A2BC6_oC18_65_j_a_2ij;3;9;65;5;7;oC18;a,b/a,c/a,y2,y3,y4,y5;-;Dy2Ni1Sn6;Dy2Ni1Sn6 (ICSD #165783)");
    vproto.push_back("A3B2C7_oC24_65_aj_i_c2ij;3;12;65;5;8;oC24;a,b/a,c/a,y3,y4,y5,y6,y7;-/-;Ce3Ni2Sn7La3Ni2Sn7;Ce3Ni2Sn7 (ICSD #102239)/La3Ni2Sn7 (ICSD #160913)");
    vproto.push_back("A2B2C_oI10_71_f_h_a;3;5;71;5;5;oI10;a,b/a,c/a,x2,y3;-/-;Ce2Ni2Sn1/Ca2Cu2Ga1;Ce2Ni2Sn1 (ICSD #55495)/Ca2Cu2Ga1 (ICSD #58885)");
    vproto.push_back("A2B2C3_oI14_71_e_f_af;3;7;71;5;6;oI14;a,b/a,c/a,x2,x3,x4;-;Al2Sn2Sr3;Al2Sn2Sr3 (ICSD #9564)");
    vproto.push_back("A5B2C_oI16_71_an_e_d;3;8;71;5;6;oI16;a,b/a,c/a,x3,x4,y4;-;Al5Ni2Pr1;Al5Ni2Pr1 (ICSD #58046)");
    vproto.push_back("A4B3C4_oI22_71_n_af_eg;3;11;71;5;8;oI22;a,b/a,c/a,x2,x3,y4,x5,y5;-/-;Ag8Ce6Sn8/Ag8Pr6Sn8;Ag8Ce6Sn8 (ICSD #55826)/Ag8Pr6Sn8 (ICSD #55827)");
    vproto.push_back("A4B3C4_oI22_71_n_af_eh;3;11;71;5;8;oI22;a,b/a,c/a,x2,x3,y4,x5,y5;-;Ag4Dy3Sn4;Ag4Dy3Sn4 (ICSD #156968)");
    vproto.push_back("A3B4C4_oI22_71_af_eh_n;3;11;71;5;8;oI22;a,b/a,c/a,x2,x3,y4,x5,y5;-;La3Pd4Zn4;La3Pd4Zn4 (ICSD #182774)");
    vproto.push_back("A6B6C_oI26_71_hk_efg_a;3;13;71;5;7;oI26;a,b/a,c/a,x2,x3,y4,y5;-;Fe6Ga6Sc1;Fe6Ga6Sc1 (ICSD #103465)");
    vproto.push_back("A3B3C8_oI28_71_cg_be_2n;3;14;71;5;9;oI28;a,b/a,c/a,x3,y4,x5,y5,x6,y6;-/-;Ag3Er3Ga8/Au3Er3Ga8;Ag3Er3Ga8 (ICSD #605107)/Au3Er3Ga8 (ICSD #611824)");
    vproto.push_back("A3B10C_oI28_71_be_g2n_c;3;14;71;5;9;oI28;a,b/a,c/a,x3,y4,x5,y5,x6,y6;-;Ba3Hg10In1;Ba3Hg10In1 (ICSD #290305)");
    vproto.push_back("A6B5C3_oI28_71_gn_cn_be;3;14;71;5;9;oI28;a,b/a,c/a,x3,y4,x5,y5,x6,y6;-;Cu6Sn5Yb3;Cu6Sn5Yb3 (ICSD #413485)");
    vproto.push_back("A9B3C2_oI28_71_c2n_be_g;3;14;71;5;9;oI28;a,b/a,c/a,x3,y4,x5,y5,x6,y6;-;Ga9Ho3Pd2;Ga9Ho3Pd2 (ICSD #634397)");
    vproto.push_back("A8B3C3_oI28_71_2n_be_cg;3;14;71;5;9;oI28;a,b/a,c/a,x3,y4,x5,y5,x6,y6;-;Ga8Ho3Ir3;Ga8Ho3Ir3 (ICSD #634379)");
    vproto.push_back("A2B2C_oI20_72_j_j_a;3;10;72;5;7;oI20;a,b/a,c/a,x2,y2,x3,y3;-;Ba2Bi2Zn1;Ba2Bi2Zn1 (ICSD #421424)");
    vproto.push_back("ABC3_tI10_107_a_a_ab;3;5;107;5;6;tI10;a,c/a,z1,z2,z3,z4;-;Ba1Ni1Sn3;Ba1Ni1Sn3 (ICSD #58662)");
    vproto.push_back("A3BC_tI10_107_ab_a_a;3;5;107;5;6;tI10;a,c/a,z1,z2,z3,z4;-;Al3Cu1Pr1;Al3Cu1Pr1 (ICSD #290390)");
    vproto.push_back("A3BC_tI10_119_bf_a_c;3;5;119;5;3;tI10;a,c/a,z4;-;In3Mg1Sr1;In3Mg1Sr1 (ICSD #249592)");
    vproto.push_back("AB4C2_tI28_120_d_i_e;3;14;120;5;6;tI28;a,c/a,x2,x3,y3,z3;-;Ce1Ni4Sn2;Ce1Ni4Sn2 (ICSD #102238)");
    vproto.push_back("A4B6C_tI22_121_i_ci_a;3;11;121;5;6;tI22;a,c/a,x3,z3,x4,z4;-;Ru4Sn6Y1;Ru4Sn6Y1 (ICSD #54354)");
    vproto.push_back("AB2C_tP4_123_b_h_a;3;4;123;5;3;tP4;a,c/a,z3;-;Cd1Pt2Zn1;Cd1Pt2Zn1 (ICSD #102057)");
    vproto.push_back("AB4C_tP6_123_b_i_a;3;6;123;5;3;tP6;a,c/a,z3;-;Co1Ga4Hf1;Co1Ga4Hf1 (ICSD #623081)");
    vproto.push_back("ABC5_tP7_123_a_b_ci;3;7;123;5;3;tP7;a,c/a,z4;-;Ce1Co1In5;Ce1Co1In5 (ICSD #102108)");
    vproto.push_back("AB5C_tP7_123_a_ci_b;3;7;123;5;3;tP7;a,c/a,z4;-;Co1In5Tb1;Co1In5Tb1 (ICSD #623944)");
    vproto.push_back("AB6C_tP8_123_a_hi_b;3;8;123;5;4;tP8;a,c/a,z3,z4;-;Ce1Ga6Pd1;Ce1Ga6Pd1 (ICSD #240161)");
    vproto.push_back("AB8C2_tP11_123_a_ehi_g;3;11;123;5;5;tP11;a,c/a,z3,z4,z5;-/-;Co1Ga8Ho2/Fe1Ga8Ho2;Co1Ga8Ho2 (ICSD #42426)/Fe1Ga8Ho2 (ICSD #180132)");
    vproto.push_back("A2BC2_tP10_127_g_a_h;3;10;127;5;4;tP10;a,c/a,x2,x3;-;Ni2Sn1Zr2;Ni2Sn1Zr2 (ICSD #54303)");
    vproto.push_back("AB2C2_tP10_127_a_g_h;3;10;127;5;4;tP10;a,c/a,x2,x3;-;Al1Cu2Re2;Al1Cu2Re2 (ICSD #57706)");
    vproto.push_back("ABC_tP6_129_b_c_c;3;6;129;5;4;tP6;a,c/a,z2,z3;-;Mn1Sn1Sr1;Mn1Sn1Sr1 (ICSD #66951)");
    vproto.push_back("ABC_tP6_129_c_c_a;3;6;129;5;4;tP6;a,c/a,z2,z3;-;Bi1K1Mn1;Bi1K1Mn1 (ICSD #601586)");
    vproto.push_back("ABC_tP6_129_a_c_c;3;6;129;5;4;tP6;a,c/a,z2,z3;-;Ba1Mg1Pb1;Ba1Mg1Pb1 (ICSD #615958)");
    vproto.push_back("A2BC_tP8_129_bc_c_a;3;8;129;5;4;tP8;a,c/a,z3,z4;-/-;Sn2Tb1Zn1/Bi2La1Li1;Sn2Tb1Zn1 (ICSD #163425)/Bi2La1Li1 (ICSD #415728)");
    vproto.push_back("AB2C2_tP10_129_c_ac_bc;3;10;129;5;5;tP10;a,c/a,z3,z4,z5;-;Nd1Ni2Sn2;Nd1Ni2Sn2 (ICSD #160053)");
    vproto.push_back("AB2C2_tP10_129_c_bc_ac;3;10;129;5;5;tP10;a,c/a,z3,z4,z5;-;Eu1Sn2Zn2;Eu1Sn2Zn2 (ICSD #162267)");
    vproto.push_back("A2BC2_tP10_129_bc_c_ac;3;10;129;5;5;tP10;a,c/a,z3,z4,z5;-;Bi2Ce1Ni2;Bi2Ce1Ni2 (ICSD #616562)");
    vproto.push_back("A2BC4_tP14_129_ac_b_j;3;14;129;5;5;tP14;a,c/a,z3,x4,z4;-;Ce2Ru1Zn4;Ce2Ru1Zn4 (ICSD #418547)");
    vproto.push_back("AB2C_tI8_139_b_d_a;3;4;139;5;2;tI8;a,c/a;-;Li1Pd2Tl1;Li1Pd2Tl1 (ICSD #54365)");
    vproto.push_back("A2BC2_tI10_139_d_a_e;3;5;139;5;3;tI10;a,c/a,z3;-/-;Al2Ce1Ga2/Al2Ca1Zn2;Al2Ce1Ga2 (ICSD #55789)/Al2Ca1Zn2 (ICSD #57550)");
    vproto.push_back("ABC_tI12_139_c_e_e;3;6;139;5;4;tI12;a,c/a,z2,z3;-;Mg1Pr1Sn1;Mg1Pr1Sn1 (ICSD #182479)");
    vproto.push_back("ABC_tI12_139_e_c_e;3;6;139;5;4;tI12;a,c/a,z2,z3;-/-;Dy1Mg1Sn1/Au1Cu1Sn1;Dy1Mg1Sn1 (ICSD #183524)/Au1Cu1Sn1 (ICSD #611768)");
    vproto.push_back("A4B2C_tI14_139_h_d_a;3;7;139;5;3;tI14;a,c/a,x3;-;Al4Mo2Yb1;Al4Mo2Yb1 (ICSD #456)");
    vproto.push_back("AB4C2_tI14_139_a_h_d;3;7;139;5;3;tI14;a,c/a,x3;-;Er1Ga4Ti2;Er1Ga4Ti2 (ICSD #630579)");
    vproto.push_back("A5B2C_tI16_139_bg_e_a;3;8;139;5;4;tI16;a,c/a,z3,z4;-;Al5Ni2Zr1;Al5Ni2Zr1 (ICSD #58084)");
    vproto.push_back("A2BC_tI16_139_ce_e_d;3;8;139;5;4;tI16;a,c/a,z3,z4;-;Bi2Sr1Zn1;Bi2Sr1Zn1 (ICSD #41924)");
    vproto.push_back("AB2C_tI16_139_e_ce_d;3;8;139;5;4;tI16;a,c/a,z3,z4;-/-;Ba1Bi2Cd1/Ba1Bi2Zn1;Ba1Bi2Cd1 (ICSD #58635)/Ba1Bi2Zn1 (ICSD #58638)");
    vproto.push_back("A2BC_tI16_139_ce_d_e;3;8;139;5;4;tI16;a,c/a,z3,z4;-;Bi2Cd1Sr1;Bi2Cd1Sr1 (ICSD #58764)");
    vproto.push_back("A2B5C3_tI20_139_e_bg_ae;3;10;139;5;5;tI20;a,c/a,z3,z4,z5;-;K2Mg5Sn3;K2Mg5Sn3 (ICSD #421342)");
    vproto.push_back("A8BC4_tI26_139_ij_a_f;3;13;139;5;4;tI26;a,c/a,x3,x4;-;Al8Ca1Cu4;Al8Ca1Cu4 (ICSD #57539)");
    vproto.push_back("A4B8C_tI26_139_i_fj_a;3;13;139;5;4;tI26;a,c/a,x3,x4;-;Cr4Fe8Y1;Cr4Fe8Y1 (ICSD #168240)");
    vproto.push_back("A8B4C_tI26_139_ij_f_a;3;13;139;5;4;tI26;a,c/a,x3,x4;-;Al8Mn4Y1;Al8Mn4Y1 (ICSD #57997)");
    vproto.push_back("ABC4_tI24_140_a_c_l;3;12;140;5;4;tI24;a,c/a,x3,z3;-;Ir1Li1Sn4;Ir1Li1Sn4 (ICSD #172149)");
    vproto.push_back("A2BC4_tI28_140_h_a_k;3;14;140;5;5;tI28;a,c/a,x2,x3,y3;-;Bi2Mn1Ti4;Bi2Mn1Ti4 (ICSD #150145)");
    vproto.push_back("A4B2C_tI28_140_l_h_a;3;14;140;5;5;tI28;a,c/a,x2,x3,z3;-;Cu4Sn2Sr1;Cu4Sn2Sr1 (ICSD #182104)");
    vproto.push_back("AB5C2_tI32_140_a_cl_h;3;16;140;5;5;tI32;a,c/a,x3,x4,z4;-;Bi1Er5Pt2;Bi1Er5Pt2 (ICSD #107217)");
    vproto.push_back("A3B4C_tI32_140_ah_l_c;3;16;140;5;5;tI32;a,c/a,x3,x4,z4;-;Bi3In4Pb1;Bi3In4Pb1 (ICSD #616725)");
    vproto.push_back("AB5C2_tI32_140_a_bk_h;3;16;140;5;5;tI32;a,c/a,x3,x4,y4;-;Ga1Nb5Sn2;Ga1Nb5Sn2 (ICSD #103842)");
    vproto.push_back("A2BC5_tI32_140_h_a_cl;3;16;140;5;5;tI32;a,c/a,x3,x4,z4;-;Au2Bi1Tb5;Au2Bi1Tb5 (ICSD #156956)");
    vproto.push_back("AB2C5_tI32_140_a_h_cl;3;16;140;5;5;tI32;a,c/a,x3,x4,z4;-;Bi1Co2Ho5;Bi1Co2Ho5 (ICSD #161658)");
    vproto.push_back("A4B7C3_hR14_146_ab_a2b_b;3;14;146;5;16;hR14;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Ga4Pd7Zn3;Ga4Pd7Zn3 (ICSD #103911)");
    vproto.push_back("ABC_hP9_156_3a_b2c_2bc;3;9;156;5;11;hP9;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8,z9;-;Ca1Li1Sn1;Ca1Li1Sn1 (ICSD #58911)");
    vproto.push_back("A5BC4_hR10_160_5a_a_4a;3;10;160;5;12;hR10;a,c/a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10;-;Li5Na1Sn4;Li5Na1Sn4 (ICSD #12142)");
    vproto.push_back("A2BC_hR4_166_c_b_a;3;4;166;5;3;hR4;a,c/a,x3;-;Cu2Ga1Sr1;Cu2Ga1Sr1 (ICSD #102938)");
    vproto.push_back("AB2C4_hR7_166_a_c_2c;3;7;166;5;5;hR7;a,c/a,x2,x3,x4;-;Au1Ni2Sn4;Au1Ni2Sn4 (ICSD #150127)");
    vproto.push_back("AB4C6_hR11_166_a_2c_h;3;11;166;5;6;hR11;a,c/a,x2,x3,x4,z4;-;Gd1Sn4Ti6;Gd1Sn4Ti6 (ICSD #166193)");
    vproto.push_back("A2BC_hR12_166_h_bc_ac;3;12;166;5;6;hR12;a,c/a,x3,x4,x5,z5;-;Al2Cu1Yb1;Al2Cu1Yb1 (ICSD #604213)");
    vproto.push_back("A7B3C2_hR12_166_bh_ac_c;3;12;166;5;6;hR12;a,c/a,x3,x4,x5,z5;-;Al7Ca3Cu2;Al7Ca3Cu2 (ICSD #57538)");
    vproto.push_back("AB2C9_hR12_166_a_c_bch;3;12;166;5;6;hR12;a,c/a,x3,x4,x5,z5;-;La1Mg2Ni9;La1Mg2Ni9 (ICSD #55614)");
    vproto.push_back("A9B2C_hR12_166_eh_c_a;3;12;166;5;5;hR12;a,c/a,x2,x4,z4;-;Ni3Ti0.67Zr0.33;Ni3Ti0.67Zr0.33 (ICSD #646978)");
    vproto.push_back("AB2C_hR12_166_bc_h_ac;3;12;166;5;6;hR12;a,c/a,x3,x4,x5,z5;-;Ag1Al2Pr1;Ag1Al2Pr1 (ICSD #604688)");
    vproto.push_back("A2B7C3_hR12_166_c_bh_ac;3;12;166;5;6;hR12;a,c/a,x3,x4,x5,z5;-;Ag2Al7Ca3;Ag2Al7Ca3 (ICSD #104173)");
    vproto.push_back("A7B4C2_hR13_166_ah_2c_c;3;13;166;5;7;hR13;a,c/a,x2,x3,x4,x5,z5;-;Au7Rb4Sn2;Au7Rb4Sn2 (ICSD #58581)");
    vproto.push_back("A3B2C2_hR14_166_3c_2c_abc;3;14;166;5;8;hR14;a,c/a,x3,x4,x5,x6,x7,x8;-;Au3Sn2Yb2;Au3Sn2Yb2 (ICSD #710044)");
    vproto.push_back("ABC_hP6_186_b_a_b;3;6;186;5;5;hP6;a,c/a,z1,z2,z3;-/-;Au1Pr1Sn1/Bi1Li1Zn1;Au1Pr1Sn1 (ICSD #54997)/Bi1Li1Zn1 (ICSD #100115)");
    vproto.push_back("ABC_hP6_186_a_b_b;3;6;186;5;5;hP6;a,c/a,z1,z2,z3;-;Ca1Sn1Zn1;Ca1Sn1Zn1 (ICSD #106358)");
    vproto.push_back("A4BC_hP12_186_ac_b_b;3;12;186;5;7;hP12;a,c/a,z1,z2,z3,x4,z4;-;Cu4In1Mn1;Cu4In1Mn1 (ICSD #424277)");
    vproto.push_back("A3BC4_hP16_186_a2b_b_2a2b;3;16;186;5;10;hP16;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8;-;Au3Li1Sn4;Au3Li1Sn4 (ICSD #412207)");
    vproto.push_back("ABC_hP6_187_i_h_ab;3;6;187;5;4;hP6;a,c/a,z3,z4;-;Au1Ga1Zr1;Au1Ga1Zr1 (ICSD #156264)");
    vproto.push_back("A4BC_hP6_187_hi_b_a;3;6;187;5;4;hP6;a,c/a,z3,z4;-;Ga4Li1Y1;Ga4Li1Y1 (ICSD #98666)");
    vproto.push_back("ABC_hP9_187_dh_ai_eg;3;9;187;5;5;hP9;a,c/a,z4,z5,z6;-;Ca1Li1Pb1;Ca1Li1Pb1 (ICSD #409533)");
    vproto.push_back("ABC_hP9_187_eh_ci_bg;3;9;187;5;5;hP9;a,c/a,z4,z5,z6;-;Li1Pb1Yb1;Li1Pb1Yb1 (ICSD #409534)");
    vproto.push_back("A5B5C_hP11_187_hk_cgi_a;3;11;187;5;6;hP11;a,c/a,z3,z4,z5,x6;-;Al5Ba5Pb1;Al5Ba5Pb1 (ICSD #416337)");
    vproto.push_back("A4B6C_hP11_187_ck_jk_a;3;11;187;5;5;hP11;a,c/a,x3,x4,x5;-;Au4In6K1;Au4In6K1 (ICSD #249520)");
    vproto.push_back("A3B2C3_hP8_189_f_c_g;3;8;189;5;4;hP8;a,c/a,x2,x3;-;In3Rh2Ti3;In3Rh2Ti3 (ICSD #410967)");
    vproto.push_back("ABC_hP9_189_g_f_ad;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Bi1Dy1Rh1;Bi1Dy1Rh1 (ICSD #51845)");
    vproto.push_back("AB2C6_hP9_189_b_c_fg;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Co1Ga2Zr6;Co1Ga2Zr6 (ICSD #20876)");
    vproto.push_back("A5B3C_hP9_189_cf_g_b;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Ga5Lu3Ni1;Ga5Lu3Ni1 (ICSD #634575)");
    vproto.push_back("ABC_hP9_189_bc_g_f;3;9;189;5;4;hP9;a,c/a,x3,x4;-/-;Ni1Pr1Zn1/Pd1Tb1Tl1;Ni1Pr1Zn1 (ICSD #54295)/Pd1Tb1Tl1 (ICSD #105713)");
    vproto.push_back("A2BC6_hP9_189_c_b_fg;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Bi2Fe1Ho6;Bi2Fe1Ho6 (ICSD #96253)");
    vproto.push_back("ABC_hP9_189_g_bc_f;3;9;189;5;4;hP9;a,c/a,x3,x4;-/-;Pb1Pd1Y1/Hf1Ir1Sn1;Pb1Pd1Y1 (ICSD #54315)/Hf1Ir1Sn1 (ICSD #638586)");
    vproto.push_back("ABC_hP9_189_f_bc_g;3;9;189;5;4;hP9;a,c/a,x3,x4;-/-;Mg1Sn1Yb1/Ga1Ni1Ti1;Mg1Sn1Yb1 (ICSD #54344)/Ga1Ni1Ti1 (ICSD #103885)");
    vproto.push_back("ABC_hP9_189_f_ad_g;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Dy1In1Mg1;Dy1In1Mg1 (ICSD #55102)");
    vproto.push_back("ABC_hP9_189_g_f_bc;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Dy1In1Ni1;Dy1In1Ni1 (ICSD #55846)");
    vproto.push_back("ABC_hP9_189_f_g_bc;3;9;189;5;4;hP9;a,c/a,x3,x4;-/-/-;Al1Dy1Ni1/Al1Hf1Pt1/Al1Nd1Ni1;Al1Dy1Ni1 (ICSD #107416)/Al1Hf1Pt1 (ICSD #608141)/Al1Nd1Ni1 (ICSD #608770)");
    vproto.push_back("ABC_hP9_189_ad_g_f;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Al1Mg1Y1;Al1Mg1Y1 (ICSD #160908)");
    vproto.push_back("ABC_hP9_189_bc_f_g;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Cu1In1Tb1;Cu1In1Tb1 (ICSD #628137)");
    vproto.push_back("A3B2C_hP6_191_g_c_a;3;6;191;5;2;hP6;a,c/a;-/-;Ag3Al2La1/Co3Ga2Tb1;Ag3Al2La1 (ICSD #57329)/Co3Ga2Tb1 (ICSD #102452)");
    vproto.push_back("A2BC3_hP6_191_c_a_g;3;6;191;5;2;hP6;a,c/a;-/-;Cu2Dy1In3/Al2Ce1Pt3;Cu2Dy1In3 (ICSD #108369)/Al2Ce1Pt3 (ICSD #658142)");
    vproto.push_back("A9BC2_hP12_191_fm_a_c;3;12;191;5;3;hP12;a,c/a,x4;-;Al9Ba1Fe2;Al9Ba1Fe2 (ICSD #57518)");
    vproto.push_back("A6BC6_hP13_191_i_a_cde;3;13;191;5;4;hP13;a,c/a,z4,z5;-;Mn6Sc1Sn6;Mn6Sc1Sn6 (ICSD #54273)");
    vproto.push_back("A6B6C_hP13_191_i_cde_a;3;13;191;5;4;hP13;a,c/a,z4,z5;-;Mn6Sn6Tb1;Mn6Sn6Tb1 (ICSD #54277)");
    vproto.push_back("A3B3C2_hP16_193_g_g_d;3;16;193;5;4;hP16;a,c/a,x2,x3;-;Ga3Hf3Nb2;Ga3Hf3Nb2 (ICSD #103733)");
    vproto.push_back("ABC_hP6_194_d_a_c;3;6;194;5;2;hP6;a,c/a;-;Pb1Sr1Zn1;Pb1Sr1Zn1 (ICSD #54319)");
    vproto.push_back("ABC_hP6_194_a_c_d;3;6;194;5;2;hP6;a,c/a;-/-/-;Ag1Ba1Bi1/Ce1Cu1Sn1/Ce1Tl1Zn1;Ag1Ba1Bi1 (ICSD #56978)/Ce1Cu1Sn1 (ICSD #106384)/Ce1Tl1Zn1 (ICSD #622301)");
    vproto.push_back("ABC_hP6_194_c_a_d;3;6;194;5;2;hP6;a,c/a;-/-;Bi1Ca1Cu1/Cu1Gd1Sn1;Bi1Ca1Cu1 (ICSD #57018)/Cu1Gd1Sn1 (ICSD #600994)");
    vproto.push_back("ABC_hP6_194_a_d_c;3;6;194;5;2;hP6;a,c/a;-/-;Ba1Bi1Cu1/Ca1Hg1Pb1;Ba1Bi1Cu1 (ICSD #106303)/Ca1Hg1Pb1 (ICSD #106345)");
    vproto.push_back("AB2C_hP8_194_b_f_c;3;8;194;5;3;hP8;a,c/a,z3;-;Al1Pt2Zr1;Al1Pt2Zr1 (ICSD #58139)");
    vproto.push_back("AB2C_hP8_194_a_f_c;3;8;194;5;3;hP8;a,c/a,z3;-;In1Pt2Y1;In1Pt2Y1 (ICSD #59505)");
    vproto.push_back("AB2C_hP8_194_c_f_a;3;8;194;5;3;hP8;a,c/a,z3;-/-;Er1Pt2Sn1/Ce1Cu2Mg1;Er1Pt2Sn1 (ICSD #103292)/Ce1Cu2Mg1 (ICSD #157110)");
    vproto.push_back("A2BC_hP8_194_f_a_c;3;8;194;5;3;hP8;a,c/a,z3;-;Cu2Li1Sn1;Cu2Li1Sn1 (ICSD #150602)");
    vproto.push_back("A3BC2_hP12_194_h_a_f;3;12;194;5;4;hP12;a,c/a,z2,x3;-;Al3Ru1Sc2;Al3Ru1Sc2 (ICSD #58159)");
    vproto.push_back("ABC_hP12_194_f_f_ab;3;12;194;5;4;hP12;a,c/a,z3,z4;-;Ga1Sn1Sr1;Ga1Sn1Sr1 (ICSD #66003)");
    vproto.push_back("ABC_hP12_194_ab_f_f;3;12;194;5;4;hP12;a,c/a,z3,z4;-/-;Gd1Sn1Zn1/La1Sn1Zn1;Gd1Sn1Zn1 (ICSD #152621)/La1Sn1Zn1 (ICSD #152623)");
    vproto.push_back("ABC_hP12_194_f_ab_f;3;12;194;5;4;hP12;a,c/a,z3,z4;-;Sn1Y1Zn1;Sn1Y1Zn1 (ICSD #152622)");
    vproto.push_back("A3B3C2_hP16_194_h_af_f;3;16;194;5;5;hP16;a,c/a,z2,z3,x4;-;Al2Ba2Ga1.33;Al2Ba2Ga1.33 (ICSD #415931)");
    vproto.push_back("AB2C_cP16_215_e_bce_ad;3;16;215;5;3;cP16;a,x5,x6;-;Ir1Li2Mg1;Ir1Li2Mg1 (ICSD #180119)");
    vproto.push_back("A3B4C4_cI22_217_b_c_c;3;11;217;5;3;cI22;a,x2,x3;-;Ca3Ga4Ni4;Ca3Ga4Ni4 (ICSD #58899)");
    vproto.push_back("A3BC_cP15_221_ag_c_d;3;15;221;5;2;cP15;a,x4;-;Ga3Os1Tb1;Ga3Os1Tb1 (ICSD #412144)");
    vproto.push_back("AB3C_cP15_221_d_ag_c;3;15;221;5;2;cP15;a,x4;-;Dy1Ga3Ru1;Dy1Ga3Ru1 (ICSD #629734)");
    // done now produce

    // FROM PROTO LIST
    for(uint i=0;i<vproto.size();i++) {
      vproto_label.push_back("");
      vproto_nspecies.push_back(0);
      vproto_natoms.push_back(0);
      vproto_spacegroup.push_back(0);
      vproto_nunderscores.push_back(0);
      vproto_nparameters.push_back(0);
      vproto_Pearson_symbol.push_back("");
      vproto_params.push_back("");
      vproto_Strukturbericht.push_back("");
      vproto_prototype.push_back("");
      vproto_dialect.push_back("");

      anrl::vproto2tokens(vproto.at(i),
          vproto_label.at(i),vproto_nspecies.at(i),vproto_natoms.at(i),vproto_spacegroup.at(i),
          vproto_nunderscores.at(i),vproto_nparameters.at(i),vproto_Pearson_symbol.at(i),
          vproto_params.at(i),vproto_Strukturbericht.at(i),vproto_prototype.at(i),vproto_dialect.at(i));
    }

    return vproto.size();
  }
}

// ***************************************************************************
namespace anrl { 
  vector<string> getANRLParameters(string _anrl_label, string library, int choice, bool keep_original_lattice_parameter){
    //choice indicates choice of parameter set
    //all: grabs all
    //"": looks for one
    //part1: set from part1
    //part2: set from part2
    //part1:0
    //part1:1
    //part2:0
    //part2:1
    //...

    string function_name = "anrl::getANRLParameters()";
    stringstream message;

    vector<string> tokens;
    vector<string> vparameters;

    string anrl_label = _anrl_label;

    // check if number suffix, i.e., predefined structure, with atomic volume scaling
    string number_id = "";
    if(aurostd::substring2bool(anrl_label,"-")){
      aurostd::string2tokens(_anrl_label,tokens,"-");
      anrl_label = tokens[0];
      number_id = tokens[1];
    }

    // add option to keep original scaling
    //bool keep_original_lattice_parameter = false;

    // A WORD TO THE WISE: When adding new preset parameter values for a given 
    // label, add them to the end.  Otherwise, you will mess up the connection
    // to the previously defined preset labels (e.g., -001, -002, etc.)

    // ---------------------------------------------------------------------------
    // Part 1
    // ---------------------------------------------------------------------------
    if(library=="all" || library == "" || library=="part1"){
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_aP12_1_4a_8a"){
        vparameters.push_back("5.417,1.0,1.0,90.0,90.0,90.0,0.001,0.002,0.003,0.4966,0.0001,0.5036,0.5001,0.502,0.0011,-0.0006,0.5013,0.5038,0.3857,0.3832,0.384,0.1149,0.6114,0.8846,0.8854,0.1157,0.6143,0.6153,0.8865,0.1141,0.6151,0.6132,0.6137,0.8854,0.3818,0.1149,0.1147,0.8856,0.3841,0.3857,0.1161,0.8842");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_aP16_1_4a_4a_8a"){
        vparameters.push_back("6.554,1.00061031431,1.92662496186,100.43475,100.46074,107.53,0.3267,0.582,0.177,0.565,-0.0132,0.4424,0.5217,0.3883,0.6767,-0.0744,0.6254,-0.0574,0.0338,0.0476,0.2599,0.0831,0.6072,0.4974,-0.0131,0.0949,0.7583,0.5449,0.1443,-0.0022,-0.0211,0.5213,0.2073,0.2907,0.5956,-0.0183,-0.0616,0.0602,0.4998,0.5068,-0.0175,0.2448,0.4596,0.0397,0.708,0.5326,0.352,0.4818,0.0,0.0,0.0,-0.078,0.569,0.7448");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_aP6_2_2i_i"){
        vparameters.push_back("4.56,1.54824561404,1.62280701754,80.2,106.96667,98.2,0.557,0.73,0.165,0.82,0.803,0.695,0.397,0.639,0.463");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_aP4_2_aci"){
        vparameters.push_back("3.307,2.24130631993,0.844572119746,89.06,85.15,85.7,0.572,0.259,0.433");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP12_3_bc3e_2e"){
        vparameters.push_back("4.1605,0.992524936907,1.78370388174,101.3752,0.15907,0.73859,0.02399,0.752,0.18927,0.38562,0.71473,0.64074,0.48963,0.20196,0.18802,0.18244,0.0,0.69651,0.38098,0.58564,0.17797");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP4_4_2a"){
        vparameters.push_back("3.104,2.42042525773,1.53350515464,92.71,0.25,0.23,0.48,0.48,0.0,0.02");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mC12_5_3c"){
        vparameters.push_back("7.42,0.578167115903,1.90026954178,92.0,0.05,0.27,0.245,0.63,0.3,0.4,0.245,0.43,0.07");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_mC10_8_ab_a_a"){
        vparameters.push_back("5.72204,0.9978207073,0.722908263486,90.498,0.5515,-0.0994,0.0,0.0,0.523,0.4492,0.288,0.2434,0.3729");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC144_9_24a_12a"){
        vparameters.push_back("18.524,0.270092852516,1.28535953358,105.82,0.5749,0.351,0.8182,0.0707,0.34,0.8476,0.7315,0.138,0.4851,0.2509,0.144,0.5152,0.4155,0.352,0.6741,-0.0873,0.352,0.6434,0.8773,0.164,-0.0787,0.416,0.168,-0.0639,0.7741,0.145,0.7538,0.2336,0.143,0.7402,0.6195,0.341,0.5847,0.0811,0.343,0.5661,-0.0034,0.011,0.6062,0.3533,0.489,0.5665,0.6498,0.005,0.6711,0.1524,0.496,0.7805,0.8636,0.499,0.7328,0.3361,0.003,0.8333,0.0052,0.493,0.7398,0.1369,0.011,-0.0732,0.4927,0.492,0.8868,0.5,0.468,0.5,0.2252,0.491,0.5898,0.2744,0.021,-0.0845,0.0507,0.041,0.5642,0.2036,0.447,0.7347,-0.0802,0.049,0.6225,0.5751,0.043,0.7955,0.4247,0.048,0.6971,0.2643,0.444,0.5386,0.8023,0.449,0.7661,0.6453,0.041,0.6027,0.8531,0.463,-0.0984,0.4493,0.466,-0.0642,0.2244,0.059,-0.0395,0.0697,0.049,0.8702");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mP4_11_e_e"){
        vparameters.push_back("2.8837,1.42393452856,1.61854561848,82.062,0.0387,0.8252,0.5887,0.7184");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_mP10_11_e_e_ef"){
        vparameters.push_back("4.63,1.20259179266,1.52203023758,110.21,0.121,0.1745,0.3531,0.7086,0.4009,0.1165,0.8544,0.5361,0.6943");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP16_11_8e"){
        vparameters.push_back("6.183,0.779880316998,1.77308749798,101.79,0.345,0.162,0.767,0.168,0.128,0.34,0.657,0.457,0.025,0.618,0.473,0.653,0.328,-0.074,0.869,0.894");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mC6_12_a_i"){
        vparameters.push_back("7.189,0.613019891501,0.705105021561,90.04,0.6879,0.2889");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mC34_12_ah3i2j"){
        vparameters.push_back("11.93871,0.876392843113,0.658278825769,129.00411,0.22,0.854,0.241,0.663,0.745,0.566,0.238,0.355,0.232,-0.037,0.333,0.35,0.586");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC16_12_g_ij"){
        vparameters.push_back("5.914,1.73047007102,1.03956712885,108.25,0.1662,0.2147,0.2263,0.2518,0.32131,0.2248");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_mC14_12_a2i_i"){
        vparameters.push_back("9.188,0.430343926861,0.705158902917,97.56,0.14286,0.42857,0.28571,0.85714,0.42857,0.28571");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mC4_12_i"){
        vparameters.push_back("5.403,0.635387747548,0.940033314825,132.32,0.106,0.173");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_mP12_13_e_a_2g"){
        vparameters.push_back("8.95,0.500335195531,1.63360893855,145.35,0.5182,0.2986,0.0278,0.0003,0.2821,0.4045,0.2366");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP84_13_21g"){
        vparameters.push_back("9.21,0.99348534202,2.45385450597,106.1,0.30089,0.20127,0.18147,0.17387,0.03262,0.11695,0.05014,-0.05231,0.18035,-0.07589,0.78099,0.11634,0.79463,0.67872,0.1738,0.68463,0.51532,0.10402,0.56601,0.44932,0.17224,0.42424,0.27741,0.11672,0.0412,0.39067,0.07245,-0.00092,0.15881,0.04497,0.78847,0.13878,0.07346,0.7486,-0.09081,0.04464,0.53574,0.87264,0.06842,0.50833,0.63715,0.03304,0.30515,0.63715,0.06617,0.25041,0.40555,0.0442,0.146,0.38905,0.17219,0.86038,0.10055,0.17357,0.59606,0.82384,0.1694,0.41856,0.64581,0.16732,-0.05418,0.32296,0.2006");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP12_14_2e_e"){
        vparameters.push_back("5.1505,1.01186292593,1.03238520532,99.23,0.07,0.3317,0.3447,0.4496,0.7569,0.4792,0.2754,0.0395,0.2083");
        vparameters.push_back("5.4972,0.880103,1,120.47,0.8881,0.2186,0.7666,0.61,0.7024,0.7014,0.7722,0.9898,0.9889"); // Friedrich-binary-oxide (ICSD #8217) 
        vparameters.push_back("5.7038,0.78679126,0.93574459,122.6,0.106,0.212,0.209,0.401,0.703,0.299,0.239,0.979,0.026"); // Friedrich-binary-oxide (ICSD #647610) 
        vparameters.push_back("5.4128,0.90027343,0.90766332,123.7682,0.192,0.404,0.885,0.321,0.594,0.145,0.235,0,0.01");  // 004, binary metal-nitride prototype (ICSD #240759)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP32_14_8e"){
        vparameters.push_back("9.31,0.866809881847,1.38023630505,93.13333,0.437,0.185,0.084,0.246,0.273,-0.023,0.24,0.102,0.828,0.05,-0.08,0.852,0.157,0.669,-0.09,0.142,0.66,0.09,0.368,0.746,0.16,0.334,0.021,0.21");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP64_14_16e"){
        vparameters.push_back("15.018,0.979691037422,0.585231056066,93.61,0.18313,0.14063,0.03451,0.22856,0.28408,0.12262,0.35548,0.31907,-0.00548,0.47826,0.28776,0.16131,0.52853,0.14438,0.09345,0.47966,0.04033,0.27102,0.36296,-0.02818,0.15123,0.22521,0.04261,0.2343,0.09552,0.48601,0.14213,0.01298,0.58883,0.27815,-0.01931,0.71476,0.12135,0.08347,0.82945,0.18553,0.19177,0.81338,0.00963,0.3102,0.73961,0.14402,0.30834,0.59137,0.04778,0.24353,0.50553,0.23353");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5_mC28_15_f_e2f"){
        vparameters.push_back("12.786,0.387533239481,0.427968090099,97.03333,0.5727,0.106,0.311,0.077,0.0958,0.0952,0.4213,0.7127,0.0726,0.3138");
        vparameters.push_back("13.332,0.378173,0.988741,155.097,0.3306,0.892,0.5091,0.2798,0.8281,0.6701,0.3752,0.4039,0.3353,0.051"); // Friedrich-binary-oxide (ICSD #1422) , BA decoration
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC8_15_c_e"){
        vparameters.push_back("4.6837,0.730747058949,1.0950317057,120.34,0.4184");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC48_15_ae3f_2f"){
        vparameters.push_back("7.1356,1.73344918437,1.00532541062,120.34,0.1163,0.266,0.1234,0.9401,0.3114,0.1038,0.3282,0.0172,0.2117,0.4782,0.14033,0.10833,0.07227,0.50682,0.15799,0.54077");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC6D2_mC40_15_e_e_3f_f"){
        vparameters.push_back("9.79,0.901123595506,0.548518896834,105.81,0.3082,-0.0942,0.3888,0.4123,0.8659,0.1365,0.2411,0.6799,0.1468,0.4802,0.0124,0.2117,0.4057,0.7764");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_oP12_16_ag_cd_2u"){
        vparameters.push_back("5.61,1.01069518717,1.61319073084,0.2,0.26,0.125,0.74,0.8,0.63");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oP16_18_ab_3c"){
        vparameters.push_back("8.32,1.15865384615,0.579326923077,0.0,0.0,0.25,0.25,0.0,0.25,0.5,0.5,0.124,0.309,0.382");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_19_2a_a"){
        vparameters.push_back("7.764,0.909582689335,0.558088614116,0.185,0.07,0.465,0.055,0.765,-0.008,0.884,-0.011,0.391");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC24_20_abc_c"){
        vparameters.push_back("8.74,0.576659038902,0.942791762014,0.3336,0.4403,0.2453,0.1971,0.2713,0.33154,0.03589,0.81143");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP2_25_b_a"){
        vparameters.push_back("2.8102,1.87104120703,1.0769696107,0.0,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP24_28_acd_2c3d"){
        vparameters.push_back("16.54,0.533252720677,0.269649334946,0.0,0.319,0.014,0.018,0.042,0.617,0.042,0.624,0.334,0.5,0.503,0.301,0.042,0.632,0.636,0.5,0.619,0.036,0.5");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oP16_31_a_ab_2ab"){
        vparameters.push_back("7.43,0.869448183042,0.831763122476,0.8268,0.0,0.1514,0.4983,0.8226,0.6454,0.1436,0.1166,0.2466,0.3255,-0.0134,0.2598,0.3364,0.6184");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_33_a_a"){
        vparameters.push_back("5.2857,1.11007056776,0.659950432298,0.1996,0.5867,0.2506,0.002,0.2003,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oP32_33_a_3a_4a"){
        vparameters.push_back("9.11,1.01866081229,1.1613611416,0.2187,0.4807,0.2031,0.4418,0.2052,0.0015,0.4488,0.1967,0.4146,0.1422,0.9176,0.2246,0.191,0.2506,0.2228,0.3424,0.5361,0.0415,0.0069,0.5876,0.2212,0.3355,0.546,0.3761");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC12_36_2a_a"){
        vparameters.push_back("4.624,1.46820934256,2.69139273356,0.333,0.0,0.061,0.134,0.395,0.366");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC8_38_e_a_b"){
        vparameters.push_back("3.875,1.17470967742,1.59019354839,0.0,0.6144,0.155,0.2914");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC12_38_de_ab"){
        vparameters.push_back("4.684,1.81084543126,1.0269000854,0.06,0.5,0.17,0.56,0.17,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_oC20_41_a_2b"){
        vparameters.push_back("6.388,1.00485284909,1.7778647464,0.0,0.673,0.327,0.376,0.827,0.673,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oC24_41_2a_2b"){
        vparameters.push_back("6.478,1.0,1.87635072553,0.01,0.238,0.342,0.158,0.125,0.25,0.25,-0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oF72_43_ab_3b"){
        vparameters.push_back("11.66,1.91595197256,0.58833619211,0.0,0.125,0.13889,0.0,0.02222,0.08056,0.18333,0.15278,-0.01389,-0.18333,0.0625,0.125,0.27778");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oI4_44_a_b"){
        vparameters.push_back("4.92,0.973577235772,0.535569105691,0.0,0.425");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C7D_oP13_47_t_aq_eqrs_h"){
        vparameters.push_back("3.8187,1.01691675177,3.05567339671,0.3554,0.1579,0.3771,0.3788,0.18445");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP4_51_e_f"){
        vparameters.push_back("4.7549,0.661969757513,1.0209678437,0.8125,0.3125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_oP20_56_ce_e"){
        vparameters.push_back("4.911,2.53797597231,1.10201588271,0.029,0.147,0.058,0.861,0.044,0.128,0.179");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD_oP16_57_d_c_d_d"){
        vparameters.push_back("6.707,0.997614432682,1.13627553303,0.208,0.7704,0.2871,0.889,0.4154,0.605,0.1087");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_57_d_d"){
        vparameters.push_back("6.09556,0.900425883758,0.850291031505,0.8593,0.0628,0.255,0.0096");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP6_58_a_g"){
        vparameters.push_back("6.24,1.03044871795,0.673076923077,0.275,0.325");
        vparameters.push_back("4.704,0.917942176871,0.601615646259,0.66667,0.25");
        vparameters.push_back("4.4446,1.22049228277,0.761913333033,0.2004,0.3787");
        vparameters.push_back("4.025,1.1965466,0.67115528,0.8682,0.5945");  // 004, binary metal-nitride prototype (ICSD #160620)
        vparameters.push_back("5.013,1.2288051,0.58567724,0.414,0.068");  // 005, binary metal-nitride prototype (ICSD #166465)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP4_59_a_b"){
        vparameters.push_back("3.15,1.29841269841,2.20634920635,0.051,0.277");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP6_59_a_a_a"){
        vparameters.push_back("5.68,0.700704225352,1.01056338028,0.1499,0.4237,0.6255");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oP8_59_bf_a"){
        vparameters.push_back("5.162,0.842115459124,0.877760557923,0.67125,0.329,0.505,0.174");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP16_61_c_c"){
        vparameters.push_back("6.471,1.27538247566,1.31757070005,0.136,0.072,0.108,0.456,0.119,0.872");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP24_61_2c_c"){
        vparameters.push_back("9.174,0.375953782429,0.560061042075,0.0095,0.1491,0.1835,0.2314,0.111,0.5366,0.1289,0.0972,0.8628");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_oP20_62_3c_2c"){
        vparameters.push_back("11.282,0.339443361106,0.994947704308,0.2922,0.19181,0.4504,0.877,0.6246,0.5611,-0.02937,0.17398,0.64939,-0.03603");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_oP20_62_c_cd_a"){
        vparameters.push_back("5.4224,1.41099881971,0.996661994689,0.4877,-0.0084,0.0313,0.0586,0.288,0.537,0.213");
        vparameters.push_back("6.2152,1.4341775,1.0315356,0.9878,0.4505,0.3956,0.5401,0.6976,0.0554,0.2"); //Friedrich-ternary-oxide (ICSD #97463)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_oP20_62_2cd_c"){
        vparameters.push_back("5.464,0.810395314788,1.36749633968,0.22451,0.65626,0.55801,0.6466,0.05131,0.36362,0.13079,0.0579,0.06543");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP16_62_c_2c_c"){
        vparameters.push_back("6.018,0.630741110003,2.4086075108,0.2522,0.8276,0.6221,0.095,0.8706,0.8244,0.226,0.06333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_62_2c_c"){
        vparameters.push_back("4.918,0.7600650671,1.44550630338,0.038,0.282,0.674,0.562,0.202,0.611");
        vparameters.push_back("12.735,0.468237141735,0.339615233608,0.733,0.125,0.508,0.722,0.874,0.447");
        vparameters.push_back("7.6204,0.595008136056,1.1869718125,0.125,0.4217,0.0202,0.837,0.2377,0.0959");
        vparameters.push_back("3.875,1.64232258065,1.89496774194,0.004,0.758,0.24,0.07,0.24,0.39"); //DX20201019 - moved from Part 2 to here
        vparameters.push_back("9.579,0.29710826,0.43334377,0.174,0.545,0.094,0.259,0.891,0.1");  // 005, binary metal-nitride prototype (ICSD #187446)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_62_c_c"){
        vparameters.push_back("10.42,0.349328214971,0.411708253359,0.375,0.333,0.139,0.389");
        vparameters.push_back("5.24160,0.606723137973,1.12622100122,0.0056,0.1952,0.1879,0.5696");
        vparameters.push_back("5.495,0.536123748863,0.737579617834,0.125,0.69,-0.18,0.125");
        vparameters.push_back("11.18,0.356171735242,0.387209302326,0.3507,0.0201,0.61937,0.3806");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oP16_62_c_cd"){
        vparameters.push_back("5.09,1.3257367387,0.888605108055,0.39,0.05,0.036,0.852,0.186,0.063,0.328");
        vparameters.push_back("7.09,1.0603667,0.67616361,0.376,0.943,0.026,0.085,0.838,0.928,0.384");  // 002, binary metal-nitride prototype (ICSD #165990)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B7_oP40_62_cd_3c2d"){
        vparameters.push_back("4.526,1.54882898807,2.68272205038,0.4594,0.5629,0.0579,0.6261,0.2501,0.2063,0.2619,0.4165,0.0288,0.0291,0.3428,0.0565,0.0642,0.8119,0.2509,0.0657,0.0218");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oP8_62_2c"){
        vparameters.push_back("6.663,0.708839861924,0.73345339937,0.464,0.292,0.181,0.658");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oC16_63_c_2c_c"){
        vparameters.push_back("3.577,4.56863293263,1.09538719597,0.06109,-0.0558,0.1792,0.33096");
        vparameters.push_back("3.5814,5.0438934,1.3735132,0.95974,0.60233,0.78501,0.19337");  // 002, ternary metal-nitride prototype (ICSD #85528)
        vparameters.push_back("4.57,4.2778993,0.98905908,0.826,0.249,0.958,0.609");  // 003, metallic ternary prototype (ICSD #103110)
        vparameters.push_back("4.4247,3.7013583,0.98573915,0.697,0.5638,0.2498,0.9044");  // 004, metallic ternary prototype (ICSD #240096)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC12_63_2c_c"){
        vparameters.push_back("3.73,3.94638069705,0.983914209115,0.061,0.75,0.396");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC8_63_c_c"){
        vparameters.push_back("2.9782,2.64253575985,0.985360284736,0.436,0.14525");
        vparameters.push_back("1.0.0185797325,1.41421356236,0.707106781182,0.625,0.125"); // Lederer-13
        vparameters.push_back("3.31,1.4441088,1.1389728,0.3174,0.6888");  // 003, binary metal-nitride prototype (ICSD #163951)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oC4_63_c"){
        vparameters.push_back("2.8444,2.06331739558,1.73379271551,0.10228");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oC8_64_f"){
        vparameters.push_back("4.523,1.69378730931,1.0002210922,0.1549,0.081");
        vparameters.push_back("3.3136,3.16211974891,1.32070859488,0.10168,0.08056");
        vparameters.push_back("7.11906,0.654575182679,1.37596817557,0.15485,0.1175");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oC80_64_efg_efg_df"){
        vparameters.push_back("10.922,0.866233290606,0.682933528658,0.84657,0.0946,0.9271,0.5886,0.276,-0.0792,0.2314,0.27981,-0.0113,0.1278,0.3415,0.2438,0.1245,0.175,0.2231");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC8_65_j_g"){
        vparameters.push_back("5.971,1.1314687657,0.468263272484,0.28,0.22");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_oC16_65_ah_bej"){
        vparameters.push_back("8.031,0.926410160628,0.491595069107,0.25,0.225");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oC8_65_a_bf"){
        vparameters.push_back("5.82068,1.35259626023,0.493507631411");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oF8_69_a_b"){
        vparameters.push_back("6.08,0.903782894737,0.851973684211");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oF8_70_a"){
        vparameters.push_back("3.1587,1.82613100326,3.21714629436");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oF24_70_e_a"){
        vparameters.push_back("8.2671,0.580614725841,1.0342804611,0.4615");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oF128_70_4h"){
        vparameters.push_back("10.4646,1.22947843205,2.33988876785,0.14415,0.04732,0.0486,0.29277,0.2269,0.25406,0.21598,0.28022,0.32618,0.21405,0.15761,0.37947");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oI6_71_a_i"){
        vparameters.push_back("3.144,0.994910941476,2.44179389313,0.339");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oI6_71_a_g"){
        vparameters.push_back("2.75984,2.9999963766,1.4241115427,0.35333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oI12_72_j_a"){
        vparameters.push_back("9.583,0.585829072316,0.578837524783,0.1182,0.2088");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tI12_82_c_g_a"){
        vparameters.push_back("4.3404,1.53216293429,0.256,0.2566,0.3722");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tI14_82_bc_a_g"){
        vparameters.push_back("5.55,1.85585585586,0.26,0.25,0.13");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP16_84_cej_k"){
        vparameters.push_back("6.429,1.02830922383,0.46779,0.25713,0.19361,0.30754,0.22904");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B5_tI18_87_h_ah"){
        vparameters.push_back("10.164,0.37111373475,0.2797,-0.0589,0.3752,0.6856");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_tI10_87_a_h"){
        vparameters.push_back("5.72,0.623076923077,0.4,0.8");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP12_92_b_a"){
        vparameters.push_back("4.957,1.39001412144,0.3047,0.2381,0.1109,0.1826");
        vparameters.push_back("4.7985,1.584391,0.036776,0.153876,0.270876,0.195876"); // Friedrich-binary-oxide (ICSD #161691)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP36_96_3b_ab"){
        vparameters.push_back("7.464,1.15487674169,0.41,0.445,0.132,0.4,0.117,0.123,0.296,0.344,0.297,0.143,0.326,0.12,0.248");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP12_96_ab"){
        vparameters.push_back("5.51889,1.25999974633,0.0849,0.1752,0.3792,0.2742");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_tP5_99_bc_a_b"){
        vparameters.push_back("4.046,1.02308452793,0.0,0.8973,0.4517,0.3785");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tP8_113_a_ce"){
        vparameters.push_back("6.871,0.606622034638,0.206,0.1797,0.476");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4D_tI16_121_d_a_i_b"){
        vparameters.push_back("5.46,1.96428571429,0.245,0.132");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tI16_122_a_b_d"){
        vparameters.push_back("5.289,1.97069389299,0.2574");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C_tP7_123_b_ci_a"){
        vparameters.push_back("4.207,1.61516520086,0.312");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tP4_123_a_ce"){
        vparameters.push_back("4.158,0.864357864358");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP2_123_a_d"){
        vparameters.push_back("2.8,1.31071428571");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP4_123_d_a_f"){
        vparameters.push_back("3.8611,0.828649866618");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP10_127_g_ah"){
        vparameters.push_back("7.3364,0.530232811733,0.3841,0.1821");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD_tP8_129_c_b_a_c"){
        vparameters.push_back("3.6736,2.60540069686,0.6793,0.2246");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP4_129_ac"){
        vparameters.push_back("4.897,0.69185215438,0.375");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_129_c_a_c"){
        vparameters.push_back("4.11,1.76301703163,0.6497,0.2058");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP6_129_ac_c"){
        vparameters.push_back("4.0006,1.52584612308,0.27,0.7");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP4_129_a_c"){
        vparameters.push_back("3.9645,1.26008323874,0.2368");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP4_129_c_c"){
        vparameters.push_back("3.107,1.90505310589,0.1,0.65");
        vparameters.push_back("1.0,2.82842712472,0.875,0.375"); // Lederer-42
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP4_131_c_e"){
        vparameters.push_back("4.9073,1.24500234345");
        vparameters.push_back("2.979,1.7582746"); //Friedrich-binary-oxide (ICSD #26598)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP50_134_b2m2n"){
        vparameters.push_back("8.74,0.575514874142,0.0048,0.1685,0.1305,0.628,0.1695,0.5228,0.1635,0.0753,0.3383,0.1485");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP30_136_bf2ij"){
        vparameters.push_back("10.59,0.532011331445,0.1033,0.3667,0.0383,0.5608,0.2354,0.3183,0.27");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP8_136_g_f"){
        vparameters.push_back("4.75,0.576842105263,0.31,0.336");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP6_136_f_a"){
        vparameters.push_back("4.5922,0.644005052045,0.30496");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="sigma_tP30_136_bf2ij" || anrl_label == "A_tP30_136_bf2ij"){ //DX20200929 - added "correct" label version
        vparameters.push_back("8.7966,0.518177477662,0.39864,0.13122,0.46349,0.06609,0.73933,0.18267,0.25202");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP4_136_f"){
        vparameters.push_back("3.957,1.29112964367,0.098");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP16_138_j"){
        vparameters.push_back("8.56,0.714953271028,0.375,-0.083,0.857");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tI16_139_cde_e"){
        vparameters.push_back("3.9993,4.3215062636,0.37498,0.11886");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tI4_139_e"){
        vparameters.push_back("3.34916,1.94217355994,0.819");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_tI14_139_a_e_ce"){
        vparameters.push_back("3.7817,3.50337149959,0.36075,0.1824");
        vparameters.push_back("4.6113,3.24412,0.645,0.848"); //Friedrich-ternary-oxide (ICSD #157402)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B_tI26_139_fij_a"){
        vparameters.push_back("8.47,0.584415584416,0.361,0.278");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tI2_139_a"){
        //DX 20191218 [this is the a' parameter (fct) vs the a parameter (bct)] vparameters.push_back("4.6002,1.07523585931");
        vparameters.push_back("3.25283,1.52061313499"); //DX 20191218 [CORRECT PARAMETERS]
        vparameters.push_back("3.932,0.823499491353");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tI8_139_h"){
        vparameters.push_back("4.33184,0.574102459925,0.17916");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tI8_139_bd_a"){
        vparameters.push_back("3.8537,2.22744375535");
        vparameters.push_back("1.0,2.0"); // Lederer-46
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI6_139_a_e"){
        vparameters.push_back("3.2064,2.44754241517,0.3353");
        vparameters.push_back("1.0,4.24264068707,0.3333333333"); // Lederer-44
        vparameters.push_back("4.4371,1.7910121,0.596"); // Friedrich-binary-oxide (ICSD #24248)
        vparameters.push_back("4.7583,1.66103,0.4025"); // Friedrich-binary-oxide (ICSD #38245)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B5_tI18_139_i_ah"){
        vparameters.push_back("8.91,0.361391694725,0.328,0.348");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_tI10_139_de_a"){
        vparameters.push_back("4.53,2.45033112583,0.38");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B_tI18_139_hi_a"){
        vparameters.push_back("8.312,0.468840230991,0.333,0.327");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI6_139_d_a"){
        vparameters.push_back("4.1,1.22682926829");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI12_140_h_a"){
        vparameters.push_back("6.04,0.804635761589,0.158");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI16_140_b_ah"){
        vparameters.push_back("6.017,1.44241316271,0.231");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI16_140_ab_h"){
        vparameters.push_back("8.03,0.87297633873,0.179");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_tI24_141_h_b_a"){
        vparameters.push_back("6.6042,0.905423821205,0.066,0.1951");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tI4_141_a"){
        vparameters.push_back("5.8318,0.545611989437");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_tI28_141_ad_h"){
        vparameters.push_back("5.765,1.63781439722,0.0278,0.2589");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI12_141_e_a"){
        vparameters.push_back("3.785,2.51360634082,0.08306");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI16_141_e_e"){
        vparameters.push_back("3.108,5.45045045045,0.227,0.071");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI24_141_2e_e"){
        vparameters.push_back("4.046,6.28917449333,0.125,0.289,-0.051");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI8_141_a_b"){
        vparameters.push_back("3.325,3.42255639098");
        vparameters.push_back("1.0,1.99999999997"); // Lederer-53
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tI80_141_ceh_3h"){
        vparameters.push_back("7.5937,4.26037373086,0.2044,0.5201,0.3324,0.516,0.2547,0.494,0.0859,0.4667,0.4164");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_tI96_142_e_ab_2g"){
        vparameters.push_back("10.914,1.77396005131,0.0375,0.2482,0.3197,-0.0867,0.0923,0.1117,0.0025");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_147_g_ad"){
        vparameters.push_back("7.636,0.369264012572,0.25,0.33333,0.0,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR16_148_cf_cf"){
        vparameters.push_back("6.29713,1.8633345667,0.11546,0.21,0.10706,0.81289,0.19519,0.1848,0.6754,0.3468");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hR8_148_c_f"){
        vparameters.push_back("7.49626,2.75900649124,0.33333,0.088,0.755,0.421");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR26_148_b2f_a2f"){
        vparameters.push_back("15.659,0.335334312536,0.054,0.346,0.098,0.754,0.15699,0.6,0.555,0.84401,0.599,0.252,0.65501,0.098");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hR10_148_c_f_c"){
        vparameters.push_back("5.0884,2.76815894977,0.35537,0.1464,0.22174,0.56249,0.95095");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_150_ef_bd"){
        vparameters.push_back("5.85,0.589743589744,0.875,0.26,0.6");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP24_151_3c_2a"){
        vparameters.push_back("6.017,2.87518697025,0.8889,0.5556,0.8889,0.1111,0.0731,0.5556,0.4444,0.0731,0.2222,0.77778,0.0731");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_152_c_a"){
        vparameters.push_back("4.914,1.10012210012,0.4699,0.413,0.2668,0.214");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP3_152_a"){
        vparameters.push_back("4.3662,1.13453346159,0.2254");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_154_a_b"){
        vparameters.push_back("4.145,2.29095295537,0.7198,0.4889");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hR8_155_c_de"){
        vparameters.push_back("4.91608,2.53341483458,0.237,0.43,0.07");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hR5_155_e_c"){
        vparameters.push_back("5.73296,1.24097324942,0.2521,0.2449");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR6_160_b_b"){
        vparameters.push_back("9.619,0.327466472606,0.00019,0.26362,0.7288,0.39161");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR6_160_3a_3a"){
        vparameters.push_back("3.01791,7.34847294982,0.0,0.22222,0.77778,0.08333,0.30556,0.86111");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hR10_161_a_a_b"){
        vparameters.push_back("5.2542,2.64091583876,0.2875,0.0128,0.74643,0.14093,0.36263");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP9_162_ad_k"){
        vparameters.push_back("4.917,0.929021761237,0.325,0.272");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2CD2_hP36_163_h_i_bf_i"){
        vparameters.push_back("7.384,2.37716684724,0.01,0.833,0.33333,0.03833,0.141,0.03167,0.365,0.083");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP5_164_ad_d"){
        vparameters.push_back("4.0282,1.21409066084,0.648,0.149");
        vparameters.push_back("4.1971,1.2861976,0.31689,0.76948");  // 002, binary metal-nitride prototype (ICSD #162795)
        vparameters.push_back("3.8572,1.7764181,0.33541,0.76913");  // 003, binary metal-nitride prototype (ICSD #162796)
        vparameters.push_back("3.938,1.2496191,0.606,0.231");  // 004, binary metal-nitride prototype (ICSD #169727)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP3_164_a_d"){
        vparameters.push_back("4.24,1.61320754717,0.252");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP24_165_adg_f"){
        vparameters.push_back("6.308,1.03994927077,0.167,0.666,0.356,0.028,0.096");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR2_166_a_b"){
        vparameters.push_back("3.13,4.78594249201");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hR2_166_c"){
        vparameters.push_back("3.7595,2.7815666977,0.22754");
        vparameters.push_back("2.456,4.08957654723,0.16667");
        vparameters.push_back("3.289,3.42991790818,0.0543");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hR1_166_a"){
        vparameters.push_back("5.07846,0.968139947937");
        vparameters.push_back("3.45741,1.92728082582");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B6_hR13_166_ah_3c"){
        vparameters.push_back("4.757,5.4319949548,0.167,0.346,0.448,0.09,0.59001");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hR3_166_ac"){
        vparameters.push_back("3.62036,7.25049442597,0.22222");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hR5_166_c_ac"){
        vparameters.push_back("4.36914,6.96313919902,0.399,0.208");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_hR7_166_a2c_c"){
        vparameters.push_back("3.011,6.9511790103,0.186,0.33333,0.075");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hR12_166_2h"){
        vparameters.push_back("4.908,2.56022616137,0.0104,0.65729,0.2206,0.6323");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hR4_166_a_b_c"){
        vparameters.push_back("3.5561,5.44557239673,0.2667");
        vparameters.push_back("1.0,4.89897948553,0.25"); //Lederer-59
        vparameters.push_back("3.2903652,5.0765915,0.2375"); //Friedrich-ternary-oxide (ICSD #28288)
        vparameters.push_back("2.8496891,5.6601982,0.899"); //Friedrich-ternary-oxide (ICSD #246912)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hR105_166_bc9h4i"){
        vparameters.push_back("10.96,2.17974452555,0.3848,0.3843,0.21309,0.4895,0.21780,0.3873,0.56899,0.1991,0.50609,0.1983,0.68740,0.1032,0.49209,0.9933,0.66980,0.1008,0.83740,0.0025,0.16801,0.3622,0.58109,0.0976,0.3765,0.68261,0.2024,0.1673,0.55209,0.8921,0.1777,0.3473,0.0033");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B_hR7_166_g_a"){
        vparameters.push_back("4.33304,3.13251204697,0.16667");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hR10_167_a_b_e"){
        vparameters.push_back("5.285,2.62039735099,0.85756666666667");
        vparameters.push_back("4.988,3.42040898156,0.5067");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hR10_167_c_e"){
        vparameters.push_back("4.7607,2.72957758313,0.35216,0.5561");
        vparameters.push_back("4.95878,2.82762,0.654,0.9381"); // Friedrich-binary-oxide (ICSD #94768)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP18_180_fi_bd"){
        vparameters.push_back("5.198,2.54136206233,0.163,0.1141");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP9_180_d_j"){
        vparameters.push_back("4.42758,1.43826876081,0.16559");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_180_j_c"){
        vparameters.push_back("4.9977,1.09252256038,0.2072");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP8_182_c_g"){
        vparameters.push_back("4.8507,0.866967654153,0.3249");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP4_186_ab"){
        vparameters.push_back("2.47,2.75303643725,0.0,0.07143");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP8_186_ab_ab"){
        vparameters.push_back("3.08051,3.27374363336,0.18784,0.0,0.43671,0.24982");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_186_b_b"){
        vparameters.push_back("3.8227,1.63776911607,0.3748,0.0");
        vparameters.push_back("3.143,1.823099,0,0.375");  // 002, binary metal-nitride prototype (ICSD #153888)
        vparameters.push_back("4.084,1.4299706,0.5,0.91");  // 003, binary metal-nitride prototype (ICSD #162195)
        vparameters.push_back("2.965,2.0856661,0,0.3748");  // 004, binary metal-nitride prototype (ICSD #168369)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP12_186_a2b_a2b"){
        vparameters.push_back("3.08129,4.90695780014,0.1254,0.0,0.29215,-0.0415,0.16675,0.8335");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3C_hP18_186_2a3b_2ab_b"){
        vparameters.push_back("3.281,6.57726302956,0.155,0.345,0.0,0.248,0.045,0.261,0.455,0.367,0.137");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_186_b_a"){
        vparameters.push_back("2.51,2.66932270916,0.0,0.05");
        vparameters.push_back("3.172,1.6314628,0,0.38");  // 002, binary metal-nitride prototype (ICSD #159250)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP3_187_a_d_f"){
        vparameters.push_back("4.535,1.0769570011");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP2_187_d_a"){
        vparameters.push_back("2.9065,0.975950455875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_189_fg_bc"){
        vparameters.push_back("5.877,0.584822188191,0.256,0.589");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_hP6_191_a_h_b"){
        vparameters.push_back("3.04436,2.20489035462,0.2413");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5_hP6_191_a_cg"){
        vparameters.push_back("5.405,0.773913043478");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP1_191_a"){
        vparameters.push_back("3.2062,0.931195808122");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP4_191_bc_a"){
        vparameters.push_back("3.6576,1.05902777778");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP3_191_a_d"){
        vparameters.push_back("3.005,1.08276206323");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP6_191_h_e"){
        vparameters.push_back("4.237,1.71040830776,0.306,0.16");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_191_f_ad"){
        vparameters.push_back("5.279,0.806914188293");
        vparameters.push_back("5.1912,0.56069888");  // 002, binary metal-nitride prototype (ICSD #25659)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP8_194_ad_f"){
        vparameters.push_back("3.64,3.37362637363,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP6_194_h"){
        vparameters.push_back("4.40445,0.568892824303,0.44799");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP12_194_af_bf"){
        vparameters.push_back("3.01,4.85382059801,0.166,0.583");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP4_194_ac"){
        vparameters.push_back("3.77,3.2175066313");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP8_194_c_bf"){
        vparameters.push_back("5.088,1.76533018868,-0.083");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP6_194_b_f"){
        vparameters.push_back("4.895,1.58324821246,0.045");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_194_c_d"){
        vparameters.push_back("2.50399,2.66023426611");
        vparameters.push_back("3.03,1.2508251");  // 002, binary metal-nitride prototype (ICSD #163950)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP8_194_d_a_f"){
        vparameters.push_back("2.86,4.48251748252,0.086");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_194_h_c"){
        vparameters.push_back("5.295,0.802077431539,0.8392");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP4_194_bc"){
        vparameters.push_back("2.464,2.72362012987");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP6_194_c_f"){
        vparameters.push_back("3.161,3.8895919013,0.6275");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_hP14_194_abdf_f"){
        vparameters.push_back("2.982,4.651240778,0.528,0.139");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP12_194_f_ah"){
        vparameters.push_back("5.223,1.64005360904,0.06286,0.83048");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_c_d_a"){
        vparameters.push_back("2.752,2.56468023256");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP4_194_f"){
        vparameters.push_back("2.508,1.66786283892,0.05995");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP6_194_c_ad"){
        vparameters.push_back("4.186,1.22527472527");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_hP16_194_c_af_ef"){
        vparameters.push_back("2.988,7.82195448461,0.1543,0.605,0.0539");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP2_194_c"){
        vparameters.push_back("3.2093,1.62359393014");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP24_194_ef_fgh"){
        vparameters.push_back("4.824,3.28067993367,0.04598,0.84417,0.12514,0.16429");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP12_194_df_ce"){
        vparameters.push_back("3.976,4.12022132797,0.0637,0.10724");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_194_c_a"){
        vparameters.push_back("3.619,1.39375518099");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP12_194_cg_f"){
        vparameters.push_back("5.052,1.63697545527,0.062");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_cI40_197_cde_c"){
        vparameters.push_back("8.4295,0.1668,0.3345,0.6476,0.7484");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_cP12_198_a_a_a"){
        vparameters.push_back("5.881,-0.024,0.39,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cP16_198_b_a"){
        vparameters.push_back("5.1305,0.2107,0.3689,0.2671,0.1159");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP8_198_2a"){
        vparameters.push_back("5.65,0.0699,-0.0378");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cP8_198_a_a"){
        vparameters.push_back("5.63,-0.042,0.067");
        vparameters.push_back("4.48688,0.13652,0.8424");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cI16_199_a_a"){
        vparameters.push_back("6.3557,0.294,0.0347");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB32C48_cI162_204_a_2efg_2gh"){
        vparameters.push_back("14.16,0.8203,0.5998,0.1836,0.2942,0.8806,0.0908,0.8499,0.1748,0.6993,0.686,0.0969,0.332");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cI32_204_g_c"){
        vparameters.push_back("7.58,0.3431,0.8497");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B_cI26_204_g_a"){
        vparameters.push_back("7.58,0.184,0.691");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP8_205_c"){
        vparameters.push_back("5.65,0.05569");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cP16_205_c_c"){
        vparameters.push_back("6.4162,0.1527,0.6297");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_cP12_205_a_c"){
        vparameters.push_back("5.417,0.3851");
        vparameters.push_back("5.9431,0.427"); //DX20200703 - binary oxide for R. Friedrich (ICSD #87178)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C6_cI80_206_b_d_e"){
        vparameters.push_back("9.4,-0.0344,0.338,0.1,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cI16_206_c"){
        vparameters.push_back("4.11971,0.1001");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP20_213_cd"){
        vparameters.push_back("6.315,0.06361,0.20224");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C_cP8_215_d_e_a"){
        vparameters.push_back("5.3912,0.2372");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_cP5_215_a_e"){
        vparameters.push_back("3.878,0.265");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_cP8_215_a_c_e"){
        vparameters.push_back("5.28,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5_cF24_216_a_ce"){
        vparameters.push_back("6.1,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_cF12_216_b_c_a"){
        vparameters.push_back("6.24");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cF8_216_c_a"){
        vparameters.push_back("5.4093");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_cI10_217_c_a"){
        vparameters.push_back("5.45858,0.165");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cI58_217_ac2g"){
        vparameters.push_back("8.911,0.31787,-0.08958,0.28194,0.64294,0.03457");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B8_cI52_217_ce_cg"){
        vparameters.push_back("8.8664,0.32774,0.10781,0.64421,0.68844,0.03674");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cI16_220_c"){
        vparameters.push_back("5.2716,0.049");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_cI40_220_d_c"){
        vparameters.push_back("8.135,0.0492,0.2896");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cP2_221_b_a"){
        vparameters.push_back("4.07925");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cP6_221_c_d"){
        vparameters.push_back("4.2101");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_cP5_221_a_c_b"){
        vparameters.push_back("3.795");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB27CD3_cP32_221_a_dij_b_c"){
        vparameters.push_back("7.04,0.245,0.26");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_cP4_221_a_c"){
        vparameters.push_back("3.7402");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP1_221_a"){
        vparameters.push_back("3.34");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB11_cP36_221_c_agij"){
        vparameters.push_back("9.6,0.345,0.225,0.115");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB11CD3_cP16_221_a_dg_b_c"){
        vparameters.push_back("5.74,0.245");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cP4_221_d_a"){
        vparameters.push_back("3.734");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B_cP7_221_f_a"){
        vparameters.push_back("4.145,0.2117");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cP8_223_c_a"){
        vparameters.push_back("4.556");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP46_223_dik"){
        vparameters.push_back("10.355,0.1837,0.1172,0.3077");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cP6_224_b_a"){
        vparameters.push_back("4.267");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B_cF32_225_bd_a"){
        vparameters.push_back("9.45");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_cF16_225_a_bc"){
        vparameters.push_back("5.853");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B16C7_cF128_225_acd_2f_be"){
        vparameters.push_back("11.48,0.25,0.875,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B_cF52_225_i_a"){
        vparameters.push_back("7.468,0.666");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_cF12_225_a_c"){
        vparameters.push_back("5.4631");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B23_cF116_225_e_acfh"){
        vparameters.push_back("10.65,0.2765,0.6191,0.6699");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_cF16_225_a_c_b"){
        vparameters.push_back("5.95");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cF4_225_a"){
        vparameters.push_back("3.61491");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB18C8_cF108_225_a_eh_f"){
        vparameters.push_back("10.56,0.325,0.65833,0.66");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cF8_225_a_b"){
        vparameters.push_back("5.63931");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cF24_227_c_a"){
        vparameters.push_back("7.166");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_cF96_227_e_cf"){
        vparameters.push_back("11.278,0.215,0.44");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cF16_227_a_b"){
        vparameters.push_back("7.483");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cF136_227_aeg"){
        vparameters.push_back("14.864,0.2624,0.1824,0.3701");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cF24_227_d_a"){
        vparameters.push_back("7.02");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cF8_227_a"){
        vparameters.push_back("3.55");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_cF56_227_d_a_e"){
        vparameters.push_back("8.0832,0.7376");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_cF48_227_c_e"){
        vparameters.push_back("8.6,0.245");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C3_cF112_227_c_de_f"){
        vparameters.push_back("11.087,0.7047,0.323");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cI2_229_a"){
        vparameters.push_back("3.155");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cI8_229_b_a"){
        vparameters.push_back("2.984");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_cI14_229_c_b"){
        vparameters.push_back("6.226");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7_cI54_229_e_afh"){
        vparameters.push_back("11.618,0.6862,0.1704,0.6503");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB12C3_cI32_229_a_h_b"){
        vparameters.push_back("7.04,0.7625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C3_cI16_229_a_c_b"){
        vparameters.push_back("5.74");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_cI112_230_af_g"){
        vparameters.push_back("11.411,0.0,0.625");
      }
    }
    // ---------------------------------------------------------------------------
    // Part 2
    // ---------------------------------------------------------------------------
    if(library=="all" || library == "" || library=="part2"){
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_aP6_2_aei_i"){
        vparameters.push_back("2.7804,1.00438785786,1.53100273342,78.28,76.53,70.42,0.259,0.415,0.663,0.192,0.183,0.255");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B5_mP13_6_a7b_3a2b"){
        vparameters.push_back("6.5369054321,0.490874879533,1.43786810261,109.592,0.5119,0.7172,0.4546,0.4285,0.044,0.5911,0.37,0.0053,0.3887,0.2195,0.6277,0.0001,0.1818,0.4663,0.7456,0.5256,0.1826,0.7901,0.7824,0.2724,0.0,0.0,0.8188,0.8026,0.0123,0.2051");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mP4_6_2b_2a"){
        vparameters.push_back("3.580975428,1.00027925162,1.00167550965,90.04,0.0,0.0,0.518,0.507,0.026,0.501,0.529,0.027");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP12_7_4a_2a"){
        vparameters.push_back("5.0942,0.627360527659,1.04603274312,90.38,0.498,-0.062,0.472,-0.023,0.574,0.143,0.777,-0.052,0.799,0.271,0.151,0.261,-0.001,0.808,0.36,0.494,0.649,0.658");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP18_7_6a_3a"){
        vparameters.push_back("6.5499839723,1.91328244277,1.22717557252,127.75,0.6753,0.4828,0.3393,0.4711,0.6498,0.3067,0.4174,0.3433,0.0893,0.2883,0.1551,0.3432,0.1472,0.1112,0.5582,0.0,0.0746,0.0,0.5693,0.07852,0.0933,0.1184,0.41567,0.2938,0.8473,0.27155,0.6656");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_mP16_7_6a_2a"){
        vparameters.push_back("5.2778002048,0.97690325515,1.45210125431,91.762,0.5044,0.292,0.01,0.5764,0.215,0.586,0.0,0.209,0.0,0.0864,0.29,0.58,0.2874,0.0717,0.287,0.7924,0.4201,0.301,0.2874,0.014,0.0012,0.7994,0.528,0.078");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B2_mP22_7_9a_2a"){
        vparameters.push_back("6.4165085102,0.999298672153,1.36910105356,93.39,0.4944,0.2369,0.0135,0.2965,0.2623,0.5568,0.4955,0.4453,0.2912,0.2817,0.0557,0.2805,0.8067,0.5491,0.0506,0.0,0.0309,0.0,0.6774,0.1539,0.7442,0.096,0.6343,0.3258,0.8638,0.2389,0.2526,0.6359,0.1217,0.4513,0.156,0.3761,0.11377");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_mC32_9_5a_3a"){
        vparameters.push_back("8.12077,0.718445418353,1.12797801194,115.809,0.009,-0.003,0.269,0.129,0.341,0.45,0.37,0.119,0.066,0.142,0.351,0.147,0.356,0.135,0.348,0.0,0.5182,0.0,0.136,0.2,0.309,0.365,0.2924,0.196");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC16_9_a_3a"){
        vparameters.push_back("3.5367,2.67777872028,0.986823875364,93.018,0.50691,0.35551,0.49997,0.28659,0.26502,0.27913,0.57076,0.06256,0.41996,0.42173,0.07037,0.56036");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP6_10_mn_bg"){
        vparameters.push_back("4.012,0.819541375872,2.93668993021,97.03,0.843,0.126,0.558,0.644");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mP16_10_mn_3m3n"){
        vparameters.push_back("3.678,0.71098423056,1.23817292007,91.0,0.263,0.339,0.08,0.059,0.795,0.341,0.438,-0.069,0.763,0.161,0.705,0.841,0.58,0.441,0.062,0.431");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mP8_10_ac_eh_mn"){
        vparameters.push_back("5.1177434753,0.86107854631,1.45123094959,90.021,0.6089,0.24179,0.1277,0.24913");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mP6_10_en_am"){
        vparameters.push_back("5.1700416367,0.61508704062,1.49709864605,104.5,0.234,0.66,0.263,0.336");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP8_10_2m2n"){
        vparameters.push_back("4.7302,0.527461840937,0.86332501797,106.1,0.1175,0.6746,0.5344,0.3333,0.1131,0.8977,0.4209,0.1319");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C2_mC22_12_aij_h_i"){
        vparameters.push_back("6.65,1.29563909774,0.704661654135,102.2,0.30503,0.38654,0.22171,0.22108,-0.08762,0.23655,0.15499,0.71826");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mC16_12_4i"){
        vparameters.push_back("9.089,0.274617669711,0.451534822313,96.96,-0.0572,0.1206,0.4419,0.3467,0.7858,-0.0594,0.2715,0.4149");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP12_13_2g_ef"){
        vparameters.push_back("5.6255,0.61198115723,1.23780997245,127.44,0.1808,-0.004,0.155,0.346,0.225,0.345,0.273,0.573");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP6_14_e_a"){
        vparameters.push_back("5.5496,0.695689779444,1.15512829753,107.151,0.255,0.2573,0.3141");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B8_mP120_14_14e_16e"){
        vparameters.push_back("7.5889,0.766725085322,3.55545599494,106.136,0.25252,0.5751,0.3402,0.76924,0.888,0.47034,0.13489,0.3866,0.33222,0.8218,0.8146,0.42747,0.12506,0.2343,0.29191,0.78423,0.9464,0.38285,0.23276,0.2673,0.25885,0.69189,0.1533,0.38027,0.34844,0.4547,0.26591,0.63874,0.2281,0.42248,0.35819,0.6069,0.3061,0.67752,0.0969,0.46714,0.26777,0.7386,0.3843,0.8119,0.7461,0.51889,0.0602,0.3617,0.3547,0.8843,0.6723,0.4287,0.0438,0.1068,0.2871,0.822,0.8945,0.354,0.2273,0.1621,0.2315,0.6653,0.243,0.3497,0.4219,0.4796,0.2431,0.5754,0.3699,0.4209,0.4385,0.7352,0.3104,0.6409,0.1506,0.496,0.9409,0.7672,0.5381,0.7891,0.5835,0.5099,0.7334,0.7952,0.5403,0.2081,0.8842,0.3711,0.3975,0.7664,0.4019,0.2077,0.6717,0.4087");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC16_15_e_cf"){
        vparameters.push_back("3.282,2.64137720902,0.970749542962,91.9,0.86,0.577,0.068,0.169");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mC24_15_2e2f"){
        vparameters.push_back("4.939,0.569143551326,0.838023891476,142.47,0.1012,0.3684,0.226,0.0672,0.2464,0.3443,0.1958,0.2227");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_17_abe_e"){
        vparameters.push_back("7.0499703347,1.11347517732,0.614184397158,0.893,0.878,0.379,0.225,0.522,0.202,0.275,0.022");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oP16_19_a_3a"){
        vparameters.push_back("5.668,0.758997882851,0.528757939308,0.584,0.123,0.027,0.31,0.159,0.417,0.257,0.073,0.603,0.983,0.124,0.227");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oC6_21_a_k"){
        vparameters.push_back("3.3982513,1.3943496174,1.40170688639,0.268");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_oF40_22_fi_ad_gh"){
        vparameters.push_back("6.4793068924,1.3977127159,1.54737394473,0.3134,0.3627,0.1144,0.0719");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oF8_22_a_c"){
        vparameters.push_back("5.5400291632,0.990433212996,0.937725631773");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oI32_23_ij2k_k"){
        vparameters.push_back("5.82463,1.24369101557,1.32254065924,0.04851,0.45153,0.75475,0.49255,0.20405,0.4421,0.23223,0.28616,0.76005,0.8216,0.36488");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B2C12D2E_oI50_23_bcfk_i_3k_j_a"){
        vparameters.push_back("10.767,0.502554100492,1.4969815176,0.2511,0.3298,0.1693,0.2465,0.0107,0.1695,0.1308,0.2443,0.0826,0.3792,0.7558,0.0801,0.1294,0.7488,0.2546");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oI16_23_ab_i_k"){
        vparameters.push_back("5.3999384419,1.15740740741,2.00555555555,0.28,0.25,0.2,0.115");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_oI12_23_a_b_k"){
        vparameters.push_back("5.2501580231,1.06666666665,1.7219047619,0.21,0.2,0.115");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB7CD2_oI44_24_a_b3d_c_ac"){
        vparameters.push_back("7.0501914381,1.035035461,1.41546099291,0.7511,0.2496,0.1361,-0.0003,0.5002,0.7501,0.2213,0.3356,0.5662,0.0681,0.1362,-0.0648,0.0709,0.1385");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_26_abc_ab"){
        vparameters.push_back("4.6806,0.627034995513,1.05710806307,0.455,0.858,0.179,0.623,0.048,0.545,0.375,0.355,0.751,0.119,0.213");
        vparameters.push_back("5.0725824349,0.881353258932,1.48474034935,0.12,0.0,0.2484,0.461,0.246,0.289,0.3781,0.0862,0.253,0.652,0.12");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B_oP24_26_3a3b2c_ab"){
        vparameters.push_back("6.465013742,1.07113689096,1.87440061873,0.1628,0.2907,0.112,0.0,0.1452,0.8256,0.4859,0.1227,0.1146,0.0096,0.3286,0.1358,0.1565,0.2905,0.7214,0.2818,0.2489,0.2854,0.3947,0.249,0.0963,0.5436");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B4C16D_oP108_27_abcd4e_4e_16e_e"){
        vparameters.push_back("13.0409854622,0.999877309088,0.703110981603,0.0,0.005,0.03,0.05,0.373,0.379,0.296,0.122,0.125,0.249,0.369,0.129,0.28,0.122,0.376,0.26,0.058,0.751,0.449,0.237,0.574,0.086,0.257,0.012,0.554,0.483,0.247,0.028,0.729,0.183,0.412,0.162,0.684,0.26,0.23,0.164,0.679,0.651,0.299,0.6217,0.4,0.254,0.238,0.433,0.091,0.44,0.445,0.395,0.455,0.245,0.09,0.304,0.404,0.057,0.126,0.093,0.054,0.102,0.243,0.402,0.336,0.099,0.466,0.124,0.392,0.467,0.154,0.103,0.253,0.191,0.047,0.385,0.425,0.052,0.111,0.41,0.74634,0.2513,0.28218");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_29_2a_a"){
        vparameters.push_back("5.2594682584,0.963498098863,0.965209125484,0.639,0.068,0.0,0.771,0.537,0.106,0.53,0.267,0.356");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP12_29_a_2a"){
        vparameters.push_back("5.4179557747,1.0,1.0,0.5049,0.2419,0.0,0.615,0.135,0.3834,0.615,0.635,0.1134");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP12_29_a_a_a"){
        vparameters.push_back("5.2594682584,0.963498098863,0.965209125484,0.61885,0.63065,0.11668,0.50496,0.24091,0.0,0.61734,0.13129,0.37996");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3C15_oP46_30_a2c_bc_a7c"){
        vparameters.push_back("5.4630061801,1.00183049606,3.84605528098,0.5,0.107,0.523,0.2442,0.513,0.019,0.3675,0.026,0.016,0.1074,0.003,0.02,0.1949,0.05,0.074,0.097,0.687,0.22,0.0844,0.254,0.226,0.413,0.52,0.105,0.505,0.718,0.31,0.305,0.763,0.271,0.695,0.752,0.289");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oP20_30_2a_c_3c"){
        vparameters.push_back("4.480982758,1.71211783083,3.19772372237,0.8543,0.5,0.3011,0.721,0.4096,0.3607,0.7305,0.1733,0.5847,0.8629,0.0496,0.6302,0.8523,0.2991");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A13B2C2_oP34_32_a6c_c_c"){
        vparameters.push_back("8.4261784988,1.1450866366,0.701103726564,0.0,0.808,0.55,0.8946,0.614,0.7595,-0.0255,0.8197,0.8044,0.642,0.5367,0.6481,0.3894,0.7689,-0.0703,0.2917,0.8392,0.688,0.2793,0.67666,0.60773,0.0844,0.8643,0.8235,0.4147");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_oP40_33_4a_6a"){
        vparameters.push_back("4.8437,1.71975968784,1.84873134174,0.6787,0.8416,0.0,0.1846,0.3432,0.7868,0.8115,0.6489,0.6972,0.6677,0.4696,0.9993,0.329,0.8313,0.8927,0.0248,0.4908,0.6292,0.4717,0.6647,0.6381,0.5145,0.6728,0.1212,0.8608,0.3301,0.8662,0.336,0.4992,0.9");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8C_oP22_34_c_4c_a"){
        vparameters.push_back("6.2740008991,2.04032515143,1.38284985656,0.5,0.605,0.8135,0.499,0.7349,0.6796,0.009,0.743,-0.0925,0.3064,0.2388,0.5935,0.2196,0.7131,0.6459,0.513");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP6_34_a_c"){
        vparameters.push_back("5.8327827022,1.12083390481,0.548158688784,0.5,0.6881,0.8565,0.0097");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C2_oC22_35_a_ab3e_e"){
        vparameters.push_back("4.534758023,1.1408839779,2.72580110498,0.0,0.5838,-0.042,0.189,0.5461,0.0961,0.122,0.2982,0.1399,0.1866,0.0038");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC8_36_a_a"){
        vparameters.push_back("5.825,0.945115879828,0.922403433476,0.25,0.0,0.081,0.83");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C2_oC36_37_d_c2d_d"){
        vparameters.push_back("5.8071602545,2.51110728429,0.821939039086,0.0,0.654,0.0584,0.0469,0.6705,0.0718,0.471,0.0932,0.1377,0.4004,0.1552,0.14836,0.0571");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_oC40_39_2d_2c2d"){
        vparameters.push_back("7.47831,2.30292940517,0.749515599113,0.3449,0.5,0.0058,0.7654,0.173,0.5398,0.3889,0.5843,0.6154,0.3917,0.28104,0.16619,0.0206,0.10455,0.6073,0.0144");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9BC_oC44_39_3c3d_a_c"){
        vparameters.push_back("6.2123847492,1.90215711527,2.62509658728,0.5,0.2673,0.7493,0.6176,0.8879,0.6246,0.6059,0.6191,0.7476,0.2266,0.0857,0.2459,0.1789,0.5932,0.0649,0.18,0.5948,0.4281");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oC16_40_a_2b_b"){
        vparameters.push_back("6.4580033647,1.33199132859,1.69634561784,0.0,0.3467,0.1773,0.8661,0.3301,0.7281,0.0093");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oC16_40_b_3b"){
        vparameters.push_back("4.3850022361,5.92335058951,0.997331752147,0.83109,0.0,0.70417,0.0021,0.56971,0.4966,-0.07002,0.4978");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B3_oF52_42_2abce_ab"){
        vparameters.push_back("7.4494846573,1.70036689767,1.04688136975,0.76,0.27,0.0,0.31,0.06,0.79,0.6,0.17,0.11,0.07");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oF8_42_a_a"){
        vparameters.push_back("2.5000573158,1.33999999997,1.73599999999,0.0,0.333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_oI20_45_c_b_c"){
        vparameters.push_back("5.9678340183,1.97721179624,0.981568364609,0.25,0.2729,0.629,0.493,0.2296,0.8629,0.474");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oI36_46_ac_bc_3b"){
        vparameters.push_back("6.9969353511,1.54780620265,0.898527940538,0.25,0.5253,0.0054,0.7207,0.7706,-0.0021,-0.0823,0.2996,0.7963,0.5295,0.6236,0.1199,0.006,0.3325,0.4952");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8CD_oP24_48_k_2m_d_b"){
        vparameters.push_back("6.3302356554,1.0,1.50710900473,0.4978,0.56,0.629,0.114,0.642,0.558,0.601");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_oP14_49_dehq_ab"){
        vparameters.push_back("3.6705354001,1.69539132804,2.12544314151,0.002,0.681");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C8D_oP24_49_g_q_2qr_e"){
        vparameters.push_back("6.3302356554,1.0,1.50710900473,0.5184,0.7996,0.258,-0.0697,0.3652,0.6307,0.2617,0.1943,0.8286");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_oP28_50_ij_ac_ijm"){
        vparameters.push_back("5.5348069961,2.26684733515,0.987895212291,0.8726,0.445,0.3867,-0.076,0.5,0.743,0.226");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC2_oP48_50_3m_m_2m"){
        vparameters.push_back("11.0940438033,1.50045069408,0.472480620154,0.0515,0.8152,0.238,0.602,0.5567,0.315,0.604,0.8447,0.116,0.61286,0.66266,0.2344,0.11964,0.66839,0.2479,0.63183,0.00412,0.2439");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP24_52_2e_cd"){
        vparameters.push_back("7.2235008877,1.34576036547,1.32098013428,0.3175,0.6759,0.8271,0.1762,0.5576,0.5093,0.0419,0.8142");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_oP20_52_de_cd"){
        vparameters.push_back("15.5832022616,0.434585119311,0.426069989123,0.443,0.4294,0.0001,0.6539,0.064,0.0788");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP16_53_h_e_gh"){
        vparameters.push_back("11.6904618594,0.282132455361,0.78890717994,0.79514,0.3193,0.1198,0.3549,0.224,0.7493");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oP20_53_e_g_hi"){
        vparameters.push_back("14.3630679002,0.312469539787,0.535821207264,0.6826,0.2856,0.3458,0.2708,0.6247,0.6057,0.3575");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oP20_54_e_d_cf"){
        vparameters.push_back("5.3467489374,0.956627452436,1.81710213776,0.8667,0.6417,0.8902,0.4055,0.2686,0.5503");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP24_55_2g2h_gh"){
        vparameters.push_back("10.1600191136,1.45275590551,0.366929133855,0.0628,0.4022,0.1014,0.1118,0.2024,0.2667,0.226,0.0384,0.3532,0.2953,0.4192,0.1378");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_oP16_55_ch_agh"){
        vparameters.push_back("5.4199981729,1.90405904061,0.730627306278,0.348,0.22,0.112,0.152,0.17,0.393");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oP16_55_2g2h"){
        vparameters.push_back("7.7886,0.613101199189,0.320442698303,0.6731,-0.037,0.8435,0.8087,-0.0454,0.8613,0.5704,0.8926");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP6_58_g_a"){
        vparameters.push_back("3.7572,2.89952624295,0.89063664431,0.3326,0.63309");
        vparameters.push_back("2.7,1.8259259,1.5296296,0.87249,0.40386");  // 002, binary metal-nitride prototype (ICSD #157283)
        vparameters.push_back("3.778,1.2916887,0.84912652,0.852,0.584");  // 003, binary metal-nitride prototype (ICSD #166463)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP6_59_a_b_a"){
        vparameters.push_back("3.301,1.14298697364,2.39612238716,0.32961,-0.04795,0.89243");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_oP20_60_d_cd"){
        vparameters.push_back("8.4598035458,0.706855791961,0.724586288418,0.547,0.394,0.75,0.033,0.348,0.611,0.396");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oP32_60_3d_d"){
        vparameters.push_back("7.3397836195,1.05524748967,1.03197678378,0.5016,0.7205,0.0322,0.2167,0.7591,0.2582,0.2197,0.5016,0.013,0.248,0.783,0.0291");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B8_oP120_60_7d_8d"){
        vparameters.push_back("13.85,0.824548736462,0.53357400722,0.382,0.3932,0.5783,0.2764,0.38,0.5411,0.2264,0.4631,0.4411,0.1289,0.4507,0.4069,0.0794,0.3554,0.4717,0.1276,0.2723,0.571,0.2251,0.2844,0.6054,0.2616,0.5329,0.393,0.094,0.5116,0.3344,0.0088,0.3466,0.4468,0.0917,0.2026,0.6185,0.2593,0.2232,0.6779,0.3972,0.4784,0.595,0.4192,0.3624,0.4724,0.4002,0.3492,0.6895");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP48_61_3c_3c"){
        vparameters.push_back("6.914,1.08128435059,1.38313566676,0.1297,0.5762,0.40803,0.1235,0.6328,0.54518,0.0057,0.4432,0.36289,0.2172,0.6275,0.346,0.2068,0.7225,0.5756,0.0095,0.4051,0.2704");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_oP20_62_2c_3c"){
        vparameters.push_back("5.5329,0.511305102207,2.07339731425,0.1008,0.2055,0.2432,-0.0464,0.0157,0.4015,0.1808,0.7737,0.8691,-0.0688");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_oP28_62_ac_2cd_c"){
        vparameters.push_back("10.193,0.586382811734,0.466202295693,0.2774,-0.0085,0.0913,0.7657,0.4474,0.2215,0.094,0.4262,0.1628,0.0331,0.2777");
        vparameters.push_back("9.3538,0.698497,0.550536,0.134,0.4967,0.1506,0.0677,0.3502,0.7211,0.3161,0.0161,0.6155,0.5647,0.8583"); // Friedrich-ternary-oxide (ICSD #16298)
      }
      // ---------------------------------------------------------------------------
      //DX20201019 [OBSOLETE - moved up ] if(anrl_label=="A2B_oP12_62_2c_c"){
      //DX20201019 [OBSOLETE - moved up ]   vparameters.push_back("3.875,1.64232258065,1.89496774194,0.004,0.758,0.24,0.07,0.24,0.39");
      //DX20201019 [OBSOLETE - moved up ] }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oP16_62_cd_c"){
        vparameters.push_back("6.5982,1.11416750023,0.727789397108,0.011,0.415,0.369,0.555,0.174,0.053,0.856");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3_oP24_62_c_d_cd"){
        vparameters.push_back("6.231,1.78414379714,1.03787514043,0.123,0.0823,0.2579,0.4127,0.1366,0.087,0.5853,0.267,0.0846,-0.088");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oP16_62_c_3c"){
        vparameters.push_back("13.855,0.266791771923,0.28601948755,0.398,0.425,0.074,-0.026,0.414,-0.066,0.276,0.49");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oP24_62_c_2cd_c"){
        vparameters.push_back("8.884,0.614362899595,0.805155335434,0.8154,0.3419,0.5878,0.6062,0.3192,0.5515,0.437,0.6914,0.4186,0.4702,0.819");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_62_c_c"){
        vparameters.push_back("5.454,0.609644297763,1.10542720939,0.2005,0.5741,0.0058,0.1993");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_oC24_63_e_c_cg"){
        vparameters.push_back("9.049,1.21770361366,0.600176815118,0.6699,0.1191,0.8502,0.2174,0.3859");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h"){
        vparameters.push_back("10.11895,1.73974572461,4.16742843872,0.38998,0.0432,0.5987,0.1275,0.13992,0.17412,0.7278,0.20648,0.19723,0.62632,0.1388,0.02884,0.39065,0.0988,0.06788,0.55569,0.59876,0.1328,0.28079,0.54761,0.0678,0.6912,0.07797,0.09224,0.45337,0.1675,0.43146,0.0091,0.18607,0.20422,0.3338,0.3774,0.18828,0.32898,0.174,0.19344,0.2014,0.10123,0.30731,0.03499,0.20655,0.30495,0.49596,0.12681,0.19578,0.33018,0.026,0.31697,0.03483,0.04522,0.2812,0.17213,0.16828,0.2807,0.35936,0.09067");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B_oC28_63_efg_c"){
        vparameters.push_back("7.5551,0.860266574896,1.17435904224,0.45686,0.32602,0.13917,0.10039,0.31768,0.28622");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_oC20_63_a_cf_c"){
        vparameters.push_back("2.456,3.27442996743,2.48086319218,0.077,0.747,0.631,-0.064");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oC24_63_a_fg_c"){
        vparameters.push_back("5.182,1.52315708221,1.25549980702,0.37,0.25,0.06,0.25,0.47");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oC24_63_c_fg_c"){
        vparameters.push_back("6.995,0.892780557541,0.999714081487,0.6524,0.15556,0.7025,-0.0819,0.1699,0.0162");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC24_64_2f_f"){
        vparameters.push_back("2.9986,1.44314013206,2.534682852,0.372,0.27,0.6,0.092,0.371,0.615");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_oC28_66_l_kl_a"){
        vparameters.push_back("6.2700590867,1.72567783094,1.73046251993,0.833,0.005,0.268,0.737,0.42");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oC64_66_gi2lm_2l"){
        vparameters.push_back("8.157,1.00294225818,0.59212945936,0.54552,0.3287,0.39296,0.16145,0.33837,0.89039,0.24136,0.07877,0.42337,0.74171,0.33285,0.66798,0.24948");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oC64_66_kl2m_bdl"){
        vparameters.push_back("8.735,2.32364052662,1.67842014883,0.1826,-0.0318,0.1994,0.327,0.1716,0.2894,0.451,0.1302,0.1133,0.3773,0.3708");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_67_ag_b_g"){
        vparameters.push_back("5.0644067238,1.60320657112,1.02379259962,0.3305,0.8198");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_67_b_g_ag"){
        vparameters.push_back("5.2729860526,1.00606865161,1.82912952778,0.765,0.3403");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC8_67_a_g"){
        vparameters.push_back("5.32495,0.997010300566,1.02896740814,0.26686");
        vparameters.push_back("5.6124272786,0.999376380873,0.8895303257,0.7642");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_oC20_68_a_i"){
        vparameters.push_back("6.44222,1.77662048176,0.991771470083,0.3313,0.1233,0.0801");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oF48_70_f_fg"){
        vparameters.push_back("4.2082,3.45504015969,1.7326647973,0.2495,0.54337,0.29445");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_oI14_71_gh_cg"){
        vparameters.push_back("3.29,4.25531914894,0.951367781155,0.375,0.18,0.444");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oI12_71_h_j_g"){
        vparameters.push_back("3.438,3.4554973822,1.37434554974,0.212,0.1232,0.235");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD3_oI48_73_d_e_e_ef"){
        vparameters.push_back("5.7750411261,1.02926406926,3.53056277057,0.63427,0.3734,0.18032,0.311,0.1349,0.6124,0.0967");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oI12_74_h_e"){
        vparameters.push_back("5.158858099,1.56976744188,1.70096899227,-0.047,0.56,0.663");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_oI20_74_beh_e"){
        vparameters.push_back("4.39,1.42369020501,3.12528473804,-0.111,0.111,-0.033,0.314");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C12D4_tP76_75_2a2b_2d_12d_4d"){
        vparameters.push_back("9.8880614494,0.922431229769,0.3333,0.0,0.5003,0.1666,0.167,0.348,0.1666,0.333,0.152,0.0,0.208,0.152,0.8333,0.458,0.152,0.8333,0.292,0.348,0.6666,0.042,0.348,0.6666,0.208,0.152,0.5,0.458,0.152,0.5,0.292,0.348,0.3333,0.042,0.348,0.3333,0.208,0.152,0.1666,0.458,0.152,0.1666,0.292,0.348,0.0,0.042,0.348,0.0,0.167,0.348,0.8333,0.333,0.152,0.6666,0.167,0.348,0.5,0.333,0.152,0.3333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tP16_76_2a_a_a"){
        vparameters.push_back("3.9810644174,3.85606631503,0.143,0.173,0.1347,0.353,0.141,0.2027,0.3461,0.3475,0.5937,0.1519,0.1578,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B7_tP40_76_3a_7a"){
        vparameters.push_back("9.046,1.84766747734,0.7435,0.3852,0.0,0.4169,0.733,0.8359,0.026,0.8404,-0.0086,0.79,0.6,0.8105,0.443,0.095,-0.053,0.106,0.473,0.8928,0.357,0.024,0.061,0.629,0.794,0.017,-0.002,0.341,0.7049,0.011,0.29,0.84");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B6CD7_tP64_77_2d_6d_d_ab6d"){
        vparameters.push_back("7.6159522234,1.07480314961,0.0,0.035,0.1167,0.1203,0.418,0.392,0.3817,0.132,0.362,0.099,0.309,0.136,0.393,0.161,0.192,0.355,0.442,0.312,0.133,0.03,0.028,0.14,0.14,0.47,0.35,0.32,0.2491,0.7483,0.273,0.2602,0.0188,0.013,0.2316,0.4834,0.184,0.1644,0.2577,0.529,0.3406,0.2365,0.012,0.0182,0.2119,0.275,0.4808,0.3008,0.268");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP48_77_8d_4d"){
        vparameters.push_back("13.47,0.304380103935,0.193,0.14,0.163,0.091,0.059,0.837,0.693,0.14,0.163,0.591,0.059,0.837,0.193,0.64,0.163,0.193,0.64,0.837,0.693,0.64,0.837,0.591,0.559,0.163,0.11,0.14,0.0,0.61,0.14,0.0,0.11,0.64,0.0,0.61,0.64,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C2_tP88_78_4a_14a_4a"){
        vparameters.push_back("7.1088964505,3.60337042297,0.14561,0.23414,0.69696,0.22057,0.48276,0.79885,0.18673,0.10525,0.29124,0.24906,0.35443,0.19001,0.4586,0.3299,0.02876,0.0936,0.3899,0.1424,0.4187,0.2147,0.16739,0.1085,0.2261,0.23449,0.2999,0.2626,0.32782,0.0777,0.3681,0.7514,0.2725,0.3751,0.65898,0.0556,0.3427,0.02051,0.3251,0.3037,0.52047,0.4811,0.0546,0.59334,0.0026,0.0009,0.31632,0.3813,0.3334,0.82115,0.2817,0.0576,0.71776,0.168,0.0503,-0.08181,0.26434,0.22687,0.42638,0.02701,0.34822,0.57318,0.15408,0.39584,-0.07135,0.37469,0.25769,0.06464");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tI20_79_c_2a_c"){
        vparameters.push_back("8.6489709127,0.842525147414,0.0,0.4896,0.337,0.164,0.2196,0.1519,0.1578,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI48_80_2b_4b"){
        vparameters.push_back("9.693,0.617455896007,0.2621,0.5076,0.0299,0.2455,0.4909,0.4804,0.3974,0.1497,0.0077,0.1102,0.3642,-0.0098,0.6086,0.3609,0.5064,0.65,0.1038,0.2484");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tP12_81_adg_2h"){
        vparameters.push_back("5.3391020438,1.87980670176,0.25,0.2739,0.234,0.128,0.2289,0.23,0.6373");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tI32_82_3g_g"){
        vparameters.push_back("8.954,0.495711413893,0.0775,0.1117,0.2391,0.3649,0.0321,0.9765,0.1689,0.22,0.7524,0.2862,0.0487,0.4807");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_tP10_83_adk_j"){
        vparameters.push_back("6.2840064744,0.638128580522,0.375,0.191,0.109,0.314");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP30_85_ab2g_cg"){
        vparameters.push_back("11.6179902668,0.614391461525,0.6517,0.5428,0.6612,0.5963,0.6531,0.541,0.1258,0.5856,0.1045,0.2524");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tP32_86_g_3g"){
        vparameters.push_back("9.997999333,0.499899979999,0.0439,0.20812,0.5354,0.11009,0.22151,0.0295,0.14275,0.66613,0.7153,0.53342,0.06957,0.7593");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_tI20_88_f_a"){
        vparameters.push_back("6.407994829,2.01685393259,0.147,0.017,0.298");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI96_88_2f_4f"){
        vparameters.push_back("13.66,0.436603221083,0.1155,0.1249,0.4746,0.1356,0.125,0.0267,-0.0134,0.1262,-0.0046,-0.0251,0.1252,0.5,0.2739,0.1245,-0.0002,0.2631,0.1241,0.5043");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A17BC4D_tP184_89_17p_p_4p_io"){
        vparameters.push_back("16.2503,0.81957871547,0.86772,0.10962,0.1709,0.485,0.7893,0.2334,0.4627,0.7103,0.2146,0.4383,0.6125,0.2885,0.416,0.5631,0.354,0.4242,0.6323,0.3206,0.4579,0.7209,0.234,0.226,0.661,0.3111,0.2275,0.696,0.3148,0.2616,0.794,0.236,0.287,0.818,0.179,0.265,0.733,-0.0696,0.3435,-0.084,-0.078,0.2512,-0.097,-0.004,0.378,0.549,-0.01,0.309,0.593,-0.003,0.238,0.558,-0.006,0.165,0.6,0.2676,0.3463,0.6943,0.1965,0.4869,0.8803,0.096,0.4889,0.7614,0.1143,0.375,0.0177,-0.0146,0.3768,0.8623");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C13D_tP40_90_g_d_cef2g_c"){
        vparameters.push_back("7.3739946979,1.45226471387,-0.0654,0.7769,0.83072,0.6476,0.7846,0.648,0.751,0.001,0.7502,0.5354,0.32757,0.7289,0.8811,0.2635");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C17D4E_tP54_90_a_g_c4g_g_c"){
        vparameters.push_back("9.5601062765,0.748953974895,0.2035,-0.023,0.7639,0.5152,0.407,0.6525,0.8664,0.4465,0.8786,0.838,0.2708,0.6666,0.8964,0.0943,0.6838,0.656,0.2468,0.7207,0.8114,0.2609");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP24_91_d_d_d"){
        vparameters.push_back("3.7620082462,6.71079213192,0.303,0.202,0.019,0.296,0.189,0.08,0.2975,0.1983,0.1795");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB32CD4E8_tP184_93_i_16p_af_2p_4p"){
        vparameters.push_back("11.461,2.92112381119,0.62279,0.369,0.14,0.5403,0.251,0.129,0.6173,0.138,0.647,0.7882,0.215,0.58,0.8645,0.138,0.517,0.591,0.132,0.58,0.5541,0.244,0.605,0.536,0.345,0.554,0.5495,0.343,0.485,0.5842,0.237,0.463,0.6031,0.008,0.366,0.6554,0.077,0.371,0.6892,0.082,0.273,0.7151,0.023,0.17,0.7031,-0.04,0.165,0.6678,-0.05,0.265,0.6408,0.2272,0.1019,0.5636,0.2601,0.5802,0.812,0.1021,0.205,0.5436,0.1936,-0.0676,0.5558,0.404,0.6758,0.8049,0.2809,0.418,0.7912");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A14B3C5_tP44_94_c3g_ad_bg"){
        vparameters.push_back("7.3451315192,1.41592920353,0.8225,-0.0016,0.7086,-0.0024,0.6281,-0.0287,0.6575,0.6242,0.54,0.75,0.5448,0.7312,0.7348,0.7403");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C_tP18_94_eg_c_a"){
        vparameters.push_back("4.6863209101,1.96124874637,0.6623,0.7093,0.684,0.707,0.6579");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP24_95_d_d_d"){
        vparameters.push_back("3.7620082462,6.71079213192,0.303,0.202,-0.019,0.296,0.189,-0.08,0.2975,0.1983,0.8205");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8CD_tI24_97_d_k_a_b"){
        vparameters.push_back("5.4068544677,1.92010356944,0.1697,0.3128,0.1237");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C2_tI44_97_e_2k_cd"){
        vparameters.push_back("9.5317026755,1.33879580768,0.1553,0.0449,0.284,0.3693,0.312,0.1212,0.1191");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI12_98_f_a"){
        vparameters.push_back("7.953376649,0.587954231111,0.44");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8C2D_tP26_100_c_abcd_c_a"){
        vparameters.push_back("8.527,0.611047261639,0.7904,0.4646,0.3707,0.32701,0.0,0.1259,0.7949,0.1282,0.4871,0.2924,0.5772,0.3571");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B11C6_tP40_100_ac_bc2d_cd"){
        vparameters.push_back("10.1311200179,0.47843253381,0.0,0.0672,0.18111,0.0154,0.6536,0.6912,0.6174,0.0432,0.5814,0.678,0.6402,0.7288,0.574,0.1742,0.7098,0.5782,0.5334");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B7C2_tP32_101_bde_ade_d"){
        vparameters.push_back("9.8510697809,0.697188102729,0.0,0.4692,0.17136,0.73945,0.2254,0.3402,0.30926,0.0086,0.2384,0.5244,0.2259,0.0352,0.3449,0.0281");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP20_102_2c_b2c"){
        vparameters.push_back("8.3289849893,0.909833113223,0.5,0.604,0.439,0.623,0.031,0.795,0.725,0.848,0.251");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_tP10_103_a_d"){
        vparameters.push_back("6.5509768136,1.04518394138,0.0,0.144,0.3276,0.242");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B5C4_tP28_104_ac_ac_c"){
        vparameters.push_back("10.6225961282,0.84830508475,0.5,0.8821,0.8116,0.6057,0.3261,0.60942,0.80921,0.00978,0.8116,-0.072,0.1681");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C4_tP22_104_a_2ac_c"){
        vparameters.push_back("9.3940153509,0.981690440703,0.786,0.5,0.0649,0.8297,0.6458,0.286,0.6491,0.8588,0.036");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tP20_105_f_ac_2e"){
        vparameters.push_back("7.7858653925,1.11276650398,0.0,0.0241,0.3351,0.3664,0.1632,0.6068,0.3453,0.229,0.2558");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC3D_tP64_106_3c_c_3c_c"){
        vparameters.push_back("10.8312004139,0.492105992062,0.836,0.614,0.2,0.531,0.845,0.11,0.858,0.825,0.08,0.5929,0.62036,0.0343,0.8146,0.6194,0.0688,0.5681,0.8298,0.1322,0.8669,-0.0994,0.0636,0.69285,-0.05677,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B7_tI24_107_ac_abd"){
        vparameters.push_back("7.6400197048,0.760471204184,0.0,0.056,0.04,0.22,0.0,0.243,0.29");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI4_107_a_a"){
        vparameters.push_back("3.5440505103,1.57477426639,0.0,0.427");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_tI32_108_ac_a2c"){
        vparameters.push_back("8.0549870847,1.94761018001,0.007,0.75,0.109,0.257,0.676,0.114,0.676,0.4");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tI12_109_a_a_a"){
        vparameters.push_back("4.2490694941,3.42174629325,0.081,0.666,0.5");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI8_109_a_a"){
        vparameters.push_back("3.4517145504,3.38383984705,0.416,0.0"); //DX 20181220 - changed from "0.5416,0.5" to "0.416,0.0"
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC8_tI176_110_2b_b_8b"){
        vparameters.push_back("13.62003,0.668135092213,0.5491,0.8362,0.259,0.8065,0.6386,0.3704,0.1998,0.0869,0.0,0.5687,-0.0955,0.3098,0.6012,0.7807,0.2082,0.5033,0.7931,0.3482,0.6413,0.5083,0.4159,0.8436,0.6054,0.2752,0.8374,0.7122,0.3904,0.7286,0.645,0.3515,0.815,0.5899,0.4761");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP12_111_2n_adf"){
        vparameters.push_back("5.1219931862,1.02616165562,0.205,0.28,0.301,0.622");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP8_111_n_n"){
        vparameters.push_back("4.1305540686,0.998959140196,0.2522,0.7473,0.24415,0.24404");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tP12_112_b_n_e"){
        vparameters.push_back("5.4410982776,1.85862633019,0.2334,0.2761,0.6295");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC7D2_tP24_113_e_a_cef_e"){
        vparameters.push_back("7.8338,0.639306594501,0.8201,0.8324,0.4935,0.6407,0.7471,0.6396,0.0642,0.0798,0.1862,0.7856");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP32_114_3e_e"){
        vparameters.push_back("9.6362543341,0.547945205485,0.6,0.743,0.809,0.836,0.555,0.604,0.881,0.899,0.844,0.5125,0.7273,0.563");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_tP10_114_e_a"){
        vparameters.push_back("5.2323591487,1.07923706139,0.626,0.768,0.846");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP5_115_g_ag"){
        vparameters.push_back("3.3269188443,1.84881274424,0.253,0.6308");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tP12_115_j_egi"){
        vparameters.push_back("8.7863400452,0.701854022742,0.0288,0.0315,0.26398,0.24918,0.24945");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP20_116_bci_fj"){
        vparameters.push_back("6.1720115185,1.606448477,0.177,0.625,0.655,0.216,0.582");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP20_117_i_adgh"){
        vparameters.push_back("7.7289660931,0.727131582349,0.73,0.73,0.75,0.52,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP16_118_ei_f"){
        vparameters.push_back("6.9983025398,1.03510852635,0.237,0.15,0.343,0.149,0.509");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_tP32_118_g2i_aceh"){
        vparameters.push_back("5.8229835854,2.43860552978,0.6709,0.675,0.5861,0.73,0.85,0.5515,0.84,0.8,0.15");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tI24_119_b2i_af"){
        vparameters.push_back("6.315,2.37529691211,0.372,0.2068,0.2229,0.3067,0.3917");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI4_119_c_a"){
        vparameters.push_back("5.4790101504,0.558496075934");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC2_tI28_120_i_d_e"){
        vparameters.push_back("8.8470588481,0.924381146154,0.856,0.6452,0.6575,0.0851");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC4D_tP10_123_gh_a_i_d"){
        vparameters.push_back("3.8757,3.38106664603,0.3336,0.1193,0.2246");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tP12_124_a_m_c"){
        vparameters.push_back("6.1884808485,0.816448537725,0.162,0.662");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_tP10_124_a_m"){
        vparameters.push_back("6.4989671731,1.05200800123,0.1425,0.3361");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_tP10_125_m_a"){
        vparameters.push_back("6.6398746049,0.899096385539,0.425,0.255");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_tP12_125_a_b_m"){
        vparameters.push_back("6.3759876428,1.30630489336,0.3822,0.2163");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tP28_126_cd_e_k"){
        vparameters.push_back("7.4920241479,1.58609183128,0.11808,0.0865,0.5924,0.1254");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_tP20_127_ehj_g"){
        vparameters.push_back("7.256,0.56684123484,0.2,0.31,0.1,0.2,0.04");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C_tP18_128_eh_d_a"){
        vparameters.push_back("7.057532571,1.41383170154,0.2523,0.2217,0.2511");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C_tP40_128_egi_h_e"){
        vparameters.push_back("6.336,2.34690656566,0.366,0.2008,0.165,0.278,0.088,0.198,0.42,0.1");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tP28_130_f_c_g"){
        vparameters.push_back("8.5103337343,0.683196239716,0.58,0.5815,0.045,0.136,0.597");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_tP32_130_cg_cf"){
        vparameters.push_back("8.465,1.94329592439,0.2271,0.0095,0.1482,0.57997,0.07997,0.10688");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C4D_tP18_132_e_i_o_d"){
        vparameters.push_back("5.6046,2.34700067801,0.2369,0.26316,0.34803");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C_tP16_132_d_io_a"){
        vparameters.push_back("5.4229923801,1.46597824081,0.3,0.2,0.333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tP32_133_h_i2j"){
        vparameters.push_back("9.3810096033,0.497068542799,-0.0329,0.8986,0.658,0.0472");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP24_135_gh_h"){
        vparameters.push_back("8.3218,0.607332548247,0.36248,-0.05789,0.17358,0.13396,0.20929");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_tP28_135_gh_h_d"){
        vparameters.push_back("8.4909023894,0.697208809336,0.169,0.114,0.386,0.167,0.175");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP40_137_cdf_3g"){
        vparameters.push_back("8.097,1.41410398913,0.0,0.011,0.511,0.533,0.147,0.467,0.864,0.5,0.603");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP6_137_d_a"){
        vparameters.push_back("3.64008007,1.4478021978,0.565");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC4_tP18_137_g_b_g"){
        vparameters.push_back("5.0593101041,1.39612571654,0.08,0.1,0.503,0.384");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tP6_137_a_d"){
        vparameters.push_back("4.3675,2.8551803091,0.389");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP12_138_bi"){
        vparameters.push_back("3.388,1.77420306966,0.086,0.107");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI8_139_e_e"){
        vparameters.push_back("4.4795,2.43451278044,0.3356,0.119");
        vparameters.push_back("1.0,2.82842712475,0.125,0.625"); // Lederer-47
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_tI32_140_ah_bk"){
        vparameters.push_back("9.64,0.515560165975,0.17,0.074,0.223");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_tI32_140_ah_cl"){
        vparameters.push_back("5.46,1.91575091575,0.625,0.166,0.15");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI12_141_e_a"){
        vparameters.push_back("4.126,3.47697527872,0.2915");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tI16_142_f"){
        vparameters.push_back("8.5939,0.420984651904,0.1405");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B14C3_hP21_143_bd_ac4d_d"){
        vparameters.push_back("7.3813879247,0.611841213925,0.0,0.305,0.497,0.453,0.09,0.302,0.087,0.605,0.536,0.253,0.059,0.583,0.141,0.429,0.08,0.444,0.296,0.071,0.2215,0.2757,0.806");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B6C_hP11_143_bd_2d_a"){
        vparameters.push_back("6.9672687469,0.526119402977,0.0004,-0.0007,0.8181,0.1915,0.4998,0.5475,0.48,0.0003,0.1859,0.799,0.5002");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP12_143_cd_ab2d"){
        vparameters.push_back("6.5003859369,0.944615384626,0.0,0.5,0.25,0.0548,0.2679,0.25,0.33333,0.16667,0.5,0.0,0.5,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_hP15_144_4a_a"){
        vparameters.push_back("6.2151204722,1.25245374096,0.4904,0.2194,0.2268,0.2226,0.4873,0.1142,0.0775,0.0012,0.0,0.6097,0.0014,0.0018,0.3178,0.0008,0.5062");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_144_a_a"){
        vparameters.push_back("4.0452497415,2.30951792334,0.33,0.16,0.197,0.13,0.32,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C3DE7_hP48_145_2a_3a_3a_a_7a"){
        vparameters.push_back("6.7260126433,2.23669342849,0.567,0.317,0.168,0.432,0.752,0.169,0.6161,-0.0009,0.0,-0.0011,0.6284,0.0045,0.3706,0.371,-0.0043,-0.003,0.266,0.01,0.271,0.002,-0.005,0.73,0.734,-0.02,0.001,0.297,0.165,0.664,0.348,0.09,0.668,0.333,0.238,0.353,0.273,0.167,0.321,0.654,-0.092,0.336,0.666,0.761,0.646,-0.079,0.165,0.203,0.202,0.831");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_hR5_146_b_a_a"){
        vparameters.push_back("6.8858547087,1.22475238897,0.47,0.0,0.49,-0.002,-0.14399");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hR10_146_2a_2a_2b"){
        vparameters.push_back("6.2648898445,3.16041500397,0.3911,0.7256,0.1132,0.0,0.1454,-0.1886,0.474,0.6151,-0.0128,0.3247");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_hR42_148_2f_4f_f"){
        vparameters.push_back("12.4401,0.661634552777,0.20546,-0.56835,0.60961,0.87265,0.10211,-0.72078,0.66221,-0.37042,-0.03997,0.41696,-0.2507,0.0834,0.70548,0.00275,0.03736,0.37614,-0.32772,0.70785,0.53815,-0.23404,-0.05466");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hR18_148_2f_f"){
        vparameters.push_back("13.0471,0.659280606418,0.7407,-0.76315,0.0238,1.14609,0.38403,-0.59217,0.08714,0.06386,0.3263");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP24_149_acgi_3l"){
        vparameters.push_back("5.1418156677,2.78268310709,0.33333,0.33333,0.0,0.33333,0.421,0.33333,0.33333,0.246,0.0,0.33333,0.088");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP24_153_3c_2b"){
        vparameters.push_back("6.017,2.87518697025,0.1111,0.4444,0.1111,0.2222,0.09357,0.4444,0.8888,0.09357,0.77778,0.55558,0.09357");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP9_154_bc"){
        vparameters.push_back("6.9082,0.616557134999,0.876,0.23,0.534,0.051");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP9_156_b2c_3a2bc"){
        vparameters.push_back("4.239807232,4.83608490569,0.0,0.66667,0.33333,0.75,0.16667,0.5,0.08333,0.41667,0.83333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP12_156_2ab3c_2ab3c"){
        vparameters.push_back("4.2499813346,4.90823529409,0.375,0.70833,0.5,0.83333,0.04167,0.16667,0.45833,0.79167,0.125,0.33333,0.66667,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_156_ab_ab"){
        vparameters.push_back("4.2794836776,1.67515774714,0.636,0.0,0.894,0.5");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B6C2_hP13_157_2ac_2c_b"){
        vparameters.push_back("5.939,1.08233709379,0.264,0.736,0.022,0.5,0.522,0.603,0.215,0.397,0.829");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_158_d_a"){
        vparameters.push_back("6.12,0.924509803922,0.0,0.318,0.027,0.237");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP20_159_bc_2c"){
        vparameters.push_back("7.7488577892,0.813266227893,0.0,0.337,0.154,0.021,0.451,0.061,0.287,0.149,0.282,0.083");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_hP28_159_ab2c_2c"){
        vparameters.push_back("7.7479899913,0.724961280333,0.0,0.3649,0.0424,0.3891,0.0408,0.3169,0.3198,0.2712,0.0821,0.5089,0.3172,0.1712,0.2563,0.0274");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C7D_hP26_159_b_ac_a2c_b"){
        vparameters.push_back("6.2653109789,1.6324735851,0.0,0.3057,0.0613,0.4379,0.3425,0.1575,0.2472,0.0015,0.5149,0.3102,0.3347,0.1157,0.0618");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hR4_160_b_a"){
        vparameters.push_back("4.405,0.610442678774,0.0,0.52073,-0.02217");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B5_hR26_160_a3bc_a3b"){
        vparameters.push_back("12.71838,0.624030733474,0.672,0.194,0.654,0.01201,0.349,0.58199,0.722,0.35601,1.003,-0.20599,0.998,-0.66,0.355,1.00601,1.033,-0.339,0.28799");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hR3_160_a_a_a"){
        vparameters.push_back("6.1703,0.949908432329,0.0,0.79356,0.25763");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR10_160_5a_5a"){
        vparameters.push_back("3.09,12.2653721683,0.0,0.13333,0.4,0.6,0.86667,0.05,0.18333,0.45,0.65,-0.08333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP5_164_d_ad"){
        vparameters.push_back("3.9381,1.55813717275,0.2467,0.647");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP9_164_bd_c2d"){
        vparameters.push_back("2.89,7.90657439446,0.0607,0.154,0.27263,0.39403");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP4_164_a_b_d"){
        vparameters.push_back("4.0482,1.26777333136,0.271");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP24_165_bdg_f"){
        vparameters.push_back("7.07,1.00919377652,0.17,0.38,0.69,0.07,0.08");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_hR7_166_2c_ac"){
        vparameters.push_back("3.335,7.48635682159,0.29422,0.12967,0.2168");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hR6_166_c_c_c"){
        vparameters.push_back("3.8548,7.95397945419,0.1159,0.3017,0.3815");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hR10_167_b_e_a"){
        vparameters.push_back("5.4577,2.40134122433,0.6946");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hR24_167_e_e_2e"){
        vparameters.push_back("12.76,0.575235109718,1.1389,0.8113,1.0343,0.3584");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B13C4_hP57_168_d_c6d_2d"){
        vparameters.push_back("15.9361426085,0.244226907628,0.0,0.4965,0.179,0.555,0.035,0.187,0.007,0.349,0.034,0.018,0.168,0.388,0.01,0.195,0.569,0.02,0.093,0.454,0.539,0.273,0.095,0.555,0.2634,0.09,0.086,0.0789,0.4368,0.071");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_hP72_168_2d_8d_2d"){
        vparameters.push_back("13.759939558,0.609738372094,0.4498,0.113,0.6253,0.1278,0.466,0.1246,0.4154,0.2053,0.0456,0.2066,0.4189,0.5576,0.4218,0.0894,0.8248,0.1514,0.4907,0.3249,0.3746,0.0106,0.0948,0.0083,0.3667,0.4879,0.5693,0.1616,0.0357,0.1505,0.5625,0.5943,0.4451,0.1173,0.0,0.1291,0.4589,0.4993");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP30_169_2a_3a"){
        vparameters.push_back("6.4300116544,2.78071539658,0.013,0.3579,0.12736,0.3339,0.3226,0.29886,0.3347,0.0,0.004,0.0119,0.3343,0.0,0.338,0.0064,0.33823");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP30_170_2a_3a"){
        vparameters.push_back("6.4300116544,2.78071539658,0.013,0.3579,0.87264,0.3339,0.3226,0.70114,0.3347,0.0,-0.004,0.0119,0.3343,0.0,0.338,0.0064,0.66177");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B2C_hP39_171_5c_c_a"){
        vparameters.push_back("6.3199906634,3.05221518986,0.0,-0.004,0.253,0.23633,0.015,0.248,0.11034,0.681,0.192,0.16733,0.448,0.188,0.29533,0.449,0.25,0.048,0.373,0.065,0.50467");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B2C_hP39_172_5c_c_a"){
        vparameters.push_back("6.3199906634,3.05221518986,0.0,-0.004,0.253,0.76367,0.015,0.248,0.88966,0.681,0.192,0.83267,0.448,0.188,0.70467,0.449,0.25,-0.048,0.373,0.065,0.49533");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_173_c_b"){
        vparameters.push_back("7.1329719992,1.03939436422,0.0,0.0337,0.3475,0.146");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_hP14_173_bc_c"){
        vparameters.push_back("7.603038022,0.382612126795,0.0,0.3284,0.0313,0.05,0.2314,0.4063,0.013");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B7C2_hP21_174_2j2k_ajk_cf"){
        vparameters.push_back("9.0004021308,0.399102242177,0.4309,0.3719,0.1189,0.2772,0.4163,0.1204,0.0495,0.4359,0.2232,0.124,0.2889,0.4096");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP12_174_cj_fk_aj"){
        vparameters.push_back("10.7303215747,0.395153774459,0.30167,0.15433,0.03467,0.51733,0.14967,0.31433");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B7C6_hP21_175_ck_aj_k"){
        vparameters.push_back("9.5057625379,0.329093950194,0.36335,0.08627,0.0661,0.221,0.15405,0.51373");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP36_175_jk_jk_jk"){
        vparameters.push_back("11.5799622371,0.317895264089,0.255,0.0598,0.1323,0.416,0.3483,0.0762,0.2334,0.5727,0.4279,0.0715,0.1367,0.5046");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP10_176_h_bc"){
        vparameters.push_back("7.8700439404,0.50025412961,0.3847,0.0915");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C_hP14_176_h_h_c"){
        vparameters.push_back("9.3500327107,0.451336898402,0.1493,0.1701,0.357,0.0462");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_176_h_c"){
        vparameters.push_back("7.4429335392,0.580545478976,0.375,0.083");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP36_177_j2lm_n"){
        vparameters.push_back("12.7835,0.291064262526,0.61855,0.39242,0.79257,0.44445,0.52169,0.86952,0.16458");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP24_178_b_ac"){
        vparameters.push_back("5.14898393,3.15789473684,0.8361,0.7601,0.5338,0.3099,-0.0053");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP6_178_a"){
        vparameters.push_back("2.355,4.43566878981,0.461");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP24_179_b_ac"){
        vparameters.push_back("5.14898393,3.15789473684,0.8361,0.7601,0.5338,0.3099,0.0053");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_181_j_c"){
        vparameters.push_back("4.9977,1.09252256038,0.2072");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP3_183_a_a_a"){
        vparameters.push_back("3.3908401495,1.49568037743,0.608,0.0,0.226");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_183_c_ab"){
        vparameters.push_back("5.3175214551,0.83247442072,0.0,0.513,0.01");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_hP72_184_d_4d_d"){
        vparameters.push_back("13.7178848276,0.616168537688,0.45652,0.12053,0.0,0.42,0.2069,0.071,0.1224,0.4519,0.2982,0.3629,0.0019,0.0649,0.1538,0.5737,0.0602,0.12298,0.4525,0.12746");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_hP30_185_cd_c_ab"){
        vparameters.push_back("11.795090915,0.502416278086,0.0,0.377,0.1598,0.2396,0.6647,0.1706,0.1732,0.5056,0.1148");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP24_185_ab2c_c"){
        vparameters.push_back("6.9593,1.02639633296,0.3213,0.1998,0.2806,0.0765,0.3761,0.4246,0.3322,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_185_c_a"){
        vparameters.push_back("6.1197165237,0.924509803925,0.0,0.305,0.245");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP24_185_c_ab2c"){
        vparameters.push_back("8.7838,1.02449964708,0.2684,0.2311,0.3321,0.25,0.3153,0.5863,0.3518,-0.0769");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B7_hP20_186_c_b2c"){
        vparameters.push_back("9.85,0.624365482234,0.06,0.815,0.31,0.126,0.25,0.544,0.31");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP4_187_e_fh"){
        vparameters.push_back("2.8065,2.53411722786,0.198");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_hP10_188_k_c_a"){
        vparameters.push_back("7.2864263258,0.928891999832,0.0041,0.32685");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB9C4_hP28_188_e_kl_ak"){
        vparameters.push_back("6.4953629976,1.43896355826,0.07103,0.48306,0.12023,0.75436,0.22923,0.00127,0.6032");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8BC3D6_hP18_189_bfh_a_g_i"){
        vparameters.push_back("6.62,1.20241691843,0.403,0.444,0.231,0.75,0.222");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9BC3D5_hP18_189_fi_a_g_bh"){
        vparameters.push_back("6.6,1.19696969697,0.378,0.43,0.266,0.755,0.236");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP18_190_gh_bf"){
        vparameters.push_back("7.946,0.789076264787,0.0225,0.294,0.612,0.01");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_hP16_190_bdh_g"){
        vparameters.push_back("6.9236970109,1.22634969237,0.3313,0.0628,0.6682");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP24_190_i_afh"){
        vparameters.push_back("5.9699820408,1.96984924623,0.52,0.6683,0.6653,0.3786,0.3233,0.623");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C18D6_hP58_192_c_f_lm_l"){
        vparameters.push_back("9.214,0.997829390059,0.3103,0.2369,0.3876,0.1159,0.4985,0.1456,0.1453");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP72_192_m_j2kl"){
        vparameters.push_back("13.7705812079,0.608458538779,0.6373,0.7895,0.4221,0.6688,0.5445,0.1221,0.4551,0.686");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_hP16_193_dg_g"){
        vparameters.push_back("6.9104160691,0.696671490596,0.2358,0.5992");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP16_194_gh_ac"){
        vparameters.push_back("5.096,1.6295133438,-0.16667");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_hP28_194_ahk_ch"){
        vparameters.push_back("7.656,0.991771159875,0.533,0.872,0.196,0.439");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B3C_hP26_194_hk_h_a"){
        vparameters.push_back("7.513,1.03087980833,0.458,0.12,0.201,-0.067");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12BC4_cP34_195_2j_ab_2e"){
        vparameters.push_back("8.0357772599,0.2493,0.7507,0.14225,0.5,0.35579,0.0002,0.14286,0.35831");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B2C_cF60_196_h_bc_a"){
        vparameters.push_back("9.9799933334,0.0,0.0625,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_cP20_198_a_a_b"){
        vparameters.push_back("6.57,0.417,0.064,0.303,0.592,0.5");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B11_cP39_200_f_aghij"){
        vparameters.push_back("8.5520223662,0.18,0.34,0.265,0.278,0.157,0.257");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_cP60_201_be_fh_g"){
        vparameters.push_back("9.5599167841,0.3889,0.6111,0.0972,0.0389,0.0972,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B6C_cF104_202_h_h_c"){
        vparameters.push_back("10.6100296668,0.5827,0.6359,0.638,0.72");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cF240_202_h2i"){
        vparameters.push_back("14.26,0.249,0.052,0.105,0.085,0.22,0.185,0.052,0.165");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BCD3E6_cF208_203_e_c_d_f_g"){
        vparameters.push_back("13.9898,0.2826,-0.099,0.2257,0.2665,0.3531");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C6D16E_cF232_203_e_d_f_eg_a"){
        vparameters.push_back("13.9038,0.28207,0.06362,0.34379,0.26626,0.22529,0.35333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C16_cF160_203_a_bc_eg"){
        vparameters.push_back("16.6600284675,0.20516,0.01201,0.111,0.42978");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C6_cP264_205_2d_ab2c2d_6d"){
        vparameters.push_back("15.263,0.2561,0.375,0.2526,0.0133,0.0197,0.2444,0.2335,0.0046,0.1386,0.3763,0.1272,0.38,0.3838,0.1209,0.2777,0.1241,0.0103,0.4835,0.1315,0.2536,0.2664,0.2841,0.1049,0.235,0.4047,0.2921,0.3491,-0.0385,-0.1074,0.1509,-0.0104,-0.0242");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP240_205_10d"){
        vparameters.push_back("14.04078,0.2294,-0.0325,0.101,0.2467,-0.054,0.0061,0.2081,0.0646,0.1289,0.2066,0.8599,-0.036,0.171,-0.0963,0.159,0.2236,0.1122,-0.0371,0.2439,0.0192,-0.0636,0.2053,0.1349,0.0616,0.1503,0.7983,0.0202,0.1323,0.8207,0.1186");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C2_cI96_206_c_e_ad"){
        vparameters.push_back("9.46,0.115,0.205,0.16,0.382,0.11");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A17B15_cP64_207_acfk_eij"){
        vparameters.push_back("10.6058825779,0.2422,0.2622,0.3319,0.2701,0.142,0.1539,0.3498");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cP16_208_j_b"){
        vparameters.push_back("6.31,0.184");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2CD6E_cP64_208_m_ad_b_m_c"){
        vparameters.push_back("10.3701312618,0.057,0.25,0.23,0.235,0.25,0.458");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A24BC_cF104_209_j_a_b"){
        vparameters.push_back("7.7099775082,0.043,0.109,0.165");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B36CD12_cF488_210_h_3h_a_fg"){
        vparameters.push_back("16.4321054599,0.12536,0.51382,0.55423,0.58845,0.50025,0.5008,0.59,0.3557,0.5123,0.5411,0.1708,0.5347,0.2125,0.5888");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B6C_cF608_210_4h_2h_e"){
        vparameters.push_back("16.4321054599,0.3771,0.0157,0.2009,0.1224,0.5287,0.1425,0.0068,0.0949,0.2596,0.2225,0.6117,0.1785,0.0224,0.0928,0.2607,0.1616,0.0194,0.0793,0.3503");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cI72_211_hi_i"){
        vparameters.push_back("9.68882,0.37338,0.15866,0.38235");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cP12_212_c_a"){
        vparameters.push_back("6.54,0.428");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C_cI56_214_g_h_a"){
        vparameters.push_back("12.31504,0.108,0.384");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC2_cI48_214_f_a_e"){
        vparameters.push_back("10.38,0.266,0.365");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B9_cP52_215_ei_3efgi"){
        vparameters.push_back("8.7068,0.1157,0.8296,0.3253,0.6066,0.3534,0.8549,0.8113,0.5332,0.3153,0.0322");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD_cF16_216_c_d_b_a"){
        vparameters.push_back("6.465");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C_cP16_218_c_e_a"){
        vparameters.push_back("6.0258147002,0.6486");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7BC3D13_cF192_219_de_b_c_ah"){
        vparameters.push_back("12.0986,0.0808,0.0987,0.0214,0.1821");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A15B4_cI76_220_ae_c"){
        vparameters.push_back("9.718,-0.042,0.12,0.16,-0.04");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_cI28_220_c_a"){
        vparameters.push_back("8.6,0.08333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C6_cP33_221_cd_ag_fh"){
        vparameters.push_back("7.624,0.3,0.2,0.2");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3C16_cP96_222_ce_d_fi"){
        vparameters.push_back("11.1199004528,0.0,0.125,0.084,0.166,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A23B6_cF116_225_bd2f_e"){
        vparameters.push_back("12.523,0.203,0.178,0.378");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C_cF36_225_e_c_a"){
        vparameters.push_back("9.725,0.24");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB13_cF112_226_a_bi"){
        vparameters.push_back("12.2836,0.1806,0.1192");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C7_cF88_227_c_d_af"){
        vparameters.push_back("10.2663,0.4157");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_cF56_227_ad_e"){
        vparameters.push_back("8.0835,0.2642");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5BCD6_cF416_228_eg_c_b_h"){
        vparameters.push_back("22.2649775014,0.19,0.425,0.05,0.18,0.28");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B_cF224_228_h_c"){
        vparameters.push_back("15.4800297791,0.043,0.138,0.278");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B10_cI52_229_e_fh"){
        vparameters.push_back("8.9822,0.3538,0.13505,0.3045");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_cI10_229_c_a"){
        vparameters.push_back("6.186");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B3_cI40_229_df_e"){
        vparameters.push_back("8.735,0.342,0.156");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C12D3_cI160_230_a_c_h_d"){
        vparameters.push_back("11.4597,0.3471,0.4664,0.0512");
      }
    }
    // ---------------------------------------------------------------------------
    // miscellaneous structures
    // ---------------------------------------------------------------------------
    if(library=="all" || library == "" || library=="misc"){
      // SQS structures
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_aP16_2_4i_4i"){
        vparameters.push_back("5.4244498076,1.04446593575,1.78376517004,72.976133815,87.0786711125,79.9750121379,0.375,0.75,0.375,0.875,0.0,0.375,0.375,0.5,0.875,0.625,0.75,0.625,0.875,0.5,0.375,0.125,0.75,0.125,0.875,0.75,0.875,0.375,-0.0,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B11_mP16_6_2abc_2a3b3c"){
        vparameters.push_back("6.5421326204,1.0,1.0,90.0,0.5,-0.0,0.5,0.5,0.0,-0.0,0.0,0.5,0.0,0.5,0.0,-0.0,0.5,-0.0,0.5,0.5,0.25,0.25,0.25,0.75,0.25,0.25,0.75,0.25,0.75,0.25,0.25,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC32_8_4a_12a"){
        vparameters.push_back("13.8779590179,0.333333333333,0.666666666667,109.471220634,-0.0,0.5,0.5,0.75,0.375,0.9375,0.125,0.8125,0.875,0.1875,0.625,0.0625,0.75,0.375,0.875,0.6875,-0.0,-0.0,0.5,0.25,0.625,0.5625,0.75,0.875,0.25,0.125,0.375,0.4375,0.125,0.3125,0.25,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC32_8_4a_4a4b"){
        vparameters.push_back("9.2519726786,1.0,0.707106781188,90.0,0.0,0.0,0.0,0.5,0.25,0.25,0.5,0.5,0.25,0.75,0.5,0.0,0.75,0.75,0.75,0.25,0.75,0.25,0.0,0.75,0.25,0.5,0.0,0.25,0.75,0.0,0.25,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B13_oC32_38_ac_a2bcdef"){
        vparameters.push_back("6.5421326204,1.41421356237,1.41421356237,0.5,0.0,0.5,0.0,0.75,0.25,0.75,0.75,0.25,0.25,0.25,0.25,0.75,0.75,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_oC32_38_abce_abcdf"){
        vparameters.push_back("6.5421326204,1.41421356237,1.41421356237,0.5,-0.0,-0.0,0.5,0.25,0.75,0.25,0.25,0.75,0.25,0.75,0.25,0.25,0.75,-0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB7_hR16_166_c_c2h"){
        vparameters.push_back("9.2519726786,1.22474487138,0.875,0.625,1.625,0.1250000001,1.125,1.6249999999");
      }
      // kesterite
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BCD4_tI16_82_ac_b_d_g"){
        vparameters.push_back("5.427,2.00313248572,0.7434,0.256,0.6278");
      }
      // misc structures (Y. Lederer)
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC8_12_a_di"){
        vparameters.push_back("1.0,0.301511344577,0.522232967868,100.024987862,0.75,0.25");
        vparameters.push_back("5.627,0.58983473,0.88484095,72.6,0.6048,0.2603");  // 002, binary metal-nitride prototype (ICSD #34675)
        vparameters.push_back("5.328,0.60041291,0.8871997,77.38,0.3756,0.7456");  // 003, binary metal-nitride prototype (ICSD #181559)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC8_12_i_i"){
        vparameters.push_back("1.0,0.301511344581,0.522232967872,100.024987862,0.625,0.875,0.875,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_mC16_12_a_di_2i"){
        vparameters.push_back("1.0,0.301511344577,0.522232967868,100.024987862,0.75,0.25,0.125,0.375,0.625,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC16_12_i_i_adi"){
        vparameters.push_back("1.0,0.301511344581,0.522232967872,100.024987862,0.625,0.875,0.875,0.625,0.75,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oP4_47_a_ct"){
        vparameters.push_back("1.0,1.41421356239,2.0,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oP4_47_cr_a"){
        vparameters.push_back("1.0,1.41421356238,2.82842712476,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oP8_47_a_ct_egs"){
        vparameters.push_back("1.0,1.41421356239,2.0,0.25,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_oP8_47_eq_g_bdt"){
        vparameters.push_back("1.0,1.41421356238,2.82842712476,0.75,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP4_51_e_e"){
        vparameters.push_back("1.0,0.707106781182,2.0,0.125,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP8_51_e_e_2f"){
        vparameters.push_back("1.0,0.707106781182,2.0,0.125,0.625,0.375,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP4_59_a_a"){
        vparameters.push_back("1.0,1.41421356238,2.0,0.875,0.375");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP8_59_a_a_2b"){
        vparameters.push_back("1.0,1.41421356238,2.0,0.875,0.375,0.875,0.375");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_63_c_c_g"){
        vparameters.push_back("1.0,1.41421356236,0.707106781182,0.625,0.125,0.25,0.375");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3_oC12_65_a_i_cj"){
        vparameters.push_back("1.0,2.99999999998,0.707106781176,0.3333333333,0.8333333333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_oC16_65_ai_b_q"){
        vparameters.push_back("1.0,1.99999999999,0.499999999994,0.75,0.75,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_oC16_65_bj_a_eh"){
        vparameters.push_back("1.0,1.41421356236,0.707106781182,0.75,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oC16_65_a_bf_hi"){
        vparameters.push_back("1.0,1.41421356239,0.5,0.75,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oC6_65_a_i"){
        vparameters.push_back("1.0,2.99999999998,0.707106781176,0.3333333333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oC8_65_ai_b"){
        vparameters.push_back("1.0,1.99999999999,0.499999999994,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oC8_65_bj_a"){
        vparameters.push_back("1.0,1.41421356236,0.707106781182,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC8_65_i_i"){
        vparameters.push_back("1.0,1.99999999999,0.499999999994,0.875,0.375");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_65_i_i_fh"){
        vparameters.push_back("1.0,1.99999999999,0.499999999994,0.75,0.875,0.375");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3_oI12_71_a_e_df"){
        vparameters.push_back("1.0,0.471404520797,0.333333333338,0.6666666667,0.6666666667");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP2_123_a_b"){
        vparameters.push_back("1.0,2.00000000002");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP2_123_a_c"){
        vparameters.push_back("1.0,0.707106781182");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP3_123_g_a"){
        vparameters.push_back("1.0,3.00000000002,0.6666666667");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP4_123_abc_d"){
        vparameters.push_back("1.0,1.41421356238");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP4_123_ag_b"){
        vparameters.push_back("1.0,4.00000000002,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP4_123_cf_a"){
        vparameters.push_back("1.0,0.499999999994");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tP4_123_a_bh"){
        vparameters.push_back("1.0,2.82842712477,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP4_123_a_b_h"){
        vparameters.push_back("1.0,2.00000000002,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP4_123_a_c_e"){
        vparameters.push_back("1.0,0.707106781182");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP4_123_a_d_bc"){
        vparameters.push_back("1.0,1.41421356238");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP4_123_g_g"){
        vparameters.push_back("1.0,4.00000000002,0.125,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_tP6_123_g_b_ch"){
        vparameters.push_back("1.0,3.00000000002,0.1666666667,0.6666666667");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_tP8_123_abc_d_i"){
        vparameters.push_back("1.0,1.41421356238,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_tP8_123_ag_b_2h"){
        vparameters.push_back("1.0,4.00000000002,0.25,0.625,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_tP8_123_cf_a_k"){
        vparameters.push_back("1.0,0.499999999994,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_tP8_123_a_bh_cdg"){
        vparameters.push_back("1.0,2.82842712477,0.25,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP8_123_h_h_abg"){
        vparameters.push_back("1.0,4.00000000002,0.25,0.625,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP8_129_c_c_2c"){
        vparameters.push_back("1.0,2.82842712472,0.875,0.375,0.125,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tI8_139_ae_b"){
        vparameters.push_back("1.0,2.82842712475,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3_tI12_139_a_e_be"){
        vparameters.push_back("1.0,4.24264068707,0.3333333333,0.8333333333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_tI16_139_ae_b_g"){
        vparameters.push_back("1.0,2.82842712475,0.25,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_tI16_139_a_bd_ce"){
        vparameters.push_back("1.0,2.0,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tI16_139_e_e_cd"){
        vparameters.push_back("1.0,2.82842712475,0.125,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tI16_141_a_b_e"){
        vparameters.push_back("1.0,1.99999999997,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP3_164_d_a"){
        vparameters.push_back("1.0,2.44948974278,0.3333333333");
        vparameters.push_back("1.0,1.22474487139,0.3333333333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_hP6_164_d_a_bd"){
        vparameters.push_back("1.0,2.44948974278,0.3333333333,0.8333333333");
        vparameters.push_back("4.704,1.6001276,0.742,0.342");  // 002, metallic ternary prototype (ICSD #616769)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_hP6_164_d_b_ad"){
        vparameters.push_back("1.0,1.22474487139,0.8333333333,0.3333333333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hR4_166_bc_a"){
        vparameters.push_back("1.0,4.89897948557,0.25");
        vparameters.push_back("3.6459845,4.1752552,0.42402");  // 002, binary metal-nitride prototype (ICSD #1144)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hR4_166_a_bc"){
        vparameters.push_back("1.0,9.79795897119,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR4_166_c_c"){
        vparameters.push_back("1.0,9.79795897107,0.625,0.125");
        vparameters.push_back("1.0,4.89897948557,0.375,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_hR8_166_bc_a_2c"){
        vparameters.push_back("1.0,4.89897948557,0.25,0.625,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_hR8_166_a_bc_2c"){
        vparameters.push_back("1.0,9.79795897119,0.75,0.375,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hR8_166_c_c_abc"){
        vparameters.push_back("1.0,9.79795897107,0.625,0.125,0.75");
        vparameters.push_back("1.0,4.89897948557,0.375,0.875,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_cP8_221_a_c_bd"){
        vparameters.push_back("1.0");
      }
      // -------------------------------------------------------------------------
      // oxide prototypes (from R. Friedrich)
      // -------------------------------------------------------------------------
      // binaries 
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC20_12_b2i_c2i"){
        vparameters.push_back("9.7915972616,0.443468950753,0.626873661672,72.47,0.825,0.343,0.821,0.833,0.66,0.17,0.658,0.669");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_mC20_12_2i_3i"){
        vparameters.push_back("12.8749583858,0.24865727853,0.999198732944,152.527888576,0.7041,0.7946,0.65536,0.31402,0.9453,0.1098,0.3899,0.5632,0.7607,0.2566");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_mC32_12_5i_3i"){
        vparameters.push_back("9.7349670426,0.385643337648,1.40884015867,133.94703436,0.0136,0.2506,0.8276,0.8736,0.5586,0.6116,0.6158,0.9368,0.267,0.574,0.442,0.2461,0.9228,0.0529,0.2983,0.4337");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_mP16_14_3e_e"){
        vparameters.push_back("8.3369865333,0.56178500085,0.838076298071,144.831245943,0.958,0.961,0.756,0.442,0.046,0.243,0.255,0.673,0.243,0.2563,0.7319,0.9681");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_mP20_14_2e_3e"){
        vparameters.push_back("7.592734201,1.09598618914,0.998035732622,133.736291174,0.524,0.816,0.158,0.04,0.959,0.263,0.778,0.702,0.071,0.235,0.948,0.11,0.27,0.969,0.758");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC8_15_a_e"){
        vparameters.push_back("4.6509294963,0.729773531807,1.35357481762,127.166266704,0.333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_mC20_15_2f_e"){
        vparameters.push_back("7.7087709852,0.481394604967,0.920353982304,63.4,0.7399,0.379,0.462,0.19,0.885,0.524,0.425");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5_oP28_19_2a_5a"){
        vparameters.push_back("4.4859196717,1.82663207955,1.86900129701,0.3784,0.3475,0.4015,0.384,0.9674,0.7836,0.0951,0.3166,0.5403,0.8671,0.6061,0.7453,0.3371,0.4845,0.8324,0.7112,0.2743,0.749,0.55,0.4939,0.5257");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP24_33_4a_2a"){
        vparameters.push_back("5.6320569384,0.884841795441,2.16335540838,0.6753,0.1608,0.593,0.6454,0.8379,0.906,0.834,0.699,0.6957,0.9256,0.188,0.8093,0.9786,0.9654,0.0,0.6272,0.9994,0.7473");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oC16_40_b_a2b"){
        vparameters.push_back("5.743,1.48999,0.833885,0.3841,0.90324,0.5,0.7323,0.3755,0.8922,0.8284");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_oP14_59_a2e_e"){
        vparameters.push_back("3.2077,3.2300714,1.2255822,0.999,0.3957,0.469,0.5689,0.997,0.39882,0.1083");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_62_c_c"){ //DX20200703 - there are other structures with this label in part1 and part2
        vparameters.push_back("6.9613,0.83485843,0.5325873,0.1136,0.7544,0.3592,0.4045");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC16_64_e_f"){
        vparameters.push_back("8.7597573116,1.03356994972,0.940698018333,0.6554,0.9104,0.0725");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_tP28_135_gh_dh"){
        vparameters.push_back("8.7553254034,0.751693002255,0.33,0.4,0.1,0.14,0.165");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP15_144_2a_3a"){
        vparameters.push_back("3.9844538913,1.92300578035,0.15,0.61,0.02,0.77,0.18,0.26,0.15,0.95,0.0,0.79,0.33,0.07,0.23,0.72,0.56");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hR3_166_c_a"){
        vparameters.push_back("5.73015,4.46193,0.744");
        vparameters.push_back("3.6378513,5.1624334,0.268");  // 002, binary metal-nitride prototype (ICSD #22231)
        vparameters.push_back("4.0290164,5.5658926,0.26647");  // 003, binary metal-nitride prototype (ICSD #409851)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hR6_166_c_2c"){
        vparameters.push_back("3.5572963728,10.7622298067,0.083,0.377,0.21");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP12_189_fg_eh"){
        vparameters.push_back("7.5487706231,0.71987757732,0.668,0.295,0.632,0.168");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP8_194_ac_f"){
        vparameters.push_back("3.7639841644,2.4272070374,0.6497");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_cI80_199_a2b_2c"){
        vparameters.push_back("10.7844999416,0.5,0.229,0.708,0.625,0.125,0.875,0.625,0.875,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_cF80_227_f_e"){
        vparameters.push_back("10.7531,0.24,0.895");
      }
      // -------------------------------------------------------------------------
      // ternaries
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B4C_aP18_2_4i_4i_i"){
        vparameters.push_back("6.63976352,1.00001361035,1.50518591105,71.844867398,80.9249357443,67.4411683666,0.7435,0.7649,0.0197,0.2323,0.4828,0.5617,0.2681,0.0455,0.3887,0.236,0.5962,0.1682,0.9572,0.3436,0.1926,0.9557,0.8502,0.3514,0.6442,0.1527,0.0884,0.5375,0.2907,0.3741,0.7148,0.2162,0.2479");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C2_aP22_2_2i_7i_2i"){
        vparameters.push_back("7.3384982604,0.986316404633,0.950710549193,96.331,116.158,86.387,0.82535,0.34153,0.25581,0.23479,0.93732,0.27034,0.5623,0.6988,0.9539,0.6304,0.3997,0.4639,0.6025,0.9461,0.341,0.9238,0.6849,0.3595,0.8578,0.0087,0.1317,0.173,0.2822,0.2504,0.232,0.6233,0.107,0.34536,0.44585,0.27538,0.73214,0.83509,0.19625");
        vparameters.push_back("10.4825,0.46040544,0.50745528,81.42,79.6425,102.8254,0.9163,0.4905,0.7907,0.58313,0.5771,0.21504,0.2925,0.1055,0.6186,0.0477,0.7911,0.917,0.2349,0.8135,0.2409,0.0666,0.3044,0.5772,0.5412,0.7648,0.8863,0.7877,0.6364,0.0493,0.5614,0.2891,0.553,0.36066,0.96292,0.34888,0.1386,0.1407,0.786");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_aP30_2_3i_9i_3i"){
        vparameters.push_back("8.440104917,0.923591309392,0.891430518055,90.055,95.217,103.426,0.19831,0.42266,0.7606,0.20241,0.92919,0.76401,0.50333,0.7504,0.52691,0.3034,0.4616,0.4628,0.3014,0.9385,0.4641,0.5705,0.7688,0.1988,0.9832,0.3739,0.2655,0.9819,0.8677,0.2648,0.4018,0.7266,0.8296,0.2183,0.1785,0.2254,0.2713,0.8704,0.0938,0.2735,0.5126,0.0931,0.1851,0.3875,0.2684,0.1849,0.9542,0.2691,0.3973,0.7236,0.0561");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C_aP32_2_4i_10i_2i"){
        vparameters.push_back("7.3087649311,1.14261433201,0.821778898093,90.89,101.58,105.99,0.3463,0.7028,0.4634,0.2875,0.6959,0.9481,0.0961,0.3757,0.6278,0.1278,0.9318,0.1815,0.1019,0.1479,0.1134,0.1351,0.6991,0.197,0.2641,0.4568,0.9447,0.2832,0.9314,0.9431,0.1337,0.6232,0.6318,0.2765,0.4421,0.4144,0.2948,0.9435,0.4818,0.4953,0.2777,0.2314,0.1057,0.1465,0.6744,0.5003,0.2278,0.7656,0.2959,0.0562,0.72,0.2808,0.3359,0.1702");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_mP28_4_4a_8a_2a"){
        vparameters.push_back("9.0939150001,1.04504737499,0.608139562228,92.875,0.807,0.829,0.272,0.204,0.832,0.263,0.302,0.661,0.774,0.321,0.171,0.214,0.7387,0.1707,0.568,0.2717,0.3202,0.4218,0.6785,0.4914,0.642,0.7422,0.9861,0.0989,0.8129,0.3034,0.0693,0.7893,0.6712,0.078,0.5051,0.2076,0.129,0.0106,0.8783,0.3432,0.5032,-0.0,0.747,0.975,0.9939,0.7889");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_mC60_5_ab8c_ab2c_3c"){
        vparameters.push_back("13.2962734468,0.579323216751,0.88233384728,68.42,0.126,0.583,0.395,0.925,0.458,0.177,0.139,0.393,0.968,0.991,0.37,0.658,0.858,0.777,0.165,0.139,0.451,0.315,0.635,0.89,0.089,0.511,0.381,0.841,0.35,0.287,0.857,0.63,0.384,0.71,0.999,0.37,0.798,0.498,0.9127,0.0,0.2494,0.2545,0.022,0.2458,0.5872,0.0207,0.2516");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5C_mC54_8_3a3b_9a3b_3a"){
        vparameters.push_back("13.993228409,0.578095627303,0.95071704775,134.2304365,0.0445,0.3334,0.9668,0.6832,0.5007,0.0005,0.2994,0.4493,0.9984,0.4969,0.6404,0.2828,0.8665,0.803,0.7033,0.5434,0.7612,0.1446,0.2517,0.8741,0.3501,0.7145,0.9192,0.057,0.0021,0.9993,0.5651,0.3477,0.4285,0.6427,0.1992,0.2336,0.6793,0.3112,0.2665,0.3313,0.7507,0.2478,0.0016,0.4657,0.8234,0.2857,0.5968,0.6856,0.0533,0.035,0.3071,0.7119");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C3_mP24_11_2e_7e_3e"){
        vparameters.push_back("9.7906148128,0.416420361249,0.938259441704,101.57,0.595,0.682,0.154,0.508,0.195,0.221,0.473,0.14,0.645,0.438,0.885,0.314,0.745,0.997,0.313,0.791,0.031,0.905,0.2806,0.0278,0.673,0.2467,0.9811,0.142");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C2_mC18_12_a_3i_i"){
        vparameters.push_back("9.7742615012,0.36510934394,0.69960238568,75.2,0.6434,0.8818,0.7535,0.3934,0.9126,0.7219,0.7658,0.6691");
        vparameters.push_back("9.1042,0.37740823,0.9940467,137.3359,0.4616,0.1133,0.3157,0.2764,0.8722,0.5652,0.53846,0.3436");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_mC48_12_gi_hi_2i3j"){
        vparameters.push_back("9.7995229272,0.904117589797,0.6838314027,73.04,0.3216,0.249,0.2996,0.3569,0.2291,0.9043,0.3587,0.0391,0.2983,0.6449,0.5415,0.8467,0.696,0.3587,0.6439,0.6088,0.6337,0.6552,0.9717");
        vparameters.push_back("8.8789,0.912805,0.986699,133.694,0.7148,0.1906,0.3547,0.6391,0.8447,0.1555,0.09,0.4376,0.5669,0.9156,0.73601,0.8485,0.23951,0.4244,0.3344,0.0772,0.4173,0.8466,0.5815");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_mP12_13_f_2g_e"){
        vparameters.push_back("6.912898268,0.837946008451,0.729410860817,136.193216009,0.1825,0.327,0.22,0.894,0.281,0.256,0.622,0.855");
        vparameters.push_back("6.5113,0.84124522,0.72976211,136.867,0.82,0.653,0.78,0.11,0.82,0.74,0.38,0.44");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_mP40_14_2e_6e_2e"){
        vparameters.push_back("9.8340572172,0.944454808597,0.553242905653,103.024236752,0.2509,0.3482,0.0296,0.2556,0.9865,0.037,0.8687,0.659,0.6887,0.1227,0.5,0.7987,0.1052,0.724,0.4962,0.3748,0.159,0.2478,0.6341,0.016,0.2471,0.6018,0.303,0.1428,0.0441,0.6597,0.7519,0.5525,0.1628,0.3184");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_mC24_15_2e_af_e"){
        vparameters.push_back("6.0875108932,1.66963707964,1.11173475337,123.707427092,0.1695,0.483,0.8419,0.0299,0.6684,0.5221");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_mC40_15_2e_3f_f"){
        vparameters.push_back("10.7421147996,0.897270659587,0.978465510971,147.31111283,0.45571,0.83777,0.78937,0.64875,0.16687,0.71802,0.50232,0.32358,0.43573,0.74259,0.03818,0.719,0.66021,0.26127");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_mC48_15_aef_3f_2e"){
        vparameters.push_back("5.8157729384,1.73589080063,2.00719678852,109.149281309,0.8352,0.16749,0.49971,0.9884,0.6705,0.0001,0.75396,0.48644,0.8628,0.21836,0.66556,0.86281,0.25504,0.34412,0.86546");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC7_mC48_15_2f_e_e3f"){
        vparameters.push_back("13.2432547747,0.690228819766,0.966427058766,155.232510527,0.9409,0.2269,0.8889,0.6633,0.303,0.3712,0.8094,0.241,0.2004,0.6988,0.5659,0.2806,0.9947,0.1491,0.6373,0.8064,0.5797");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_mC60_15_cf_e4f_ef"){
        vparameters.push_back("12.5843824384,0.579976673828,0.889919034503,68.755,0.6408,0.16,0.587,0.2597,0.0011,0.6275,0.0955,0.4035,0.729,0.1009,0.1137,0.6109,0.8071,0.2522,0.5497,0.5564,0.1061,0.6278,0.047,0.2555");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP16_33_a_a_2a"){
        vparameters.push_back("5.9547396971,1.2544327646,0.947018118432,0.93178,0.87342,-0.0,0.5759,0.8771,0.503,0.9576,0.9125,0.656,0.6123,0.839,0.9169");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_oC24_36_b_ab_a"){
        vparameters.push_back("10.2764819278,0.576935497153,0.499471075467,0.108,0.155,0.167,0.504,0.33,0.83,0.0,0.353,0.81,0.594");
        vparameters.push_back("11.8784,0.57718211,0.4611648,0.077,0.895,0.166,0.563,0.666,0.839,0,0.63,0.786,0.5");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C_oP32_58_eg_3gh_g"){
        vparameters.push_back("7.7975902487,1.01536533044,0.715127549044,0.2419,0.8702,0.3614,0.922,0.1368,0.574,0.6392,0.8969,0.6001,0.754,0.7495,0.7697,0.8676,0.2387");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_oC28_63_c_ac_fg"){
        vparameters.push_back("6.2543235982,1.56175972927,1.21827411167,0.3472,0.7,0.75,0.5694,0.7778,0.45");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_oC32_63_c_c2f_f"){
        vparameters.push_back("3.8117041359,2.59759679576,2.66755674233,0.1943,0.7775,0.9542,0.886,0.687,0.9362,0.8659,0.4343");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_tI24_88_a_b_f"){
        vparameters.push_back("5.7543,2.2969084,0.8831,0.5138,0.2973");
        vparameters.push_back("5.7062,2.20627,0.2638,0.1105,0.4568");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC2_tP28_91_2d_b_ac"){
        vparameters.push_back("6.0709632596,1.40149875103,0.233,0.246,0.246,0.975,0.73,0.25,0.513,0.264,0.231");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_hP56_173_2b2c_ac_b5c"){
        vparameters.push_back("11.1609261164,0.841147428044,0.25,0.95,0.56,0.74,0.163,0.341,0.056,0.152,0.325,0.45,0.507,0.009,0.253,0.181,0.0,0.978,0.69,0.005,0.058,0.498,0.181,0.004,0.175,0.498,0.993,0.1161,0.32,0.255");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_cF56_227_b_c_e"){
        vparameters.push_back("9.8805440428,0.2617");
      }
      // -------------------------------------------------------------------------
      // nitride prototypes (from R. Friedrich)
      // -------------------------------------------------------------------------
      // binaries
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_cI80_206_ad_e"){
        vparameters.push_back("9.8892799651,0.7716,0.1259,0.6002,0.6475");
      }
      // -------------------------------------------------------------------------
      // ternaries
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP12_62_c_c_c"){
        vparameters.push_back("9.95851159,0.433951127369,0.653641836852,0.15617,0.97349,0.04761,0.42551,0.40749,0.75982");
        vparameters.push_back("7.1586,0.48988629,0.70044422,0.0948,0.4719,0.1442,0.9797,0.3792,0.731");  // 002, ternary metal-nitride prototype (ICSD #93259)
        vparameters.push_back("7.86,0.59287532,1.1119593,0.482,0.317,0.352,0.939,0.732,0.886");  // 003, metallic ternary prototype (ICSD #42757)
        vparameters.push_back("7.146,0.62244612,1.0720683,0.8048,0.9093,0.6959,0.589,0.012,0.2975");  // 004, metallic ternary prototype (ICSD #54301)
        vparameters.push_back("7.199,0.62911515,1.0247257,0.859,0.435,0.017,0.7973,0.204,0.371");  // 005, metallic ternary prototype (ICSD #54415)
        vparameters.push_back("7.3628,0.62931222,1.0731651,0.0105,0.7963,0.7969,0.4187,0.1866,0.4123");  // 006, metallic ternary prototype (ICSD #54460)
        vparameters.push_back("6.897,0.6415833,1.1204872,0.8555,0.0637,0.7261,0.3829,0.9679,0.6799");  // 007, metallic ternary prototype (ICSD #54938)
        vparameters.push_back("7.77,0.61132561,1.0862291,0.731,0.901,0.494,0.296,0.387,0.917");  // 008, metallic ternary prototype (ICSD #58762)
        vparameters.push_back("6.625,0.65403774,1.1375094,0.198,0.431,0.462,0.715,0.7839,0.408");  // 009, metallic ternary prototype (ICSD #105338)
        vparameters.push_back("7.156,0.6352711,1.0444382,0.839,0.439,0.974,0.794,0.17,0.426");  // 010, metallic ternary prototype (ICSD #106459)
        vparameters.push_back("7.671,0.60643984,0.98970147,0.478,0.809,0.688,0.452,0.815,0.07");  // 011, metallic ternary prototype (ICSD #108571)
        vparameters.push_back("7.7672,0.61037439,0.94261767,0.9694,0.3161,0.176,0.9257,0.3232,0.6007");  // 012, metallic ternary prototype (ICSD #157921)
        vparameters.push_back("7.6,0.61628947,0.97888158,0.1794,0.924,0.9737,0.3099,0.3181,0.5972");  // 013, metallic ternary prototype (ICSD #157923)
        vparameters.push_back("6.723,0.62977837,1.1518667,0.012,0.304,0.677,0.551,0.794,0.903");  // 014, metallic ternary prototype (ICSD #630573)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_tI10_139_e_e_a"){
        vparameters.push_back("4.0191421247,3.53238454024,0.664,0.8545");
        vparameters.push_back("3.8568,3.3538166,0.1449,0.3409");  // 002, ternary metal-nitride prototype (ICSD #80376)
        vparameters.push_back("4.152,3.1442678,0.6557,0.8589");  // 003, ternary metal-nitride prototype (ICSD #80377)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC2_hP6_191_f_a_d"){
        vparameters.push_back("5.3369,0.57059342");
      }

      // -------------------------------------------------------------------------
      // metal-nitride prototypes (from DX)
      // ---------------------------------------------------------------------------
      // binaries 
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP6_6_3ab_ab"){
        vparameters.push_back("4.86,0.58148148,0.84053498,99.49,0.192,0.064,0.825,0.654,0.62,0.359,0.201,0.425,0.353,0.805,0.778,0.028");  // 001, binary metal-nitride prototype (ICSD #187448)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC12_8_4a_2a"){
        vparameters.push_back("4.964,0.57836422,1.9284851,107.49,0.252,0.798,0.162,0.643,0.008,0.376,0.407,0.065,0.999,0.936,0.749,0.505");  // 001, binary metal-nitride prototype (ICSD #187447)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_mC16_8_6a_b"){
        vparameters.push_back("9.0123,0.83308367,0.70781044,134.9427,0.86936937,0.5,0.12962963,0.5,0.63036963,0.76073926,0.36862937,0.23926074,0,0.4995005,0.5004995,0.5004995,0.5,0.75,0");  // 001, binary metal-nitride prototype (ICSD #155169)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP6_11_2e_e"){
        vparameters.push_back("5.1164,0.56651161,0.82880541,113.029,0.29511678,0.65490933,0.17069835,0.88921147,0.71680688,0.70665402");  // 001, binary metal-nitride prototype (ICSD #187444)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6_mP14_11_e_6e"){
        vparameters.push_back("9.63,0.45794393,0.56386293,99.34,0.7819,0.8285,0.093,0.86,0.108,0.646,0.112,0.437,0.416,0.221,0.54,0.25,0.667,0.298");  // 001, binary metal-nitride prototype (ICSD #14244)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_mC8_12_di_a"){
        vparameters.push_back("6.2233,0.58666624,0.85740684,71.556,0.0939,0.2687");  // 001, binary metal-nitride prototype (ICSD #29370)
        vparameters.push_back("6.3105,0.57963711,0.86886934,68.788,0.0836,0.2705");  // 002, binary metal-nitride prototype (ICSD #29376)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC12_12_2i_i"){
        vparameters.push_back("6.817,0.41587208,0.84264339,83.9924,0.373,0.864,0.902,0.47,0.675,0.769");  // 001, binary metal-nitride prototype (ICSD #187441)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC16_12_2i_2i"){
        vparameters.push_back("13.46,0.2830312,0.49946508,85.287,0.7559,0.7486,0.9779,0.0833,0.846,0.3926,0.6441,0.092");  // 001, binary metal-nitride prototype (ICSD #411555)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_mC30_12_a4i_3i"){
        vparameters.push_back("15.0498,0.25136547,0.62120427,79.46,0.62567,0.72225,0.328,0.96482,0.28973,0.62051,0.96826,0.66055,0.14105,0.51172,0.18256,0.86977,0.46936,0.81225");  // 001, binary metal-nitride prototype (ICSD #162794)
        vparameters.push_back("14.095,0.25902802,0.63213906,81.748,0.328,0.965,0.299,0.629,0.631,0.715,0.974,0.658,0.141,0.508,0.18,0.868,0.465,0.807");  // 002, binary metal-nitride prototype (ICSD #169726)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP6_13_g_e"){
        vparameters.push_back("5.053,0.58321789,0.98495943,129.9799,0.809,0.679,0.698,0.37");  // 001, binary metal-nitride prototype (ICSD #187451)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mP12_14_e_2e"){
        vparameters.push_back("5.6053,0.86029294,0.85947229,126.0785,0.2336,0,0.0131,0.1841,0.4,0.8837,0.3261,0.5888,0.169");  // 001, binary metal-nitride prototype (ICSD #160623)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mP16_14_e_3e"){
        vparameters.push_back("5.9607,2.1122855,0.60051672,103.253,0.2956,0.03505,0.3855,0.949,0.0987,0.19,0.914,0.1852,0.297,0.865,0.2681,0.362");  // 001, binary metal-nitride prototype (ICSD #98661)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mC12_15_e_f"){
        vparameters.push_back("7.1712,0.61281236,1.2249554,127.232,0.5516,0.491,0.607,0.032");  // 001, binary metal-nitride prototype (ICSD #280681)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC16_36_2a_2a"){
        vparameters.push_back("2.59,4.6100386,1.7722008,0.86,0,0.77,0.82,0.47,0.82,0.15,0.93");  // 001, binary metal-nitride prototype (ICSD #167514)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_53_h_h"){
        vparameters.push_back("2.79,1.344086,1.5197133,0.35,0.7,0.86,0.89");  // 001, binary metal-nitride prototype (ICSD #167513)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_57_e_c"){
        vparameters.push_back("2.901,1.7087211,2.7083764,0.297,0.792,0.08,0.838");  // 001, binary metal-nitride prototype (ICSD #187450)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_oC32_63_c2f_cf"){
        vparameters.push_back("3.893,2.6365271,2.6365271,0.2373,0.70243,0.047,0.8812,0.3097,0.9268,0.13429,0.4399");  // 001, binary metal-nitride prototype (ICSD #16253)
        vparameters.push_back("3.8862,2.6277083,2.6407287,0.23678,0.8029,0.04701,0.88051,0.30862,0.92622,0.13455,0.44094");  // 002, binary metal-nitride prototype (ICSD #66533)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC12_65_gj_e"){
        vparameters.push_back("5.701,1.3513419,0.49587792,0.117,0.31");  // 001, binary metal-nitride prototype (ICSD #187452)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oI4_71_b_c"){
        vparameters.push_back("2.8796,1.479268,1.0460133");  // 001, binary metal-nitride prototype (ICSD #53146)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oI12_71_2e_f"){
        vparameters.push_back("13.0613,0.28741396,0.215415,0.45057355,0.20799023,0.84176015");  // 001, binary metal-nitride prototype (ICSD #187449)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oI16_71_abe_n"){
        vparameters.push_back("10.7912,0.41118689,0.2889484,0.751,0.12321,0.6466");  // 001, binary metal-nitride prototype (ICSD #423831)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oI16_72_b_cj"){
        vparameters.push_back("5.94,0.94107744,1.0185185,0.855,0.145");  // 001, binary metal-nitride prototype (ICSD #27135)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B4_tI18_87_ah_h"){
        vparameters.push_back("6.873,0.62534556,0.6,0.2,0.9,0.3");  // 001, binary metal-nitride prototype (ICSD #26251)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI32_88_d_cf"){
        vparameters.push_back("8.653,0.64648099,0.173,0.827,0.625");  // 001, binary metal-nitride prototype (ICSD #30633)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP2_123_d_a"){
        vparameters.push_back("4.129,0.7808186");  // 001, binary metal-nitride prototype (ICSD #168645)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP4_123_ag_d"){
        vparameters.push_back("4.497,0.82432733,0.304");  // 001, binary metal-nitride prototype (ICSD #16963)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP8_123_egh_ab"){
        vparameters.push_back("3.745,2,0.25,0.25");  // 001, binary metal-nitride prototype (ICSD #180237)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP6_127_g_b"){
        vparameters.push_back("4.408,0.60528584,0.61441863");  // 001, binary metal-nitride prototype (ICSD #187442)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP4_131_e_c"){
        vparameters.push_back("2.704,2.2525888");  // 001, binary metal-nitride prototype (ICSD #187709)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI4_139_a_b"){
        vparameters.push_back("3.009,1.3921569");  // 001, binary metal-nitride prototype (ICSD #106932)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_tI10_139_ae_e"){
        vparameters.push_back("4.2,2.8878571,0.667,0.8412");  // 001, binary metal-nitride prototype (ICSD #71638)
        vparameters.push_back("2.974,4.0773369,0.33341,0.1593");  // 002, binary metal-nitride prototype (ICSD #84202)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_tI14_139_ad_ce"){
        vparameters.push_back("4.382,1.9698768,0.2521");  // 001, binary metal-nitride prototype (ICSD #76389)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B_tI18_139_deh_a"){
        vparameters.push_back("5.72,1.0996503,0.293,0.757");  // 001, binary metal-nitride prototype (ICSD #41953)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI16_140_a_dh"){
        vparameters.push_back("6.106,1.1578775,0.635");  // 001, binary metal-nitride prototype (ICSD #24007)
        vparameters.push_back("6.5412,1.2368984,0.374");  // 002, binary metal-nitride prototype (ICSD #25008)
        vparameters.push_back("5.52,1.009058,0.84");  // 003, binary metal-nitride prototype (ICSD #183201)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI16_140_a_ch"){
        vparameters.push_back("6.1624,1.0959853,0.4359");  // 001, binary metal-nitride prototype (ICSD #187018)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI16_140_c_ah"){
        vparameters.push_back("6.1188,1.1606034,0.36406");  // 001, binary metal-nitride prototype (ICSD #1145)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI12_141_a_e"){
        vparameters.push_back("4.14,2.1268116,0.64");  // 001, binary metal-nitride prototype (ICSD #23403)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hR4_160_3a_a"){
        vparameters.push_back("3.6458651,4.1802255,0.428,0.5,0.583,0");  // 001, binary metal-nitride prototype (ICSD #644523)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hR5_160_3a_2a"){
        vparameters.push_back("2.8499718,5.1435171,0.14085,0.68246,0.91112,0.79584,0.29163");  // 001, binary metal-nitride prototype (ICSD #185490)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP16_162_ek_ci"){
        vparameters.push_back("5.665,0.97440424,0.25,0.83333,0.494,0.75");  // 001, binary metal-nitride prototype (ICSD #43559)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP5_164_ac_d"){
        vparameters.push_back("2.82,2.2163121,0.393,0.802");  // 001, binary metal-nitride prototype (ICSD #182700)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP16_164_ci_di"){
        vparameters.push_back("5.72,0.97902098,0.25,0,0.511,0.25,0.83333,0");  // 001, binary metal-nitride prototype (ICSD #60168)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP16_186_ac_bc"){
        vparameters.push_back("5.745,0.97859008,0.25,0,0.51,0.25,0.833,0");  // 001, binary metal-nitride prototype (ICSD #76280)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP3_187_a_h"){
        vparameters.push_back("3.2689,1.1635107,0.687679");  // 001, binary metal-nitride prototype (ICSD #290427)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP3_191_e_a"){
        vparameters.push_back("2.829,1.7575115,0.384");  // 001, binary metal-nitride prototype (ICSD #260545)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP6_191_af_d"){
        vparameters.push_back("5.263,0.5166255");  // 001, binary metal-nitride prototype (ICSD #182351)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_191_ad_f"){
        vparameters.push_back("5.186,0.56170459");  // 001, binary metal-nitride prototype (ICSD #76457)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_191_bc_ad"){
        vparameters.push_back("3.98,1.2084925");  // 001, binary metal-nitride prototype (ICSD #161079)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_193_g_b"){
        vparameters.push_back("7.6418,0.9225706,0.2721");  // 001, binary metal-nitride prototype (ICSD #77730)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_194_a_c"){
        vparameters.push_back("2.46,1.5859756");  // 001, binary metal-nitride prototype (ICSD #183157)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP6_194_c_e"){
        vparameters.push_back("3.2647,2.3412565,0.5930479");  // 001, binary metal-nitride prototype (ICSD #290428)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP10_194_bf_ac"){
        vparameters.push_back("2.8413,3.4114666,0.075");  // 001, binary metal-nitride prototype (ICSD #25656)
        vparameters.push_back("3.8911,3.3888361,0.5862");  // 002, binary metal-nitride prototype (ICSD #162797)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP10_194_cf_f"){
        vparameters.push_back("2.89,5.2892734,0.07,0.84");  // 001, binary metal-nitride prototype (ICSD #186207)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP16_194_bh_ag"){
        vparameters.push_back("5.72,0.96153846,0.52");  // 001, binary metal-nitride prototype (ICSD #106926)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_cI32_204_c_g"){
        vparameters.push_back("6.592,0.6933,0.0944");  // 001, binary metal-nitride prototype (ICSD #162105)
        vparameters.push_back("6.289,0.8855,0.3784");  // 002, binary metal-nitride prototype (ICSD #162106)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_cP5_221_a_bc"){
        vparameters.push_back("3.74");  // binary metal-nitride prototype (ICSD #76403)
      }
      // -------------------------------------------------------------------------
      // ternaries
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C5_aP8_1_a_2a_5a"){
        vparameters.push_back("6.57,0.6173516,0.56754947,90.52,90.53,91.16,0.788,0,0.759,0.791,0.497,0.745,0.291,0.505,0.238,0.958,0.555,0.254,0.124,0.447,0.732,0.285,0.998,0.127,0.621,0.465,0.246,0.459,0.536,0.74");  // 001, ternary metal-nitride prototype (ICSD #92316)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C5_mC16_5_a_c_a2bc"){
        vparameters.push_back("6.731,0.88307829,0.94829892,88.82,0.093,0.631,0.084,0.617,0.73,0.85,0.763,0.715,0.848,0.232");  // 001, ternary metal-nitride prototype (ICSD #92315)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C3_mP8_6_2ab_ab_a2b"){
        vparameters.push_back("6.32,0.58164557,0.60981013,90.31,0.972,0.462,0.283,0,0.297,0.5,0.638,0.566,0.819,0,0.814,0.502,0.473,0.395,0.143,0.617");  // 001, ternary metal-nitride prototype (ICSD #92312)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C_mP8_6_3a2b_ab_b"){
        vparameters.push_back("6.6606,0.57952737,0.64258475,91.4279,0.363,0.521,0.67,0.005,0.026,0.447,0.688,0.505,0.221,0,0.535,0.61,0.207,0.498,0.872,0.659");  // 001, ternary metal-nitride prototype (ICSD #92314)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_mC20_12_2i_i_2i"){
        vparameters.push_back("10.93,0.45379689,0.86336688,141.6047,0.726,0.7657,0.0848,0.4066,0.5029,0.8751,0.7185,0.1958,0.4256,0.613");  // 001, ternary metal-nitride prototype (ICSD #72389)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C12_mC30_12_a_i_6i"){
        vparameters.push_back("14.272,0.26534473,0.62268778,87.17,0.347,0.6745,0.3979,0.0033,0.3277,0.0835,0.2593,0.1569,0,0.7378,0.068,0.6555,0.133,0.5702");  // 001, ternary metal-nitride prototype (ICSD #31297)
        vparameters.push_back("14.254,0.2655395,0.62319349,87.1,0.847,0.674,0.398,0.003,0.328,0.084,0.259,0.157,0,0.738,0.068,0.656,0.133,0.57");  // 002, ternary metal-nitride prototype (ICSD #659621)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_mP12_14_e_e_e"){
        vparameters.push_back("4.54,1.0484581,1.2790749,124.9,0.91998,0.39152,0.608,0.67587,0.54156,0.91072,0.84932,0.19521,0.8007");  // 001, ternary metal-nitride prototype (ICSD #402341)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_63_2c_c_c"){
        vparameters.push_back("3.532,5.7202718,1.4068516,0.78124,0.62006,0.96137,0.19906");  // 001, ternary metal-nitride prototype (ICSD #96228)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oC20_63_c_a_cf"){
        vparameters.push_back("3.369,3.4128822,2.6662214,0.2526,0.9554,0.3719,0.958");  // 001, ternary metal-nitride prototype (ICSD #29521)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC3_oC28_63_cg_c_cg"){
        vparameters.push_back("8.503,1.2094555,0.59179113,0.8921,0.3055,0.1243,0.2843,0.1174,0.6918,0.8714");  // 001, ternary metal-nitride prototype (ICSD #40205)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC8_65_i_c_b"){
        vparameters.push_back("3.698,1.7212006,1.2439156,0.166");  // 001, ternary metal-nitride prototype (ICSD #92307)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C3_oC24_65_a2ij_i_cj"){
        vparameters.push_back("4.9,5.0531224,0.99591837,0.1982,0.40004,0.1027,0.0998,0.29317");  // 001, ternary metal-nitride prototype (ICSD #98178)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_oI12_71_e_a_de"){
        vparameters.push_back("12.4601,0.30095264,0.2764424,0.35332,0.1638");  // 001, ternary metal-nitride prototype (ICSD #50579)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_99_ab_ab_ab"){
        vparameters.push_back("3.0285,2.430246,0.76,0.05,0.35,0.58,0.3,0");  // 001, ternary metal-nitride prototype (ICSD #23779)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C3_tP8_115_abc_g_de"){
        vparameters.push_back("3.8929,2.021886,0.23552817,0.098179023");  // 001, ternary metal-nitride prototype (ICSD #92311)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP3_123_c_a_b"){
        vparameters.push_back("2.83,1.3120141");  // 001, ternary metal-nitride prototype (ICSD #53505)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C_tP8_123_ai_bc_d"){
        vparameters.push_back("3.965,1.3881463,0.75");  // 001, ternary metal-nitride prototype (ICSD #92313)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_129_c_c_c"){
        vparameters.push_back("3.5701,2.1170836,0.65036,0.08393,0.3303");  // 001, ternary metal-nitride prototype (ICSD #2027)
        vparameters.push_back("3.081,2.5316456,0.9,0.617,0.335");  // 002, ternary metal-nitride prototype (ICSD #100437)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP8_129_c_c_bc"){
        vparameters.push_back("4.1279,2.0304755,0.8479,0.4142,0.168");  // 001, ternary metal-nitride prototype (ICSD #50994)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP8_129_a_c_bc"){
        vparameters.push_back("3.895,1.5697047,0.221,0.789");  // 001, ternary metal-nitride prototype (ICSD #92309)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_131_e_c_b"){
        vparameters.push_back("3.5809,1.9574967");  // 001, ternary metal-nitride prototype (ICSD #69044)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_131_b_c_e"){
        vparameters.push_back("3.924,1.8055556");  // 001, ternary metal-nitride prototype (ICSD #87414)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tI8_139_d_b_a"){
        vparameters.push_back("3.653,1.460991");  // 001, ternary metal-nitride prototype (ICSD #92305)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B6C2_tI22_139_ae_eg_e"){
        vparameters.push_back("3.905,5.1766965,0.1805,0.6949,0.6023,0.5871");  // 001, ternary metal-nitride prototype (ICSD #98477)
        vparameters.push_back("4.0698,4.9520861,0.82019,0.3019,0.39895,0.9153");  // 002, ternary metal-nitride prototype (ICSD #411473)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_tI20_140_h_h_a"){
        vparameters.push_back("5.6213,1.3433192,0.872,0.333");  // 001, ternary metal-nitride prototype (ICSD #413356)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tI20_140_h_a_h"){
        vparameters.push_back("5.5615,1.2370224,0.3794,0.8312");  // 001, ternary metal-nitride prototype (ICSD #413357)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tI20_140_a_h_h"){
        vparameters.push_back("5.6646,1.4818699,0.372,0.836");  // 001, ternary metal-nitride prototype (ICSD #415304)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hR4_160_a_2a_a"){
        vparameters.push_back("2.8560976,5.4641029,0.8251,0.2629,0.4068,0");  // 001, ternary metal-nitride prototype (ICSD #84639)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_hP5_164_d_d_a"){
        vparameters.push_back("3.282,1.6636197,0.39,0.77");  // 001, ternary metal-nitride prototype (ICSD #16231)
        vparameters.push_back("3.622,1.7556599,0.6314,0.2797");  // 002, ternary metal-nitride prototype (ICSD #410826)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_hP5_164_a_d_d"){
        vparameters.push_back("3.557,1.5451223,0.65,0.27");  // 001, ternary metal-nitride prototype (ICSD #34003)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hR4_166_a_c_b"){
        vparameters.push_back("3.1359829,5.5605955,0.881");  // 001, ternary metal-nitride prototype (ICSD #71136)
        vparameters.push_back("4.3180071,3.9173242,0.189");  // 002, metallic ternary prototype (ICSD #615828)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC3_hP14_176_h_c_h"){
        vparameters.push_back("8.0141,0.69977914,0.35834333,0.086966667,0.12945233,0.68896467");  // 001, ternary metal-nitride prototype (ICSD #36502)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_186_b_b_a"){
        vparameters.push_back("3.3306,3.2702516,0,0.598,0.2466");  // 001, ternary metal-nitride prototype (ICSD #172471)
        vparameters.push_back("4.3972,1.3683253,0,0.777,0.245");  // 002, metallic ternary prototype (ICSD #54657)
        vparameters.push_back("4.5825,1.7151555,0,0.2398,0.7854");  // 003, metallic ternary prototype (ICSD #416143)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_186_a_ab_b"){
        vparameters.push_back("2.9213,3.7507274,0,0.2,0.4,0.8");  // 001, ternary metal-nitride prototype (ICSD #80029)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP3_187_a_c_d"){
        vparameters.push_back("3.64,0.93956044");  // 001, ternary metal-nitride prototype (ICSD #247028)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3C3_hP11_189_dg_g_f"){
        vparameters.push_back("6.475,0.54903475,0.3557,0.7281,0.3537");  // 001, ternary metal-nitride prototype (ICSD #411152)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_hP4_191_c_a_b"){
        vparameters.push_back("3.65,1.260274");  // 001, ternary metal-nitride prototype (ICSD #92308)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP4_191_b_a_c"){
        vparameters.push_back("4,1.05");  // 001, ternary metal-nitride prototype (ICSD #92310)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP8_194_c_a_f"){
        vparameters.push_back("3.6506,3.4680053,0.595");  // 001, ternary metal-nitride prototype (ICSD #74791)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_194_a_f_b"){
        vparameters.push_back("2.8724,3.8201504,0.128");  // 001, ternary metal-nitride prototype (ICSD #75971)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP8_194_a_c_f"){
        vparameters.push_back("2.9108,3.6245362,0.6384");  // 001, ternary metal-nitride prototype (ICSD #185913)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_hP10_194_h_a_c"){
        vparameters.push_back("8.4414,0.82707845,0.8561");  // 001, ternary metal-nitride prototype (ICSD #67497)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_hP10_194_h_c_a"){
        vparameters.push_back("7.6111,0.87756041,0.83953333");  // 001, ternary metal-nitride prototype (ICSD #152055)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_hP16_194_c_ae_2f"){
        vparameters.push_back("2.795,8.5144902,0.6078,0.0543,0.6591");  // 001, ternary metal-nitride prototype (ICSD #181351)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_cF32_227_b_c_a"){
        vparameters.push_back("8.74");  // ternary metal-nitride prototype (ICSD #72546)
      }
      
      // -------------------------------------------------------------------------
      // metal prototypes (from DX) //DX20201028
      // ---------------------------------------------------------------------------
      // ternaries
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC16_8_b_2a_2ab"){
        vparameters.push_back("7.171,1.2100126,0.7774369,71.74,0.12,0,0.641,0.978,0.961,0.538,0.551,0.528,0.281,0.224,0.42,0.417,0.756,0.876");  // 001, metallic ternary prototype (ICSD #58830)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C3_mP9_10_m_2n_am"){
        vparameters.push_back("8.225,0.45483283,0.91671733,109.65,0.89515,0.35185,0.4983,0.7077,0.1823,0.39514,0.2748,0.90637");  // 001, metallic ternary prototype (ICSD #417035)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C4_mP11_10_m_a2m_2n"){
        vparameters.push_back("9.046,0.36325448,0.82367897,94.53,0.6346,0.5766,0.8578,0.343,0.6676,0.9576,0.1541,0.3154,0.5619,0.295");  // 001, metallic ternary prototype (ICSD #55578)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_mP10_11_e_2e_2e"){
        vparameters.push_back("10.783,0.45015302,0.45618103,91.3,0.7587,0.2575,0.4999,0.7506,0.1184,0.2145,0,0.6962,0.3728,0.2462");  // 001, metallic ternary prototype (ICSD #416299)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_mP10_11_2e_e_2e"){
        vparameters.push_back("11.3765,0.40427196,0.40718147,90.542,0.13236,0.2261,0.49959,0.7565,0.75464,0.2577,0.3606,0.2566,0.001,0.7307");  // 001, metallic ternary prototype (ICSD #424108)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_mP14_11_2e_4e_e"){
        vparameters.push_back("9.591,0.47361068,0.79799812,107.838,0.05956,0.27432,0.57956,0.27679,0.8773,0.90331,0.57123,0.89518,0.36827,0.46211,0.87129,0.49589,0.21178,0.8337");  // 001, metallic ternary prototype (ICSD #261042)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_mC10_12_i_a_i"){
        vparameters.push_back("11.1,0.38801802,0.95268468,154.4287,0.5942,0.8555,0.9238,0.5702");  // 001, metallic ternary prototype (ICSD #182050)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_mC12_12_i_i_i"){
        vparameters.push_back("11.098,0.3973689,0.98442963,154.8009,0.2238,0.8102,0.8046,0.675,0.3393,0.1465");  // 001, metallic ternary prototype (ICSD #20632)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C2_mC14_12_ai_i_i"){
        vparameters.push_back("11.019,0.38767583,0.94896089,136.7589,0.30129,0.44343,0.42164,0.79411,0.8673,0.22347");  // 001, metallic ternary prototype (ICSD #410985)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C_mC18_12_3i_i_a"){
        vparameters.push_back("10.077,0.54314776,0.92448149,81.93,0.5297,0.2753,0.7074,0.9295,0.8264,0.5908,0.85526,0.24903");  // 001, metallic ternary prototype (ICSD #260159)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C4_mC22_12_i_a2i_2i"){
        vparameters.push_back("16.9369,0.27129522,0.97586926,153.6114,0.38824,0.6115,0.04121,0.23355,0.40758,0.34101,0.71885,0.06875,0.2201,0.643");  // 001, metallic ternary prototype (ICSD #419134)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C8_mC30_12_ai_hi_2ij"){
        vparameters.push_back("18.549,0.36039679,0.36616529,81.66,0.213,0.67755,0.6913,0.791,0.13,0.8708,0.4541,0.51517,0.202,0.14945,0.26752,0.19789");  // 001, metallic ternary prototype (ICSD #240016)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C7_mC30_12_3i_i_a3i"){
        vparameters.push_back("14.257,0.32012345,0.86806481,86.021,0.63732,0.49687,0.46213,0.20062,0.19481,0.20687,0.83052,0.11118,0.81799,0.32444,0.01039,0.38194,0.65953,0.01844");  // 001, metallic ternary prototype (ICSD #171243)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_mC20_15_e_f_f"){
        vparameters.push_back("10.487,0.5754744,0.79202823,76.32,0.1345,0.63313,0.36433,0.99699,0.85144,0.3604,0.64786");  // 001, metallic ternary prototype (ICSD #391432)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_mC28_15_f_2f_e"){
        vparameters.push_back("19.0003,0.35246812,0.35005237,86.065,0.7786,0.2981,0.0987,0.1999,0.42255,0.20721,0.77201,0.34791,0.51916,0.01422");  // 001, metallic ternary prototype (ICSD #261117)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_20_b_a_c"){
        vparameters.push_back("8.38,1.2601432,0.75536993,0.7698,0.504,0.351,0.673,0.513");  // 001, metallic ternary prototype (ICSD #1156)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP12_33_a_a_a"){
        vparameters.push_back("7.5356,1.0569563,0.62430331,0.51211,0.79877,0.75,0.29931,0.41498,0.7559,0.18658,0.0868,0.7533");  // 001, metallic ternary prototype (ICSD #106418)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC12_36_a_a_a"){
        vparameters.push_back("4.8891,1.720603,1.6021354,0.66671918,0.026102028,0.99999382,0.23452562,0.66654623,0.45737235");  // 001, metallic ternary prototype (ICSD #55819)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC32_36_b_2a_2ab"){
        vparameters.push_back("8.681,0.82363783,1.2135699,0.638,0.011,0.12,0,0.082,0.731,0.673,0.736,0.7763,0.622,0.29,0.756,0.619,0.563");  // 001, metallic ternary prototype (ICSD #56278)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC18_38_ad_be_bd"){
        vparameters.push_back("4.0307,3.2806709,1.8867442,0.74976771,0.40946388,0.0012150248,0.8747892,0.1251434,0.66729312,0.00026107832,0.70531985,0.29437221");  // 001, metallic ternary prototype (ICSD #59430)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC18_38_bd_be_ad"){
        vparameters.push_back("4.1297,3.1724823,1.8520958,0.2464,0,0.5914,0.83335,0.50005,0.3768,0.3768,0.2043,0.2043");  // 001, metallic ternary prototype (ICSD #260992)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC18_38_be_ad_bd"){
        vparameters.push_back("3.9497,3.3839026,1.9576677,0.25494797,0.59120572,1.0483556e-05,0.12745096,0.87256621,0.33335379,2.1695924e-05,0.79560917,0.20433001");  // 001, metallic ternary prototype (ICSD #415506)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10BC3_oP14_47_2q2rs_h_bt"){
        vparameters.push_back("4.261,1.0073926,3.6026754,0.3187,0.1603,0.0819,0.4156,0.4159,0.255");  // 001, metallic ternary prototype (ICSD #20663)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP8_51_f_e_be"){
        vparameters.push_back("5.254,0.78625809,1.3644842,0.2824,0.8886,0.656");  // 001, metallic ternary prototype (ICSD #106800)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oP8_51_ae_f_f"){
        vparameters.push_back("5.09,0.79980354,1.3709234,0.315,0.118,0.683");  // 001, metallic ternary prototype (ICSD #623100)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oP8_51_af_f_e"){
        vparameters.push_back("4.998,0.80712285,1.412565,0.346,0.103,0.715");  // 001, metallic ternary prototype (ICSD #658675)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_oP12_51_afj_e_e"){
        vparameters.push_back("7.7097,0.52203069,0.91056721,0.17610349,0.61101001,0.95660968,0.92790085,0.71070497");  // 001, metallic ternary prototype (ICSD #9986)
        vparameters.push_back("7.6602,0.5295162,0.90252213,0.61759,0.18412,0.95685,0.93045,0.70183");  // 002, metallic ternary prototype (ICSD #55598)
        vparameters.push_back("8.637,0.48917448,0.86036818,0.19372,0.60198,0.93439,0.93656,0.68597");  // 003, metallic ternary prototype (ICSD #411393)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C_oP12_51_ak_ef_f"){
        vparameters.push_back("5.66,1.1909894,1.5980565,0.3741,0.1328,0.57444,0.20949,0.78478");  // 001, metallic ternary prototype (ICSD #249924)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oP12_51_e_afj_e"){
        vparameters.push_back("8.676,0.48744813,0.85892116,0.5998,0.19219,0.93475,0.93556,0.68673");  // 001, metallic ternary prototype (ICSD #410891)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP12_57_c_d_d"){
        vparameters.push_back("5.75,1.3547826,0.9773913,0.104,0.6363,0.3994,0.8726,0.0244");  // 001, metallic ternary prototype (ICSD #370036)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP16_57_d_2d_c"){
        vparameters.push_back("6.804,1.96796,0.95517343,0.093,0.6601,0.8945,0.841,0.4885,0.611,0.1664");  // 001, metallic ternary prototype (ICSD #107616)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP8_59_a_e_b"){
        vparameters.push_back("4.216,1.2018501,1.0661765,0.333,0.333,0.5,0.833");  // 001, metallic ternary prototype (ICSD #628580)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C_oP14_59_ae_ae_b"){
        vparameters.push_back("4.5526,1.7043008,1.9876554,0.04449,0.36673,0.2694,0.52383,0.58763,0.55436,0.88559");  // 001, metallic ternary prototype (ICSD #245679)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oP16_59_b_ae_ef"){
        vparameters.push_back("6.938,1.8080138,0.94018449,0.4367,0.68751,0.9941,0.1706,0.62693,0.33466,0.52489,0.95941");  // 001, metallic ternary prototype (ICSD #107444)
        vparameters.push_back("7.333,1.8847675,0.93113323,0.4473,0.6873,0.99594,0.1735,0.63365,0.33473,0.52812,0.962");  // 002, metallic ternary prototype (ICSD #107448)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP16_62_c_d_c"){
        vparameters.push_back("7.946,1.0187516,0.69808709,0.0332,0.89792,0.37781,0.11849,0.82793,0.9443,0.59466");  // 001, metallic ternary prototype (ICSD #425509)
        vparameters.push_back("7.38,0.95799458,0.75880759,0.379,0.652,0.033,0.35,0.823,0.947,0.088");  // 002, metallic ternary prototype (ICSD #602790)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP16_62_c_c_d"){
        vparameters.push_back("7.292,0.96941854,0.76851344,0.033,0.35,0.379,0.652,0.823,0.947,0.088");  // 001, metallic ternary prototype (ICSD #656966)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC12_63_c_c_c"){
        vparameters.push_back("3.6684,4.8307164,1.086441,0.7978,0.95664,0.58754");  // 001, metallic ternary prototype (ICSD #166874)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC12_63_c_c_a"){
        vparameters.push_back("4.2145,2.2335271,1.2607189,0.5972146,0.29158002");  // 001, metallic ternary prototype (ICSD #424600)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_63_f_c_c"){
        vparameters.push_back("4.08,2.4877451,1.7303922,0.213,0.939,0.348,0.956");  // 001, metallic ternary prototype (ICSD #57645)
        vparameters.push_back("4,2.3075,1.785,0.222,0.928,0.356,0.944");  // 002, metallic ternary prototype (ICSD #57693)
        vparameters.push_back("4.077,2.4724062,1.6997792,0.061,0.787,0.348,0.044");  // 003, metallic ternary prototype (ICSD #57759)
        vparameters.push_back("4.037,2.3661135,1.6019321,0.796,0.085,0.364,0.05");  // 004, metallic ternary prototype (ICSD #103875)
        vparameters.push_back("4.1824,2.448355,1.9307096,0.21983,0.93483,0.35459,0.9494");  // 005, metallic ternary prototype (ICSD #152048)
        vparameters.push_back("4.3642,2.589478,1.8801613,0.77214,0.05365,0.66459,0.94045");  // 006, metallic ternary prototype (ICSD #425493)
        vparameters.push_back("4.223,2.6381719,1.6663509,0.561,0.287,0.152,0.956");  // 007, metallic ternary prototype (ICSD #658140)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_63_a_c_g"){
        vparameters.push_back("7.7183,0.90569167,0.70617882,0.6485,0.819,0.3032");  // 001, metallic ternary prototype (ICSD #261785)
        vparameters.push_back("9.1362,0.8928548,0.62475646,0.3044,0.8289,0.7637");  // 002, metallic ternary prototype (ICSD #261786)
        vparameters.push_back("9.2835,0.8182582,0.63188453,0.3153,0.8185,0.6763");  // 003, metallic ternary prototype (ICSD #261788)
        vparameters.push_back("10.632,0.74473288,0.62046652,0.72594,0.8269,0.328");  // 004, metallic ternary prototype (ICSD #380341)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_63_c_c_2c"){
        vparameters.push_back("4.75,4.1263158,0.97894737,0.609,0.824,0.25,0.957");  // 001, metallic ternary prototype (ICSD #58647)
        vparameters.push_back("4.3923,3.6935091,0.98661294,0.6994,0.038,0.5648,0.2504");  // 002, metallic ternary prototype (ICSD #240098)
        vparameters.push_back("4.499,3.9319849,1,0.894,0.685,0.55,0.25");  // 003, metallic ternary prototype (ICSD #621687)
        vparameters.push_back("4.46,3.9506726,0.99080717,0.327,0.109,0.448,0.749");  // 004, metallic ternary prototype (ICSD #628742)
        vparameters.push_back("4.675,3.7411765,0.98524064,0.107,0.315,0.452,0.749");  // 005, metallic ternary prototype (ICSD #641754)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oC16_63_c_f_c"){
        vparameters.push_back("4.2266,2.4376567,1.9779492,0.55971,0.27339,0.1543,0.9448");  // 001, metallic ternary prototype (ICSD #96152)
        vparameters.push_back("4.621,2.3142177,1.6132872,0.9244,0.2069,0.6388,0.0506");  // 002, metallic ternary prototype (ICSD #150132)
        vparameters.push_back("4.433,2.5971126,1.8181818,0.94579,0.21864,0.66381,0.05497");  // 003, metallic ternary prototype (ICSD #414434)
        vparameters.push_back("4.397,2.2967933,1.5010234,0.589,0.3,0.13,0.94");  // 004, metallic ternary prototype (ICSD #600150)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_63_g_c_c"){
        vparameters.push_back("7.097,1.5019022,0.6299845,0.6747,0.13117,0.7127,0.90293");  // 001, metallic ternary prototype (ICSD #99139)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_63_c_c_f"){
        vparameters.push_back("4.4744,2.4125693,1.8175174,0.43302,0.71694,0.14403,0.44878");  // 001, metallic ternary prototype (ICSD #172230)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oC20_63_g_2c_a"){
        vparameters.push_back("7.0102,2.1704231,0.98316738,0.584,0.205,0.816,0.888");  // 001, metallic ternary prototype (ICSD #658701)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_oC24_63_acf_c_c"){
        vparameters.push_back("4.056,3.7302761,1.6346154,0.557,0.379,0.729,0.814,0.946");  // 001, metallic ternary prototype (ICSD #57760)
        vparameters.push_back("4.075,3.7398773,1.6220859,0.557,0.729,0.379,0.814,0.946");  // 002, metallic ternary prototype (ICSD #58060)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C2_oC24_63_c_cg_e"){
        vparameters.push_back("9.765,0.80747568,0.52892985,0.329,0.948,0.202,0.194,0.708");  // 001, metallic ternary prototype (ICSD #58898)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC3_oC28_63_cf_a_cf"){
        vparameters.push_back("4.098,2.4670571,3.1747194,0.3088,0.9983,0.5714,0.9114,0.289,0.8959");  // 001, metallic ternary prototype (ICSD #10044)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC32_63_ac_f_2cf"){
        vparameters.push_back("4.639,3.6861393,2.0773874,0.7976,0.5921,0.4277,0.3515,0.4931,0.1961,0.4111");  // 001, metallic ternary prototype (ICSD #410732)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oC10_65_j_i_a"){
        vparameters.push_back("4.0172,3.4774221,0.9045106,0.1986,0.3599");  // 001, metallic ternary prototype (ICSD #54612)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_oC12_65_i_cf_a"){
        vparameters.push_back("5.3308,1.7232498,0.69906956,0.66645548");  // 001, metallic ternary prototype (ICSD #107864)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C2_oC14_65_aj_j_i"){
        vparameters.push_back("3.9788,5.3289937,0.91944807,0.17641,0.4198,0.28698");  // 001, metallic ternary prototype (ICSD #240761)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_65_aci_j_i"){
        vparameters.push_back("4.192,4.1898855,0.98592557,0.288,0.416,0.136");  // 001, metallic ternary prototype (ICSD #103850)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC6_oC18_65_j_a_2ij"){
        vparameters.push_back("4.313,5.2124275,1.0159286,0.2046,0.4311,0.3222,0.0697");  // 001, metallic ternary prototype (ICSD #165783)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C7_oC24_65_aj_i_c2ij"){
        vparameters.push_back("4.565,5.9824754,1,0.875,0.593,0.7813,0.6847,0.9107");  // 001, metallic ternary prototype (ICSD #102239)
        vparameters.push_back("4.579,5.7896921,1,0.8718,0.592,0.7812,0.6845,0.9087");  // 002, metallic ternary prototype (ICSD #160913)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oI10_71_f_h_a"){
        vparameters.push_back("8.5967,0.66765154,0.51107983,0.7969,0.2835");  // 001, metallic ternary prototype (ICSD #55495)
        vparameters.push_back("8.994,0.60973983,0.47520569,0.1961,0.2537");  // 002, metallic ternary prototype (ICSD #58885)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C3_oI14_71_e_f_af"){
        vparameters.push_back("19.795,0.24976004,0.22556201,0.4351,0.8584,0.6832");  // 001, metallic ternary prototype (ICSD #9564)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C_oI16_71_an_e_d"){
        vparameters.push_back("9.562,0.73457436,0.41612633,0.758,0.354,0.31");  // 001, metallic ternary prototype (ICSD #58046)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3C4_oI22_71_n_af_eg"){
        vparameters.push_back("16.5484,0.51916802,0.27410505,0.80940162,0.81788308,0.65938423,0.61795128,0.69904145");  // 001, metallic ternary prototype (ICSD #55826)
        vparameters.push_back("15.489,0.47588611,0.30221447,0.204,0.152,0.31,0.337,0.322");  // 002, metallic ternary prototype (ICSD #55827)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3C4_oI22_71_n_af_eh"){
        vparameters.push_back("15.1225,0.48197719,0.29965945,0.2162,0.8716,0.3059,0.3285,0.6944");  // 001, metallic ternary prototype (ICSD #156968)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C4_oI22_71_af_eh_n"){
        vparameters.push_back("14.407,0.51613799,0.2911779,0.78461,0.87104,0.6954,0.67309,0.6906");  // 001, metallic ternary prototype (ICSD #182774)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B6C_oI26_71_hk_efg_a"){
        vparameters.push_back("8.613,0.97886915,0.58156275,0.655,0.665,0.336,0.757");  // 001, metallic ternary prototype (ICSD #103465)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C8_oI28_71_cg_be_2n"){
        vparameters.push_back("12.707,0.7480916,0.33721571,0.201,0.181,0.369,0.282,0.153,0.367");  // 001, metallic ternary prototype (ICSD #605107)
        vparameters.push_back("12.373,0.80029096,0.34009537,0.799,0.181,0.631,0.282,0.847,0.367");  // 002, metallic ternary prototype (ICSD #611824)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B10C_oI28_71_be_g2n_c"){
        vparameters.push_back("15.007,0.6622243,0.34248018,0.20187,0.8153,0.35777,0.71484,0.15945,0.64223");  // 001, metallic ternary prototype (ICSD #290305)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B5C3_oI28_71_gn_cn_be"){
        vparameters.push_back("12.827,0.76666407,0.34029781,0.19179,0.19303,0.3783,0.28322,0.16931,0.34743");  // 001, metallic ternary prototype (ICSD #413485)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B3C2_oI28_71_c2n_be_g"){
        vparameters.push_back("12.442,0.75952419,0.34078122,0.799,0.819,0.631,0.718,0.847,0.633");  // 001, metallic ternary prototype (ICSD #634397)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B3C3_oI28_71_2n_be_cg"){
        vparameters.push_back("12.1749,0.81103746,0.34336216,0.82977124,0.78536152,0.62439285,0.71382786,0.83444541,0.62824919");  // 001, metallic ternary prototype (ICSD #634379)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oI20_72_j_j_a"){
        vparameters.push_back("14.325,0.50820244,0.49486911,0.86413,0.70673,0.61373,0.71129");  // 001, metallic ternary prototype (ICSD #421424)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_tI10_107_a_a_ab"){
        vparameters.push_back("4.82,2.2676349,0,0.6654,0.4241,0.7494");  // 001, metallic ternary prototype (ICSD #58662)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_tI10_107_ab_a_a"){
        vparameters.push_back("4.2357,2.5098095,0.405,0.635,0,0.748");  // 001, metallic ternary prototype (ICSD #290390)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_tI10_119_bf_a_c"){
        vparameters.push_back("4.6895,2.6930376,0.8608");  // 001, metallic ternary prototype (ICSD #249592)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_tI28_120_d_i_e"){
        vparameters.push_back("7.753,1.014962,0.856,0.648,0.625,0.417");  // 001, metallic ternary prototype (ICSD #102238)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B6C_tI22_121_i_ci_a"){
        vparameters.push_back("6.84,1.4307018,0.3274,0.9228,0.3237,0.2067");  // 001, metallic ternary prototype (ICSD #54354)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_tP4_123_b_h_a"){
        vparameters.push_back("2.9,2.4965517,0.23");  // 001, metallic ternary prototype (ICSD #102057)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tP6_123_b_i_a"){
        vparameters.push_back("4.103,1.5781136,0.3");  // 001, metallic ternary prototype (ICSD #623081)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC5_tP7_123_a_b_ci"){
        vparameters.push_back("4.601,1.6387742,0.3103");  // 001, metallic ternary prototype (ICSD #102108)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C_tP7_123_a_ci_b"){
        vparameters.push_back("4.549,1.6322269,0.805");  // 001, metallic ternary prototype (ICSD #623944)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C_tP8_123_a_hi_b"){
        vparameters.push_back("4.35,1.8211494,0.15149,0.32935");  // 001, metallic ternary prototype (ICSD #240161)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C2_tP11_123_a_ehi_g"){
        vparameters.push_back("4.217,2.6013754,0.3093,0.2933,0.1189");  // 001, metallic ternary prototype (ICSD #42426)
        vparameters.push_back("4.23,2.3841608,0.3104,0.2925,0.1138");  // 002, metallic ternary prototype (ICSD #180132)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tP10_127_g_a_h"){
        vparameters.push_back("7.0659,0.48323639,0.374,0.169");  // 001, metallic ternary prototype (ICSD #54303)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tP10_127_a_g_h"){
        vparameters.push_back("5.716,0.55563331,0.389,0.181");  // 001, metallic ternary prototype (ICSD #57706)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_129_b_c_c"){
        vparameters.push_back("4.7819,1.6426943,0.72227037,0.16905509");  // 001, metallic ternary prototype (ICSD #66951)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_129_c_c_a"){
        vparameters.push_back("4.741,1.7527948,0.201,0.65");  // 001, metallic ternary prototype (ICSD #601586)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_129_a_c_c"){
        vparameters.push_back("4.9,1.6857143,0.273,0.72");  // 001, metallic ternary prototype (ICSD #615958)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tP8_129_bc_c_a"){
        vparameters.push_back("4.3365,2.2772743,0.1843,0.7494");  // 001, metallic ternary prototype (ICSD #163425)
        vparameters.push_back("4.5206,2.4320223,0.16198,0.73829");  // 002, metallic ternary prototype (ICSD #415728)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tP10_129_c_ac_bc"){
        vparameters.push_back("4.401,2.2842536,0.2446,0.6257,0.8776");  // 001, metallic ternary prototype (ICSD #160053)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tP10_129_c_bc_ac"){
        vparameters.push_back("4.528,2.5110424,0.75271,0.13526,0.3628");  // 001, metallic ternary prototype (ICSD #162267)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tP10_129_bc_c_ac"){
        vparameters.push_back("4.546,2.1240651,0.135,0.76,0.38");  // 001, metallic ternary prototype (ICSD #616562)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tP14_129_ac_b_j"){
        vparameters.push_back("7.196,0.72290161,0.7755,0.43525,0.29");  // 001, metallic ternary prototype (ICSD #418547)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_tI8_139_b_d_a"){
        vparameters.push_back("4.137,1.7689147");  // 001, metallic ternary prototype (ICSD #54365)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tI10_139_d_a_e"){
        vparameters.push_back("4.242,2.6692834,0.38820402");  // 001, metallic ternary prototype (ICSD #55789)
        vparameters.push_back("4.127,2.8279622,0.6073");  // 002, metallic ternary prototype (ICSD #57550)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tI12_139_c_e_e"){
        vparameters.push_back("4.4914,3.6117024,0.66644,0.8675");  // 001, metallic ternary prototype (ICSD #182479)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tI12_139_e_c_e"){
        vparameters.push_back("4.359,3.6237669,0.333,0.135");  // 001, metallic ternary prototype (ICSD #183524)
        vparameters.push_back("4.073,3.1723545,0.131,0.341");  // 002, metallic ternary prototype (ICSD #611768)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_tI14_139_h_d_a"){
        vparameters.push_back("6.717,0.79082924,0.303");  // 001, metallic ternary prototype (ICSD #456)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_tI14_139_a_h_d"){
        vparameters.push_back("7.606,0.71916908,0.3");  // 001, metallic ternary prototype (ICSD #630579)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C_tI16_139_bg_e_a"){
        vparameters.push_back("4.023,3.5893612,0.762,0.851");  // 001, metallic ternary prototype (ICSD #58084)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tI16_139_ce_e_d"){
        vparameters.push_back("4.64,4.7327586,0.3248,0.1233");  // 001, metallic ternary prototype (ICSD #41924)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_tI16_139_e_ce_d"){
        vparameters.push_back("4.768,4.9496644,0.1158,0.3296");  // 001, metallic ternary prototype (ICSD #58635)
        vparameters.push_back("4.846,4.5356995,0.8717,0.6859");  // 002, metallic ternary prototype (ICSD #58638)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tI16_139_ce_d_e"){
        vparameters.push_back("4.635,4.9363538,0.6636,0.8862");  // 001, metallic ternary prototype (ICSD #58764)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C3_tI20_139_e_bg_ae"){
        vparameters.push_back("4.9461,4.7316067,0.19934,0.64091,0.57874");  // 001, metallic ternary prototype (ICSD #421342)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8BC4_tI26_139_ij_a_f"){
        vparameters.push_back("8.823,0.58381503,0.6526,0.7244");  // 001, metallic ternary prototype (ICSD #57539)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B8C_tI26_139_i_fj_a"){
        vparameters.push_back("8.301,0.56161908,0.3594,0.7234");  // 001, metallic ternary prototype (ICSD #168240)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B4C_tI26_139_ij_f_a"){
        vparameters.push_back("8.86,0.5744921,0.365,0.175");  // 001, metallic ternary prototype (ICSD #57997)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_tI24_140_a_c_l"){
        vparameters.push_back("6.5562,1.721424,0.1575,0.867");  // 001, metallic ternary prototype (ICSD #172149)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tI28_140_h_a_k"){
        vparameters.push_back("10.491,0.47450195,0.861,0.8027,0.0827");  // 001, metallic ternary prototype (ICSD #150145)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_tI28_140_l_h_a"){
        vparameters.push_back("8.175,0.95425076,0.33913,0.36404,0.66158");  // 001, metallic ternary prototype (ICSD #182104)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_tI32_140_a_cl_h"){
        vparameters.push_back("7.7012,1.7607516,0.36124,0.65581,0.63793");  // 001, metallic ternary prototype (ICSD #107217)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C_tI32_140_ah_l_c"){
        vparameters.push_back("8.565,1.4824285,0.85,0.15,0.34");  // 001, metallic ternary prototype (ICSD #616725)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_tI32_140_a_bk_h"){
        vparameters.push_back("10.586,0.48904213,0.832,0.7796,0.0724");  // 001, metallic ternary prototype (ICSD #103842)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC5_tI32_140_h_a_cl"){
        vparameters.push_back("7.9856,1.6805249,0.1336,0.6569,0.8577");  // 001, metallic ternary prototype (ICSD #156956)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C5_tI32_140_a_h_cl"){
        vparameters.push_back("7.4816,1.855726,0.8881,0.6658,0.8611");  // 001, metallic ternary prototype (ICSD #161658)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B7C3_hR14_146_ab_a2b_b"){
        vparameters.push_back("11.53099,0.4552076,0.944,0.46,0.352,0.533,0.121,0.408,0.149,0.972,0.617,0.864,0.026,0.629,0.461,0.903");  // 001, metallic ternary prototype (ICSD #103911)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_156_3a_b2c_2bc"){
        vparameters.push_back("4.94,2.2064777,0.334,0.668,0,0.454,0.146,0.853,0.214,0.788,0.5");  // 001, metallic ternary prototype (ICSD #58911)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5BC4_hR10_160_5a_a_4a"){
        vparameters.push_back("4.7099904,6.7346082,0.94,0.862,0.409,0.061,0.145,0.501,0.7787,0.6863,0.2327,0.3215");  // 001, metallic ternary prototype (ICSD #12142)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_hR4_166_c_b_a"){
        vparameters.push_back("4.2709972,3.7286311,0.1888");  // 001, metallic ternary prototype (ICSD #102938)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_hR7_166_a_c_2c"){
        vparameters.push_back("4.2259957,6.2863209,0.90091,0.79965,0.61204");  // 001, metallic ternary prototype (ICSD #150127)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C6_hR11_166_a_2c_h"){
        vparameters.push_back("5.9937387,3.892206,0.66786054,0.87208523,0.39951361,0.88682193");  // 001, metallic ternary prototype (ICSD #166193)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_hR12_166_h_bc_ac"){
        vparameters.push_back("5.4710243,4.6349773,0.335,0.142,0.593,0.054");  // 001, metallic ternary prototype (ICSD #604213)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B3C2_hR12_166_bh_ac_c"){
        vparameters.push_back("5.6140311,4.6050915,0.1426,0.335,0.5787,0.0745");  // 001, metallic ternary prototype (ICSD #57538)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C9_hR12_166_a_c_bch"){
        vparameters.push_back("4.9234811,4.8473628,0.146,0.334,0.5872,0.0827");  // 001, metallic ternary prototype (ICSD #55614)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B2C_hR12_166_eh_c_a"){
        vparameters.push_back("5.1579625,3.6593664,0.222,0.278,0.777");  // 001, metallic ternary prototype (ICSD #646978)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hR12_166_bc_h_ac"){
        vparameters.push_back("5.6320019,4.8020254,0.333,0.14,0.586,0.074");  // 001, metallic ternary prototype (ICSD #604688)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C3_hR12_166_c_bh_ac"){
        vparameters.push_back("5.6139959,4.7716445,0.6644,0.8571,0.419,0.908");  // 001, metallic ternary prototype (ICSD #104173)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B4C2_hR13_166_ah_2c_c"){
        vparameters.push_back("6.8009851,4.2773178,0.3152,0.1878,0.5493,0.2133,0.8046");  // 001, metallic ternary prototype (ICSD #58581)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C2_hR14_166_3c_2c_abc"){
        vparameters.push_back("4.7683821,9.8351827,0.93576,0.79728,0.64381,0.87437,0.70884,0.57697");  // 001, metallic ternary prototype (ICSD #710044)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_186_b_a_b"){
        vparameters.push_back("4.7199,1.6142715,0.75,0.0593,0.473");  // 001, metallic ternary prototype (ICSD #54997)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_186_b_a_b"){
        vparameters.push_back("4.576,1.6121066,0,0.382,0");  // 002, metallic ternary prototype (ICSD #100115)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_186_a_b_b"){
        vparameters.push_back("4.655,1.6386681,0.247,0.467,0.0518");  // 001, metallic ternary prototype (ICSD #106358)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_hP12_186_ac_b_b"){
        vparameters.push_back("4.9708,1.5973284,0,0.9386,0.5644,0.83026667,0.7503");  // 001, metallic ternary prototype (ICSD #424277)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_hP16_186_a2b_b_2a2b"){
        vparameters.push_back("4.4831,4.5854431,0.16473,0.02937,0.29581,0.98426,0.84576,0.711,0.12513,0.40817");  // 001, metallic ternary prototype (ICSD #412207)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_187_i_h_ab"){
        vparameters.push_back("4.3442,1.5409512,0.7853,0.7161");  // 001, metallic ternary prototype (ICSD #156264)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_hP6_187_hi_b_a"){
        vparameters.push_back("4.3168,1.5918273,0.30839,0.22249");  // 001, metallic ternary prototype (ICSD #98666)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_187_dh_ai_eg"){
        vparameters.push_back("4.933,2.2278532,0.6533,0.8342,0.708");  // 001, metallic ternary prototype (ICSD #409533)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_187_eh_ci_bg"){
        vparameters.push_back("4.918,2.2206181,0.8334,0.696,0.6544");  // 001, metallic ternary prototype (ICSD #409534)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B5C_hP11_187_hk_cgi_a"){
        vparameters.push_back("6.0692,1.8322019,0.32828,0.322,0.2,0.5122");  // 001, metallic ternary prototype (ICSD #416337)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B6C_hP11_187_ck_jk_a"){
        vparameters.push_back("8.0581,0.55113488,0.53424,0.79441,0.19539");  // 001, metallic ternary prototype (ICSD #249520)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C3_hP8_189_f_c_g"){
        vparameters.push_back("7.2873,0.42088291,0.76243633,0.41574233");  // 001, metallic ternary prototype (ICSD #410967)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_g_f_ad"){
        vparameters.push_back("7.579,0.51801029,0.60163,0.26746");  // 001, metallic ternary prototype (ICSD #51845)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C6_hP9_189_b_c_fg"){
        vparameters.push_back("7.801,0.43276503,0.249,0.592");  // 001, metallic ternary prototype (ICSD #20876)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3C_hP9_189_cf_g_b"){
        vparameters.push_back("6.759,0.61236869,0.24,0.578");  // 001, metallic ternary prototype (ICSD #634575)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_bc_g_f"){
        vparameters.push_back("7.1941,0.54156879,0.755,0.415");  // 001, metallic ternary prototype (ICSD #54295)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC6_hP9_189_c_b_fg"){
        vparameters.push_back("8.213,0.50651406,0.235,0.601");  // 001, metallic ternary prototype (ICSD #96253)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_g_bc_f"){
        vparameters.push_back("7.717,0.49941687,0.597,0.26");  // 001, metallic ternary prototype (ICSD #54315)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_f_bc_g"){
        vparameters.push_back("7.637,0.61332984,0.75,0.4277");  // 001, metallic ternary prototype (ICSD #54344)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_f_ad_g"){
        vparameters.push_back("7.4992,0.62021282,0.5659,0.243");  // 001, metallic ternary prototype (ICSD #55102)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_g_f_bc"){
        vparameters.push_back("7.4174,0.50893844,0.2522,0.5889");  // 001, metallic ternary prototype (ICSD #55846)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_f_bc_g"){
        vparameters.push_back("6.798,0.46587232,0.273,0.615");  // 002, metallic ternary prototype (ICSD #103885)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_bc_g_f"){
        vparameters.push_back("7.751,0.48651787,0.245,0.585");  // 002, metallic ternary prototype (ICSD #105713)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_f_g_bc"){
        vparameters.push_back("7.0015,0.54982504,0.793,0.446");  // 001, metallic ternary prototype (ICSD #107416)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_ad_g_f"){
        vparameters.push_back("7.345,0.6250919,0.56643333,0.24173333");  // 001, metallic ternary prototype (ICSD #160908)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_f_g_bc"){
        vparameters.push_back("7.0941,0.50138848,0.76,0.6667");  // 002, metallic ternary prototype (ICSD #608141)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_f_g_bc"){
        vparameters.push_back("7.016,0.5791049,0.24,0.3333");  // 003, metallic ternary prototype (ICSD #608770)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_bc_f_g"){
        vparameters.push_back("7.45,0.53288591,0.75506667,0.51506667");  // 001, metallic ternary prototype (ICSD #628137)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_g_bc_f"){
        vparameters.push_back("7.286,0.49725501,0.25,0.3333");  // 002, metallic ternary prototype (ICSD #638586)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C_hP6_191_g_c_a"){
        vparameters.push_back("5.565,0.81239892");  // 001, metallic ternary prototype (ICSD #57329)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C_hP6_191_g_c_a"){
        vparameters.push_back("8.813,0.46102349");  // 002, metallic ternary prototype (ICSD #102452)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_hP6_191_c_a_g"){
        vparameters.push_back("7.44,0.52956989");  // 001, metallic ternary prototype (ICSD #108369)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_hP6_191_c_a_g"){
        vparameters.push_back("5.581,0.70757929");  // 002, metallic ternary prototype (ICSD #658142)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9BC2_hP12_191_fm_a_c"){
        vparameters.push_back("8.04,0.48383085,0.212");  // 001, metallic ternary prototype (ICSD #57518)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC6_hP13_191_i_a_cde"){
        vparameters.push_back("5.48,1.6394161,0.6683,0.7529");  // 001, metallic ternary prototype (ICSD #54273)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B6C_hP13_191_i_cde_a"){
        vparameters.push_back("5.5301,1.6316522,0.168,0.746");  // 001, metallic ternary prototype (ICSD #54277)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C2_hP16_193_g_g_d"){
        vparameters.push_back("7.9396,0.68974508,0.60359025,0.24573816");  // 001, metallic ternary prototype (ICSD #103733)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_d_a_c"){
        vparameters.push_back("4.713,1.8949714");  // 001, metallic ternary prototype (ICSD #54319)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_a_c_d"){
        vparameters.push_back("4.953,1.8465576");  // 001, metallic ternary prototype (ICSD #56978)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_c_a_d"){
        vparameters.push_back("4.553,1.7665276");  // 001, metallic ternary prototype (ICSD #57018)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_a_d_c"){
        vparameters.push_back("4.726,1.9970377");  // 001, metallic ternary prototype (ICSD #106303)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_a_d_c"){
        vparameters.push_back("4.909,1.500713");  // 002, metallic ternary prototype (ICSD #106345)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_a_c_d"){
        vparameters.push_back("4.583,1.7161248");  // 002, metallic ternary prototype (ICSD #106384)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_c_d_a"){
        vparameters.push_back("4.922,1.4691182");  // 002, metallic ternary prototype (ICSD #108545)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_c_a_d"){
        vparameters.push_back("4.5813,1.6472617");  // 002, metallic ternary prototype (ICSD #600994)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_a_c_d"){
        vparameters.push_back("4.783,1.6048505");  // 003, metallic ternary prototype (ICSD #622301)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_194_b_f_c"){
        vparameters.push_back("4.267,2.0180455,0.58");  // 001, metallic ternary prototype (ICSD #58139)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_194_a_f_c"){
        vparameters.push_back("4.539,1.9649703,0.592");  // 001, metallic ternary prototype (ICSD #59505)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_194_c_f_a"){
        vparameters.push_back("4.498,1.9873277,0.59");  // 001, metallic ternary prototype (ICSD #103292)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_hP8_194_f_a_c"){
        vparameters.push_back("4.303,1.7748083,0.583");  // 001, metallic ternary prototype (ICSD #150602)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_194_c_f_a"){
        vparameters.push_back("4.6571,1.8583024,0.5887");  // 002, metallic ternary prototype (ICSD #157110)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC2_hP12_194_h_a_f"){
        vparameters.push_back("5.27,1.626186,0.063,0.833");  // 001, metallic ternary prototype (ICSD #58159)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP12_194_f_f_ab"){
        vparameters.push_back("4.552,4.1001757,0.1617,0.6162");  // 001, metallic ternary prototype (ICSD #66003)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP12_194_ab_f_f"){
        vparameters.push_back("4.4951,3.5496874,0.11481,0.65681");  // 001, metallic ternary prototype (ICSD #152621)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP12_194_f_ab_f"){
        vparameters.push_back("4.468,3.5505819,0.61414,0.15824");  // 001, metallic ternary prototype (ICSD #152622)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP12_194_ab_f_f"){
        vparameters.push_back("4.588,3.7149085,0.6154,0.1507");  // 002, metallic ternary prototype (ICSD #152623)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C2_hP16_194_h_af_f"){
        vparameters.push_back("5.9892,2.4310425,0.12318,0.612,0.84546667");  // 001, metallic ternary prototype (ICSD #415931)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_cP16_215_e_bce_ad"){
        vparameters.push_back("6.315,0.7511,0.2331");  // 001, metallic ternary prototype (ICSD #180119)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C4_cI22_217_b_c_c"){
        vparameters.push_back("7.41,0.2031,0.366");  // 001, metallic ternary prototype (ICSD #58899)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_cP15_221_ag_c_d"){
        vparameters.push_back("6.408,0.29015");  // 001, metallic ternary prototype (ICSD #412144)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_cP15_221_d_ag_c"){
        vparameters.push_back("6.384,0.711");  // 001, metallic ternary prototype (ICSD #629734)
      }

    }
    if(library=="" && !keep_original_lattice_parameter){
      //DX 20190314 - loop over parameters - START
      for(uint p=0;p<vparameters.size();p++){
        aurostd::string2tokens(vparameters[p],tokens,",");
        tokens[0]="-1";
        vparameters[p]=aurostd::joinWDelimiter(tokens,",");
      }
      //DX 20190314 - loop over parameters - END
    }
    if(library != "all"){
      if(library=="" && choice!=-1){
        if((uint)choice<vparameters.size()){
          vector<string> tmp; tmp.push_back(vparameters[choice]);
          vparameters.clear();
          vparameters=tmp;
        }
        else{
          message << "anrl::getANRLParameters(): ERROR - " << anrl_label << " does not have more than " << vparameters.size() << " choice(s):";
          for(uint i=0;i<vparameters.size();i++){
            message << "  " << anrl_label << "-" << std::setw(3) << std::setfill('0') << i+1 << " : " << vparameters[i] << endl;
          }
          throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name, message, _VALUE_RANGE_); //DX 20191118 - exit to throw
        }
      }
      else if((library=="" && choice==-1) || (vparameters.size() && (library=="part1" || library=="part2" || library=="misc"))){
        // -------------------------------------------------------------------------
        // we do not automatically assign the enumeration suffix to the label
        // (even if there there is currently only one prototype with a label)
        // otherwise, the command without the suffix may break in the future
        // (cannot auto-assign when more than one choice, e.g., -001 or -002)
        message << "anrl::getANRLParameters(): ERROR - " << anrl_label << " has " << vparameters.size() << " preset parameter set(s): " << endl;
        for(uint i=0;i<vparameters.size();i++){
          message << "  " << anrl_label << "-" << std::setw(3) << std::setfill('0') << i+1 << " : " << vparameters[i] << endl;
        }
        message << "Rerun command and specify the parameters or the preset suffix, e.g., aflow --proto=" << anrl_label << "-" << std::setw(3) << std::setfill('0') << 1; //DX 20190826 - changed "./aflow" to "aflow"
        throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name, message, _VALUE_ERROR_); //DX 20191118 - exit to throw
      }
    }
    return vparameters;  
  }
}

#endif // AFLOW_REMOVE_GREP // _AFLOW_ANRL_LIST_CPP // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************

