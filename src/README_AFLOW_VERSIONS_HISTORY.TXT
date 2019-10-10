3.1.227 - 2019/11/01
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.227.tar.xz
        -- Fixed bug in xDOSCAR reader error detection. (ME)
        -- Additional web output for entry page applets. (ME)
3.1.226 - 2019/10/01
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.226.tar.xz
        -- Modification on the allocation of web entries, through web cash of AUID. Speed up of 20-200 times in serving aflowlib entry pages (SC)
        -- Bug fixes in producing html/json for icsd# entries (SC)
        -- Modified aflow.in generation in AEL and AGL to use common routines in aflow_avasp.cpp and aflow_modules.cpp (CT)  
        -- Write AGL enthalpy as a function of pressure in an AFLOW readable format for future phase diagram application (CT)
        -- Added --quiet option to LIB2LIB routine (CT)
        -- Sets the VASP symmetry to off by default for AEL calculations (CT)
        -- Added projected phonon DOS calculations. (ME)
        -- Linear tetrahedron method does not reduce upon construction anymore (necessary for projected DOS). (ME)
        -- Cleaned up code and added more comments to QMesh and LTMethod files. (ME)
        -- Made moveToBZ inside QMesh more robust. (ME)
        -- Replaced all instances of exit in APL with throw APLRuntimeError. (ME)
        -- Sped up reduction of q-mesh grid. (ME)
        -- Added safeguards for DOSCAR readers for broken DOSCAR files. (ME)
        -- Fixed plotter bugs for HTQC prototype names. (ME)
        -- Added findInList function to AUROSTD. (ME)
        -- Generalized BringInCell method. (ME/DX)
        -- created overloads for BringInCell method (double/xvector/_atom/xstructure), along with "InPlace" variants, i.e., const input (DX)
        -- replaced all mod_one/mod_one_xvec/mod_one_atom for new BringInCell method (DX)
        -- replaced old BringInCell functions with new version (DX)
        -- fixed indenting in symmetry functions (DX)
        -- replaced all SYM::CrossPro() with aurostd::vector_product() (DX)
        -- replaced all SYM::normalize() with explicity division of aurostd::modulus() (DX)
        -- replaced SYM::extract_row() with aurostd::operator() (DX)
        -- fabs() -> aurostd::abs() (CO)
        -- changed abinit writer lattice from column-based to row-based (CO)
        -- updated some class definitions (free()/copy()) (CO)
        (AFLOW: aconvasp.cpp, aconvasp_main.cpp, aflow.cpp, aflow.h, aflowlib.h, aflowlib_libraries.cpp, aflowlib_webapp_entry.js, aflowlib_web_interface.cpp, aflowrc.cpp, apennsy_gndstate.cpp, apennsy_main.cpp, apennsy_vaspin.cpp, avasp.cpp, chull.cpp, compare_structure_function.cpp, compare_structure.h, contrib_cormac.cpp, contrib_wahyu.cpp, data.cpp, init.cpp, Makefile, modules.cpp, oaims.cpp, ovasp.cpp, pflow_funcs.cpp, pflow_print.cpp, plotter.cpp, pocc.cpp, pocc.h, pocc_old.cpp, README_AFLOW_ACONVASP.TXT, README_AFLOW_SYM.TXT, spacegroup.cpp, surface.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_functions.cpp, symmetry_spacegroup.h, symmetry_spacegroup_ITC_library.cpp, xatom.cpp, xclasses.cpp, xproto.cpp)
        (AEL: ael_elastic_fit.cpp, ael_elasticity.cpp, ael_elasticity.h, ael_get_stress.cpp)
        (AGL: agl_debye.cpp, agl_debye.h, agl_electronic.cpp, agl_eqn_state.cpp, agl_get_ev.cpp, agl_hugoniot.cpp, agl_polynomial.cpp, agl_rungibbs.cpp, agl_thermal.cpp)
        (APL: apl_atomic_disp.cpp, apl_dm_pdos_save.cpp, apl_doscalc.cpp, apl_function_fitting.cpp, apl_group_velocity.cpp, apl.h, apl_hsqpoints.cpp, apl_kphonons.cpp, apl_ltet.cpp, apl_pdisc.cpp, apl_phoncalc.cpp, apl_qmesh.cpp, apl_shellhandle.cpp, apl_supercell.cpp, qha3phonons_eos.cpp, qha_dm_pdos_save.cpp, qha_energies.cpp, qha_eoscalc.cpp, qha_gruneisen.cpp, qha_operations.cpp, README_AFLOW_APL.TXT, scqha_eos.cpp, scqha_gruneisen.cpp)
        (AUROSTD: aurostd.h, main.cpp, xmatrix.cpp, xmatrix.h)
3.1.225 - 2019/07/26
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.225.tar.xz
        -- fixing bug with stability criterion calculation in --output=web, need structure comparison for correct stability criterion calculation (WS)
        -- changed output name for --output=web to aflow_ALLOY_hull_web.json (WS)
        -- integrated N+1_enthalpy_gain and stability criterion to --output=web (WS)
        -- patched web json for quaternary chull (WS)
        -- added CCE functionality (RF)
        -- added CCE readme (RF)
        -- added a user input cutoff energy to the GFA code (DF)
        -- Added vector version of RemoveComments. (ME)
        -- RemoveControlCharactersFromStringstream now also removes carriage returns because they break string additions. (ME)
        -- SPIN_REMOVE_RELAX_2 did not remove spin after the second relaxation; this is fixed now. (ME)
        -- Removed bug that prevented certain NBANDS errors from being fixed while an NBANDS warning was present in vasp.out (ME)
        -- APL and AAPL now read forces from qmvasp file. (ME)
        -- AAPL_KPPRA added to aflow_modules. (ME)
        -- Minor APL code clean-ups. (ME)
        -- Added new plotter functions for band structures, phonons, and thermal properties. (ME)
        -- Updated numbering scheme for partial DOS from being zero-based to being one-based. (ME)
        -- Changed xDOSCAR, xEIGENVAL, and xKPOINTS to be compatible with new plotter. (ME)
        -- Added functionality to APL that outputs phonon DOS and dispersion in
        -- VASP format for new plotter functions. (ME)
        -- Fixed bug for the path builder for path with - and | as separators. (ME)
        -- Fixed xDOSCAR reader to be more general. (ME)
        -- patched change_tolerance() bug for double& tolerance return (DX)
        -- renamed minimum distance functions, more descriptive (minimizeDistanceCartesianMethod() and minimizeDistanceFractionalMethod()) (DX)
        -- reduced number of inputs for minimum distance functions (DX)
        -- added functionality to read VASP5 format in CHGCAR (DX)
        -- cleaned species/atom names in xstructure2json; remove pseudopotential information, oxidation numbers, etc. (DX)
        -- added boolean to ignore minimum distance check in foldAtomsInCell() (DX)
        -- added overloads for mapping functions (DX) 
        -- fixed bug in in xstructure2json in LIB2RAW runs; store corresponding structure from relaxations (DX) 
        -- read species for xstructure2json in LIB2RAW runs (DX) 
        -- added default permutation designation to ANRL prototypes for --proto and --aflow_proto, e.g., ".AB" for "AB_cF8_225_a_b.AB" (DX) 
        -- added functionality to read "_symmetry_space_group_name_H-M" in CIF reader (DX)
        -- added GetSpaceGroupNumber() to get the space group number from the Hermann-Mauguin symbol (DX)
        -- added RemoveCharacterFromBack(), RemoveCharacterFromFront(), and RemoveCharacterFromBackandFront() to AUROSTD; useful for parsing CIF fields (DX)
        -- added function to compare Wyckoff position strings regardless of order, necessary for matching equations in CIF to ITC standard (DX)
        -- fixed bug in reading symops from CIF; check for symop_id (DX)
        -- fixed bug in "_space_group_name_Hall" reader in CIF (DX)
        -- converted exits to xerrors in GetSpaceGroup*() functions (DX)
        -- added DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF to aflow.rc (CO)
        -- added functions to find pointLineIntersection() and linePlaneIntersection() (CO)
        -- added function to find 3D rotation between two 3D vectors (CO)
        -- remapped const utype& tol -> utype tol in function inputs (CO)
        -- remapped const bool& -> bool in function inputs (CO)
        -- added check in generalHouseHolderQRDecomposition() that Q*R=A (CO)
        -- changed gcd() to GCD() (CO)
        -- added least common multiple function LCM() (CO)
        -- added --slab_test (CO)
        -- added setctau/setftau functions to _sym_op (CO)
        -- added ReplaceAtoms() to xstructure to safely RemoveAtom() then AddAtom() (CO)
        -- modified foldAtomsInCell() to remove fold_in_only complexity, this was simplified by GenerateGridAtoms() (CO)
        -- added check to foldAtomsInCell() that min_dist does not change (CO)
        -- added CreateSlab_RigidRotation(), CreateSlab_SurfaceLattice(), and associated functions (many overloads) to create slabs for surface calculations: the two functions should yield about the same structure but by two different methods (CO)
        -- reformulate basis in CreateSlab_RigidRotation() so k-points grid is created appropriately (K1xK2x1) (CO)
        -- added new modules for calculating the generalized stacking fault energy and cleavage energy (CO)
        -- patched bug in CHULL pdf output for missing column-types in compiling aflow_CsIPb_hull.tex (CO)
        -- added cout as default for ostream& oss in header files (CO)
        -- converted Normal2HKLPlane (CO)
        -- reduce total layers creation by calculating how many layers are in a cell for that particular direction (CO)
        -- added XHOST.WEB_MODE to reduce cluttered AFLOW-Online output (CO)
        -- removed pocc progress bar in web_mode (CO)
        -- added FPOSDistance() for skew considerations (CO)
        -- reduced dims considered for GetSuperCell(), RadiusSphereLattice and LatticeDimensionSphere() are enough (CO)
        -- added general getPeaks() function to AUROSTD (CO)
        -- calculated X-ray amplitude on-the-fly, no need to pass it around (CO)
        -- moved double lambda to end of inputs and default to 1.5418 Copper K-alpha (CO)
        -- patched RemoveControlCodeCharactersFromString() for carriage return (CO)
        -- added width_maximum to getPeaks() to resolve arbitrarily defined peaks (CO)
        -- made x-ray analysis code more modular (CO)
        -- added --plot_xray_file option (CO)
        -- patched moving_average function (!= convolution with box_filter) (CO)
        -- defaulted lambda=Cu_Ka for XRD functions (CO)
        -- added compareVecElements(), compareXVecElements(), and class compareVecElement(ind) to sort vector<vector<utype> > and vector<xvector<utype> > and retire classes ids_cmp and hkl_cmp (CO)
        -- changed signal -> signal_input so it does not get confused with std::signal (CO)
        -- added --force_generic_title to PLOT_XRAY (CO)
        -- removed ids_cmp and hkl_cmp classes (CO)
        -- added filetype and vector_reduction_type enums (CO)
        -- patched aflags.Directory() early in kbin::KBIN_Main() (CO)
        -- corrected GetCompoundAttenuationLenght -> GetCompoundAttenuationLength and toogle -> toggle (CO)
        -- extended PARTCAR functionality for VASP5 POSCARs, including the optional Selective Dynamics tag (CO)
        -- patched 'make check' to return an error, instead of a string, if encountered (CO)
        -- patched 'S' -> 'P' conversion in --pocc_params (CO)
        -- incorporated some new warnings from VASP6 (CO)
        -- patched stringElements2VectorElements() to handle two types of compound designation: composition_string and pp_string (CO)
        -- added new suffixes to pseudopotential list (CO)
        -- address capital letters in new pp suffixes directly (CO)
        -- added faster (in place) routines for CleanStringASCII(), RemoveNumbers(), VASP_PseudoPotential_CleanName() (CO)
        -- added AFLOWIN_FLAG::NOAUTOPP (--noautopp) functionality to --aflow_proto (CO)
        -- converted flag PFLOW::LOAD_ENTRIES_ONLY_ALPHABETICAL to PFLOW::LOAD_ENTRIES_NON_ALPHABETICAL so the default is always OFF (CO)
        -- added 'make clean_autogen' to delete files that are automatically downloaded/generated by aflow during compilation (CO)
        -- set explicit double->int conversion for roundDouble() in chull (CO)
        (AFLOW: aconvasp.cpp, aconvasp_main.cpp, aflow.cpp, aflow.h, aflowlib_libraries.cpp, aflowlib_web_interface.cpp, aflowrc.cpp, avasp.cpp, cce.cpp, cce.h, chull.cpp, chull.h, compare_structure_function.cpp, contrib_shidong_main.cpp, data.cpp, gfa.cpp, gfa.h, init.cpp, ivasp.cpp, kbin.cpp, kvasp.cpp, Makefile, mix.cpp, modules.cpp, ovasp.cpp, pflow_funcs.cpp, pflow.h, pflow_print.cpp, plotter.cpp, pocc.cpp, README_AFLOW_ACONVASP.TXT, README_AFLOW_AFLOWRC.TXT, README_AFLOW_CCE.TXT, README_AFLOW_GFA.TXT, surface.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_functions.cpp, symmetry_spacegroup.h, symmetry_spacegroup_ITC_library.cpp, xatom.cpp, xclasses.cpp, xelement.cpp, xelement.h, xproto.cpp)
        (AGL: agl_electronic.cpp)
        (APL: aapl_cluster.cpp, aapl_tcond.cpp, apl_dirphoncalc.cpp, apl_doscalc.cpp, apl.h, apl_kphonons.cpp, apl_lrphoncalc.cpp, apl_pathbuilder.cpp, apl_pdisc.cpp, apl_phoncalc.cpp, apl_qmesh.cpp, apl_supercell.cpp, apl_thermalpc.cpp, README_AFLOW_APL.TXT)
        (AUROSTD: aurostd.h, boot.cpp, main.cpp, xmatrix.cpp, xmatrix.h, xscalar.cpp, xscalar.h, xvector.cpp, xvector.h)
3.1.224 - 2019/05/20
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.224.tar.xz
        -- Fixed bug that didn't allow users to override the standard lattice phonon dispersion path. (ME)
        -- Replaced anharmonic IFC option struct with xoption. (ME)
        -- Simplified APL class structure: IReciprocalPointGrid, MonkhorstPackMesh, and UniformMesh are consolidated into QMesh. IDOSCalculator, DOSCalculator, RootSamplingMethod, and LinearTetrahedronMethod are now described only by the DOSCalculator and LTMethod classes. (ME)
        -- Fixed QHA file name bugs for phonon dispersion and DOS. (ME)
        -- Fixed bug that caused AAPL to read APL calculations even though AAPL calculations are still missing. (ME)
        -- Fixed seg fault bug when DC=OFF and TP=ON. (ME)
        -- Fixed bug where DOS parameters are not properly set when DOS=OFF and TP=ON. (ME)
        -- Added KPPRA option for AAPL. (ME)
        -- Increased number of iterations for the AAPL BTE solver. (ME)
        -- fixed symmetry bug: if AddAtom() removes too many atoms and alters stoichiometry, then change to better tolerance (DX) 
        -- fixed comparison bug: ordering of stoichiometry for material-type comparisons (DX)
        -- changed name of prototyping flag (DX)
        -- added --structure2json command line functionality (DX)
        -- cleaned/reduced comparison command line calls (DX)
        -- added to comparison code: JSON output, print to screen only, keep non-matching info, store comparison logs, store family properties, and speed up comparison grouping function (DX)
        -- removed atomic library for communicating between threads, obsolete (DX)
        -- removed unused line in SYM::PBC() function (DX)
        -- added types information in CIF reader, otherwise aflowSG fails (DX)
        -- alphabetize elements in CIF reader function (DX)
        -- fixed --machine names (DX)
        -- fixed removal of control code characters in aflow.in; only overwrites if invalid chars are detected and saves old version to aflow.in_old (DX)
        -- changed internal indexing of xvectors in CHULL/GFA code from 0 to 1 to match rest of code (CO)
        -- fixed bug in web-json writer for AFLOW-CHULL Online (CO)
        -- added function for defining 3D rotation matrix between two vectors (CO)
        -- added functions to convert xvector<double> to xvector<int> and vice versa, useful for hkl (CO)
        -- added functionality for generalized stacking fault calculation, alpha testing only (CO)
        -- added functions for HKLPlane -> direct normal vector and vice versa (CO)
        -- added function to find distance between HKL planes (CO)
        -- added function to get distance between images of atoms along particular direction (CO)
        -- for many sym functions, removed pointers for read-only bools and doubles, added const where possible (CO)
        -- updated xstructure::Rotate() (CO)
        -- in xclasses, stringstream.clear() -> stringstream.str("") (CO)
        -- removed '#' from sg specification in POSCAR, it is interpreted as comment and removed (CO)
        -- bug fix with AFLOWDATA path in automatic aflow.in generator (CO)
        -- added force_strict_pc2scMap to GetSuperCell() constraining the map to the true primitive cell and not equivalent atoms (CO)
        -- added general function to fix string for latex output (CO)
        -- added XHOST.vflag_control.flag("DIRECTORY_CLEAN") for accurate logging of current directory (CO)
        -- major restructuring of pocc code, consolidating classes, renaming variables, etc. (CO)
        -- fixed old pocc to run with any compression extension (CO)
        -- eliminating atom.print_RHT (CO)
        -- patches for xStream() (CO)
        -- fixed webapp_entry for displaying first bader charge isosurface online with help from Bob Hanson (CO)
        -- added convolution function for signal processing (CO)
        -- added moving average and box/gaussian filters (CO)
        -- added mean/stddev functions to xvector (CO)
        -- added plotter for XRD (CO)
        -- added analysis for peaks with moving average (CO)
        (AFLOW: aconvasp_main.cpp, aflow.cpp, aflow.h, aflowlib_libraries.cpp, aflowlib_webapp_entry.js, aflowrc.cpp, avasp.cpp, chull.cpp, chull.h, compare_structure.cpp, compare_structure_function.cpp, compare_structure.h, gfa.cpp, init.cpp, kbin.cpp, kvasp.cpp, lattice.cpp, Makefile, modules.cpp, pflow_funcs.cpp, pflow.h, pflow_print.cpp, pocc.cpp, pocc.h, pocc_old.cpp, README_AFLOW_ACONVASP.TXT, README_AFLOW_AFLOWRC.TXT, symmetry.cpp, symmetry_spacegroup.cpp, xatom.cpp, xclasses.cpp, xproto.cpp)
        (ANRL: aflow_anrl.cpp, ANRL_CPPS_20180710/A10B2C_hP39_171_5c_c_a.cpp, ANRL_CPPS_20180710/A10B2C_hP39_172_5c_c_a.cpp, ANRL_CPPS_20180710/A10B3_oF52_42_2abce_ab.cpp, ANRL_CPPS_20180710/A12B2C_cF60_196_h_bc_a.cpp, ANRL_CPPS_20180710/A12B36CD12_cF488_196_2h_6h_ac_fgh.cpp, ANRL_CPPS_20180710/A12B36CD12_cF488_210_h_3h_a_fg.cpp, ANRL_CPPS_20180710/A12B6C_cF608_210_4h_2h_e.cpp, ANRL_CPPS_20180710/A12B7C2_hP21_174_2j2k_ajk_cf.cpp, ANRL_CPPS_20180710/A12BC4_cP34_195_2j_ab_2e.cpp, ANRL_CPPS_20180710/A12B_cF52_225_i_a.cpp, ANRL_CPPS_20180710/A12B_cI26_204_g_a.cpp, ANRL_CPPS_20180710/A12B_tI26_139_fij_a.cpp, ANRL_CPPS_20180710/A13B2C2_oP34_32_a6c_c_c.cpp, ANRL_CPPS_20180710/A14B3C5_tP44_94_c3g_ad_bg.cpp, ANRL_CPPS_20180710/A15B4_cI76_220_ae_c.cpp, ANRL_CPPS_20180710/A17B15_cP64_207_acfk_eij.cpp, ANRL_CPPS_20180710/A17BC4D_tP184_89_17p_p_4p_io.cpp, ANRL_CPPS_20180710/A23B6_cF116_225_bd2f_e.cpp, ANRL_CPPS_20180710/A24BC_cF104_209_j_a_b.cpp, ANRL_CPPS_20180710/A2B11_cP39_200_f_aghij.cpp, ANRL_CPPS_20180710/A2B13C4_hP57_168_d_c6d_2d.cpp, ANRL_CPPS_20180710/A2B2C4D_tP18_132_e_i_o_d.cpp, ANRL_CPPS_20180710/A2B2C7_cF88_227_c_d_af.cpp, ANRL_CPPS_20180710/A2B2C_oC80_64_efg_efg_df.cpp, ANRL_CPPS_20180710/A2B3C12D3_cI160_230_a_c_h_d.cpp, ANRL_CPPS_20180710/A2B3C18D6_hP58_192_c_f_lm_l.cpp, ANRL_CPPS_20180710/A2B3C3DE7_hP48_145_2a_3a_3a_a_7a.cpp, ANRL_CPPS_20180710/A2B3C6_cP264_205_2d_ab2c2d_6d.cpp, ANRL_CPPS_20180710/A2B3C6_cP33_221_cd_ag_fh.cpp, ANRL_CPPS_20180710/A2B3C7D_oP13_47_t_aq_eqrs_h.cpp, ANRL_CPPS_20180710/A2B3_hP20_159_bc_2c.cpp, ANRL_CPPS_20180710/A2B3_hP30_169_2a_3a.cpp, ANRL_CPPS_20180710/A2B3_hP30_170_2a_3a.cpp, ANRL_CPPS_20180710/A2B3_hP5_164_d_ad.cpp, ANRL_CPPS_20180710/A2B3_hR10_167_c_e.cpp, ANRL_CPPS_20180710/A2B3_hR5_166_c_ac.cpp, ANRL_CPPS_20180710/A2B3_oC40_39_2d_2c2d.cpp, ANRL_CPPS_20180710/A2B3_oP20_60_d_cd.cpp, ANRL_CPPS_20180710/A2B3_oP20_62_2c_3c.cpp, ANRL_CPPS_20180710/A2B3_oP40_33_4a_6a.cpp, ANRL_CPPS_20180710/A2B3_tI80_141_ceh_3h.cpp, ANRL_CPPS_20180710/A2B3_tP10_127_g_ah.cpp, ANRL_CPPS_20180710/A2B3_tP20_102_2c_b2c.cpp, ANRL_CPPS_20180710/A2B3_tP20_116_bci_fj.cpp, ANRL_CPPS_20180710/A2B3_tP20_117_i_adgh.cpp, ANRL_CPPS_20180710/A2B3_tP40_137_cdf_3g.cpp, ANRL_CPPS_20180710/A2B3_tP5_115_g_ag.cpp, ANRL_CPPS_20180710/A2B4C_hR42_148_2f_4f_f.cpp, ANRL_CPPS_20180710/A2B4C_oC28_66_l_kl_a.cpp, ANRL_CPPS_20180710/A2B4C_oP28_62_ac_2cd_c.cpp, ANRL_CPPS_20180710/A2B5C2_oC36_37_d_c2d_d.cpp, ANRL_CPPS_20180710/A2B5_mC28_15_f_e2f.cpp, ANRL_CPPS_20180710/A2B6CD7_tP64_77_2d_6d_d_ab6d.cpp, ANRL_CPPS_20180710/A2B7C2_tP88_78_4a_14a_4a.cpp, ANRL_CPPS_20180710/A2B7_cI54_229_e_afh.cpp, ANRL_CPPS_20180710/A2B8C2D_tP26_100_c_abcd_c_a.cpp, ANRL_CPPS_20180710/A2B8CD_oP24_48_k_2m_d_b.cpp, ANRL_CPPS_20180710/A2B8CD_tI24_97_d_k_a_b.cpp, ANRL_CPPS_20180710/A2B8C_oP22_34_c_4c_a.cpp, ANRL_CPPS_20180710/A2B_aP6_2_2i_i.cpp, ANRL_CPPS_20180710/A2B_aP6_2_aei_i.cpp, ANRL_CPPS_20180710/A2BC2_oF40_22_fi_ad_gh.cpp, ANRL_CPPS_20180710/A2BC2_oI20_45_c_b_c.cpp, ANRL_CPPS_20180710/A2BC2_tI20_79_c_2a_c.cpp, ANRL_CPPS_20180710/A2BC2_tP20_105_f_ac_2e.cpp, ANRL_CPPS_20180710/A2BC3_oC24_63_e_c_cg.cpp, ANRL_CPPS_20180710/A2BC4_cF56_227_d_a_e.cpp, ANRL_CPPS_20180710/A2BC4D_tI16_121_d_a_i_b.cpp, ANRL_CPPS_20180710/A2BC4_oP28_50_ij_ac_ijm.cpp, ANRL_CPPS_20180710/A2BC4_tI14_82_bc_a_g.cpp, ANRL_CPPS_20180710/A2BC4_tP28_126_cd_e_k.cpp, ANRL_CPPS_20180710/A2BC4_tP28_130_f_c_g.cpp, ANRL_CPPS_20180710/A2BC7D2_tP24_113_e_a_cef_e.cpp, ANRL_CPPS_20180710/A2BC8_tI176_110_2b_b_8b.cpp, ANRL_CPPS_20180710/A2BCD3E6_cF208_203_e_c_d_f_g.cpp, ANRL_CPPS_20180710/A2BCD4_tI16_82_ac_b_d_g.cpp, ANRL_CPPS_20180710/A2B_cF24_227_c_a.cpp, ANRL_CPPS_20180710/A2B_cF24_227_d_a.cpp, ANRL_CPPS_20180710/A2B_cI72_211_hi_i.cpp, ANRL_CPPS_20180710/A2BC_oC16_67_ag_b_g.cpp, ANRL_CPPS_20180710/A2BC_oC8_38_e_a_b.cpp, ANRL_CPPS_20180710/A2B_cP12_212_c_a.cpp, ANRL_CPPS_20180710/A2B_cP6_224_b_a.cpp, ANRL_CPPS_20180710/A2BC_tP16_76_2a_a_a.cpp, ANRL_CPPS_20180710/A2B_hP12_194_cg_f.cpp, ANRL_CPPS_20180710/A2B_hP18_180_fi_bd.cpp, ANRL_CPPS_20180710/A2B_hP18_190_gh_bf.cpp, ANRL_CPPS_20180710/A2B_hP36_177_j2lm_n.cpp, ANRL_CPPS_20180710/A2B_hP6_191_h_e.cpp, ANRL_CPPS_20180710/A2B_hP9_147_g_ad.cpp, ANRL_CPPS_20180710/A2B_hP9_150_ef_bd.cpp, ANRL_CPPS_20180710/A2B_hP9_152_c_a.cpp, ANRL_CPPS_20180710/A2B_hP9_180_j_c.cpp, ANRL_CPPS_20180710/A2B_hP9_181_j_c.cpp, ANRL_CPPS_20180710/A2B_hP9_189_fg_bc.cpp, ANRL_CPPS_20180710/A2B_hR18_148_2f_f.cpp, ANRL_CPPS_20180710/A2B_mC144_9_24a_12a.cpp, ANRL_CPPS_20180710/A2B_mC48_15_ae3f_2f.cpp, ANRL_CPPS_20180710/A2B_mP12_13_2g_ef.cpp, ANRL_CPPS_20180710/A2B_mP12_14_2e_e.cpp, ANRL_CPPS_20180710/A2B_mP12_3_bc3e_2e.cpp, ANRL_CPPS_20180710/A2B_mP12_7_4a_2a.cpp, ANRL_CPPS_20180710/A2B_mP18_7_6a_3a.cpp, ANRL_CPPS_20180710/A2B_mP6_10_mn_bg.cpp, ANRL_CPPS_20180710/A2B_mP6_14_e_a.cpp, ANRL_CPPS_20180710/A2B_oC12_36_2a_a.cpp, ANRL_CPPS_20180710/A2B_oC12_38_de_ab.cpp, ANRL_CPPS_20180710/A2B_oC12_63_2c_c.cpp, ANRL_CPPS_20180710/A2B_oC24_20_abc_c.cpp, ANRL_CPPS_20180710/A2B_oC24_64_2f_f.cpp, ANRL_CPPS_20180710/A2B_oF24_70_e_a.cpp, ANRL_CPPS_20180710/A2B_oI12_72_j_a.cpp, ANRL_CPPS_20180710/A2B_oI12_74_h_e.cpp, ANRL_CPPS_20180710/A2B_oP12_17_abe_e.cpp, ANRL_CPPS_20180710/A2B_oP12_19_2a_a.cpp, ANRL_CPPS_20180710/A2B_oP12_26_abc_ab.cpp, ANRL_CPPS_20180710/A2B_oP12_29_2a_a.cpp, ANRL_CPPS_20180710/A2B_oP12_62_2c_c.cpp, ANRL_CPPS_20180710/A2B_oP24_52_2e_cd.cpp, ANRL_CPPS_20180710/A2B_oP24_55_2g2h_gh.cpp, ANRL_CPPS_20180710/A2B_oP24_61_2c_c.cpp, ANRL_CPPS_20180710/A2B_oP6_58_g_a.cpp, ANRL_CPPS_20180710/A2B_tI12_140_h_a.cpp, ANRL_CPPS_20180710/A2B_tI12_141_e_a.cpp, ANRL_CPPS_20180710/A2B_tI12_98_f_a.cpp, ANRL_CPPS_20180710/A2B_tI24_141_2e_e.cpp, ANRL_CPPS_20180710/A2B_tI6_139_d_a.cpp, ANRL_CPPS_20180710/A2B_tP12_111_2n_adf.cpp, ANRL_CPPS_20180710/A2B_tP12_92_b_a.cpp, ANRL_CPPS_20180710/A2B_tP24_135_gh_h.cpp, ANRL_CPPS_20180710/A2B_tP30_85_ab2g_cg.cpp, ANRL_CPPS_20180710/A2B_tP36_96_3b_ab.cpp, ANRL_CPPS_20180710/A2B_tP48_77_8d_4d.cpp, ANRL_CPPS_20180710/A2B_tP6_129_ac_c.cpp, ANRL_CPPS_20180710/A2B_tP6_136_f_a.cpp, ANRL_CPPS_20180710/A2B_tP6_137_d_a.cpp, ANRL_CPPS_20180710/A3B10_cI52_229_e_fh.cpp, ANRL_CPPS_20180710/A3B11C6_tP40_100_ac_bc2d_cd.cpp, ANRL_CPPS_20180710/A3B13_oC32_38_ac_a2bcdef.cpp, ANRL_CPPS_20180710/A3B2_cI40_220_d_c.cpp, ANRL_CPPS_20180710/A3B2_hP10_176_h_bc.cpp, ANRL_CPPS_20180710/A3B2_hP10_176_h_bd.cpp, ANRL_CPPS_20180710/A3B2_hP5_164_ad_d.cpp, ANRL_CPPS_20180710/A3B2_hR5_155_e_c.cpp, ANRL_CPPS_20180710/A3B2_oP20_52_de_cd.cpp, ANRL_CPPS_20180710/A3B2_oP20_56_ce_e.cpp, ANRL_CPPS_20180710/A3B2_oP20_62_3c_2c.cpp, ANRL_CPPS_20180710/A3B2_tP10_83_adk_j.cpp, ANRL_CPPS_20180710/A3B3C_cI56_214_g_h_a.cpp, ANRL_CPPS_20180710/A3B3C_hP14_176_h_h_c.cpp, ANRL_CPPS_20180710/A3B3C_hP14_176_h_h_d.cpp, ANRL_CPPS_20180710/A3B4C_cP16_218_c_e_a.cpp, ANRL_CPPS_20180710/A3B4C_cP8_215_d_e_a.cpp, ANRL_CPPS_20180710/A3B4_cF56_227_ad_e.cpp, ANRL_CPPS_20180710/A3B4_tI28_141_ad_h.cpp, ANRL_CPPS_20180710/A3B5_oC16_65_ah_bej.cpp, ANRL_CPPS_20180710/A3B5_oC32_38_abce_abcdf.cpp, ANRL_CPPS_20180710/A3B5_oP16_55_ch_agh.cpp, ANRL_CPPS_20180710/A3B5_tI32_108_ac_a2c.cpp, ANRL_CPPS_20180710/A3B5_tI32_140_ah_bk.cpp, ANRL_CPPS_20180710/A3B5_tI32_140_ah_cl.cpp, ANRL_CPPS_20180710/A3B7_hP20_186_c_b2c.cpp, ANRL_CPPS_20180710/A3B7_oP40_62_cd_3c2d.cpp, ANRL_CPPS_20180710/A3B7_tP40_76_3a_7a.cpp, ANRL_CPPS_20180710/A3BC2_cI48_214_f_a_e.cpp, ANRL_CPPS_20180710/A3BC2_oP48_50_3m_m_2m.cpp, ANRL_CPPS_20180710/A3BC3D_tP64_106_3c_c_3c_c.cpp, ANRL_CPPS_20180710/A3BC_hP10_188_k_a_e.cpp, ANRL_CPPS_20180710/A3BC_hP10_188_k_c_a.cpp, ANRL_CPPS_20180710/A3BC_hP30_185_cd_c_ab.cpp, ANRL_CPPS_20180710/A3BC_hR5_146_b_a_a.cpp, ANRL_CPPS_20180710/A3B_cI32_204_g_c.cpp, ANRL_CPPS_20180710/A3B_cI8_229_b_a.cpp, ANRL_CPPS_20180710/A3BC_mC10_8_ab_a_a.cpp, ANRL_CPPS_20180710/A3B_cP16_198_b_a.cpp, ANRL_CPPS_20180710/A3B_cP16_208_j_b.cpp, ANRL_CPPS_20180710/A3B_cP4_221_d_a.cpp, ANRL_CPPS_20180710/A3B_cP8_223_c_a.cpp, ANRL_CPPS_20180710/A3BC_tP5_99_bc_a_b.cpp, ANRL_CPPS_20180710/A3B_hP16_194_gh_ac.cpp, ANRL_CPPS_20180710/A3B_hP24_151_3c_2a.cpp, ANRL_CPPS_20180710/A3B_hP24_153_3c_2b.cpp, ANRL_CPPS_20180710/A3B_hP24_165_adg_f.cpp, ANRL_CPPS_20180710/A3B_hP24_165_bdg_f.cpp, ANRL_CPPS_20180710/A3B_hP24_185_ab2c_c.cpp, ANRL_CPPS_20180710/A3B_hP4_191_bc_a.cpp, ANRL_CPPS_20180710/A3B_hP8_158_d_a.cpp, ANRL_CPPS_20180710/A3B_hP8_173_c_b.cpp, ANRL_CPPS_20180710/A3B_hP8_176_h_c.cpp, ANRL_CPPS_20180710/A3B_hP8_176_h_d.cpp, ANRL_CPPS_20180710/A3B_hP8_185_c_a.cpp, ANRL_CPPS_20180710/A3B_hP8_194_h_c.cpp, ANRL_CPPS_20180710/A3B_hR4_160_b_a.cpp, ANRL_CPPS_20180710/A3B_mP16_7_6a_2a.cpp, ANRL_CPPS_20180710/A3B_oC64_66_gi2lm_2l.cpp, ANRL_CPPS_20180710/A3B_oC64_66_kl2m_bdl.cpp, ANRL_CPPS_20180710/A3B_oI32_23_ij2k_k.cpp, ANRL_CPPS_20180710/A3B_oP16_62_cd_c.cpp, ANRL_CPPS_20180710/A3B_oP32_60_3d_d.cpp, ANRL_CPPS_20180710/A3B_oP8_59_bf_a.cpp, ANRL_CPPS_20180710/A3B_tI16_139_cde_e.cpp, ANRL_CPPS_20180710/A3B_tI24_119_b2i_af.cpp, ANRL_CPPS_20180710/A3B_tI32_82_3g_g.cpp, ANRL_CPPS_20180710/A3B_tI8_139_bd_a.cpp, ANRL_CPPS_20180710/A3B_tP16_118_ei_f.cpp, ANRL_CPPS_20180710/A3B_tP32_114_3e_e.cpp, ANRL_CPPS_20180710/A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h.cpp, ANRL_CPPS_20180710/A4B14C3_hP21_143_bd_ac4d_d.cpp, ANRL_CPPS_20180710/A4B2C13D_tP40_90_g_d_cef2g_c.cpp, ANRL_CPPS_20180710/A4B2C6D16E_cF232_203_e_d_f_eg_a.cpp, ANRL_CPPS_20180710/A4B2C_tP28_135_gh_h_d.cpp, ANRL_CPPS_20180710/A4B3_cI112_230_af_g.cpp, ANRL_CPPS_20180710/A4B3_cI14_229_c_b.cpp, ANRL_CPPS_20180710/A4B3_cI28_220_c_a.cpp, ANRL_CPPS_20180710/A4B3_hP14_173_bc_c.cpp, ANRL_CPPS_20180710/A4B3_hP28_159_ab2c_2c.cpp, ANRL_CPPS_20180710/A4B3_hR7_166_2c_ac.cpp, ANRL_CPPS_20180710/A4B3_oI14_71_gh_cg.cpp, ANRL_CPPS_20180710/A4B5_tI18_139_i_ah.cpp, ANRL_CPPS_20180710/A4B5_tI18_87_h_ah.cpp, ANRL_CPPS_20180710/A4B6C_hP11_143_bd_2d_a.cpp, ANRL_CPPS_20180710/A4B9_cP52_215_ei_3efgi.cpp, ANRL_CPPS_20180710/A4BC2_tI28_120_i_d_e.cpp, ANRL_CPPS_20180710/A4BC4D_tP10_123_gh_a_i_d.cpp, ANRL_CPPS_20180710/A4BC4_tP18_137_g_b_g.cpp, ANRL_CPPS_20180710/A4B_cI10_217_c_a.cpp, ANRL_CPPS_20180710/A4B_cI10_229_c_a.cpp, ANRL_CPPS_20180710/A4B_cI40_197_cde_c.cpp, ANRL_CPPS_20180710/A4BC_tI24_141_h_b_a.cpp, ANRL_CPPS_20180710/A4B_hP15_144_4a_a.cpp, ANRL_CPPS_20180710/A4B_oI20_74_beh_e.cpp, ANRL_CPPS_20180710/A4B_oP20_62_2cd_c.cpp, ANRL_CPPS_20180710/A4B_tI10_139_de_a.cpp, ANRL_CPPS_20180710/A4B_tI20_88_f_a.cpp, ANRL_CPPS_20180710/A4B_tP10_114_e_a.cpp, ANRL_CPPS_20180710/A4B_tP10_125_m_a.cpp, ANRL_CPPS_20180710/A4B_tP20_127_ehj_g.cpp, ANRL_CPPS_20180710/A5B11_mP16_6_2abc_2a3b3c.cpp, ANRL_CPPS_20180710/A5B2_hP14_194_abdf_f.cpp, ANRL_CPPS_20180710/A5B2_hP28_194_ahk_ch.cpp, ANRL_CPPS_20180710/A5B2_hR7_166_a2c_c.cpp, ANRL_CPPS_20180710/A5B2_mC14_12_a2i_i.cpp, ANRL_CPPS_20180710/A5B2_oP14_49_dehq_ab.cpp, ANRL_CPPS_20180710/A5B3C15_oP46_30_a2c_bc_a7c.cpp, ANRL_CPPS_20180710/A5B3C16_cP96_222_ce_d_fi.cpp, ANRL_CPPS_20180710/A5B3C_hP18_186_2a3b_2ab_b.cpp, ANRL_CPPS_20180710/A5B3_hP16_190_bdh_g.cpp, ANRL_CPPS_20180710/A5B3_hP16_193_dg_g.cpp, ANRL_CPPS_20180710/A5B3_mC32_9_5a_3a.cpp, ANRL_CPPS_20180710/A5B3_tP32_118_g2i_aceh.cpp, ANRL_CPPS_20180710/A5B3_tP32_130_cg_cf.cpp, ANRL_CPPS_20180710/A5B5C4_tP28_104_ac_ac_c.cpp, ANRL_CPPS_20180710/A5B6C2_hP13_157_2ac_2c_b.cpp, ANRL_CPPS_20180710/A5B7_tI24_107_ac_abd.cpp, ANRL_CPPS_20180710/A5B8_cI52_217_ce_cg.cpp, ANRL_CPPS_20180710/A5BCD6_cF416_228_eg_c_b_h.cpp, ANRL_CPPS_20180710/A5B_oP24_26_3a3b2c_ab.cpp, ANRL_CPPS_20180710/A6B23_cF116_225_e_acfh.cpp, ANRL_CPPS_20180710/A6B2C_cF36_225_e_c_a.cpp, ANRL_CPPS_20180710/A6B2CD6E_cP64_208_m_ad_b_m_c.cpp, ANRL_CPPS_20180710/A6B2C_tP18_128_eh_d_a.cpp, ANRL_CPPS_20180710/A6B2C_tP18_128_eh_d_b.cpp, ANRL_CPPS_20180710/A6B2C_tP18_94_eg_c_a.cpp, ANRL_CPPS_20180710/A6B4C16D_oP108_27_abcd4e_4e_16e_e.cpp, ANRL_CPPS_20180710/A6B6C_cF104_202_h_h_c.cpp, ANRL_CPPS_20180710/A6B_cF224_228_h_c.cpp, ANRL_CPPS_20180710/A6B_cP7_221_f_a.cpp, ANRL_CPPS_20180710/A6B_hR7_166_g_a.cpp, ANRL_CPPS_20180710/A6B_oC28_63_efg_c.cpp, ANRL_CPPS_20180710/A7B2C2_mC22_12_aij_h_i.cpp, ANRL_CPPS_20180710/A7B2C_tP40_128_egi_h_e.cpp, ANRL_CPPS_20180710/A7B3_cI40_229_df_e.cpp, ANRL_CPPS_20180710/A7B6_hR13_166_ah_3c.cpp, ANRL_CPPS_20180710/A7B7C2_tP32_101_bde_ade_d.cpp, ANRL_CPPS_20180710/A7B8_mP120_14_14e_16e.cpp, ANRL_CPPS_20180710/A7B8_oP120_60_7d_8d.cpp, ANRL_CPPS_20180710/A7BC3D13_cF192_219_de_b_c_ah.cpp, ANRL_CPPS_20180710/A7B_cF32_225_bd_a.cpp, ANRL_CPPS_20180710/A8B2C12D2E_oI50_23_bcfk_i_3k_j_a.cpp, ANRL_CPPS_20180710/A8B5_hR26_160_a3bc_a3b.cpp, ANRL_CPPS_20180710/A8B5_mP13_6_a7b_3a2b.cpp, ANRL_CPPS_20180710/A8B7C6_hP21_175_ck_aj_k.cpp, ANRL_CPPS_20180710/A8BC3D6_hP18_189_bfh_a_g_i.cpp, ANRL_CPPS_20180710/A8B_tI18_139_hi_a.cpp, ANRL_CPPS_20180710/A9B16C7_cF128_225_acd_2f_be.cpp, ANRL_CPPS_20180710/A9B2_mP22_7_9a_2a.cpp, ANRL_CPPS_20180710/A9B3C_hP26_194_hk_h_a.cpp, ANRL_CPPS_20180710/A9BC3D5_hP18_189_fi_a_g_bh.cpp, ANRL_CPPS_20180710/A9BC_oC44_39_3c3d_a_c.cpp, ANRL_CPPS_20180710/A_aP4_2_aci.cpp, ANRL_CPPS_20180710/AB11CD3_cP16_221_a_dg_b_c.cpp, ANRL_CPPS_20180710/AB11_cP36_221_c_agij.cpp, ANRL_CPPS_20180710/AB12C3_cI32_229_a_h_b.cpp, ANRL_CPPS_20180710/AB13_cF112_226_a_bi.cpp, ANRL_CPPS_20180710/AB18C8_cF108_225_a_eh_f.cpp, ANRL_CPPS_20180710/AB27CD3_cP32_221_a_dij_b_c.cpp, ANRL_CPPS_20180710/AB2_aP12_1_4a_8a.cpp, ANRL_CPPS_20180710/AB2C12D4_tP76_75_2a2b_2d_12d_4d.cpp, ANRL_CPPS_20180710/AB2C3_oP24_62_c_d_cd.cpp, ANRL_CPPS_20180710/AB2C4_tI14_139_a_e_ce.cpp, ANRL_CPPS_20180710/AB2C8D_oP24_49_g_q_2qr_e.cpp, ANRL_CPPS_20180710/AB2C_cF16_225_a_c_b.cpp, ANRL_CPPS_20180710/AB2CD2_hP36_163_h_i_bf_i.cpp, ANRL_CPPS_20180710/AB2_cF12_225_a_c.cpp, ANRL_CPPS_20180710/AB2_cF48_227_c_e.cpp, ANRL_CPPS_20180710/AB2_cF96_227_e_cf.cpp, ANRL_CPPS_20180710/AB2C_oC16_40_a_2b_b.cpp, ANRL_CPPS_20180710/AB2C_oC16_63_c_2c_c.cpp, ANRL_CPPS_20180710/AB2C_oP16_62_c_2c_c.cpp, ANRL_CPPS_20180710/AB2_cP12_205_a_c.cpp, ANRL_CPPS_20180710/AB2_hP12_143_cd_ab2d.cpp, ANRL_CPPS_20180710/AB2_hP12_194_f_ah.cpp, ANRL_CPPS_20180710/AB2_hP24_194_ef_fgh.cpp, ANRL_CPPS_20180710/AB2_hP3_164_a_d.cpp, ANRL_CPPS_20180710/AB2_hP3_191_a_d.cpp, ANRL_CPPS_20180710/AB2_hP6_194_b_f.cpp, ANRL_CPPS_20180710/AB2_hP6_194_c_ad.cpp, ANRL_CPPS_20180710/AB2_hP6_194_c_f.cpp, ANRL_CPPS_20180710/AB2_hP72_192_m_j2kl.cpp, ANRL_CPPS_20180710/AB2_hP9_156_b2c_3a2bc.cpp, ANRL_CPPS_20180710/AB2_hP9_162_ad_k.cpp, ANRL_CPPS_20180710/AB2_hP9_164_bd_c2d.cpp, ANRL_CPPS_20180710/AB2_hP9_180_d_j.cpp, ANRL_CPPS_20180710/AB2_mC6_12_a_i.cpp, ANRL_CPPS_20180710/AB2_oC24_41_2a_2b.cpp, ANRL_CPPS_20180710/AB2_oC6_21_a_k.cpp, ANRL_CPPS_20180710/AB2_oF48_70_f_fg.cpp, ANRL_CPPS_20180710/AB2_oF72_43_ab_3b.cpp, ANRL_CPPS_20180710/AB2_oI6_71_a_g.cpp, ANRL_CPPS_20180710/AB2_oI6_71_a_i.cpp, ANRL_CPPS_20180710/AB2_oP12_29_a_2a.cpp, ANRL_CPPS_20180710/AB2_oP24_28_acd_2c3d.cpp, ANRL_CPPS_20180710/AB2_oP6_34_a_c.cpp, ANRL_CPPS_20180710/AB2_oP6_58_a_g.cpp, ANRL_CPPS_20180710/AB2_tI48_80_2b_4b.cpp, ANRL_CPPS_20180710/AB2_tI6_139_a_e.cpp, ANRL_CPPS_20180710/AB2_tI96_88_2f_4f.cpp, ANRL_CPPS_20180710/AB2_tP12_115_j_egi.cpp, ANRL_CPPS_20180710/AB2_tP12_81_adg_2h.cpp, ANRL_CPPS_20180710/AB2_tP6_137_a_d.cpp, ANRL_CPPS_20180710/AB32C48_cI162_204_a_2efg_2gh.cpp, ANRL_CPPS_20180710/AB32CD4E8_tP184_93_i_16p_af_2p_4p.cpp, ANRL_CPPS_20180710/AB3C16_cF160_203_a_bc_eg.cpp, ANRL_CPPS_20180710/AB3C16_cF160_203_b_ad_eg.cpp, ANRL_CPPS_20180710/AB3C2_cI96_206_c_e_ad.cpp, ANRL_CPPS_20180710/AB3C3_cF112_227_c_de_f.cpp, ANRL_CPPS_20180710/AB3C4_cP8_215_a_c_e.cpp, ANRL_CPPS_20180710/AB3C4_hP16_194_c_af_ef.cpp, ANRL_CPPS_20180710/AB3C4_oP16_31_a_ab_2ab.cpp, ANRL_CPPS_20180710/AB3C4_oP32_33_a_3a_4a.cpp, ANRL_CPPS_20180710/AB3C6_cI80_206_a_d_e.cpp, ANRL_CPPS_20180710/AB3C_cP5_221_a_c_b.cpp, ANRL_CPPS_20180710/AB3C_cP60_201_be_fh_g.cpp, ANRL_CPPS_20180710/AB3C_cP60_201_ce_fh_g.cpp, ANRL_CPPS_20180710/AB3_cF16_225_a_bc.cpp, ANRL_CPPS_20180710/AB3C_hR10_148_c_f_c.cpp, ANRL_CPPS_20180710/AB3C_hR10_167_b_e_a.cpp, ANRL_CPPS_20180710/AB3C_oC20_63_a_cf_c.cpp, ANRL_CPPS_20180710/AB3C_oP20_62_c_cd_a.cpp, ANRL_CPPS_20180710/AB3_cP4_221_a_c.cpp, ANRL_CPPS_20180710/AB3_hP24_149_acgi_3l.cpp, ANRL_CPPS_20180710/AB3_hP24_178_b_ac.cpp, ANRL_CPPS_20180710/AB3_hP24_179_b_ac.cpp, ANRL_CPPS_20180710/AB3_hP24_185_c_ab2c.cpp, ANRL_CPPS_20180710/AB3_hP4_187_e_fh.cpp, ANRL_CPPS_20180710/AB3_hP8_182_c_g.cpp, ANRL_CPPS_20180710/AB3_hP8_194_c_bf.cpp, ANRL_CPPS_20180710/AB3_hR8_148_c_f.cpp, ANRL_CPPS_20180710/AB3_hR8_155_c_de.cpp, ANRL_CPPS_20180710/AB3_mC16_12_g_ij.cpp, ANRL_CPPS_20180710/AB3_mC16_15_e_cf.cpp, ANRL_CPPS_20180710/AB3_mC16_9_a_3a.cpp, ANRL_CPPS_20180710/AB3_mC32_8_4a_12a.cpp, ANRL_CPPS_20180710/AB3_mC32_8_4a_4a4b.cpp, ANRL_CPPS_20180710/AB3_mP16_10_mn_3m3n.cpp, ANRL_CPPS_20180710/AB3_oC16_40_b_3b.cpp, ANRL_CPPS_20180710/AB3_oC8_65_a_bf.cpp, ANRL_CPPS_20180710/AB3_oP16_18_ab_3c.cpp, ANRL_CPPS_20180710/AB3_oP16_19_a_3a.cpp, ANRL_CPPS_20180710/AB3_oP16_62_c_3c.cpp, ANRL_CPPS_20180710/AB3_oP16_62_c_cd.cpp, ANRL_CPPS_20180710/AB3_tI16_140_b_ah.cpp, ANRL_CPPS_20180710/AB3_tP32_133_h_i2j.cpp, ANRL_CPPS_20180710/AB3_tP32_86_g_3g.cpp, ANRL_CPPS_20180710/AB3_tP4_123_a_ce.cpp, ANRL_CPPS_20180710/AB3_tP8_113_a_ce.cpp, ANRL_CPPS_20180710/AB4C17D4E_tP54_90_a_g_c4g_g_c.cpp, ANRL_CPPS_20180710/AB4C3_cI16_229_a_c_b.cpp, ANRL_CPPS_20180710/AB4C7D_hP26_159_b_ac_a2c_b.cpp, ANRL_CPPS_20180710/AB4C_hP6_191_a_h_b.cpp, ANRL_CPPS_20180710/AB4C_hP72_168_2d_8d_2d.cpp, ANRL_CPPS_20180710/AB4C_hP72_184_d_4d_d.cpp, ANRL_CPPS_20180710/AB4C_oC24_63_a_fg_c.cpp, ANRL_CPPS_20180710/AB4C_oC24_63_c_fg_c.cpp, ANRL_CPPS_20180710/AB4C_oP24_62_c_2cd_c.cpp, ANRL_CPPS_20180710/AB4_cP5_215_a_e.cpp, ANRL_CPPS_20180710/AB4C_tI12_82_c_g_a.cpp, ANRL_CPPS_20180710/AB4C_tP12_112_b_n_e.cpp, ANRL_CPPS_20180710/AB4C_tP12_124_a_m_c.cpp, ANRL_CPPS_20180710/AB4_oC20_41_a_2b.cpp, ANRL_CPPS_20180710/AB4_oC20_68_a_i.cpp, ANRL_CPPS_20180710/AB4_tI10_87_a_h.cpp, ANRL_CPPS_20180710/AB4_tP10_103_a_d.cpp, ANRL_CPPS_20180710/AB4_tP10_124_a_m.cpp, ANRL_CPPS_20180710/AB5_cF24_216_a_ce.cpp, ANRL_CPPS_20180710/AB5C_tP7_123_b_ci_a.cpp, ANRL_CPPS_20180710/AB5_hP6_191_a_cg.cpp, ANRL_CPPS_20180710/AB6C4_tP22_104_a_2ac_c.cpp, ANRL_CPPS_20180710/AB6C_tP16_132_d_io_a.cpp, ANRL_CPPS_20180710/AB7CD2_oI44_24_a_b3d_c_ac.cpp, ANRL_CPPS_20180710/AB7_hR16_166_c_c2h.cpp, ANRL_CPPS_20180710/AB8C2_oC22_35_a_ab3e_e.cpp, ANRL_CPPS_20180710/AB8C2_tI44_97_e_2k_cd.cpp, ANRL_CPPS_20180710/AB9C4_hP28_188_e_kl_ak.cpp, ANRL_CPPS_20180710/AB_aP16_2_4i_4i.cpp, ANRL_CPPS_20180710/ABC2_aP16_1_4a_4a_8a.cpp, ANRL_CPPS_20180710/ABC2_hP4_164_a_b_d.cpp, ANRL_CPPS_20180710/ABC2_hP8_194_d_a_f.cpp, ANRL_CPPS_20180710/ABC2_hR24_167_e_e_2e.cpp, ANRL_CPPS_20180710/ABC2_hR4_166_a_b_c.cpp, ANRL_CPPS_20180710/ABC2_mP8_10_ac_eh_mn.cpp, ANRL_CPPS_20180710/ABC2_oC16_67_b_g_ag.cpp, ANRL_CPPS_20180710/ABC2_oI16_23_ab_i_k.cpp, ANRL_CPPS_20180710/ABC2_oP16_53_h_e_gh.cpp, ANRL_CPPS_20180710/ABC2_tI16_122_a_b_d.cpp, ANRL_CPPS_20180710/ABC2_tP4_123_d_a_f.cpp, ANRL_CPPS_20180710/ABC3_cP20_198_a_a_b.cpp, ANRL_CPPS_20180710/ABC3_hR10_146_2a_2a_2b.cpp, ANRL_CPPS_20180710/ABC3_hR10_161_a_a_b.cpp, ANRL_CPPS_20180710/ABC3_hR10_167_a_b_e.cpp, ANRL_CPPS_20180710/ABC3_mP10_11_e_e_ef.cpp, ANRL_CPPS_20180710/ABC3_oP20_30_2a_c_3c.cpp, ANRL_CPPS_20180710/ABC3_oP20_53_e_g_hi.cpp, ANRL_CPPS_20180710/ABC3_oP20_54_e_d_cf.cpp, ANRL_CPPS_20180710/ABC4_mP12_13_e_a_2g.cpp, ANRL_CPPS_20180710/ABC4_oI12_23_a_b_k.cpp, ANRL_CPPS_20180710/ABC4_oP12_16_ag_cd_2u.cpp, ANRL_CPPS_20180710/ABC4_tI96_142_e_ab_2g.cpp, ANRL_CPPS_20180710/ABC4_tP12_125_a_b_m.cpp, ANRL_CPPS_20180710/ABC6D2_mC40_15_e_e_3f_f.cpp, ANRL_CPPS_20180710/ABC_cF12_216_b_c_a.cpp, ANRL_CPPS_20180710/ABC_cP12_198_a_a_a.cpp, ANRL_CPPS_20180710/ABCD3_oI48_73_d_e_e_ef.cpp, ANRL_CPPS_20180710/ABCD_cF16_216_c_d_b_a.cpp, ANRL_CPPS_20180710/ABCD_oP16_57_d_c_d_d.cpp, ANRL_CPPS_20180710/ABCD_tP8_129_c_b_a_c.cpp, ANRL_CPPS_20180710/AB_cF16_227_a_b.cpp, ANRL_CPPS_20180710/AB_cF8_216_c_a.cpp, ANRL_CPPS_20180710/AB_cF8_225_a_b.cpp, ANRL_CPPS_20180710/ABC_hP12_174_cj_fk_aj.cpp, ANRL_CPPS_20180710/ABC_hP3_183_a_a_a.cpp, ANRL_CPPS_20180710/ABC_hP3_187_a_d_f.cpp, ANRL_CPPS_20180710/ABC_hP36_175_jk_jk_jk.cpp, ANRL_CPPS_20180710/ABC_hP6_194_c_d_a.cpp, ANRL_CPPS_20180710/ABC_hR3_160_a_a_a.cpp, ANRL_CPPS_20180710/ABC_hR6_166_c_c_c.cpp, ANRL_CPPS_20180710/AB_cI16_199_a_a.cpp, ANRL_CPPS_20180710/ABC_oI12_71_h_j_g.cpp, ANRL_CPPS_20180710/ABC_oI36_46_ac_bc_3b.cpp, ANRL_CPPS_20180710/ABC_oP12_29_a_a_a.cpp, ANRL_CPPS_20180710/ABC_oP6_59_a_a_a.cpp, ANRL_CPPS_20180710/ABC_oP6_59_a_b_a.cpp, ANRL_CPPS_20180710/AB_cP16_205_c_c.cpp, ANRL_CPPS_20180710/AB_cP2_221_b_a.cpp, ANRL_CPPS_20180710/AB_cP6_221_c_d.cpp, ANRL_CPPS_20180710/AB_cP8_198_a_a.cpp, ANRL_CPPS_20180710/ABC_tI12_109_a_a_a.cpp, ANRL_CPPS_20180710/ABC_tP24_91_d_d_d.cpp, ANRL_CPPS_20180710/ABC_tP24_95_d_d_d.cpp, ANRL_CPPS_20180710/ABC_tP6_129_c_a_c.cpp, ANRL_CPPS_20180710/AB_hP12_156_2ab3c_2ab3c.cpp, ANRL_CPPS_20180710/AB_hP12_186_a2b_a2b.cpp, ANRL_CPPS_20180710/AB_hP12_194_af_bf.cpp, ANRL_CPPS_20180710/AB_hP12_194_df_ce.cpp, ANRL_CPPS_20180710/AB_hP2_187_d_a.cpp, ANRL_CPPS_20180710/AB_hP24_190_i_afh.cpp, ANRL_CPPS_20180710/AB_hP4_156_ab_ab.cpp, ANRL_CPPS_20180710/AB_hP4_156_ac_ac.cpp, ANRL_CPPS_20180710/AB_hP4_186_b_a.cpp, ANRL_CPPS_20180710/AB_hP4_186_b_b.cpp, ANRL_CPPS_20180710/AB_hP4_194_c_a.cpp, ANRL_CPPS_20180710/AB_hP4_194_c_d.cpp, ANRL_CPPS_20180710/AB_hP6_144_a_a.cpp, ANRL_CPPS_20180710/AB_hP6_154_a_b.cpp, ANRL_CPPS_20180710/AB_hP6_183_c_ab.cpp, ANRL_CPPS_20180710/AB_hP6_191_f_ad.cpp, ANRL_CPPS_20180710/AB_hP8_186_ab_ab.cpp, ANRL_CPPS_20180710/AB_hP8_194_ad_f.cpp, ANRL_CPPS_20180710/AB_hR10_160_5a_5a.cpp, ANRL_CPPS_20180710/AB_hR16_148_cf_cf.cpp, ANRL_CPPS_20180710/AB_hR2_166_a_b.cpp, ANRL_CPPS_20180710/AB_hR26_148_b2f_a2f.cpp, ANRL_CPPS_20180710/AB_hR6_160_3a_3a.cpp, ANRL_CPPS_20180710/AB_hR6_160_b_b.cpp, ANRL_CPPS_20180710/AB_mC8_15_c_e.cpp, ANRL_CPPS_20180710/AB_mP4_11_e_e.cpp, ANRL_CPPS_20180710/AB_mP4_6_2b_2a.cpp, ANRL_CPPS_20180710/AB_mP6_10_en_am.cpp, ANRL_CPPS_20180710/AB_oC8_36_a_a.cpp, ANRL_CPPS_20180710/AB_oC8_63_c_c.cpp, ANRL_CPPS_20180710/AB_oC8_65_j_g.cpp, ANRL_CPPS_20180710/AB_oC8_67_a_g.cpp, ANRL_CPPS_20180710/AB_oF8_22_a_c.cpp, ANRL_CPPS_20180710/AB_oF8_42_a_a.cpp, ANRL_CPPS_20180710/AB_oF8_69_a_b.cpp, ANRL_CPPS_20180710/AB_oI4_44_a_b.cpp, ANRL_CPPS_20180710/AB_oP16_61_c_c.cpp, ANRL_CPPS_20180710/AB_oP2_25_b_a.cpp, ANRL_CPPS_20180710/AB_oP4_51_e_f.cpp, ANRL_CPPS_20180710/AB_oP4_59_a_b.cpp, ANRL_CPPS_20180710/AB_oP48_61_3c_3c.cpp, ANRL_CPPS_20180710/AB_oP8_33_a_a.cpp, ANRL_CPPS_20180710/AB_oP8_57_d_d.cpp, ANRL_CPPS_20180710/AB_oP8_62_c_c.cpp, ANRL_CPPS_20180710/AB_tI16_140_ab_h.cpp, ANRL_CPPS_20180710/AB_tI16_141_e_e.cpp, ANRL_CPPS_20180710/AB_tI4_107_a_a.cpp, ANRL_CPPS_20180710/AB_tI4_119_c_a.cpp, ANRL_CPPS_20180710/AB_tI8_109_a_a.cpp, ANRL_CPPS_20180710/AB_tI8_139_e_e.cpp, ANRL_CPPS_20180710/AB_tI8_141_a_b.cpp, ANRL_CPPS_20180710/AB_tP16_84_cej_k.cpp, ANRL_CPPS_20180710/AB_tP2_123_a_d.cpp, ANRL_CPPS_20180710/AB_tP4_129_a_c.cpp, ANRL_CPPS_20180710/AB_tP4_129_c_c.cpp, ANRL_CPPS_20180710/AB_tP4_131_c_e.cpp, ANRL_CPPS_20180710/AB_tP8_111_n_n.cpp, ANRL_CPPS_20180710/AB_tP8_136_g_f.cpp, ANRL_CPPS_20180710/A_cF136_227_aeg.cpp, ANRL_CPPS_20180710/A_cF240_202_h2i.cpp, ANRL_CPPS_20180710/A_cF4_225_a.cpp, ANRL_CPPS_20180710/A_cF8_227_a.cpp, ANRL_CPPS_20180710/A_cI16_206_c.cpp, ANRL_CPPS_20180710/A_cI16_220_c.cpp, ANRL_CPPS_20180710/A_cI2_229_a.cpp, ANRL_CPPS_20180710/A_cI58_217_ac2g.cpp, ANRL_CPPS_20180710/A_cP1_221_a.cpp, ANRL_CPPS_20180710/A_cP20_213_cd.cpp, ANRL_CPPS_20180710/A_cP240_205_10d.cpp, ANRL_CPPS_20180710/A_cP46_223_dik.cpp, ANRL_CPPS_20180710/A_cP8_198_2a.cpp, ANRL_CPPS_20180710/A_cP8_205_c.cpp, ANRL_CPPS_20180710/A_hP1_191_a.cpp, ANRL_CPPS_20180710/A_hP2_194_c.cpp, ANRL_CPPS_20180710/A_hP3_152_a.cpp, ANRL_CPPS_20180710/A_hP4_186_ab.cpp, ANRL_CPPS_20180710/A_hP4_194_ac.cpp, ANRL_CPPS_20180710/A_hP4_194_bc.cpp, ANRL_CPPS_20180710/A_hP4_194_f.cpp, ANRL_CPPS_20180710/A_hP6_178_a.cpp, ANRL_CPPS_20180710/A_hP6_194_h.cpp, ANRL_CPPS_20180710/A_hP9_154_bc.cpp, ANRL_CPPS_20180710/A_hR105_166_bc9h4i.cpp, ANRL_CPPS_20180710/A_hR1_166_a.cpp, ANRL_CPPS_20180710/A_hR12_166_2h.cpp, ANRL_CPPS_20180710/A_hR2_166_c.cpp, ANRL_CPPS_20180710/A_hR3_166_ac.cpp, ANRL_CPPS_20180710/A_mC12_5_3c.cpp, ANRL_CPPS_20180710/A_mC16_12_4i.cpp, ANRL_CPPS_20180710/A_mC24_15_2e2f.cpp, ANRL_CPPS_20180710/A_mC34_12_ah3i2j.cpp, ANRL_CPPS_20180710/A_mC4_12_i.cpp, ANRL_CPPS_20180710/A_mP16_11_8e.cpp, ANRL_CPPS_20180710/A_mP32_14_8e.cpp, ANRL_CPPS_20180710/A_mP4_4_2a.cpp, ANRL_CPPS_20180710/A_mP64_14_16e.cpp, ANRL_CPPS_20180710/A_mP8_10_2m2n.cpp, ANRL_CPPS_20180710/A_mP84_13_21g.cpp, ANRL_CPPS_20180710/A_oC4_63_c.cpp, ANRL_CPPS_20180710/A_oC8_64_f.cpp, ANRL_CPPS_20180710/A_oF128_70_4h.cpp, ANRL_CPPS_20180710/A_oF8_70_a.cpp, ANRL_CPPS_20180710/A_oP16_55_2g2h.cpp, ANRL_CPPS_20180710/A_oP8_62_2c.cpp, ANRL_CPPS_20180710/A_tI16_142_f.cpp, ANRL_CPPS_20180710/A_tI2_139_a.cpp, ANRL_CPPS_20180710/A_tI4_139_e.cpp, ANRL_CPPS_20180710/A_tI4_141_a.cpp, ANRL_CPPS_20180710/A_tI8_139_h.cpp, ANRL_CPPS_20180710/A_tP12_138_bi.cpp, ANRL_CPPS_20180710/A_tP12_96_ab.cpp, ANRL_CPPS_20180710/A_tP16_138_j.cpp, ANRL_CPPS_20180710/A_tP30_136_bf2ij.cpp, ANRL_CPPS_20180710/A_tP4_129_ac.cpp, ANRL_CPPS_20180710/A_tP4_136_f.cpp, ANRL_CPPS_20180710/A_tP50_134_b2m2n.cpp, ANRL_CPPS_20180710/sigma_tP30_136_bf2ij.cpp)
        (APL: aapl_ifcs.cpp, aapl_setup.cpp, aapl_tcond.cpp, apl_atomic_disp.cpp, apl_dirphoncalc.cpp, apl_doscalc.cpp, apl_group_velocity.cpp, apl.h, apl_kphonons.cpp, apl_lrphoncalc.cpp, apl_ltet.cpp, apl_pathbuilder.cpp, apl_phoncalc.cpp, apl_qmesh.cpp, apl_supercell.cpp, apl_thermalpc.cpp, qha_energies.cpp, qha_eoscalc.cpp, qha_gruneisen.cpp, README_AFLOW_APL.TXT, scqha_gruneisen.cpp, scqha_T_freqs.cpp)
        (AUROSTD: aurostd.h, boot.cpp, main.cpp, xmatrix.cpp, xmatrix.h, xscalar.cpp, xscalar.h, xvector.cpp, xvector.h)
        Deleted files: APL/apl_mpmesh.cpp, APL/apl_uniform_mesh.cpp, APL/apl_ltetdos.cpp, APL/apl_rsmdos.cpp
3.1.223 - 2019/04/18
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.223.tar.xz
        -- fixed logics for filename extended characters (SC)
        -- created compliance in aurostd for extended characters and added aurostd::LinkFile (SC)
        -- aurostd::RenameFile => aurostd::file2file (SC)
        -- aurostd::MoveFile => aurostd::file2directory (SC)
        -- --use_tmpfs=XXXXX (for rerouting /tmp directories, useful for [il]logical mapping) (SC)
         (AFLOW: aflowlib_libraries.cpp aflow_ivasp.cpp aurostd.h aurostd_main.cpp)
3.1.222 - 2019/04/01
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.222.tar.xz
        -- added GFA code+README (DF)
        -- added Cygwin support (ME)
        -- fixed bugs in the phonon property plotter (ME)
        -- added defaults/options for DoD MUSTANG machine (--machine=dod_mustang); also added in aflowrc (DX)
        -- remove pseudopotential information for comparing materials (DX)
        -- fixed comparison bugs (e.g., duplicate count, load from URL, printing properties, BrinInCell() for supercell expansion, etc.) (DX)
        -- fixed --aflow_proto bug when generating POCC structures with ANRL parameters (DX)
        -- patched GetDistMatrix() header declaration (CO)
        -- patched x-ray analysis for POSCARs with no atom names (CO)
        -- silenced cematrix::InverseMatrix (CO)
        -- minor patches for clang (CO)
        -- patched fancy print (colors) for POCC AFLOW-Online web mode (CO)
        -- added --scrub=LIBS for lib2raw (SC)
        -- added ProgressBar(std::ostream& oss,string prelim,uint j,uint jmax,bool VERBOSE_PERCENTAGE,bool VERBOSE_ROLLER,bool VERBOSE_CURSOR) in aurostd*
        (AFLOW: aconvasp.cpp, aconvasp_main.cpp, aflow.cpp, aflow.h, aflowrc.cpp, avasp.cpp, bader.cpp, compare_structure.cpp, compare_structure_function.cpp, data.cpp, gfa.cpp, gfa.h, init.cpp, ivasp.cpp, kbin.cpp, kvasp.cpp, Makefile, pflow_funcs.cpp, pflow.h, pflow_print.cpp, pocc.cpp, README_AFLOW_ACONVASP.TXT, README_AFLOW_GFA.TXT, README_AFLOW.TXT, README_CONTRIBS.TXT, symmetry_spacegroup_ITC_library.cpp, xclasses.cpp)
        (AUROSTD: aurostd.h, boot.cpp, xmatrix.cpp)
3.1.221 - 2019/03/20
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.221.tar.xz
        -- changed the compiler for Mac OS to clang++ as g++ does not find the standard libraries on some systems. (ME)
        -- fixed bug in VASP_Produce_POTCAR for POTCAR paths (ME)
        -- APL now switches RELAX to off with --generate_aflowin_only. Also does not check for the vasp binary anymore if POLAR is ON, which would break the code if no vasp binary is present. (ME)
        -- fixed ANRL setting for symmetry cell choice (DX)
        -- added fast supercell function for quick expansion (DX)
        -- speed increase for compare functions (DX)
        -- fixed tolerance scan issue introduced in V3.1.220 after removing global symmetry variables (DX)
        -- added missing directory string in GetSpaceGroupLabel() function (DX)
        -- fixed primitive reduction corner case; should use PBC() instead of BringInCell() for distance vectors (DX)
        -- added ANRL directory name for --aflow_proto: check if ANRL proto matches to entry in library and return ANRL suffix, otherwise add parameter values to directory name (DX)
        -- added FileMESSAGE argument to compare functions (DX)
        -- patched species re-decoration for parent structure, volumes should be consistent across different decorations (CO)
        -- sort unique structures by HNF matrix/site configuration indices so order is always fixed (CO)
        -- removed file writing for command-line POCC commands, speeds up AFLOW-Online (CO)
        -- propagate flags into POCC structures for command-line control (CO)
        -- fixed AFLOW_CHULL_JUPYTER subdirectory creation (CO)
        (AFLOW: aconvasp_main.cpp, aflow.h, aflowlib_libraries.cpp, avasp.cpp, chull.cpp, chull.h, compare_structure.cpp, compare_structure_function.cpp, compare_structure.h, ivasp.cpp, Makefile, pflow.h, pocc.cpp, pocc.h, symmetry_spacegroup.cpp, xatom.cpp, xclasses.cpp)
        (ANRL: aflow_anrl.cpp, ANRL_CPPS_20180710/list.cpp)
        (APL: apl_kphonons.cpp)
3.1.220 - 2019/03/11
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.220.tar.xz
        -- added jupyter|jupyter2|jupyter3 functionality to chull (MB)
        -- added functionality to generate prototypes in CIF format (--cif flag); for --proto and --aflow_proto (DX)
        -- added functionality to generate prototypes in ABCCAR format (--abccar flag); for --proto and --aflow_proto (DX)
        -- added original crystal keywords to aflowlib entry (volume_cell_orig, volume_atom_orig, density_orig, crystal_family_orig, crystal_system_orig, point_group_Hermann_Mauguin_orig, point_group_Schoenflies_orig, point_group_orbifold_orig, point_group_type_orig, point_group_order_orig, point_group_structure_orig, Bravais_lattice_lattice_type_orig, Bravais_lattice_lattice_system_orig, Bravais_superlattice_lattice_type_orig, Bravais_superlattice_lattice_variation_type_orig, Bravais_superlattice_lattice_system_orig, Pearson_symbol_superlattice_orig, reciprocal_geometry_orig, reciprocal_volume_cell_orig, reciprocal_lattice_type_orig, reciprocal_lattice_variation_type_orig, Wyckoff_letters_orig, Wyckoff_multiplicities_orig, Wyckoff_site_symmetries_orig) (DX)
        -- fixed type for point_group_order in JSON; should be int, not string (DX)  
        -- added character check when loading aflow.in, to prevent null bytes (DX) 
        -- added function to remove null bytes in aflow.in file; rewrites cleaned aflow.in (DX) 
        -- added defaults/options for DoD ONYX, GORDON, COPPER, GAFFNEY, and KOEHR machines (--machine=dod_onyx,dod_gordon,dod_copper,dod_gaffney,dod_koehr); also added in aflowrc (DX)
        -- added ANRL space group setting option; monoclinic: unique axis-b, rhombohedral: rhombohedral setting, centrosymmetric: origin centered on inversion (DX)
        -- added functionality to convert a structure into an ANRL designation, i.e., label, parameter list, and parameter values (DX)
        -- added ANRL keywords to aflowlib entry (anrl_label_orig, anrl_parameter_list_orig, anrl_parameter_values_orig, anrl_label_relax, anrl_parameter_list_relax, anrl_parameter_values_relax) (DX)
        -- added functions to get the Wyckoff equations for future symbolic notation (DX)
        -- added ANRL preset values for prototypes in getANRLParameters() (DX)
        -- added ANRL prototypes to README_PROTO.TXT along with unique permutation information (DX)
        -- updated headers for each AFLOW prototype in README_LIBRARY_HTQC.TXT to include the stoichometry, Pearson symbol, space group, and Wyckoff letters (DX)
        -- added functionality to search AFLOW prototypes (--proto_labels) via stoichiometry, space group, arity, etc. (DX)
        -- added AFLUX command line functionality (DX+FR)
        -- added AFLUX helper functions to run inside AFLOW and extract properties into vectors of properties (DX)
        -- added functionality to --wyccar, print letters/multiplicities/site symmetries as a string (DX)
        -- added function to perform symmetry on vector of xstructures (DX)
        -- added function to get enantiomophs of space group, if any exist (DX)
        -- added function to return list of particular element classes, e.g., alkali, transition metals, metals, non-metals, etc. (DX) 
        -- added function to convert compound name into a stoichiometry (DX) 
        -- updated comparison code; added new functionality to compare compounds to AFLOW database, compare compounds to AFLOW prototypes, identify unique permutations, etc. (DX)
        -- added options to ignore symmetry during structure comparison (DX)
        -- improved comparison code output, i.e., more information and cleaner (DX)
        -- improved multithreading of comparison code (DX)
        -- refactor symmetry code to remove extern variables, otherwise it breaks multithreaded symmetry analyses (DX)
        -- refactor symmetry code to remove global variable _SYM_TOL_, otherwise it breaks multithreaded symmetry analyses (DX)
        -- fixed lattice transformation in GetLatticeType(); account for unit cell orientation change using GetPrimitive() (DX)
        -- fixed bug in CIF writer; set VASP version booleans to false (DX)
        -- updated "make check" hash since Wyckoff letters are now in title line; changed from #174c76b2b2928dcdf2f3b39069a8b59 to #ffb9c5681045fb80f391e9a931f21d1 (DX)   
        -- fixed typo in AEL/AGL json functions (CT)
        -- fixed zerostate for AAPL (ME)
        -- small bug fixes for AAPL and APL (ME)
        -- made changes to write(A)APL functions (ME)
        -- added keywords to aflowlib.out/json: title, ldau_type, ldau_l, ldau_u, ldau_j, kpoints_relax, kpoints_static, kpoints_bands_path, kpoints_bands_nkpts (ME)
        -- fixed type declaration bug in aurostd::StringStreamSubst (ME)
        -- fixed bug in QHA with uninitialized Booleans (ME)
        -- accelerated reading of forces from vasprun.xml files for APL and AAPL (ME)
        -- AFLOW not produces POSCAR in the format required by the VASP binary instead of just taking the format in the aflow.in file (ME)
        -- added functionality to read CHGCAR and WAVECAR files between relaxations (ME + Rico Friedrich)
        -- removed duplicate APL functions - apl_hroutines.cpp is now obsolete (ME)
        -- auto-sort by inequivalent atoms in APL (CO + Xiaoyu Wang from UBuffalo)
        -- added sortAtomsEquivalent() to pre-APL/POCC and CONVERT-sprim, -sconv, -niggli, -minkowski for prospective APL calculations: better to sort before relaxations and not have to sort again in the future (CO)
        -- patched sortAtoms*() to include basis and ensure relative order (CO)
        -- added switches for symmetrization of distortions (DISTORTION_SYMMETRIZE) and considering iatoms only (DISTORTION_INEQUIVONLY) in APL (CO)
        -- populate forceConstants matrix in iatoms-sorted-agnostic fashion as we move away from this dependence (CO)
        -- clear pgroupk and pgroupk_xtal symmetries before re-calculating - bug in APL (CO)
        -- patching pc2scpMap in GetSuperCell() (CO)
        -- removed exit's from balanceChemicalEquations(), now throwing exceptions (CO)
        -- overloaded StringsAlphabetic (CO)
        -- set up Makefile with auto-alerts to AFLOW-Forum (CO)
        -- added ClearSymmetry() for _atom (CO)
        -- patched robust structure comparison in POCC (CO)
        -- patched --hnf/--hnfcell routines to work with new POCC (AFLOW-Online) and added --pocc_count_total and --pocc_count_unique (CO)
        -- pre-sort sites for --proto with POCC (CO)
        -- added 'S'+'P' combo functionality to --proto with POCC (CO)
        -- patching robust structure comparison analysis within POCC (CO)
        -- added counts to chull txt and json outputs (CO)
        -- patched logo+count centering on chull doc (CO)
        -- removed dependency of chull .tex on tabu, which breaks in TeX Live 2019 (CO)
        -- fixed image alignment in chull pdf doc (CO)
        -- added options for --aflow_proto, including --run_relax, --run_relax_static, --run_relax_static_bands, --run_static, --run_static_bands, --relax_count (CO)
        (AFLOW: aconvasp_main.cpp, aflow.h, aflowlib.h, aflowlib_libraries.cpp, aflowlib_webapp_bands.js, aflowlib_web_interface.cpp, aflowrc.cpp, avasp.cpp, bader.cpp, chull.cpp, chull.h, chull_jupyter.json, chull_jupyter_plotter.py, chull_jupyter_requirements.txt, chull_python.py, compare_structure.cpp, compare_structure_function.cpp, compare_structure.h, init.cpp, ivasp.cpp, kaims.cpp, kalien.cpp, kbin.cpp, kvasp.cpp, lattice.cpp, Makefile, matlab.cpp, modules.cpp, ovasp.cpp, pflow_funcs.cpp, pflow.h, pocc.cpp, pocc.h, pocc_old.cpp, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT, README_AFLOW_COMPARE.TXT, README_AFLOW.TXT, README_LIBRARY_HTQC_BORIDES.TXT, README_LIBRARY_HTQC_CARBIDES.TXT, README_LIBRARY_HTQC.TXT, README_PROTO.TXT, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_functions.cpp, symmetry_spacegroup.h, symmetry_spacegroup_ITC_library.cpp, xatom.cpp, xclasses.cpp, xproto.cpp, xproto_gus.cpp)
        (ANRL: aflow_anrl.cpp, ANRL_CPPS_20180710/list.cpp, README_AFLOW_ANRL.TXT)
        (APL: aapl_cluster.cpp, aapl_ifcs.cpp, aapl_setup.cpp, aapl_tcond.cpp, apl_dirphoncalc.cpp, apl.h, apl_kphonons.cpp, apl_lrphoncalc.cpp, apl_pathbuilder.cpp, apl_phoncalc.cpp, apl_shellhandle.cpp, apl_supercell.cpp, README_AFLOW_APL.TXT)
        (AUROSTD: aurostd.h, main.cpp)
3.1.219 - 2019/01/20
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.219.tar.xz
        -- small modifications for --aflow_proto (--generate_aflowin_only) that enable aflow.in generation on any mac (CO)
        -- fixed JSON output chull (CO)
        (AFLOW: aconvasp_main.cpp, avasp.cpp, chull.cpp, Makefile)
3.1.218 - 2019/01/14
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.218.tar.xz
        -- added defaults/options for CMU EULER machine (--machine=cmu_euler); also added in aflowrc (DX)
        -- modified AEL/AGL aflow.in keys (CT)
        -- aurostd::xoption: fixed bug in addattachedscheme, added additional Boolean functionality in options2entry. (ME)
        -- Added debug options to AAPL. (ME)
        -- Improved APL and AAPL defaults. (ME)
        -- Integrated APL and AAPL into the standardized aflow.in creator for VASP calculations. (ME)
        -- Updated the APL README. (ME)
        -- Fixed bugs for APL and AAPL output file names. (ME)
        -- Updated the PREC=PHONONS INCAR parameters. (ME)
        -- Added RELAX feature to APL: structures can now be relaxed specifically for phonon calculations inside the APL environment. (ME)
        -- Added features to handle NCPUS=MAX. (ME)
        -- *NEW* automatic aflow.in generator for APL/AAPL/POCC (beta) (CO+ME)
        -- integrated POCC into the standardized aflow.in creator for VASP calculations (ME)
        -- added pocc params and tol to --proto (CO)
        -- added aflags to APL so sym analysis works with -D (CO)
        -- added SYMMETRIZE=OFF option for APL distortions (CO)
        -- added --aflow_proto options: --bader, --spin_remove_relax_1 _2, --kscheme _static, --kppra _static, --relax_count=XX, --run_relax_static, --run_relax_static_bands (CO)
        -- patched PrototypeLibrariesSpeciesNumber() to handle ANY ICSD (CO)
        -- added -001 for select ANRL prototypes (CO)
        -- added pseudopotential information to --aflow_proto directory + system name (CO)
        -- fixed m_initialized in chull/pocc (CO)
        -- added chull plotting functionality: if unstable, figure out how far above hull to plot automatically (CO)
        -- modified plot axes titles in chull (CO)
        -- added more fixes for "Reciprocal lattice and k-lattice belong to different class of lattices", off until further testing (CO)
        -- moved pocc settings to aflow.rc (CO)
        -- now write out all derivative structures, unique derivative structures, hnf matrices, and site configurations in pocc (CO)
        -- fixed pocc partial occupancy optimizer table settings (CO)
        -- added AEL/AGL settings to aflow.in generator (CO)
        (AFLOW: aconvasp_main.cpp, aflow.h, aflowlib.h, aflowlib_libraries.cpp, aflowlib_web_interface.cpp, aflowrc.cpp, avasp.cpp, chull.cpp, contrib_shidong.cpp, ifrozsl.cpp, init.cpp, ivasp.cpp, kaims.cpp, kbin.cpp, kvasp.cpp, Makefile, modules.cpp, neighbours.cpp, pflow.h, pocc.cpp, pocc.h, pocc_old.cpp, pthreads.cpp, README_AFLOW_ACONVASP.TXT, README_AFLOW_AFLOWRC.TXT, symmetry.cpp, xatom.cpp, xclasses.cpp, xproto.cpp)
        (ANRL: aflow_anrl.cpp)
        (AEL: ael_elasticity.cpp, README_AFLOW_AEL.TXT)
        (AGL: agl_debye.cpp, README_AFLOW_AGL.TXT)
        (APL: aapl_setup.cpp, aapl_tcond.cpp, apl_dirphoncalc.cpp, apl_doscalc.cpp, apl.h, apl_hsqpoints.cpp, apl_kphonons.cpp, apl_lrphoncalc.cpp, apl_pdisc.cpp, apl_phoncalc.cpp, apl_supercell.cpp, apl_thermalpc.cpp, qha_aflowin_creator.cpp, qha_eosrun.cpp, README_AFLOW_APL.TXT)
        (AUROSTD: aurostd.h, boot.cpp, main.cpp, xerror.cpp, xoption.cpp)
3.1.217 - 2018/12/13
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.217.tar.xz
        -- Added LIB2LIB function to automatically perform AEL and AGL post processing when running LIB2RAW (CT)
        -- Updated LIB2RAW and webpage generation functions to incorporate additional AEL and AGL properties (CT)
        -- Added conversion of pressure extracted from OUTCAR file from kB to GPa in AEL and AGL to keep units consistent (CT)
        -- added species scaling to volume for ANRL prototypes, if --params=-1.0,... (DX)
        -- added ANRL prototype for the kesterite structure (DX)
        -- improved CIF reader, i.e., account for different formats and partial occupation (DX)
        -- fixed bug in kpath determination (primitive vs original lattice) (DX)
        (AFLOW: aconvasp_main.cpp, anrl.cpp, xatom.cpp, AFLOW_README_PROTO.TXT, aflow.h, aconvasp_main.cpp, ael_elastic_fit.cpp, ael_elasticity.cpp, agl_debye.cpp, agl_electronic.cpp, init.cpp, xclasses.cpp, aflowlib.h, aflowlib_libraries.cpp, aflowlib_web_interface.cpp, Makefile, README_PROTO.TXT)
        (ANRL: ANRL_CPPS_20180710/aflow_anrl_A2BCD4_tI16_82_ac_b_d_g.cpp, ANRL_CPPS_20180710/aflow_anrl_list.cpp)
3.1.216 - 2018/12/05
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.216.tar.xz
        -- fixed std in AAPL (ME)
        -- added seven SQS structures to ANRL prototypes (DX)
        -- Write Gibbs free energy as a function of pressure and temperature in AFLOW readable format for future phase diagram application (CT)
        -- Add option to suppress extrapolation in Hugoniot relation calculation (CT)
        -- Functions to extract electronic properties as a function of pressure from AGL data (CT)
        -- Calculate and write additional elastic properties from AEL (Young's modulus, Pugh's modulus ratio) (CT)
        -- Flattened JSON structure used for writing and storing elastic stiffness and compliance tensors (CT)
        -- fixed anisotropy spelling (CO)
        -- added ISMEAR/SIGMA VASP option (CO)
        -- rerouted GetAtomVolume() and GetAtomMass() for properly cleaning pp (CO)
        -- added N+1 analysis to chull (CO)
        -- removed _pc from QH_ENERGIES() class (CO)
        (AUROSTD: AUROSTD/aurostd.h, main.cpp, xcombos.cpp, xcombos.h)
        (ANRL: ANRL_CPPS_20180710/aflow_anrl_A3B13_oC32_38_ac_a2bcdef.cpp, ANRL_CPPS_20180710/aflow_anrl_A3B5_oC32_38_abce_abcdf.cpp, ANRL_CPPS_20180710/aflow_anrl_A5B11_mP16_6_2abc_2a3b3c.cpp, ANRL_CPPS_20180710/aflow_anrl_AB3_mC32_8_4a_12a.cpp, ANRL_CPPS_20180710/aflow_anrl_AB3_mC32_8_4a_4a4b.cpp, ANRL_CPPS_20180710/aflow_anrl_AB7_hR16_166_c_c2h.cpp, ANRL_CPPS_20180710/aflow_anrl_AB_aP16_2_4i_4i.cpp, ANRL_CPPS_20180710/aflow_anrl_list.cpp)
        (AFLOW: ael_elastic_fit.cpp, ael_elasticity.cpp, ael_elasticity.h, agl_debye.cpp, agl_debye.h, agl_electronic.cpp, agl_hugoniot.cpp, anrl.cpp, avasp.cpp, chull.cpp, chull.h, aflow.cpp, aflow.h, ivasp.cpp, kaims.cpp, kvasp.cpp, aflowlib.h, aflowlib_libraries.cpp, aflowlib_web_interface.cpp, pocc.cpp, pocc.h, xatom.cpp, xclasses.cpp, xproto.cpp, xproto_gus.cpp, APL/aapl_cluster.cpp, APL/aapl_ifcs.cpp, APL/aapl_tcond.cpp, APL/apl.h, APL/apl_kphonons.cpp, APL/qha_energies.cpp, Makefile, README_AFLOW_AEL.TXT, README_AFLOW_AGL.TXT, README_AFLOW.TXT, README_CONTRIBS.TXT)
3.1.215 - 2018/12/04
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.215.tar.xz
        -- added machinery for vaiud (auid bytes for the AUID directory cashed) (SC)
	      aflowlib.h aflowlib_libraries.cpp aflowlib_web_interface.cpp aflow.h aflow_aflowrc.cpp aflow_init.cpp
3.1.214 - 2018/11/30
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.214.tar.xz
        -- added auid.out and auid.json to RAW/WEB (SC)
3.1.213 - 2018/11/27
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.213.tar.xz
        -- added LIB0 (a bunch of h and cpp) (SC)
3.1.212 - 2018/11/12
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.212.tar.xz
        -- in --bzd command, added option to transform high-symmetry kpaths to input lattice representation with [--transform2original] (DX)
        -- in --bzd command, added option to print transformation matrices between input lattice and AFLOW standard lattice [--print_transformation_matrix] (DX)
        -- store both the coordinate system and unit cell (rigid rotation) transformations in xstructure (DX)
        -- created _kpoints class; easier to transform and print kpoint information (DX)
        (AFLOW: aflow.h, aconvasp_main.cpp, lattice.cpp, pflow.h, xatom.cpp, README_AFLOW_ACONVASP.TXT)
3.1.211 - 2018/10/19
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.211.tar.xz
        -- added Pearson coefficient to element properties (ME)
        -- fixed some constants in xscalar (ME)
        -- Rewrote AAPL to be faster, require less DFT calculations, and include four-phonon processes. (ME)
        -- Added more functionality to xvector and xmatrix for complex numbers. (ME)
        -- Added ability to update progress bar using percentages instead of indices (ME)
        -- Redesigned xtensor to be more lightweight and faster. (ME)
        -- Fixed bug in aurostd::xcombos::reset() (ME)
        -- Fixed minor typos (ME)
        -- Edited APL readme for grammar and clarity, and added new AAPL features (ME)
        -- Added APL and AAPL parameters to the aflow.rc (ME)
        -- Reformatted exception readme to be more consistent with other readmes and fixed typos (ME)
        -- updated get_datetime_formatted() and StringsAlphabetic() (CO)
        -- new chull date format avoids time stamp, too much resolution (CO)
        -- added png output option for chull and resolution option to aflowrc (CO)
        -- switched to tight tol spacegroups (vsg2) in chull (CO)
        -- centralized checking sign of distances in chull, then flipped sign of distance to hull to be positive by default (CO)
        -- funneled all points of facet through addVertex() (CO)
        -- fixed facet content tolerance scaling with increasing dimensionality, important for 6D hulls (CO)
        -- added default initializations to a bunch of variables in chull to avoid spurious warnings (CO)
        -- added polymorph case to decomposition phases and coefficients (CO)
        -- allow for eq_g_states to be calculated on the fly if not already calculated (CO)
        -- intelligently avoid printing unknown (or unset) stability criterion (CO)
        -- fixed unsorted directories in aflow_compare_structure (CO)
        -- now check for negative coefficients when balancing chemical equations, means we missed the facet (CO)
        -- correctBadDatabase() now checks for unaries too (as per apennsy) (CO)
        -- skipping ".old" prototypes (CO)
        -- fixed shidong warnings with strncat (CO)
        (AUROSTD: boot.cpp, AUROSTD/aurostd.h, main.cpp, xcombos.cpp, xcombos.h, xcomplex.cpp, xcomplex.h, xmatrix.cpp, xmatrix.h, xtensor.cpp, xtensor.h, xvector.cpp, xvector.h)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, chull.cpp, chull.h, compare_structure.cpp, contrib_shidong_auxiliary.cpp, aflowlib_web_interface.cpp, pflow_funcs.cpp, APL/aapl_cluster.cpp, APL/aapl_ifcs.cpp, APL/aapl_setup.cpp, APL/aapl_tcond.cpp, APL/apl_atomic_disp.cpp, APL/apl_dirphoncalc.cpp, APL/apl_dm_pdos_save.cpp, APL/apl_doscalc.cpp, APL/apl_group_velocity.cpp, APL/apl_gsa.cpp, APL/apl.h, APL/apl_hsqpoints.cpp, APL/apl_kphonons.cpp, APL/apl_logger.cpp, APL/apl_lrphoncalc.cpp, APL/apl_pdisc.cpp, APL/apl_phoncalc.cpp, APL/apl_supercell.cpp, APL/apl_thermalpc.cpp, APL/qha3phonons_eos.cpp, APL/qha_aflowin_creator.cpp, APL/qha_dm_pdos_save.cpp, APL/qha_eosrun.cpp, APL/qha_gruneisen.cpp, APL/scqha_gruneisen.cpp, APL/scqha_T_freqs.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_APL.TXT, README_AFLOW_CHULL.TXT, README_AFLOW_EXCEPTIONS.TXT)
3.1.210 - 2018/10/01
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.210.tar.xz
        -- changed 8 ANRL prototype labels (part 2) to match with Wyckoff positions in reference literature (DX)
        -- fixed prototype names for A4BC4D_tP10_123_gh_a_i_d and AB_hP6_144_a_a in aflow_anrl_list.cpp (DX) 
        -- added Strukturbericht designations to aflow_anrl_list.cpp (DX)
        -- fixed minimum enumerated Wyckoff letter determination; more robust (DX)
        -- cleaned workflow for Wyckoff functions (DX)
        -- added more debugging messages to symmetry functions (DX)
        (AFLOW: ANRL/aflow_anrl_list.cpp ANRL/aflow_anrl_A12B36CD12_cF488_210_h_3h_a_fg.cpp ANRL/aflow_anrl_A3B3C_hP14_176_h_h_c.cpp ANRL/aflow_anrl_A3B_hP8_176_h_c.cpp ANRL/aflow_anrl_AB3C_cP60_201_be_fh_g.cpp ANRL/aflow_anrl_A3B2_hP10_176_h_bc.cpp ANRL/aflow_anrl_A3BC_hP10_188_k_c_a.cpp ANRL/aflow_anrl_AB3C16_cF160_203_a_bc_eg.cpp ANRL/aflow_anrl_AB_hP4_156_ab_ab.cpp Makefile anrl.cpp symmetry_spacegroup.cpp symmetry_spacegroup.h symmetry_spacegroup_functions.cpp)
3.1.209 - 2018/08/31
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.209.tar.xz
        -- added all origin choice 2 possibilities to ITC space group list (DX)
        -- added all rhombohedral setting possibilities to ITC space group list (DX)
        -- added monoclinic unique axis choices (b or c) to ITC space group list (DX)
        -- added SYM::findRhombohedralSetting() in space group functions to distinguish between hexagonal and rhombohedral settings (DX)
        -- added other settings for Hall space group symbol (settings 1 or 2, H or R for rhombohedral, and unique axis-b or -c for monoclinic systems) (DX)
        -- added CIF reader; can read/expand CIFs with only representative Wyckoff position specified in settings 1 or 2 (H or R for rhombohedral systems and unique axis-b or -c for monoclinic systems) (DX)
        -- improved CIF writer; default functionality calculates the symmetry and prints the representative Wyckoff positions (DX)
        -- added symmetry tolerance option to CIF writer via --cif[=<tolerance>] (default: tight; other options: loose, <number>) (DX)
        -- added space group setting option to CIF writer via [--setting=1| =2] (default: 1) (DX)
        -- added [--no_symmetry] option to CIF writer, which returns CIF as space group 1 (DX)
        -- added space group setting option to --wyccar via [--setting=1| =2] (default: 1) (DX)
        -- added space group setting option to --sgdata and --edata via [--setting=1| =2] (default: 1) (DX)
        -- added [--no_scan] and [--magmom] options to --wyccar and moved function to pflow::WYCCAR() (DX)
        -- added more debugging messages to space group functions (ConventionalCell(), SpaceGroup_ITC(), etc.) (DX)
        -- changed ANRL prototype A6B2C_tP18_128_eh_d_b to A6B2C_tP18_128_eh_d_a (b vs a Wyckoff letter); consistent with reference (DX)
        -- fixed bug in minimumDistanceVector(); should return xvector<double> not double (DX)
        (AFLOW: ANRL/aflow_anrl_list.cpp, ANRL/aflow_anrl_A6B2C_tP18_128_eh_d_a.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_SYM.TXT, aflow.h, aconvasp_main.cpp, anrl.cpp, pflow.h , pflow_print.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup.h, symmetry_spacegroup_ITC_library.cpp, symmetry_spacegroup_functions.cpp, xatom.cpp)
3.1.208 - 2018/08/27
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.208.tar.xz
        -- new small banner style in chull (CO)
        -- perform thermo loop with static run in aflowlib_libraries (CO)
        -- added xoption to edata and sgdata functions; stores all data (DX)
        -- added more keywords to aflowlib entry (crystal_family, crystal_system, point_group_Hermann_Mauguin, point_group_Schoenflies, point_group_orbifold, point_group_type, point_group_order, point_group_structure, Bravais_lattice_lattice_type, Bravais_lattice_lattice_system, Bravais_superlattice_lattice_type, Bravais_superlattice_lattice_variation_type, Bravais_superlattice_lattice_system, Pearson_symbol_superlattice, reciprocal_geometry, reciprocal_volume_cell, reciprocal_lattice_type, reciprocal_lattice_variation_type, Wyckoff_letters, Wyckoff_multiplicities, Wyckoff_site_symmetries) (DX)
        -- new aflowlib keywords extracted from xoption and integrated into aflowlib.out/.json (DX)
        -- prepared website function to read in new keywords once database is populated; currently commented out (DX)
        -- fixed typo in sgdata; Shoenflies is now Schoenflies (DX)
        -- added functions to create Wyckoff strings in aflowlib entry; ExtractWyckoffLettersString(), ExtractWyckoffMultiplicitiesString(), and ExtractWyckoffSiteSymmetriesString() (DX)
        (AUROSTD: xoption.cpp)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, pflow.h, pflow_print.cpp, symmetry_spacegroup.h, symmetry_spacegroup_functions.cpp, aflowlib.h, aflowlib_libraries.cpp, aflowlib_web_interface.cpp)
3.1.207 - 2018/08/19
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.207.tar.xz
        -- integrated new html for entry page (JPO)
        -- fixed xaxis of bands plotter for entry page (PC)
        -- Improved distribution of APL calculations over threads. (ME)
        -- Added function to SYM that returns the minimum distance vector. (ME)
        -- Added docstring to xcombos. (ME)
        -- added xvector::normalizeSumToOne() (CO)
        -- changed decomposition reaction to atomic concentrations (CO)
        -- added fractional_compound to chull output (CO)
        -- replaced AFLOWLogicError() and AFLOWRuntimeError() with xerror() (CO)
        -- added chull plot ICSD labels mode (CO)
        -- revamped small banner setting (CO)
        -- added verbose output if skipping entries above/below half hull (CO)
        -- added --output=png --png_resolution=500 options to chull (CO)
        -- added --keep=gpl to phonons gnuplot script (CO)
        (AUROSTD: boot.cpp, AUROSTD/aurostd.h, xcombos.cpp, xscalar.cpp, xscalar.h, xvector.cpp, xvector.h)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, avasp.cpp, bader.cpp, chull.cpp, chull.h, aflow.cpp, aflow.h, aflowlib_webapp_bands.js, aflowlib_web_interface.cpp, pflow_funcs.cpp, pflow.h, pocc.cpp, symmetry.cpp, APL/aapl_tcond.cpp, APL/apl_atomic_disp.cpp, APL/apl_dm_pdos_save.cpp, APL/apl_doscalc.cpp, APL/apl_group_velocity.cpp, APL/apl.h, APL/apl_hroutines.cpp, APL/apl_pdisc.cpp, APL/qha3phonons_eos.cpp, APL/qha_dm_pdos_save.cpp, APL/qha_gruneisen.cpp, APL/scqha_eos.cpp, APL/scqha_gruneisen.cpp, APL/scqha_T_freqs.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_AFLOWRC.TXT, README_AFLOW_CHULL.TXT)
3.1.206 - 2018/08/08
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.206.tar.xz
        -- Fix JVXL (SC)
3.1.205 - 2018/07/27
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.205.tar.xz
        -- Improved destructor for xtensor (ME)
        -- Added function to test if xvector is a zero vector (ME)
        -- fixed parameter list for ANRL prototype (A_hR105_166_bc9h4i); z2 to x2 (DX)
  	-- more rigorous check for atoms within the new cell in GetSuperCell(); check periodic images of atoms (DX)
  	-- fixed "over-reduced" issue in primitivization routine in space group function (DX)
  	-- added roundoff to axis and SU2 matrix in JSON output of symmetry elements (DX)
  	-- fixed function for printing fractions in general Wyckoff positions (DX)
        -- created --generate_aflowin_only (CT)
        -- QHA3P and SCQHA conflict resolved (PN)
        -- reorganize QHA modes (PN)
        -- added four QHA3P options (PN)
        -- reorganize QHA, SCQHA, and QHA3P options (PN)
        -- filename changed and accordingly modified in QHA-README (PN)
        -- replaced some QHA functions with aurostd (PN)
        -- fixed slab AddAtom() function (CO+DU)
        -- added --readme=aflowrc (CO)
        -- fixed --readme=xaflow for local configuration (CO)
        -- fixed a few warnings for beta (CO)
        -- fixed AMIX/BMIX typo (CO)
        -- added defaults/options for DOD CONRAD machine (--machine=dod_conrad); also added in aflowrc (DX)
        (AUROSTD: xscalar.cpp, xscalar.h, xtensor.cpp, xvector.cpp, xvector.h)
        (AFLOW: aflowrc.cpp, avasp.cpp, chull.cpp, aflow.cpp, data.cpp, aflow.h, init.cpp, ivasp.cpp, kbin.cpp, kvasp.cpp, aflowlib_webapp_entry.js, aflowlib_web_interface.cpp, pocc_old.cpp, surface.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_functions.cpp, xatom.cpp, xclasses.cpp, ANRL_CPPS_20180710/anrl_list.cpp, APL/apl_group_velocity.cpp, APL/apl.h, APL/apl_kphonons.cpp, APL/qha3phonons_eos.cpp, APL/qha_aflowin_creator.cpp, APL/qha_energies.cpp, APL/qha_eos.cpp, APL/qha_gruneisen.cpp, APL/scqha_eos.cpp, APL/scqha_gruneisen.cpp, APL/scqha_T_freqs.cpp, Makefile, README_AFLOW_AEL.TXT, README_AFLOW_AFLOWRC.TXT, README_AFLOW_AGL.TXT, README_AFLOW_POCC.TXT, README_AFLOW_QHA_SCQHA_QHA3P.TXT, README_AFLOW.TXT, README_AFLOW_XAFLOW.TXT)
3.1.204 - 2018/07/12
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.204.tar.xz
        -- updated QHA and added QHA3P and SCQHA functionality (PN)
        -- Restrucured xcombos. Added enumerations. (ME) [AUROSTD/aflow_xcombos.cpp, AUROSTD/aflow_xcombos.h]
        -- Introduced exception handlier class aurostd::xerror. (ME)
        -- Introduced xtensor class for tensors of arbitrary dimension. Other xtensor classes are obsolete now. (ME).
        -- Changed all current xtensor3 instances to the new xtensor format. (ME)
        -- Added Kronecker product to xmatrix. (ME)
        -- Option to write and use AEL data at lowest finite pressure where the material is elastically stable (CT)
        -- Option to specify a separate set of finite pressures for AEL calculations than are used for AGL post-processing (CT)
        -- Integrated workflow option for using finite pressure Poisson ratio in AGL calculations (CT)
        -- Write nominal target pressure and calculated external pressure for AEL calculations in aflow.ael.out file (CT)
        -- Fixed keyword in JSON output file (CT)
        -- Increase grid for AGL DOSCAR (CT)
        -- cleaning up webapp_bands.js and added more mouse functionality (PC)
        -- added citation information to entry page (PC)
        -- added ANRL prototypes from part 2 of library (302 prototypes) (DX)
        -- added option to print symbolic math representation of ANRL prototypes (--add_equations or --equations_only) in aims and vasp formats (DX)
        -- updated DOI for ANRL part 1 and added arXiv for part 2 to title line of each prototype (DX) 
        -- fixed bug with --vasp keyword, i.e., need to check if used with --proto command (DX) 
        -- removed unused variables in standard lattice function (DX)
        -- added applyCombo() to xcombos (CO)
        -- fixed wget *.xz issue if *.xz already exists (CO)
        -- fixed apl2agr to handle xz compression format (CO)
        -- fixed --aflow_proto empty BZ issue, occurs when structure is downloaded from online (CO)
        -- fixed reading in forces from aims.out (CO)
        -- force species input for aims structure (CO)
        -- added checks for broken API (CO)
        -- added GENERATOR to aflow.in generation (CO)
        -- added auid to chull PDF output (CO)
        -- fixed bader num_each_type issue (relax2 vs. static primitivization) (CO)
        -- fixed command line aflow_qmvasp generation (CO)
        -- fixed search for ./aflow_data issue upon initial installation (CO)
        -- PARTCAR now handles '+' in write out (CO)
        -- fixed AddAtom() to handle different occupations (CO)
        -- added combination parameters to chull (CO)
        -- fixed aflowrc load issue (remove spaces between quotes and comment) (CO)
        (AUROSTD: boot.cpp, AUROSTD/aurostd.cpp, AUROSTD/aurostd.h, main.cpp, xcombos.cpp, xcombos.h, xerror.cpp, xerror.h, xmatrix.cpp, xmatrix.h, xtensor.cpp, xtensor.h)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, anrl.cpp, avasp.cpp, bader.cpp, bader.h, chull.cpp, chull.h, aflow.cpp, data.cpp, aflow.h, init.cpp, ivasp.cpp, kbin.cpp, lattice.cpp, aflowlib_libraries.cpp, aflowlib_webapp_bands.js, aflowlib_web_interface.cpp, oaims.cpp, ovasp.cpp, pflow.h, pocc.cpp, pocc.h, surface.cpp, xatom.cpp, xclasses.cpp, xproto.cpp, ANRL_CPPS_20180710, APL/apl_atomic_disp.cpp, APL/apl_doscalc.cpp, APL/apl_group_velocity.cpp, APL/apl.h, APL/apl_hroutines.cpp, APL/apl_hsqpoints.cpp, APL/apl_kphonons.cpp, APL/apl_ltetdos.cpp, APL/apl_mpmesh.cpp, APL/apl_pdisc.cpp, APL/apl_phoncalc.cpp, APL/apl_thermalpc.cpp, APL/apl_uniform_mesh.cpp, APL/qha3phonons_eos.cpp, APL/qha_aflowin_creator.cpp, APL/qha_dm_pdos_save.cpp, APL/qha_energies.cpp, APL/qha_eoscalc.cpp, APL/qha_eos.cpp, APL/qha_gruneisen.cpp, APL/scqha_eos.cpp, APL/scqha_gruneisen.cpp, APL/scqha_T_freqs.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_AEL.TXT, README_AFLOW_AGL.TXT, README_AFLOW_ANRL.TXT, README_AFLOW_APL.TXT, README_AFLOW_CHULL.TXT, README_AFLOW_EXCEPTIONS.TXT, README_AFLOW_POCC.TXT, README_AFLOW_QHA_SCQHA_QHA3P.TXT, README_AFLOW.TXT, README_AFLOW_XAFLOW.TXT, README_CONTRIBS.TXT)
3.1.203 - 2018/06/15
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.203.tar.xz
        -- fixed tolerance scan issue in edata (try new tolerance from PrintSGData() on GetLatticeType() routine before changing tolerance) (DX)
        -- fixed primitivization routine in aflowSG functions (Minkowski/Niggli to fix left-handed candidate lattices and recheck moduli after Minkowski/Niggli) (DX)
        -- check all possible generator choices to match to ITC convention before changing tolerance (DX)
        -- added more rigorous check of tetragonal symmetry operations for determining conventional cell (DX)
        -- added generator information for P1 symmetry systems (DX)
        -- added check of Cartesian distance before removing fractional copies in GetPrimitive() (DX)
        (AFLOW: lattice.cpp, pflow_print.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_ITC.cpp, xatom.cpp) 
3.1.202 - 2018/06/07
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.202.tar.xz
        -- dropping EXTRA to AFLOW3_FREE/EXTRA on local machine or wget (SC)
3.1.201 - 2018/06/04
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.201.tar.xz
        -- more space saving in lib2raw (linking OUTCAR.relax instead of copying) SC
3.1.200 - 2018/05/27
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.200.tar.xz
        -- fixed structure rescaling issue for space group determination (affecting LIB4) (DX)
        -- print geometry file location for errors/debug (DX)
        -- added more debugging messages along with file location in symmetry functions (DX) 
        -- initialize variables for -O3 in symmetry and structure comparison routines (DX)
        -- speed increase for minimumDistance() function (DX)
        -- speed increase for primitivization routine in aflowSG functions (DX)
        -- fixed Wyckoff position typo for space group 62 (8d not 8e) (DX)
        -- fixed Wyckoff position typo for space group 89 (4i not 2i) (DX)
        (AFLOW: aflow.h, aconvasp_main.cpp, compare_structure.cpp, compare_structure_function.cpp, pflow_print.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_ITC_library.cpp, symmetry_spacegroup_functions.cpp, xatom.cpp, xproto.cpp, aflowlib_libraries.cpp)
3.1.199 - 2018/05/27
       Location: http://materials.duke.edu/AFLOW/aflow.3.1.199.tar.xz
	      -- compressing aflow.***.json/out in LIB2RAW and linking them (SC)
        -- removing the brainy/useless stuff about file compression in aurostd (SC)
	      -- fighting sloppyness in aflow_libraries about aflow_pgroup[x][_xtal] (SC)
	      -- fixed aflow_libraries about compress, delete and link files (SC)
	      -- fixed inconsitencies in aflow_convasp_main compress (SC)
	      -- fixed as bunch of inconsistencies in compressing and XHOST.command (SC)
	      -- more inconsistencies RAW-WEB in aflow_libraries  (SC)
	      -- fixed removal useless files in LIB2RAW (SC)
	      -- added BZ2XZ engine to aurostd_main (SC)
	      -- added GZ2XZ engine to aurostd_main (SC)
	      -- added ZIP2ZIP engine to aurostd_main (SC)
3.1.198 - 2018/05/24
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.198.tar.xz
        -- fixed bader extension finder for --lib2raw (CO)
        (AFLOW: bader.cpp, bader.h, aflowlib_libraries.cpp, Makefile)
3.1.197 - 2018/05/24
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.197.tar.xz
        -- updates to entry page bands plotter (PC)
        -- added some helpful comments for getGeneralNormal(), CMdet(), and ZVAL in bader (CO)
        -- check if /www directory exists for jmol display on entry page (CO)
        (AUROSTD: xmatrix.cpp, xvector.cpp)
        (AFLOW: bader.cpp, chull.cpp, aflowlib_webapp_bands.js, aflowlib_web_interface.cpp, Makefile)
3.1.196 - 2018/05/21
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.196.tar.xz
        -- fixed entry page property line wrapping issues (PC)
        -- added button for spin-polarized band structure selection in webapp (PC)
        -- fixed precision inconsistency printing xstr.json (FK)
        -- added some -O3 compatibility (FK)
        -- citation added to aflow_aapl_pairs.cpp (ME)
        -- fixed bug in aflow_kvasp.cpp that prevented the creation of the primitive cell structure for APL, AAPL, and QHA calculations (ME)
        -- implemented combinations with repetitions with sequence taken into account (ME)
        -- moved the code to check for input and output files in APL into separate functions to make them available for AAPL, remove duplicate code, and make the code more readable (ME)
        -- removed the ENCUT and EDIFF tags from AAPL input files as they may result in lower cut-off energies and higher energy differences (ME)
        -- fixed gcc8 issue in AGL (CO)
        -- added simplex content and hypercollinearity properties to convex hull data (CO)
        -- decoupled internal links between graph2doc and withindoc (CO)
        -- fixed site error calculation in POCC for vacancies (CO)
        -- added eyes/ones xmatrix constructors (CO)
        -- added Cayley-Menger determinant to xmatrix (CO)
        -- included cstdlib in aflow_data.cpp for compilation on qrats (CO)
        (AUROSTD: boot.cpp, xcombos.cpp, xcombos.h, xmatrix.cpp, xmatrix.h, xvector.cpp)
        (AFLOW: aconvasp_main.cpp, apennsy_main.cpp, chull.cpp, chull.h, contrib_shidong_auxiliary.cpp, data.cpp, aflow.h, kvasp.cpp, aflowlib.h, aflowlib_libraries.cpp, aflowlib_webapp_bands.js, aflowlib_web_interface.cpp, pflow_funcs.cpp, pflow.h, pocc.cpp, pocc.h, pocc_old.cpp, pocc_old.h, xatom.cpp, APL/aapl_tensor.cpp, APL/apl_dirphoncalc.cpp, APL/apl.h, APL/apl_phoncalc.cpp, APL/apl_thermalpc.cpp, Makefile)
3.1.195 - 2018/05/21
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.195.tar.xz
        -- beta of release with XZ (SC)
3.1.194 - 2018/05/16
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.194.tar.xz
        -- small bug fixes for g++/gcc 7 and 8 (SC)
        -- preparing for xz compression (SC) 
        194 xatom.cpp  aurostd.h aurostd_main.cpp pthread.cpp README_SCRIPTING README_AFLOW init.cpp aflow.cpp aconvasp_main.cpp
        194c fix aconvasp_main.cpp
	194d avasp.cpp ael_elasticity.cpp ael_get_stress.cpp agl_debye.cpp agl_get_ev.cpp contrib_cormac.cpp contrib_junkai_phasediag.cpp
	194e avasp.cpp xatom.cpp  aurostd_main libraries.cpp ael_elasticity.cpp ael_get_stress.cpp agl_debye.cpp agl_get_ev.cpp contrib_cormac.cpp aconvasp_main.cpp pthread.cpp init.cpp aflow.cpp (heavy)
	194f good bye bzip2 (SC)
	194g aflow_contrib_wahyu.cpp 
	194h aflow.cpp aflow_estructure.cpp aflow_ifrozsl.cpp aflow_ivasp.cpp
	194i aflow_pthreads.cpp
	194j aflow_aflowrc.cpp
	194l BASE64 for pseudopotentials and aflow_data.cpp
	194m aflow_kbin.cpp aflow_kvasp.cpp
	194n aflowlib_web_interface.cpp aflow_matlab_funcs.cpp (EXTRA/MATLAB/plotband.m) aflow_ovasp.cpp aflow_pocc_edos.cpp
	shortened aflow_xatom.cpp working on ZVAL POMASS and removing all EXT stuff. Tests work.  Fixed even further aflow_pthread.cpp aflow_xproto_gus.cpp
3.1.193 - 2018/05/11
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.193.tar.xz
        -- extending xoptions push pop (SC)
3.1.192 - 2018/05/10
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.192.tar.xz
        -- extending xoptions push pop (SC)
3.1.191 - 2018/05/08
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.191.tar.xz
        -- rationalized orthogonality search in xmatrix (SC)
        aurostd_xmatrix.h aurostd_xmatrix.cpp aflow.h
3.1.190 - 2018/05/08
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.190.tar.xz
        -- fixed units in AGL output (CT)
        -- fixed permutation vector initialization in APL (required for compiling on DoD CONRAD machine) (DX)
        -- fixed bug in tolerance scan; was only scanning in one direction (DX) 
        -- created faster minimum cartesian distance calculator for skewed cells; fewer duplicate operations (DX) 
        -- fixed conflict between CUT_RAD and CUT_SHELL in AAPL (ME)
        -- properly added directory to bader error output (CO)
        -- fixed slow down with loadEntries() (stringElements2VectorElements()) (CO)
        -- fixed AIMS read-in issue with '#' comments (CO)
        -- added guard around BANDSDATA_JSON() for html generation (CO)
        -- fixed static compile settings (CO)
        -- fixed _sym_op.basis_map_calculated issue in ApplyAtom() for AAPL (CO)
        -- added xStream class for logging workflow updates (CO)
        -- fixed stability criterion vs. (Delta H) fonts in chull report (CO)
        -- substantial clean-up of classes in POCC in anticipation for AVASP_MakeSingleAFLOWIN integration (CO)
        -- added avasp function in preparation for AVASP_MakeSingleAFLOWIN integration (CO)
        -- added try/catches for easy debugging (CO)
        (AUROSTD: main.cpp)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, agl_debye.cpp, avasp.cpp, bader.cpp, chull.cpp, chull.h, contrib_kesong_ipocc.cpp, aflow.cpp, aflow.h, ivasp.cpp, kbin.cpp, kvasp.cpp, aflowlib_web_interface.cpp, pflow.h, pocc.cpp, pocc.h, symmetry.cpp, xatom.cpp, xclasses.cpp, APL/aapl_pairs.cpp, APL/apl_kphonons.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT)
3.1.189 - 2018/05/04 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.189.tar.xz
        Moved these defaults from aflow.h to aflow_aflowrc.cpp for user tuning (SC)
        #define DEFAULT_AFLOW_PRESCRIPT_OUT            string("aflow.prescript.out")  
        #define DEFAULT_AFLOW_PRESCRIPT_COMMAND        string("aflow.prescript.command")  
        #define DEFAULT_AFLOW_POSTSCRIPT_OUT           string("aflow.postscript.out")  
        #define DEFAULT_AFLOW_POSTSCRIPT_COMMAND       string("aflow.postscript.command") 
        #define DEFAULT_AFLOW_PGROUP_OUT               string("aflow.pgroup.out")
        #define DEFAULT_AFLOW_PGROUP_XTAL_OUT          string("aflow.pgroup_xtal.out")
        #define DEFAULT_AFLOW_PGROUPK_OUT              string("aflow.pgroupk.out")
        #define DEFAULT_AFLOW_PGROUPK_XTAL_OUT         string("aflow.pgroupk_xtal.out")
        #define DEFAULT_AFLOW_FGROUP_OUT               string("aflow.fgroup.out")
        #define DEFAULT_AFLOW_SGROUP_OUT               string("aflow.sgroup.out")
        #define DEFAULT_AFLOW_AGROUP_OUT               string("aflow.agroup.out")
        #define DEFAULT_AFLOW_IATOMS_OUT               string("aflow.iatoms.out")
        #define DEFAULT_AFLOW_PGROUP_JSON              string("aflow.pgroup.json")      
        #define DEFAULT_AFLOW_PGROUP_XTAL_JSON         string("aflow.pgroup_xtal.json") 
        #define DEFAULT_AFLOW_PGROUPK_JSON             string("aflow.pgroupk.json")    
        #define DEFAULT_AFLOW_PGROUPK_XTAL_JSON        string("aflow.pgroupk_xtal.json")
        #define DEFAULT_AFLOW_FGROUP_JSON              string("aflow.fgroup.json")   
        #define DEFAULT_AFLOW_SGROUP_JSON              string("aflow.sgroup.json")  
        #define DEFAULT_AFLOW_AGROUP_JSON              string("aflow.agroup.json")    
        #define DEFAULT_AFLOW_IATOMS_JSON              string("aflow.iatoms.json")   
        #define DEFAULT_AFLOW_ICAGES_OUT               string("aflow.icages.out")
        #define DEFAULT_AFLOW_SURFACE_OUT              string("aflow.surface.out")
        #define DEFAULT_AFLOW_QMVASP_OUT               string("aflow.qmvasp.out")
        #define DEFAULT_AFLOW_ERVASP_OUT               string("aflow.error.out")
        #define DEFAULT_AFLOW_IMMISCIBILITY_OUT        string("aflow.immiscibility.out")
        #define DEFAULT_AFLOW_MEMORY_OUT               string("aflow.memory.out")
        #define DEFAULT_AFLOW_FROZSL_INPUT_OUT         string("aflow.frozsl_input.out")
        #define DEFAULT_AFLOW_FROZSL_POSCAR_OUT        string("aflow.frozsl_poscar.out")
        #define DEFAULT_AFLOW_FROZSL_MODES_OUT         string("aflow.frozsl_energies.out")
        #define DEFAULT_AFLOW_FROZSL_EIGEN_OUT         string("aflow.frozsl_eigen.out")
        #define DEFAULT_AFLOW_END_OUT                  string("aflow.end.out")
3.1.188- 2018/05/03 -
       Location: http://materials.duke.edu/AFLOW/aflow.3.1.188.tar.xz
       Fixing kvasp.cpp (SC)
3.1.187- 2018/04/26 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.187.tar.xz
        -- added missing scaling factor information to lattice for AIMS xstructure output (DX) 
        -- fixed bug in edata lattice type/variation determination; accounts for lattices that do not reflect crystal symmetry (DX) 
        -- fixed bug in edata reciprocal lattice type/variation determination; accounts for reciprocal lattices transformed from a lattice that does not reflect crystal symmetry (DX) 
        -- added directory (pwd) information to LDEBUG/ERROR messages for symmetry functions (DX)
        (AFLOW: aconvasp_main.cpp, lattice.cpp, symmetry.cpp, xatom.cpp) 
3.1.186- 2018/04/25 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.186.tar.xz
        -- fixed bug in tolerance scan; was only scanning in one direction (DX) 
        -- created faster minimum cartesian distance calculator for skewed cells; fewer duplicate operations (DX) 
        (AFLOW: symmetry_spacegroup.cpp, symmetry.cpp) 
3.1.185- 2018/04/24 -
       Location: http://materials.duke.edu/AFLOW/aflow.3.1.185.tar.xz
       Fixing kvasp.cpp for dying jobs (SC)
3.1.184- 2018/04/23 -
       Location: http://materials.duke.edu/AFLOW/aflow.3.1.184.tar.xz
       Fixing ivasp.cpp for AFLOW_PSEUDOPOTENTIALS.TXT AFLOW_PSEUDOPOTENTIALS_LIST.TXT (SC)
       Tuning aflow_aflowrc.cpp for HYPERTHREADING in eos,draco,cobra,hydra (SC)
3.1.183- 2018/04/20 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.183.tar.xz
        Fixing ovasp for METAGGA/isKINETIC (SC)
3.1.182- 2018/04/19 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.182.tar.xz
        AFLOW_PSEUDOPOTENTIALS.TXT AFLOW_PSEUDOPOTENTIALS_LIST.TXT (SC)
3.1.181- 2018/04/17 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.181.tar.xz
        Creating defaults in aflow_aflowrc.cpp for HYPERTHREADING in eos,draco,cobra,hydra (SC)
3.1.180- 2018/04/16 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.180.tar.xz
        Creating defaults in aflow_aflowrc.cpp for NCPUS in eos,draco,cobra,hydra (SC)
3.1.179 - 2018/04/13 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.179.tar.xz
        -- fixed PP settings for SCAN (RF)
        -- fixed warning in APIget() (CO)
        -- added relative stability criterion and latex-formatted sg's to chull properties list (CO)
        -- added control.in and geometry.in file name specification in aflowrc (CO)
        -- --generate now applies for aims output as well (CO)
        -- force aflow.in generation even for MODE=AIMS (CO)
        -- fixed k-point mismatch for non-primitive APL runs (CO)
        -- added MINATOMS_RESTRICTED tag for APL (CO)
        -- added functionality for user-defined path in phonon dispersion, specify coords/labels in aflow.in (CO)
        -- swapped out quser for qflow, but maintain backwards compatibility (CO)
        (AUROSTD: xoption.cpp, xoption.cpp, xoption.h, xoption.h, xvector.cpp, xvector.cpp, xvector.h, xvector.h)
        (AFLOW: aconvasp_main.cpp, aconvasp_main.cpp, aflowrc.cpp, aflowrc.cpp, avasp.cpp, avasp.cpp, chull.cpp, chull.cpp, chull.h, chull.h, compare_structure.cpp, compare_structure.cpp, contrib_kesong_hnfcell.cpp, contrib_kesong_hnfcell.cpp, contrib_kesong_ipocc.cpp, contrib_kesong_ipocc.cpp, contrib_kesong_pocc_basic.cpp, contrib_kesong_pocc_basic.cpp, aflow.cpp, aflow.cpp, aflow.h, aflow.h, init.cpp, init.cpp, kbin.cpp, kbin.cpp, kvasp.cpp, kvasp.cpp, aflowlib.h, aflowlib.h, aflowlib_web_interface.cpp, aflowlib_web_interface.cpp, pocc.cpp, pocc.h, poccupation_forcefield.cpp, poccupation_forcefield.cpp, symmetry.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_functions.cpp, symmetry_spacegroup_functions.cpp, symmetry_spacegroup.h, symmetry_spacegroup.h, xatom.cpp, xatom.cpp, xclasses.cpp, xclasses.cpp, APL/aapl_tcond.cpp, APL/aapl_tcond.cpp, APL/aapl_tensor.cpp, APL/aapl_tensor.cpp, APL/apl_dirphoncalc.cpp, APL/apl_dirphoncalc.cpp, APL/apl.h, APL/apl.h, APL/apl_kphonons.cpp, APL/apl_kphonons.cpp, APL/apl_lrphoncalc.cpp, APL/apl_lrphoncalc.cpp, APL/apl_pathbuilder.cpp, APL/apl_pathbuilder.cpp, APL/apl_pdisc.cpp, APL/apl_pdisc.cpp, APL/apl_phoncalc.cpp, APL/apl_phoncalc.cpp, APL/apl_supercell.cpp, APL/apl_supercell.cpp, APL/qha_eosrun.cpp, APL/qha_eosrun.cpp, Makefile, Makefile, README_AFLOW_APL.TXT, README_AFLOW_APL.TXT, README_AFLOW.TXT, README_AFLOW.TXT)
3.1.178 - 2018/04/13 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.178.tar.xz
        Moved these defaults from aflow.h to aflow_aflowrc.cpp for user tuning (SC)
        #define DEFAULT_FILE_AFLOWLIB_ENTRY_OUT        string("aflowlib.out")
        #define DEFAULT_FILE_AFLOWLIB_ENTRY_JSON       string("aflowlib.json")
        #define DEFAULT_FILE_EDATA_ORIG_OUT            string("edata.orig.out")
        #define DEFAULT_FILE_EDATA_RELAX_OUT           string("edata.relax.out")
        #define DEFAULT_FILE_EDATA_BANDS_OUT           string("edata.bands.out")
        #define DEFAULT_FILE_DATA_ORIG_OUT             string("data.orig.out")
        #define DEFAULT_FILE_DATA_RELAX_OUT            string("data.relax.out")
        #define DEFAULT_FILE_DATA_BANDS_OUT            string("data.bands.out")
        #define DEFAULT_FILE_EDATA_ORIG_JSON           string("edata.orig.json")
        #define DEFAULT_FILE_EDATA_RELAX_JSON          string("edata.relax.json")
        #define DEFAULT_FILE_EDATA_BANDS_JSON          string("edata.bands.json")
        #define DEFAULT_FILE_DATA_ORIG_JSON            string("data.orig.json")
        #define DEFAULT_FILE_DATA_RELAX_JSON           string("data.relax.json")
        #define DEFAULT_FILE_DATA_BANDS_JSON           string("data.bands.json")
        #define DEFAULT_FILE_TIME_OUT                  string("time")
        #define DEFAULT_FILE_SPACEGROUP1_OUT           string("SpaceGroup")
        #define DEFAULT_FILE_SPACEGROUP2_OUT           string("SpaceGroup2")
        #define DEFAULT_FILE_VOLDISTPARAMS_OUT         string("VOLDISTParams")
        #define DEFAULT_FILE_VOLDISTEVOLUTION_OUT      string("VOLDISTEvolution")
3.1.177 - 2018/04/06 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.177.tar.xz
        Bug fixes on autopseudootentialsl (SC)
3.1.176 - 2018/04/06 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.176.tar.xz
        Working PAW_PBE_KIN and PAW_LDA_KIN autopseudopotential (SC)
3.1.175 - 2018/04/06 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.175.tar.xz
        Working on mpcdf-cobra (SC)
3.1.174 - 2018/03/29
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.174.tar.xz    
        -- added bands app to entry page (GG)
        -- fixed scaling factor type in xstructure2json() (DX)
        -- account for Wyckoff positions in represented as fractions (hex/rhl) in minimum enumerated Wyckoff search (DX)
        -- improved minimum enumerated Wyckoff search (consider combinations of origin shifts) (DX)
        -- fixed typo in xmatrix2json function (DX)
        -- limited outlier detection to binaries only, statistics not globally favorable for ternaries and above YET (CO)
        -- added explanatory comments in aflowrc for chull settings (CO)
        -- added statistics output to chull logger (CO)
        -- fixed Greek letter issues with Helvetica fonts in chull (CO)
        -- further decoupled ChullFacet() initialization from ConvexHull() with initialize() and addVertex() (CO)
        -- added user-defined dft_type restrictions in aflowrc (CO)
        -- added LIB1 to loadEntries() default (CO)
        -- fixed stringElements2VectorElements() bugs with LIB1 colons (CO)
        -- added plotting option for iso-max latent heat (CO)
        -- added kJ/mol axis for formation enthalpy hulls (CO)
        -- added logos to hull illustrations, generally fixed header/footers (CO)
        -- added stability criterion analysis to default routine (CO)
        -- fixed sign of decomposition energy/stability criterion (CO)
        -- added equivalent ICSD structures analysis to default routine (CO)
        -- fixed --readme vs. --readme=XX (CO)
        -- fixed superfluous output in --poscar2aflowin (CO)
        -- fixed bug in aconvasp's Ewald summation function as pointed out by Wei Xie (CO)
        -- changed generic xstructure title to include cleannames in case there is pp info (CO)
        (AUROSTD: AUROSTD/aurostd.h, main.cpp, xscalar.h, xvector.cpp)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, chull.cpp, chull.h, init.cpp, aflowlib.h, aflowlib_web_interface.cpp, pflow_funcs.cpp, pflow.h, symmetry_spacegroup.cpp, webapp_bands.js, webapp_entry.js, xatom.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT)
3.1.173 - 2018/03/12
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.173.tar.xz    
        -- Create LIB7/LIB8/LIB9 framework (SC)
3.1.172 - 2018/02/27
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.172.tar.xz    
        -- Write data in JSON format in AEL/AGL (CT)
        -- Option to turn off VASP symmetry in AEL/AGL (CT)
        -- Write vibrational free energy and vibrational entropy at 300K in aflow.agl.out in AEL/AGL (CT)
        -- Write enthalpy (H = E + pV) in AEL/AGL (CT)
        -- Option to skip pressure/temperature points where no minimum Gibbs free energy is found, instead of truncating pressure/temperature range in AEL/AGL (CT)
        -- fixed tolerance scan counter for lib2raw runs (no longer static variable) (DX)
        -- edata speed increase, only calculate up to pgroup_xtal for lattice, superlattice, and reciprocal lattice (DX)
        -- fixed typo in Wyckoff position "b" for space group #160 (DX)
        -- fixed --kppra command line bug, divided by natoms twice (CO)
        -- added debug output to nanoparticle (CO)
        -- fixed --readme=chull empty string error (CO)
        -- automated plot_unstable with z_filter_cutoff in chull (CO)
        -- fixed plotting ranges with z_filter_cutoff in chull (CO)
        (AFLOW: aconvasp_main.cpp, ael_elastic_fit.cpp, ael_elasticity.cpp, ael_elasticity.h, ael_get_stress.cpp, agl_debye.cpp, agl_debye.h, agl_eqn_state.cpp, agl_get_ev.cpp, agl_hugoniot.cpp, agl_polynomial.cpp, agl_rungibbs.cpp, chull.cpp, data.cpp, aflow.h, lattice.cpp, aflowlib_libraries.cpp, pflow.h, pflow_print.cpp, symmetry.cpp, symmetry_spacegroup.cpp, symmetry_spacegroup_ITC_library.cpp, xatom.cpp, Makefile, README_AFLOW_AEL.TXT, README_AFLOW_AGL.TXT, README_AFLOW_SYM.TXT)
3.1.171 - 2018/02/21
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.171.tar.xz    
        -- fixed multiple degeneracy prints in POCC (CO)
        -- added directory information to logging functions in bader/chull (CO)
        -- added --destination=|--path= flag for output of chull (CO)
        -- use full path (pwd) for working directory in AFLOW-SYM functions (DX)
        -- fixed directory flag for AFLOW-SYM functions (DX)
        -- fixed typo in site symmetry of  Wyckoff position "a" for space group #109 (DX)
        (AUROSTD: aurostd.h, main.cpp)
        (AFLOW: aconvasp_main.cpp, bader.cpp, bader.h, chull.cpp, chull.h, poccupation_edos.cpp, symmetry_spacegroup_ITC_library.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT)
3.1.170 - 2018/02/18
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.170.tar.xz    
        -- fixed chmod in LIB2RAW (CO)
        -- fixed DOS-extraction bug (exit() vs. return false) in LIB2RAW for POCC+AEL/AGL runs (CO)
        -- fixed joinWDelimiter() xvector<int> bug: removed erroneous delimiter at the end (CO)
        -- fixed DOSDATA_JSON() + BANDSDATA_JSON() amalgamation: removed erroneous wrapping brackets around DOS (CO)
        -- fixed misleading logging message in APL on primitivization of input (CO)
        (AUROSTD: main.cpp)
        (AFLOW: estructure.cpp, aflow.h, aflowlib_libraries.cpp, aflowlib_web_interface.cpp, poccupation_edos.cpp, APL/apl_supercell.cpp, Makefile)
3.1.169 - 2018/02/16
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.169.tar.xz    
        -- added aflow_proto functionality:  --relax_type=IONS, --module=APL, --apl_supercell=3x3x3, --no_volume_adjustment (CO)
        -- NCPUS in APL now respects parent aflow.in (default MAX, then looks at parent aflow.in, otherwise overrides with --np=XX) (CO)
        -- fixed spacegroup bug in apl post hibernation (CO)
        -- added --print and --screen_only options to chull for python integration (CO)
        -- removed extraneous ytick lines from 3D hull (CO)
        -- increased spacing between axes and ticklabels on hulls (CO)
        -- fixed bugs in entropic_temperature hull visualization (CO)
        -- fixed apl2agr/subst make commands (CO)
        -- fixed Niggli tolerance; more robust and use _ZERO_TOL_ (DX)
        -- fixed bug in Niggl step 6 (DX)
        -- fixed comparison of SU(2) to exp(theta*su(2)); parentheses issue (DX)
        -- fixed Quantum Espresso output to include lattice scaling factor (DX)
        -- fixed Quantum Espresso celldm units (Bohr) (DX)
        -- fixed Quantum Espresso alat flag to multiply cell parameters by celldm(1) or A parameter (DX)
        -- fixed tolerance issue with PrintData for lib2raw runs (DX)
        (AFLOW: aconvasp_main.cpp, avasp.cpp, chull.cpp, aflow.h, pflow_print.cpp, symmetry.cpp, xatom.cpp, APL/aapl_tcond.cpp, APL/apl_dm_pdos_save.cpp, APL/apl_doscalc.cpp, APL/apl.h, APL/apl_pdisc.cpp, APL/apl_phoncalc.cpp, APL/qha_gruneisen.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT)
3.1.168 - 2018/02/08
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.168.tar.xz
        -- fixed structure comparison bug (overwritten match)  (DX)
        (AFLOW: aflow_compare_structure.h, aflow_compare_structure_function.cpp)
3.1.167 - 2018/02/06
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.167.tar.xz
        -- speed increase for structure comparison (DX)
        -- fixed bugs for structure comparison (DX)
        -- added Wyckoff position analysis to group similar structures in directory comparison (DX)
        -- added grouped Wyckoff position information to the JSON output (DX)
        -- added logger/status for directory comparisons (DX)
        -- added functionality to read .bz2 geometry files (DX)
        -- added --no_scale option for comparisons (required for chull) (CO/DX)
        (AFLOW: aconvasp_main.cpp, chull.cpp, symmetry_spacegroup.cpp, compare_structure.cpp, compare_structure_function.cpp, compare_structure.h) (DX)
3.1.166 - 2018/02/02
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.165.tar.xz
        -- fixed broken --slab commmand (CO)
        (AFLOW: slab.cpp)
3.1.165 - 2018/01/30
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.165.tar.xz
        -- integrated .aflow.rc settings into APL/QHA/AEL/AGL (CO)
        -- added LORBIT=10 to relaxation INCAR (spinD) (CO)
        -- improved magnetic properties extraction, pull from relax first, overwrite with static (CO)
        -- added safety hull coordinates assignment (CO)
        -- new default for load entries (no load xstructures) (CO)
        (AFLOW: aconvasp_main.cpp, ael_get_stress.cpp, agl_get_ev.cpp, chull.cpp, ivasp.cpp, aflowlib_libraries.cpp, APL/apl_phoncalc.cpp, APL/qha_eosrun.cpp, Makefile) (CO)
3.1.164 - 2018/01/25
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.164.tar.xz
        -- fixed minor bugs in Quantum Espresso reader
        -- added ibrav options to Quantum Espresso reader (ibrav2lattice function)
        -- added celldm and a, b, c, cosAB, cosAC, cosBC readers for Quantum Espresso
        -- added Bohr units reader for Quantum Espresso
        (AFLOW: aconvasp_main.cpp, pflow.h, xatom.cpp) (DX)
3.1.163 - 2018/01/23
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.163.tar.xz
        -- added --dist2hull=0.25,0.25 option for chull, which provides the value of the convex hull surface at the specified coordinate/concentration (CO)
        -- added separate readme for chull (--readme=chull) (CO)
        -- customized avasp.cpp (--proto) with options from .aflow.rc, --mpi, and --np (CO)
        -- fixed QHA/AAPL aflow.in keyword bug (CO)
        -- added ANRL modifiers to directory for --proto (Ex: label:ANRL=param1,param2,...) (DX)
        (AFLOW: aconvasp_main.cpp, aflowrc.cpp, avasp.cpp, chull.cpp, aflow.cpp, aflow.h, init.cpp, kbin.cpp, Makefile, README_AFLOW_ACONVASP.TXT, README_AFLOW_CHULL.TXT, README_AFLOW.TXT) (CO/DX)
3.1.162 - 2018/01/19
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.162.tar.xz
        -- added SU(2) complex matrix and su(2) generator coefficients representation of symmetry elements 
        -- added SU(2) and su(2) information to .out and .json symmetry files        
        -- extended functionality for complex matrices and vectors (trace, exponential, ostream, etc.)
        -- added xcomplex2json to represent complex numbers in json
        -- added pgroupk_xtal (dual of crystal point group operations) standalone function
        -- extended --proto for ANRL prototypes (to generate aflow.in)
        -- updated README_AFLOW_SYM.TXT
        (AUROSTD: aurostd.h, boot.cpp, main.cpp, xmatrix.cpp, xmatrix.h)
        (AFLOW: README_AFLOW_ACONVASP.TXT, README_AFLOW_SYM.TXT, aflow.h, aconvasp_main.cpp, avasp.cpp, symmetry.cpp, xatom.cpp) (DX)
3.1.161 - 2018/01/17
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.161.tar.xz
        -- bool xscalar and xvector gccV7 warning fixes (SC)
        -- double xscalar and xmatrix clang warning fixes (CO)
        -- uint chull gccV7 warning fix (CO)
        (AUROSTD: xscalar.h, xvector.h, xmatrix.h)
        (AFLOW: chull.cpp)
        Added README_AFLOW_SYM.TXT (init.cpp aflow.h aflow.cpp data.cpp) (DX/SC)
3.1.160 - 2018/01/13
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.160.tar.xz
        convex hull mods:
          - keep order of points (m_points) same as input from user/database
          - size will only differ if there are artificial points (they go at the end of m_points)
          - moved point banning to structurePoints() (i.e., m_coord_groups will only contain validated points)
          - only print decomposition phases with non-zero coefficients (tie-line)
          - fixed bug in initial facet.f_outside_set calculation
          - fixed bug in using hull centroid vs. hull reference to align normals
          - fixed bug in vcoords input (shiftlrows)
        (AFLOW: chull.cpp, chull.h)
3.1.159 - 2018/01/10 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.159.tar.xz
        MAJOR UPDATE - convex hull code rewritten, includes: (CO)
          - complete generalization of input (need not be thermodynamic hull, can be ANY data)
          - new major classes for naries, alloys, and coordinate (stoichiometry) groups
          - for hulls with stoichiometric coordinates, the hull and properties are calculated in steps of increasing dimensionality
          - automatic calculation of equivalent ground state structures, near-equivalent (by symmetry) ground state structures, equilibrium phases (by mixture)
          - automatic removal of equivalent ground states for calculation of stability criterion
          - automatic detection/removal of outliers (IQR/MAD)
          - moved all major plotting options to .aflow.rc
          - report any bugs to aflow.org/forum
        added xcombo class for combinations/permutations (CO)
        integrated roundoff/precision/format into stream2stream/utype2tring and variants (xvecDouble2String(),xmatDouble2String()) (CO)
        added another date format for chull (CO)
        added tol to isinteger (CO)
        added +=, -=, gcd (reduced composition), shiftlrows(), isCollinear(), getCentroid(), getGeneralAngles() (similar to Euler angles), getGeneralNormal() (n-dimensional normal), and getQuartiles()/getMAD() (outlier analysis) for xvector (CO)
        added general Householder QR decomposition (CO)
        smarter overloading of loadEntries()/loadXstructures() (CO)
        updated logger() to take directory (aflags) (CO)
        fixed Wyckoff position (a) for space group 117 (DX)
        fixed monoclinic space group determination with multiple unique axis choices (DX)
        (AUROSTD: Makefile,aurostd.cpp,aurostd.h,boot.cpp,main.cpp,xcombos.cpp,xcombos.h,xmatrix.cpp,xmatrix.h,xoption.cpp,xoption.h,xscalar.cpp,xscalar.h,xvector.cpp,xvector.h) (CO, patched SC)
	(AFLOW: Makefile,aflow.h,aconvasp_main.cpp,aflowrc.cpp,anrl.cpp,chull.cpp,chull.h,compare_structure.cpp,compare_structure.h,compare_structure_function.cpp,pflow.h,pflow_funcs.cpp,xatom.cpp) (CO)
	(AFLOWLIB: aflowlib.h,aflowlib_libraries.cpp,aflowlib_web_interface.cpp) (CO)
	(SYM: symmetry_spacegroup.cpp,symmetry_spacegroup_ITC_library.cpp) (CO)
3.1.158 - 2018/01/09 -
        Location: http://materials.duke.edu/AFLOW/aflow.3.1.158.tar.xz
        Tentative distribution GNU (SC)
3.1.157 - 2017/12/19 -
	Working on mpcdf-draco (SC) 
        Replaced pgroupk with pgroupk_xtal for resolving IBZ in APL/AAPL (CO)
        Fixed pgroupk_xtal flag and wyckoff position for sg133 (DX)
        (aflowrc.cpp,APL/apl_supercell.cpp,APL/apl_mpmesh.cpp,APL/aapl_tcond.cpp,lattice.cpp,symmetry_spacegroup_ITC_library.cpp)
3.1.156 - 2017/12/18 -
        Working on mpcdf-hydra (SC)
3.1.155 - 2017/12/13 -
        updated web component of ANRL CPP files to include the prototype generator and Jmol visualization
        fixed Wyckoff position for SG #171
        improvement for monoclinic space group determination (consider alternative unique axis choices, if possible)
        added more to debug output for the space group routine
        (ANRL/,ANRL_CPPS_20171213/,aconvasp_main.cpp,symmetry_spacegroup.cpp,symmetry_spacegroup_ITC_library.cpp,symmetry_spacegroup_functions.cpp,xatom.cpp) (DX)
      
3.1.154 - 2017/12/07 -
        added pgroupk_xtal (dual of crystal point group operations)
        added magnetic symmetry analysis (crystal + spin) to edata 
        added non-collinear magnetic symmetry analysis (crystal + spin) to 
        edata,sgdata,aflowSG,aflowSYM,fgroup,pgroup_xtal,pgroupk_xtal,agroup,and sgroup
        fixed ApplyAtom (tolerance issue)
        fixed AFLOW-SYM printing bug (--screen_only for Python environment)
        changed point group matrix comparisons (uses Uf; exact)
        added SO(3) generator expansion coefficients onto Lx, Ly, and Lz
        (README_ACONVASP.TXT,aflow.h,aconvasp_main.cpp,ovasp.cpp,pflow.h,symmetry.cpp,symmetry_spacegroup.cpp,symmetry_spacegroup_functions.cpp,xatom.cpp,lattice.cpp,aflowlib_libraries.cpp) (DX)       
 
3.1.153 - 2017/11/23 -
        Fixing Makefile for GNU version (SC) (fixed Library_ICSD..) (SC)
	Updating directories from AFLOW2 to AFLOW3 (SC)
	Updating directories for findsym platon frozsl (SC)
	
3.1.152 - 2017/11/20 -
        Fixing Makefile for GNU version (SC) (fixed DATA_CVs..)
3.1.151 - 2017/11/17 -
        Fixing Makefile for GNU version (SC) (fixed DATA_CVs..)
3.1.150 - 2017/11/06 -
        NFS cache-cleaning HACK
        (kbin.cpp) (CO)
3.1.149 - 2017/11/06 -
        fixed --xplug
        (libraries.cpp) (SC)
3.1.148 - 2017/10/27 -
        fixed --edata (txt + json print outs)
        (aconvasp_main.cpp,pflow.h,pflow_print.cpp,aflowlib_libraries.cpp) (CO/DX)
3.1.147 - 2017/10/25 -
	added skew test for structures (sym_eps),
	incorporated full symmetry analysis in more functions,
	improved --lib2raw_local,
	improved dos/bands-2json functionality,
	actually fixed --delta_kpoints,
	made min_dist (nearest neighbor distance) function more robust/faster,
	fixed nbondxx units,
	added jsons for sym/structure/bandsdata/dosdata to lib2raw,
	incorporated scale (ReScale()) in more functions
	(aurostd.h,boot.cpp,main.cpp,xmatrix.cpp,xmatrix.h,Makefile,README_AFLOW_ACONVASP.TXT,aflow.h,aconvasp_main.cpp,bader.cpp,defects.cpp,estructure.cpp,init.cpp,ovasp.cpp,pflow.h,pflow_print.cpp,surface.cpp,symmetry.cpp,symmetry_spacegroup.cpp,xatom.cpp,xclasses.cpp,xproto.cpp,aflowlib.h,aflowlib_libraries.cpp) (CO/DX)
3.1.146 - 2017/10/16 -
	rewritten bandgap code,
	added extensive FHI-AIMS support (APL),
	modified loadEntries() for GFA code,
	added new xvector operations (*=,/=),
	fixed delta_kpoints function,
	fixed bands/dos-2json functions,
	added local lib2raw
	(aapl_tensor.cpp,apl.h,apl_dirphoncalc.cpp,apl_gsa.cpp,apl_kphonons.cpp,apl_lrphoncalc.cpp,apl_phoncalc.cpp,qha_eosrun.cpp,aurostd.h,boot.cpp,main.cpp,xmatrix.cpp,xvector.cpp,xvector.h,xmatrix.h,Makefile,README_AFLOW_ACONVASP.TXT,aflow.h,aconvasp_main.cpp,aflowrc.cpp,chull.cpp,kbin.cpp,ovasp.cpp,pflow.h,poccupation_edos.cpp,xatom.cpp,xclasses.cpp,aflowlib_libraries.cpp,apl_supercell.cpp,iaims.cpp,kaims.cpp,oaims.cpp,estructure.cpp,ivasp.cpp,aflowlib.h) (CO)
3.1.145 - 2017/10/12 -
	fixed tolerance scan range (symmetry.cpp) (DX) 
3.1.144 - 2017/10/03 -
	fixed bandgap type bug,
	fixed precision mistmatch in bands data when plotting band structure,
	changed misleading metric_tensor function name,
	fixed false positive message in structure_comparison,
	fixed APL bugs with reading aflow.in,
	avoid SYMPREC vs. PREC confusion
	(ovasp.cpp,estructure.cpp,xatom.cpp,aconvasp_main.cpp,pflow.h,aflow.h,aurostd.h,aflowlib_libraries.cpp,poccupation_edos.cpp,Makefile,ivasp.cpp,apl_phoncalc.cpp,qha_eosrun.cpp,apl.h,apl_kphonons.cpp,compare_structure.cpp,symmetry_spacegroup_functions.cpp,symmetry_spacegroup.h) (CO/DX)
3.1.143 - 2017/09/28 -
	added magnetic symmetry analysis (crystal + spin) to sgdata,
	aflowSG,aflowSYM,agroup,fgroup,pgroup_xtal,agroup,and sgroup; speed increase for space group routine,
	alter tolerance scan range
	(README_AFLOW_ACONVASP.TXT,aflow.h,aconvasp.cpp,aconvasp_main.cpp,lattice.cpp,pflow.h,symmetry.cpp,symmetry_spacegroup.cpp,symmetry_spacegroup_functions.cpp,xatom.cpp) (DX)  
3.1.142 - 2017/09/14 -
	added mpcdf_eos machine,
	fixed proper nouns in symmetry output
	(aflow.cpp,kbin.cpp,aflowrc.cpp,init.cpp,kvasp.cpp,README_AFLOW.TXT,Makefile,pflow_print.cpp) (CO/DX)
3.1.141 - 2017/09/12 -
	fixed sgdata/wyccar issue,
	added space group functions to SYM namespace,
	more detailed warning messages
	(symmetry_spacegroup.h,symmetry_spacegroup.cpp,symmetry_spacegroup_function.cpp,symmetry_spacegroup_ITC_library.cpp,aconvasp_main.cpp,pflow_print.cpp,symmetry.cpp) (DX) 
3.1.140 - 2017/09/05 -
	added sgdata function (along with options for json output),
	edata updated (along with options for json output),
	pgroupk speed up,
	added covalent radii info,
	dos/band structure gnuplot updates
	(main.cpp,README_AFLOW_ACONVASP.TXT,aflow.h,aconvasp_main.cpp,lattice.cpp,pflow.h,pflow_print.cpp,symmetry.cpp,symmetry_spacegroup.cpp,symmetry_spacegroup.h,symmetry_spacegroup_functions.cpp,xatom.cpp,estructure.cpp,init.cpp) (DX/CO) 
3.1.139 - 2017/08/29 -
	fixed edata issue (agroups/iatoms were not calculated) for library runs
	(README_AFLOW_ACONVASP.TXT,aconvasp_main.cpp,lattice.cpp,xatom.cpp,aflow.h) (DX)
3.1.138 - 2017/08/25 -
	fully functional AAPL,
	Pearson symbol/lattice type speed up,
	json serializers for symmetry groups,
	tolerance scan implemented for separate symmetry groups,
	chull latex fix
	(chull.cpp,aconvasp_main.cpp,aflow.h,lattice.cpp,pflow.h,symmetry.cpp,xatom.cpp,xclasses.cpp,aapl_pairs.cpp,aapl_tcond.cpp,aapl_tensor.cpp,apl.h,apl_kphonons.cpp,aurostd.h,main.cpp) (JJPR/DX/CO)
3.1.137 - 2017/08/04 -
	incorporated aflow.Xgroup.json,
	json serializers for bands/dos,
	print to screen options for AFLOW-SYM API,
	xvector/xmatrix lrows/urows/lcols/ucols bug fix,
	APL symmetry output append to ofstream,
	removed validation for sgroup (NOT NEEDED)
	(aconvasp_main.cpp,aflow.h,pflow.h,symmetry.cpp,xatom.cpp,README_AFLOW_ACONVASP.TXT,main.cpp,aurostd.h,xvector.cpp,xmatrix.cpp,apl_supercell.cpp,Makefile,estructure.cpp) (DX/EG/CO)
3.1.136 - 2017/08/01 -
	speed up getFullSymBasis,
	fixed some Wyckoff positions,
	fixed printing in apl (spacegroup vs. lattice type)
	(apl_pdisc.cpp,apl_phoncalc.cpp,aflow.h,symmetry.cpp,symmetry_spacegroup.cpp,symmetry_spacegroup.h,symmetry_spacegroup_ITC_library.cpp,symmetry_spacegroup_functions.cpp) (DX/CO)
3.1.135 - 2017/07/25 -
	fixed xOUTCAR parsing issue with lattice vectors (vasp bug)
	(aflow.h,ovasp.cpp) (CO)
3.1.134 - 2017/07/19 -
	fixed symmetry for library runs (angle/vol tolerance),
	added directory output,
	removed relaxed data from orig structure section in entry page
	(aflowlib_web_interface.cpp,Makefile,aconvasp_main.cpp,aflow.h,lattice.cpp,aflowlib_libraries.cpp,pflow_print.cpp,symmetry.cpp,symmetry_spacegroup.cpp,xatom.cpp,xproto.cpp) (DX/GG/CO)
3.1.133 - 2017/07/19 -
	AEL/AGL: Hugoniots,
	AEL vs pressure (CT) 
3.1.132 - 2017/07/14 -
	added full_sym output to RAW
	(aflowlib_libraries.cpp) (DX/CO)
3.1.131 - 2017/07/14 -
	allow global aflowrc in /etc/aflow.conf plus bug fixes in aflowrc (SC)
3.1.130 - 2017/07/13 -
	bug fixes in aflowrc (SC)
3.1.129 - 2017/07/13 -
	bug fixes init.cpp aflowrc.cpp xvector.cpp (SC/CO)
3.1.128 - 2017/07/12 -
	more ~/.aflowrc definitions and porting to MPI-DE (SC)
3.1.127 - 2017/07/11 -
	more ~/.aflowrc aflow.h ivasp ovasp avasp kbin init VASP_POTCAR_* AFLOW_LIBRARY AFLOW_PROJECT (SC)
3.1.126 - 2017/07/10 -
	more ~/.aflowrc name_dynamics aflow.h  DEFAULT_DOS estructure.cpp poccupation_edos (SC)
3.1.125 - 2017/07/10 -
	more ~/.aflowrc aflow.h  kbin.cpp apl_phoncalc.cpp apl_kphonons.cpp (SC)
3.1.124 - 2017/07/06 -
	quaternion representation additions,
	aflowSYM bug fixes (HfV2)
	(lattice.cpp,xatom.cpp,aflow.h,symmetry.cpp,Makefile) (GG/DX/CO)
3.1.123 - 2017/07/06 -
	adding flags for METAGGA and IVDW (xclasses.cpp) (SC)
3.1.122 - 2017/07/04 -
	aflowrc bug fixes,
	clean names (SC)
3.1.121 - 2017/07/04 -
	aflowrc bug fixes,
	clean names,
	clean METAGGA_SET,
	IVDW_SET (SC)
3.1.120 - 2017/07/02 -
	more ~/.aflowrc
	aflow.h init.cpp aflowrc.cpp kbin kvasp README* VASP_FORCE* PLATON/FINDSYM (SC)
3.1.119 - 2017/07/02 -
	more ~/.aflowrc
	aflow.h init.cpp aflowrc.cpp kbin kvasp README* MACHINE DEPENDENT MPI* (SC)
3.1.118 - 2017/07/02 -
	more ~/.aflowrc
	aflow.h init.cpp aflowrc.cpp kbin kvasp README* GENERIC MPI* (SC)
3.1.117 - 2017/07/02 -
	more ~/.aflowrc
	aflow.h init.cpp aflowrc.cpp README*.TXT: AFLOW_CORE_TEMPERATURE,VASP_PREC (SC)
3.1.116 - 2017/07/01 -
	starting ~/.aflowrc
	aflow.h init.cpp aflowrc.cpp xoptions.h.cpp avasp.cpp apennsy lattice README_AFLOW_ACONVASP.TXT (SC)
3.1.115 - 2017/06/29 -
	XVASP_INCAR_PREPARE_GENERIC aflow.h ivasp.cpp (SC)
3.1.114 - 2017/06/29 -
	added new Jmol js scripts,
	fixed DISMAG issue with APL,
	checkMEMORY() mostly removed from messaging
	(apl_dirphoncalc.cpp,apl_kphonons.cpp,aflowlib_web_interface.cpp,entry.cpp,entry.js,kvasp.cpp,kbin.cpp,init.cpp,avasp.cpp,aflow.h,README_AFLOW_POCC.TXT,README_AFLOW_APL.TXT,README_AFLOW_AGL.TXT,README_AFLOW_AEL.TXT,README_AFLOW.TXT,Makefile) (CO)
3.1.113 - 2017/06/28 -
	XVASP_INCAR_PREPARE_GENERIC aflow.h ivasp.cpp (SC)
3.1.112 - 2017/06/28 -
	added EDIFFG --ediffg
	aflow.h avasp aconvasp_main ivasp kvasp xclasses README_AFLOW_CONVASP README_AFLOW (SC)
3.1.111 - 2017/06/13 -
	SQLITE integration
	(avasp.cpp ow_chull,aflowlib_web_interface,aflowlib_libraries,kbin,apl_phoncalc,README_AFLOW_VERSIONS_HISTORY.TXT,pflow_print,pflow_funcs,main,aurostd.h,bader,Makefile,aflow,SQLITE) (CO)
3.1.110 - 2017/06/13 -
	fixed chull invalid inputs,
	removed exit(0) in loadLIBX,
	force output for mpi settings in APL,
	fixed settings in CIF writer,
	removed aflow.end.out from CompressDirectory(),
	major fixes to CHGCAR/AECCAR readers for bader analysis (CO)
3.1.109 - 2017/06/06 -
	fixed some chull issues (CO)
3.1.108 - 2017/06/06 -
	individualized QHA/AAPL workflows,
	set reliable convex hull standard to 200 entries (binaries),
	updated ill-calculated systems in the database,
	incorporated sort for loadEntries(),
	fixed DPM=AUTO in APL,
	fixed MCL systems
	(APL/aapl_tensor,apl.h,apl_kphonons,apl_lrphoncalc,apl_phoncalc,qha_eosrun,apl_dirphoncalc,Makefile,README_AFLOW.TXT,README_AFLOW_ACONVASP.TXT,README_AFLOW_AEL.TXT,README_AFLOW_AGL.TXT,README_AFLOW_
APL.TXT,README_AFLOW_POCC.TXT,README_AFLOW_VERSIONS_HISTORY.TXT,aflow,aflow.h,avasp,init,kbin,kvasp,xclasses,aconvasp_main,chull,chull.h,symmetry_spacegroup,symmetry) (CO)
3.1.107 - 2017/06/05 -
	online pubs management curtarolo/oses/toher (SC)
3.1.106 - 2017/06/02 -
	added IVDW
	(aflow.h avasp xclasses aconvasp_main ivasp kvasp) (SC)
3.1.105 - 2017/06/01 -
	added METAGGA
	(aflow.h avasp xclasses aconvasp_main ivasp kvasp) (SC)
3.1.104 - 2017/05/30 -
	set reliable convex hull standard to 200 entries (CO)
3.1.103 - 2017/05/30 -
	added pgroup_xtal propagation in supercell creation,
	also fixed ".nfs*" issues in CompressDirectory() (CO)
3.1.102 - 2017/05/29 -
	issues with symmetry in the aflowlib.out calculation (CO)
3.1.101 - 2017/05/26 -
	cosmetics and added PaddedCENTER (SC)
3.1.100 - 2017/05/26 -
	changed numbering style + cosmetics(SC)
