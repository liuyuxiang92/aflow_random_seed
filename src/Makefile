# ***************************************************************************
# *                                                                         *
# *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
# *                                                                         *
# ***************************************************************************
#
#  Copyright 2003-2019 - Stefano Curtarolo - AFLOW.ORG consortium
#
#  This file is part of AFLOW software.
#
#  AFLOW is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ***************************************************************************

VERSNUMBER=3.1.227

# MAKEFILE FOR AFLOW3
# FOR FROZSL -lmkl_lapack95_ilp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential

#system variables
UNAME=$(shell uname)

# for aflow
# https://stackoverflow.com/questions/714100/os-detecting-makefile
ifeq ($(UNAME),Darwin)
CPP=clang++
CC=clang
ARCH=-e ARCH=-D_MACOSX_
else
CPP=g++
#-O2 -frepo
# for sqllite
CC=gcc
ARCH=
endif

#CPP=g++-fsf-7 
#CPP=g++-fsf-4.9 -std=gnu++11
# -std=gnu++11
# -O3
#CPP=g++ -std=c++0x-4.7 -std=c++0x
# -std=gnu++0x
#CCFLAGS= -DCYGWIN -Wall -W -fno-inline
CCFLAGS=-Wall -W -fno-inline
CCFLAGS_MT=$(CCFLAGS) -O2
#CCFLAGS_GSA=$(CCFLAGS) -O1	#CO OBSOLETE 170601
INCLUDE = -I ./
#INCLUDE = -I ./  -I /opt/intel/mkl100/include

# -fno-inline -fno-implicit-inline-templates-frepo
OPTS=
#-O1
#-O3
# -march=opteron
TODAY=-DTODAY=\"$(TIME)\"
VERS=-DAFLOW_VERSION=\"$(VERSNUMBER)\" $(TODAY)

# -DHOSTNAME=\"`hostname`\" -DYYY=\"2003-2012\"

#BIN_MARYLOU=/fslhome/fslcollab8/group/bin
BIN_MARYLOU=/fslhome/fslcollab8/bin
#BIN_RANGER=/share/home/00470/tg457283/bin
BIN=~/bin
BIN_KRAKEN=/nics/a/proj/aflow/bin
BIN_OHAD=/home/aflow/bin

EXTRA_FREE=../AFLOW3_FREE/EXTRA
ICSD_LIB=$(EXTRA_FREE)/ICSD
NIST_LIB=/common/NIST/2014
AFLOW_LIB=/common/AFLOW/LIBS
#NIST_LIB=/common/NIST/2010
PROJECT_DIRECTORY_GNDSTATE=/common/LIB2

# -static
OPTSA=-O3
# -static
DAY=$(shell date +%y)$(shell date +%m)$(shell date +%d)
#TIME="20$(shell date +%y)-$(shell date +%m)-$(shell date +%d) $(shell date +%H):$(shell date +%M):$(shell date +%S)"
TIME="20$(shell date +%y)-$(shell date +%m)-$(shell date +%d)"
LIBS=-lpthread
# -static
#CO START 170614 easy integration
LIBS_AFLOW=-ldl $(LIBS)
LIBS_STATIC=-static
LIBS_AFLOW_STATIC=$(LIBS_AFLOW) $(LIBS_STATIC)
#CO END 170614 easy integration
# -lglut -lGLU -lGL -D_HOST_=\"`hostname`\"
#AFLOW_VERSION := $(shell cat aflow.cpp | /home/auro/bin/getvalue AFLOW_VERSION)-$(TIME)
#AFLOW_VERSIONM := $(shell cat Makefile | /home/auro/bin/getvalue MAKEFILE)-$(TIME)
#DEBUG=-DXSTR_DEBUG
DEBUG=
BACKUP=/common/AFLOW/SRC/

# Library MKL dynamic link
#LIB_MKL= -L /opt/intel/mkl100/lib/32 -lmkl_lapack -lmkl_core -lmkl_intel -lmkl_sequential
# Library MKL static link
#LIB_MKL= -L /opt/intel/mkl100/lib/32/ /opt/intel/mkl100/lib/32/libmkl_lapack.a \
# /opt/intel/mkl100/lib/32/libmkl_core.a \
# /opt/intel/mkl100/lib/32/libguide.a

AUROSTD_CPP = ./AUROSTD/aurostd.cpp ./AUROSTD/aurostd_argv.cpp ./AUROSTD/aurostd_xoption.cpp ./AUROSTD/aurostd_xscalar.cpp ./AUROSTD/aurostd_xvector.cpp ./AUROSTD/aurostd_xmatrix.cpp ./AUROSTD/aurostd_xcomplex.cpp ./AUROSTD/aurostd_xtensor.cpp ./AUROSTD/aurostd_xcombos.cpp ./AUROSTD/aurostd_xerror.cpp ./AUROSTD/aurostd_xrandom.cpp ./AUROSTD/aurostd_xtensor_template.cpp ./AUROSTD/aurostd_main.cpp ./AUROSTD/aurostd_boot.cpp ./AUROSTD/aurostd_crc64.cpp
AUROSTD_HPP = ./AUROSTD/aurostd.h 	./AUROSTD/aurostd_argv.h   ./AUROSTD/aurostd_xoption.h   ./AUROSTD/aurostd_xscalar.h   ./AUROSTD/aurostd_xvector.h   ./AUROSTD/aurostd_xmatrix.h   ./AUROSTD/aurostd_xcomplex.h   ./AUROSTD/aurostd_xtensor.h   ./AUROSTD/aurostd_xcombos.h   ./AUROSTD/aurostd_xerror.h   ./AUROSTD/aurostd_xrandom.h 
#AUROSTD_OBJ = $(AUROSTD_CPP:.cpp=.o)
#AUROSTD_OBJ = aurostd.o
AUROSTD_OBJ	= ./AUROSTD/aurostd.o

#AUROSTD_CPP = aurostd.cpp
#AUROSTD_HPP = aurostd.h

AFLOW_HPP = aflow.h aflow_pflow.h aflow_agl_debye.h aflow_apennsy.h aflow_chull.h aflow_pocc.h aflow_pocc_old.h aflow_cce.h aflow_xgndstate.h aflow_xphases.h aflow_xvaspin.h APL/aflow_apl.h APL/aflow_qha_operations.h
AFLOWLIB_HPP = aflowlib.h

CALCULATED_CPP = aflow_data_calculated.cpp
CALCULATED_OBJ = $(CALCULATED_CPP:.cpp=.o)

AFLOW_CPP = aflow.cpp
AFLOW_CLASS_CPP = aflow_xclasses.cpp 
AFLOW_FUNCS_CPP = aflow_aflowrc.cpp aflow_init.cpp aflow_xatom.cpp aflow_xelement.cpp aflow_kbin.cpp aflow_ialien.cpp aflow_kalien.cpp aflow_kvasp.cpp aflow_ivasp.cpp aflow_avasp.cpp aflow_ovasp.cpp aflow_kaims.cpp aflow_iaims.cpp aflow_oaims.cpp aflow_lattice.cpp aflow_spacegroup.cpp aflow_mix.cpp aflow_nomix.cpp aflow_qsub.cpp aflow_pthreads.cpp aflow_symmetry.cpp aflow_neighbours.cpp aflow_surface.cpp aflow_defects.cpp aflow_wyckoff.cpp aflow_xproto.cpp aflow_xproto_gus.cpp aflow_matlab.cpp aflow_matlab_funcs.cpp aflow_gnuplot_funcs.cpp aflow_ifrozsl.cpp aflow_modules.cpp
# aflow_ifrozsl.cpp
PFLOW_FUNCS_CPP = aflow_aconvasp_main.cpp aflow_pflow_funcs.cpp aflow_pflow_print.cpp aflow_aconvasp.cpp 
ESTRUCTURE_FUNCS_CPP = aflow_estructure.cpp aflow_plotter.cpp
#[OBSOLETE] [KESONG] POCC_FUNCS_OLD_CPP = aflow_poccupation_params.cpp aflow_poccupation_edos.cpp aflow_poccupation_forcefield.cpp 
POCC_FUNCS_OLD_CPP = aflow_pocc_old.cpp
#CO 180409
POCC_FUNCS_CPP = aflow_pocc.cpp
APENNSY_FUNCS_CPP = aflow_apennsy_gndstate.cpp aflow_apennsy_structures.cpp aflow_apennsy.cpp aflow_apennsy_phases.cpp aflow_apennsy_main.cpp aflow_apennsy_vaspin.cpp
LIBRARIES_FUNCS_CPP = aflowlib_libraries.cpp 
MIX_EXPERIMENTS_CPP = aflow_mix_experiments.cpp
CONTRIB_FUNCS_SW_CPP = aflow_contrib_shidong_main.cpp aflow_contrib_shidong.cpp aflow_contrib_shidong_auxiliary.cpp aflow_contrib_shidong_cluster_expansion.cpp     
#[OBSOLETE] [JUNKAI] CONTRIB_FUNCS_JX_CPP = aflow_contrib_junkai_cproto.cpp aflow_contrib_junkai_protocheck.cpp aflow_contrib_junkai_primitivetrans.cpp aflow_contrib_junkai_basic.cpp aflow_contrib_junkai_magnet.cpp aflow_contrib_junkai_phasediag.cpp
#[OBSOLETE] [JUNKAI] CONTRIB_FUNCS_JX_CPP = aflow_contrib_junkai_basic.cpp aflow_contrib_junkai_phasediag.cpp
#[OBSOLETE] [KESONG] CONTRIB_FUNCS_KY_CPP = aflow_contrib_kesong_main.cpp aflow_contrib_kesong_ipocc.cpp aflow_contrib_kesong_hnfcell.cpp aflow_contrib_kesong_pocc_basic.cpp aflow_contrib_kesong_multienum.cpp
CONTRIB_FUNCS_RT_CPP = aflow_symmetry_spacegroup.cpp aflow_symmetry_spacegroup_functions.cpp aflow_symmetry_spacegroup_ITC_library.cpp
CONTRIB_FUNCS_WS_CPP = aflow_contrib_wahyu.cpp
CONTRIB_FUNCS_SC_CPP = 
#aflow_contrib_stefano_htsecurity.cpp
CONTRIB_FUNCS_KR_CPP = aflow_bader.cpp
CONTRIB_FUNCS_CT_CPP = aflow_contrib_cormac.cpp
AFLOWLIB_CPP = aflowlib_web_interface.cpp aflowlib_web_outreach.cpp 
AFLOWLIB_WEBAPP_CPP = aflowlib_webapp_entry.cpp aflowlib_webapp_bands.cpp
AGL_CPP = aflow_agl_debye.cpp aflow_agl_eqn_state.cpp aflow_agl_get_ev.cpp aflow_agl_plot_write.cpp aflow_agl_polynomial.cpp aflow_agl_rungibbs.cpp aflow_agl_thermal.cpp aflow_agl_hugoniot.cpp aflow_agl_electronic.cpp
AEL_CPP = aflow_ael_elasticity.cpp aflow_ael_get_stress.cpp aflow_ael_elastic_fit.cpp 
APL_CPP = APL/aflow_apl_kphonons.cpp APL/aflow_qha_operations.cpp APL/aflow_apl_function_fitting.cpp APL/aflow_apl_dirphoncalc.cpp APL/aflow_apl_qmesh.cpp APL/aflow_apl_ltet.cpp APL/aflow_apl_doscalc.cpp APL/aflow_apl_gsa.cpp APL/aflow_apl_logger.cpp APL/aflow_apl_lrphoncalc.cpp APL/aflow_apl_pathbuilder.cpp APL/aflow_apl_pdisc.cpp APL/aflow_apl_phoncalc.cpp APL/aflow_apl_shellhandle.cpp APL/aflow_apl_supercell.cpp APL/aflow_apl_thermalpc.cpp APL/aflow_aapl_cluster.cpp APL/aflow_aapl_ifcs.cpp APL/aflow_aapl_setup.cpp APL/aflow_aapl_tcond.cpp APL/aflow_qha_aflowin_creator.cpp APL/aflow_qha_dm_pdos_save.cpp APL/aflow_apl_hsqpoints.cpp APL/aflow_qha_gruneisen.cpp APL/aflow_apl_group_velocity.cpp APL/aflow_apl_atomic_disp.cpp APL/aflow_scqha_gruneisen.cpp APL/aflow_scqha_eos.cpp APL/aflow_qha_energies.cpp APL/aflow_qha3phonons_eos.cpp APL/aflow_qha_eos.cpp APL/aflow_scqha_T_freqs.cpp
ANRL_CPP = aflow_anrl.cpp ANRL/aflow_anrl_list.cpp ANRL/aflow_anrl_A10B2C_hP39_171_5c_c_a.cpp ANRL/aflow_anrl_A10B2C_hP39_172_5c_c_a.cpp ANRL/aflow_anrl_A10B3_oF52_42_2abce_ab.cpp ANRL/aflow_anrl_A12B2C_cF60_196_h_bc_a.cpp ANRL/aflow_anrl_A12B36CD12_cF488_210_h_3h_a_fg.cpp ANRL/aflow_anrl_A12B6C_cF608_210_4h_2h_e.cpp ANRL/aflow_anrl_A12B7C2_hP21_174_2j2k_ajk_cf.cpp ANRL/aflow_anrl_A12BC4_cP34_195_2j_ab_2e.cpp ANRL/aflow_anrl_A12B_cF52_225_i_a.cpp ANRL/aflow_anrl_A12B_cI26_204_g_a.cpp ANRL/aflow_anrl_A12B_tI26_139_fij_a.cpp ANRL/aflow_anrl_A13B2C2_oP34_32_a6c_c_c.cpp ANRL/aflow_anrl_A14B3C5_tP44_94_c3g_ad_bg.cpp ANRL/aflow_anrl_A15B4_cI76_220_ae_c.cpp ANRL/aflow_anrl_A17B15_cP64_207_acfk_eij.cpp ANRL/aflow_anrl_A17BC4D_tP184_89_17p_p_4p_io.cpp ANRL/aflow_anrl_A23B6_cF116_225_bd2f_e.cpp ANRL/aflow_anrl_A24BC_cF104_209_j_a_b.cpp ANRL/aflow_anrl_A2B11_cP39_200_f_aghij.cpp ANRL/aflow_anrl_A2B13C4_hP57_168_d_c6d_2d.cpp ANRL/aflow_anrl_A2B2C4D_tP18_132_e_i_o_d.cpp ANRL/aflow_anrl_A2B2C7_cF88_227_c_d_af.cpp ANRL/aflow_anrl_A2B2C_oC80_64_efg_efg_df.cpp ANRL/aflow_anrl_A2B3C12D3_cI160_230_a_c_h_d.cpp ANRL/aflow_anrl_A2B3C18D6_hP58_192_c_f_lm_l.cpp ANRL/aflow_anrl_A2B3C3DE7_hP48_145_2a_3a_3a_a_7a.cpp ANRL/aflow_anrl_A2B3C6_cP264_205_2d_ab2c2d_6d.cpp ANRL/aflow_anrl_A2B3C6_cP33_221_cd_ag_fh.cpp ANRL/aflow_anrl_A2B3C7D_oP13_47_t_aq_eqrs_h.cpp ANRL/aflow_anrl_A2B3_hP20_159_bc_2c.cpp ANRL/aflow_anrl_A2B3_hP30_169_2a_3a.cpp ANRL/aflow_anrl_A2B3_hP30_170_2a_3a.cpp ANRL/aflow_anrl_A2B3_hP5_164_d_ad.cpp ANRL/aflow_anrl_A2B3_hR10_167_c_e.cpp ANRL/aflow_anrl_A2B3_hR5_166_c_ac.cpp ANRL/aflow_anrl_A2B3_oC40_39_2d_2c2d.cpp ANRL/aflow_anrl_A2B3_oP20_60_d_cd.cpp ANRL/aflow_anrl_A2B3_oP20_62_2c_3c.cpp ANRL/aflow_anrl_A2B3_oP40_33_4a_6a.cpp ANRL/aflow_anrl_A2B3_tI80_141_ceh_3h.cpp ANRL/aflow_anrl_A2B3_tP10_127_g_ah.cpp ANRL/aflow_anrl_A2B3_tP20_102_2c_b2c.cpp ANRL/aflow_anrl_A2B3_tP20_116_bci_fj.cpp ANRL/aflow_anrl_A2B3_tP20_117_i_adgh.cpp ANRL/aflow_anrl_A2B3_tP40_137_cdf_3g.cpp ANRL/aflow_anrl_A2B3_tP5_115_g_ag.cpp ANRL/aflow_anrl_A2B4C_hR42_148_2f_4f_f.cpp ANRL/aflow_anrl_A2B4C_oC28_66_l_kl_a.cpp ANRL/aflow_anrl_A2B4C_oP28_62_ac_2cd_c.cpp ANRL/aflow_anrl_A2B5C2_oC36_37_d_c2d_d.cpp ANRL/aflow_anrl_A2B5_mC28_15_f_e2f.cpp ANRL/aflow_anrl_A2B6CD7_tP64_77_2d_6d_d_ab6d.cpp ANRL/aflow_anrl_A2B7C2_tP88_78_4a_14a_4a.cpp ANRL/aflow_anrl_A2B7_cI54_229_e_afh.cpp ANRL/aflow_anrl_A2B8C2D_tP26_100_c_abcd_c_a.cpp ANRL/aflow_anrl_A2B8CD_oP24_48_k_2m_d_b.cpp ANRL/aflow_anrl_A2B8CD_tI24_97_d_k_a_b.cpp ANRL/aflow_anrl_A2B8C_oP22_34_c_4c_a.cpp ANRL/aflow_anrl_A2B_aP6_2_2i_i.cpp ANRL/aflow_anrl_A2B_aP6_2_aei_i.cpp ANRL/aflow_anrl_A2BC2_oF40_22_fi_ad_gh.cpp ANRL/aflow_anrl_A2BC2_oI20_45_c_b_c.cpp ANRL/aflow_anrl_A2BC2_tI20_79_c_2a_c.cpp ANRL/aflow_anrl_A2BC2_tP20_105_f_ac_2e.cpp ANRL/aflow_anrl_A2BC3_oC24_63_e_c_cg.cpp ANRL/aflow_anrl_A2BC4_cF56_227_d_a_e.cpp ANRL/aflow_anrl_A2BC4D_tI16_121_d_a_i_b.cpp ANRL/aflow_anrl_A2BC4_oP28_50_ij_ac_ijm.cpp ANRL/aflow_anrl_A2BC4_tI14_82_bc_a_g.cpp ANRL/aflow_anrl_A2BC4_tP28_126_cd_e_k.cpp ANRL/aflow_anrl_A2BC4_tP28_130_f_c_g.cpp ANRL/aflow_anrl_A2BC7D2_tP24_113_e_a_cef_e.cpp ANRL/aflow_anrl_A2BC8_tI176_110_2b_b_8b.cpp ANRL/aflow_anrl_A2BCD3E6_cF208_203_e_c_d_f_g.cpp ANRL/aflow_anrl_A2BCD4_tI16_82_ac_b_d_g.cpp ANRL/aflow_anrl_A2B_cF24_227_c_a.cpp ANRL/aflow_anrl_A2B_cF24_227_d_a.cpp ANRL/aflow_anrl_A2B_cI72_211_hi_i.cpp ANRL/aflow_anrl_A2BC_oC16_67_ag_b_g.cpp ANRL/aflow_anrl_A2BC_oC8_38_e_a_b.cpp ANRL/aflow_anrl_A2B_cP12_212_c_a.cpp ANRL/aflow_anrl_A2B_cP6_224_b_a.cpp ANRL/aflow_anrl_A2BC_tP16_76_2a_a_a.cpp ANRL/aflow_anrl_A2B_hP12_194_cg_f.cpp ANRL/aflow_anrl_A2B_hP18_180_fi_bd.cpp ANRL/aflow_anrl_A2B_hP18_190_gh_bf.cpp ANRL/aflow_anrl_A2B_hP36_177_j2lm_n.cpp ANRL/aflow_anrl_A2B_hP6_191_h_e.cpp ANRL/aflow_anrl_A2B_hP9_147_g_ad.cpp ANRL/aflow_anrl_A2B_hP9_150_ef_bd.cpp ANRL/aflow_anrl_A2B_hP9_152_c_a.cpp ANRL/aflow_anrl_A2B_hP9_180_j_c.cpp ANRL/aflow_anrl_A2B_hP9_181_j_c.cpp ANRL/aflow_anrl_A2B_hP9_189_fg_bc.cpp ANRL/aflow_anrl_A2B_hR18_148_2f_f.cpp ANRL/aflow_anrl_A2B_mC144_9_24a_12a.cpp ANRL/aflow_anrl_A2B_mC48_15_ae3f_2f.cpp ANRL/aflow_anrl_A2B_mP12_13_2g_ef.cpp ANRL/aflow_anrl_A2B_mP12_14_2e_e.cpp ANRL/aflow_anrl_A2B_mP12_3_bc3e_2e.cpp ANRL/aflow_anrl_A2B_mP12_7_4a_2a.cpp ANRL/aflow_anrl_A2B_mP18_7_6a_3a.cpp ANRL/aflow_anrl_A2B_mP6_10_mn_bg.cpp ANRL/aflow_anrl_A2B_mP6_14_e_a.cpp ANRL/aflow_anrl_A2B_oC12_36_2a_a.cpp ANRL/aflow_anrl_A2B_oC12_38_de_ab.cpp ANRL/aflow_anrl_A2B_oC12_63_2c_c.cpp ANRL/aflow_anrl_A2B_oC24_20_abc_c.cpp ANRL/aflow_anrl_A2B_oC24_64_2f_f.cpp ANRL/aflow_anrl_A2B_oF24_70_e_a.cpp ANRL/aflow_anrl_A2B_oI12_72_j_a.cpp ANRL/aflow_anrl_A2B_oI12_74_h_e.cpp ANRL/aflow_anrl_A2B_oP12_17_abe_e.cpp ANRL/aflow_anrl_A2B_oP12_19_2a_a.cpp ANRL/aflow_anrl_A2B_oP12_26_abc_ab.cpp ANRL/aflow_anrl_A2B_oP12_29_2a_a.cpp ANRL/aflow_anrl_A2B_oP12_62_2c_c.cpp ANRL/aflow_anrl_A2B_oP24_52_2e_cd.cpp ANRL/aflow_anrl_A2B_oP24_55_2g2h_gh.cpp ANRL/aflow_anrl_A2B_oP24_61_2c_c.cpp ANRL/aflow_anrl_A2B_oP6_58_g_a.cpp ANRL/aflow_anrl_A2B_tI12_140_h_a.cpp ANRL/aflow_anrl_A2B_tI12_141_e_a.cpp ANRL/aflow_anrl_A2B_tI12_98_f_a.cpp ANRL/aflow_anrl_A2B_tI24_141_2e_e.cpp ANRL/aflow_anrl_A2B_tI6_139_d_a.cpp ANRL/aflow_anrl_A2B_tP12_111_2n_adf.cpp ANRL/aflow_anrl_A2B_tP12_92_b_a.cpp ANRL/aflow_anrl_A2B_tP24_135_gh_h.cpp ANRL/aflow_anrl_A2B_tP30_85_ab2g_cg.cpp ANRL/aflow_anrl_A2B_tP36_96_3b_ab.cpp ANRL/aflow_anrl_A2B_tP48_77_8d_4d.cpp ANRL/aflow_anrl_A2B_tP6_129_ac_c.cpp ANRL/aflow_anrl_A2B_tP6_136_f_a.cpp ANRL/aflow_anrl_A2B_tP6_137_d_a.cpp ANRL/aflow_anrl_A3B10_cI52_229_e_fh.cpp ANRL/aflow_anrl_A3B11C6_tP40_100_ac_bc2d_cd.cpp ANRL/aflow_anrl_A3B13_oC32_38_ac_a2bcdef.cpp ANRL/aflow_anrl_A3B2_cI40_220_d_c.cpp ANRL/aflow_anrl_A3B2_hP10_176_h_bc.cpp ANRL/aflow_anrl_A3B2_hP5_164_ad_d.cpp ANRL/aflow_anrl_A3B2_hR5_155_e_c.cpp ANRL/aflow_anrl_A3B2_oP20_52_de_cd.cpp ANRL/aflow_anrl_A3B2_oP20_56_ce_e.cpp ANRL/aflow_anrl_A3B2_oP20_62_3c_2c.cpp ANRL/aflow_anrl_A3B2_tP10_83_adk_j.cpp ANRL/aflow_anrl_A3B3C_cI56_214_g_h_a.cpp ANRL/aflow_anrl_A3B3C_hP14_176_h_h_c.cpp ANRL/aflow_anrl_A3B4C_cP16_218_c_e_a.cpp ANRL/aflow_anrl_A3B4C_cP8_215_d_e_a.cpp ANRL/aflow_anrl_A3B4_cF56_227_ad_e.cpp ANRL/aflow_anrl_A3B4_tI28_141_ad_h.cpp ANRL/aflow_anrl_A3B5_oC16_65_ah_bej.cpp ANRL/aflow_anrl_A3B5_oC32_38_abce_abcdf.cpp ANRL/aflow_anrl_A3B5_oP16_55_ch_agh.cpp ANRL/aflow_anrl_A3B5_tI32_108_ac_a2c.cpp ANRL/aflow_anrl_A3B5_tI32_140_ah_bk.cpp ANRL/aflow_anrl_A3B5_tI32_140_ah_cl.cpp ANRL/aflow_anrl_A3B7_hP20_186_c_b2c.cpp ANRL/aflow_anrl_A3B7_oP40_62_cd_3c2d.cpp ANRL/aflow_anrl_A3B7_tP40_76_3a_7a.cpp ANRL/aflow_anrl_A3BC2_cI48_214_f_a_e.cpp ANRL/aflow_anrl_A3BC2_oP48_50_3m_m_2m.cpp ANRL/aflow_anrl_A3BC3D_tP64_106_3c_c_3c_c.cpp ANRL/aflow_anrl_A3BC_hP10_188_k_c_a.cpp ANRL/aflow_anrl_A3BC_hP30_185_cd_c_ab.cpp ANRL/aflow_anrl_A3BC_hR5_146_b_a_a.cpp ANRL/aflow_anrl_A3B_cI32_204_g_c.cpp ANRL/aflow_anrl_A3B_cI8_229_b_a.cpp ANRL/aflow_anrl_A3BC_mC10_8_ab_a_a.cpp ANRL/aflow_anrl_A3B_cP16_198_b_a.cpp ANRL/aflow_anrl_A3B_cP16_208_j_b.cpp ANRL/aflow_anrl_A3B_cP4_221_d_a.cpp ANRL/aflow_anrl_A3B_cP8_223_c_a.cpp ANRL/aflow_anrl_A3BC_tP5_99_bc_a_b.cpp ANRL/aflow_anrl_A3B_hP16_194_gh_ac.cpp ANRL/aflow_anrl_A3B_hP24_151_3c_2a.cpp ANRL/aflow_anrl_A3B_hP24_153_3c_2b.cpp ANRL/aflow_anrl_A3B_hP24_165_adg_f.cpp ANRL/aflow_anrl_A3B_hP24_165_bdg_f.cpp ANRL/aflow_anrl_A3B_hP24_185_ab2c_c.cpp ANRL/aflow_anrl_A3B_hP4_191_bc_a.cpp ANRL/aflow_anrl_A3B_hP8_158_d_a.cpp ANRL/aflow_anrl_A3B_hP8_173_c_b.cpp ANRL/aflow_anrl_A3B_hP8_176_h_c.cpp ANRL/aflow_anrl_A3B_hP8_185_c_a.cpp ANRL/aflow_anrl_A3B_hP8_194_h_c.cpp ANRL/aflow_anrl_A3B_hR4_160_b_a.cpp ANRL/aflow_anrl_A3B_mP16_7_6a_2a.cpp ANRL/aflow_anrl_A3B_oC64_66_gi2lm_2l.cpp ANRL/aflow_anrl_A3B_oC64_66_kl2m_bdl.cpp ANRL/aflow_anrl_A3B_oI32_23_ij2k_k.cpp ANRL/aflow_anrl_A3B_oP16_62_cd_c.cpp ANRL/aflow_anrl_A3B_oP32_60_3d_d.cpp ANRL/aflow_anrl_A3B_oP8_59_bf_a.cpp ANRL/aflow_anrl_A3B_tI16_139_cde_e.cpp ANRL/aflow_anrl_A3B_tI24_119_b2i_af.cpp ANRL/aflow_anrl_A3B_tI32_82_3g_g.cpp ANRL/aflow_anrl_A3B_tI8_139_bd_a.cpp ANRL/aflow_anrl_A3B_tP16_118_ei_f.cpp ANRL/aflow_anrl_A3B_tP32_114_3e_e.cpp ANRL/aflow_anrl_A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h.cpp ANRL/aflow_anrl_A4B14C3_hP21_143_bd_ac4d_d.cpp ANRL/aflow_anrl_A4B2C13D_tP40_90_g_d_cef2g_c.cpp ANRL/aflow_anrl_A4B2C6D16E_cF232_203_e_d_f_eg_a.cpp ANRL/aflow_anrl_A4B2C_tP28_135_gh_h_d.cpp ANRL/aflow_anrl_A4B3_cI112_230_af_g.cpp ANRL/aflow_anrl_A4B3_cI14_229_c_b.cpp ANRL/aflow_anrl_A4B3_cI28_220_c_a.cpp ANRL/aflow_anrl_A4B3_hP14_173_bc_c.cpp ANRL/aflow_anrl_A4B3_hP28_159_ab2c_2c.cpp ANRL/aflow_anrl_A4B3_hR7_166_2c_ac.cpp ANRL/aflow_anrl_A4B3_oI14_71_gh_cg.cpp ANRL/aflow_anrl_A4B5_tI18_139_i_ah.cpp ANRL/aflow_anrl_A4B5_tI18_87_h_ah.cpp ANRL/aflow_anrl_A4B6C_hP11_143_bd_2d_a.cpp ANRL/aflow_anrl_A4B9_cP52_215_ei_3efgi.cpp ANRL/aflow_anrl_A4BC2_tI28_120_i_d_e.cpp ANRL/aflow_anrl_A4BC4D_tP10_123_gh_a_i_d.cpp ANRL/aflow_anrl_A4BC4_tP18_137_g_b_g.cpp ANRL/aflow_anrl_A4B_cI10_217_c_a.cpp ANRL/aflow_anrl_A4B_cI10_229_c_a.cpp ANRL/aflow_anrl_A4B_cI40_197_cde_c.cpp ANRL/aflow_anrl_A4BC_tI24_141_h_b_a.cpp ANRL/aflow_anrl_A4B_hP15_144_4a_a.cpp ANRL/aflow_anrl_A4B_oI20_74_beh_e.cpp ANRL/aflow_anrl_A4B_oP20_62_2cd_c.cpp ANRL/aflow_anrl_A4B_tI10_139_de_a.cpp ANRL/aflow_anrl_A4B_tI20_88_f_a.cpp ANRL/aflow_anrl_A4B_tP10_114_e_a.cpp ANRL/aflow_anrl_A4B_tP10_125_m_a.cpp ANRL/aflow_anrl_A4B_tP20_127_ehj_g.cpp ANRL/aflow_anrl_A5B11_mP16_6_2abc_2a3b3c.cpp ANRL/aflow_anrl_A5B2_hP14_194_abdf_f.cpp ANRL/aflow_anrl_A5B2_hP28_194_ahk_ch.cpp ANRL/aflow_anrl_A5B2_hR7_166_a2c_c.cpp ANRL/aflow_anrl_A5B2_mC14_12_a2i_i.cpp ANRL/aflow_anrl_A5B2_oP14_49_dehq_ab.cpp ANRL/aflow_anrl_A5B3C15_oP46_30_a2c_bc_a7c.cpp ANRL/aflow_anrl_A5B3C16_cP96_222_ce_d_fi.cpp ANRL/aflow_anrl_A5B3C_hP18_186_2a3b_2ab_b.cpp ANRL/aflow_anrl_A5B3_hP16_190_bdh_g.cpp ANRL/aflow_anrl_A5B3_hP16_193_dg_g.cpp ANRL/aflow_anrl_A5B3_mC32_9_5a_3a.cpp ANRL/aflow_anrl_A5B3_tP32_118_g2i_aceh.cpp ANRL/aflow_anrl_A5B3_tP32_130_cg_cf.cpp ANRL/aflow_anrl_A5B5C4_tP28_104_ac_ac_c.cpp ANRL/aflow_anrl_A5B6C2_hP13_157_2ac_2c_b.cpp ANRL/aflow_anrl_A5B7_tI24_107_ac_abd.cpp ANRL/aflow_anrl_A5B8_cI52_217_ce_cg.cpp ANRL/aflow_anrl_A5BCD6_cF416_228_eg_c_b_h.cpp ANRL/aflow_anrl_A5B_oP24_26_3a3b2c_ab.cpp ANRL/aflow_anrl_A6B23_cF116_225_e_acfh.cpp ANRL/aflow_anrl_A6B2C_cF36_225_e_c_a.cpp ANRL/aflow_anrl_A6B2CD6E_cP64_208_m_ad_b_m_c.cpp ANRL/aflow_anrl_A6B2C_tP18_128_eh_d_a.cpp ANRL/aflow_anrl_A6B2C_tP18_94_eg_c_a.cpp ANRL/aflow_anrl_A6B4C16D_oP108_27_abcd4e_4e_16e_e.cpp ANRL/aflow_anrl_A6B6C_cF104_202_h_h_c.cpp ANRL/aflow_anrl_A6B_cF224_228_h_c.cpp ANRL/aflow_anrl_A6B_cP7_221_f_a.cpp ANRL/aflow_anrl_A6B_hR7_166_g_a.cpp ANRL/aflow_anrl_A6B_oC28_63_efg_c.cpp ANRL/aflow_anrl_A7B2C2_mC22_12_aij_h_i.cpp ANRL/aflow_anrl_A7B2C_tP40_128_egi_h_e.cpp ANRL/aflow_anrl_A7B3_cI40_229_df_e.cpp ANRL/aflow_anrl_A7B6_hR13_166_ah_3c.cpp ANRL/aflow_anrl_A7B7C2_tP32_101_bde_ade_d.cpp ANRL/aflow_anrl_A7B8_mP120_14_14e_16e.cpp ANRL/aflow_anrl_A7B8_oP120_60_7d_8d.cpp ANRL/aflow_anrl_A7BC3D13_cF192_219_de_b_c_ah.cpp ANRL/aflow_anrl_A7B_cF32_225_bd_a.cpp ANRL/aflow_anrl_A8B2C12D2E_oI50_23_bcfk_i_3k_j_a.cpp ANRL/aflow_anrl_A8B5_hR26_160_a3bc_a3b.cpp ANRL/aflow_anrl_A8B5_mP13_6_a7b_3a2b.cpp ANRL/aflow_anrl_A8B7C6_hP21_175_ck_aj_k.cpp ANRL/aflow_anrl_A8BC3D6_hP18_189_bfh_a_g_i.cpp ANRL/aflow_anrl_A8B_tI18_139_hi_a.cpp ANRL/aflow_anrl_A9B16C7_cF128_225_acd_2f_be.cpp ANRL/aflow_anrl_A9B2_mP22_7_9a_2a.cpp ANRL/aflow_anrl_A9B3C_hP26_194_hk_h_a.cpp ANRL/aflow_anrl_A9BC3D5_hP18_189_fi_a_g_bh.cpp ANRL/aflow_anrl_A9BC_oC44_39_3c3d_a_c.cpp ANRL/aflow_anrl_A_aP4_2_aci.cpp ANRL/aflow_anrl_AB11CD3_cP16_221_a_dg_b_c.cpp ANRL/aflow_anrl_AB11_cP36_221_c_agij.cpp ANRL/aflow_anrl_AB12C3_cI32_229_a_h_b.cpp ANRL/aflow_anrl_AB13_cF112_226_a_bi.cpp ANRL/aflow_anrl_AB18C8_cF108_225_a_eh_f.cpp ANRL/aflow_anrl_AB27CD3_cP32_221_a_dij_b_c.cpp ANRL/aflow_anrl_AB2_aP12_1_4a_8a.cpp ANRL/aflow_anrl_AB2C12D4_tP76_75_2a2b_2d_12d_4d.cpp ANRL/aflow_anrl_AB2C3_oP24_62_c_d_cd.cpp ANRL/aflow_anrl_AB2C4_tI14_139_a_e_ce.cpp ANRL/aflow_anrl_AB2C8D_oP24_49_g_q_2qr_e.cpp ANRL/aflow_anrl_AB2C_cF16_225_a_c_b.cpp ANRL/aflow_anrl_AB2CD2_hP36_163_h_i_bf_i.cpp ANRL/aflow_anrl_AB2_cF12_225_a_c.cpp ANRL/aflow_anrl_AB2_cF48_227_c_e.cpp ANRL/aflow_anrl_AB2_cF96_227_e_cf.cpp ANRL/aflow_anrl_AB2C_oC16_40_a_2b_b.cpp ANRL/aflow_anrl_AB2C_oC16_63_c_2c_c.cpp ANRL/aflow_anrl_AB2C_oP16_62_c_2c_c.cpp ANRL/aflow_anrl_AB2_cP12_205_a_c.cpp ANRL/aflow_anrl_AB2_hP12_143_cd_ab2d.cpp ANRL/aflow_anrl_AB2_hP12_194_f_ah.cpp ANRL/aflow_anrl_AB2_hP24_194_ef_fgh.cpp ANRL/aflow_anrl_AB2_hP3_164_a_d.cpp ANRL/aflow_anrl_AB2_hP3_191_a_d.cpp ANRL/aflow_anrl_AB2_hP6_194_b_f.cpp ANRL/aflow_anrl_AB2_hP6_194_c_ad.cpp ANRL/aflow_anrl_AB2_hP6_194_c_f.cpp ANRL/aflow_anrl_AB2_hP72_192_m_j2kl.cpp ANRL/aflow_anrl_AB2_hP9_156_b2c_3a2bc.cpp ANRL/aflow_anrl_AB2_hP9_162_ad_k.cpp ANRL/aflow_anrl_AB2_hP9_164_bd_c2d.cpp ANRL/aflow_anrl_AB2_hP9_180_d_j.cpp ANRL/aflow_anrl_AB2_mC6_12_a_i.cpp ANRL/aflow_anrl_AB2_oC24_41_2a_2b.cpp ANRL/aflow_anrl_AB2_oC6_21_a_k.cpp ANRL/aflow_anrl_AB2_oF48_70_f_fg.cpp ANRL/aflow_anrl_AB2_oF72_43_ab_3b.cpp ANRL/aflow_anrl_AB2_oI6_71_a_g.cpp ANRL/aflow_anrl_AB2_oI6_71_a_i.cpp ANRL/aflow_anrl_AB2_oP12_29_a_2a.cpp ANRL/aflow_anrl_AB2_oP24_28_acd_2c3d.cpp ANRL/aflow_anrl_AB2_oP6_34_a_c.cpp ANRL/aflow_anrl_AB2_oP6_58_a_g.cpp ANRL/aflow_anrl_AB2_tI48_80_2b_4b.cpp ANRL/aflow_anrl_AB2_tI6_139_a_e.cpp ANRL/aflow_anrl_AB2_tI96_88_2f_4f.cpp ANRL/aflow_anrl_AB2_tP12_115_j_egi.cpp ANRL/aflow_anrl_AB2_tP12_81_adg_2h.cpp ANRL/aflow_anrl_AB2_tP6_137_a_d.cpp ANRL/aflow_anrl_AB32C48_cI162_204_a_2efg_2gh.cpp ANRL/aflow_anrl_AB32CD4E8_tP184_93_i_16p_af_2p_4p.cpp ANRL/aflow_anrl_AB3C16_cF160_203_a_bc_eg.cpp ANRL/aflow_anrl_AB3C2_cI96_206_c_e_ad.cpp ANRL/aflow_anrl_AB3C3_cF112_227_c_de_f.cpp ANRL/aflow_anrl_AB3C4_cP8_215_a_c_e.cpp ANRL/aflow_anrl_AB3C4_hP16_194_c_af_ef.cpp ANRL/aflow_anrl_AB3C4_oP16_31_a_ab_2ab.cpp ANRL/aflow_anrl_AB3C4_oP32_33_a_3a_4a.cpp ANRL/aflow_anrl_AB3C6_cI80_206_a_d_e.cpp ANRL/aflow_anrl_AB3C_cP5_221_a_c_b.cpp ANRL/aflow_anrl_AB3C_cP60_201_be_fh_g.cpp ANRL/aflow_anrl_AB3_cF16_225_a_bc.cpp ANRL/aflow_anrl_AB3C_hR10_148_c_f_c.cpp ANRL/aflow_anrl_AB3C_hR10_167_b_e_a.cpp ANRL/aflow_anrl_AB3C_oC20_63_a_cf_c.cpp ANRL/aflow_anrl_AB3C_oP20_62_c_cd_a.cpp ANRL/aflow_anrl_AB3_cP4_221_a_c.cpp ANRL/aflow_anrl_AB3_hP24_149_acgi_3l.cpp ANRL/aflow_anrl_AB3_hP24_178_b_ac.cpp ANRL/aflow_anrl_AB3_hP24_179_b_ac.cpp ANRL/aflow_anrl_AB3_hP24_185_c_ab2c.cpp ANRL/aflow_anrl_AB3_hP4_187_e_fh.cpp ANRL/aflow_anrl_AB3_hP8_182_c_g.cpp ANRL/aflow_anrl_AB3_hP8_194_c_bf.cpp ANRL/aflow_anrl_AB3_hR8_148_c_f.cpp ANRL/aflow_anrl_AB3_hR8_155_c_de.cpp ANRL/aflow_anrl_AB3_mC16_12_g_ij.cpp ANRL/aflow_anrl_AB3_mC16_15_e_cf.cpp ANRL/aflow_anrl_AB3_mC16_9_a_3a.cpp ANRL/aflow_anrl_AB3_mC32_8_4a_12a.cpp ANRL/aflow_anrl_AB3_mC32_8_4a_4a4b.cpp ANRL/aflow_anrl_AB3_mP16_10_mn_3m3n.cpp ANRL/aflow_anrl_AB3_oC16_40_b_3b.cpp ANRL/aflow_anrl_AB3_oC8_65_a_bf.cpp ANRL/aflow_anrl_AB3_oP16_18_ab_3c.cpp ANRL/aflow_anrl_AB3_oP16_19_a_3a.cpp ANRL/aflow_anrl_AB3_oP16_62_c_3c.cpp ANRL/aflow_anrl_AB3_oP16_62_c_cd.cpp ANRL/aflow_anrl_AB3_tI16_140_b_ah.cpp ANRL/aflow_anrl_AB3_tP32_133_h_i2j.cpp ANRL/aflow_anrl_AB3_tP32_86_g_3g.cpp ANRL/aflow_anrl_AB3_tP4_123_a_ce.cpp ANRL/aflow_anrl_AB3_tP8_113_a_ce.cpp ANRL/aflow_anrl_AB4C17D4E_tP54_90_a_g_c4g_g_c.cpp ANRL/aflow_anrl_AB4C3_cI16_229_a_c_b.cpp ANRL/aflow_anrl_AB4C7D_hP26_159_b_ac_a2c_b.cpp ANRL/aflow_anrl_AB4C_hP6_191_a_h_b.cpp ANRL/aflow_anrl_AB4C_hP72_168_2d_8d_2d.cpp ANRL/aflow_anrl_AB4C_hP72_184_d_4d_d.cpp ANRL/aflow_anrl_AB4C_oC24_63_a_fg_c.cpp ANRL/aflow_anrl_AB4C_oC24_63_c_fg_c.cpp ANRL/aflow_anrl_AB4C_oP24_62_c_2cd_c.cpp ANRL/aflow_anrl_AB4_cP5_215_a_e.cpp ANRL/aflow_anrl_AB4C_tI12_82_c_g_a.cpp ANRL/aflow_anrl_AB4C_tP12_112_b_n_e.cpp ANRL/aflow_anrl_AB4C_tP12_124_a_m_c.cpp ANRL/aflow_anrl_AB4_oC20_41_a_2b.cpp ANRL/aflow_anrl_AB4_oC20_68_a_i.cpp ANRL/aflow_anrl_AB4_tI10_87_a_h.cpp ANRL/aflow_anrl_AB4_tP10_103_a_d.cpp ANRL/aflow_anrl_AB4_tP10_124_a_m.cpp ANRL/aflow_anrl_AB5_cF24_216_a_ce.cpp ANRL/aflow_anrl_AB5C_tP7_123_b_ci_a.cpp ANRL/aflow_anrl_AB5_hP6_191_a_cg.cpp ANRL/aflow_anrl_AB6C4_tP22_104_a_2ac_c.cpp ANRL/aflow_anrl_AB6C_tP16_132_d_io_a.cpp ANRL/aflow_anrl_AB7CD2_oI44_24_a_b3d_c_ac.cpp ANRL/aflow_anrl_AB7_hR16_166_c_c2h.cpp ANRL/aflow_anrl_AB8C2_oC22_35_a_ab3e_e.cpp ANRL/aflow_anrl_AB8C2_tI44_97_e_2k_cd.cpp ANRL/aflow_anrl_AB9C4_hP28_188_e_kl_ak.cpp ANRL/aflow_anrl_AB_aP16_2_4i_4i.cpp ANRL/aflow_anrl_ABC2_aP16_1_4a_4a_8a.cpp ANRL/aflow_anrl_ABC2_hP4_164_a_b_d.cpp ANRL/aflow_anrl_ABC2_hP8_194_d_a_f.cpp ANRL/aflow_anrl_ABC2_hR24_167_e_e_2e.cpp ANRL/aflow_anrl_ABC2_hR4_166_a_b_c.cpp ANRL/aflow_anrl_ABC2_mP8_10_ac_eh_mn.cpp ANRL/aflow_anrl_ABC2_oC16_67_b_g_ag.cpp ANRL/aflow_anrl_ABC2_oI16_23_ab_i_k.cpp ANRL/aflow_anrl_ABC2_oP16_53_h_e_gh.cpp ANRL/aflow_anrl_ABC2_tI16_122_a_b_d.cpp ANRL/aflow_anrl_ABC2_tP4_123_d_a_f.cpp ANRL/aflow_anrl_ABC3_cP20_198_a_a_b.cpp ANRL/aflow_anrl_ABC3_hR10_146_2a_2a_2b.cpp ANRL/aflow_anrl_ABC3_hR10_161_a_a_b.cpp ANRL/aflow_anrl_ABC3_hR10_167_a_b_e.cpp ANRL/aflow_anrl_ABC3_mP10_11_e_e_ef.cpp ANRL/aflow_anrl_ABC3_oP20_30_2a_c_3c.cpp ANRL/aflow_anrl_ABC3_oP20_53_e_g_hi.cpp ANRL/aflow_anrl_ABC3_oP20_54_e_d_cf.cpp ANRL/aflow_anrl_ABC4_mP12_13_e_a_2g.cpp ANRL/aflow_anrl_ABC4_oI12_23_a_b_k.cpp ANRL/aflow_anrl_ABC4_oP12_16_ag_cd_2u.cpp ANRL/aflow_anrl_ABC4_tI96_142_e_ab_2g.cpp ANRL/aflow_anrl_ABC4_tP12_125_a_b_m.cpp ANRL/aflow_anrl_ABC6D2_mC40_15_e_e_3f_f.cpp ANRL/aflow_anrl_ABC_cF12_216_b_c_a.cpp ANRL/aflow_anrl_ABC_cP12_198_a_a_a.cpp ANRL/aflow_anrl_ABCD3_oI48_73_d_e_e_ef.cpp ANRL/aflow_anrl_ABCD_cF16_216_c_d_b_a.cpp ANRL/aflow_anrl_ABCD_oP16_57_d_c_d_d.cpp ANRL/aflow_anrl_ABCD_tP8_129_c_b_a_c.cpp ANRL/aflow_anrl_AB_cF16_227_a_b.cpp ANRL/aflow_anrl_AB_cF8_216_c_a.cpp ANRL/aflow_anrl_AB_cF8_225_a_b.cpp ANRL/aflow_anrl_ABC_hP12_174_cj_fk_aj.cpp ANRL/aflow_anrl_ABC_hP3_183_a_a_a.cpp ANRL/aflow_anrl_ABC_hP3_187_a_d_f.cpp ANRL/aflow_anrl_ABC_hP36_175_jk_jk_jk.cpp ANRL/aflow_anrl_ABC_hP6_194_c_d_a.cpp ANRL/aflow_anrl_ABC_hR3_160_a_a_a.cpp ANRL/aflow_anrl_ABC_hR6_166_c_c_c.cpp ANRL/aflow_anrl_AB_cI16_199_a_a.cpp ANRL/aflow_anrl_ABC_oI12_71_h_j_g.cpp ANRL/aflow_anrl_ABC_oI36_46_ac_bc_3b.cpp ANRL/aflow_anrl_ABC_oP12_29_a_a_a.cpp ANRL/aflow_anrl_ABC_oP6_59_a_a_a.cpp ANRL/aflow_anrl_ABC_oP6_59_a_b_a.cpp ANRL/aflow_anrl_AB_cP16_205_c_c.cpp ANRL/aflow_anrl_AB_cP2_221_b_a.cpp ANRL/aflow_anrl_AB_cP6_221_c_d.cpp ANRL/aflow_anrl_AB_cP8_198_a_a.cpp ANRL/aflow_anrl_ABC_tI12_109_a_a_a.cpp ANRL/aflow_anrl_ABC_tP24_91_d_d_d.cpp ANRL/aflow_anrl_ABC_tP24_95_d_d_d.cpp ANRL/aflow_anrl_ABC_tP6_129_c_a_c.cpp ANRL/aflow_anrl_AB_hP12_156_2ab3c_2ab3c.cpp ANRL/aflow_anrl_AB_hP12_186_a2b_a2b.cpp ANRL/aflow_anrl_AB_hP12_194_af_bf.cpp ANRL/aflow_anrl_AB_hP12_194_df_ce.cpp ANRL/aflow_anrl_AB_hP2_187_d_a.cpp ANRL/aflow_anrl_AB_hP24_190_i_afh.cpp ANRL/aflow_anrl_AB_hP4_156_ab_ab.cpp ANRL/aflow_anrl_AB_hP4_186_b_a.cpp ANRL/aflow_anrl_AB_hP4_186_b_b.cpp ANRL/aflow_anrl_AB_hP4_194_c_a.cpp ANRL/aflow_anrl_AB_hP4_194_c_d.cpp ANRL/aflow_anrl_AB_hP6_144_a_a.cpp ANRL/aflow_anrl_AB_hP6_154_a_b.cpp ANRL/aflow_anrl_AB_hP6_183_c_ab.cpp ANRL/aflow_anrl_AB_hP6_191_f_ad.cpp ANRL/aflow_anrl_AB_hP8_186_ab_ab.cpp ANRL/aflow_anrl_AB_hP8_194_ad_f.cpp ANRL/aflow_anrl_AB_hR10_160_5a_5a.cpp ANRL/aflow_anrl_AB_hR16_148_cf_cf.cpp ANRL/aflow_anrl_AB_hR2_166_a_b.cpp ANRL/aflow_anrl_AB_hR26_148_b2f_a2f.cpp ANRL/aflow_anrl_AB_hR6_160_3a_3a.cpp ANRL/aflow_anrl_AB_hR6_160_b_b.cpp ANRL/aflow_anrl_AB_mC8_15_c_e.cpp ANRL/aflow_anrl_AB_mP4_11_e_e.cpp ANRL/aflow_anrl_AB_mP4_6_2b_2a.cpp ANRL/aflow_anrl_AB_mP6_10_en_am.cpp ANRL/aflow_anrl_AB_oC8_36_a_a.cpp ANRL/aflow_anrl_AB_oC8_63_c_c.cpp ANRL/aflow_anrl_AB_oC8_65_j_g.cpp ANRL/aflow_anrl_AB_oC8_67_a_g.cpp ANRL/aflow_anrl_AB_oF8_22_a_c.cpp ANRL/aflow_anrl_AB_oF8_42_a_a.cpp ANRL/aflow_anrl_AB_oF8_69_a_b.cpp ANRL/aflow_anrl_AB_oI4_44_a_b.cpp ANRL/aflow_anrl_AB_oP16_61_c_c.cpp ANRL/aflow_anrl_AB_oP2_25_b_a.cpp ANRL/aflow_anrl_AB_oP4_51_e_f.cpp ANRL/aflow_anrl_AB_oP4_59_a_b.cpp ANRL/aflow_anrl_AB_oP48_61_3c_3c.cpp ANRL/aflow_anrl_AB_oP8_33_a_a.cpp ANRL/aflow_anrl_AB_oP8_57_d_d.cpp ANRL/aflow_anrl_AB_oP8_62_c_c.cpp ANRL/aflow_anrl_AB_tI16_140_ab_h.cpp ANRL/aflow_anrl_AB_tI16_141_e_e.cpp ANRL/aflow_anrl_AB_tI4_107_a_a.cpp ANRL/aflow_anrl_AB_tI4_119_c_a.cpp ANRL/aflow_anrl_AB_tI8_109_a_a.cpp ANRL/aflow_anrl_AB_tI8_139_e_e.cpp ANRL/aflow_anrl_AB_tI8_141_a_b.cpp ANRL/aflow_anrl_AB_tP16_84_cej_k.cpp ANRL/aflow_anrl_AB_tP2_123_a_d.cpp ANRL/aflow_anrl_AB_tP4_129_a_c.cpp ANRL/aflow_anrl_AB_tP4_129_c_c.cpp ANRL/aflow_anrl_AB_tP4_131_c_e.cpp ANRL/aflow_anrl_AB_tP8_111_n_n.cpp ANRL/aflow_anrl_AB_tP8_136_g_f.cpp ANRL/aflow_anrl_A_cF136_227_aeg.cpp ANRL/aflow_anrl_A_cF240_202_h2i.cpp ANRL/aflow_anrl_A_cF4_225_a.cpp ANRL/aflow_anrl_A_cF8_227_a.cpp ANRL/aflow_anrl_A_cI16_206_c.cpp ANRL/aflow_anrl_A_cI16_220_c.cpp ANRL/aflow_anrl_A_cI2_229_a.cpp ANRL/aflow_anrl_A_cI58_217_ac2g.cpp ANRL/aflow_anrl_A_cP1_221_a.cpp ANRL/aflow_anrl_A_cP20_213_cd.cpp ANRL/aflow_anrl_A_cP240_205_10d.cpp ANRL/aflow_anrl_A_cP46_223_dik.cpp ANRL/aflow_anrl_A_cP8_198_2a.cpp ANRL/aflow_anrl_A_cP8_205_c.cpp ANRL/aflow_anrl_A_hP1_191_a.cpp ANRL/aflow_anrl_A_hP2_194_c.cpp ANRL/aflow_anrl_A_hP3_152_a.cpp ANRL/aflow_anrl_A_hP4_186_ab.cpp ANRL/aflow_anrl_A_hP4_194_ac.cpp ANRL/aflow_anrl_A_hP4_194_bc.cpp ANRL/aflow_anrl_A_hP4_194_f.cpp ANRL/aflow_anrl_A_hP6_178_a.cpp ANRL/aflow_anrl_A_hP6_194_h.cpp ANRL/aflow_anrl_A_hP9_154_bc.cpp ANRL/aflow_anrl_A_hR105_166_bc9h4i.cpp ANRL/aflow_anrl_A_hR1_166_a.cpp ANRL/aflow_anrl_A_hR12_166_2h.cpp ANRL/aflow_anrl_A_hR2_166_c.cpp ANRL/aflow_anrl_A_hR3_166_ac.cpp ANRL/aflow_anrl_A_mC12_5_3c.cpp ANRL/aflow_anrl_A_mC16_12_4i.cpp ANRL/aflow_anrl_A_mC24_15_2e2f.cpp ANRL/aflow_anrl_A_mC34_12_ah3i2j.cpp ANRL/aflow_anrl_A_mC4_12_i.cpp ANRL/aflow_anrl_A_mP16_11_8e.cpp ANRL/aflow_anrl_A_mP32_14_8e.cpp ANRL/aflow_anrl_A_mP4_4_2a.cpp ANRL/aflow_anrl_A_mP64_14_16e.cpp ANRL/aflow_anrl_A_mP8_10_2m2n.cpp ANRL/aflow_anrl_A_mP84_13_21g.cpp ANRL/aflow_anrl_A_oC4_63_c.cpp ANRL/aflow_anrl_A_oC8_64_f.cpp ANRL/aflow_anrl_A_oF128_70_4h.cpp ANRL/aflow_anrl_A_oF8_70_a.cpp ANRL/aflow_anrl_A_oP16_55_2g2h.cpp ANRL/aflow_anrl_A_oP8_62_2c.cpp ANRL/aflow_anrl_A_tI16_142_f.cpp ANRL/aflow_anrl_A_tI2_139_a.cpp ANRL/aflow_anrl_A_tI4_139_e.cpp ANRL/aflow_anrl_A_tI4_141_a.cpp ANRL/aflow_anrl_A_tI8_139_h.cpp ANRL/aflow_anrl_A_tP12_138_bi.cpp ANRL/aflow_anrl_A_tP12_96_ab.cpp ANRL/aflow_anrl_A_tP16_138_j.cpp ANRL/aflow_anrl_A_tP30_136_bf2ij.cpp ANRL/aflow_anrl_A_tP4_129_ac.cpp ANRL/aflow_anrl_A_tP4_136_f.cpp ANRL/aflow_anrl_A_tP50_134_b2m2n.cpp ANRL/aflow_anrl_sigma_tP30_136_bf2ij.cpp 
#[CO191217 - has duplicates]ANRL/aflow_anrl_AB2_aP12_1_4a_8a.cpp ANRL/aflow_anrl_ABC2_aP16_1_4a_4a_8a.cpp ANRL/aflow_anrl_A2B_aP6_2_2i_i.cpp ANRL/aflow_anrl_A_aP4_2_aci.cpp ANRL/aflow_anrl_A2B_mP12_3_bc3e_2e.cpp ANRL/aflow_anrl_A_mP4_4_2a.cpp ANRL/aflow_anrl_A_mC12_5_3c.cpp ANRL/aflow_anrl_A3BC_mC10_8_ab_a_a.cpp ANRL/aflow_anrl_A2B_mC144_9_24a_12a.cpp ANRL/aflow_anrl_AB_mP4_11_e_e.cpp ANRL/aflow_anrl_ABC3_mP10_11_e_e_ef.cpp ANRL/aflow_anrl_A_mP16_11_8e.cpp ANRL/aflow_anrl_AB2_mC6_12_a_i.cpp ANRL/aflow_anrl_A_mC34_12_ah3i2j.cpp ANRL/aflow_anrl_AB3_mC16_12_g_ij.cpp ANRL/aflow_anrl_A5B2_mC14_12_a2i_i.cpp ANRL/aflow_anrl_A_mC4_12_i.cpp ANRL/aflow_anrl_ABC4_mP12_13_e_a_2g.cpp ANRL/aflow_anrl_A_mP84_13_21g.cpp ANRL/aflow_anrl_A2B_mP12_14_2e_e.cpp ANRL/aflow_anrl_A_mP32_14_8e.cpp ANRL/aflow_anrl_A_mP64_14_16e.cpp ANRL/aflow_anrl_A2B5_mC28_15_f_e2f.cpp ANRL/aflow_anrl_AB_mC8_15_c_e.cpp ANRL/aflow_anrl_A2B_mC48_15_ae3f_2f.cpp ANRL/aflow_anrl_ABC6D2_mC40_15_e_e_3f_f.cpp ANRL/aflow_anrl_ABC4_oP12_16_ag_cd_2u.cpp ANRL/aflow_anrl_AB3_oP16_18_ab_3c.cpp ANRL/aflow_anrl_A2B_oP12_19_2a_a.cpp ANRL/aflow_anrl_A2B_oC24_20_abc_c.cpp ANRL/aflow_anrl_AB_oP2_25_b_a.cpp ANRL/aflow_anrl_AB2_oP24_28_acd_2c3d.cpp ANRL/aflow_anrl_AB3C4_oP16_31_a_ab_2ab.cpp ANRL/aflow_anrl_AB_oP8_33_a_a.cpp ANRL/aflow_anrl_AB3C4_oP32_33_a_3a_4a.cpp ANRL/aflow_anrl_A2B_oC12_36_2a_a.cpp ANRL/aflow_anrl_A2BC_oC8_38_e_a_b.cpp ANRL/aflow_anrl_A2B_oC12_38_de_ab.cpp ANRL/aflow_anrl_AB4_oC20_41_a_2b.cpp ANRL/aflow_anrl_AB2_oC24_41_2a_2b.cpp ANRL/aflow_anrl_AB2_oF72_43_ab_3b.cpp ANRL/aflow_anrl_AB_oI4_44_a_b.cpp ANRL/aflow_anrl_A2B3C7D_oP13_47_t_aq_eqrs_h.cpp ANRL/aflow_anrl_AB_oP4_51_e_f.cpp ANRL/aflow_anrl_A3B2_oP20_56_ce_e.cpp ANRL/aflow_anrl_ABCD_oP16_57_d_c_d_d.cpp ANRL/aflow_anrl_AB_oP8_57_d_d.cpp ANRL/aflow_anrl_AB2_oP6_58_a_g.cpp ANRL/aflow_anrl_AB2_oP6_58_a_g.cpp ANRL/aflow_anrl_AB2_oP6_58_a_g.cpp ANRL/aflow_anrl_AB_oP4_59_a_b.cpp ANRL/aflow_anrl_ABC_oP6_59_a_a_a.cpp ANRL/aflow_anrl_A3B_oP8_59_bf_a.cpp ANRL/aflow_anrl_AB_oP16_61_c_c.cpp ANRL/aflow_anrl_A2B_oP24_61_2c_c.cpp ANRL/aflow_anrl_A3B2_oP20_62_3c_2c.cpp ANRL/aflow_anrl_AB3C_oP20_62_c_cd_a.cpp ANRL/aflow_anrl_A4B_oP20_62_2cd_c.cpp ANRL/aflow_anrl_AB2C_oP16_62_c_2c_c.cpp ANRL/aflow_anrl_A2B_oP12_62_2c_c.cpp ANRL/aflow_anrl_A2B_oP12_62_2c_c.cpp ANRL/aflow_anrl_A2B_oP12_62_2c_c.cpp ANRL/aflow_anrl_AB_oP8_62_c_c.cpp ANRL/aflow_anrl_AB_oP8_62_c_c.cpp ANRL/aflow_anrl_AB3_oP16_62_c_cd.cpp ANRL/aflow_anrl_A3B7_oP40_62_cd_3c2d.cpp ANRL/aflow_anrl_A_oP8_62_2c.cpp ANRL/aflow_anrl_AB_oP8_62_c_c.cpp ANRL/aflow_anrl_AB_oP8_62_c_c.cpp ANRL/aflow_anrl_AB2C_oC16_63_c_2c_c.cpp ANRL/aflow_anrl_A2B_oC12_63_2c_c.cpp ANRL/aflow_anrl_AB_oC8_63_c_c.cpp ANRL/aflow_anrl_A_oC4_63_c.cpp ANRL/aflow_anrl_A_oC8_64_f.cpp ANRL/aflow_anrl_A2B2C_oC80_64_efg_efg_df.cpp ANRL/aflow_anrl_A_oC8_64_f.cpp ANRL/aflow_anrl_A_oC8_64_f.cpp ANRL/aflow_anrl_AB_oC8_65_j_g.cpp ANRL/aflow_anrl_A3B5_oC16_65_ah_bej.cpp ANRL/aflow_anrl_AB3_oC8_65_a_bf.cpp ANRL/aflow_anrl_AB_oF8_69_a_b.cpp ANRL/aflow_anrl_A_oF8_70_a.cpp ANRL/aflow_anrl_A2B_oF24_70_e_a.cpp ANRL/aflow_anrl_A_oF128_70_4h.cpp ANRL/aflow_anrl_AB2_oI6_71_a_i.cpp ANRL/aflow_anrl_AB2_oI6_71_a_g.cpp ANRL/aflow_anrl_A2B_oI12_72_j_a.cpp ANRL/aflow_anrl_AB4C_tI12_82_c_g_a.cpp ANRL/aflow_anrl_A2BC4_tI14_82_bc_a_g.cpp ANRL/aflow_anrl_AB_tP16_84_cej_k.cpp ANRL/aflow_anrl_A4B5_tI18_87_h_ah.cpp ANRL/aflow_anrl_AB4_tI10_87_a_h.cpp ANRL/aflow_anrl_A2B_tP12_92_b_a.cpp ANRL/aflow_anrl_A2B_tP36_96_3b_ab.cpp ANRL/aflow_anrl_A_tP12_96_ab.cpp ANRL/aflow_anrl_A3BC_tP5_99_bc_a_b.cpp ANRL/aflow_anrl_AB3_tP8_113_a_ce.cpp ANRL/aflow_anrl_A2BC4D_tI16_121_d_a_i_b.cpp ANRL/aflow_anrl_ABC2_tI16_122_a_b_d.cpp ANRL/aflow_anrl_AB5C_tP7_123_b_ci_a.cpp ANRL/aflow_anrl_AB3_tP4_123_a_ce.cpp ANRL/aflow_anrl_AB_tP2_123_a_d.cpp ANRL/aflow_anrl_ABC2_tP4_123_d_a_f.cpp ANRL/aflow_anrl_A2B3_tP10_127_g_ah.cpp ANRL/aflow_anrl_ABCD_tP8_129_c_b_a_c.cpp ANRL/aflow_anrl_A_tP4_129_ac.cpp ANRL/aflow_anrl_ABC_tP6_129_c_a_c.cpp ANRL/aflow_anrl_A2B_tP6_129_ac_c.cpp ANRL/aflow_anrl_AB_tP4_129_a_c.cpp ANRL/aflow_anrl_AB_tP4_129_c_c.cpp ANRL/aflow_anrl_AB_tP4_131_c_e.cpp ANRL/aflow_anrl_A_tP50_134_b2m2n.cpp ANRL/aflow_anrl_A_tP30_136_bf2ij.cpp ANRL/aflow_anrl_AB_tP8_136_g_f.cpp ANRL/aflow_anrl_A2B_tP6_136_f_a.cpp ANRL/aflow_anrl_sigma_tP30_136_bf2ij.cpp ANRL/aflow_anrl_A_tP4_136_f.cpp ANRL/aflow_anrl_A_tP16_138_j.cpp ANRL/aflow_anrl_A3B_tI16_139_cde_e.cpp ANRL/aflow_anrl_A_tI4_139_e.cpp ANRL/aflow_anrl_AB2C4_tI14_139_a_e_ce.cpp ANRL/aflow_anrl_A12B_tI26_139_fij_a.cpp ANRL/aflow_anrl_A_tI2_139_a.cpp ANRL/aflow_anrl_A_tI8_139_h.cpp ANRL/aflow_anrl_A3B_tI8_139_bd_a.cpp ANRL/aflow_anrl_AB2_tI6_139_a_e.cpp ANRL/aflow_anrl_A4B5_tI18_139_i_ah.cpp ANRL/aflow_anrl_A4B_tI10_139_de_a.cpp ANRL/aflow_anrl_A8B_tI18_139_hi_a.cpp ANRL/aflow_anrl_A2B_tI6_139_d_a.cpp ANRL/aflow_anrl_A_tI2_139_a.cpp ANRL/aflow_anrl_A2B_tI12_140_h_a.cpp ANRL/aflow_anrl_AB3_tI16_140_b_ah.cpp ANRL/aflow_anrl_AB_tI16_140_ab_h.cpp ANRL/aflow_anrl_A4BC_tI24_141_h_b_a.cpp ANRL/aflow_anrl_A_tI4_141_a.cpp ANRL/aflow_anrl_A3B4_tI28_141_ad_h.cpp ANRL/aflow_anrl_A2B_tI12_141_e_a.cpp ANRL/aflow_anrl_AB_tI16_141_e_e.cpp ANRL/aflow_anrl_A2B_tI24_141_2e_e.cpp ANRL/aflow_anrl_AB_tI8_141_a_b.cpp ANRL/aflow_anrl_A2B3_tI80_141_ceh_3h.cpp ANRL/aflow_anrl_ABC4_tI96_142_e_ab_2g.cpp ANRL/aflow_anrl_A2B_hP9_147_g_ad.cpp ANRL/aflow_anrl_AB_hR16_148_cf_cf.cpp ANRL/aflow_anrl_AB3_hR8_148_c_f.cpp ANRL/aflow_anrl_AB_hR26_148_b2f_a2f.cpp ANRL/aflow_anrl_AB3C_hR10_148_c_f_c.cpp ANRL/aflow_anrl_A2B_hP9_150_ef_bd.cpp ANRL/aflow_anrl_A3B_hP24_151_3c_2a.cpp ANRL/aflow_anrl_A2B_hP9_152_c_a.cpp ANRL/aflow_anrl_A_hP3_152_a.cpp ANRL/aflow_anrl_AB_hP6_154_a_b.cpp ANRL/aflow_anrl_AB3_hR8_155_c_de.cpp ANRL/aflow_anrl_A3B2_hR5_155_e_c.cpp ANRL/aflow_anrl_AB_hR6_160_b_b.cpp ANRL/aflow_anrl_AB_hR6_160_3a_3a.cpp ANRL/aflow_anrl_ABC3_hR10_161_a_a_b.cpp ANRL/aflow_anrl_AB2_hP9_162_ad_k.cpp ANRL/aflow_anrl_AB2CD2_hP36_163_h_i_bf_i.cpp ANRL/aflow_anrl_A3B2_hP5_164_ad_d.cpp ANRL/aflow_anrl_AB2_hP3_164_a_d.cpp ANRL/aflow_anrl_A3B_hP24_165_adg_f.cpp ANRL/aflow_anrl_AB_hR2_166_a_b.cpp ANRL/aflow_anrl_A_hR2_166_c.cpp ANRL/aflow_anrl_A_hR1_166_a.cpp ANRL/aflow_anrl_A7B6_hR13_166_ah_3c.cpp ANRL/aflow_anrl_A_hR3_166_ac.cpp ANRL/aflow_anrl_A2B3_hR5_166_c_ac.cpp ANRL/aflow_anrl_A_hR1_166_a.cpp ANRL/aflow_anrl_A5B2_hR7_166_a2c_c.cpp ANRL/aflow_anrl_A_hR2_166_c.cpp ANRL/aflow_anrl_A_hR12_166_2h.cpp ANRL/aflow_anrl_ABC2_hR4_166_a_b_c.cpp ANRL/aflow_anrl_A_hR2_166_c.cpp ANRL/aflow_anrl_A_hR105_166_bc9h4i.cpp ANRL/aflow_anrl_A6B_hR7_166_g_a.cpp ANRL/aflow_anrl_ABC3_hR10_167_a_b_e.cpp ANRL/aflow_anrl_ABC3_hR10_167_a_b_e.cpp ANRL/aflow_anrl_A2B3_hR10_167_c_e.cpp ANRL/aflow_anrl_A2B_hP18_180_fi_bd.cpp ANRL/aflow_anrl_AB2_hP9_180_d_j.cpp ANRL/aflow_anrl_A2B_hP9_180_j_c.cpp ANRL/aflow_anrl_AB3_hP8_182_c_g.cpp ANRL/aflow_anrl_A_hP4_186_ab.cpp ANRL/aflow_anrl_AB_hP8_186_ab_ab.cpp ANRL/aflow_anrl_AB_hP4_186_b_b.cpp ANRL/aflow_anrl_AB_hP12_186_a2b_a2b.cpp ANRL/aflow_anrl_A5B3C_hP18_186_2a3b_2ab_b.cpp ANRL/aflow_anrl_AB_hP4_186_b_a.cpp ANRL/aflow_anrl_ABC_hP3_187_a_d_f.cpp ANRL/aflow_anrl_AB_hP2_187_d_a.cpp ANRL/aflow_anrl_A2B_hP9_189_fg_bc.cpp ANRL/aflow_anrl_AB4C_hP6_191_a_h_b.cpp ANRL/aflow_anrl_AB5_hP6_191_a_cg.cpp ANRL/aflow_anrl_A_hP1_191_a.cpp ANRL/aflow_anrl_A3B_hP4_191_bc_a.cpp ANRL/aflow_anrl_AB2_hP3_191_a_d.cpp ANRL/aflow_anrl_A2B_hP6_191_h_e.cpp ANRL/aflow_anrl_AB_hP6_191_f_ad.cpp ANRL/aflow_anrl_AB_hP8_194_ad_f.cpp ANRL/aflow_anrl_A_hP6_194_h.cpp ANRL/aflow_anrl_AB_hP12_194_af_bf.cpp ANRL/aflow_anrl_A_hP4_194_ac.cpp ANRL/aflow_anrl_AB3_hP8_194_c_bf.cpp ANRL/aflow_anrl_AB2_hP6_194_b_f.cpp ANRL/aflow_anrl_AB_hP4_194_c_d.cpp ANRL/aflow_anrl_ABC2_hP8_194_d_a_f.cpp ANRL/aflow_anrl_A3B_hP8_194_h_c.cpp ANRL/aflow_anrl_A_hP4_194_bc.cpp ANRL/aflow_anrl_AB2_hP6_194_c_f.cpp ANRL/aflow_anrl_A5B2_hP14_194_abdf_f.cpp ANRL/aflow_anrl_AB2_hP12_194_f_ah.cpp ANRL/aflow_anrl_ABC_hP6_194_c_d_a.cpp ANRL/aflow_anrl_A_hP4_194_f.cpp ANRL/aflow_anrl_AB2_hP6_194_c_ad.cpp ANRL/aflow_anrl_AB3C4_hP16_194_c_af_ef.cpp ANRL/aflow_anrl_A_hP2_194_c.cpp ANRL/aflow_anrl_AB2_hP24_194_ef_fgh.cpp ANRL/aflow_anrl_AB_hP12_194_df_ce.cpp ANRL/aflow_anrl_AB_hP4_194_c_a.cpp ANRL/aflow_anrl_A2B_hP12_194_cg_f.cpp ANRL/aflow_anrl_A4B_cI40_197_cde_c.cpp ANRL/aflow_anrl_ABC_cP12_198_a_a_a.cpp ANRL/aflow_anrl_A3B_cP16_198_b_a.cpp ANRL/aflow_anrl_A_cP8_198_2a.cpp ANRL/aflow_anrl_AB_cP8_198_a_a.cpp ANRL/aflow_anrl_AB_cP8_198_a_a.cpp ANRL/aflow_anrl_AB_cI16_199_a_a.cpp ANRL/aflow_anrl_AB32C48_cI162_204_a_2efg_2gh.cpp ANRL/aflow_anrl_A3B_cI32_204_g_c.cpp ANRL/aflow_anrl_A12B_cI26_204_g_a.cpp ANRL/aflow_anrl_A_cP8_205_c.cpp ANRL/aflow_anrl_AB_cP16_205_c_c.cpp ANRL/aflow_anrl_AB2_cP12_205_a_c.cpp ANRL/aflow_anrl_AB3C6_cI80_206_a_d_e.cpp ANRL/aflow_anrl_A_cI16_206_c.cpp ANRL/aflow_anrl_A_cP20_213_cd.cpp ANRL/aflow_anrl_A3B4C_cP8_215_d_e_a.cpp ANRL/aflow_anrl_AB4_cP5_215_a_e.cpp ANRL/aflow_anrl_AB3C4_cP8_215_a_c_e.cpp ANRL/aflow_anrl_AB5_cF24_216_a_ce.cpp ANRL/aflow_anrl_ABC_cF12_216_b_c_a.cpp ANRL/aflow_anrl_AB_cF8_216_c_a.cpp ANRL/aflow_anrl_A4B_cI10_217_c_a.cpp ANRL/aflow_anrl_A_cI58_217_ac2g.cpp ANRL/aflow_anrl_A5B8_cI52_217_ce_cg.cpp ANRL/aflow_anrl_A_cI16_220_c.cpp ANRL/aflow_anrl_A3B2_cI40_220_d_c.cpp ANRL/aflow_anrl_AB_cP2_221_b_a.cpp ANRL/aflow_anrl_AB_cP6_221_c_d.cpp ANRL/aflow_anrl_AB3C_cP5_221_a_c_b.cpp ANRL/aflow_anrl_AB27CD3_cP32_221_a_dij_b_c.cpp ANRL/aflow_anrl_AB3_cP4_221_a_c.cpp ANRL/aflow_anrl_A_cP1_221_a.cpp ANRL/aflow_anrl_AB11_cP36_221_c_agij.cpp ANRL/aflow_anrl_AB11CD3_cP16_221_a_dg_b_c.cpp ANRL/aflow_anrl_A3B_cP4_221_d_a.cpp ANRL/aflow_anrl_A6B_cP7_221_f_a.cpp ANRL/aflow_anrl_A3B_cP8_223_c_a.cpp ANRL/aflow_anrl_A_cP46_223_dik.cpp ANRL/aflow_anrl_A2B_cP6_224_b_a.cpp ANRL/aflow_anrl_A7B_cF32_225_bd_a.cpp ANRL/aflow_anrl_AB3_cF16_225_a_bc.cpp ANRL/aflow_anrl_A9B16C7_cF128_225_acd_2f_be.cpp ANRL/aflow_anrl_A12B_cF52_225_i_a.cpp ANRL/aflow_anrl_AB2_cF12_225_a_c.cpp ANRL/aflow_anrl_A6B23_cF116_225_e_acfh.cpp ANRL/aflow_anrl_AB2C_cF16_225_a_c_b.cpp ANRL/aflow_anrl_A_cF4_225_a.cpp ANRL/aflow_anrl_AB18C8_cF108_225_a_eh_f.cpp ANRL/aflow_anrl_AB_cF8_225_a_b.cpp ANRL/aflow_anrl_A2B_cF24_227_c_a.cpp ANRL/aflow_anrl_AB2_cF96_227_e_cf.cpp ANRL/aflow_anrl_AB_cF16_227_a_b.cpp ANRL/aflow_anrl_A_cF136_227_aeg.cpp ANRL/aflow_anrl_A2B_cF24_227_d_a.cpp ANRL/aflow_anrl_A_cF8_227_a.cpp ANRL/aflow_anrl_A2BC4_cF56_227_d_a_e.cpp ANRL/aflow_anrl_AB2_cF48_227_c_e.cpp ANRL/aflow_anrl_AB3C3_cF112_227_c_de_f.cpp ANRL/aflow_anrl_A_cI2_229_a.cpp ANRL/aflow_anrl_A3B_cI8_229_b_a.cpp ANRL/aflow_anrl_A4B3_cI14_229_c_b.cpp ANRL/aflow_anrl_A2B7_cI54_229_e_afh.cpp ANRL/aflow_anrl_AB12C3_cI32_229_a_h_b.cpp ANRL/aflow_anrl_AB4C3_cI16_229_a_c_b.cpp ANRL/aflow_anrl_A4B3_cI112_230_af_g.cpp ANRL/aflow_anrl_A2B_aP6_2_aei_i.cpp ANRL/aflow_anrl_A8B5_mP13_6_a7b_3a2b.cpp ANRL/aflow_anrl_AB_mP4_6_2b_2a.cpp ANRL/aflow_anrl_A2B_mP12_7_4a_2a.cpp ANRL/aflow_anrl_A2B_mP18_7_6a_3a.cpp ANRL/aflow_anrl_A3B_mP16_7_6a_2a.cpp ANRL/aflow_anrl_A9B2_mP22_7_9a_2a.cpp ANRL/aflow_anrl_A5B3_mC32_9_5a_3a.cpp ANRL/aflow_anrl_AB3_mC16_9_a_3a.cpp ANRL/aflow_anrl_A2B_mP6_10_mn_bg.cpp ANRL/aflow_anrl_AB3_mP16_10_mn_3m3n.cpp ANRL/aflow_anrl_ABC2_mP8_10_ac_eh_mn.cpp ANRL/aflow_anrl_AB_mP6_10_en_am.cpp ANRL/aflow_anrl_A_mP8_10_2m2n.cpp ANRL/aflow_anrl_A7B2C2_mC22_12_aij_h_i.cpp ANRL/aflow_anrl_A_mC16_12_4i.cpp ANRL/aflow_anrl_A2B_mP12_13_2g_ef.cpp ANRL/aflow_anrl_A2B_mP6_14_e_a.cpp ANRL/aflow_anrl_A7B8_mP120_14_14e_16e.cpp ANRL/aflow_anrl_AB3_mC16_15_e_cf.cpp ANRL/aflow_anrl_A_mC24_15_2e2f.cpp ANRL/aflow_anrl_A2B_oP12_17_abe_e.cpp ANRL/aflow_anrl_AB3_oP16_19_a_3a.cpp ANRL/aflow_anrl_AB2_oC6_21_a_k.cpp ANRL/aflow_anrl_A2BC2_oF40_22_fi_ad_gh.cpp ANRL/aflow_anrl_AB_oF8_22_a_c.cpp ANRL/aflow_anrl_A3B_oI32_23_ij2k_k.cpp ANRL/aflow_anrl_A8B2C12D2E_oI50_23_bcfk_i_3k_j_a.cpp ANRL/aflow_anrl_ABC2_oI16_23_ab_i_k.cpp ANRL/aflow_anrl_ABC4_oI12_23_a_b_k.cpp ANRL/aflow_anrl_AB7CD2_oI44_24_a_b3d_c_ac.cpp ANRL/aflow_anrl_A2B_oP12_26_abc_ab.cpp ANRL/aflow_anrl_A2B_oP12_26_abc_ab.cpp ANRL/aflow_anrl_A5B_oP24_26_3a3b2c_ab.cpp ANRL/aflow_anrl_A6B4C16D_oP108_27_abcd4e_4e_16e_e.cpp ANRL/aflow_anrl_A2B_oP12_29_2a_a.cpp ANRL/aflow_anrl_AB2_oP12_29_a_2a.cpp ANRL/aflow_anrl_ABC_oP12_29_a_a_a.cpp ANRL/aflow_anrl_A5B3C15_oP46_30_a2c_bc_a7c.cpp ANRL/aflow_anrl_ABC3_oP20_30_2a_c_3c.cpp ANRL/aflow_anrl_A13B2C2_oP34_32_a6c_c_c.cpp ANRL/aflow_anrl_A2B3_oP40_33_4a_6a.cpp ANRL/aflow_anrl_A2B8C_oP22_34_c_4c_a.cpp ANRL/aflow_anrl_AB2_oP6_34_a_c.cpp ANRL/aflow_anrl_AB8C2_oC22_35_a_ab3e_e.cpp ANRL/aflow_anrl_AB_oC8_36_a_a.cpp ANRL/aflow_anrl_A2B5C2_oC36_37_d_c2d_d.cpp ANRL/aflow_anrl_A2B3_oC40_39_2d_2c2d.cpp ANRL/aflow_anrl_A9BC_oC44_39_3c3d_a_c.cpp ANRL/aflow_anrl_AB2C_oC16_40_a_2b_b.cpp ANRL/aflow_anrl_AB3_oC16_40_b_3b.cpp ANRL/aflow_anrl_A10B3_oF52_42_2abce_ab.cpp ANRL/aflow_anrl_AB_oF8_42_a_a.cpp ANRL/aflow_anrl_A2BC2_oI20_45_c_b_c.cpp ANRL/aflow_anrl_ABC_oI36_46_ac_bc_3b.cpp ANRL/aflow_anrl_A2B8CD_oP24_48_k_2m_d_b.cpp ANRL/aflow_anrl_A5B2_oP14_49_dehq_ab.cpp ANRL/aflow_anrl_AB2C8D_oP24_49_g_q_2qr_e.cpp ANRL/aflow_anrl_A2BC4_oP28_50_ij_ac_ijm.cpp ANRL/aflow_anrl_A3BC2_oP48_50_3m_m_2m.cpp ANRL/aflow_anrl_A2B_oP24_52_2e_cd.cpp ANRL/aflow_anrl_A3B2_oP20_52_de_cd.cpp ANRL/aflow_anrl_ABC2_oP16_53_h_e_gh.cpp ANRL/aflow_anrl_ABC3_oP20_53_e_g_hi.cpp ANRL/aflow_anrl_ABC3_oP20_54_e_d_cf.cpp ANRL/aflow_anrl_A2B_oP24_55_2g2h_gh.cpp ANRL/aflow_anrl_A3B5_oP16_55_ch_agh.cpp ANRL/aflow_anrl_A_oP16_55_2g2h.cpp ANRL/aflow_anrl_A2B_oP6_58_g_a.cpp ANRL/aflow_anrl_ABC_oP6_59_a_b_a.cpp ANRL/aflow_anrl_A2B3_oP20_60_d_cd.cpp ANRL/aflow_anrl_A3B_oP32_60_3d_d.cpp ANRL/aflow_anrl_A7B8_oP120_60_7d_8d.cpp ANRL/aflow_anrl_AB_oP48_61_3c_3c.cpp ANRL/aflow_anrl_A2B3_oP20_62_2c_3c.cpp ANRL/aflow_anrl_A2B4C_oP28_62_ac_2cd_c.cpp ANRL/aflow_anrl_A2B_oP12_62_2c_c.cpp ANRL/aflow_anrl_A3B_oP16_62_cd_c.cpp ANRL/aflow_anrl_AB2C3_oP24_62_c_d_cd.cpp ANRL/aflow_anrl_AB3_oP16_62_c_3c.cpp ANRL/aflow_anrl_AB4C_oP24_62_c_2cd_c.cpp ANRL/aflow_anrl_AB_oP8_62_c_c.cpp ANRL/aflow_anrl_A2BC3_oC24_63_e_c_cg.cpp ANRL/aflow_anrl_A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h.cpp ANRL/aflow_anrl_A6B_oC28_63_efg_c.cpp ANRL/aflow_anrl_AB3C_oC20_63_a_cf_c.cpp ANRL/aflow_anrl_AB4C_oC24_63_a_fg_c.cpp ANRL/aflow_anrl_AB4C_oC24_63_c_fg_c.cpp ANRL/aflow_anrl_A2B_oC24_64_2f_f.cpp ANRL/aflow_anrl_A2B4C_oC28_66_l_kl_a.cpp ANRL/aflow_anrl_A3B_oC64_66_gi2lm_2l.cpp ANRL/aflow_anrl_A3B_oC64_66_kl2m_bdl.cpp ANRL/aflow_anrl_A2BC_oC16_67_ag_b_g.cpp ANRL/aflow_anrl_ABC2_oC16_67_b_g_ag.cpp ANRL/aflow_anrl_AB_oC8_67_a_g.cpp ANRL/aflow_anrl_AB_oC8_67_a_g.cpp ANRL/aflow_anrl_AB4_oC20_68_a_i.cpp ANRL/aflow_anrl_AB2_oF48_70_f_fg.cpp ANRL/aflow_anrl_A4B3_oI14_71_gh_cg.cpp ANRL/aflow_anrl_ABC_oI12_71_h_j_g.cpp ANRL/aflow_anrl_ABCD3_oI48_73_d_e_e_ef.cpp ANRL/aflow_anrl_A2B_oI12_74_h_e.cpp ANRL/aflow_anrl_A4B_oI20_74_beh_e.cpp ANRL/aflow_anrl_AB2C12D4_tP76_75_2a2b_2d_12d_4d.cpp ANRL/aflow_anrl_A2BC_tP16_76_2a_a_a.cpp ANRL/aflow_anrl_A3B7_tP40_76_3a_7a.cpp ANRL/aflow_anrl_A2B6CD7_tP64_77_2d_6d_d_ab6d.cpp ANRL/aflow_anrl_A2B_tP48_77_8d_4d.cpp ANRL/aflow_anrl_A2B7C2_tP88_78_4a_14a_4a.cpp ANRL/aflow_anrl_A2BC2_tI20_79_c_2a_c.cpp ANRL/aflow_anrl_AB2_tI48_80_2b_4b.cpp ANRL/aflow_anrl_AB2_tP12_81_adg_2h.cpp ANRL/aflow_anrl_A3B_tI32_82_3g_g.cpp ANRL/aflow_anrl_A3B2_tP10_83_adk_j.cpp ANRL/aflow_anrl_A2B_tP30_85_ab2g_cg.cpp ANRL/aflow_anrl_AB3_tP32_86_g_3g.cpp ANRL/aflow_anrl_A4B_tI20_88_f_a.cpp ANRL/aflow_anrl_AB2_tI96_88_2f_4f.cpp ANRL/aflow_anrl_A17BC4D_tP184_89_17p_p_4p_io.cpp ANRL/aflow_anrl_A4B2C13D_tP40_90_g_d_cef2g_c.cpp ANRL/aflow_anrl_AB4C17D4E_tP54_90_a_g_c4g_g_c.cpp ANRL/aflow_anrl_ABC_tP24_91_d_d_d.cpp ANRL/aflow_anrl_AB32CD4E8_tP184_93_i_16p_af_2p_4p.cpp ANRL/aflow_anrl_A14B3C5_tP44_94_c3g_ad_bg.cpp ANRL/aflow_anrl_A6B2C_tP18_94_eg_c_a.cpp ANRL/aflow_anrl_ABC_tP24_95_d_d_d.cpp ANRL/aflow_anrl_A2B8CD_tI24_97_d_k_a_b.cpp ANRL/aflow_anrl_AB8C2_tI44_97_e_2k_cd.cpp ANRL/aflow_anrl_A2B_tI12_98_f_a.cpp ANRL/aflow_anrl_A2B8C2D_tP26_100_c_abcd_c_a.cpp ANRL/aflow_anrl_A3B11C6_tP40_100_ac_bc2d_cd.cpp ANRL/aflow_anrl_A7B7C2_tP32_101_bde_ade_d.cpp ANRL/aflow_anrl_A2B3_tP20_102_2c_b2c.cpp ANRL/aflow_anrl_AB4_tP10_103_a_d.cpp ANRL/aflow_anrl_A5B5C4_tP28_104_ac_ac_c.cpp ANRL/aflow_anrl_AB6C4_tP22_104_a_2ac_c.cpp ANRL/aflow_anrl_A2BC2_tP20_105_f_ac_2e.cpp ANRL/aflow_anrl_A3BC3D_tP64_106_3c_c_3c_c.cpp ANRL/aflow_anrl_A5B7_tI24_107_ac_abd.cpp ANRL/aflow_anrl_AB_tI4_107_a_a.cpp ANRL/aflow_anrl_A3B5_tI32_108_ac_a2c.cpp ANRL/aflow_anrl_ABC_tI12_109_a_a_a.cpp ANRL/aflow_anrl_AB_tI8_109_a_a.cpp ANRL/aflow_anrl_A2BC8_tI176_110_2b_b_8b.cpp ANRL/aflow_anrl_A2B_tP12_111_2n_adf.cpp ANRL/aflow_anrl_AB_tP8_111_n_n.cpp ANRL/aflow_anrl_AB4C_tP12_112_b_n_e.cpp ANRL/aflow_anrl_A2BC7D2_tP24_113_e_a_cef_e.cpp ANRL/aflow_anrl_A3B_tP32_114_3e_e.cpp ANRL/aflow_anrl_A4B_tP10_114_e_a.cpp ANRL/aflow_anrl_A2B3_tP5_115_g_ag.cpp ANRL/aflow_anrl_AB2_tP12_115_j_egi.cpp ANRL/aflow_anrl_A2B3_tP20_116_bci_fj.cpp ANRL/aflow_anrl_A2B3_tP20_117_i_adgh.cpp ANRL/aflow_anrl_A3B_tP16_118_ei_f.cpp ANRL/aflow_anrl_A5B3_tP32_118_g2i_aceh.cpp ANRL/aflow_anrl_A3B_tI24_119_b2i_af.cpp ANRL/aflow_anrl_AB_tI4_119_c_a.cpp ANRL/aflow_anrl_A4BC2_tI28_120_i_d_e.cpp ANRL/aflow_anrl_A4BC4D_tP10_123_gh_a_i_d.cpp ANRL/aflow_anrl_AB4C_tP12_124_a_m_c.cpp ANRL/aflow_anrl_AB4_tP10_124_a_m.cpp ANRL/aflow_anrl_A4B_tP10_125_m_a.cpp ANRL/aflow_anrl_ABC4_tP12_125_a_b_m.cpp ANRL/aflow_anrl_A2BC4_tP28_126_cd_e_k.cpp ANRL/aflow_anrl_A4B_tP20_127_ehj_g.cpp ANRL/aflow_anrl_A6B2C_tP18_128_eh_d_a.cpp ANRL/aflow_anrl_A7B2C_tP40_128_egi_h_e.cpp ANRL/aflow_anrl_A2BC4_tP28_130_f_c_g.cpp ANRL/aflow_anrl_A5B3_tP32_130_cg_cf.cpp ANRL/aflow_anrl_A2B2C4D_tP18_132_e_i_o_d.cpp ANRL/aflow_anrl_AB6C_tP16_132_d_io_a.cpp ANRL/aflow_anrl_AB3_tP32_133_h_i2j.cpp ANRL/aflow_anrl_A2B_tP24_135_gh_h.cpp ANRL/aflow_anrl_A4B2C_tP28_135_gh_h_d.cpp ANRL/aflow_anrl_A2B3_tP40_137_cdf_3g.cpp ANRL/aflow_anrl_A2B_tP6_137_d_a.cpp ANRL/aflow_anrl_A4BC4_tP18_137_g_b_g.cpp ANRL/aflow_anrl_AB2_tP6_137_a_d.cpp ANRL/aflow_anrl_A_tP12_138_bi.cpp ANRL/aflow_anrl_AB_tI8_139_e_e.cpp ANRL/aflow_anrl_A3B5_tI32_140_ah_bk.cpp ANRL/aflow_anrl_A3B5_tI32_140_ah_cl.cpp ANRL/aflow_anrl_A2B_tI12_141_e_a.cpp ANRL/aflow_anrl_A_tI16_142_f.cpp ANRL/aflow_anrl_A4B14C3_hP21_143_bd_ac4d_d.cpp ANRL/aflow_anrl_A4B6C_hP11_143_bd_2d_a.cpp ANRL/aflow_anrl_AB2_hP12_143_cd_ab2d.cpp ANRL/aflow_anrl_A4B_hP15_144_4a_a.cpp ANRL/aflow_anrl_AB_hP6_144_a_a.cpp ANRL/aflow_anrl_A2B3C3DE7_hP48_145_2a_3a_3a_a_7a.cpp ANRL/aflow_anrl_A3BC_hR5_146_b_a_a.cpp ANRL/aflow_anrl_ABC3_hR10_146_2a_2a_2b.cpp ANRL/aflow_anrl_A2B4C_hR42_148_2f_4f_f.cpp ANRL/aflow_anrl_A2B_hR18_148_2f_f.cpp ANRL/aflow_anrl_AB3_hP24_149_acgi_3l.cpp ANRL/aflow_anrl_A3B_hP24_153_3c_2b.cpp ANRL/aflow_anrl_A_hP9_154_bc.cpp ANRL/aflow_anrl_AB2_hP9_156_b2c_3a2bc.cpp ANRL/aflow_anrl_AB_hP12_156_2ab3c_2ab3c.cpp ANRL/aflow_anrl_AB_hP4_156_ab_ab.cpp ANRL/aflow_anrl_A5B6C2_hP13_157_2ac_2c_b.cpp ANRL/aflow_anrl_A3B_hP8_158_d_a.cpp ANRL/aflow_anrl_A2B3_hP20_159_bc_2c.cpp ANRL/aflow_anrl_A4B3_hP28_159_ab2c_2c.cpp ANRL/aflow_anrl_AB4C7D_hP26_159_b_ac_a2c_b.cpp ANRL/aflow_anrl_A3B_hR4_160_b_a.cpp ANRL/aflow_anrl_A8B5_hR26_160_a3bc_a3b.cpp ANRL/aflow_anrl_ABC_hR3_160_a_a_a.cpp ANRL/aflow_anrl_AB_hR10_160_5a_5a.cpp ANRL/aflow_anrl_A2B3_hP5_164_d_ad.cpp ANRL/aflow_anrl_AB2_hP9_164_bd_c2d.cpp ANRL/aflow_anrl_ABC2_hP4_164_a_b_d.cpp ANRL/aflow_anrl_A3B_hP24_165_bdg_f.cpp ANRL/aflow_anrl_A4B3_hR7_166_2c_ac.cpp ANRL/aflow_anrl_ABC_hR6_166_c_c_c.cpp ANRL/aflow_anrl_AB3C_hR10_167_b_e_a.cpp ANRL/aflow_anrl_ABC2_hR24_167_e_e_2e.cpp ANRL/aflow_anrl_A2B13C4_hP57_168_d_c6d_2d.cpp ANRL/aflow_anrl_AB4C_hP72_168_2d_8d_2d.cpp ANRL/aflow_anrl_A2B3_hP30_169_2a_3a.cpp ANRL/aflow_anrl_A2B3_hP30_170_2a_3a.cpp ANRL/aflow_anrl_A10B2C_hP39_171_5c_c_a.cpp ANRL/aflow_anrl_A10B2C_hP39_172_5c_c_a.cpp ANRL/aflow_anrl_A3B_hP8_173_c_b.cpp ANRL/aflow_anrl_A4B3_hP14_173_bc_c.cpp ANRL/aflow_anrl_A12B7C2_hP21_174_2j2k_ajk_cf.cpp ANRL/aflow_anrl_ABC_hP12_174_cj_fk_aj.cpp ANRL/aflow_anrl_A8B7C6_hP21_175_ck_aj_k.cpp ANRL/aflow_anrl_ABC_hP36_175_jk_jk_jk.cpp ANRL/aflow_anrl_A3B2_hP10_176_h_bc.cpp ANRL/aflow_anrl_A3B3C_hP14_176_h_h_c.cpp ANRL/aflow_anrl_A3B_hP8_176_h_c.cpp ANRL/aflow_anrl_A2B_hP36_177_j2lm_n.cpp ANRL/aflow_anrl_AB3_hP24_178_b_ac.cpp ANRL/aflow_anrl_A_hP6_178_a.cpp ANRL/aflow_anrl_AB3_hP24_179_b_ac.cpp ANRL/aflow_anrl_A2B_hP9_181_j_c.cpp ANRL/aflow_anrl_ABC_hP3_183_a_a_a.cpp ANRL/aflow_anrl_AB_hP6_183_c_ab.cpp ANRL/aflow_anrl_AB4C_hP72_184_d_4d_d.cpp ANRL/aflow_anrl_A3BC_hP30_185_cd_c_ab.cpp ANRL/aflow_anrl_A3B_hP24_185_ab2c_c.cpp ANRL/aflow_anrl_A3B_hP8_185_c_a.cpp ANRL/aflow_anrl_AB3_hP24_185_c_ab2c.cpp ANRL/aflow_anrl_A3B7_hP20_186_c_b2c.cpp ANRL/aflow_anrl_AB3_hP4_187_e_fh.cpp ANRL/aflow_anrl_A3BC_hP10_188_k_c_a.cpp ANRL/aflow_anrl_AB9C4_hP28_188_e_kl_ak.cpp ANRL/aflow_anrl_A8BC3D6_hP18_189_bfh_a_g_i.cpp ANRL/aflow_anrl_A9BC3D5_hP18_189_fi_a_g_bh.cpp ANRL/aflow_anrl_A2B_hP18_190_gh_bf.cpp ANRL/aflow_anrl_A5B3_hP16_190_bdh_g.cpp ANRL/aflow_anrl_AB_hP24_190_i_afh.cpp ANRL/aflow_anrl_A2B3C18D6_hP58_192_c_f_lm_l.cpp ANRL/aflow_anrl_AB2_hP72_192_m_j2kl.cpp ANRL/aflow_anrl_A5B3_hP16_193_dg_g.cpp ANRL/aflow_anrl_A3B_hP16_194_gh_ac.cpp ANRL/aflow_anrl_A5B2_hP28_194_ahk_ch.cpp ANRL/aflow_anrl_A9B3C_hP26_194_hk_h_a.cpp ANRL/aflow_anrl_A12BC4_cP34_195_2j_ab_2e.cpp ANRL/aflow_anrl_A12B2C_cF60_196_h_bc_a.cpp ANRL/aflow_anrl_ABC3_cP20_198_a_a_b.cpp ANRL/aflow_anrl_A2B11_cP39_200_f_aghij.cpp ANRL/aflow_anrl_AB3C_cP60_201_be_fh_g.cpp ANRL/aflow_anrl_A6B6C_cF104_202_h_h_c.cpp ANRL/aflow_anrl_A_cF240_202_h2i.cpp ANRL/aflow_anrl_A2BCD3E6_cF208_203_e_c_d_f_g.cpp ANRL/aflow_anrl_A4B2C6D16E_cF232_203_e_d_f_eg_a.cpp ANRL/aflow_anrl_AB3C16_cF160_203_a_bc_eg.cpp ANRL/aflow_anrl_A2B3C6_cP264_205_2d_ab2c2d_6d.cpp ANRL/aflow_anrl_A_cP240_205_10d.cpp ANRL/aflow_anrl_AB3C2_cI96_206_c_e_ad.cpp ANRL/aflow_anrl_A17B15_cP64_207_acfk_eij.cpp ANRL/aflow_anrl_A3B_cP16_208_j_b.cpp ANRL/aflow_anrl_A6B2CD6E_cP64_208_m_ad_b_m_c.cpp ANRL/aflow_anrl_A24BC_cF104_209_j_a_b.cpp ANRL/aflow_anrl_A12B36CD12_cF488_210_h_3h_a_fg.cpp ANRL/aflow_anrl_A12B6C_cF608_210_4h_2h_e.cpp ANRL/aflow_anrl_A2B_cI72_211_hi_i.cpp ANRL/aflow_anrl_A2B_cP12_212_c_a.cpp ANRL/aflow_anrl_A3B3C_cI56_214_g_h_a.cpp ANRL/aflow_anrl_A3BC2_cI48_214_f_a_e.cpp ANRL/aflow_anrl_A4B9_cP52_215_ei_3efgi.cpp ANRL/aflow_anrl_ABCD_cF16_216_c_d_b_a.cpp ANRL/aflow_anrl_A3B4C_cP16_218_c_e_a.cpp ANRL/aflow_anrl_A7BC3D13_cF192_219_de_b_c_ah.cpp ANRL/aflow_anrl_A15B4_cI76_220_ae_c.cpp ANRL/aflow_anrl_A4B3_cI28_220_c_a.cpp ANRL/aflow_anrl_A2B3C6_cP33_221_cd_ag_fh.cpp ANRL/aflow_anrl_A5B3C16_cP96_222_ce_d_fi.cpp ANRL/aflow_anrl_A23B6_cF116_225_bd2f_e.cpp ANRL/aflow_anrl_A6B2C_cF36_225_e_c_a.cpp ANRL/aflow_anrl_AB13_cF112_226_a_bi.cpp ANRL/aflow_anrl_A2B2C7_cF88_227_c_d_af.cpp ANRL/aflow_anrl_A3B4_cF56_227_ad_e.cpp ANRL/aflow_anrl_A5BCD6_cF416_228_eg_c_b_h.cpp ANRL/aflow_anrl_A6B_cF224_228_h_c.cpp ANRL/aflow_anrl_A3B10_cI52_229_e_fh.cpp ANRL/aflow_anrl_A4B_cI10_229_c_a.cpp ANRL/aflow_anrl_A7B3_cI40_229_df_e.cpp ANRL/aflow_anrl_A2B3C12D3_cI160_230_a_c_h_d.cpp ANRL/aflow_anrl_AB_aP16_2_4i_4i.cpp ANRL/aflow_anrl_A5B11_mP16_6_2abc_2a3b3c.cpp ANRL/aflow_anrl_AB3_mC32_8_4a_12a.cpp ANRL/aflow_anrl_AB3_mC32_8_4a_4a4b.cpp ANRL/aflow_anrl_A3B13_oC32_38_ac_a2bcdef.cpp ANRL/aflow_anrl_A3B5_oC32_38_abce_abcdf.cpp ANRL/aflow_anrl_AB7_hR16_166_c_c2h.cpp ANRL/aflow_anrl_A2BCD4_tI16_82_ac_b_d_g.cpp 
COMPARE_FUNCS_DH_CPP = aflow_compare_structure.cpp aflow_compare_structure_function.cpp
CHULL_FUNCS_CPP = aflow_chull.cpp
#MB 190301 - added jupyter notebook files for chull
CHULL_JUPYTER_CPP=aflow_chull_jupyter.cpp aflow_chull_jupyter_plotter.cpp aflow_chull_python.cpp aflow_chull_jupyter_requirements.cpp
GFA_FUNCS_CPP = aflow_gfa.cpp
#RF 190603
CCE_CPP = aflow_cce.cpp
#CO 170614
MAKEFILE_CPP = aflow_makefile.cpp
SQLITE_C = ./SQLITE/sqlite3.c
#README_CPP= aflow_readme.cpp
DATA_CPP=aflow_data_readme.cpp aflow_data_calculated.cpp aflow_data_htqc.cpp aflow_data_libraries.cpp aflow_data_stokes.cpp aflow_data_nist.cpp aflow_data_latex.cpp aflow_data_aflow_potcars.cpp aflow_data.cpp
#aflow_data aflow_data.o aflow_data_aflow_potcars.o aflow_data_calculated.o aflow_data_htqc.o aflow_data_latex.o aflow_data_libraries.o aflow_data_nist.o aflow_data_readme.o aflow_data_stokes.o

AFLOW_OBJ              = $(AFLOW_CPP:.cpp=.o)
AFLOW_CLASS_OBJ        = $(AFLOW_CLASS_CPP:.cpp=.o)
AFLOW_FUNCS_OBJ        = $(AFLOW_FUNCS_CPP:.cpp=.o)
PFLOW_FUNCS_OBJ        = $(PFLOW_FUNCS_CPP:.cpp=.o)
ESTRUCTURE_FUNCS_OBJ   = $(ESTRUCTURE_FUNCS_CPP:.cpp=.o)
POCC_FUNCS_OLD_OBJ     = $(POCC_FUNCS_OLD_CPP:.cpp=.o)
POCC_OBJ               = $(POCC_FUNCS_CPP:.cpp=.o) #CO 180409
APENNSY_FUNCS_OBJ      = $(APENNSY_FUNCS_CPP:.cpp=.o)
LIBRARIES_FUNCS_OBJ    = $(LIBRARIES_FUNCS_CPP:.cpp=.o)
MIX_EXPERIMENTS_OBJ    = $(MIX_EXPERIMENTS_CPP:.cpp=.o)
CONTRIB_FUNCS_OBJ      = $(CONTRIB_FUNCS_SW_CPP:.cpp=.o) $(CONTRIB_FUNCS_JX_CPP:.cpp=.o) $(CONTRIB_FUNCS_KY_CPP:.cpp=.o) $(CONTRIB_FUNCS_RT_CPP:.cpp=.o) $(CONTRIB_FUNCS_WS_CPP:.cpp=.o) $(CONTRIB_FUNCS_SC_CPP:.cpp=.o) $(CONTRIB_FUNCS_KR_CPP:.cpp=.o) $(CONTRIB_FUNCS_CT_CPP:.cpp=.o) 
#README_OBJ 						 = $(README_CPP:.cpp=.o)
AFLOWLIB_OBJ           = $(AFLOWLIB_CPP:.cpp=.o)
AGL_OBJ                = $(AGL_CPP:.cpp=.o)
AEL_OBJ                = $(AEL_CPP:.cpp=.o)
APL_OBJ                = $(APL_CPP:.cpp=.o)
ANRL_OBJ               = $(ANRL_CPP:.cpp=.o)
COMPARE_FUNCS_DH_OBJ   = $(COMPARE_FUNCS_DH_CPP:.cpp=.o)
CHULL_OBJ              = $(CHULL_FUNCS_CPP:.cpp=.o)
GFA_OBJ                = $(GFA_FUNCS_CPP:.cpp=.o)
#RF 190603
CCE_OBJ                = $(CCE_CPP:.cpp=.o)
#CO 170614
MAKEFILE_OBJ           = $(MAKEFILE_CPP:.cpp=.o)
SQLITE_OBJ             = $(SQLITE_C:.c=.o)

AFLOW_OBJS = $(AFLOW_OBJ) $(AUROSTD_OBJ) $(AFLOW_CLASS_OBJ) $(AFLOW_FUNCS_OBJ) $(PFLOW_FUNCS_OBJ) $(ESTRUCTURE_FUNCS_OBJ) $(POCC_FUNCS_OLD_OBJ) $(POCC_OBJ) $(APENNSY_FUNCS_OBJ) $(LIBRARIES_FUNCS_OBJ) $(MIX_EXPERIMENTS_OBJ) $(CONTRIB_FUNCS_OBJ) $(AFLOWLIB_OBJ) $(AGL_OBJ) $(AEL_OBJ) $(APL_OBJ) $(ANRL_OBJ) $(COMPARE_FUNCS_DH_OBJ) $(CHULL_OBJ) $(GFA_OBJ) $(CCE_OBJ) $(MAKEFILE_OBJ) $(SQLITE_OBJ)
AFLOW_DEP = $(AFLOW_OBJS) $(AFLOWLIB_WEBAPP_CPP) $(CHULL_JUPYTER_CPP)

all: aflow aflow_data 
	$(MAKE) check
aflow: $(AFLOW_DEP) Makefile 
#	$(MAKE) -C ../AFLOW2_FREE/EXTRA/APL2AGR/ && cd ../../
#	cp -f aflow_apl*.o  aflow_aapl*.o  aflow_qha*.o ./APL/
#	cp -f aflow_anrl*.o ./ANRL/
#	cp sqlite3.o ./SQLITE/
#	$(CPP) $(LIBS_AFLOW) -o aflow *.o $(LIB_MKL) $(LIBS_AFLOW)
#	$(CPP) $(LIBS_AFLOW) -o aflow `ls *.o | grep -v aflow_data` $(LIB_MKL) $(LIBS_AFLOW)
	@echo "creating single aflow binary"
	@$(CPP) $(LIBS_AFLOW) -o aflow $(AFLOW_OBJS) $(LIB_MKL) $(LIBS_AFLOW)
	ln -sf aflow aflowd
	ln -sf aflow apennsy

$(SQLITE_OBJ): $(SQLITE_C)	#CO 1706147
	#$(CC) -DSQLITE_ENABLE_JSON1 -DSQLITE_OMIT_LOAD_EXTENSION -o $@ -c $^ $(LIBS)	#-DSQLITE_OMIT_LOAD_EXTENSION doesn't really work amalgamation
	#$(CC) -DSQLITE_ENABLE_JSON1 -o $@ -c $^ $(LIBS)
	#$(CC) -DSQLITE_ENABLE_JSON1 -c $^ $(LIBS)	#CO190401 - $(LIBS) -lpthreads not needed, issues with clang++ (mac)
	$(CC) -DSQLITE_ENABLE_JSON1 -c $^ -o $@
	#cp sqlite3.o ./SQLITE/

.PHONY: extra extra_clean
extra:	apl2agr subst findsym frozsl platon
apl2agr:
	echo "MAKING APL2AGR"
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR realclean; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/apl2agr ] && [ -d ../AFLOW3_NONFREE/BIN ]; then cp $(EXTRA_FREE)/APL2AGR/apl2agr ../AFLOW3_NONFREE/BIN; fi
subst:
	echo "MAKING SUBST"
	if [ -f $(EXTRA_FREE)/SUBST/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/SUBST clean; fi
	if [ -f $(EXTRA_FREE)/SUBST/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/SUBST; fi
	if [ -f $(EXTRA_FREE)/SUBST/subst ] && [ -d ../AFLOW3_NONFREE/BIN ]; then cp $(EXTRA_FREE)/SUBST/subst ../AFLOW3_NONFREE/BIN; fi
findsym:
	echo "MAKING FINDSYM"
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym ]; then cp ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym ../AFLOW3_NONFREE/BIN; fi
frozsl:
	echo "MAKING FROZSL"
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl ]; then cp ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl ../AFLOW3_NONFREE/BIN; fi
platon:
	echo "MAKING PLATON"
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/PLATON/51108 clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/PLATON/51108; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/platon ]; then cp ../AFLOW3_NONFREE/EXTRA/PLATON/51108/platon ../AFLOW3_NONFREE/BIN; fi
enum:
	echo "MAKING ENUM"
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS; fi
extra_clean:
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR realclean; fi
	if [ -f $(EXTRA_FREE)/SUBST/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/SUBST clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/PLATON/51108 clean; fi
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi

.PHONY: mac
mac:
#ME190313 - Change compiler to clang because g++ fails to compile on some Macs
#	[CO191217 - fix above for the OS]$(MAKE) CC=clang CPP=clang++ aflow_data
#	[CO191217 - fix above for the OS]$(MAKE) CC=clang CPP=clang++ -e ARCH=-D_MACOSX_ -j 4
	$(MAKE) aflow_data
	$(MAKE) -e ARCH=-D_MACOSX_ -j 4
#	ln -sf aflow aflowd
#	ln -sf aflow apennsy
	$(MAKE) check
#mac_static: mac
#	tar Jcfvp aflow_$(VERSNUMBER)_maxOSX_10.X.X.tar.xz aflow aflow_data

.PHONY: static_no_pdf
static_no_pdf: aflow aflow_data
	rm -rf $(SQLITE_OBJ)
	#$(CC) -DSQLITE_ENABLE_JSON1 -c $(SQLITE_C) $(LIBS) $(LIBS_STATIC)	#CO190401 - $(LIBS) -lpthreads not needed, issues with clang++ (mac)
	$(CC) -DSQLITE_ENABLE_JSON1 -c $(SQLITE_C) $(LIBS_STATIC)
	cp sqlite3.o ./SQLITE/
	rm -f aflow
#	$(CPP) $(LIBS_AFLOW_STATIC) -o aflow *.o $(LIB_MKL) $(LIBS_AFLOW_STATIC)
#	$(CPP) $(LIBS_AFLOW_STATIC) -o aflow `ls *.o | grep -v aflow_data` $(LIB_MKL) $(LIBS_AFLOW_STATIC)
	@echo "creating single aflow binary"
	@$(CPP) $(LIBS_AFLOW_STATIC) -o aflow $(AFLOW_OBJS) $(LIB_MKL) $(LIBS_AFLOW_STATIC)
	rm -f aflow_data
#	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -o aflow_data aflow_data.cpp $(LIBS_STATIC)
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -o aflow_data aflow_data.o aflow_data_*.o $(LIBS_STATIC)
#static: static_no_pdf
#	tar Jcfvp /www/AFLOW/CURRENT/aflow_$(VERSNUMBER)_amd64_debian_static.tar.xz aflow aflow_data	

$(AFLOW_OBJ): $(AFLOW_CPP) aflow_test.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $(AFLOW_CPP) -c -o $@
#$(AUROSTD_OBJ): aurostd.cpp aurostd.h Makefile
#	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) aurostd.cpp -c -o $(AUROSTD_OBJ)
#aurostd.cpp:
#	cd ./AUROSTD && ppp4cpp.py
#aurostd.h:
#	cd ./AUROSTD && ppp4cpp.py

$(AUROSTD_OBJ): $(AUROSTD_HPP) $(AUROSTD_CPP) Makefile
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) ./AUROSTD/aurostd.cpp -c -o $@
	touch $@
.cpp.o: .cpp Makefile $(AFLOW_HPP) $(AUROSTD_HPP) $(AUROSTD_CPP)
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@
	touch $@
aflow_init.o: aflow_init.cpp Makefile $(AFLOW_HPP) $(AUROSTD_HPP) $(AUROSTD_CPP)
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@
aflow_apennsy_phases.o: aflow_apennsy_phases.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $< -c -o $@
#aflow_mix_nomix_table.o: aflow_mix_nomix_table.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_mix_nomix_table.cpp -c
#CO 170622 - START
aflowlib_webapp_entry.cpp: aflowlib_webapp_entry.js
	rm -rf $@
	echo "// $@ automatic generated" >> $@
	echo "std::string AFLOW_WEBAPP_ENTRY_JS=\"\\" >> $@
	cat $^ | perl -p -e 's/\r//g' | perl -p -e 's/\\u/\\\\u/g' | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
#CO 170622 - END
#CO 180305 - START
aflowlib_webapp_bands.cpp: aflowlib_webapp_bands.js
	rm -rf $@
	echo "// $@ automatic generated" >> $@
	echo "std::string AFLOW_WEBAPP_BANDS_JS=\"\\" >> $@
	cat $^ | perl -p -e 's/\r//g' | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
#CO 180305 - END
aflow_chull.o: aflow_chull.cpp Makefile $(AFLOW_HPP) $(AUROSTD_HPP) $(AUROSTD_CPP) aflow_chull_jupyter.cpp aflow_chull_jupyter_plotter.cpp aflow_chull_jupyter_requirements.cpp aflow_chull_python.cpp
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@
	touch $@
#MB 190301 - START
aflow_chull_jupyter.cpp: aflow_chull_jupyter.json
	rm -rf $@
	echo "// $@ automatic generated" >> $@
	echo "std::string AFLOW_CHULL_JUPYTER_JSON=\"\\" >> $@
	cat $^ | perl -p -e 's/\r//g' | perl -p -e 's/\\n/\\\\n/g' | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@	
	echo "\";" >> $@
#MB 190301 - END
#MB 190301 - START
aflow_chull_jupyter_plotter.cpp: aflow_chull_jupyter_plotter.py
	rm -rf $@
	echo "// $@ automatic generated" >> $@
	echo "std::string AFLOW_CHULL_JUPYTER_PLOTTER_PY=\"\\" >> $@
	cat $^ | perl -p -e 's/\r//g' | perl -p -e 's/\\d/\\\\d/g' | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
#MB 190301 - END
#MB 190301 - START
aflow_chull_python.cpp: aflow_chull_python.py
	rm -rf $@
	echo "// $@ automatic generated" >> $@
	echo "std::string AFLOW_CHULL_PYTHON_PY=\"\\" >> $@
	cat $^ | perl -p -e 's/\r//g' | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
#MB 190301 - END
#MB 190301 - START
aflow_chull_jupyter_requirements.cpp: aflow_chull_jupyter_requirements.txt
	rm -rf $@
	echo "// $@ automatic generated" >> $@
	echo "std::string AFLOW_CHULL_JUPYTER_REQUIREMENTS_TXT=\"\\" >> $@
	cat $^ | perl -p -e 's/\r//g' | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
#MB 190301 - END
aflowlib_web_interface.o: aflowlib_web_interface.cpp $(AFLOWLIB_WEBAPP_CPP)
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@
aflow_mix_nomix_table.cpp: Makefile $(AFLOW_HPP) $(AUROSTD_HPP)
	rm -rf nomix.tmp
	find $(PROJECT_DIRECTORY_GNDSTATE)/NOMIX/ -name aflow.in | sed "s/\/common\/$(PROJECT_DIRECTORY_GNDSTATE)\/NOMIX\/LIB2\/LIB\///g" | sed "s/aflow.in//g" >> nomix.tmp
	rm -rf $@
	echo "// $@ automatic generated " >> $@
	echo "#include <string>" >> $@
	echo "std::string NomixTable=\"\\" >> $@
	cat nomix.tmp | sort | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	echo "\";" >> $@
	chmod 755 $@
	rm -rf nomix.tmp
$(MIX_EXPERIMENTS_OBJ): $(MIX_EXPERIMENTS_CPP) $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $(MIX_EXPERIMENTS_CPP) -c -o $@
aflow_xproto.o: aflow_xproto.cpp aflow_xproto_library_default.cpp Makefile
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@
aflow_library_icsd.dat.xz: $(ICSD_LIB)/README_LIBRARY_ICSD1.TXT $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT $(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	$(MAKE) aflow -j 16
	rm -rf aflow_library_icsd.dat
	echo "[AFLOW] aflow_library_icsd.dat  " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD1.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD1.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD1.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD2.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD2.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD3.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD3.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD4.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD4.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD5.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD5.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD6.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD6.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD7.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD7.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD8.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD8.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD9.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD9.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD9.TXT]STOP " >> aflow_library_icsd.dat
	cat aflow_library_icsd.dat |sed "s/WYCK_ICSD/WYC_ICSD/g" |sed "s/  / /g"|sed "s/  / /g"|sed "s/  / /g"|sed "s/  / /g"|sed "s/  / /g" | grep -v "//" > aflow_library_icsd.dat2
	mv aflow_library_icsd.dat2 aflow_library_icsd.dat
	xz -9vf aflow_library_icsd.dat
	cp -f $@ ../AFLOW3_NONFREE/LIBS
	cp -f $@ /common/AFLOW/LIBS
	rm -f aflow_data_*
	$(MAKE) aflow_data
#	rm -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt*
	./aflow --protos_icsd > $(EXTRA_FREE)/NIST/ICSD_List_2014.txt
	rm -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz
	xz -9vf $(EXTRA_FREE)/NIST/ICSD_List_2014.txt
$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT:
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD1.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT:
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD2.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT:
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD3.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT:
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD4.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT:
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD5.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT:
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD6.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT:
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD7.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT:
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD8.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT:
	rm -rf $@
	echo "// ***************************************************************************">>$@
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD9.TXT.xz >>$@
	echo "// ***************************************************************************">>$@
	echo " STRUCTURE ___PROTO_END___">>$@
	echo " END">>$@
	echo "// ***************************************************************************">>$@
############################################################################################################
aflow_data_calculated.cpp: Makefile
	rm -rf $@
	if [ -f $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz ]; then echo "// $@ automatic generated from CALCULATED/vAUID.dat" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz ]; then echo "#include <string>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz ]; then echo "std::string vAUID=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz ]; then xzcat $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz | grep -v "//" | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz ]; then echo "// $@ automatic generated from CALCULATED/vAURL.dat" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz ]; then echo "#include <string>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz ]; then echo "std::string vAURL=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz ]; then xzcat $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz | grep -v "//" | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz ]; then echo "// $@ automatic generated from CALCULATED/vLOOP.dat" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz ]; then echo "#include <string>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz ]; then echo "std::string vLOOP=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz ]; then xzcat $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz | grep -v "//" | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz ]; then echo "\";" >> $@ ; fi
	if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	touch $@
aflow_data_htqc.cpp: README_LIBRARY_HTQC.TXT Makefile
	rm -f $@
	echo "// $@ automatic generated from $<" >> $@
	echo "#include <string>" >> $@
	echo "std::string Library_HTQC=\"\\" >> $@
	cat $< | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	cat README_LIBRARY_HTQC_BORIDES.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	cat README_LIBRARY_HTQC_CARBIDES.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	echo "STRUCTURE ___PROTO_END___"  | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	echo " END"  | perl -p -e 's/\n/ \\n\\\n/g' >> $@
	echo "\";" >> $@
#aflow_data_binary.cpp:  Makefile
#	rm -rf aflow_data_binary.cpp
#	echo "// aflow_data_binary.cpp automatic generated from PROTO/Binary.exp" >> aflow_data_binary.cpp
#	echo "std::string AFLOW_BinaryRead=\"\\" >> aflow_data_binary.cpp
#	cat $(EXTRA_FREE)/PROTO/Binary.exp | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_binary.cpp
#	echo "\";" >> aflow_data_binary.cpp
#	echo "std::string AFLOW_Binary_Angle_Read=\"\\" >> aflow_data_binary.cpp
#	cat $(EXTRA_FREE)/PROTO/Binary_angle.exp | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_binary.cpp
#	echo "\";" >> aflow_data_binary.cpp
aflow_data_libraries.cpp: Makefile
	rm -rf $@
	echo "// $@ automatic generated from ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_libX.dat" >> $@
	echo "#include <string>" >> $@
	@echo "LIB0"
	if [ -f ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib0.dat.xz ]; then echo "std::string aflowlib_lib0=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib0.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; else echo "std::string aflowlib_lib0=\"\";" >> $@; fi
	@echo "LIB1"
	if [ -f ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib1.dat.xz ]; then echo "std::string aflowlib_lib1=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib1.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; else echo "std::string aflowlib_lib1=\"\";" >> $@; fi
	@echo "LIB2"
	if [ -f ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib2.dat.xz ]; then echo "std::string aflowlib_lib2=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib2.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; else echo "std::string aflowlib_lib2=\"\";" >> $@; fi
	@echo "LIB3"
	if [ -f ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib3.dat.xz ]; then echo "std::string aflowlib_lib3=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib3.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; else echo "std::string aflowlib_lib3=\"\";" >> $@; fi
	@echo "LIB4"
	if [ -f ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib4.dat.xz ]; then echo "std::string aflowlib_lib4=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib4.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; else echo "std::string aflowlib_lib4=\"\";" >> $@; fi
	@echo "LIB5"
	if [ -f ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib5.dat.xz ]; then echo "std::string aflowlib_lib5=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib5.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; else echo "std::string aflowlib_lib5=\"\";" >> $@; fi
	@echo "LIB6"
	if [ -f ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib6.dat.xz ]; then echo "std::string aflowlib_lib6=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib6.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; else echo "std::string aflowlib_lib6=\"\";" >> $@; fi
	@echo "LIB7"
	if [ -f ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib7.dat.xz ]; then echo "std::string aflowlib_lib7=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib7.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; else echo "std::string aflowlib_lib7=\"\";" >> $@; fi
	@echo "LIB8"
	if [ -f ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib8.dat.xz ]; then echo "std::string aflowlib_lib8=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib8.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; else echo "std::string aflowlib_lib8=\"\";" >> $@; fi
	@echo "LIB9"
	if [ -f ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib9.dat.xz ]; then echo "std::string aflowlib_lib9=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib9.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; else echo "std::string aflowlib_lib9=\"\";" >> $@; fi
	@echo "ICSD"
	if [ -f ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_icsd.dat.xz ]; then echo "std::string aflowlib_icsd=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_icsd.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; else echo "std::string aflowlib_icsd=\"\";" >> $@; fi
#[CO191217 - OBSOLETE]#	echo "std::string aflowlib_lib0=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib0.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@	
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib0=\"\";" >> $@
#[CO191217 - OBSOLETE]	echo "LIB1"
#[CO191217 - OBSOLETE]#	echo "std::string aflowlib_lib1=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib1.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@	
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib1=\"\";" >> $@
#[CO191217 - OBSOLETE]	echo "LIB2"
#[CO191217 - OBSOLETE]#	echo "std::string aflowlib_lib2=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib2.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@	
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib2=\"\";" >> $@
#[CO191217 - OBSOLETE]	echo "LIB3"
#[CO191217 - OBSOLETE]#	echo "std::string aflowlib_lib3=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib3.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib3=\"\";" >> $@
#[CO191217 - OBSOLETE]	echo "LIB4"
#[CO191217 - OBSOLETE]#	echo "std::string aflowlib_lib4=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib4.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib4=\"\";" >> $@
#[CO191217 - OBSOLETE]	echo "LIB5"
#[CO191217 - OBSOLETE]#	echo "std::string aflowlib_lib5=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib5.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib5=\"\";" >> $@
#[CO191217 - OBSOLETE]	echo "LIB6"
#[CO191217 - OBSOLETE]#	echo "std::string aflowlib_lib6=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib6.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib6=\"\";" >> $@
#[CO191217 - OBSOLETE]	echo "LIB7"
#[CO191217 - OBSOLETE]#	echo "std::string aflowlib_lib7=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib7.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib7=\"\";" >> $@
#[CO191217 - OBSOLETE]	echo "LIB8"
#[CO191217 - OBSOLETE]#	echo "std::string aflowlib_lib8=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib8.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib8=\"\";" >> $@
#[CO191217 - OBSOLETE]	echo "LIB9"
#[CO191217 - OBSOLETE]#	echo "std::string aflowlib_lib9=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib9.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib9=\"\";" >> $@
#[CO191217 - OBSOLETE]	echo "ICSD"
#[CO191217 - OBSOLETE]#	echo "std::string aflowlib_icsd=\"\\" >> $@ && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_icsd.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> $@ && echo "\";" >> $@
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_icsd=\"\";" >> $@
	touch $@
aflow_data_stokes.cpp: Makefile 
	rm -rf $@
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "// $@ automatic generated " >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "#include <string>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "std::string FINDSYM_data_space_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then echo "std::string FINDSYM_data_wyckoff_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then echo "std::string FROZSL_data_space_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then echo "std::string FROZSL_data_wyckoff_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then echo "std::string FROZSL_data_images_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g'  | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then echo "std::string FROZSL_data_irreps_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then echo "std::string FROZSL_data_isotropy_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then echo "std::string FROZSL_data_little_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then echo "std::string FROZSL_symmetry2_dat=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then echo "std::string FROZSL_const_dat=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then echo "std::string FROZSL_phvaspsetup_AFLOW=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then xzcat $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz | perl -p -e 's/\#/_ASCII_23_/g' | perl -p -e 's/\"/_ASCII_22_/g' | perl -p -e 's/\\/_ASCII_5C_/g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then echo "std::string FROZSL_phvaspsetup_POSCAR=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then xzcat $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz | perl -p -e 's/\#/_ASCII_23_/g' | perl -p -e 's/\"/_ASCII_22_/g' | perl -p -e 's/\\/_ASCII_5C_/g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then echo "\";" >> $@ ; fi
	if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	touch $@
aflow_data_nist.cpp: Makefile 
	rm -rf $@
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "// $@ automatic generated from NIST" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "#include <string>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "std::string ElectronStoppingPower_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then echo "std::string PhotonCrossSection_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then echo "std::string PhotonStoppingPower_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then echo "std::string ICSD_List_txt=\"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then echo "\";" >> $@ ; fi
	if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	touch $@
aflow_data_aflow_potcars.cpp: Makefile 
	rm -rf $@
	echo "// $@ automatic generated from AFLOW_PSEUDOPOTENTIALS.TXT" >> $@
	echo "#include <string>" >> $@
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then echo "bool AFLOW_PSEUDOPOTENTIALS=1;" >> $@ ; else echo "bool AFLOW_PSEUDOPOTENTIALS=0;" >> $@ ; fi
	echo "std::string AFLOW_PSEUDOPOTENTIALS_TXT=\"\\" >> $@
	#if[ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT_XZ_BASE64.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT_XZ_BASE64.xz | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	echo "\";" >> $@
	echo "std::string AFLOW_PSEUDOPOTENTIALS_LIST_TXT=\"\\" >> $@
	#if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz | grep "POTCAR=" | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS_LIST.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS_LIST.TXT.xz | grep "POTCAR=" | perl -p -e 's/\n/\\n\\\n/g' >> $@ ; fi
	echo "\";" >> $@
aflow_data_latex.cpp: Makefile
	rm -rf $@
	echo "// $@ automatic generated from LATEX" >> $@
	echo "#include <string>" >> $@
	echo "std::string f144468a7ccc2d3a72ba44000715efdb=\"\\" >> $@
	if [ -f ../AFLOW3_NONFREE/EXTRA/LATEX/f144468a7ccc2d3a72ba44000715efdb ]; then cat ../AFLOW3_NONFREE/EXTRA/LATEX/f144468a7ccc2d3a72ba44000715efdb* | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g'  | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	echo "\";" >> $@
	echo "std::string d0f1b0e47f178ae627a388d3bf65d2d2=\"\\" >> $@
	if [ -f ../AFLOW3_NONFREE/EXTRA/LATEX/d0f1b0e47f178ae627a388d3bf65d2d2 ]; then cat ../AFLOW3_NONFREE/EXTRA/LATEX/d0f1b0e47f178ae627a388d3bf65d2d2 | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g'  | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	echo "\";" >> $@
	echo "std::string decf00ca3ad2fe494eea8e543e929068=\"\\" >> $@
	if [ -f ../AFLOW3_NONFREE/EXTRA/LATEX/decf00ca3ad2fe494eea8e543e929068 ]; then cat ../AFLOW3_NONFREE/EXTRA/LATEX/decf00ca3ad2fe494eea8e543e929068 | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g'  | perl -p -e 's/\n/ \\n\\\n/g' >> $@ ; fi
	echo "\";" >> $@
aflow_data_readme.cpp: README_AFLOW_LICENSE_GPL3.TXT README_AFLOW_ACONVASP.TXT README_AFLOW_APL.TXT README_AFLOW_QHA_SCQHA_QHA3P.TXT README_AFLOW_AGL.TXT README_AFLOW_AEL.TXT README_AFLOW_POCC.TXT README_AFLOW_ANRL.TXT README_AFLOW_COMPARE.TXT README_AFLOW_GFA.TXT README_AFLOW.TXT README_PROTO.TXT README_AFLOW_APENNSY.TXT README_AFLOW_FROZSL.TXT README_AFLOW_SCRIPTING.TXT README_AFLOW_SYM.TXT README_AFLOW_CCE.TXT README_AFLOW_CHULL.TXT README_AFLOW_XAFLOW.TXT README_AFLOW_AFLOWRC.TXT README_AFLOW_EXCEPTIONS.TXT README_CONTRIBS.TXT README_AFLOW_VERSIONS_HISTORY.TXT Makefile
	rm -rf $@
	echo "// $@ automatic generated" >> $@
	echo "#include <string>" >> $@
	echo "std::string README_AFLOW_LICENSE_GPL3_TXT=\"\\" >> $@
	cat README_AFLOW_LICENSE_GPL3.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_TXT=\"\\" >> $@
	cat README_AFLOW.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_VERSIONS_HISTORY_TXT=\"\\" >> $@
	cat README_AFLOW_VERSIONS_HISTORY.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_PFLOW_TXT=\"\\" >> $@
	cat README_AFLOW_ACONVASP.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_APENNSY_TXT=\"\\" >> $@
	cat README_AFLOW_APENNSY.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_SCRIPTING_TXT=\"\\" >> $@
	cat README_AFLOW_SCRIPTING.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_FROZSL_TXT=\"\\" >> $@
	cat README_AFLOW_FROZSL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_POCC_TXT=\"\\" >> $@
	cat README_AFLOW_POCC.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_APL_TXT=\"\\" >> $@
	cat README_AFLOW_APL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_QHA_SCQHA_QHA3P_TXT=\"\\" >> $@
	cat README_AFLOW_QHA_SCQHA_QHA3P.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_AGL_TXT=\"\\" >> $@
	cat README_AFLOW_AGL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_AEL_TXT=\"\\" >> $@
	cat README_AFLOW_AEL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_ANRL_TXT=\"\\" >> $@
	cat README_AFLOW_ANRL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_COMPARE_TXT=\"\\" >> $@
	cat README_AFLOW_COMPARE.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_GFA_TXT=\"\\" >> $@
	cat README_AFLOW_GFA.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_SYM_TXT=\"\\" >> $@
	cat README_AFLOW_SYM.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_CCE_TXT=\"\\" >> $@
	cat README_AFLOW_CCE.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_CHULL_TXT=\"\\" >> $@
	cat README_AFLOW_CHULL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_EXCEPTIONS_TXT=\"\\" >> $@
	cat README_AFLOW_EXCEPTIONS.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_PROTO_TXT=\"\\" >> $@
	cat README_PROTO.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_XAFLOW_TXT=\"\\" >> $@
	cat README_AFLOW_XAFLOW.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
	echo "std::string README_AFLOW_AFLOWRC_TXT=\"\\" >> $@
	cat README_AFLOW_AFLOWRC.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
aflow_data_extra.cpp: Makefile aflow_library_icsd.dat.xz
	rm -rf $@
	echo "// $@ automatic generated from $(NIST_LIB)/X-ary.icsd" >> $@
	echo "#include <string>" >> aflow_data_libraries.cpp
#CO191217 - the ls command will fail if $(NIST_LIB)/1-ary.icsd.xz is not present: mimicking previous behavior
	@ls $(NIST_LIB)/1-ary.icsd.xz > /dev/null
	if [ -f $(NIST_LIB)/1-ary.icsd.xz ]; then echo "std::string icsd_1_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/1-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/2-ary.icsd.xz ]; then echo "std::string icsd_2_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/2-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/3-ary.icsd.xz ]; then echo "std::string icsd_3_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/3-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/4-ary.icsd.xz ]; then echo "std::string icsd_4_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/4-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/5-ary.icsd.xz ]; then echo "std::string icsd_5_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/5-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/6-ary.icsd.xz ]; then echo "std::string icsd_6_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/6-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/7-ary.icsd.xz ]; then echo "std::string icsd_7_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/7-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/8-ary.icsd.xz ]; then echo "std::string icsd_8_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/8-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/9-ary.icsd.xz ]; then echo "std::string icsd_9_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/9-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/10-ary.icsd.xz ]; then echo "std::string icsd_10_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/10-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/11-ary.icsd.xz ]; then echo "std::string icsd_11_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/11-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/12-ary.icsd.xz ]; then echo "std::string icsd_12_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/12-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/13-ary.icsd.xz ]; then echo "std::string icsd_13_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/13-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/14-ary.icsd.xz ]; then echo "std::string icsd_14_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/14-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
	if [ -f $(NIST_LIB)/15-ary.icsd.xz ]; then echo "std::string icsd_15_ary=\"\\" >> $@ && xzcat $(NIST_LIB)/15-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@; fi
#	echo "std::string Library_ICSD=\"\\" >> $@ && xzcat  ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> $@ && echo "\";" >> $@
#[CO191217 - OBSOLETE]	$(MAKE) aflow_data_libraries.cpp
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib0=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib0.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib1=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib1.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib2=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib2.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib3=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib3.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib4=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib4.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib5=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib5.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib6=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib6.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib7=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib7.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib8=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib8.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_lib9=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib9.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
#[CO191217 - OBSOLETE]	echo "std::string aflowlib_icsd=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_icsd.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	touch $@
#aflow_data: aflow_data_readme.cpp aflow_data_calculated.cpp aflow_data_htqc.cpp aflow_data_libraries.cpp aflow_data_stokes.cpp aflow_data_nist.cpp aflow_data_latex.cpp aflow_data_aflow_potcars.cpp aflow_data.cpp
aflow_data: aflow_data_readme.o aflow_data_calculated.o aflow_data_htqc.o aflow_data_libraries.o aflow_data_stokes.o aflow_data_nist.o aflow_data_latex.o aflow_data_aflow_potcars.o aflow_data.cpp
	rm -f aflow_data_extra.cpp
	echo "" >> aflow_data_extra.cpp
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "#define defined_Library_ICSD" >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "std::string Library_ICSD=\"\\" >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then xzcat  ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "\";" >> aflow_data_extra.cpp; fi
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -c aflow_data.cpp -o aflow_data.o
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -o aflow_data aflow_data.o aflow_data_*.o 
	echo aflow_data aflow_data.o aflow_data_*.o
.PHONY: amir
amir:	aflow_data_readme.cpp aflow_data_calculated.cpp aflow_data_htqc.cpp aflow_data_libraries.cpp aflow_data_stokes.cpp aflow_data_nist.cpp aflow_data_latex.cpp aflow_data_aflow_potcars.cpp aflow_data.cpp aflow_data_extra.cpp aflow_data_libraries.cpp
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -DICSD -o aflow_data aflow_data.cpp
############################################################################################################
aflow_phonons.o: aflow_phonons.cpp $(AFLOW_HPP) $(AUROSTD_HPP)  Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@
aflow_matlab_funcs.cpp: Makefile
	rm -rf $@
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "// $@ automatic generated from MATLAB/*" >>  $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "#include <sstream>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "// $@ automatic generated from MATLAB/param.m" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "std::string MATLAB_FUNCS_param(void){" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "   std::stringstream strstream; strstream << \"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then cat $(EXTRA_FREE)/MATLAB/param.m | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "\" << std::endl; return strstream.str();};" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "// $@ automatic generated from MATLAB/plotband.m" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "// std::string MATLAB_FUNCS_plotband(std::string DIRECTORY,std::string OPTION1){" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "std::string MATLAB_FUNCS_plotband(void){" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "   std::stringstream strstream; strstream << \"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then cat $(EXTRA_FREE)/MATLAB/plotband.m | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' | perl -p -e 's/MATLAB_FUNCS_DIRECTORY/\" << DIRECTORY << \"/g' | perl -p -e 's/DOS_SCALE_MODE/\" << OPTION1 << \"/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "\" << std::endl; return strstream.str();};" >> $@ ; fi
	if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	touch $@
aflow_gnuplot_funcs.cpp: Makefile
	rm -rf $@
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "// $@ automatic generated from GNUPLOT/*" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "#include <sstream>" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "// $@ automatic generated from GNUPLOT/plotbz.sh" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "//std::string GNUPLOT_FUNCS_plotbz(std::string DIRECTORY,std::string OPTION1){" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "std::string GNUPLOT_FUNCS_plotbz(void){" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "   std::stringstream strstream; strstream << \"\\" >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then cat $(EXTRA_FREE)/GNUPLOT/plotbz.sh | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' | perl -p -e 's/GNUPLOT_FUNCS_DIRECTORY/\" << DIRECTORY << \"/g' | perl -p -e 's/GNUPLOT_FUNCS_OPTION1/\" << OPTION1 << \"/g' >> $@ ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "\" << std::endl; return strstream.str();};" >> $@ ; fi
	if [ ! -f $@ ]; then rm -rf $@.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/$@.xz && xz -dv $@.xz ; fi
	touch $@
aurostd_xtemplates.o: aurostd_xtemplates.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@
#
###### HTRESOURCES #################################################
aflow_contrib_stefano_htcurriculum.o: $(AFLOWLIB_HPP) $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -c -o $@
aflow_contrib_stefano_htsecurity.o: aflow_contrib_stefano_htsecurity.cpp $(AFLOWLIB_HPP) $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@
aflowlib_web_outreach.o: aflowlib_web_outreach.cpp $(AFLOWLIB_HPP) $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@

###### SHIDONG #################################################
aflow_contrib_shidong_main.o: aflow_contrib_shidong_main.cpp aflow_contrib_shidong_main.h aflow_contrib_shidong.h aflow_contrib_shidong_cluster_expansion.h aflow_contrib_shidong_auxiliary.h 
aflow_contrib_shidong.o: aflow_contrib_shidong.cpp aflow_contrib_shidong.h aflow_contrib_shidong_auxiliary.h aflow_contrib_shidong_cluster_expansion.h 
aflow_contrib_shidong_auxiliary.o: aflow_contrib_shidong_auxiliary.cpp aflow_contrib_shidong_auxiliary.h aflow_contrib_shidong_cluster_expansion.h
aflow_contrib_shidong_cexstructure.o: aflow_contrib_shidong_cluster_expansion.cpp aflow_contrib_shidong_cluster_expansion.h
aflow_contrib_shidong_mc.o: aflow_contrib_shidong_mc.cpp aflow_contrib_shidong_mc.h aflow_contrib_shidong_statistic.h aflow_contrib_shidong_cluster_expansion.h 
aflow_contrib_shidong_statistic.o: aflow_contrib_shidong_statistic.cpp aflow_contrib_shidong_statistic.h

#####COMPARE###################################################
aflow_compare_structure_function.o: aflow_compare_structure_function.cpp Makefile $(AFLOW_HPP) $(AUROSTD_HPP)
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" ${INCLUDE} ${CFLAGS} -std=gnu++0x $< -c -o $@
aflow_compare_structure.o: aflow_compare_structure.cpp Makefile $(AFLOW_HPP) $(AUROSTD_HPP)
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" ${INCLUDE} ${CFLAGS} -std=gnu++0x $< -c -o $@

#[OBSOLETE] [JUNKAI]###### JUNKAI #################################################
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_cproto.o:aflow_contrib_junkai_cproto.cpp aflow_contrib_junkai_basic.h
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_primitivetrans.o:aflow_contrib_junkai_primitivetrans.cpp aflow_contrib_junkai_basic.h
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_basic.o:aflow_contrib_junkai_basic.cpp aflow_contrib_junkai_basic.h
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_protocheck.o:aflow_contrib_junkai_protocheck.cpp aflow_contrib_junkai_basic.h
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_magnet.o:aflow_contrib_junkai_magnet.cpp aflow_contrib_junkai_basic.h
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_phasediag.o:aflow_contrib_junkai_phasediag.cpp aflow_contrib_junkai_basic.h
###### ACONVASP #################################################
aflow_aconvasp_main.o: aflow_aconvasp_main.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@
aflow_apennsy.o: aflow_apennsy.cpp $(AFLOW_HPP) $(AUROSTD_HPP) aflow_apennsy.h Makefile
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@
aflow_apennsy_main.o: aflow_apennsy_main.cpp $(AFLOW_HPP) $(AUROSTD_HPP) aflow_apennsy.h Makefile
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $< -c -o $@
aflow_wyckoff.o: aflow_wyckoff.cpp $(AFLOW_HPP) $(AUROSTD_HPP)
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" -D_WICKOFF_PRE_ $< -c -o $@
####### APL STUFF ####################################################
#DX and CO - START
APL/aflow_apl_doscalc.o: APL/aflow_apl_doscalc.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x $< -c -o $@
	touch $@
APL/aflow_apl_pdisc.o: APL/aflow_apl_pdisc.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x $< -c -o $@
	touch $@
APL/aflow_aapl_tcond.o: APL/aflow_aapl_tcond.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x $< -c -o $@
	touch $@
APL/aflow_qha_dm_pdos_save.o: APL/aflow_qha_dm_pdos_save.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x $< -c -o $@
	touch $@
APL/aflow_apl_group_velocity.o: APL/aflow_apl_group_velocity.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x $< -c -o $@
	touch $@
APL/aflow_apl_atomic_disp.o: APL/aflow_apl_atomic_disp.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x $< -c -o $@
	touch $@
APL/aflow_qha_gruneisen.o: APL/aflow_qha_gruneisen.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x $< -c -o $@
	touch $@
APL/aflow_scqha_gruneisen.o: APL/aflow_scqha_gruneisen.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x $< -c -o $@
	touch $@
APL/aflow_scqha_eos.o: APL/aflow_scqha_eos.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x $< -c -o $@
	touch $@
APL/aflow_qha3phonons_eos.o: APL/aflow_qha3phonons_eos.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x $< -c -o $@
	touch $@
APL/aflow_scqha_T_freqs.o: APL/aflow_scqha_T_freqs.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x $< -c -o $@
	touch $@
#DX - CO - END
#APL/aflow_apl_gsa.o: APL/aflow_apl_gsa.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
##	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
#	${CPP} $(VERS) -D_AFLOW_FILE_NAME_=\""$<"\" $(INCLUDE) $(CCFLAGS_GSA) $< -c -o $@
#	touch $@

################################################################################################
.PHONY: install install_auro install_root 
install: all
	cp -f aflow $(BIN)/aflow.$(VERSNUMBER)
	cp -f aflow_data $(BIN)/
	cp -f $(BIN)/aflow.$(VERSNUMBER) $(BIN)/aflow
	ls -las $(BIN)/aflow*
	$(MAKE) check
install_auro: install_user
	./aflow --version | grep AFLOW >> $(BIN)/etc_aflow_status.txt
	chmod 644 $(BIN)/etc_aflow_status.txt
	./aflow --calculated | grep Calcs >> $(BIN)/etc_aflowlib_status.txt
	chmod 644 $(BIN)/etc_aflowlib_status.txt
	./aflow --apool > $(BIN)/etc_aflowlib_apool.txt
	chmod 644 $(BIN)/etc_aflowlib_apool.txt
install_root: all
	cp -f aflow /usr/local/bin/aflow.$(VERSNUMBER)
	cp -f aflow_data /usr/local/bin/aflow_data
	cp -f /usr/local/bin/aflow.$(VERSNUMBER) /usr/local/bin/aflow
	ls -las /usr/local/bin/aflow*
	$(MAKE) check

.PHONY: send_version_updates_email
#CO190131 - START
send_version_updates_email:
	rm -f version_updates_email.txt
	$(MAKE) version_updates_email.txt
	[ -s version_updates_email.txt ]
	cat version_updates_email.txt | mail -r aflow@materials.duke.edu -s "AFLOW Code Update: v$(VERSNUMBER)" aflow+int+737+2740356048620877160@groups.io
	rm -f version_updates_email.txt
#CO190131 - STOP
#CO190218 - START
version_updates_email.txt: README_AFLOW_VERSIONS_HISTORY.TXT
	rm -f $(@:.txt=_tmp.txt)
	@file_content=""; \
	found_content=false; \
	while IFS= read -r line; do \
		if [ "$$found_content" = true ]; then \
			if [ -n "$$(echo $$line | grep "Location" )" ] && [ -z "$$(echo $$line | grep "$(VERSNUMBER)" )" ]; then \
				found_content=false; \
				break; \
			fi; \
			file_content="$${file_content}\n$${line}"; \
		elif [ -n "$$( echo "$$line" | grep "$(VERSNUMBER)" )" ]; then \
			found_content=true; \
			file_content="$$line"; \
		fi; \
	done < $<; \
	line_count=$$(echo $$file_content | wc -l); \
	IFS=; \
	file_content=$$(echo $$file_content | head -n $$(( $$line_count - 1 ))); \
	if [ -z "$$file_content" ]; then echo "ERROR - No relevant version content found. Exiting."; exit 1; fi; \
	echo $$file_content >> $(@:.txt=_tmp.txt);
	rm -f $@;
	echo "Dear AFLOW user," >> $@
	echo "" >> $@
	echo "A new version of the AFLOW C++ code is now available (v$(VERSNUMBER))." >> $@
	echo "Update details:" >> $@
	echo "" >> $@
	cat $(@:.txt=_tmp.txt) >> $@; rm -f $(@:.txt=_tmp.txt)
	echo "" >> $@
	echo "Thanks," >> $@
	echo "AFLOW Dev" >> $@
#CO190218 - STOP

.PHONY: niet
niet:	
	dropbox stop
	$(MAKE) clean
	$(MAKE) icsdclean
	$(MAKE) extra_clean
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	$(MAKE) web
	$(MAKE) backup
	$(MAKE) icsdclean
	$(MAKE) -j 32
	rm -fv aflow_data_calculated.cpp* && $(MAKE) aflow_data_calculated.cpp && cp -f aflow_data_calculated.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv /www/AFLOW/AFLOW3_FREE/aflow_data_calculated.cpp
	rm -fv aflow_data_stokes.cpp* && $(MAKE) aflow_data_stokes.cpp && cp -f aflow_data_stokes.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv /www/AFLOW/AFLOW3_FREE/aflow_data_stokes.cpp
	rm -fv aflow_data_nist.cpp* && $(MAKE) aflow_data_nist.cpp && cp -f aflow_data_nist.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv /www/AFLOW/AFLOW3_FREE/aflow_data_nist.cpp
	rm -fv aflow_matlab_funcs.cpp* && $(MAKE) aflow_matlab_funcs.cpp && cp -f aflow_matlab_funcs.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv /www/AFLOW/AFLOW3_FREE/aflow_matlab_funcs.cpp
	rm -fv aflow_gnuplot_funcs.cpp* && $(MAKE) aflow_gnuplot_funcs.cpp && cp -f aflow_gnuplot_funcs.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv /www/AFLOW/AFLOW3_FREE/aflow_gnuplot_funcs.cpp
	dropbox start

.PHONY: backup
backup: realclean clean icsdclean 
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean ; fi
	mkdir -p aflow.$(VERSNUMBER)
#	mkdir -p aflow.$(VERSNUMBER)/SOURCES
	cp -drvf a*cpp a*.h *.js *.json *.py *.txt README* Makefile APL ANRL* AUROSTD SQLITE aflow.$(VERSNUMBER)
# no EXTRA in SRC	if [ -d $(EXTRA_FREE) ]; then cp -drvf $(EXTRA_FREE) aflow.$(VERSNUMBER) ; fi
	if [ -d SOURCE ]; then cp -drvf SOURCE aflow.$(VERSNUMBER) ; fi
	tar Jcfvp aflow.$(VERSNUMBER).tar.xz aflow.$(VERSNUMBER)
	rm -rfv $(BACKUP)/aflow.$(VERSNUMBER) $(BACKUP)/aflow.$(VERSNUMBER).tar.xz
	if [ -d /www/AFLOW ]; then cp -f aflow.$(VERSNUMBER).tar.xz /www/AFLOW ; fi
	mv aflow.$(VERSNUMBER).tar.xz $(BACKUP)
	rm -rfv aflow.$(VERSNUMBER)

.PHONY: www web
www:	web
web:	realclean
	rm -rvf aflow3.tar.xz aflow_test.tar.xz aflow3 ~/aflow3.tar.xz ~/aflow_test.tar.xz ~/aflow3
	tar Jcfvp ~/aflow3.tar.xz Makefile *cpp *.h *.js *.json *.py *.txt READ* APL ANRL* AUROSTD SQLITE
# no EXTRA IN SRC	if [ -d $(EXTRA_FREE) ]; then tar Jcfvp ~/aflow3.tar.xz Makefile *cpp *.h *.js *.json *.py *.txt READ* APL ANRL* AUROSTD SQLITE $(EXTRA_FREE) ; fi
# no other SOURCES in SRC	if [ -d SOURCES ]; then tar Jcfvp ~/aflow3.tar.xz Makefile *cpp *.h *.js *.json *.py *.txt READ* APL ANRL* AUROSTD SQLITE $(EXTRA_FREE) SOURCES ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_AEL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_AGL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_ANRL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_APENNSY.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_APL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_QHA_SCQHA_QHA3P.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_ACONVASP.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_COMPARE.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_GFA.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_FROZSL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_LICENSE_GPL3.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_POCC.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_SCRIPTING.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_SYM.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_CCE.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_CHULL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_EXCEPTIONS.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_VERSIONS_HISTORY.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_XAFLOW.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_AFLOWRC.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_CONTRIBS.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_PROTO.TXT /www/AFLOW ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl ]; then cp ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl /www/AFLOW ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl_init ]; then cp ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl_init /www/AFLOW ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym ]; then cp ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym /www/AFLOW ; fi
	if [ -d /www/auro ]; then mv -f  ~/aflow3.tar.xz /www/auro ; fi

.PHONY: wget distroget patch_makefile1 patch_makefile2
wget:	realclean
	rm -rf aflow3.tar.xz* aflow_test.tar.xz*
	wget http://materials.duke.edu/auro/aflow3.tar.xz
	tar xfvp aflow3.tar.xz
	rm -rf aflow3.tar.xz 
distroget: wget
	rm -rf TESTS PLATON* FINDSYM* FROZSL*
	cat Makefile | sed "s/Makefile /Makefile /g" | sed "s/\*\.a//g" > M
	mv M Makefile
patch_makefile1: Makefile
	cat Makefile | sed "s/ \-Wno\-unused\-but\-set\-variable//g" > M	
	mv M Makefile
patch_makefile2: Makefile
	cat Makefile | sed "s/\-O1//g" | sed "s/\-O2//g" | sed "s/\-O3//g" > M	
	mv M Makefile

.PHONY: check proto_check
check:
	$(MAKE) proto_check

.PHONY: check_full gcd_check smith_check proto_check slab_check pocc_check aflowSG_check aflowSYM_check aflow_xtal_match_check
check_full:
	$(MAKE) gcd_check
	$(MAKE) smith_check
	$(MAKE) proto_check
	$(MAKE) slab_check
	$(MAKE) pocc_check
	$(MAKE) aflowSG_check
	$(MAKE) aflowSYM_check
	$(MAKE) aflow_xtal_match_check
proto_check:
	@echo "Performing 'aflow --proto' test"
	@if [ -z `which md5sum` ]; then echo "Cannot find binary md5sum - check skipped."; else \
		PROTO4SUM=$$(./aflow --proto=4:Bi:Tc | ./aflow --sprim | md5sum | cut -d ' ' -f1) && \
		[ "$${PROTO4SUM#*ffb9c5681045fb80f391e9a931f21d1}" != "$$PROTO4SUM" ] && \
		echo "'aflow --proto' check pass!" && exit 0 || \
		echo "'aflow --proto' test failed. Check compilation of aflow and aflow_data." && exit 1; \
	fi
#ME190313 - only execute check if md5sum is present
#	-if [ -z `which md5sum` ]; then echo "Cannot find binary md5sum - check skipped."; else ./aflow --proto=4:Bi:Tc | ./aflow --sprim | md5sum; echo "#ffb9c5681045fb80f391e9a931f21d1"; fi
#	./aflow --proto=4:Bi:Tc | ./aflow --sprim | md5sum
#	#ffb9c5681045fb80f391e9a931f21d1
#	./aflow --checki | grep FAIL
slab_check:
	./aflow --slab_test
gcd_check:
	./aflow --gcd_test
smith_check:
	./aflow --smith_test
pocc_check:
	@echo "Running POCC test 1" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:Ti --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 2" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 3" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 4" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 5" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 6" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 7" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 8" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 9" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 10" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 11" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 12" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 13" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 14" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 15" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 16" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 17" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 18" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 19" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 20" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 21" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 22" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 23" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 24" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 25" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 26" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 27" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 28" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 29" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 30" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 31" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 32" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 33" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 34" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 35" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 36" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 37" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 38" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 39" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 40" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 41" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 42" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 43" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 44" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 45" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 46" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 47" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 48" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 49" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 50" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 51" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 52" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 53" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 54" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 55" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 56" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Ta:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 57" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 58" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Nb:Ta:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 59" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Nb:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 60" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Hf:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 61" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Hf:Mo:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 62" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "POCC tests passed"
aflowSG_check:
	@echo "Running aflowSG test 1" && ANSWER="Pm-3m #221" && RESULT=$$(./aflow --proto=AB_cP2_221_b_a:Cl:Cs | ./aflow --aflowSG) && \
		if [ "$$RESULT" != "$$ANSWER" ]; then "aflowSG did not return $$ANSWER for './aflow --proto=AB_cP2_221_b_a:Cl:Cs | ./aflow --aflowSG', found $$RESULT instead" && exit 1; fi
	@echo "aflowSG test 1 passed!"
aflowSYM_check:
	@echo "Running aflowSYM test 1" && ANSWER_PGROUPS=8 && ANSWER_PGROUPKS=8 && ANSWER_FGROUPS=4 && ANSWER_PGROUP_XTALS=4 && ANSWER_PGROUPK_XTALS=4 && ANSWER_AGROUPS=24 && ANSWER_SGROUPS=108 && RESULT=$$(./aflow --proto=ABC2_mP8_10_ac_eh_mn-001:Au:Ag:Te | ./aflow --aflowSYM --screen_only) && COUNT_PGROUPS=$$(echo $$RESULT | awk -F'pgroup ' '{print NF-1}') && COUNT_PGROUPKS=$$(echo $$RESULT | awk -F'pgroupk ' '{print NF-1}') && COUNT_FGROUPS=$$(echo $$RESULT | awk -F'fgroup ' '{print NF-1}') && COUNT_PGROUP_XTALS=$$(echo $$RESULT | awk -F'pgroup_xtal ' '{print NF-1}') && COUNT_PGROUPK_XTALS=$$(echo $$RESULT | awk -F'pgroupk_xtal ' '{print NF-1}') && COUNT_AGROUPS=$$(echo $$RESULT | awk -F'agroup ' '{print NF-1}') && COUNT_SGROUPS=$$(echo $$RESULT | awk -F'sgroup ' '{print NF-1}') && \
		if [ "$$ANSWER_PGROUPS" -ne "$$COUNT_PGROUPS" ] || [ "$$ANSWER_PGROUPKS" -ne "$$COUNT_PGROUPKS" ] || [ "$$ANSWER_FGROUPS" -ne "$$COUNT_FGROUPS" ] || [ "$$ANSWER_PGROUP_XTALS" -ne "$$COUNT_PGROUP_XTALS" ] || [ "$$ANSWER_PGROUPK_XTALS" -ne "$$COUNT_PGROUPK_XTALS" ] || [ "$$ANSWER_AGROUPS" -ne "$$COUNT_AGROUPS" ] || [ "$$ANSWER_SGROUPS" -ne "$$COUNT_SGROUPS" ]; then "aflowSYM did not return one of the correct group cardinalities for './aflow --proto=ABC2_mP8_10_ac_eh_mn-001:Au:Ag:Te | ./aflow --aflowSYM --screen_only'; (ANSWER==FOUND) $$ANSWER_PGROUPS==$$COUNT_PGROUPS, $$ANSWER_PGROUPKS==$$COUNT_PGROUPKS, $$ANSWER_FGROUPS==$$COUNT_FGROUPS, $$ANSWER_PGROUP_XTALS==$$COUNT_PGROUP_XTALS, $$ANSWER_PGROUPK_XTALS==$$COUNT_PGROUPK_XTALS, $$ANSWER_AGROUPS==$$COUNT_AGROUPS, $$ANSWER_SGROUPS==$$COUNT_SGROUPS"  && exit 1; fi
	@echo "aflowSYM test 1 passed!"
aflow_xtal_match_check:
	@echo "Running aflow_xtal_match test 1" && ANSWER=3 && RESULT=$$(./aflow --proto=A2B3_hR10_167_c_e:Al:O --params=-1,2.72957758313,0.36,0.5561 | ./aflow --compare2protos --print=txt --quiet) && NUM_MATCHING_PROTOS=$$(echo $$RESULT | awk -F'duplicate_compounds=' '{print $$2}' | awk '{print $$1}') && \
		if [ "$$NUM_MATCHING_PROTOS" -ne "$$ANSWER" ]; then "aflow_xtal_match did not return $$ANSWER for './aflow --proto=A2B3_hR10_167_c_e:Al:O --params=-1,2.72957758313,0.36,0.5561 | ./aflow --compare2protos --print=txt --quiet', found $$NUM_MATCHING_PROTOS instead" && exit 1; fi
	@echo "aflow_xtal_match test 1 passed!"

.PHONY: clean realclean clean_autogen
clean:
	rm -rf *~ *.o *.o-* *.a *.rpo \#*
	rm -rf ./AUROSTD/*~ ./AUROSTD/*.o ./AUROSTD/*.o-* ./AUROSTD/*.a ./AUROSTD/*.rpo ./AUROSTD/\#*
	rm -rf ./APL/*~ ./APL/*.o ./APL/*.o-* ./APL/*.a ./APL/*.rpo ./APL/\#*
	rm -rf ./SQLITE/*~ ./SQLITE/*.o ./SQLITE/*.o-* ./SQLITE/*.a ./SQLITE/*.rpo ./SQLITE/\#*
	rm -rf ./EXTRA/GUS/*~ ./EXTRA/GUS/*.o ./EXTRA/GUS/*.o-* ./EXTRA/GUS/*.a ./EXTRA/GUS/*.rpo ./EXTRA/GUS/\#*
	rm -rf `find . -name *~`
	rm -rf aflow aflow_data aflowd aflow2 aflowd2 aflow1 aflowd1 aflow.exe aconvasp apennsy aflow_convasp* aflow_pennsy
	rm -rf *.eps *.pdf *ps *aux *dvi *log *m *tex 
	rm -rf nohup* core* VASP_TEST* GRND_TEST* xa xb xc xd xe xf aflow.out*
	rm -rf *exe aurostd_xtensor_extra.cpp aflow_data_*cpp xaflowauto boo* xlist* _xdat 
	rm -rfv `find . -name "LOCK*" | grep -v TESTS | grep -v SOURCES` *.out aflow.rasmol.xyz* aflow_potcar.cpp aflow_xpotcar.cpp aflow3.tar.xz ./AFLOWDATA ./LIBRARY* *html *toc *m *.jpg *.gif *.png
	rm -rfv ../AFLOW3_NONFREE/BIN/*
realclean: clean
	$(MAKE) clean_autogen
	rm -rf aurostd.cpp aurostd.h *conflic* 
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM realclean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then rm -fv ../AFLOW3_NONFREE/EXTRA/FROZSL/*conflic*; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR realclean; fi
clean_autogen:
	rm -rf $(AFLOWLIB_WEBAPP_CP) $(CHULL_JUPYTER_CPP) 
	rm -rf aflow_data_readme.cpp aflow_readme_aflow.cpp aflow_readme_aconvasp.cpp aflow_readme_apennsy.cpp aflow_readme_proto.cpp aflow_xproto_icsd_lib.cpp aflow_matlab_funcs.cpp aflow_gnuplot_funcs.cpp aflow_data_*cpp *.dat

# ****************************************************************************
# icsd
.PHONY: icsdclean icsd icsd_test
icsdclean: clean realclean
	echo "dummy" > $(ICSD_LIB)//dummy.txt
	rm $(ICSD_LIB)//*
icsd:	$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
icsd_test: aflow
	cat $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT | grep PROTOTYPE | sed "s/\[/\n/g" | grep PROTOTYPE | sed "s/PROTOTYPE/\.\/aflow --proto_icsd /g" > xgo
	sh ./xgo > /dev/null

# ****************************************************************************
# PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY
.PHONY: energy psenergy MIX mix NOMIX nomix EXPERIMENTS experiments TABLE table xaflowauto boo WEB APOOL_OLD shull vasp
energy:	all
	./aflow --energy > apennsy.tex
	latex apennsy.tex | grep -v Overfull | grep -v cmtt | grep -v type | grep -v "\[\]" | grep "\["
psenergy: all
	./aflow --psenergy > psapennsy.tex
	latex psapennsy.tex
	dvips psapennsy.dvi -o psapennsy.eps
	ps2pdf psapennsy.eps
	cp psapennsy*ps EPS
	cp psapennsy*pdf EPS
MIX:	aflow
	./aflow --libraryX --mix > aflow_nomix.$(TIME).cpp
	chmod 755 aflow_nomix.$(TIME).cpp
	rm -rf aflow_nomix.cpp
	ln -sf aflow_nomix.$(TIME).cpp aflow_nomix.cpp
	grep -c UNKNOWN aflow_nomix.2*
mix:	MIX
NOMIX:	aflow
	cat aflow_nomix.cpp | grep NOMIX | sed "s/\/\//\n/g" | grep calcs | sed "s/\=/\n/g" | grep calcs | sed "s/ calcs/xxxx/g" | sed "s/ /cd \.\//g" | sed "s/xxxx/ \&\& aflow --multi \&\& cd \.\.\//g" > xnomix
	scp xnomix fslcollab8@mary:LIBS/LIB2/LIB
	rm -rf xnomix
nomix: NOMIX
EXPERIMENTS:aflow
	./aflow --libraryX --experiments > aflow_mix_experiments.$(TIME).cpp
experiments:EXPERIMENTS
TABLE:aflow
	./aflow --libraryX --table > aflow_table.$(TIME).tex
table:TABLE
xaflowauto: aflow
	./aflow --libraryX --xaflowauto > xaflowauto
	scp xaflowauto fslcollab8@marylou5.byu.edu:LIBS
boo: aflow
	rm -rf boo
	./aflow --calculated=icsd --remove | sed "s/\.\//\.\/ICSD\//g" >> boo
	./aflow --calculated=lib0 --remove | sed "s/\.\//\.\/LIB0\/LIB\//g" >> boo
	./aflow --calculated=lib1 --remove | sed "s/\.\//\.\/LIB1\/LIB\//g" >> boo
	./aflow --calculated=lib2 --remove | sed "s/\.\//\.\/LIB2\/LIB\//g" >> boo
	./aflow --calculated=lib3 --remove | sed "s/\.\//\.\/LIB3\/LIB\//g" >> boo
	./aflow --calculated=lib4 --remove | sed "s/\.\//\.\/LIB4\/LIB\//g" >> boo
	./aflow --calculated=lib5 --remove | sed "s/\.\//\.\/LIB5\/LIB\//g" >> boo
	./aflow --calculated=lib6 --remove | sed "s/\.\//\.\/LIB6\/LIB\//g" >> boo
	./aflow --calculated=lib7 --remove | sed "s/\.\//\.\/LIB7\/LIB\//g" >> boo
	./aflow --calculated=lib8 --remove | sed "s/\.\//\.\/LIB8\/LIB\//g" >> boo
	./aflow --calculated=lib9 --remove | sed "s/\.\//\.\/LIB9\/LIB\//g" >> boo
	#scp boo aflow@lonsdale.tchpc.tcd.ie:LIB3/
	scp boo fslcollab8@mary:LIBS/
WEB:	aflow
	rm -rf *eps *pdf apennsy.m *html
	#	./aflow --libraryX --QUIET --print=pdf --hull > apennsy.m
	./aflow --libraryX --QUIET --print=gif --print=jpg --hull --oss=cout > apennsy.m
	./aflow --libraryX --QUIET --print=html > apennsy.html
	/usr/local/bin/matlab -nodesktop -r apennsy
	mv *jpg *gif /www/AFLOW
	rm -rf *eps
APOOL_OLD: aflow aflow_apennsy_main.cpp .cpp
	rm -rf /www/php/apool.html
	cat /www/php/apool_pre.html >> /www/php/apool.html
#	./aflow --apool_test >> /www/php/apool.html
	./aflow --apool >> /www/php/apool.html
	cat /www/php/apool_post.html >> /www/php/apool.html	
	./aflow --apool | grep available | sed "s/ checked//g" | sed "s/></XXXX\n/g" | grep radio | sed "s/<input type=\"radio\" name=\"alloy\" value=//g" | sed "s/XXXX/,/g" > /home/auro/work/AWRAPPER/awrapper_tables_public.cpp
	#
	#
#	rm -rf /www/php/apool2.html
#	cat /www/php/apool_pre.html >> /www/php/apool2.html
#	./aflow --apool_private >> /www/php/apool2.html
#	cat /www/php/apool_post.html >> /www/php/apool2.html	
#	#
shull:	all
	./aflow --library2 --shull --oss=cout > apennsys.m
vasp:	all
	./aflow --VASPIN | sed "s/\t/ /g" | sed "s/  / /g" | sed "s/  / /g" > structures.tex
