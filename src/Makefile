# ***************************************************************************
# *                                                                         *
# *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
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

VERSNUMBER=3.2.3

# MAKEFILE FOR AFLOW3
# FOR FROZSL -lmkl_lapack95_ilp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential

#system variables
# https://stackoverflow.com/questions/714100/os-detecting-makefile
UNAME=$(shell uname)

# for aflow
CPP=g++
#-O2 -frepo
# for sqllite
CC=gcc
ARCH=
ifeq ($(UNAME),Darwin)
CPP=clang++
CC=clang
ARCH=-e -D_MACOSX_
endif


#CPP=g++-fsf-7 
#CPP=g++-fsf-4.9 -std=gnu++11
# -std=gnu++11
# -O3
#CPP=g++ -std=c++0x-4.7 -std=c++0x
# -std=gnu++0x
#CCFLAGS= -DCYGWIN -Wall -W -fno-inline
CCFLAGS=-Wall -W -fno-inline
#[CO191217 - put in OPTS_MT]CCFLAGS_MT=$(CCFLAGS) -O2
#CCFLAGS_GSA=$(CCFLAGS) -O1	#CO OBSOLETE 170601
INCLUDE = -I ./
#INCLUDE = -I ./  -I /opt/intel/mkl100/include

# -fno-inline -fno-implicit-inline-templates-frepo
OPTS=
#-O1
#-O3
# -march=opteron
OPTS_MT=-O2 -std=gnu++0x
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

# ME190823 - SQLite options
# See https://www.sqlite.org/compile.html for compile-time options
# The link recommends -DSQLITE_THREADSAFE=0, but that would disable
# all threading capabilities
SQLITE_OPTS=-DSQLITE_THREADSAFE=1

# Library MKL dynamic link
#LIB_MKL= -L /opt/intel/mkl100/lib/32 -lmkl_lapack -lmkl_core -lmkl_intel -lmkl_sequential
# Library MKL static link
#LIB_MKL= -L /opt/intel/mkl100/lib/32/ /opt/intel/mkl100/lib/32/libmkl_lapack.a \
# /opt/intel/mkl100/lib/32/libmkl_core.a \
# /opt/intel/mkl100/lib/32/libguide.a

#CO20170614
#ME191125 -- Updated SQLite to version 3.30.1 (2019-10-11, https://www.sqlite.org/2019/sqlite-amalgamation-3300100.zip)
SQLITE_C = ./SQLITE/sqlite3.c

#README_SRC= aflow_readme.cpp

DATA_SRC=aflow_data_readme.cpp aflow_data_calculated.cpp aflow_data_htqc.cpp aflow_data_libraries.cpp aflow_data_stokes.cpp aflow_data_nist.cpp aflow_data_latex.cpp aflow_data_aflow_potcars.cpp aflow_data.cpp
DATA_OBJ                = $(DATA_FUNCS_SRC:.cpp=.o)

SQLITE_OBJ             = $(SQLITE_C:.c=.o)	#CO20170614

CHULL_JUPYTER_OBJS_ALL=aflow_chull_jupyter.cpp aflow_chull_jupyter_plotter.cpp aflow_chull_python.cpp aflow_chull_jupyter_requirements.cpp
OBJS_ALL= aflow.o $(AUROSTD_OBJ) $(AFLOW_CLASS_OBJ) $(AFLOW_FUNCS_OBJ) $(PFLOW_FUNCS_OBJ) $(ESTRUCTURE_FUNCS_OBJ) $(POCCUPATION_FUNCS_OBJ) $(APENNSY_FUNCS_OBJ) $(AUROSTD) $(LIBRARIES_FUNCS_OBJ) aflow_mix_experiments.o $(CONTRIB_FUNCS_OBJ) aflowlib_webapp_entry.cpp aflowlib_webapp_bands.cpp $(CHULL_JUPYTER_OBJS_ALL) $(AFLOWLIB_OBJ) $(AGL_OBJ) $(AEL_OBJ) $(APL_OBJ) $(ANRL_OBJ) $(COMPARE_FUNCS_DH_OBJ) $(CHULL_OBJ) $(GFA_OBJ) $(POCC_OBJ) $(CCE_OBJ) $(MAKEFILE_OBJ) $(SQLITE_OBJ)

.PHONY: check send_version_updates_email

#CO20200508 - all MUST come first
all: aflow aflow_data 
	$(MAKE) check

include Makefile.aflow

aflow: $(SQLITE_OBJ) $(DEPS_ALL) Makefile 
#	$(MAKE) -C ../AFLOW2_FREE/EXTRA/APL2AGR/ && cd ../../
#	cp -f aflow_apl*.o  aflow_aapl*.o  aflow_qha*.o ./APL/
#	cp -f aflow_anrl*.o ./ANRL/
#	cp sqlite3.o ./SQLITE/
#	$(CPP) $(LIBS_AFLOW) -o aflow *.o $(LIB_MKL) $(LIBS_AFLOW)
#	$(CPP) $(LIBS_AFLOW) -o aflow `ls *.o | grep -v aflow_data` $(LIB_MKL) $(LIBS_AFLOW)
	$(CPP) $(LIBS_AFLOW) -o aflow aflow_matlab_funcs.cpp aflow_gnuplot_funcs.cpp $(SQLITE_OBJ) $(OBJS_ALL) $(LIB_MKL) $(LIBS_AFLOW)
	ln -sf aflow aflowd
	ln -sf aflow apennsy
static_no_pdf: aflow aflow_data
	rm -rf $(SQLITE_OBJ)
	#$(CC) -DSQLITE_ENABLE_JSON1 -c $(SQLITE_C) $(LIBS) $(LIBS_STATIC)	#CO190401 - $(LIBS) -lpthreads not needed, issues with clang++ (mac)
	$(CC) $(SQLITE_OPTS) -c $(SQLITE_C) $(LIBS_STATIC)
#	cp sqlite3.o ./SQLITE/
	rm -f aflow
#	$(CPP) $(LIBS_AFLOW_STATIC) -o aflow *.o $(LIB_MKL) $(LIBS_AFLOW_STATIC)
	$(CPP) $(LIBS_AFLOW_STATIC) -o aflow `ls *.o | grep -v aflow_data` $(LIB_MKL) $(LIBS_AFLOW_STATIC)
	rm -f aflow_data
#	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -o aflow_data aflow_data.cpp $(LIBS_STATIC)
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -o aflow_data aflow_data.o aflow_data_*.o $(LIBS_STATIC)
#static: static_no_pdf
#	tar Jcfvp /www/AFLOW/CURRENT/aflow_$(VERSNUMBER)_amd64_debian_static.tar.xz aflow aflow_data	
$(SQLITE_OBJ): $(SQLITE_C)	#CO201706147
	$(CC) $(SQLITE_OPTS) -c $^ -o $@
#	$(CC) -DSQLITE_ENABLE_JSON1 -DSQLITE_OMIT_LOAD_EXTENSION -o $@ -c $^ $(LIBS)	#-DSQLITE_OMIT_LOAD_EXTENSION doesn't really work amalgamation
#	$(CC) -DSQLITE_ENABLE_JSON1 -o $@ -c $^ $(LIBS)
#	$(CC) -DSQLITE_ENABLE_JSON1 -c $^ $(LIBS)	#CO190401 - $(LIBS) -lpthreads not needed, issues with clang++ (mac)
#	$(CC) $(SQLITE_OPTS) -c $^
#	cp sqlite3.o ./SQLITE/
extra:	findsym frozsl apl2agr platon subst
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
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS; fi
extra_clean:
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/PLATON/51108 clean; fi
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR realclean; fi
#ME190313 - Change compiler to clang because g++ fails to compile on some Macs
mac:
	$(MAKE) CC=clang CPP=clang++ aflow_data
	$(MAKE) CC=clang CPP=clang++ -e ARCH=-D_MACOSX_ -j 4
	ln -sf aflow aflowd
	ln -sf aflow apennsy
	$(MAKE) check
#mac_static: mac
#	tar Jcfvp aflow_$(VERSNUMBER)_maxOSX_10.X.X.tar.xz aflow aflow_data

.cpp.o: .cpp Makefile $(AFLOW_HPP) $(AUROSTD_HPP) $(AUROSTD_CPP)
	$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\""$(@:.o=.cpp)"\" $(INCLUDE) $(CCFLAGS) $(OPTS) $(@:.o=.cpp) -c
	touch $@
#CO20170622 - START
aflowlib_webapp_entry.cpp: aflowlib_webapp_entry.js
	rm -rf $@
	echo "// $@ automatic generated" >> $@
	echo "std::string AFLOW_WEBAPP_ENTRY_JS=\"\\" >> $@
	cat $^ | perl -p -e 's/\r//g' | perl -p -e 's/\\u/\\\\u/g' | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
#CO20170622 - END
#CO20180305 - START
aflowlib_webapp_bands.cpp: aflowlib_webapp_bands.js
	rm -rf $@
	echo "// $@ automatic generated" >> $@
	echo "std::string AFLOW_WEBAPP_BANDS_JS=\"\\" >> $@
	cat $^ | perl -p -e 's/\r//g' | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
#CO20180305 - END
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
aflow_mix_nomix_table.cpp: Makefile $(AFLOW_HPP) $(AUROSTD_HPP)
	rm -rf nomix.tmp
	find $(PROJECT_DIRECTORY_GNDSTATE)/NOMIX/ -name aflow.in | sed "s/\/common\/$(PROJECT_DIRECTORY_GNDSTATE)\/NOMIX\/LIB2\/LIB\///g" | sed "s/aflow.in//g" >> nomix.tmp
	rm -rf aflow_mix_nomix_table.cpp
	echo "// aflow_mix_nomix_table.cpp automatic generated " >> aflow_mix_nomix_table.cpp
	echo "#include <string>" >> aflow_mix_nomix_table.cpp
	echo "std::string NomixTable=\"\\" >> aflow_mix_nomix_table.cpp
	cat nomix.tmp | sort | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_mix_nomix_table.cpp
	echo "\";" >> aflow_mix_nomix_table.cpp
	chmod 755 aflow_mix_nomix_table.cpp
	rm -rf nomix.tmp
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
	cp -f aflow_library_icsd.dat.xz ../AFLOW3_NONFREE/LIBS
	cp -f aflow_library_icsd.dat.xz /common/AFLOW/LIBS
	rm -f aflow_data_*
	$(MAKE) aflow_data
#	rm -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt*
	./aflow --protos_icsd > $(EXTRA_FREE)/NIST/ICSD_List_2014.txt
	rm -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz
	xz -9vf $(EXTRA_FREE)/NIST/ICSD_List_2014.txt
$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD1.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD2.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD3.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD4.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD5.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD6.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD7.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD8.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD9.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
############################################################################################################
aflow_data_calculated.cpp: Makefile
	rm -rf aflow_data_calculated.cpp
	# vLIBS	
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then echo "// aflow_data_calculated.cpp automatic generated from CALCULATED/vLIBS.dat" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then echo "#include <string>" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then echo "std::string vAUID,vAURL,vLOOP;" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then echo "std::string vLIBS=\"\\" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then xzcat $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz | grep -v "//" | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLIBS.dat.xz ]; then echo "\";" >> aflow_data_calculated.cpp ; fi
	# done
	-if [ ! -f aflow_data_calculated.cpp ]; then rm -rf aflow_data_calculated.cpp.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/aflow_data_calculated.cpp.xz && xz -dv aflow_data_calculated.cpp.xz ; fi
	-if [ ! -f aflow_data_calculated.cpp ]; then rm -rf aflow_data_calculated.cpp.xz*; wget http://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/aflow_data_calculated.cpp.xz && xz -dv aflow_data_calculated.cpp.xz ; fi
	if [ ! -f aflow_data_calculated.cpp ]; then exit 1; fi
	touch $@
aflow_data_htqc.cpp: README_LIBRARY_HTQC.TXT Makefile
	rm -f aflow_data_htqc.cpp
	echo "// aflow_data_htqc.cpp automatic generated from README_LIBRARY_HTQC.TXT" >> aflow_data_htqc.cpp
	echo "#include <string>" >> aflow_data_htqc.cpp
	echo "std::string Library_HTQC=\"\\" >> aflow_data_htqc.cpp
	cat README_LIBRARY_HTQC.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_htqc.cpp
	cat README_LIBRARY_HTQC_BORIDES.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_htqc.cpp
	cat README_LIBRARY_HTQC_CARBIDES.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_htqc.cpp
	echo "STRUCTURE ___PROTO_END___"  | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_htqc.cpp
	echo " END"  | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_htqc.cpp
	echo "\";" >> aflow_data_htqc.cpp
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
	rm -rf aflow_data_libraries.cpp
	echo "// aflow_data_libraries.cpp automatic generated from ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_libX.dat" >> aflow_data_libraries.cpp
	echo "#include <string>" >> aflow_data_libraries.cpp
#	echo "std::string aflowlib_icsd=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_icsd.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_icsd=\"\";" >> aflow_data_libraries.cpp
	echo "LIB0"
#	echo "std::string aflowlib_lib0=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib0.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib0=\"\";" >> aflow_data_libraries.cpp
	echo "LIB1"
#	echo "std::string aflowlib_lib1=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib1.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib1=\"\";" >> aflow_data_libraries.cpp
	echo "LIB2"
#	echo "std::string aflowlib_lib2=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib2.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib2=\"\";" >> aflow_data_libraries.cpp
	echo "LIB3"
#	echo "std::string aflowlib_lib3=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib3.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib3=\"\";" >> aflow_data_libraries.cpp
	echo "LIB4"
#	echo "std::string aflowlib_lib4=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib4.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib4=\"\";" >> aflow_data_libraries.cpp
	echo "LIB5"
#	echo "std::string aflowlib_lib5=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib5.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib5=\"\";" >> aflow_data_libraries.cpp
	echo "LIB6"
#	echo "std::string aflowlib_lib6=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib6.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib6=\"\";" >> aflow_data_libraries.cpp
	echo "LIB7"
#	echo "std::string aflowlib_lib7=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib7.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib7=\"\";" >> aflow_data_libraries.cpp
	echo "LIB8"
#	echo "std::string aflowlib_lib8=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib8.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib8=\"\";" >> aflow_data_libraries.cpp
	echo "LIB9"
#	echo "std::string aflowlib_lib9=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib9.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib9=\"\";" >> aflow_data_libraries.cpp
aflow_data_stokes.cpp: Makefile 
	rm -rf aflow_data_stokes.cpp
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "// aflow_data_stokes.cpp automatic generated " >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "#include <string>" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "std::string FINDSYM_data_space_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then echo "std::string FINDSYM_data_wyckoff_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then echo "std::string FROZSL_data_space_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then echo "std::string FROZSL_data_wyckoff_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then echo "std::string FROZSL_data_images_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g'  | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then echo "std::string FROZSL_data_irreps_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then echo "std::string FROZSL_data_isotropy_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then echo "std::string FROZSL_data_little_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then echo "std::string FROZSL_symmetry2_dat=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then echo "std::string FROZSL_const_dat=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then echo "std::string FROZSL_phvaspsetup_AFLOW=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then xzcat $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz | perl -p -e 's/\#/_ASCII_23_/g' | perl -p -e 's/\"/_ASCII_22_/g' | perl -p -e 's/\\/_ASCII_5C_/g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then echo "std::string FROZSL_phvaspsetup_POSCAR=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then xzcat $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz | perl -p -e 's/\#/_ASCII_23_/g' | perl -p -e 's/\"/_ASCII_22_/g' | perl -p -e 's/\\/_ASCII_5C_/g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	-if [ ! -f aflow_data_stokes.cpp ]; then rm -rf aflow_data_stokes.cpp.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/aflow_data_stokes.cpp.xz && xz -dv aflow_data_stokes.cpp.xz ; fi
	-if [ ! -f aflow_data_stokes.cpp ]; then rm -rf aflow_data_stokes.cpp.xz*; wget http://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/aflow_data_stokes.cpp.xz && xz -dv aflow_data_stokes.cpp.xz ; fi
	if [ ! -f aflow_data_stokes.cpp ]; then exit 1; fi
	touch $@
aflow_data_nist.cpp: Makefile 
	rm -rf aflow_data_nist.cpp
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "// aflow_data_nist.cpp automatic generated from NIST" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "#include <string>" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "std::string ElectronStoppingPower_txt=\"\\" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "\";" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then echo "std::string PhotonCrossSection_txt=\"\\" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then echo "\";" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then echo "std::string PhotonStoppingPower_txt=\"\\" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then echo "\";" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then echo "std::string ICSD_List_txt=\"\\" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then echo "\";" >> aflow_data_nist.cpp ; fi
	-if [ ! -f aflow_data_nist.cpp ]; then rm -rf aflow_data_nist.cpp.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/aflow_data_nist.cpp.xz && xz -dv aflow_data_nist.cpp.xz ; fi
	-if [ ! -f aflow_data_nist.cpp ]; then rm -rf aflow_data_nist.cpp.xz*; wget http://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/aflow_data_nist.cpp.xz && xz -dv aflow_data_nist.cpp.xz ; fi
	if [ ! -f aflow_data_nist.cpp ]; then exit 1; fi
	touch $@
aflow_data_aflow_potcars.cpp: Makefile 
	rm -rf aflow_data_aflow_potcars.cpp
	echo "// aflow_data_aflow_potcars.cpp automatic generated from AFLOW_PSEUDOPOTENTIALS.TXT" >> aflow_data_aflow_potcars.cpp
	echo "#include <string>" >> aflow_data_aflow_potcars.cpp
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then echo "bool AFLOW_PSEUDOPOTENTIALS=1;" >> aflow_data_aflow_potcars.cpp ; else echo "bool AFLOW_PSEUDOPOTENTIALS=0;" >> aflow_data_aflow_potcars.cpp ; fi
	echo "std::string AFLOW_PSEUDOPOTENTIALS_TXT=\"\\" >> aflow_data_aflow_potcars.cpp
	#if[ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_aflow_potcars.cpp ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT_XZ_BASE64.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT_XZ_BASE64.xz | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_aflow_potcars.cpp ; fi
	echo "\";" >> aflow_data_aflow_potcars.cpp
	echo "std::string AFLOW_PSEUDOPOTENTIALS_LIST_TXT=\"\\" >> aflow_data_aflow_potcars.cpp
	#if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz | grep "POTCAR=" | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_aflow_potcars.cpp ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS_LIST.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS_LIST.TXT.xz | grep "POTCAR=" | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_aflow_potcars.cpp ; fi
	echo "\";" >> aflow_data_aflow_potcars.cpp
aflow_data_latex.cpp: Makefile
	rm -rf aflow_data_latex.cpp
	echo "// aflow_data_latex.cpp automatic generated from LATEX" >> aflow_data_latex.cpp
	echo "#include <string>" >> aflow_data_latex.cpp
	echo "std::string f144468a7ccc2d3a72ba44000715efdb=\"\\" >> aflow_data_latex.cpp
	if [ -f ../AFLOW3_NONFREE/EXTRA/LATEX/f144468a7ccc2d3a72ba44000715efdb ]; then cat ../AFLOW3_NONFREE/EXTRA/LATEX/f144468a7ccc2d3a72ba44000715efdb* | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g'  | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_latex.cpp ; fi
	echo "\";" >> aflow_data_latex.cpp
	echo "std::string d0f1b0e47f178ae627a388d3bf65d2d2=\"\\" >> aflow_data_latex.cpp
	if [ -f ../AFLOW3_NONFREE/EXTRA/LATEX/d0f1b0e47f178ae627a388d3bf65d2d2 ]; then cat ../AFLOW3_NONFREE/EXTRA/LATEX/d0f1b0e47f178ae627a388d3bf65d2d2 | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g'  | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_latex.cpp ; fi
	echo "\";" >> aflow_data_latex.cpp
	echo "std::string decf00ca3ad2fe494eea8e543e929068=\"\\" >> aflow_data_latex.cpp
	if [ -f ../AFLOW3_NONFREE/EXTRA/LATEX/decf00ca3ad2fe494eea8e543e929068 ]; then cat ../AFLOW3_NONFREE/EXTRA/LATEX/decf00ca3ad2fe494eea8e543e929068 | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g'  | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_latex.cpp ; fi
	echo "\";" >> aflow_data_latex.cpp
aflow_data_readme.cpp: README_AFLOW_LICENSE_GPL3.TXT README_AFLOW_ACONVASP.TXT README_AFLOW_APL.TXT README_AFLOW_QHA_SCQHA_QHA3P.TXT README_AFLOW_AGL.TXT README_AFLOW_AEL.TXT README_AFLOW_POCC.TXT README_AFLOW_ANRL.TXT README_AFLOW_COMPARE.TXT README_AFLOW_GFA.TXT README_AFLOW.TXT README_PROTO.TXT README_AFLOW_APENNSY.TXT README_AFLOW_FROZSL.TXT README_AFLOW_SCRIPTING.TXT README_AFLOW_SYM.TXT README_AFLOW_CCE.TXT README_AFLOW_CHULL.TXT README_AFLOW_XAFLOW.TXT README_AFLOW_AFLOWRC.TXT README_AFLOW_EXCEPTIONS.TXT README_CONTRIBS.TXT README_AFLOW_VERSIONS_HISTORY.TXT Makefile
	rm -rf aflow_data_readme.cpp
	echo "// aflow_data_readme.cpp automatic generated" >> aflow_data_readme.cpp
	echo "#include <string>" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_LICENSE_GPL3_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_LICENSE_GPL3.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_VERSIONS_HISTORY_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_VERSIONS_HISTORY.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_PFLOW_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_ACONVASP.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_APENNSY_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_APENNSY.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_SCRIPTING_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_SCRIPTING.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_FROZSL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_FROZSL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_POCC_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_POCC.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_APL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_APL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_QHA_SCQHA_QHA3P_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_QHA_SCQHA_QHA3P.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_AGL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_AGL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_AEL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_AEL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_ANRL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_ANRL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_COMPARE_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_COMPARE.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_GFA_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_GFA.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_SYM_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_SYM.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_CCE_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_CCE.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_CHULL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_CHULL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_EXCEPTIONS_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_EXCEPTIONS.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_PROTO_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_PROTO.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_XAFLOW_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_XAFLOW.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_AFLOWRC_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_AFLOWRC.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
aflow_data_extra.cpp: Makefile aflow_library_icsd.dat.xz
	rm -rf aflow_data_extra.cpp aflow_data_libraries.cpp
	echo "#include <string>" >> aflow_data_libraries.cpp
	echo "std::string icsd_1_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/1-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_2_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/2-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_3_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/3-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_4_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/4-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_5_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/5-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_6_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/6-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_7_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/7-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_8_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/8-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_9_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/9-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_10_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/10-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string icsd_11_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/11-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string icsd_12_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/12-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string icsd_13_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/13-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string icsd_14_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/14-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string icsd_15_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/15-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
#	echo "std::string Library_ICSD=\"\\" >> aflow_data_extra.cpp && xzcat  ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string aflowlib_lib0=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib0.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib1=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib1.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib2=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib2.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib3=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib3.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib4=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib4.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib5=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib5.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib6=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib6.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib7=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib7.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib8=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib8.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib9=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib9.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_icsd=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_icsd.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
#aflow_data: aflow_data_readme.cpp aflow_data_calculated.cpp aflow_data_htqc.cpp aflow_data_libraries.cpp aflow_data_stokes.cpp aflow_data_nist.cpp aflow_data_latex.cpp aflow_data_aflow_potcars.cpp aflow_data.cpp
aflow_data: aflow_data_readme.o aflow_data_calculated.o aflow_data_htqc.o aflow_data_libraries.o aflow_data_stokes.o aflow_data_nist.o aflow_data_latex.o aflow_data_aflow_potcars.o aflow_data.cpp
	rm -f aflow_data_extra.cpp
	echo "" >> aflow_data_extra.cpp
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "#define defined_Library_ICSD" >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "std::string Library_ICSD=\"\\" >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then xzcat  ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "\";" >> aflow_data_extra.cpp; fi
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -c aflow_data.cpp
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -o aflow_data aflow_data.o aflow_data_*.o 
amir:	aflow_data_readme.cpp aflow_data_calculated.cpp aflow_data_htqc.cpp aflow_data_libraries.cpp aflow_data_stokes.cpp aflow_data_nist.cpp aflow_data_latex.cpp aflow_data_aflow_potcars.cpp aflow_data.cpp aflow_data_extra.cpp
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -DICSD -o aflow_data aflow_data.cpp
############################################################################################################
aflow_phonons.o: aflow_phonons.cpp $(AFLOW_HPP) $(AUROSTD_HPP)  Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_phonons.cpp -c
aflow_matlab_funcs.cpp: Makefile
	rm -rf aflow_matlab_funcs.cpp
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "// aflow_matlab_funcs.cpp automatic generated from MATLAB/*" >>  aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "#include <sstream>" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "// aflow_matlab_funcs.cpp automatic generated from MATLAB/param.m" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "std::string MATLAB_FUNCS_param(void){" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "   std::stringstream strstream; strstream << \"\\" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then cat $(EXTRA_FREE)/MATLAB/param.m | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "\" << std::endl; return strstream.str();};" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "// aflow_matlab_funcs.cpp automatic generated from MATLAB/plotband.m" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "// std::string MATLAB_FUNCS_plotband(std::string DIRECTORY,std::string OPTION1){" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "std::string MATLAB_FUNCS_plotband(void){" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "   std::stringstream strstream; strstream << \"\\" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then cat $(EXTRA_FREE)/MATLAB/plotband.m | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' | perl -p -e 's/MATLAB_FUNCS_DIRECTORY/\" << DIRECTORY << \"/g' | perl -p -e 's/DOS_SCALE_MODE/\" << OPTION1 << \"/g' >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "\" << std::endl; return strstream.str();};" >> aflow_matlab_funcs.cpp ; fi
	-if [ ! -f aflow_matlab_funcs.cpp ]; then rm -rf aflow_matlab_funcs.cpp.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/aflow_matlab_funcs.cpp.xz && xz -dv aflow_matlab_funcs.cpp.xz ; fi
	-if [ ! -f aflow_matlab_funcs.cpp ]; then rm -rf aflow_matlab_funcs.cpp.xz*; wget http://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/aflow_matlab_funcs.cpp.xz && xz -dv aflow_matlab_funcs.cpp.xz ; fi
	if [ ! -f aflow_matlab_funcs.cpp ]; then exit 1; fi
	touch $@
aflow_gnuplot_funcs.cpp: Makefile
	rm -rf aflow_gnuplot_funcs.cpp
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "// aflow_gnuplot_funcs.cpp automatic generated from GNUPLOT/*" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "#include <sstream>" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "// aflow_gnuplot_funcs.cpp automatic generated from GNUPLOT/plotbz.sh" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "//std::string GNUPLOT_FUNCS_plotbz(std::string DIRECTORY,std::string OPTION1){" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "std::string GNUPLOT_FUNCS_plotbz(void){" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "   std::stringstream strstream; strstream << \"\\" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then cat $(EXTRA_FREE)/GNUPLOT/plotbz.sh | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' | perl -p -e 's/GNUPLOT_FUNCS_DIRECTORY/\" << DIRECTORY << \"/g' | perl -p -e 's/GNUPLOT_FUNCS_OPTION1/\" << OPTION1 << \"/g' >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "\" << std::endl; return strstream.str();};" >> aflow_gnuplot_funcs.cpp ; fi
	-if [ ! -f aflow_gnuplot_funcs.cpp ]; then rm -rf aflow_gnuplot_funcs.cpp.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/aflow_gnuplot_funcs.cpp.xz && xz -dv aflow_gnuplot_funcs.cpp.xz ; fi
	-if [ ! -f aflow_gnuplot_funcs.cpp ]; then rm -rf aflow_gnuplot_funcs.cpp.xz*; wget http://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/aflow_gnuplot_funcs.cpp.xz && xz -dv aflow_gnuplot_funcs.cpp.xz ; fi
	if [ ! -f aflow_gnuplot_funcs.cpp ]; then exit 1; fi
	touch $@
aurostd_xtemplates.o: aurostd_xtemplates.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aurostd_xtemplates.cpp -c
#
###### HTRESOURCES #################################################
aflow_contrib_stefano_htcurriculum.o: $(AFLOWLIB_HPP) $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS)  -c
aflow_contrib_stefano_htsecurity.o: aflow_contrib_stefano_htsecurity.cpp $(AFLOWLIB_HPP) $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_contrib_stefano_htsecurity.cpp -c

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
send_version_updates_email:
	#CO190131 - START
	rm -f version_updates_email.txt
	$(MAKE) version_updates_email.txt
	[ -s version_updates_email.txt ]
	cat version_updates_email.txt | mail -r aflow@materials.duke.edu -s "AFLOW Code Update: v$(VERSNUMBER)" aflow+int+737+2740356048620877160@groups.io
	rm -f version_updates_email.txt
	#CO190131 - STOP
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
	done < README_AFLOW_VERSIONS_HISTORY.TXT; \
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
wget:	realclean
	rm -rf aflow3.tar.xz* aflow_test.tar.xz*
	-if [ ! -f aflow3.tar.xz ]; then wget http://materials.duke.edu/auro/aflow3.tar.xz; fi
	-if [ ! -f aflow3.tar.xz ]; then wget http://aflowlib.duke.edu/AFLOW/AFLOW3_FREE/aflow3.tar.xz; fi
	if [ ! -f aflow3.tar.xz ]; then exit 1; fi
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
check:
	$(MAKE) check_proto
check_full:
	#[WAITING FOR ICSD LIB2RAW REDO]$(MAKE) check_Egap
	$(MAKE) check_gcd
	$(MAKE) check_smith
	$(MAKE) check_proto
	$(MAKE) check_slab
	$(MAKE) check_chull
	$(MAKE) check_aflowSG
	$(MAKE) check_AFLOWSYM
	$(MAKE) check_aflow_xtal_match
	$(MAKE) check_aflow_gcc_compatibility
	$(MAKE) check_pocc
check_proto:
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
check_Egap:
	./aflow --Egap_test
check_gcd:
	./aflow --gcd_test
check_smith:
	./aflow --smith_test
check_slab:
	./aflow --slab_test
check_chull:
	@echo "Running CHULL test 1 - stability criterion of Mn2PdPt" && ANSWER=32 && SC=$$(./aflow --chull --alloy=MnPdPt --sc=aflow:fb9eaa58604ce774 --screen_only --o=t | cut -d' ' -f2) && \
		if [ -n "$${SC##$$ANSWER*}" ]; then echo "CHULL did not return a stability criterion of ~$$ANSWER meV/atom for Mn2PdPt './aflow --chull --alloy=MnPdPt --sc=aflow:fb9eaa58604ce774 --screen_only --o=t', found $$SC instead" && exit 1; fi
	@echo "CHULL tests passed"
check_pocc:
	@echo "Running POCC test 1 - unique structures count of C5HfMoNbTaTi" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:Ti --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 2 - unique structures count of C5HfMoNbTaV" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 3 - unique structures count of C5HfMoNbTaW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 4 - unique structures count of C5HfMoNbTaZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 5 - unique structures count of C5HfMoNbTiV" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 6 - unique structures count of C5HfMoNbTiW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 7 - unique structures count of C5HfMoNbTiZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 8 - unique structures count of C5HfMoNbVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 9 - unique structures count of C5HfMoNbVZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 10 - unique structures count of C5HfMoNbWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 11 - unique structures count of C5HfMoTaTiV" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 12 - unique structures count of C5HfMoTaTiW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 13 - unique structures count of C5HfMoTaTiZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 14 - unique structures count of C5HfMoTaVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 15 - unique structures count of C5HfMoTaVZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 16 - unique structures count of C5HfMoTaWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 17 - unique structures count of C5HfMoTiVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 18 - unique structures count of C5HfMoTiVZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 19 - unique structures count of C5HfMoTiWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 20 - unique structures count of C5HfMoVWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 21 - unique structures count of C5HfNbTaTiV" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 22 - unique structures count of C5HfNbTaTiW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 23 - unique structures count of C5HfNbTaTiZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 24 - unique structures count of C5HfNbTaVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 25 - unique structures count of C5HfNbTaVZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 26 - unique structures count of C5HfNbTaWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 27 - unique structures count of C5HfNbTiVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 28 - unique structures count of C5HfNbTiVZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 29 - unique structures count of C5HfNbTiWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 30 - unique structures count of C5HfNbVWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 31 - unique structures count of C5HfTaTiVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 32 - unique structures count of C5HfTaTiVZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 33 - unique structures count of C5HfTaTiWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 34 - unique structures count of C5HfTaVWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 35 - unique structures count of C5HfTiVWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 36 - unique structures count of C5MoNbTaTiV" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:Ti:V --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 37 - unique structures count of C5MoNbTaTiW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 38 - unique structures count of C5MoNbTaTiZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:Ti:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 39 - unique structures count of C5MoNbTaVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 40 - unique structures count of C5MoNbTaVZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 41 - unique structures count of C5MoNbTaWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 42 - unique structures count of C5MoNbTiVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 43 - unique structures count of C5MoNbTiVZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 44 - unique structures count of C5MoNbTiWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 45 - unique structures count of C5MoNbVWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 46 - unique structures count of C5MoTaTiVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 47 - unique structures count of C5MoTaTiVZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 48 - unique structures count of C5MoTaTiWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 49 - unique structures count of C5MoTaVWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ta:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 50 - unique structures count of C5MoTiVWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 51 - unique structures count of C5NbTaTiVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:Ti:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 52 - unique structures count of C5NbTaTiVZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:Ti:V:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 53 - unique structures count of C5NbTaTiWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 54 - unique structures count of C5NbTaVWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ta:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 55 - unique structures count of C5NbTiVWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Nb:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 56 - unique structures count of C5TaTiVWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Ta:Ti:V:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 57 - unique structures count of C5CrMoTaVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 58 - unique structures count of C5CrMoNbTaW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Nb:Ta:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 59 - unique structures count of C5CrMoNbVW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Nb:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 60 - unique structures count of C5CrHfTaWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Hf:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 61 - unique structures count of C5CrHfMoTiW" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Hf:Mo:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "Running POCC test 62 - unique structures count of C5CrMoTiWZr" && ANSWER=49 && COUNT=$$(./aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique | grep 'Structure bin' | tail -1 | awk -F 'bin ' '{print $$2}' | cut -d ':' -f1) && \
		if [ "$$COUNT" -ne "$$ANSWER" ]; then echo "POCC did not return $$ANSWER unique structures for './aflow --proto=AB_cF8_225_a_b.AB:C:Cr:Mo:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ./aflow --pocc_count_unique', found $$COUNT instead" && exit 1; fi
	@echo "POCC tests passed"
check_aflowSG:
	@echo "Running aflowSG test 1" && ANSWER="Pm-3m #221" && RESULT=$$(./aflow --proto=AB_cP2_221_b_a:Cl:Cs | ./aflow --aflowSG) && \
		if [ "$$RESULT" != "$$ANSWER" ]; then "aflowSG did not return $$ANSWER for './aflow --proto=AB_cP2_221_b_a:Cl:Cs | ./aflow --aflowSG', found $$RESULT instead" && exit 1; fi
	@echo "aflowSG test 1 passed!"
check_AFLOWSYM:
	@echo "Running aflowSYM test 1" && ANSWER_PGROUPS=8 && ANSWER_PGROUPKS=8 && ANSWER_FGROUPS=4 && ANSWER_PGROUP_XTALS=4 && ANSWER_PGROUPK_XTALS=4 && ANSWER_AGROUPS=24 && ANSWER_SGROUPS=108 && RESULT=$$(./aflow --proto=ABC2_mP8_10_ac_eh_mn-001:Au:Ag:Te | ./aflow --aflowSYM --screen_only) && COUNT_PGROUPS=$$(echo $$RESULT | awk -F'pgroup ' '{print NF-1}') && COUNT_PGROUPKS=$$(echo $$RESULT | awk -F'pgroupk ' '{print NF-1}') && COUNT_FGROUPS=$$(echo $$RESULT | awk -F'fgroup ' '{print NF-1}') && COUNT_PGROUP_XTALS=$$(echo $$RESULT | awk -F'pgroup_xtal ' '{print NF-1}') && COUNT_PGROUPK_XTALS=$$(echo $$RESULT | awk -F'pgroupk_xtal ' '{print NF-1}') && COUNT_AGROUPS=$$(echo $$RESULT | awk -F'agroup ' '{print NF-1}') && COUNT_SGROUPS=$$(echo $$RESULT | awk -F'sgroup ' '{print NF-1}') && \
		if [ "$$ANSWER_PGROUPS" -ne "$$COUNT_PGROUPS" ] || [ "$$ANSWER_PGROUPKS" -ne "$$COUNT_PGROUPKS" ] || [ "$$ANSWER_FGROUPS" -ne "$$COUNT_FGROUPS" ] || [ "$$ANSWER_PGROUP_XTALS" -ne "$$COUNT_PGROUP_XTALS" ] || [ "$$ANSWER_PGROUPK_XTALS" -ne "$$COUNT_PGROUPK_XTALS" ] || [ "$$ANSWER_AGROUPS" -ne "$$COUNT_AGROUPS" ] || [ "$$ANSWER_SGROUPS" -ne "$$COUNT_SGROUPS" ]; then "aflowSYM did not return one of the correct group cardinalities for './aflow --proto=ABC2_mP8_10_ac_eh_mn-001:Au:Ag:Te | ./aflow --aflowSYM --screen_only'; (ANSWER==FOUND) $$ANSWER_PGROUPS==$$COUNT_PGROUPS, $$ANSWER_PGROUPKS==$$COUNT_PGROUPKS, $$ANSWER_FGROUPS==$$COUNT_FGROUPS, $$ANSWER_PGROUP_XTALS==$$COUNT_PGROUP_XTALS, $$ANSWER_PGROUPK_XTALS==$$COUNT_PGROUPK_XTALS, $$ANSWER_AGROUPS==$$COUNT_AGROUPS, $$ANSWER_SGROUPS==$$COUNT_SGROUPS"  && exit 1; fi
	@echo "aflowSYM test 1 passed!"
check_aflow_xtal_match:
	@echo "Running aflow_xtal_match test 1" && ANSWER=1 && RESULT=$$(./aflow --proto=A2B3_hR10_167_c_e:Al:O --params=-1,2.72957758313,0.36,0.5561 | ./aflow --compare2protos --print=txt --catalog=anrl --quiet) && NUM_MATCHING_PROTOS=$$(echo $$RESULT | awk -F'duplicate_compounds=' '{print $$2}' | awk '{print $$1}') && \
		if [ "$$NUM_MATCHING_PROTOS" -ne "$$ANSWER" ]; then "aflow_xtal_match did not return $$ANSWER for './aflow --proto=A2B3_hR10_167_c_e:Al:O --params=-1,2.72957758313,0.36,0.5561 | ./aflow --compare2protos --print=txt --catalog=anrl --quiet', found $$NUM_MATCHING_PROTOS instead" && exit 1; fi
	@echo "aflow_xtal_match test 1 passed!"
check_aflow_gcc_compatibility:
	@echo "Running aflow_gcc_compatibility test 1: no 'nullptr'" && ANSWER=0 && NUMBER_NULLPTRS=$$(grep -P '^(?=[\s]*+[^//])[^//]*(nullptr)' *.cpp | wc -l) && \
		if [ "$$NUMBER_NULLPTRS" -ne "$$ANSWER" ]; then echo "Found $$NUMBER_NULLPTRS in the source code; this is not backwards compatible with older GCC version. Please swap with 'NULL' and recompile." && exit; fi
	@echo "aflow_gcc_compatibility test 1 passed!"
www:	web
clean:
	rm -rf *~ *.o *.o-* *.rpo AUROSTD/*~ AUROSTD/*.o AUROSTD/*.rpo ./APL/*~ ./APL/*o ./APL/*a ./ANRL/*~ ./ANRL/*o ./ANRL/*a ./ANRL*/*~ ./ANRL*/*o ./ANRL*/*a ./SQLITE/*o aflow aflow_data aflowd aflow2 aflowd2 aflow1 aflowd1 aflow.exe aconvasp apennsy aflow_convasp* *.eps *.pdf *ps *aux *dvi *log *m *tex nohup* core* VASP_TEST* GRND_TEST* xa xb xc xd xe xf aflow.out* core* *~ \#* *.o aflow_pennsy *.eps *.pdf *ps *aux *dvi *log *m *tex *exe aurostd_xtensor_extra.cpp aflow_data_*cpp xaflowauto boo* xlist* _xdat `find . -name *~`
	rm -rfv `find . -name "LOCK*" | grep -v TESTS | grep -v SOURCES` *.out aflow.rasmol.xyz* aflow_potcar.cpp aflow_xpotcar.cpp aflow3.tar.xz ./AFLOWDATA ./LIBRARY* *html *toc *m *.jpg *.gif *.png
	rm -rfv ./APL/*~ ./ANRL/*~ ./ANRL*/*~ ./SQLITE/*~ ../AFLOW3_NONFREE/BIN/* ./EXTRA/GUS/*~ ./EXTRA/GUS/*.o 
realclean: clean
	$(MAKE) clean_autogen
	rm -rf aurostd.cpp aurostd.h *conflic* 
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM realclean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then rm -fv ../AFLOW3_NONFREE/EXTRA/FROZSL/*conflic*; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR realclean; fi
dataclean:
	rm -f aflow_data aflow_data_*
data: aflow_data
clean_autogen:
	rm -rf aflowlib_webapp_entry.cpp aflowlib_webapp_bands.cpp aflow_chull_jupyter.cpp aflow_chull_jupyter_plotter.cpp aflow_chull_python.cpp aflow_chull_jupyter_requirements.cpp aflow_data_readme.cpp aflow_readme_aflow.cpp aflow_readme_aconvasp.cpp aflow_readme_apennsy.cpp aflow_readme_proto.cpp aflow_xproto_icsd_lib.cpp aflow_matlab_funcs.cpp aflow_gnuplot_funcs.cpp aflow_data_*cpp *.dat
# ****************************************************************************
# icsd
icsdclean: clean realclean
	echo "dummy" > $(ICSD_LIB)//dummy.txt
	rm $(ICSD_LIB)//*
icsd:	$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
icsd_test: aflow
	cat $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT | grep PROTOTYPE | sed "s/\[/\n/g" | grep PROTOTYPE | sed "s/PROTOTYPE/\.\/aflow --proto_icsd /g" > xgo
	sh ./xgo > /dev/null
# ****************************************************************************
# PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY
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
