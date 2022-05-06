****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2022                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                       materials.duke.edu/AFLOW/README_AFLOW_APEC.TXT
****************************************************************************************************

AFLOW-APEC README:

README written by: Simon Divilov (simon.divilov@duke.edu)

Citation info: TBA





****************************************************************************************************
GENERAL OVERVIEW

AFLOW-APEC is an automated tool to calculate phase equilibria properties of alloys for a given 
parent lattice. Alloy data is downloaded from the AFLOW database and the cluster expansion is 
performed to calculcate the energies of additional derivatice structures using the Alloy Theoretic 
Automated Toolkit [J. Phase Equilib. 23, 348 (2002)].

Using the generalized quasi-chemical approximation, the binodal curve is calculcated, based on the
LTVC method [Acta Mater. 159, 364 (2018)].

****************************************************************************************************
LIST OF AFLOW-APEC COMMANDS:

  aflow --phase_equilibria --plattice=fcc --elements=Au,Pt[,Zn] [apec_options] [--directory=[DIRECTORY]]
      Loads the alloy data from the AFLOW database and constructs the binodal curve and returns,
      the output as a TXT (default, see --format_data).
      --phase_equilibria : Necessary argument, enter the mode to calculate the phase equilibria
          properties.
      --plattice : Necessary argument, specify the parent lattice of the alloy
      --elements : Necessary argument, specify the elements present in the alloy
      [--directory=[DIRECTORY]] : Optional argument, specify the directory where to run the
          calculation.
          Default: "./"

      GENERAL OPTIONS:
      [--usage] : Returns usage commands and options.
          Default: OFF
      [--screen_only] : Output directly to screen (stdout) instead of writing to a file.
          Default: OFF
      [--image_only|--image] : Only read the output data and plot the results. Plotting
          is only implemented for binary alloys.
          Default: OFF
      [--no_plot|--noplot] : Do not plot the results after the calculcation is complete.
          Default: OFF

      BINODAL OPTIONS:
      [--max_num_atoms=|--mna=8] : Specify the maximum number of atoms in the alloy
          Default: 8
      [--conc_curve_range=|--conc_curve=0,1,1,0] : Define the concentration curve
          to evalute the equilibria properties. The format of this argument is:
          <c_1i>,<c_c1f>,<c_2i>,<c_2f>,...,<c_[K]i>,<c_[K]f> where <c_[k]i> is the 
          initial concentration and <c_[k]f> is the final concentration.
          Default: NONE
      [--conc_npts=20] : Number of points used to determine the sub-division of 
          the concentration profile.
          Default: 20
      [--temp_range=300,2000] : Define the temperature range to evaluate the
          equilibria properties. The format of this argument is: <T_i>,<T_f>
          where <T_i> and <T_f> is the initial and final temperature, respectively.
          Default: 300,2000
      [--temp_npts=100] : Number of points used to determine the sub-division of
          the temperature profile.
          Default: 100

      SPINODAL OPTIONS:
      [--spinodal] : Calculate the spinodal curve (NOT IMPLEMENTED)
          Default: OFF

      FORMAT OPTIONS:
      [--format_data=|--data_format=txt|json] : Specify the format of the output
          data.
          Default: txt
      [--format_image=|--image_format=pdf|eps|png] : Specify the format of the
          image produced by the plotting function.
          Default: pdf

For additional information contact: Simon Divilov (simon.divilov@duke.edu)

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2022                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
