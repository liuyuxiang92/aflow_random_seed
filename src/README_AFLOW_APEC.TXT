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

AFLOW-APEC is an automated thermodynamic stability analysis tool. Compound data is downloaded from the
AFLOW REST-API/AFLUX Search-API and the minimum energy surface is calculated using a modified
variant (dimensionally-iterative) of the qhull convex hull algorithm:
ACM Trans. on Mathematical Software, 22(4):469-483, Dec 1996.
Various thermodynamic properties are extracted, including off-hull (full decomposition reaction,
distances to the hull) and on-hull (equilibrium phases, equivalent ground state structures) features.
Determination of the stability criterion is also available, which has served as a powerful
screening tool for the discovery of new permanent magnets: Science Adv. 3(4), e1602241 (2017).

Other important analyses running under the hood include:
 - structure comparison analysis via Burzlaff criteria (--readme=compare)
 - outlier detection and removal (interquartile range).

Various output formats are provided, including plain text, json, PDF, and Jupyter.

****************************************************************************************************
LIST OF AFLOW-APEC COMMANDS:


