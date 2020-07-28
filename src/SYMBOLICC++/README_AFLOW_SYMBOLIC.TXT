###########################################
 File name     : README
 Description   : Describes the changes made to the SYMBOLICC++/ subdirectory 
 Creation date : 2020-07-28
 Last modified : Tue 28 Jul 2020 12:48:47 PM EDT
 Created by    : David Hicks (DX)
###########################################

All changes from the original source are comments with "//DX<date>".

Modifications from original source code:

- Wrapped all functions in this directory in namespace "symbolic". Otherwise, this can cause namespace conflicts.

- Removed importing of the standard ("std") namespace.  Otherwise, this can cause namespace conflicts.

- Removed directory hierarchy in the SYMBOLICC++/ subdirectory; otherwise this does not work with the current Makefile scheme.
  Specifically, the original hierarchy had a header redirect which cause issues with the automatic dependency tree in aflow_makefile.cpp.

- Initialization of variables to fix compiler errors/warnings.
