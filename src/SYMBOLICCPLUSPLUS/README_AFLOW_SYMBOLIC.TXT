###########################################
 File name     : README_AFLOW_SYMBOLIC.TXT
 Description   : Describes the changes made to the SYMBOLICCPLUSPLUS/ subdirectory 
 Creation date : 2020-07-28
 Last modified : Tue 28 Jul 2020 12:48:47 PM EDT
 Created by    : David Hicks (DX)
###########################################

All changes from the original source are comments with "//DX<date>".

Modifications from original source code:

- The directory name was changed from SYMBOLICC++/ to SYMBOLICCPLUSPLUS/ since the former changes the Makefile formatting (i.e., colors for targets, etc.).

- Moved function definitions from the header files into CPP (.cpp) files.
  The original source code was all header (.h) files.
  This prevents inclusion of the code into more than one CPP file (i.e., violates the "One Definition Rule").

- The original source code had a class hierarchy (i.e., symbolic/ subdirectory) and a "redirect" header file.
  For clarity (and to ensure no circular/infinite dependency issues), the redirect was removed and all functions 
  were moved to a single directory hierarchy.

- Wrapped all functions in this directory in namespace "symbolic". Otherwise, this can cause namespace conflicts.

- Removed importing of the standard ("std") namespace.  Otherwise, this can cause namespace conflicts.

- Removed directory hierarchy in the SYMBOLICCPLUSPLUS/ subdirectory; otherwise this does not work with the current Makefile scheme.
  Specifically, the original hierarchy had a header redirect which cause issues with the automatic dependency tree in aflow_makefile.cpp.

- Initialization of variables to fix compiler errors/warnings.
