  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT AND OF THE ACCOMPANYING FILE LICENSE.txt.               *
  ***************************************************************************

   AUTHORS:

       Per Christian Hansen
       Technical University of Denmark, Denmark
       Email: pcha@dtu.dk

       Jakob Sauer Jorgensen
       University of Manchester, United Kingdom
       Email: jakob.jorgensen@manchester.ac.uk

   REFERENCE:

       AIR Tools II: Algebraic Iterative Reconstruction Methods,
       Improved Implementation
       NUMERICAL ALGORITHMS, XX (XXXX), pp. XXX-XXX
       DOI: XXXXXXXXXXXXXXXXXXXXXXXXXXX.

   SOFTWARE REVISION DATE:

       V1.0, October 2017

   SOFTWARE LANGUAGE:

       MATLAB 9.0 (R2016a)


=====================================================================
SOFTWARE
=====================================================================

The AIR Tools II toolbox for MATLAB provides efficient, robust and 
flexible implementations of algebraic iterative reconstruction (AIR) 
methods for computing regularized solutions to discretizations of 
inverse problems.

The AIR Tools II toolbox is an updated and expanded version of the 
AIR Tools toolbox from 2012.

The software has been developed and tested using MATLAB version 9.0.
No other MathWorks products or toolboxes are required.


=====================================================================
HOW TO INSTALL AND CHECK THE INSTALLATION
=====================================================================

Please follow these steps:

- Extract the archive in a directory of your choice. This creates a  
  directory called AIRToolsII with all files of the toolbox.

- Start MATLAB.

- Add the directory AIRToolsII containing the toolbox as well as all 
  subdirectories to the MATLAB search path.
  - The easiest way to do this is by using the function
    AIRToolsII_setup which updates the path either temporarily (for
    the current MATLAB session only - this is the default if no input
    is given) or permanently (for the current and all future MATLAB
    sessions).
  - The alternative is to use the "Set Path" tool in the "Home" tab 
    of the main MATLAB window, by clicking "Add with Subfolders..." 
    and choosing the AIRToolsII directory.

- The software installation can be checked by running demo scripts 
  from the demos subdirectory.


=====================================================================
SOFTWARE UPDATES AND BUG FIXES
=====================================================================

In addition to the refereed version of the software published along 
with the journal paper, the AIR Tools II software will be maintained 
in the GitHub code repository

  https://github.com/jakobsj/AIRToolsII

Please check this location for software updates and bug fixes.


=====================================================================
PACKAGE
=====================================================================

The AIRToolsII package is organized into a main directory with seven 
subdirectories:

AIRToolsII             contains this README.txt file, a LICENSE.txt 
                       file, a file CHANGES.txt listing major changes 
                       since previous versions, Contents.m which 
                       provides a detailed overview of all files in 
                       the package (can be listed from within MATLAB 
                       using "help AIRToolsII"), and an installation 
                       function called AIRToolsII_setup.m.
AIRToolsII/art         contains all ART methods.
AIRToolsII/cart        contains all CART methods.
AIRToolsII/sirt        contains all SIRT methods.
AIRToolsII/train       contains methods to train parameters.
AIRToolsII/testprobs   contains tomography test problems.
AIRToolsII/demos       contains demonstration scripts illustrating 
                       important functionality.
AIRToolsII/auxil       contains auxiliary files.

See the file Contents.m for a list of the files contained in the 
toolbox.


=====================================================================
DEMONSTRATION SCRIPTS
=====================================================================

The directory AIRToolsII/demos contains 12 demo scripts, which can be
executed for illustrating the basic use and functionality of the 
toolbox.

The demo demo_astra_2d.m illustrates how to interface from the AIR 
Tools II package to external software, in this case the ASTRA 
Tomography Toolbox, which must be installed separately, see 
  http://www.astra-toolbox.com/
for instructions.


=====================================================================
DOCUMENTATION
=====================================================================

More information can be found:

     - in the paper
       AIR Tools II: Algebraic Iterative Reconstruction Methods,
       Improved Implementation
       NUMERICAL ALGORITHMS, XX (XXXX), pp. XXX-XXX
       DOI: XXXXXXXXXXXXXXXXXXXXXXXXXXX.

     - in the GitHub code repository
       https://github.com/jakobsj/AIRToolsII

