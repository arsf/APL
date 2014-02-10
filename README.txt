General Notes
-------------

The Airborne Processing Library (APL) is a toolbox for processing hyperspectral data. It has initially been created by the Airborne Research Survey Facility (ARSF) for processing the data collected from the ARSF platform. As such, to use it for non-ARSF data, may require some extensive coding. For information on software use and capabilities please see the ARSF Data Processing wiki: http://arsf-dan.nerc.ac.uk/trac/wiki/Help

Licencing information
---------------------

This software is available under a licence derived from the Non-Profit Open Software License version 3.0 (http://opensource.org/licenses/NPOSL-3.0), modifying it to exclude commercial uses of the original APL software while permitting non-profit use, retaining open access to source code and allowing contributions.  If you wish to use the original APL software commercially, please contact NERC ARSF (email: arsf@nerc.ac.uk cc'ing arsf-processing@pml.ac.uk) or, in the event of ARSF ceasing to exist, NERC or its successor organisations directly.  Note that derivative works cannot be released from these non-commercial provisions; only copies of APL where NERC has full ownership can be relicensed under other terms. Please see LICENCE.txt for full licencing details.

External packages (dependencies)
--------------------------------

Parts of the apl-suite make use of the PROJ.4 library of cartographic projections. This is required for a fully working version of the full apl-suite. If your system does not already have PROJ installed then it can be obtained from http://trac.osgeo.org/proj/wiki. Tested with PROJ.4 version 4.7.1.

Parts of the apl-suite make use of the Blitz++ library for matrix manipulation. This is required for a fully working version the full apl-suite. If your system does not already have Blitz++ it can be obtained from http://sourceforge.net/projects/blitz/. Tested with Blitz++ version 0.9 - currently incompatible with version 0.10.

Compiling on Fedora Linux
-------------------------

Basic example Makefiles (for linux and Windows versions) have been included for building (locally) the APL executables on Fedora systems. These will likely need editing to point to specific directories on your system, compilers etc.

Referencing
-----------

When acknowledging the use of APL for scientific papers, reports etc please cite the following reference:
M. A. Warren, B. H. Taylor, M. G. Grant, J. D. Shutler, Data processing of remotely sensed airborne hyperspectral data using the Airborne Processing Library (APL): Geocorrection algorithm descriptions and spatial accuracy assessment, Computers & Geosciences, Volume 64, March 2014, Pages 24-34, ISSN 0098-3004, http://dx.doi.org/10.1016/j.cageo.2013.11.006.


