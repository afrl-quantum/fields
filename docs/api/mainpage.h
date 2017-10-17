/*@HEADER
 *         fields:  A generic representation and mechanics for calculating
 *       fields/forces.
 *                 Copyright 2004-2008 Spencer Olson
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *  
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *                                                                                 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA.                                                                           .
 * 
 * Questions? Contact Spencer Olson (olsonse@umich.edu) 
 */

/** \mainpage fields


        <h3>Technical Reference</h3>
    - \f$\copyright\f$ 1998-2008 Spencer E. Olson


This package provides a generic representation of fields and forces as well as a
mechanism to calculate such in a generic fashion.

I would really appreciate any and all technical feedback, patches,
contributions.  See http://www.umich.edu/~olsonse/ for my best contact info.  



  This manual is divided in the following sections:
  - \subpage fields_intro 
  - \subpage fields_licence "License/Copying Information"
  - \subpage fields_readme
  - \subpage fields_install
  - \subpage fields_authors
*/



//-----------------------------------------------------------
/** \page fields_intro Introduction
    This library is for defining abstracted fields (from various arbitrary
    sources) and related forces.  Using C++ template programming, this library
    attempts to provide these abstracted interfaces while allowing for high
    performance computing.  In other words, this abstracted template library
    allows a computational programmer to develop various simulations using the
    inheritance and object oriented nature of C++ without the negative impact of
    using these features.
*/
    //Now you can proceed to the \ref advanced "advanced section".  



//-----------------------------------------------------------
/** \page fields_licence License/Copying Information
<p>
This software is released under the \ref fields_lgpl "LGPL v3" license except otherwise
explicitly stated in individual files included in this package. 
</p>

    - \f$\copyright\f$ 1998-2008 Spencer E. Olson
<p>
Relevant licenses:
    - \ref fields_lgpl "LGPL v3"
    - \ref gplv3 "GPL v3" 
    - \ref gplv2 "GPL v2" 
<p>

\section fields_lgpl ""
\verbinclude lgpl.txt

\section gplv3 ""
\verbinclude gpl-v3.txt

\section gplv2 ""
\verbinclude gpl-v2.txt

*/



//-----------------------------------------------------------
/** \page fields_readme README
    \verbinclude README
*/



//-----------------------------------------------------------
/** \page fields_install INSTALL
    \verbinclude INSTALL
*/



//-----------------------------------------------------------
/** \page fields_authors AUTHORS
    \verbinclude AUTHORS
*/
