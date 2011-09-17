// -*- c++ -*-
// $Id$
/*@HEADER
 *         olson-tools:  A variety of routines and algorithms that
 *      I've developed and collected over the past few years.  This collection
 *      represents tools that are most useful for scientific and numerical
 *      software.  This software is released under the LGPL license except
 *      otherwise explicitly stated in individual files included in this
 *      package.  Generally, the files in this package are copyrighted by
 *      Spencer Olson--exceptions will be noted.   
 *                 Copyright 2004-2008 Spencer E. Olson
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

/** \file
 * Scale a Force in a time dependent manner.
 *
 * Copyright 2004-2008 Spencer Olson
 *
 */



#ifndef fields_ScaleForce_h
#define fields_ScaleForce_h

#include <fields/Forces.h>
#include <fields/make_options.h>

#include <xylose/timing/Timing.h>
#include <xylose/timing/element/Exponential.h>

#include <limits>

namespace fields {

  namespace timing = xylose::timing;

  template < typename Force >
  struct ScaleForce
    : virtual BaseForce<typename Force::super0::options>, Force {
    /* TYPEDEFS */
    typedef BaseForce<typename Force::super0::options> super0;
    typedef Force F;

    /* NON-MEMBER STORAGE */
  private:
    /** Default timing element applys a unity scaling. */
    static timing::element::Exponential * mkDefaultTiming() {
      return new timing::element::Exponential(
        -std::numeric_limits<double>::infinity(), 1.0, 1.0, 1.0
      );
    }

    /* MEMBER STORAGE */
  public:
    /** Timing function for this Field scaling. */
    timing::Timing timing;


    /* MEMBER FUNCTIONS */
    /** Constructor adds default timing element to timing. */
    ScaleForce() : super0(), F(), timing() {
      /* default to having no timing effect. */
      timing.timings.push_back( mkDefaultTiming() );
      timing.set_time(0.0);
    }

    /** Assignment operator. */
    inline const ScaleForce & operator= ( const ScaleForce & that ) {
      super0::operator=(that);
      F::operator=(that);
      timing = that.timing;
      return *this;
    }

    /** Calculate acceleration. */
    inline void accel(       Vector<double,3> & a,
                       const Vector<double,3> & r,
                       const Vector<double,3> & v = V3(0.,0.,0.),
                       const double & t = 0.0,
                       const double & dt = 0.0,
                       const unsigned int & species = 0u ) const {
      F::accel(a,r,v,t,dt,species);
      a *= timing.getVal();
    }

    template < typename P >
    inline void accel(       Vector<double,3> & a,
                       const Vector<double,3> & r,
                       const Vector<double,3> & v,
                       const double & t,
                       const double & dt,
                             P & p ) const {
      F::accel(a,r,v,t,dt,p);
      a *= timing.getVal();
    }

    inline double potential( const Vector<double,3> & r,
                             const Vector<double,3> & v = V3(0,0,0),
                             const double & t = 0.0,
                             const unsigned int & species = 0u ) const {
      return timing.getVal() * F::potential(r,v,t,species);
    }

    template < typename P >
    inline double potential( const Vector<double,3> & r,
                             const Vector<double,3> & v,
                             const double & t,
                                   P & p ) const {
      return timing.getVal() * F::potential(r,v,t,p);
    }

    template < unsigned int ndim,
               typename Particle >
    inline void applyStatisticalForce(       Vector<double,ndim> & xv,
                                       const double & t,
                                       const double & dt,
                                             Particle & particle ) const {
      F::applyStatisticalForce( xv, t, timing.getVal() * dt, particle );
    }
  };

}/* namespace fields */

#endif //fields_ScaleForce_h
