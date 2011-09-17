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

#ifndef fields_Forces_h
#define fields_Forces_h

#include <fields/Fields.h>
#include <fields/indices.h>
#include <fields/make_options.h>

#include <chimp/property/mass.h>

#include <xylose/Vector.h>

namespace fields {

  using xylose::Vector;
  using xylose::V3;

  /** Wrapper class to provide 'derivs' functionality mainly for use with
   * integrating functions.  Note that this class should generally NOT be used
   * with fields::AddForce and like classes since this does not inherit
   * BaseForce virtually.
   *
   * This version allows for a velocity dependence in the force. 
   * @see xylose::integrate::RK.
   */
  template < typename Force >
  struct Derivs : public Force {
    /** Compute phase-space derivatives for integrating methods (e.g.
     * runge-kutta.
     * Because this is a template, you will have to explicity instantiate this
     * function to get a pointer to pass into rk or the like.
     */
    template < typename Particle >
    void operator() ( const Vector<double,6u> & x,
                      const double & time,
                      const double & dt,
                            Vector<double,6u> & rkF,
                      const Particle & p ) {
      const Vector<double,3> & r = V3C(x.val);
      const Vector<double,3> & v = V3C(x.val + VX);
            Vector<double,3> & a = V3C(rkF.val + VX);

      using namespace indices;
      rkF[X]  = x[VX];
      rkF[Y]  = x[VY];
      rkF[Z]  = x[VZ];
      Force::accel( a, r, v, time, dt, p );
    }
  };

  /** Base of all force classes;  This class provides access to chimp.
   *
   * All forces need to implement one or both of the accel functions.  Which one
   * they implement (or both) will determine where they can be used.  A force
   * can certainly ignore OptionalParticle if it is unnecessary, or it can also
   * assume some appropriate default value.
   * <code>inline void accel(       Vector<double,3> & a,
   *                          const Vector<double,3> & r,
   *                          const double & t,
   *                          const double & dt,
   *                          const unsigned int & species ) const = 0;
   * </code><br>
   * <code>template < typename OptionalParticle >
   *       inline void accel(       Vector<double,3> & a,
   *                          const Vector<double,3> & r,
   *                          const double & t,
   *                          const double & dt,
   *                                OptionalParticle & ) const = 0;
   * </code><br>
   * All forces need to implement one or both of the potential functions.
   * <code>template < typename OptionalParticle >
   *       inline double potential( const Vector<double,3> & r,
   *                                const double & t,
   *                                OptionalParticle & ) const = 0;
   * </code><br>
   * <code>inline double potential( const Vector<double,3> & r,
   *                                const double & t ),
   *                                const unsigned int & species ) const = 0;
   * </code><br>
  */
  template < typename _options = fields::make_options<>::type >
  struct BaseForce {
    typedef _options options;
    /** Pointer to chimp instance. */
    typename options::ChimpDB * db;

    BaseForce() : db(NULL) { }
  };

  template < typename options = fields::make_options<>::type >
  class Gravity : public virtual BaseForce<options>,
                  public BgField< Vector<double,3> > {
  public:
    typedef BaseForce<options> super0;
    typedef BgField< Vector<double,3> > super1;

    Gravity() : super0(), super1() {
      bg = V3(0.,0.,-9.81);
    }

    void accel(       Vector<double,3> & a,
                const Vector<double,3> & r = V3(0.,0.,0.),
                const Vector<double,3> & v = V3(0.,0.,0.),
                const double & t = 0.0,
                const double & dt = 0.0,
                const unsigned int & species = 0u ) const {
      a = super1::bg;
    }

    template < typename P >
    void accel(       Vector<double,3> & a,
                const Vector<double,3> & r,
                const Vector<double,3> & v,
                const double & t,
                const double & dt,
                      P & ) const {
      this->accel( a, r, v, t, dt );
    }

    /** Calculate the potential of \f$^{87}{\rm Rb}\f$ |F=1,mF=-1>.
     * Gravitational energy is referenced to (0,0,0).
     */
    double potential( const Vector<double,3> & r,
                      const Vector<double,3> & v = V3(0.,0.,0.),
                      const double & t = 0.0,
                      const unsigned int & species = 0u ) const {
      using chimp::property::mass;
      return - (*super0::db)[species].mass::value * (super1::bg * r);
    }

    template < typename P >
    double potential( const Vector<double,3> & r,
                      const Vector<double,3> & v,
                      const double & t,
                            P & p ) const {
      return this->potential( r, v, t, species(p) );
    }

    template < unsigned int ndim,
               typename Particle >
    void applyStatisticalForce(       Vector<double,ndim> & xv,
                                const double & t,
                                const double & dt,
                                      Particle & particle ) const { }
  };

  /** Adds Forces to get a total acceleration/ potential energy.
   * Note that this is only really helpful for physically disjoint forces.  It
   * will not be physically correct to add two forces due to magnetic fields for
   * example.  In this case, one must add the fields together and then compute
   * the forces with BCalcs.
   *
   * An an example use of this class, let's add the forces each due to gravity and
   * magnetic fields:
   *      chimp_instance.addParticleType("87Rb");
   *      ...
   *      typedef AddForce< BCalcs< BSrc >, Gravity > myForce;
   *      myForce force;
   *      force.db = &chimp_instance;//both forces use same chimp instance.
   *
   * Note that if F0 and F1 are from fields that use BaseField, this template
   * class will NOT cause their BaseField::delta values to be shared.
   *
   * @see BField::BCalcs.
   */
  template < class _F0, class _F1 >
  class AddForce : public virtual BaseForce<typename _F0::super0::options>,
                   public _F0, public _F1 {
  public:
    /* using _F1 here is intentional to help make sure that both forces are
     * using the same BaseForce::options. */
    typedef BaseForce<typename _F1::super0::options> super0;
    typedef _F0 F0;
    typedef _F1 F1;

    AddForce() : super0(), F0(), F1() {}

    const AddForce & operator=(const AddForce & that) {
      super0::operator=(that);
      F0::operator=(that);
      F1::operator=(that);
      return *this;
    }

    void accel(       Vector<double,3> & a,
                const Vector<double,3> & r,
                const Vector<double,3> & v = V3(0,0,0),
                const double & t = 0.0,
                const double & dt = 0.0,
                const unsigned int & species = 0u ) const {
      F0::accel(a,r,v,t,dt,species);
      Vector<double,3> a2;
      F1::accel(a2,r,v,t,dt,species);
      a +=  a2;
    }

    template < typename P >
    void accel(       Vector<double,3> & a,
                const Vector<double,3> & r,
                const Vector<double,3> & v,
                const double & t,
                const double & dt,
                      P & p ) const {
      F0::accel(a,r,v,t,dt,p);
      Vector<double,3> a2;
      F1::accel(a2,r,v,t,dt,p);
      a +=  a2;
    }

    double potential( const Vector<double,3> & r,
                      const Vector<double,3> & v = V3(0,0,0),
                      const double & t = 0.0,
                      const unsigned int & species = 0u ) const {
      return F0::potential(r,v,t,species) + F1::potential(r,v,t,species);
    }

    template < typename P >
    double potential( const Vector<double,3> & r,
                      const Vector<double,3> & v,
                      const double & t,
                            P & p ) const {
      return F0::potential(r,v,t,p) + F1::potential(r,v,t,p);
    }

    template < unsigned int ndim,
               typename Particle >
    void applyStatisticalForce(       Vector<double,ndim> & xv,
                                const double & t,
                                const double & dt,
                                      Particle & particle ) const {
      F0::applyStatisticalForce( xv, t, dt, particle );
      F1::applyStatisticalForce( xv, t, dt, particle );
    }
  };


 /** A wrapper class for a statistical force for using with RK5 driver. */
 template <class _StatisticalForce>
 class StatisticalForceRKWrapper {
  public:
    StatisticalForceRKWrapper() : statisticalForcePtr(NULL) {}
    _StatisticalForce * statisticalForcePtr;

    template < unsigned int ndim,
               typename Particle >
    void first(       Vector<double,ndim> & xv,
                const double & t,
                const double & dt_step_current,
                      double & dt_step_next,
                      Particle & particle ) {
      statisticalForcePtr
        ->applyStatisticalForce(xv, t, 0.5*dt_step_current, particle);
    }

    template < unsigned int ndim,
               typename Particle >
    void second(       Vector<double,ndim> & xv,
                 const double & t,
                 const double & dt_step_current,
                       double & dt_step_next,
                       Particle & particle ) {
      statisticalForcePtr
        ->applyStatisticalForce(xv, t, 0.5*dt_step_current, particle);
    }
  };

}/* namespace fields */

#endif // fields_Forces_h
