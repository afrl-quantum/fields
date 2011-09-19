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
#include <fields/detail/assert.h>
#include <fields/detail/NullForce.h>

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
    const typename options::ChimpDB * db;

    BaseForce( const typename options::ChimpDB * db = NULL ) : db(db) { }
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

  /** Adds several Forces to get a total acceleration/ potential energy.
   * This metafunction can only be used ONCE to aggregate forces.  To add
   * additional forces, use fields::AddForce.
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
   *      myForce force(& chimp_instance);
   *        // or init chimp_instance separately...
   *      force.db = &chimp_instance;//both forces use same chimp instance.
   *
   * Note that if F0 and F1 are from fields that use BaseField, this template
   * class will NOT cause their BaseField::delta values to be shared.
   *
   * @see fields::AddForce for repeated combining of forces.
   */
  template <
    typename _F0,
    typename _F1,
    typename _F2 = detail::NullForce<typename _F0::super0::options, 2u>,
    typename _F3 = detail::NullForce<typename _F0::super0::options, 3u>,
    typename _F4 = detail::NullForce<typename _F0::super0::options, 4u>,
    typename _F5 = detail::NullForce<typename _F0::super0::options, 5u>,
    typename _F6 = detail::NullForce<typename _F0::super0::options, 6u>,
    typename _F7 = detail::NullForce<typename _F0::super0::options, 7u>,
    typename _F8 = detail::NullForce<typename _F0::super0::options, 8u>,
    typename _F9 = detail::NullForce<typename _F0::super0::options, 9u>
  >
  class AddForces: public virtual BaseForce<typename _F0::super0::options>,
                   public _F0, public _F1, public _F2, public _F3, public _F4,
                   public _F5, public _F6, public _F7, public _F8, public _F9 {
    /* TYPEDEFS */
  public:
    typedef BaseForce<
      typename fields::detail
      ::assert< typename _F0::super0::options >
        ::template same< typename _F1::super0::options >::AND
        ::template same< typename _F2::super0::options >::AND
        ::template same< typename _F3::super0::options >::AND
        ::template same< typename _F4::super0::options >::AND
        ::template same< typename _F5::super0::options >::AND
        ::template same< typename _F6::super0::options >::AND
        ::template same< typename _F7::super0::options >::AND
        ::template same< typename _F8::super0::options >::AND
        ::template same< typename _F9::super0::options >::value
    > super0;
    typedef _F0 F0;
    typedef _F1 F1;
    typedef _F2 F2;
    typedef _F3 F3;
    typedef _F4 F4;
    typedef _F5 F5;
    typedef _F6 F6;
    typedef _F7 F7;
    typedef _F8 F8;
    typedef _F9 F9;


    /* MEMBER FUNCTIONS */
  public:
    AddForces( const typename super0::options::ChimpDB * db = NULL )
      : super0(db),
        F0(), F1(), F2(), F3(), F4(), F5(), F6(), F7(), F8(), F9() {}

    const AddForces & operator=(const AddForces & that) {
      super0::operator=(that);
      F0::operator=(that);
      F1::operator=(that);
      F2::operator=(that);
      F3::operator=(that);
      F4::operator=(that);
      F5::operator=(that);
      F6::operator=(that);
      F7::operator=(that);
      F8::operator=(that);
      F9::operator=(that);
      return *this;
    }

    void accel(       Vector<double,3> & a,
                const Vector<double,3> & r,
                const Vector<double,3> & v = V3(0.,0.,0.),
                const double & t = 0.0,
                const double & dt = 0.0,
                const unsigned int & species = 0u ) const {
      F0::accel(a,r,v,t,dt,species);
      detail::addF<F1>().accel(*this,a,r,v,t,dt,species);
      detail::addF<F2>().accel(*this,a,r,v,t,dt,species);
      detail::addF<F3>().accel(*this,a,r,v,t,dt,species);
      detail::addF<F4>().accel(*this,a,r,v,t,dt,species);
      detail::addF<F5>().accel(*this,a,r,v,t,dt,species);
      detail::addF<F6>().accel(*this,a,r,v,t,dt,species);
      detail::addF<F7>().accel(*this,a,r,v,t,dt,species);
      detail::addF<F8>().accel(*this,a,r,v,t,dt,species);
      detail::addF<F9>().accel(*this,a,r,v,t,dt,species);
    }

    template < typename P >
    void accel(       Vector<double,3> & a,
                const Vector<double,3> & r,
                const Vector<double,3> & v,
                const double & t,
                const double & dt,
                      P & p ) const {
      F0::accel(a,r,v,t,dt,p);
      detail::addF<F1>().accel(*this,a,r,v,t,dt,p);
      detail::addF<F2>().accel(*this,a,r,v,t,dt,p);
      detail::addF<F3>().accel(*this,a,r,v,t,dt,p);
      detail::addF<F4>().accel(*this,a,r,v,t,dt,p);
      detail::addF<F5>().accel(*this,a,r,v,t,dt,p);
      detail::addF<F6>().accel(*this,a,r,v,t,dt,p);
      detail::addF<F7>().accel(*this,a,r,v,t,dt,p);
      detail::addF<F8>().accel(*this,a,r,v,t,dt,p);
      detail::addF<F9>().accel(*this,a,r,v,t,dt,p);
    }

    double potential( const Vector<double,3> & r,
                      const Vector<double,3> & v = V3(0.,0.,0.),
                      const double & t = 0.0,
                      const unsigned int & species = 0u ) const {
      double retval = F0::potential(r,v,t,species);
      detail::addF<F1>().potential(*this,retval,r,v,t,species);
      detail::addF<F2>().potential(*this,retval,r,v,t,species);
      detail::addF<F3>().potential(*this,retval,r,v,t,species);
      detail::addF<F4>().potential(*this,retval,r,v,t,species);
      detail::addF<F5>().potential(*this,retval,r,v,t,species);
      detail::addF<F6>().potential(*this,retval,r,v,t,species);
      detail::addF<F7>().potential(*this,retval,r,v,t,species);
      detail::addF<F8>().potential(*this,retval,r,v,t,species);
      detail::addF<F9>().potential(*this,retval,r,v,t,species);
      return retval;
    }

    template < typename P >
    double potential( const Vector<double,3> & r,
                      const Vector<double,3> & v,
                      const double & t,
                            P & p ) const {
      double retval = F0::potential(r,v,t,p);
      detail::addF<F1>().potential(*this,retval,r,v,t,p);
      detail::addF<F2>().potential(*this,retval,r,v,t,p);
      detail::addF<F3>().potential(*this,retval,r,v,t,p);
      detail::addF<F4>().potential(*this,retval,r,v,t,p);
      detail::addF<F5>().potential(*this,retval,r,v,t,p);
      detail::addF<F6>().potential(*this,retval,r,v,t,p);
      detail::addF<F7>().potential(*this,retval,r,v,t,p);
      detail::addF<F8>().potential(*this,retval,r,v,t,p);
      detail::addF<F9>().potential(*this,retval,r,v,t,p);
      return retval;
    }

    template < unsigned int ndim,
               typename Particle >
    void applyStatisticalForce(       Vector<double,ndim> & xv,
                                const double & t,
                                const double & dt,
                                      Particle & particle ) const {
      F0::applyStatisticalForce( xv, t, dt, particle );
      F1::applyStatisticalForce( xv, t, dt, particle );
      F2::applyStatisticalForce( xv, t, dt, particle );
      F3::applyStatisticalForce( xv, t, dt, particle );
      F4::applyStatisticalForce( xv, t, dt, particle );
      F5::applyStatisticalForce( xv, t, dt, particle );
      F6::applyStatisticalForce( xv, t, dt, particle );
      F7::applyStatisticalForce( xv, t, dt, particle );
      F8::applyStatisticalForce( xv, t, dt, particle );
      F9::applyStatisticalForce( xv, t, dt, particle );
    }
  };


  /** Adds two Forces to get a total acceleration/ potential energy.
   * This metafunction can be used over and over again to aggregate forces.
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
   * @see fields::AddForces for aggregating a larger set of forces.
   */
  template < class _F0, class _F1 >
  class AddForce : public virtual BaseForce<typename _F0::super0::options>,
                   public _F0, public _F1 {
  public:
    typedef BaseForce<
      typename fields::detail
      ::assert< typename _F0::super0::options >
        ::template same< typename _F1::super0::options >::value
    > super0;
    typedef _F0 F0;
    typedef _F1 F1;

    AddForce( const typename super0::options::ChimpDB * db = NULL )
      : super0(db), F0(), F1() {}

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
