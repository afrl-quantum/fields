#define BOOST_TEST_MODULE  AddForce

#include <fields/Forces.h>
#include <fields/make_options.h>

#include <xylose/Vector.h>
#include <xylose/timing/Timing.h>
#include <xylose/timing/element/PowerLaw.h>

#include <physical/physical.h>

#include <sstream>
#include <cmath>

#include <boost/test/unit_test.hpp>

namespace {
  using fields::Gravity;
  using fields::AddForce;

  namespace timing = xylose::timing;
  using xylose::Vector;
  using xylose::V3;

  using namespace physical::units;

  template < typename options = fields::make_options<>::type >
  class Other : public virtual fields::BaseForce<options>,
                public fields::BgField< Vector<double,3> > {
  public:
    typedef fields::BaseForce<options> super0;
    typedef fields::BgField< Vector<double,3> > super1;

    Other() : super0(), super1() { bg = V3(5.,4.,0.); }

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
      return - super0::db[species].mass::value * (super1::bg * r);
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


  typedef AddForce< Gravity<>, Other<> > DualGravity;
}/* namespace (anon) */

BOOST_AUTO_TEST_CASE( added_gravities ) {
  DualGravity gravity;
}

