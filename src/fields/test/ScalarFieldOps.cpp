#define BOOST_TEST_MODULE  ScalarFieldOps

#include <fields/Fields.h>
#include <fields/indices.h>

#include <xylose/Vector.h>

#include <cmath>

#include <boost/test/unit_test.hpp>

namespace {

  using xylose::Vector;
  using xylose::V3;
  using namespace fields::indices;

  typedef fields::BgField<double> BgField;

  struct SFunctor {
    double operator() ( const Vector<double,3u> & r ) const {
      return r[X]*r[Y]*r[Z];
    }
  };

  typedef fields::FieldFunctor<SFunctor> ScalarField;

}/* namespace (anon) */

BOOST_AUTO_TEST_CASE( gradient ) {
  {
    BOOST_MESSAGE( "testing gradient for fixed scalar field" );
    BgField bgF( 42.0 );
    Vector<double,3u> gMag(-1.);
    BOOST_CHECK_EQUAL(
      fields::gradient(gMag, bgF, V3(1.,2.,3.) ),
      V3(0.,0.,0.)
    );
  }

  {
    BOOST_MESSAGE( "testing gradient for variable vector field" );
    ScalarField F;
    Vector<double,3u> gMag(-1.);
    fields::gradient(gMag, F, V3(1.,2.,3.) );
    BOOST_CHECK_CLOSE( gMag[X], 6., 1e-3 );
    BOOST_CHECK_CLOSE( gMag[Y], 3., 1e-3 );
    BOOST_CHECK_CLOSE( gMag[Z], 2., 1e-3 );
  }
}
