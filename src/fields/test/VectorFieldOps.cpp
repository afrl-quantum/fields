#define BOOST_TEST_MODULE  VectorFieldOps

#include <fields/Fields.h>
#include <fields/indices.h>

#include <xylose/Vector.h>

#include <cmath>

#include <boost/test/unit_test.hpp>

namespace {

  using xylose::Vector;
  using xylose::V3;
  using namespace fields::indices;

  typedef fields::BgField< Vector<double,3u> > BgField;

  struct VFunctor {
    void operator() (       Vector<double,3u> & F,
                      const Vector<double,3u> & r ) const {
      F[X] = r[X];
      F[Y] = r[Y]*r[Y];
      F[Z] = r[X]*r[Y]*r[Z];
    }
  };

  typedef fields::FieldFunctor<VFunctor> VectorField;

}/* namespace (anon) */

BOOST_AUTO_TEST_CASE( magnitude_of ) {
  {
    BOOST_MESSAGE( "testing magnitude_of for fixed vector field" );
    BgField bgF(V3(1.,2.,3.));
    fields::magnitude_of< BgField > mag( bgF );
    BOOST_CHECK_EQUAL( mag( V3(0.,0.,0.) ), std::sqrt(1.+ 4. + 9.) );
  }

  {
    BOOST_MESSAGE( "testing magnitude_of for variable vector field" );
    VectorField F;
    fields::magnitude_of< VectorField > mag( F );
    BOOST_CHECK_EQUAL( mag( V3(1.,2.,3.) ), std::sqrt(1.+ 16. + 36.) );
  }
}

BOOST_AUTO_TEST_CASE( gradient_of_magnitude ) {
  {
    BOOST_MESSAGE( "testing gradient_of_magnitude for fixed vector field" );
    BgField bgF(V3(1.,2.,3.));
    Vector<double,3u> gMag(-1.);
    BOOST_CHECK_EQUAL(
      fields::gradient_of_magnitude(gMag, bgF, V3(1.,2.,3.) ),
      V3(0.,0.,0.)
    );
  }

  {
    BOOST_MESSAGE( "testing gradient_of_magnitude for variable vector field" );
    VectorField F;
    Vector<double,3u> gMag(-1.);
    fields::gradient_of_magnitude(gMag, F, V3(1.,2.,3.) );
    BOOST_CHECK_CLOSE( gMag[X], 5.08234, 1e-3 );
    BOOST_CHECK_CLOSE( gMag[Y], 4.67026, 1e-3 );
    BOOST_CHECK_CLOSE( gMag[Z], 1.64833, 1e-3 );
  }
}
