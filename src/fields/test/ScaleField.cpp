#define BOOST_TEST_MODULE  ScaleField

#include <fields/Fields.h>
#include <fields/indices.h>
#include <fields/ScaleField.h>

#include <xylose/Vector.h>
#include <xylose/timing/Timing.h>
#include <xylose/timing/element/PowerLaw.h>

#include <physical/physical.h>

#include <sstream>
#include <cmath>

#include <boost/test/unit_test.hpp>

namespace {

  using fields::BgField;
  using fields::ScaleField;
  using namespace fields::indices;

  using xylose::Vector;
  using xylose::V3;
  namespace timing = xylose::timing;
  using xylose::timing::element::Base;
  using xylose::timing::element::PowerLaw;

  using namespace physical::units;

  typedef ScaleField< BgField<double> > ScaledScalarField;
  typedef ScaleField< BgField< Vector<double,3> > > ScaledVectorField;

}/* namespace (anon) */

BOOST_AUTO_TEST_CASE( scale_scalarfield ) {
  ScaledScalarField field;
  field.bg = 1.0;

  timing::TimingsVector timings;
  timings.push_back(new timing::element::PowerLaw(1.5*ms, 1.0, 0.0, 0.0));
  timings.push_back(new timing::element::PowerLaw(3.0*ms, 1.0, 0.0, 5.0));
  timings.push_back(new timing::element::PowerLaw(1.0*ms, 1.0, 0.0, 8.0));
  timings.push_back(new timing::element::PowerLaw(3.0*ms, 3.0, 3.0, 0.0));
  timings.push_back(new timing::element::PowerLaw(2.0*ms, 1.0, 0.0, 1.0));
  field.timing.timings = timings;

  const char * ans = 
    "0\t0\n"
    "0.001\t0\n"
    "0.002\t0.833333\n"
    "0.003\t2.5\n"
    "0.004\t4.16667\n"
    "0.005\t4\n"
    "0.006\t2.98611\n"
    "0.007\t2.625\n"
    "0.008\t1.26389\n"
    "0.009\t0.25\n";

  std::ostringstream out;
  const Vector<double,3u> r(0.0);
  for ( double t = 0.0; t <= 10*ms; t += ms ) {
    field.timing.set_time(t);
    out << t << '\t' << field(r) << '\n';
  }
  BOOST_CHECK_EQUAL( out.str(), ans );
}

BOOST_AUTO_TEST_CASE( scale_vectorfield ) {
  ScaledVectorField field;
  field.bg = V3(1.0,-1.0,-.5);

  timing::TimingsVector timings;
  timings.push_back(new timing::element::PowerLaw(1.5*ms, 1.0, 0.0, 0.0));
  timings.push_back(new timing::element::PowerLaw(3.0*ms, 1.0, 0.0, 5.0));
  timings.push_back(new timing::element::PowerLaw(1.0*ms, 1.0, 0.0, 8.0));
  timings.push_back(new timing::element::PowerLaw(3.0*ms, 3.0, 3.0, 0.0));
  timings.push_back(new timing::element::PowerLaw(2.0*ms, 1.0, 0.0, 1.0));
  field.timing.timings = timings;

  const char * ans =
    "0\t0\t-0\t-0\n"
    "0.001\t0\t-0\t-0\n"
    "0.002\t0.833333\t-0.833333\t-0.416667\n"
    "0.003\t2.5\t-2.5\t-1.25\n"
    "0.004\t4.16667\t-4.16667\t-2.08333\n"
    "0.005\t4\t-4\t-2\n"
    "0.006\t2.98611\t-2.98611\t-1.49306\n"
    "0.007\t2.625\t-2.625\t-1.3125\n"
    "0.008\t1.26389\t-1.26389\t-0.631944\n"
    "0.009\t0.25\t-0.25\t-0.125\n";

  std::ostringstream out;
  const Vector<double,3u> r(0.0);
  for ( double t = 0.0; t <= 10*ms; t += ms ) {
    field.timing.set_time(t);
    Vector<double,3u> F;
    field(F,r);
    out << t << '\t' << F << '\n';
  }
  BOOST_CHECK_EQUAL( out.str(), ans );
}
