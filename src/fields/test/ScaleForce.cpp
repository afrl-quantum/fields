#define BOOST_TEST_MODULE  ScaleForce

#include <fields/Forces.h>
#include <fields/indices.h>
#include <fields/ScaleForce.h>

#include <xylose/Vector.h>
#include <xylose/timing/Timing.h>
#include <xylose/timing/element/PowerLaw.h>

#include <physical/physical.h>

#include <sstream>
#include <cmath>

#include <boost/test/unit_test.hpp>

namespace {
  typedef fields::Gravity<> Gravity;
  using fields::ScaleForce;
  using namespace fields::indices;

  namespace timing = xylose::timing;
  using xylose::Vector;
  using xylose::V3;

  using namespace physical::units;

  typedef ScaleForce< Gravity > ScaledGravity;
  typedef ScaledGravity::options::ChimpDB ChimpDB;
}/* namespace (anon) */

BOOST_AUTO_TEST_CASE( scale_vectorforce ) {
  ChimpDB db;
  db.addParticleType("87Rb");
  db.initBinaryInteractions();

  ScaledGravity gravity;
  gravity.db = &db;

  timing::TimingsVector timings;
  timings.push_back(new timing::element::PowerLaw(1.5*ms, 1.0, 0.0, 0.0));
  timings.push_back(new timing::element::PowerLaw(3.0*ms, 1.0, 0.0, 5.0));
  timings.push_back(new timing::element::PowerLaw(1.0*ms, 1.0, 0.0, 8.0));
  timings.push_back(new timing::element::PowerLaw(3.0*ms, 3.0, 3.0, 0.0));
  timings.push_back(new timing::element::PowerLaw(2.0*ms, 1.0, 0.0, 1.0));
  gravity.timing.timings = timings;

  const char * ans =
    "0	0	0	-0	0\n"
    "0.001	0	0	-0	0\n"
    "0.002	0	0	-8.175	1.18102e-23\n"
    "0.003	0	0	-24.525	3.54305e-23\n"
    "0.004	0	0	-40.875	5.90508e-23\n"
    "0.005	0	0	-39.24	5.66888e-23\n"
    "0.006	0	0	-29.2938	4.23198e-23\n"
    "0.007	0	0	-25.7513	3.7202e-23\n"
    "0.008	0	0	-12.3988	1.79121e-23\n"
    "0.009	0	0	-2.4525	3.54305e-24\n";

  std::ostringstream out;
  const Vector<double,3u> r(10.0);
  for ( double t = 0.0; t <= 10*ms; t += ms ) {
    gravity.timing.set_time(t);
    Vector<double,3u> a;
    gravity.accel(a, r);
    out << t << '\t' << a << '\t' << gravity.potential(r) << '\n';
  }
  BOOST_CHECK_EQUAL( out.str(), ans );
}

