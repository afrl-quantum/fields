#define BOOST_TEST_MODULE  ScaleForce

#include <fields/Forces.h>
#include <fields/indices.h>
#include <fields/ScaleForce.h>

#include <xylose/Vector.h>
#include <xylose/timing/Timing.h>
#include <xylose/timing/element/Exponential.h>

#include <physical/physical.h>

#include <sstream>
#include <cmath>

#include <boost/test/unit_test.hpp>

namespace {
  using fields::Gravity;
  using fields::ScaleForce;
  using namespace fields::indices;

  namespace timing = xylose::timing;
  using xylose::Vector;
  using xylose::V3;

  using namespace physical::units;

  typedef ScaleForce< Gravity > ScaledGravity;
}/* namespace (anon) */

BOOST_AUTO_TEST_CASE( scale_vectorforce ) {
  ScaledGravity gravity;

  timing::TimingsVector timings;
  timings.push_back(new timing::element::Exponential(1.5*ms, 1.0, 0.0, 0.0));
  timings.push_back(new timing::element::Exponential(3.0*ms, 1.0, 0.0, 5.0));
  timings.push_back(new timing::element::Exponential(1.0*ms, 1.0, 0.0, 8.0));
  timings.push_back(new timing::element::Exponential(3.0*ms, 3.0, 3.0, 0.0));
  timings.push_back(new timing::element::Exponential(2.0*ms, 1.0, 0.0, 1.0));
  gravity.timing.timings = timings;

  const char * ans =  
    "0\t0\t0\t-0\n"
    "0.001\t0\t0\t-0\n"
    "0.002\t0\t0\t-8.175\n"
    "0.003\t0\t0\t-24.525\n"
    "0.004\t0\t0\t-40.875\n"
    "0.005\t0\t0\t-39.24\n"
    "0.006\t0\t0\t-29.2938\n"
    "0.007\t0\t0\t-25.7513\n"
    "0.008\t0\t0\t-12.3988\n"
    "0.009\t0\t0\t-2.4525\n";

  std::ostringstream out;
  const Vector<double,3u> r(0.0);
  for ( double t = 0.0; t <= 10*ms; t += ms ) {
    gravity.timing.set_time(t);
    Vector<double,3u> a;
    gravity.accel(a, r);
    out << t << '\t' << a << '\n';
  }
  BOOST_CHECK_EQUAL( out.str(), ans );
}

