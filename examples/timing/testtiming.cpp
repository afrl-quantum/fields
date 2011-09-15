#include <fields/Fields.h>
#include <fields/Forces.h>
#include <fields/ScaleField.h>
#include <fields/ScaleForce.h>

#include <xylose/Vector.h>
#include <xylose/timing/Timing.h>
#include <xylose/timing/Printer.h>
#include <xylose/timing/element/Exponential.h>

#include <physical/physical.h>

#include <fstream>

//#include <cfloat>


namespace {
  using fields::ScaleField;
  using fields::ScaleForce;
  using fields::BgField;
  using fields::Gravity;

  namespace timing = xylose::timing;
  using xylose::Vector;
  using xylose::V3;

  using namespace physical::units;

  typedef ScaleField< BgField<double> > ScaledScalarField;
  typedef ScaleField< BgField< Vector<double,3> > > ScaledVectorField;
  typedef ScaleForce< Gravity > ScaledGravity;
}/* namespace anon */


int main() {
  double t_max = 15.*ms;
  double dt    = 0.1*ms;

  ScaledGravity gravity;
  ScaledScalarField sfield;
  ScaledVectorField vfield;

  timing::TimingsVector gtimings;
  gtimings.push_back(new timing::element::Exponential(3.*ms, 1.0, 0.0, 0.0));
  gtimings.push_back(new timing::element::Exponential(3.*ms, 1.0, 0.0, 1.0));
  gtimings.push_back(new timing::element::Exponential(3.*ms, 1.0, 0.0, 1.0));
  gtimings.push_back(new timing::element::Exponential(3.*ms, 3.0, 1.0, 0.0));
  gtimings.push_back(new timing::element::Exponential(3.*ms, 1.0, 0.0, 1.0));
  gravity.timing.timings = gtimings;

  sfield.bg = 1.0;
  timing::TimingsVector stimings;
  stimings.push_back(new timing::element::Exponential(2.*ms, 1.0, 0.0, 0.0));
  stimings.push_back(new timing::element::Exponential(4.*ms, 1.0, 0.0, 5.0));
  stimings.push_back(new timing::element::Exponential(1.*ms, 1.0, 0.0, 8.0));
  stimings.push_back(new timing::element::Exponential(5.*ms, 3.0, 3.0, 0.0));
  stimings.push_back(new timing::element::Exponential(3.*ms, 1.0, 0.0, 1.0));
  sfield.timing.timings = stimings;

  vfield.bg = V3(1.0,-1.0,-.5);
  timing::TimingsVector vtimings;
  vtimings.push_back(new timing::element::Exponential(4.*ms, 1.0, 0.0, 0.0));
  vtimings.push_back(new timing::element::Exponential(2.*ms, 1.0, 0.0, 3.0));
  vtimings.push_back(new timing::element::Exponential(5.*ms, 1.0, 0.0, 7.0));
  vtimings.push_back(new timing::element::Exponential(2.*ms, 3.0, 7.0, 0.0));
  vtimings.push_back(new timing::element::Exponential(2.*ms, 1.0, 0.0, 1.0));
  vfield.timing.timings = vtimings;


  std::ofstream gout("gravity-force.dat");
  std::ofstream sout("scalar-field.dat");
  std::ofstream vout("vector-field.dat");
  Vector<double,3> r(0.0);

  for (double t = 0.0; t <= t_max ; t+=dt) {
    gravity.timing.set_time(t);
    sfield.timing.set_time(t);
    vfield.timing.set_time(t);
    Vector<double,3> a;
    gravity.accel(a,r);

    gout << t << '\t'
         << a << '\n'
         << std::flush;

    sout << t << '\t'
         << sfield(r) << '\n'
         << std::flush;

    vfield(a,r);
    vout << t << '\t'
         << a << '\n'
         << std::flush;
  }

  timing::Printer tp;
  tp.timers.push_back(&vfield.timing);
  tp.timers.push_back(&sfield.timing);
  tp.timers.push_back(&gravity.timing);
  tp.print("timing.dat", 0.0, dt, t_max);
}

