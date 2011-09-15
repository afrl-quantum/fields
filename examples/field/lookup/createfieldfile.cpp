#include "common.h"

#include <fields/createFieldFile.h>
#include <fields/indices.h>

#include <xylose/strutil.h>
#include <xylose/Vector.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>


namespace {
  using xylose::Vector;
  using xylose::V3;
  using namespace fields::indices;

  const Vector<double,3> X_MINc   = V3(-20.0*um,          -20.0*um,         -20.*um );
  const Vector<double,3> X_MAXc   = V3( 20.0*um  + 1e-12,  20.0*um  + 1e-12, 20.*um + 1e-12 );
  const Vector<double,3> X_MINs   = V3(-100.0*um,         -100.0*um,        -20.*um );
  const Vector<double,3> X_MAXs   = V3( 100.0*um + 1e-12,  100.0*um + 1e-12, 20.*um + 1e-12 );
  const Vector<double,3> dxc      = V3(5.*um,   5.*um,  5.*um);
  const Vector<double,3> dxs      = V3(20.*um, 20.*um, 20.*um);
}

int main() {
  BFieldForceTableSrc bsrc;
  addwires(bsrc);

  bsrc.mass = mass;
  bsrc.Gravity::bg[Z] = -physical::unit::gravity;
  bsrc.delta = delta_B;

  createFieldFile(bsrc,
                  X_MINc,
                  X_MAXc,
                  dxc,
                  X_MINs,
                  X_MAXs,
                  dxs,
                  FIELD_FILENAME);

  return 0;
}

