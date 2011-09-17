
#ifndef fields_examples_field_lookup_common_h
#define fields_examples_field_lookup_common_h

#include <fields/bfield.h>
#include <fields/force-lookup.h>

#include <physical/physical.h>

namespace {

  using namespace physical::constants::si;
  using namespace physical::units;

  const double I = 300.0 * Ampere;
  const double delta_B = 0.1 * microns;

  namespace fbf = fields::BField;

  /** Guide wires without the knee. */
  fbf::ThinCurrentElement wires[] = {
    /* source for wire 2. */
    fbf::ThinCurrentElement(-0.0032567, 0.0, -2.0, -0.0012567, 0.0, 4.0,I),

    /* source for wire 1. */
    fbf::ThinCurrentElement( 0.0032567, 0.0, -2.0,  0.0012567, 0.0, 4.0,I),

    /* end marker. */
    fbf::ThinCurrentElement(0,0,0,0,0,0,0)
  };

  /** Function for easily adding ThiWireBSrc elements to a ThinWireSrc. */
  template <class ThinWireBSrc>
  inline void addwires(ThinWireBSrc & bsrc) {
    for (int i = 0; fabs(wires[i].I) > 0; bsrc.currents.push_back(wires[i++]));
  }

  typedef fields::Gravity<> Gravity;
  typedef fields::AddForce<
    fbf::BCalcs< fbf::ThinWireSrc >,
    Gravity
  > BFieldForce;

  typedef BFieldForce::options::ChimpDB ChimpDB;

  typedef fields::ForceTableWrapper< BFieldForce > BFieldForceTableSrc;

  #define FIELD_FILENAME "field.dat"

  #define ERR_FILE "error.dat"

}/* namespace (anon) */

#endif // fields_examples_field_lookup_common_h
