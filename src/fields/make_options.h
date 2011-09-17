
#ifndef fields_make_options_h
#define fields_make_options_h

#include <chimp/RuntimeDB.h>
#include <chimp/make_options.h>

namespace fields {

  /** Metafunction to generate the fields::options class.
   * @param _chimp_options
   *  All options used for the Chimp::RuntimeDB class.
   *  Note that _chimp_options::setParticle< _Particle >::type will be used to
   *  define the actual ChimpDB interface to use.
   *  [Default chimp::make_options<>::type]
   * */
  template < typename _chimp_options = chimp::make_options<>::type >
  struct make_options {
    struct type {
      typedef _chimp_options                    chimp_options;
      typedef chimp::RuntimeDB< chimp_options > ChimpDB;

      /* **** Template metafunctions to change DSMC configuration. **** */

      /** Simple little sub-template class to make it easier to change the
       * extra data fields in each node. */
      template < typename T >
      struct changeChimpOptions {
        typedef typename make_options<
          T
        >::type type;
      };

    };
  };/* make_options */
}/* namespace fields */

#endif // fields_make_options_h
