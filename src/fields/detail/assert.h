
#ifndef fields_detail_assert_h
#define fields_detail_assert_h

namespace fields {
  namespace detail {

    template < typename T0 >
    struct assert {
      template < typename T1 >
      struct same {
        typedef T0 value;
        typedef assert<T0> AND;
      };
    };

  }/* namespace fields::detail */
}/* namespace fields */

#endif // fields_detail_assert_h
