
#ifndef fields_detail_NullForce_h
#define fields_detail_NullForce_h

#include <xylose/Vector.h>

namespace fields {

  namespace detail {
    using xylose::Vector;
    using xylose::V3;

    template < typename _options, unsigned int >
    struct NullForce {
      struct super0 {
        typedef _options options;
      };

      template < unsigned int ndim, typename Particle >
      void applyStatisticalForce(       Vector<double,ndim> & xv,
                                  const double & t,
                                  const double & dt,
                                        Particle & particle ) const { }
    };

    template < typename F >
    struct addF {
      void accel( const F & f,
                        Vector<double,3> & a,
                  const Vector<double,3> & r,
                  const Vector<double,3> & v,
                  const double & t,
                  const double & dt,
                  const unsigned int & species ) const {
        Vector<double,3u> a2;
        f.accel(a2,r,v,t,dt,species);
        a += a2;
      }

      template < typename P >
      void accel( const F & f,
                        Vector<double,3> & a,
                  const Vector<double,3> & r,
                  const Vector<double,3> & v,
                  const double & t,
                  const double & dt,
                        P & p ) const {
        Vector<double,3u> a2;
        f.accel(a2,r,v,t,dt,p);
        a += a2;
      }

      void potential( const F & f,
                            double & retval,
                      const Vector<double,3> & r,
                      const Vector<double,3> & v,
                      const double & t,
                      const unsigned int & species ) const {
        retval += f.potential(r,v,t,species);
      }

      template < typename P >
      void potential( const F & f,
                            double & retval,
                      const Vector<double,3> & r,
                      const Vector<double,3> & v,
                      const double & t,
                            P & p ) const {
        retval += f.potential(r,v,t,p);
      }
    };



    template < typename O, unsigned int u >
    struct addF< NullForce<O,u> > {
      void accel( const NullForce<O,u> & f,
                        Vector<double,3> & a,
                  const Vector<double,3> & r,
                  const Vector<double,3> & v,
                  const double & t,
                  const double & dt,
                  const unsigned int & species ) const { }

      template < typename P >
      void accel( const NullForce<O,u> & f,
                        Vector<double,3> & a,
                  const Vector<double,3> & r,
                  const Vector<double,3> & v,
                  const double & t,
                  const double & dt,
                        P & p ) const { }

      void potential( const NullForce<O,u> & f,
                            double & retval,
                      const Vector<double,3> & r,
                      const Vector<double,3> & v,
                      const double & t,
                      const unsigned int & species ) const { }

      template < typename P >
      void potential( const NullForce<O,u> & f,
                            double & retval,
                      const Vector<double,3> & r,
                      const Vector<double,3> & v,
                      const double & t,
                            P & p ) const { }
    };


  }/* namespace fields::detail */
}/* namespace fields */

#endif // fields_detail_NullForce_h
