// -*- c++ -*-
// $Id$
/*@HEADER
 *         olson-tools:  A variety of routines and algorithms that
 *      I've developed and collected over the past few years.  This collection
 *      represents tools that are most useful for scientific and numerical
 *      software.  This software is released under the LGPL license except
 *      otherwise explicitly stated in individual files included in this
 *      package.  Generally, the files in this package are copyrighted by
 *      Spencer Olson--exceptions will be noted.   
 *                 Copyright 2004-2008 Spencer Olson
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *  
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *                                                                                 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA.                                                                           .
 * 
 * Questions? Contact Spencer Olson (olsonse@umich.edu) 
 */

/** \file
 * Force lookup specialization of FieldLookup class.
 * @see createFieldFile.h for routines to help creating the field-lookup table
 * file.
 */

/** \example field/lookup/testfield.cpp
 * \input field/lookup/common.h
 *
 * Demonstrates the field lookup code by using a set of field/force data from
 * a previously calculated field file.  The data pertains to magnetic forces
 * on an atom in the trapping hyperfine ground-state.
 */

/** \example field/lookup/createfieldfile.cpp
 * \input field/lookup/common.h
 *
 * Demonstrates the field field writing utility that produces a file that can
 * be read in by the field-lookup code.  The data pertains to magnetic forces
 * on an atom in the trapping hyperfine ground-state.
 */


#ifndef fields_force_lookup_h
#define fields_force_lookup_h

#include <fields/field-lookup.h>

#include <xylose/except.h>

#include <string>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace fields {

  using xylose::Vector;

  /* ************** BEGIN FORCE LOOKUP SPECIALIZATION *********** */

  template <unsigned int L = 3U, unsigned n_species = 1u>
  class ForceRecord;

  namespace detail { template < unsigned int > struct dummy; }
  template < unsigned int L >
  class ForceRecord<L,0u> {
    typedef typename
      detail::dummy<L>::ZERO_SPECIES_ForceRecord_NOT_SUPPORTED
      T;
  };

  template < unsigned int L >
  class ForceRecord<L,1u> {
  public:
    ForceRecord() : a(0.0), V(0.0) {}
    Vector<double,L> a;
    double V;

    /* This is an attempt to generalize the vector_lookup function but not
     * implement it as a macro (which would probably be faster though). */

    /** For using the FieldLookup::vector_lookup routine. */
    inline Vector<double,L> & vector(const unsigned int & i) { return a; }

    /** For using the FieldLookup::vector_lookup routine. */
    inline const Vector<double,L> & vector(const unsigned int & i) const { return a; }

    /** For using the FieldLookup::scalar_lookup routine. */
    inline double & scalar(const unsigned int & i) { return V; }

    /** For using the FieldLookup::scalar_lookup routine. */
    inline const double & scalar(const unsigned int & i) const { return V; }
  };

  /** Implementation of force record for multiple species.
   * Note that this does NOT do bounds checking on species number. */
  template < unsigned int L, unsigned int N >
  class ForceRecord {
  public:
    ForceRecord() {
      std::fill( a, a+N, 0.0 );
      std::fill( V, V+N, 0.0 );
    }
    Vector<double,L> a[N];
    double V[N];

    /* This is an attempt to generalize the vector_lookup function but not
     * implement it as a macro (which would probably be faster though). */

    /** For using the FieldLookup::vector_lookup routine. */
    inline Vector<double,L> & vector(const unsigned int & i) { return a[i]; }

    /** For using the FieldLookup::vector_lookup routine. */
    inline const Vector<double,L> & vector(const unsigned int & i) const {
      return a[i];
    }

    /** For using the FieldLookup::scalar_lookup routine. */
    inline double & scalar(const unsigned int & i) { return V[i]; }

    /** For using the FieldLookup::scalar_lookup routine. */
    inline const double & scalar(const unsigned int & i) const { return V[i]; }
  };


  namespace detail {
    template < unsigned int L, unsigned int N >
    struct ForceRecordIO {
      std::istream & in( std::istream & input, ForceRecord<L,N> & fr ) const {
        for ( unsigned int i = 0u; i < N; ++i )
          input >> fr.a[i] >> fr.V[i];
        return input;
      }
      std::ostream & out( std::ostream & output,
                          const ForceRecord<L,N> & fr ) const {
        const char * presep = "";
        for ( unsigned int i = 0u; i < N; ++i ) {
          output << presep << fr.a[i] << '\t' << fr.V[i];
          presep = "\t";
        }
        return output;
      }
    };

    template < unsigned int L >
    struct ForceRecordIO<L,1u> {
      std::istream & in( std::istream & input, ForceRecord<L,1u> & fr ) const {
        return input >> fr.a >> fr.V;
      }
      std::ostream & out( std::ostream & output,
                          const ForceRecord<L,1u> & fr ) const {
        return output << fr.a << '\t' << fr.V;
      }
    };
  }/* namespace detail */

  template <unsigned int L, unsigned int N>
  inline std::istream & operator>>( std::istream & input,
                                    ForceRecord<L,N> & fr) {
    return detail::ForceRecordIO<L,N>().in( input, fr );
  }

  template <unsigned int L, unsigned int N>
  inline std::ostream & operator<<(std::ostream & output, const ForceRecord<L,N> & fr) {
    return detail::ForceRecordIO<L,N>().out( output, fr );
  }

  /** The base Lookup class of the associated FieldLookup container class.
   * If you have an axially symmetric force, you can more optimally use
   * AxiSymFieldLookup< ForceRecord<L> >.  This will both speed up the lookup
   * process as well as minimize the memory footprint.
   *
   * Currently, this lookup table elemlent only supports multi-species,
   * single-velocity, single-time data.
   *
   * @see FieldLookup for cartesian lookup table.
   * @see AxiSymFieldLookup for axially symmetric lookup table.
   */
  template <unsigned int L = 3U, class T = FieldLookup< ForceRecord<L> > >
  class ForceLookup : public T {
  public:
    typedef T super;
    inline void accel( Vector<double,3> & a,
                       const Vector<double,3> & r,
                       const Vector<double,3> & v = V3(0.,0.,0.),
                       const double & t = 0.0,
                       const double & dt = 0.0,
                       const unsigned int & species = 0u ) const {
      super::vector_lookup(a, r, species);
    }

    template < typename P >
    inline void accel( Vector<double,3> & a,
                       const Vector<double,3> & r,
                       const Vector<double,3> & v,
                       const double & t,
                       const double & dt,
                       P & p ) const {
      super::vector_lookup(a, r, species(p) );
    }

    inline double potential( const Vector<double,3> & r,
                             const Vector<double,3> & v = V3(0.,0.,0.),
                             const double & t = 0.0,
                             const unsigned int & species = 0u ) const {
      return super::scalar_lookup(r, species);
    }

    template < typename P >
    inline double potential( const Vector<double,3> & r,
                             const Vector<double,3> & v,
                             const double & t,
                             const P & p ) const {
      return super::scalar_lookup(r, species(p));
    }
  }; /* ForceLookup class */

  template < class T, unsigned int L = 3U, unsigned int n_species = 1u >
  class ForceTableWrapper : public T {
  public:
    typedef T super;
    ForceRecord<L,n_species> getRecord(const Vector<double,3> & r) const {
      ForceRecord<L,n_species> retval;
      for ( unsigned int i = 0u; i < n_species; ++i ) {
        super::accel(retval.a[i], r, 0.0, 0.0, 0.0, i );
        retval.V[i] = super::potential(r, 0.0, 0.0, i );
      }
      return retval;
    }
  };

  template < class T, unsigned int L >
  class ForceTableWrapper<T,L,1u> : public T {
  public:
    typedef T super;
    ForceRecord<L,1u> getRecord(const Vector<double,3> & r) const {
      ForceRecord<L,1u> retval;
      super::accel(retval.a, r, 0.0, 0.0, 0.0, 0u );
      retval.V = super::potential(r, 0.0, 0.0, 0u );
      return retval;
    }
  };

  /* **************   END FORCE LOOKUP SPECIALIZATION *********** */

}/* namespace fields */

#endif // fields_force_lookup_h
