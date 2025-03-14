// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
#define DUNE_LAGRANGEBASIS_INTERPOLATION_HH

#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/typeutilities.hh>

#include <dune/localfunctions/lagrange/lagrangecoefficients.hh>

namespace Dune
{

  template< template <class,unsigned int> class LP,
      unsigned int dim, class F >
  struct LagrangeInterpolationFactory;

  // LocalLagrangeInterpolation
  // --------------------------

  template< template <class,unsigned int> class LP, unsigned int dim, class F >
  class LocalLagrangeInterpolation
  {
    typedef LocalLagrangeInterpolation< LP,dim,F > This;

  public:
    typedef LP<F,dim> LagrangePointSet;
    typedef typename LagrangePointSet::Field Field;

    static const unsigned int dimension = LagrangePointSet::dimension;

  private:
    friend struct LagrangeInterpolationFactory<LP,dim,F>;
    typedef typename LagrangePointSet::LagrangePoint::Vector Point;

    /*
     * The LagrangePointSet object is created inside the
     * LagrangeInterpolationFactory::create method and
     * the object is released in the corresponding
     * LagrangeInterpolationFactory::release method.
     */
    const LagrangePointSet &lagrangePoints_;

    // The object is constructed by the LagrangeInterpolationFactory
    explicit LocalLagrangeInterpolation ( const LagrangePointSet &lagrangePoints )
      : lagrangePoints_( lagrangePoints )
    {}

    // Access to the LagrangePointSet provided to the LagrangeInterpolationFactory
    const LagrangePointSet *points () const { return &lagrangePoints_; }

  public:
    /**
     * \brief Compute the Lagrange interpolation of the function `fn` and store
     * the interpolation coefficients associated to the Lagrange basis in the
     * output vector `coefficients`.
     *
     * The interpolation simply evaluates the function `fn` in all Lagrange points
     * given by the `LagrangePointSet`.
     */
    template< class Fn, class Vector,
      decltype(std::declval<Fn>()(std::declval<Point>()),bool{}) = true,
      decltype(std::declval<Vector>().size(),bool{}) = true,
      decltype(std::declval<Vector>().resize(0u),bool{}) = true>
    void interpolate ( const Fn &fn, Vector &coefficients ) const
    {
      coefficients.resize( lagrangePoints_.size() );

      unsigned int index = 0;
      for( const auto &lp : lagrangePoints_ )
        field_cast( fn( lp.point() ), coefficients[index++] );
    }

    /**
     * \brief Evaluate the (monomial) `basis` functions in the Lagrange points
     * to build up a Vandermonde matrix for computing the Lagrange basis functions.
     *
     * The Vandermonde matrix is stored in `coefficients`. It is assumed that the
     * matrix is represented as a vector-of-vectors like data structure. The rows
     * of the matrix correspond to the Lagrange point and the columns to the
     * (monomial) basis functions.
     */
    template< class Basis, class Matrix,
      decltype(std::declval<Matrix>().rows(),bool{}) = true,
      decltype(std::declval<Matrix>().cols(),bool{}) = true,
      decltype(std::declval<Matrix>().resize(0u,0u),bool{}) = true>
    void interpolate ( const Basis &basis, Matrix &coefficients ) const
    {
      coefficients.resize( lagrangePoints_.size(), basis.size( ) );

      unsigned int index = 0;
      for( const auto &lp : lagrangePoints_ )
        basis.template evaluate< 0 >( lp.point(), coefficients[index++] );
    }

    /// \brief Return the set of Lagrange points
    const LagrangePointSet &lagrangePoints () const { return lagrangePoints_; }
  };



  // LocalLagrangeInterpolationFactory
  // ---------------------------------
  template< template <class,unsigned int> class LP,
      unsigned int dim, class F >
  struct LagrangeInterpolationFactory
  {
    typedef LagrangeCoefficientsFactory<LP,dim,F> LagrangePointSetFactory;
    typedef typename LagrangePointSetFactory::Object LagrangePointSet;

    typedef typename LagrangePointSetFactory::Key Key;
    typedef const LocalLagrangeInterpolation< LP,dim,F > Object;

    template< GeometryType::Id geometryId >
    static Object *create ( const Key &key )
    {
      const LagrangePointSet *lagrangeCoeff
        = LagrangePointSetFactory::template create< geometryId >( key );
      if ( lagrangeCoeff == 0 )
        return 0;
      else
        return new Object( *lagrangeCoeff );
    }
    template< GeometryType::Id geometryId >
    static bool supports ( const Key &key )
    {
      return true;
    }
    static void release( Object *object)
    {
      LagrangePointSetFactory::release( object->points() );
      delete object;
    }
  };

}

#endif // #ifndef DUNE_LAGRANGEBASIS_INTERPOLATION_HH
