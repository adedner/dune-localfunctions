// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGEPOINTSET_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGEPOINTSET_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /** \brief A single Lagrange point
   *
   * \tparam F Number type
   * \tparam dim Dimension of the domain
   */
  template< class F, unsigned int dim >
  class LagrangePoint
  {
    typedef LagrangePoint< F, dim > This;

    template< class, class >
    friend class LagrangePointSetImpl;

  public:
    static const int dimension = dim;

    typedef F Field;

    /** \brief The vector type used for the Lagrange point position */
    typedef FieldVector< Field, dimension > Vector;

    /** \brief Position of the Lagrange point */
    const Vector &point () const
    {
      return point_;
    }

    /** \brief Assignment to a (sub-)entity */
    const LocalKey &localKey () const
    {
      return localKey_;
    }

    const Field weight () const
    {
      return weight_;
    }

    Vector point_ = {};
    LocalKey localKey_ = {};
    Field weight_ = {};
  };

  /** \brief A set of Lagrange points on an unspecified domain
   *
   * \tparam F The number type used for weights and coordinates
   * \tparam dim The number of the domain
   */
  template< class F, unsigned int dim >
  class LagrangePointSet
  {
    typedef EmptyPointSet< F, dim > This;

  public:
    typedef F Field;

    static const unsigned int dimension = dim;

    typedef Dune::LagrangePoint< Field, dimension > LagrangePoint;

    typedef typename std::vector< LagrangePoint >::const_iterator iterator;

  protected:
    LagrangePointSet ( const std::size_t order )
      : order_( order )
    {}

  public:
    const LagrangePoint &operator[] ( const unsigned int i ) const
    {
      assert( i < size() );
      return points_[ i ];
    }

    iterator begin () const
    {
      return points_.begin();
    }

    iterator end () const
    {
      return points_.end();
    }

    const LocalKey &localKey ( const unsigned int i ) const
    {
      return (*this)[ i ].localKey();
    }

    std::size_t order () const
    {
      return order_;
    }

    std::size_t size () const
    {
      return points_.size();
    }

  protected:
    std::size_t order_;
    std::vector< LagrangePoint > points_;
  };

}

#endif // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGEPOINTSET_HH
