// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALBASIS_HH

#include <array>
#include <bitset>
#include <numeric>
#include <vector>

#include <dune/common/fmatrix.hh>

#include "../../common/localbasis.hh"

namespace Dune
{
  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief First order Brezzi-Douglas-Marini shape functions on the reference
   *        tetrahedron.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   */
  template<class D, class R>
  class BDM1Simplex3DLocalBasis
  {

  public:
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,
        R,3,Dune::FieldVector<R,3>,
        Dune::FieldMatrix<R,3,3> > Traits;

    //! \brief Standard constructor
    BDM1Simplex3DLocalBasis()
    {
    }

    /**
     * \brief Make set number s, where 0 <= s < 64
     *
     * \param s Edge orientation indicator
     */
    BDM1Simplex3DLocalBasis(std::bitset<6> s)
    {
    }

    //! \brief number of shape functions
    unsigned int size() const
    {
      return 12;
    }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in Position
     * \param out return value
     */
    inline void evaluateFunction(const typename Traits::DomainType& in,
                                 std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      typedef typename Traits::DomainFieldType Field;

      const Field& x = in[0];
      const Field& y = in[1];
      const Field& z = in[2];

      out[  0 ][ 0 ] = 6*x;
      out[  0 ][ 1 ] = 6*y;
      out[  0 ][ 2 ] = -24*x - 24*y - 18*z + 18;

      out[  1 ][ 0 ] = -18*x;
      out[  1 ][ 1 ] = 6*y;
      out[  1 ][ 2 ] = 24*x + 6*z - 6;

      out[  2 ][ 0 ] = 6*x;
      out[  2 ][ 1 ] = -18*y;
      out[  2 ][ 2 ] = 24*y + 6*z - 6;

      out[  3 ][ 0 ] = -6*x;
      out[  3 ][ 1 ] = 24*x + 18*y + 24*z - 18;
      out[  3 ][ 2 ] = -6*z;

      out[  4 ][ 0 ] = 18*x;
      out[  4 ][ 1 ] = -24*x - 6*y + 6;
      out[  4 ][ 2 ] = -6*z;

      out[  5 ][ 0 ] = -6*x;
      out[  5 ][ 1 ] = -6*y - 24*z + 6;
      out[  5 ][ 2 ] = 18*z;

      out[  6 ][ 0 ] = -18*x - 24*y - 24*z + 18;
      out[  6 ][ 1 ] = 6*y;
      out[  6 ][ 2 ] = 6*z;

      out[  7 ][ 0 ] = 6*x + 24*y - 6;
      out[  7 ][ 1 ] = -18*y;
      out[  7 ][ 2 ] = 6*z;

      out[  8 ][ 0 ] = 6*x + 24*z - 6;
      out[  8 ][ 1 ] = 6*y;
      out[  8 ][ 2 ] = -18*z;

      out[  9 ][ 0 ] = 18*x;
      out[  9 ][ 1 ] = -6*y;
      out[  9 ][ 2 ] = -6*z;

      out[ 10 ][ 0 ] = -6*x;
      out[ 10 ][ 1 ] = 18*y;
      out[ 10 ][ 2 ] = -6*z;

      out[ 11 ][ 0 ] = -6*x;
      out[ 11 ][ 1 ] = -6*y;
      out[ 11 ][ 2 ] = 18*z;

    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in Position
     * \param out return value
     */
    inline void evaluateJacobian(const typename Traits::DomainType& in,
                                 std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      out[  0 ][ 0 ] = { 6, 0, 0 };
      out[  0 ][ 1 ] = { 0, 6, 0 };
      out[  0 ][ 2 ] = { -24, -24, -18 };

      out[  1 ][ 0 ] = { -18, 0, 0 };
      out[  1 ][ 1 ] = { 0, 6, 0 };
      out[  1 ][ 2 ] = { 24, 0, 6 };

      out[  2 ][ 0 ] = { 6, 0, 0 };
      out[  2 ][ 1 ] = { 0, -18, 0 };
      out[  2 ][ 2 ] = { 0, 24, 6 };

      out[  3 ][ 0 ] = { -6, 0, 0 };
      out[  3 ][ 1 ] = { 24, 18, 24 };
      out[  3 ][ 2 ] = { 0, 0, -6 };

      out[  4 ][ 0 ] = { 18, 0, 0 };
      out[  4 ][ 1 ] = { -24, -6, 0 };
      out[  4 ][ 2 ] = { 0, 0, -6 };

      out[  5 ][ 0 ] = { -6, 0, 0 };
      out[  5 ][ 1 ] = { 0, -6, -24 };
      out[  5 ][ 2 ] = { 0, 0, 18 };

      out[  6 ][ 0 ] = { -18, -24, -24 };
      out[  6 ][ 1 ] = { 0, 6, 0 };
      out[  6 ][ 2 ] = { 0, 0, 6 };

      out[  7 ][ 0 ] = { 6, 24, 0 };
      out[  7 ][ 1 ] = { 0, -18, 0 };
      out[  7 ][ 2 ] = { 0, 0, 6 };

      out[  8 ][ 0 ] = { 6, 0, 24 };
      out[  8 ][ 1 ] = { 0, 6, 0 };
      out[  8 ][ 2 ] = { 0, 0, -18 };

      out[  9 ][ 0 ] = { 18, 0, 0 };
      out[  9 ][ 1 ] = { 0, -6, 0 };
      out[  9 ][ 2 ] = { 0, 0, -6 };

      out[ 10 ][ 0 ] = { -6, 0, 0 };
      out[ 10 ][ 1 ] = { 0, 18, 0 };
      out[ 10 ][ 2 ] = { 0, 0, -6 };

      out[ 11 ][ 0 ] = { -6, 0, 0 };
      out[ 11 ][ 1 ] = { 0, -6, 0 };
      out[ 11 ][ 2 ] = { 0, 0, 18 };

    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 3>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      }
      else if (totalOrder == 1)
      {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        switch (direction) {
        case 0:
          out[  0 ] = { 6, 0, 0 };
          out[  1 ] = { -18, 0, 0 };
          out[  2 ] = { 6, 0, 0 };
          out[  3 ] = { -6, 0, 0 };
          out[  4 ] = { 18, 0, 0 };
          out[  5 ] = { -6, 0, 0 };
          out[  6 ] = { -18, -24, 0 };
          out[  7 ] = { 6, 24, -24 };
          out[  8 ] = { 6, 0, 0 };
          out[  9 ] = { 18, 0, 24 };
          out[ 10 ] = { -6, 0, 0 };
          out[ 11 ] = { -6, 0, 0 };
          break;

        case 1:
          out[  0 ] = { 0, 6, 0 };
          out[  1 ] = { 0, 6, 0 };
          out[  2 ] = { 0, -18, 0 };
          out[  3 ] = { 24, 18, 0 };
          out[  4 ] = { -24, -6, 24 };
          out[  5 ] = { 0, -6, 0 };
          out[  6 ] = { 0, 6, -24 };
          out[  7 ] = { 0, -18, 0 };
          out[  8 ] = { 0, 6, 0 };
          out[  9 ] = { 0, -6, 0 };
          out[ 10 ] = { 0, 18, 0 };
          out[ 11 ] = { 0, -6, 0 };
          break;

        case 2:
          out[  0 ] = { -24, -24, 18 };
          out[  1 ] = { 24, 0, -18 };
          out[  2 ] = { 0, 24, 6 };
          out[  3 ] = { 0, 0, 6 };
          out[  4 ] = { 0, 0, -6 };
          out[  5 ] = { 0, 0, -6 };
          out[  6 ] = { 0, 0, 18 };
          out[  7 ] = { 0, 0, 6 };
          out[  8 ] = { 0, 0, 6 };
          out[  9 ] = { 0, 0, -18 };
          out[ 10 ] = { 0, 0, -6 };
          out[ 11 ] = { 0, 0, -6 };
          break;


        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order() const
    {
      return 1;
    }
  };
} // end namespace Dune
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALBASIS_HH
