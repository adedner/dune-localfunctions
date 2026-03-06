// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_CUBE_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_CUBE_LOCALBASIS_HH

#include <algorithm>
#include <array>
#include <bitset>
#include <numeric>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  template<class D, class R, unsigned int dim, unsigned int order>
  class BDMCubeLocalBasis
  {
    static_assert( AlwaysFalse<D>::value,
                   "`BDMCubeLocalBasis` not implemented for chosen `dim` and `order`." );
  };


#ifndef DOXYGEN
  template<class D, class R, unsigned int dim>
  class BDMCubeLocalBasis<D, R, dim, 0>
  {
    static_assert( AlwaysFalse<D>::value,
                   "`BDMCubeLocalBasis` not defined for order 0." );
  };
#endif // #ifndef DOXYGEN


  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief First order Brezzi-Douglas-Marini shape functions on the reference quadrilateral.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   */
  template<class D, class R >
  class BDMCubeLocalBasis<D, R, 2, 1>
  {
    using DomainType    = FieldVector<D, 2>;
    using RangeType     = FieldVector<R, 2>;
    using JacobianType  = FieldMatrix<R, 2, 2>;

  public:
    using Traits = LocalBasisTraits<D, 2, DomainType, R, 2, RangeType, JacobianType>;

    //! \brief Standard constructor
    BDMCubeLocalBasis ()
    {
      for (size_t i=0; i<4; i++)
        sign_[i] = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    BDMCubeLocalBasis (std::bitset<4> s)
    {
      for (size_t i=0; i<4; i++)
        sign_[i] = s[i] ? -1.0 : 1.0;
    }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 8;
    }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in Position
     * \param out return value
     */
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(8);

      out[0][0] = sign_[0]*(in[0] - 1.0);
      out[0][1] = 0.0;
      out[1][0] = 6.0*in[0]*in[1] - 3.0*in[0]-6*in[1] + 3.0;
      out[1][1] = -3.0*in[1]*in[1] + 3.0*in[1];
      out[2][0] = sign_[1]*(in[0]);
      out[2][1] = 0.0;
      out[3][0] = -6.0*in[0]*in[1] + 3.0*in[0];
      out[3][1] = 3.0*in[1]*in[1] - 3.0*in[1];
      out[4][0] = 0.0;
      out[4][1] = sign_[2]*(in[1] - 1.0);
      out[5][0] = 3.0*in[0]*in[0] - 3.0*in[0];
      out[5][1] = -6.0*in[0]*in[1] + 6.0*in[0] + 3.0*in[1] - 3.0;
      out[6][0] = 0.0;
      out[6][1] = sign_[3]*(in[1]);
      out[7][0] = -3.0*in[0]*in[0] + 3.0*in[0];
      out[7][1] = 6.0*in[0]*in[1] - 3.0*in[1];
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in Position
     * \param out return value
     */
    inline void evaluateJacobian (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(8);

      out[0][0][0] = sign_[0];
      out[0][0][1] = 0.0;
      out[0][1][0] = 0.0;
      out[0][1][1] = 0.0;

      out[1][0][0] = 6.0*in[1] - 3.0;
      out[1][0][1] = 6.0*in[0] - 6.0;
      out[1][1][0] = 0.0;
      out[1][1][1] = -6.0*in[1] + 3.0;

      out[2][0][0] = sign_[1];
      out[2][0][1] = 0.0;
      out[2][1][0] = 0.0;
      out[2][1][1] = 0.0;

      out[3][0][0] = -6.0*in[1] + 3.0;
      out[3][0][1] = -6.0*in[0];
      out[3][1][0] = 0.0;
      out[3][1][1] = 6.0*in[1] - 3.0;

      out[4][0][0] = 0.0;
      out[4][0][1] = 0.0;
      out[4][1][0] = 0.0;
      out[4][1][1] = sign_[2];

      out[5][0][0] = 6.0*in[0] - 3.0;
      out[5][0][1] = 0.0;
      out[5][1][0] = -6.0*in[1] + 6.0;
      out[5][1][1] = -6.0*in[0] + 3.0;

      out[6][0][0] = 0.0;
      out[6][0][1] = 0.0;
      out[6][1][0] = 0.0;
      out[6][1][1] = sign_[3];

      out[7][0][0] = -6.0*in[0] + 3.0;
      out[7][0][1] = 0.0;
      out[7][1][0] = 6.0*in[1];
      out[7][1][1] = 6.0*in[0] - 3.0;
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 2>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        switch (direction) {
        case 0:
          out[0][0] = sign_[0];
          out[0][1] = 0.0;

          out[1][0] = 6.0*in[1] - 3.0;
          out[1][1] = 0.0;

          out[2][0] = sign_[1];
          out[2][1] = 0.0;

          out[3][0] = -6.0*in[1] + 3.0;
          out[3][1] = 0.0;

          out[4][0] = 0.0;
          out[4][1] = 0.0;

          out[5][0] = 6.0*in[0] - 3.0;
          out[5][1] = -6.0*in[1] + 6.0;

          out[6][0] = 0.0;
          out[6][1] = 0.0;

          out[7][0] = -6.0*in[0] + 3.0;
          out[7][1] = 6.0*in[1];
          break;
        case 1:
          out[0][0] = 0.0;
          out[0][1] = 0.0;

          out[1][0] = 6.0*in[0] - 6.0;
          out[1][1] = -6.0*in[1] + 3.0;

          out[2][0] = 0.0;
          out[2][1] = 0.0;

          out[3][0] = -6.0*in[0];
          out[3][1] = 6.0*in[1] - 3.0;

          out[4][0] = 0.0;
          out[4][1] = sign_[2];

          out[5][0] = 0.0;
          out[5][1] = -6.0*in[0] + 3.0;

          out[6][0] = 0.0;
          out[6][1] = sign_[3];

          out[7][0] = 0.0;
          out[7][1] = 6.0*in[0] - 3.0;
          break;
        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return 2;
    }

  private:
    std::array<R,4> sign_;
  };

  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief Second order Brezzi-Douglas-Marini shape functions on quadrilaterals.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   */
  template<class D, class R >
  class BDMCubeLocalBasis<D, R, 2, 2>
  {
    using DomainType    = FieldVector<D, 2>;
    using RangeType     = FieldVector<R, 2>;
    using JacobianType  = FieldMatrix<R, 2, 2>;

  public:
    using Traits = LocalBasisTraits<D, 2, DomainType, R, 2, RangeType, JacobianType>;

    //! \brief Standard constructor
    BDMCubeLocalBasis()
    {
      for (size_t i=0; i<4; i++)
        sign_[i] = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    BDMCubeLocalBasis(std::bitset<4> s)
    {
      for (size_t i=0; i<4; i++)
        sign_[i] = s[i] ? -1.0 : 1.0;
    }

    //! \brief number of shape functions
    unsigned int size() const
    {
      return 14;
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

      out[0][0] = sign_[0]*(-2.25 + 5.25*in[0] + 7.5*in[1] - 7.5*in[0]*in[1] - 3.0*in[0]*in[0] - 7.5*in[1]*in[1] + 7.5*in[0]*in[1]*in[1]);
      out[0][1] = sign_[0]*(-1.25*in[1] + 3.75*in[1]*in[1] - 2.5*in[1]*in[1]*in[1]);
      out[1][0] = 3.0 - 3.0*in[0]-6.0*in[1] + 6.0*in[0]*in[1];
      out[1][1] = 0.0;
      out[2][0] = sign_[0]*(-3.75 + 3.75*in[0] + 22.5*in[1] - 22.5*in[0]*in[1] - 22.5*in[1]*in[1] + 22.5*in[0]*in[1]*in[1]);
      out[2][1] = sign_[0]*(-3.75*in[1] + 11.25*in[1]*in[1] - 7.5*in[1]*in[1]*in[1]);
      out[3][0] = sign_[1]*(-0.75*in[0] - 7.5*in[0]*in[1] + 3.0*in[0]*in[0] + 7.5*in[0]*in[1]*in[1]);
      out[3][1] = sign_[1]*(-1.25*in[1] + 3.75*in[1]*in[1] - 2.5*in[1]*in[1]*in[1]);
      out[4][0] = 3.0*in[0] - 6.0*in[0]*in[1];
      out[4][1] = 0.0;
      out[5][0] = sign_[1]*(+3.75*in[0] - 22.5*in[0]*in[1] + 22.5*in[0]*in[1]*in[1]);
      out[5][1] = sign_[1]*(-3.75*in[1] + 11.25*in[1]*in[1] - 7.5*in[1]*in[1]*in[1]);
      out[6][0] = sign_[2]*(-1.25*in[0] + 3.75*in[0]*in[0] - 2.5*in[0]*in[0]*in[0]);
      out[6][1] = sign_[2]*(-2.25 + 7.5*in[0] + 5.25*in[1] - 7.5*in[0]*in[1] - 7.5*in[0]*in[0] - 3.0*in[1]*in[1] + 7.5*in[0]*in[0]*in[1]);
      out[7][0] = 0.0;
      out[7][1] = -3.0 + 6.0*in[0] + 3.0*in[1] - 6.0*in[0]*in[1];
      out[8][0] = sign_[2]*(-3.75*in[0] + 11.25*in[0]*in[0] - 7.5*in[0]*in[0]*in[0]);
      out[8][1] = sign_[2]*(-3.75 + 22.5*in[0] + 3.75*in[1] - 22.5*in[0]*in[1] - 22.5*in[0]*in[0] + 22.5*in[0]*in[0]*in[1]);
      out[9][0] = sign_[3]*(-1.25*in[0] + 3.75*in[0]*in[0] - 2.5*in[0]*in[0]*in[0]);
      out[9][1] = sign_[3]*(-0.75*in[1] - 7.5*in[0]*in[1] + 3.0*in[1]*in[1] + 7.5*in[0]*in[0]*in[1]);
      out[10][0] = 0.0;
      out[10][1] = -3.0*in[1] + 6.0*in[0]*in[1];
      out[11][0] = sign_[3]*(-3.75*in[0] + 11.25*in[0]*in[0] - 7.5*in[0]*in[0]*in[0]);
      out[11][1] = sign_[3]*(3.75*in[1] - 22.5*in[0]*in[1] + 22.5*in[0]*in[0]*in[1]);
      out[12][0] = 6.0*in[0] - 6.0*in[0]*in[0];
      out[12][1] = 0.0;
      out[13][0] = 0.0;
      out[13][1] = 6.0*in[1] - 6.0*in[1]*in[1];
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

      out[0][0][0] = sign_[0]*(5.25 - 7.5*in[1] - 6.0*in[0] + 7.5*in[1]*in[1]);
      out[0][0][1] = sign_[0]*(7.5 - 7.5*in[0] - 15.0*in[1] + 15.0*in[0]*in[1]);
      out[0][1][0] = 0.0;
      out[0][1][1] = sign_[0]*(-1.25 + 7.5*in[1] - 7.5*in[1]*in[1]);

      out[1][0][0] = -3.0 + 6.0*in[1];
      out[1][0][1] = -6.0 + 6.0*in[0];
      out[1][1][0] = 0.0;
      out[1][1][1] = 0.0;

      out[2][0][0] = sign_[0]*(3.75 - 22.5*in[1] + 22.5*in[1]*in[1]);
      out[2][0][1] = sign_[0]*(22.5 - 22.5*in[0] - 45.0*in[1] + 45.0*in[0]*in[1]);
      out[2][1][0] = 0.0;
      out[2][1][1] = sign_[0]*(-3.75 + 22.5*in[1] - 22.5*in[1]*in[1]);

      out[3][0][0] = sign_[1]*(-0.75 - 7.5*in[1] + 6.0*in[0] + 7.5*in[1]*in[1]);
      out[3][0][1] = sign_[1]*(-7.5*in[0] + 15.0*in[0]*in[1]);
      out[3][1][0] = 0.0;
      out[3][1][1] = sign_[1]*(-1.25 + 7.5*in[1] - 7.5*in[1]*in[1]);

      out[4][0][0] = 3.0 - 6.0*in[1];
      out[4][0][1] = -6.0*in[0];
      out[4][1][0] = 0.0;
      out[4][1][1] = 0.0;

      out[5][0][0] = sign_[1]*(3.75 - 22.5*in[1] + 22.5*in[1]*in[1]);
      out[5][0][1] = sign_[1]*(-22.5*in[0] + 45.0*in[0]*in[1]);
      out[5][1][0] = 0.0;
      out[5][1][1] = sign_[1]*(-3.75 + 22.5*in[1] - 22.5*in[1]*in[1]);

      out[6][0][0] = sign_[2]*(-1.25 + 7.5*in[0] - 7.5*in[0]*in[0]);
      out[6][0][1] = 0.0;
      out[6][1][0] = sign_[2]*(7.5 - 7.5*in[1] - 15.0*in[0] + 15.0*in[0]*in[1]);
      out[6][1][1] = sign_[2]*(5.25 - 7.5*in[0]- 6.0*in[1] + 7.5*in[0]*in[0]);

      out[7][0][0] = 0.0;
      out[7][0][1] = 0.0;
      out[7][1][0] = 6.0 - 6.0*in[1];
      out[7][1][1] = 3.0 - 6.0*in[0];

      out[8][0][0] = sign_[2]*(-3.75 + 22.5*in[0] - 22.5*in[0]*in[0]);
      out[8][0][1] = 0.0;
      out[8][1][0] = sign_[2]*(22.5 - 22.5*in[1] - 45.0*in[0] + 45.0*in[0]*in[1]);
      out[8][1][1] = sign_[2]*(3.75 - 22.5*in[0] + 22.5*in[0]*in[0]);

      out[9][0][0] = sign_[3]*(-1.25 + 7.5*in[0] - 7.5*in[0]*in[0]);
      out[9][0][1] = 0.0;
      out[9][1][0] = sign_[3]*(-7.5*in[1] + 15.0*in[0]*in[1]);
      out[9][1][1] = sign_[3]*(-0.75 - 7.5*in[0] + 6.0*in[1] + 7.5*in[0]*in[0]);

      out[10][0][0] = 0.0;
      out[10][0][1] = 0.0;
      out[10][1][0] = 6.0*in[1];
      out[10][1][1] = -3.0 + 6.0*in[0];

      out[11][0][0] = sign_[3]*(-3.75 + 22.5*in[0] - 22.5*in[0]*in[0]);
      out[11][0][1] = 0.0;
      out[11][1][0] = sign_[3]*(-22.5*in[1] + 45*in[0]*in[1]);
      out[11][1][1] = sign_[3]*(3.75 - 22.5*in[0] + 22.5*in[0]*in[0]);

      out[12][0][0] = 6.0 - 12.0*in[0];
      out[12][0][1] = 0.0;
      out[12][1][0] = 0.0;
      out[12][1][1] = 0.0;

      out[13][0][0] = 0.0;
      out[13][0][1] = 0.0;
      out[13][1][0] = 0.0;
      out[13][1][1] = 6.0 - 12.0*in[1];
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 2>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        switch (direction) {
        case 0:
          out[0][0] = sign_[0]*(5.25 - 7.5*in[1] - 6.0*in[0] + 7.5*in[1]*in[1]);
          out[0][1] = 0.0;

          out[1][0] = -3.0 + 6.0*in[1];
          out[1][1] = 0.0;

          out[2][0] = sign_[0]*(3.75 - 22.5*in[1] + 22.5*in[1]*in[1]);
          out[2][1] = 0.0;

          out[3][0] = sign_[1]*(-0.75 - 7.5*in[1] + 6.0*in[0] + 7.5*in[1]*in[1]);
          out[3][1] = 0.0;

          out[4][0] = 3.0 - 6.0*in[1];
          out[4][1] = 0.0;

          out[5][0] = sign_[1]*(3.75 - 22.5*in[1] + 22.5*in[1]*in[1]);
          out[5][1] = 0.0;

          out[6][0] = sign_[2]*(-1.25 + 7.5*in[0] - 7.5*in[0]*in[0]);
          out[6][1] = sign_[2]*(7.5 - 7.5*in[1] - 15.0*in[0] + 15.0*in[0]*in[1]);

          out[7][0] = 0.0;
          out[7][1] = 6.0 - 6.0*in[1];

          out[8][0] = sign_[2]*(-3.75 + 22.5*in[0] - 22.5*in[0]*in[0]);
          out[8][1] = sign_[2]*(22.5 - 22.5*in[1] - 45.0*in[0] + 45.0*in[0]*in[1]);

          out[9][0] = sign_[3]*(-1.25 + 7.5*in[0] - 7.5*in[0]*in[0]);
          out[9][1] = sign_[3]*(-7.5*in[1] + 15.0*in[0]*in[1]);

          out[10][0] = 0.0;
          out[10][1] = 6.0*in[1];

          out[11][0] = sign_[3]*(-3.75 + 22.5*in[0] - 22.5*in[0]*in[0]);
          out[11][1] = sign_[3]*(-22.5*in[1] + 45*in[0]*in[1]);

          out[12][0] = 6.0 - 12.0*in[0];
          out[12][1] = 0.0;

          out[13][0] = 0.0;
          out[13][1] = 0.0;
          break;
        case 1:
          out[0][0] = sign_[0]*(7.5 - 7.5*in[0] - 15.0*in[1] + 15.0*in[0]*in[1]);
          out[0][1] = sign_[0]*(-1.25 + 7.5*in[1] - 7.5*in[1]*in[1]);

          out[1][0] = -6.0 + 6.0*in[0];
          out[1][1] = 0.0;

          out[2][0] = sign_[0]*(22.5 - 22.5*in[0] - 45.0*in[1] + 45.0*in[0]*in[1]);
          out[2][1] = sign_[0]*(-3.75 + 22.5*in[1] - 22.5*in[1]*in[1]);

          out[3][0] = sign_[1]*(-7.5*in[0] + 15.0*in[0]*in[1]);
          out[3][1] = sign_[1]*(-1.25 + 7.5*in[1] - 7.5*in[1]*in[1]);

          out[4][0] = -6.0*in[0];
          out[4][1] = 0.0;

          out[5][0] = sign_[1]*(-22.5*in[0] + 45.0*in[0]*in[1]);
          out[5][1] = sign_[1]*(-3.75 + 22.5*in[1] - 22.5*in[1]*in[1]);

          out[6][0] = 0.0;
          out[6][1] = sign_[2]*(5.25 - 7.5*in[0]- 6.0*in[1] + 7.5*in[0]*in[0]);

          out[7][0] = 0.0;
          out[7][1] = 3.0 - 6.0*in[0];

          out[8][0] = 0.0;
          out[8][1] = sign_[2]*(3.75 - 22.5*in[0] + 22.5*in[0]*in[0]);

          out[9][0] = 0.0;
          out[9][1] = sign_[3]*(-0.75 - 7.5*in[0] + 6.0*in[1] + 7.5*in[0]*in[0]);

          out[10][0] = 0.0;
          out[10][1] = -3.0 + 6.0*in[0];

          out[11][0] = 0.0;
          out[11][1] = sign_[3]*(3.75 - 22.5*in[0] + 22.5*in[0]*in[0]);

          out[12][0] = 0.0;
          out[12][1] = 0.0;

          out[13][0] = 0.0;
          out[13][1] = 6.0 - 12.0*in[1];
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
      return 3;
    }

  private:
    std::array<R,4> sign_;
  };


  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief First order Brezzi-Douglas-Marini shape functions on hexahedrons.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   */
  template<class D, class R >
  class BDMCubeLocalBasis<D, R, 3, 1>
  {
    using DomainType    = FieldVector<D, 3>;
    using RangeType     = FieldVector<R, 3>;
    using JacobianType  = FieldMatrix<R, 3, 3>;

  public:
    using Traits = LocalBasisTraits<D, 3, DomainType, R, 3, RangeType, JacobianType>;

    //! \brief Standard constructor
    BDMCubeLocalBasis()
    {
      for (size_t i=0; i<6; i++)
        sign_[i] = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 64
     *
     * \param s Edge orientation indicator
     */
    BDMCubeLocalBasis(std::bitset<6> s)
    {
      for (size_t i=0; i<6; i++)
        sign_[i] = s[i] ? -1.0 : 1.0;
    }

    //! \brief number of shape functions
    unsigned int size() const
    {
      return 18;
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

      out[0][0]  = sign_[0] * (in[0] - 1.0);
      out[0][1]  = 0;
      out[0][2]  = 0;
      out[1][0]  = sign_[1] * in[0];
      out[1][1]  = 0;
      out[1][2]  = 0;
      out[2][0]  = 0;
      out[2][1]  = sign_[2] * (in[1] - 1.0);
      out[2][2]  = 0;
      out[3][0]  = 0;
      out[3][1]  = sign_[3] * in[1];
      out[3][2]  = 0;
      out[4][0]  = 0;
      out[4][1]  = 0;
      out[4][2]  = sign_[4] * (in[2] - 1.0);
      out[5][0]  = 0;
      out[5][1]  = 0;
      out[5][2]  = sign_[5] * in[2];
      out[6][0]  =  6.0 * in[0] * in[1] - 3 * in[0]-6 * in[1] + 3.0;
      out[6][1]  = -3.0 * in[1] * in[1] + 3 * in[1];
      out[6][2]  =  0;
      out[7][0]  = -6.0 * in[0] * in[1] + 3 * in[0];
      out[7][1]  =  3.0 * in[1] * in[1] - 3 * in[1];
      out[7][2]  =  0;
      out[8][0]  =  3.0 * in[0] * in[0] - 3 * in[0];
      out[8][1]  = -6.0 * in[0] * in[1] + 3 * in[1]+6 * in[0]-3.0;
      out[8][2]  =  0;
      out[9][0]  = -3.0 * in[0] * in[0] + 3 * in[0];
      out[9][1]  =  6.0 * in[0] * in[1] - 3 * in[1];
      out[9][2]  =  0;
      out[10][0] = -3.0 * in[0] * in[0] + 3 * in[0];
      out[10][1] =  0;
      out[10][2] =  6.0 * in[0] * in[2]-6 * in[0]-3 * in[2] + 3.0;
      out[11][0] =  3.0 * in[0] * in[0]-3 * in[0];
      out[11][1] =  0;
      out[11][2] = -6.0 * in[0] * in[2] + 3 * in[2];
      out[12][0] = -6.0 * in[0] * in[2]+6 * in[2] + 3 * in[0]-3.0;
      out[12][1] =  0;
      out[12][2] =  3.0 * in[2] * in[2]-3 * in[2];
      out[13][0] = -3 * in[0]+6 * in[0] * in[2];
      out[13][1] =  0;
      out[13][2] = -3.0 * in[2] * in[2] + 3 * in[2];
      out[14][0] =  0;
      out[14][1] =  6.0 * in[1] * in[2]-3 * in[1]-6 * in[2] + 3.0;
      out[14][2] = -3 * in[2] * in[2] + 3 * in[2];
      out[15][0] =  0;
      out[15][1] = -6.0 * in[1] * in[2] + 3 * in[1];
      out[15][2] =  3.0 * in[2] * in[2]-3 * in[2];
      out[16][0] =  0;
      out[16][1] =  3.0 * in[1] * in[1]-3 * in[1];
      out[16][2] = -6.0 * in[1] * in[2] + 3 * in[2]+6 * in[1]-3.0;
      out[17][0] =  0;
      out[17][1] = -3.0 * in[1] * in[1] + 3 * in[1];
      out[17][2] =  6.0 * in[1] * in[2] - 3.0 * in[2];
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

      out[0][0]  = {  sign_[0],          0,          0};
      out[0][1]  = {         0,          0,          0};
      out[0][2]  = {         0,          0,          0};

      out[1][0]  = {  sign_[1],          0,          0};
      out[1][1]  = {         0,          0,          0};
      out[1][2]  = {         0,          0,          0};

      out[2][0]  = {         0,          0,          0};
      out[2][1]  = {         0,   sign_[2],          0};
      out[2][2]  = {         0,          0,          0};

      out[3][0]  = {         0,          0,          0};
      out[3][1]  = {         0,   sign_[3],          0};
      out[3][2]  = {         0,          0,          0};

      out[4][0]  = {         0,          0,          0};
      out[4][1]  = {         0,          0,          0};
      out[4][2]  = {         0,          0,   sign_[4]};

      out[5][0]  = {         0,          0,          0};
      out[5][1]  = {         0,          0,          0};
      out[5][2]  = {         0,          0,   sign_[5]};

      out[6][0]  = { 6*in[1]-3,  6*in[0]-6,          0};
      out[6][1]  = {         0, -6*in[1]+3,          0};
      out[6][2]  = {         0,          0,          0};

      out[7][0]  = {-6*in[1]+3,   -6*in[0],          0};
      out[7][1]  = {         0,  6*in[1]-3,          0};
      out[7][2]  = {         0,          0,          0};

      out[8][0]  = { 6*in[0]-3,          0,          0};
      out[8][1]  = {-6*in[1]+6, -6*in[0]+3,          0};
      out[8][2]  = {         0,          0,          0};

      out[9][0]  = {-6*in[0]+3,          0,          0};
      out[9][1]  = {   6*in[1],  6*in[0]-3,          0};
      out[9][2]  = {         0,          0,          0};

      out[10][0] = {-6*in[0]+3,          0,          0};
      out[10][1] = {         0,          0,          0};
      out[10][2] = { 6*in[2]-6,          0,  6*in[0]-3};

      out[11][0] = { 6*in[0]-3,          0,          0};
      out[11][1] = {         0,          0,          0};
      out[11][2] = {  -6*in[2],          0, -6*in[0]+3};

      out[12][0] = {-6*in[2]+3,          0, -6*in[0]+6};
      out[12][1] = {         0,          0,          0};
      out[12][2] = {         0,          0,  6*in[2]-3};

      out[13][0] = { 6*in[2]-3,          0,    6*in[0]};
      out[13][1] = {         0,          0,          0};
      out[13][2] = {         0,          0, -6*in[2]+3};

      out[14][0] = {         0,          0,          0};
      out[14][1] = {         0,  6*in[2]-3,  6*in[1]-6};
      out[14][2] = {         0,          0, -6*in[2]+3};

      out[15][0] = {         0,          0,          0};
      out[15][1] = {         0, -6*in[2]+3,   -6*in[1]};
      out[15][2] = {         0,          0,  6*in[2]-3};

      out[16][0] = {         0,          0,          0};
      out[16][1] = {         0,  6*in[1]-3,          0};
      out[16][2] = {         0, -6*in[2]+6, -6*in[1]+3};

      out[17][0] = {         0,          0,          0};
      out[17][1] = {         0, -6*in[1]+3,          0};
      out[17][2] = {         0,    6*in[2],  6*in[1]-3};
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 3>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        switch (direction) {
        case 0:
          out[0]  = {  sign_[0],          0,          0};
          out[1]  = {  sign_[1],          0,          0};
          out[2]  = {         0,          0,          0};
          out[3]  = {         0,          0,          0};
          out[4]  = {         0,          0,          0};
          out[5]  = {         0,          0,          0};
          out[6]  = { 6*in[1]-3,          0,          0};
          out[7]  = {-6*in[1]+3,          0,          0};
          out[8]  = { 6*in[0]-3, -6*in[1]+6,          0};
          out[9]  = {-6*in[0]+3,    6*in[1],          0};
          out[10] = {-6*in[0]+3,          0,  6*in[2]-6};
          out[11] = { 6*in[0]-3,          0,   -6*in[2]};
          out[12] = {-6*in[2]+3,          0,          0};
          out[13] = { 6*in[2]-3,          0,          0};
          out[14] = {         0,          0,          0};
          out[15] = {         0,          0,          0};
          out[16] = {         0,          0,          0};
          out[17] = {         0,          0,          0};
          break;
        case 1:
          out[0]  = {         0,          0,          0};
          out[1]  = {         0,          0,          0};
          out[2]  = {         0,   sign_[2],          0};
          out[3]  = {         0,   sign_[3],          0};
          out[4]  = {         0,          0,          0};
          out[5]  = {         0,          0,          0};
          out[6]  = { 6*in[0]-6, -6*in[1]+3,          0};
          out[7]  = {  -6*in[0],  6*in[1]-3,          0};
          out[8]  = {         0, -6*in[0]+3,          0};
          out[9]  = {         0,  6*in[0]-3,          0};
          out[10] = {         0,          0,          0};
          out[11] = {         0,          0,          0};
          out[12] = {         0,          0,          0};
          out[13] = {         0,          0,          0};
          out[14] = {         0,  6*in[2]-3,          0};
          out[15] = {         0, -6*in[2]+3,          0};
          out[16] = {         0,  6*in[1]-3, -6*in[2]+6};
          out[17] = {         0, -6*in[1]+3,    6*in[2]};
          break;
        case 2:
          out[0]  = {         0,          0,          0};
          out[1]  = {         0,          0,          0};
          out[2]  = {         0,          0,          0};
          out[3]  = {         0,          0,          0};
          out[4]  = {         0,          0,   sign_[4]};
          out[5]  = {         0,          0,   sign_[5]};
          out[6]  = {         0,          0,          0};
          out[7]  = {         0,          0,          0};
          out[8]  = {         0,          0,          0};
          out[9]  = {         0,          0,          0};
          out[10] = {         0,          0,  6*in[0]-3};
          out[11] = {         0,          0, -6*in[0]+3};
          out[12] = {-6*in[0]+6,          0,  6*in[2]-3};
          out[13] = {   6*in[0],          0, -6*in[2]+3};
          out[14] = {         0,  6*in[1]-6, -6*in[2]+3};
          out[15] = {         0,   -6*in[1],  6*in[2]-3};
          out[16] = {         0,          0, -6*in[1]+3};
          out[17] = {         0,          0,  6*in[1]-3};
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
      return 2;
    }

  private:
    std::array<R,6> sign_;
  };


} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_CUBE_LOCALBASIS_HH
