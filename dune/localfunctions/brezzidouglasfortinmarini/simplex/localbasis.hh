// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_SIMPLEX_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_SIMPLEX_LOCALBASIS_HH

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

#ifdef BDFM_USE_RT0_BASIS
#include <dune/localfunctions/raviartthomas/raviartthomas02d/raviartthomas02dlocalbasis.hh>
#endif

namespace Dune
{
  /**
   * \ingroup BrezziDouglasFortinMariniImpl
   * \brief Brezzi-Douglas-Fortin-Marini shape functions on a reference simplex.
   *
   * \tparam D      Type to represent the field in the domain.
   * \tparam R      Type to represent the field in the range.
   * \tparam dim    dimension of the reference element, must be >= 2.
   * \tparam order  order of the element, must be >= 1.
   *
   * \nosubgrouping
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMSimplexLocalBasis
  {
    static_assert( AlwaysFalse<D>::value,
                   "`BDFMSimplexLocalBasis` not implemented for chosen `dim` and `order`." );
  };


  /**
   * \brief Zero-th order Brezzi-Douglas-Fortin-Marini shape functions on the reference triangle.
   *
   * \note This coincides with the Raviart-Thomas element or order 0 for triangles.
   *
   * \nosubgrouping
   */
#ifdef BDFM_USE_RT0_BASIS
  template<class D, class R >
  class BDFMSimplexLocalBasis<D, R, 2, 0>
    : public RT02DLocalBasis< D, R >
  {
    using DomainType    = FieldVector<D, 2>;
    using RangeType     = FieldVector<R, 2>;
    using JacobianType  = FieldMatrix<R, 2, 2>;

  public:
    using Traits = LocalBasisTraits<D, 2, DomainType, R, 2, RangeType, JacobianType>;

    //! \brief Standard constructor
    BDFMSimplexLocalBasis () {}

    /**
     * \brief Make set number s, where 0<= s < 16
     *
     * \param s  Edge orientation indicator
     */
    BDFMSimplexLocalBasis (std::bitset<4> s) : RT02DLocalBasis< D, R >(1) {}

  };
#else
  template<class D, class R >
  class BDFMSimplexLocalBasis<D, R, 2, 0>
  {
    using DomainType    = FieldVector<D, 2>;
    using RangeType     = FieldVector<R, 2>;
    using JacobianType  = FieldMatrix<R, 2, 2>;

  public:
    using Traits = LocalBasisTraits<D, 2, DomainType, R, 2, RangeType, JacobianType>;

    //! \brief Standard constructor
    BDFMSimplexLocalBasis () {}

    /**
     * \brief Make set number s, where 0<= s < 16
     *
     * \param s  Edge orientation indicator
     */
    BDFMSimplexLocalBasis (std::bitset<3> s) :
    {
      for (auto i : range(3))
        s_[i] = s[i] ? -1 : 1;
    }

    //! \brief number of shape functions
    unsigned int size () const { return 3; }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateFunction (const DomainType& in, std::vector<RangeType>& out) const
    {
      out.resize(size());

      const auto& x = in[0];
      const auto& y = in[1];

      out[  0 ][ 0 ] = -x;
      out[  0 ][ 1 ] = 1 - y;

      out[  1 ][ 0 ] = x - 1;
      out[  1 ][ 1 ] = y;

      out[  2 ][ 0 ] = -x;
      out[  2 ][ 1 ] = -y;
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateJacobian (const DomainType& in, std::vector<JacobianType>& out) const
    {
      out.resize(size());

      out[  0 ][ 0 ] = { -1, 0 };
      out[  0 ][ 1 ] = { 0, -1 };

      out[  1 ][ 0 ] = { 1, 0 };
      out[  1 ][ 1 ] = { 0, 1 };

      out[  2 ][ 0 ] = { -1, 0 };
      out[  2 ][ 1 ] = { 0, -1 };
    }

    /**
     * \brief Evaluate all partial derivatives of all shape functions
     *
     * \param order  order the partial derivative
     * \param in     Position
     * \param out    return value
     */
    void partial (const std::array<unsigned int, 2>& order,
                  const DomainType& in,
                  std::vector<RangeType>& out) const
    {
      if (std::accumulate(order.begin(), order.end(), 0) == 0)
        evaluateFunction(in, out);
      else
      {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        switch (direction) {
        case 0:
          out[  0 ] = { -1, 0 };
          out[  1 ] = { 1, 0 };
          out[  2 ] = { -1, 0 };
          break;

        case 1:
          out[  0 ] = { 0, -1 };
          out[  1 ] = { 0, 1 };
          out[  2 ] = { 0, -1 };
          break;
        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const { return 0; }

  private:
    std::array<R, 3> s_;
  };
#endif


  /**
   * \brief First order Brezzi-Douglas-Fortin-Marini shape functions on the reference triangle.
   *
   * \nosubgrouping
   */
  template<class D, class R >
  class BDFMSimplexLocalBasis<D, R, 2, 1>
  {
    using DomainType    = FieldVector<D, 2>;
    using RangeType     = FieldVector<R, 2>;
    using JacobianType  = FieldMatrix<R, 2, 2>;

  public:
    using Traits = LocalBasisTraits<D, 2, DomainType, R, 2, RangeType, JacobianType>;

    //! \brief Standard constructor
    BDFMSimplexLocalBasis ()
    {
      std::fill(s_.begin(), s_.end(), 1);
    }

    /**
     * \brief Make set number s, where 0<= s < 16
     *
     * \param s  Edge orientation indicator
     */
    BDFMSimplexLocalBasis (std::bitset<4> s)
    {
      for (auto i : range(9))
        s_[i] = s[i] ? -1 : 1;
    }

    //! \brief number of shape functions
    unsigned int size () const { return 9; }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateFunction (const DomainType& in, std::vector<RangeType>& out) const
    {
      out.resize(size());

      const auto& x = in[0];
      const auto& y = in[1];

      out[  0 ][ 0 ] = 5*x*x - 2*x*(x + y) - x;
      out[  0 ][ 1 ] = 5*x*y - 6*x + 13*y*(x + y) - 17*y + 4;

      out[  1 ][ 0 ] = -13*x*x + 10*x*(x + y) - x;
      out[  1 ][ 1 ] = -13*x*y + 6*x - 5*y*(x + y) + 7*y - 2;

      out[  2 ][ 0 ] = 5*x*x - 18*x*(x + y) + 17*x + 6*y - 4;
      out[  2 ][ 1 ] = 5*x*y - 3*y*(x + y) + y;

      out[  3 ][ 0 ] = -13*x*x + 18*x*(x + y) - 7*x - 6*y + 2;
      out[  3 ][ 1 ] = -13*x*y + 3*y*(x + y) + y;

      out[  4 ][ 0 ] = -3*x*x - 10*x*(x + y) + 9*x;
      out[  4 ][ 1 ] = -3*x*y + 5*y*(x + y) - 3*y;

      out[  5 ][ 0 ] = 3*x*x + 2*x*(x + y) - 3*x;
      out[  5 ][ 1 ] = 3*x*y - 13*y*(x + y) + 9*y;

      out[  6 ][ 0 ] = 12*x*x - 48*x*(x + y) + 36*x;
      out[  6 ][ 1 ] = 12*x*y + 12*y*(x + y) - 12*y;

      out[  7 ][ 0 ] = -12*x*x + 24*x*(x + y) - 12*x;
      out[  7 ][ 1 ] = -12*x*y - 36*y*(x + y) + 36*y;

      out[  8 ][ 0 ] = 36*x*x - 48*x*(x + y) + 12*x;
      out[  8 ][ 1 ] = 36*x*y + 12*y*(x + y) - 12*y;
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateJacobian (const DomainType& in, std::vector<JacobianType>& out) const
    {
      out.resize(size());

      const auto& x = in[0];
      const auto& y = in[1];

      out[  0 ][ 0 ] = { 6*x - 2*y - 1, -2*x };
      out[  0 ][ 1 ] = { 18*y - 6, 18*x + 26*y - 17 };

      out[  1 ][ 0 ] = { -6*x + 10*y - 1, 10*x };
      out[  1 ][ 1 ] = { 6 - 18*y, -18*x - 10*y + 7 };

      out[  2 ][ 0 ] = { -26*x - 18*y + 17, 6 - 18*x };
      out[  2 ][ 1 ] = { 2*y, 2*x - 6*y + 1 };

      out[  3 ][ 0 ] = { 10*x + 18*y - 7, 18*x - 6 };
      out[  3 ][ 1 ] = { -10*y, -10*x + 6*y + 1 };

      out[  4 ][ 0 ] = { -26*x - 10*y + 9, -10*x };
      out[  4 ][ 1 ] = { 2*y, 2*x + 10*y - 3 };

      out[  5 ][ 0 ] = { 10*x + 2*y - 3, 2*x };
      out[  5 ][ 1 ] = { -10*y, -10*x - 26*y + 9 };

      out[  6 ][ 0 ] = { -72*x - 48*y + 36, -48*x };
      out[  6 ][ 1 ] = { 24*y, 24*x + 24*y - 12 };

      out[  7 ][ 0 ] = { 24*x + 24*y - 12, 24*x };
      out[  7 ][ 1 ] = { -48*y, -48*x - 72*y + 36 };

      out[  8 ][ 0 ] = { -24*x - 48*y + 12, -48*x };
      out[  8 ][ 1 ] = { 48*y, 48*x + 24*y - 12 };
    }

    /**
     * \brief Evaluate all partial derivatives of all shape functions
     *
     * \param order  order the partial derivative
     * \param in     Position
     * \param out    return value
     */
    void partial (const std::array<unsigned int, 2>& order,
                  const DomainType& in,
                  std::vector<RangeType>& out) const
    {
      if (std::accumulate(order.begin(), order.end(), 0) == 0)
        evaluateFunction(in, out);
      else
      {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        const auto& x = in[0];
        const auto& y = in[1];

        switch (direction) {
        case 0:
          out[  0 ] = { 6*x - 2*y - 1, -2*x };
          out[  1 ] = { -6*x + 10*y - 1, 10*x };
          out[  2 ] = { -26*x - 18*y + 17, 6 - 18*x };
          out[  3 ] = { 10*x + 18*y - 7, 18*x - 6 };
          out[  4 ] = { -26*x - 10*y + 9, -10*x };
          out[  5 ] = { 10*x + 2*y - 3, 2*x };
          out[  6 ] = { -72*x - 48*y + 36, -48*x };
          out[  7 ] = { 24*x + 24*y - 12, 24*x };
          out[  8 ] = { -24*x - 48*y + 12, -48*x };
          break;

        case 1:
          out[  0 ] = { 18*y - 6, 18*x + 26*y - 17 };
          out[  1 ] = { 6 - 18*y, -18*x - 10*y + 7 };
          out[  2 ] = { 2*y, 2*x - 6*y + 1 };
          out[  3 ] = { -10*y, -10*x + 6*y + 1 };
          out[  4 ] = { 2*y, 2*x + 10*y - 3 };
          out[  5 ] = { -10*y, -10*x - 26*y + 9 };
          out[  6 ] = { 24*y, 24*x + 24*y - 12 };
          out[  7 ] = { -48*y, -48*x - 72*y + 36 };
          out[  8 ] = { 48*y, 48*x + 24*y - 12 };
          break;
        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const { return 1; }

  private:
    std::array<R, 9> s_;
  };

  /**
   * \brief Zero-th order Brezzi-Douglas-Fortin-Marini shape functions on the reference tetrahedron.
   *
   * \nosubgrouping
   */
  template<class D, class R >
  class BDFMSimplexLocalBasis<D, R, 3, 0>
  {
    static constexpr int dim = 3;
    using DomainType    = FieldVector<D, dim>;
    using RangeType     = FieldVector<R, dim>;
    using JacobianType  = FieldMatrix<R, dim, dim>;

  public:
    using Traits = LocalBasisTraits<D, dim, DomainType, R, dim, RangeType, JacobianType>;

    //! \brief Standard constructor
    BDFMSimplexLocalBasis ()
    {
      //std::fill(s_.begin(), s_.end(), 1);
    }

    /**
     * \brief Make set number s, where 0<= s < 16
     *
     * \param s  Edge orientation indicator
     */
    BDFMSimplexLocalBasis (int s)
    {
      /*
      for (auto i : range(3))
        s_[i] = s[i] ? -1 : 1;
        */
    }

    //! \brief number of shape functions
    unsigned int size () const { return dim+1; }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateFunction (const DomainType& in, std::vector<RangeType>& out) const
    {
      out.resize(size());

      const auto& x = in[0];
      const auto& y = in[1];
      const auto& z = in[2];

      out[  0 ][ 0 ] = -2*x;
      out[  0 ][ 1 ] = -2*y;
      out[  0 ][ 2 ] = 2 - 2*z;

      out[  1 ][ 0 ] = 2*x;
      out[  1 ][ 1 ] = 2*y - 2;
      out[  1 ][ 2 ] = 2*z;

      out[  2 ][ 0 ] = 2 - 2*x;
      out[  2 ][ 1 ] = -2*y;
      out[  2 ][ 2 ] = -2*z;

      out[  3 ][ 0 ] = 2*x;
      out[  3 ][ 1 ] = 2*y;
      out[  3 ][ 2 ] = 2*z;
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateJacobian (const DomainType& in, std::vector<JacobianType>& out) const
    {
      out.resize(size());

      out[  0 ][ 0 ] = { -2, 0, 0 };
      out[  0 ][ 1 ] = { 0, -2, 0 };
      out[  0 ][ 2 ] = { 0, 0, -2 };

      out[  1 ][ 0 ] = { 2, 0, 0 };
      out[  1 ][ 1 ] = { 0, 2, 0 };
      out[  1 ][ 2 ] = { 0, 0, 2 };

      out[  2 ][ 0 ] = { -2, 0, 0 };
      out[  2 ][ 1 ] = { 0, -2, 0 };
      out[  2 ][ 2 ] = { 0, 0, -2 };

      out[  3 ][ 0 ] = { 2, 0, 0 };
      out[  3 ][ 1 ] = { 0, 2, 0 };
      out[  3 ][ 2 ] = { 0, 0, 2 };
    }

    /**
     * \brief Evaluate all partial derivatives of all shape functions
     *
     * \param order  order the partial derivative
     * \param in     Position
     * \param out    return value
     */
    void partial (const std::array<unsigned int, dim>& order,
                  const DomainType& in,
                  std::vector<RangeType>& out) const
    {
      if (std::accumulate(order.begin(), order.end(), 0) == 0)
        evaluateFunction(in, out);
      else
      {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        switch (direction) {
        case 0:
          out[  0 ] = { -2, 0, 0 };
          out[  1 ] = { 2, 0, 0 };
          out[  2 ] = { -2, 0, 0 };
          out[  3 ] = { 2, 0, 0 };
          break;

        case 1:
          out[  0 ] = { 0, -2, 0 };
          out[  1 ] = { 0, 2, 0 };
          out[  2 ] = { 0, -2, 0 };
          out[  3 ] = { 0, 2, 0 };
          break;

        case 2:
          out[  0 ] = { 0, 0, 2 };
          out[  1 ] = { 0, 0, -2 };
          out[  2 ] = { 0, 0, 2 };
          out[  3 ] = { 0, 0, -2 };
          break;

        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const { return 0; }

  private:
    //std::array<R, dim+1> s_;
  };

  /**
   * \brief First order Brezzi-Douglas-Fortin-Marini shape functions on the reference tetrahedron.
   *
   * \nosubgrouping
   */
  template<class D, class R >
  class BDFMSimplexLocalBasis<D, R, 3, 1>
  {
    static constexpr int dim = 3;
    using DomainType    = FieldVector<D, dim>;
    using RangeType     = FieldVector<R, dim>;
    using JacobianType  = FieldMatrix<R, dim, dim>;

  public:
    using Traits = LocalBasisTraits<D, dim, DomainType, R, dim, RangeType, JacobianType>;

    //! \brief Standard constructor
    BDFMSimplexLocalBasis ()
    {
      //std::fill(s_.begin(), s_.end(), 1);
    }

    /**
     * \brief Make set number s, where 0<= s < 16
     *
     * \param s  Edge orientation indicator
     */
    BDFMSimplexLocalBasis (int s)
    {
      /*
      for (auto i : range(3))
        s_[i] = s[i] ? -1 : 1;
        */
    }

    //! \brief number of shape functions
    unsigned int size () const { return 18; }

    /**
     * \brief Evaluate all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateFunction (const DomainType& in, std::vector<RangeType>& out) const
    {
      out.resize(size());

      const auto& x = in[0];
      const auto& y = in[1];
      const auto& z = in[2];

      out[  0 ][ 0 ] = 18*x*x + 18*x*y - 6*x*(x + y + z) - 6*x;
      out[  0 ][ 1 ] = 18*x*y + 18*y*y - 6*y*(x + y + z) - 6*y;
      out[  0 ][ 2 ] = 18*x*z - 24*x + 18*y*z - 24*y + 66*z*(x + y + z) - 84*z + 18;

      out[  1 ][ 0 ] = -6*x*x*0.2 + 126*x*y*0.2 - 144*x*(x + y) + 594*x*(x + y + z)*0.2 - 66*x*0.2;
      out[  1 ][ 1 ] = -6*x*y*0.2 + 126*y*y*0.2 - 54*y*(x + y + z)*0.2 - 6*y*0.2;
      out[  1 ][ 2 ] = -144*x*y - 6*x*z*0.2 + 264*x*0.2 + 126*y*z*0.2 + 144*y*0.2 - 54*z*(x + y + z) + 192*z*0.2 - 66*0.2;

      out[  2 ][ 0 ] = -234*x*x*0.2 - 366*x*y*0.2 + 144*x*(x + y) - 414*x*(x + y + z)*0.2 + 66*x*0.2;
      out[  2 ][ 1 ] = -234*x*y*0.2 - 366*y*y*0.2 + 234*y*(x + y + z)*0.2 + 6*y*0.2;
      out[  2 ][ 2 ] = 144*x*y - 234*x*z*0.2 - 144*x*0.2 - 366*y*z*0.2 - 24*y*0.2 + 18*z*(x + y + z) + 48*z*0.2 + 6*0.2;

      out[  3 ][ 0 ] = 18*x*y - 12*x*(x + y + z) + 6*x;
      out[  3 ][ 1 ] = 24*x + 18*y*y - 84*y*(x + y + z) + 84*y + 24*z - 18;
      out[  3 ][ 2 ] = 18*y*z - 12*z*(x + y + z) + 6*z;

      out[  4 ][ 0 ] = 384*x*x*0.2 + 126*x*y*0.2 - 144*x*(x + y) + 324*x*(x + y + z)*0.2 - 6*x*0.2;
      out[  4 ][ 1 ] = 384*x*y*0.2 - 24*x + 126*y*y*0.2 + 36*y*(x + y + z)*0.2 - 156*y*0.2 + 6;
      out[  4 ][ 2 ] = -144*x*y + 384*x*z*0.2 + 144*x*0.2 + 126*y*z*0.2 + 144*y*0.2 - 36*z*(x + y + z) + 42*z*0.2 - 36*0.2;

      out[  5 ][ 0 ] = -384*x*x*0.2 - 366*x*y*0.2 + 144*x*(x + y) - 264*x*(x + y + z)*0.2 + 6*x*0.2;
      out[  5 ][ 1 ] = -384*x*y*0.2 - 366*y*y*0.2 + 384*y*(x + y + z)*0.2 - 84*y*0.2 - 24*z + 6;
      out[  5 ][ 2 ] = 144*x*y - 384*x*z*0.2 - 144*x*0.2 - 366*y*z*0.2 - 144*y*0.2 + 48*z*(x + y + z) - 42*z*0.2 + 36*0.2;

      out[  6 ][ 0 ] = -18*x*x + 84*x*(x + y + z) - 84*x - 24*y - 24*z + 18;
      out[  6 ][ 1 ] = -18*x*y + 12*y*(x + y + z) - 6*y;
      out[  6 ][ 2 ] = -18*x*z + 12*z*(x + y + z) - 6*z;

      out[  7 ][ 0 ] = 234*x*x*0.2 - 24*x*y*0.2 - 144*x*(x + y) + 324*x*(x + y + z)*0.2 + 84*x*0.2 + 24*y - 6;
      out[  7 ][ 1 ] = 234*x*y*0.2 - 24*y*y*0.2 + 36*y*(x + y + z)*0.2 - 66*y*0.2;
      out[  7 ][ 2 ] = -144*x*y + 234*x*z*0.2 + 144*x*0.2 - 24*y*z*0.2 + 144*y*0.2 - 36*z*(x + y + z) + 102*z*0.2 - 36*0.2;

      out[  8 ][ 0 ] = 6*x*x*0.2 + 24*x*y*0.2 + 144*x*(x + y) - 744*x*(x + y + z)*0.2 + 156*x*0.2 + 24*z - 6;
      out[  8 ][ 1 ] = 6*x*y*0.2 + 24*y*y*0.2 - 96*y*(x + y + z)*0.2 + 66*y*0.2;
      out[  8 ][ 2 ] = 144*x*y + 6*x*z*0.2 - 144*x*0.2 + 24*y*z*0.2 - 144*y*0.2 + 24*z*(x + y + z) - 102*z*0.2 + 36*0.2;

      out[  9 ][ 0 ] = 12*x*x + 54*x*(x + y + z) - 48*x;
      out[  9 ][ 1 ] = 12*x*y - 18*y*(x + y + z) + 12*y;
      out[  9 ][ 2 ] = 12*x*z - 18*z*(x + y + z) + 12*z;

      out[ 10 ][ 0 ] = 12*x*y - 18*x*(x + y + z) + 12*x;
      out[ 10 ][ 1 ] = 12*y*y + 54*y*(x + y + z) - 48*y;
      out[ 10 ][ 2 ] = 12*y*z - 18*z*(x + y + z) + 12*z;

      out[ 11 ][ 0 ] = -12*x*x - 12*x*y - 6*x*(x + y + z) + 12*x;
      out[ 11 ][ 1 ] = -12*x*y - 12*y*y - 6*y*(x + y + z) + 12*y;
      out[ 11 ][ 2 ] = -12*x*z - 12*y*z + 66*z*(x + y + z) - 48*z;

      out[ 12 ][ 0 ] = 60*x*x - 300*x*(x + y + z) + 240*x;
      out[ 12 ][ 1 ] = 60*x*y + 60*y*(x + y + z) - 60*y;
      out[ 12 ][ 2 ] = 60*x*z + 60*z*(x + y + z) - 60*z;

      out[ 13 ][ 0 ] = 60*x*y + 60*x*(x + y + z) - 60*x;
      out[ 13 ][ 1 ] = 60*y*y - 300*y*(x + y + z) + 240*y;
      out[ 13 ][ 2 ] = 60*y*z + 60*z*(x + y + z) - 60*z;

      out[ 14 ][ 0 ] = -60*x*x - 60*x*y + 120*x*(x + y + z) - 60*x;
      out[ 14 ][ 1 ] = -60*x*y - 60*y*y + 120*y*(x + y + z) - 60*y;
      out[ 14 ][ 2 ] = -60*x*z - 60*y*z - 240*z*(x + y + z) + 240*z;

      out[ 15 ][ 0 ] = 264*x*x + 96*x*y - 720*x*(x + y) + 324*x*(x + y + z) + 24*x;
      out[ 15 ][ 1 ] = 264*x*y + 96*y*y + 36*y*(x + y + z) - 96*y;
      out[ 15 ][ 2 ] = -720*x*y + 264*x*z + 144*x + 96*y*z + 144*y - 180*z*(x + y + z) + 72*z - 36;

      out[ 16 ][ 0 ] = -84*x*x - 96*x*y + 720*x*(x + y) - 624*x*(x + y + z) + 96*x;
      out[ 16 ][ 1 ] = -84*x*y - 96*y*y + 24*y*(x + y + z) + 36*y;
      out[ 16 ][ 2 ] = 720*x*y - 84*x*z - 144*x - 96*y*z - 144*y + 240*z*(x + y + z) - 132*z + 36;

      out[ 17 ][ 0 ] = 264*x*x + 276*x*y - 720*x*(x + y) + 384*x*(x + y + z) - 36*x;
      out[ 17 ][ 1 ] = 264*x*y + 276*y*y - 264*y*(x + y + z) + 24*y;
      out[ 17 ][ 2 ] = -720*x*y + 264*x*z + 144*x + 276*y*z + 144*y - 120*z*(x + y + z) + 12*z - 36;
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param in   Position
     * \param out  return value
     */
    inline void evaluateJacobian (const DomainType& in, std::vector<JacobianType>& out) const
    {
      out.resize(size());

      const auto& x = in[0];
      const auto& y = in[1];
      const auto& z = in[2];

      out[  0 ][ 0 ] = { 24*x + 12*y - 6*z - 6, 12*x, -6*x };
      out[  0 ][ 1 ] = { 12*y, 12*x + 24*y - 6*z - 6, -6*y };
      out[  0 ][ 2 ] = { 84*z - 24, 84*z - 24, 84*x + 84*y + 132*z - 84 };

      out[  1 ][ 0 ] = { -264*x*0.2 + 594*z*0.2 - 66*0.2, 0, 594*x*0.2 };
      out[  1 ][ 1 ] = { -12*y, -12*x + 144*y*0.2 - 54*z*0.2 - 6*0.2, -54*y*0.2 };
      out[  1 ][ 2 ] = { -144*y - 276*z*0.2 + 264*0.2, -144*x - 144*z*0.2 + 144*0.2, -276*x*0.2 - 144*y*0.2 - 108*z + 192*0.2 };

      out[  2 ][ 0 ] = { 144*x*0.2 - 12*y - 414*z*0.2 + 66*0.2, -12*x, -414*x*0.2 };
      out[  2 ][ 1 ] = { 0, -264*y*0.2 + 234*z*0.2 + 6*0.2, 234*y*0.2 };
      out[  2 ][ 2 ] = { 144*y - 144*z*0.2 - 144*0.2, 144*x - 276*z*0.2 - 24*0.2, -144*x*0.2 - 276*y*0.2 + 36*z + 48*0.2 };

      out[  3 ][ 0 ] = { -24*x + 6*y - 12*z + 6, 6*x, -12*x };
      out[  3 ][ 1 ] = { 24 - 84*y, -84*x - 132*y - 84*z + 84, 24 - 84*y };
      out[  3 ][ 2 ] = { -12*z, 6*z, -12*x + 6*y - 24*z + 6 };

      out[  4 ][ 0 ] = { -24*x*0.2 - 54*y + 324*z*0.2 - 6*0.2, -54*x, 324*x*0.2 };
      out[  4 ][ 1 ] = { 84*y - 24, 84*x + 324*y*0.2 + 36*z*0.2 - 156*0.2, 36*y*0.2 };
      out[  4 ][ 2 ] = { -144*y + 204*z*0.2 + 144*0.2, -144*x - 54*z*0.2 + 144*0.2, 204*x*0.2 - 54*y*0.2 - 72*z + 42*0.2 };

      out[  5 ][ 0 ] = { 144*x*0.2 + 18*y - 264*z*0.2 + 6*0.2, 18*x, -264*x*0.2 };
      out[  5 ][ 1 ] = { 0, 36*y*0.2 + 384*z*0.2 - 84*0.2, 384*y*0.2 - 24 };
      out[  5 ][ 2 ] = { 144*y - 144*z*0.2 - 144*0.2, 144*x - 126*z*0.2 - 144*0.2, -144*x*0.2 - 126*y*0.2 + 96*z - 42*0.2 };

      out[  6 ][ 0 ] = { 132*x + 84*y + 84*z - 84, 84*x - 24, 84*x - 24 };
      out[  6 ][ 1 ] = { -6*y, -6*x + 24*y + 12*z - 6, 12*y };
      out[  6 ][ 2 ] = { -6*z, 12*z, -6*x + 12*y + 24*z - 6 };

      out[  7 ][ 0 ] = { -324*x*0.2 - 84*y + 324*z*0.2 + 84*0.2, 24 - 84*x, 324*x*0.2 };
      out[  7 ][ 1 ] = { 54*y, 54*x + 24*y*0.2 + 36*z*0.2 - 66*0.2, 36*y*0.2 };
      out[  7 ][ 2 ] = { -144*y + 54*z*0.2 + 144*0.2, -144*x - 204*z*0.2 + 144*0.2, 54*x*0.2 - 204*y*0.2 - 72*z + 102*0.2 };

      out[  8 ][ 0 ] = { -36*x*0.2 - 744*z*0.2 + 156*0.2, 0, 24 - 744*x*0.2 };
      out[  8 ][ 1 ] = { -18*y, -18*x - 144*y*0.2 - 96*z*0.2 + 66*0.2, -96*y*0.2 };
      out[  8 ][ 2 ] = { 144*y + 126*z*0.2 - 144*0.2, 144*x + 144*z*0.2 - 144*0.2, 126*x*0.2 + 144*y*0.2 + 48*z - 102*0.2 };

      out[  9 ][ 0 ] = { 132*x + 54*y + 54*z - 48, 54*x, 54*x };
      out[  9 ][ 1 ] = { -6*y, -6*x - 36*y - 18*z + 12, -18*y };
      out[  9 ][ 2 ] = { -6*z, -18*z, -6*x - 18*y - 36*z + 12 };


      out[ 10 ][ 0 ] = { -36*x - 6*y - 18*z + 12, -6*x, -18*x };
      out[ 10 ][ 1 ] = { 54*y, 54*x + 132*y + 54*z - 48, 54*y };
      out[ 10 ][ 2 ] = { -18*z, -6*z, -18*x - 6*y - 36*z + 12 };

      out[ 11 ][ 0 ] = { -36*x - 18*y - 6*z + 12, -18*x, -6*x };
      out[ 11 ][ 1 ] = { -18*y, -18*x - 36*y - 6*z + 12, -6*y };
      out[ 11 ][ 2 ] = { 54*z, 54*z, 54*x + 54*y + 132*z - 48 };

      out[ 12 ][ 0 ] = { -480*x - 300*y - 300*z + 240, -300*x, -300*x };
      out[ 12 ][ 1 ] = { 120*y, 120*x + 120*y + 60*z - 60, 60*y };
      out[ 12 ][ 2 ] = { 120*z, 60*z, 120*x + 60*y + 120*z - 60 };

      out[ 13 ][ 0 ] = { 120*x + 120*y + 60*z - 60, 120*x, 60*x };
      out[ 13 ][ 1 ] = { -300*y, -300*x - 480*y - 300*z + 240, -300*y };
      out[ 13 ][ 2 ] = { 60*z, 120*z, 60*x + 120*y + 120*z - 60 };

      out[ 14 ][ 0 ] = { 120*x + 60*y + 120*z - 60, 60*x, 120*x };
      out[ 14 ][ 1 ] = { 60*y, 60*x + 120*y + 120*z - 60, 120*y };
      out[ 14 ][ 2 ] = { -300*z, -300*z, -300*x - 300*y - 480*z + 240 };

      out[ 15 ][ 0 ] = { -264*x - 300*y + 324*z + 24, -300*x, 324*x };
      out[ 15 ][ 1 ] = { 300*y, 300*x + 264*y + 36*z - 96, 36*y };
      out[ 15 ][ 2 ] = { -720*y + 84*z + 144, -720*x - 84*z + 144, 84*x - 84*y - 360*z + 72 };

      out[ 16 ][ 0 ] = { 24*x - 624*z + 96, 0, -624*x };
      out[ 16 ][ 1 ] = { -60*y, -60*x - 144*y + 24*z + 36, 24*y };
      out[ 16 ][ 2 ] = { 720*y + 156*z - 144, 720*x + 144*z - 144, 156*x + 144*y + 480*z - 132 };

      out[ 17 ][ 0 ] = { -144*x - 60*y + 384*z - 36, -60*x, 384*x };
      out[ 17 ][ 1 ] = { 0, 24*y - 264*z + 24, -264*y };
      out[ 17 ][ 2 ] = { -720*y + 144*z + 144, -720*x + 156*z + 144, 144*x + 156*y - 240*z + 12 };

    }

    /**
     * \brief Evaluate all partial derivatives of all shape functions
     *
     * \param order  order the partial derivative
     * \param in     Position
     * \param out    return value
     */
    void partial (const std::array<unsigned int, dim>& order,
                  const DomainType& in,
                  std::vector<RangeType>& out) const
    {
      if (std::accumulate(order.begin(), order.end(), 0) == 0)
        evaluateFunction(in, out);
      else
      {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        const auto& x = in[0];
        const auto& y = in[1];
        const auto& z = in[2];

        switch (direction) {
        case 0:
          out[  0 ] = { 24*x + 12*y - 6*z - 6, 12*x, 384*x };
          out[  1 ] = { -264*x*0.2 + 594*z*0.2 - 66*0.2, 0, -6*x };
          out[  2 ] = { 144*x*0.2 - 12*y - 414*z*0.2 + 66*0.2, -12*x, 594*x*0.2 };
          out[  3 ] = { -24*x + 6*y - 12*z + 6, 6*x, -414*x*0.2 };
          out[  4 ] = { -24*x*0.2 - 54*y + 324*z*0.2 - 6*0.2, -54*x, -12*x };
          out[  5 ] = { 144*x*0.2 + 18*y - 264*z*0.2 + 6*0.2, 18*x, 324*x*0.2 };
          out[  6 ] = { 132*x + 84*y + 84*z - 84, 84*x - 24, -264*x*0.2 };
          out[  7 ] = { -324*x*0.2 - 84*y + 324*z*0.2 + 84*0.2, 24 - 84*x, 84*x - 24 };
          out[  8 ] = { -36*x*0.2 - 744*z*0.2 + 156*0.2, 0, 324*x*0.2 };
          out[  9 ] = { 132*x + 54*y + 54*z - 48, 54*x, 24 - 744*x*0.2 };
          out[ 10 ] = { -36*x - 6*y - 18*z + 12, -6*x, 54*x };
          out[ 11 ] = { -36*x - 18*y - 6*z + 12, -18*x, -18*x };
          out[ 12 ] = { -480*x - 300*y - 300*z + 240, -300*x, -6*x };
          out[ 13 ] = { 120*x + 120*y + 60*z - 60, 120*x, -300*x };
          out[ 14 ] = { 120*x + 60*y + 120*z - 60, 60*x, 60*x };
          out[ 15 ] = { -264*x - 300*y + 324*z + 24, -300*x, 120*x };
          out[ 16 ] = { 24*x - 624*z + 96, 0, 324*x };
          out[ 17 ] = { -144*x - 60*y + 384*z - 36, -60*x, -624*x };
          break;

        case 1:
          out[  0 ] = { 12*y, 12*x + 24*y - 6*z - 6, -264*y };
          out[  1 ] = { -12*y, -12*x + 144*y*0.2 - 54*z*0.2 - 6*0.2, -6*y };
          out[  2 ] = { 0, -264*y*0.2 + 234*z*0.2 + 6*0.2, -54*y*0.2 };
          out[  3 ] = { 24 - 84*y, -84*x - 132*y - 84*z + 84, 234*y*0.2 };
          out[  4 ] = { 84*y - 24, 84*x + 324*y*0.2 + 36*z*0.2 - 156*0.2, 24 - 84*y };
          out[  5 ] = { 0, 36*y*0.2 + 384*z*0.2 - 84*0.2, 36*y*0.2 };
          out[  6 ] = { -6*y, -6*x + 24*y + 12*z - 6, 384*y*0.2 - 24 };
          out[  7 ] = { 54*y, 54*x + 24*y*0.2 + 36*z*0.2 - 66*0.2, 12*y };
          out[  8 ] = { -18*y, -18*x - 144*y*0.2 - 96*z*0.2 + 66*0.2, 36*y*0.2 };
          out[  9 ] = { -6*y, -6*x - 36*y - 18*z + 12, -96*y*0.2 };
          out[ 10 ] = { 54*y, 54*x + 132*y + 54*z - 48, -18*y };
          out[ 11 ] = { -18*y, -18*x - 36*y - 6*z + 12, 54*y };
          out[ 12 ] = { 120*y, 120*x + 120*y + 60*z - 60, -6*y };
          out[ 13 ] = { -300*y, -300*x - 480*y - 300*z + 240, 60*y };
          out[ 14 ] = { 60*y, 60*x + 120*y + 120*z - 60, -300*y };
          out[ 15 ] = { 300*y, 300*x + 264*y + 36*z - 96, 120*y };
          out[ 16 ] = { -60*y, -60*x - 144*y + 24*z + 36, 36*y };
          out[ 17 ] = { 0, 24*y - 264*z + 24, 24*y };
          break;

        case 2:
          out[  0 ] = { 84*z - 24, 84*z - 24, 144*x + 156*y - 240*z + 12 };
          out[  1 ] = { -144*y - 276*z*0.2 + 264*0.2, -144*x - 144*z*0.2 + 144*0.2, 84*x + 84*y + 132*z - 84 };
          out[  2 ] = { 144*y - 144*z*0.2 - 144*0.2, 144*x - 276*z*0.2 - 24*0.2, -276*x*0.2 - 144*y*0.2 - 108*z + 192*0.2 };
          out[  3 ] = { -12*z, 6*z, -144*x*0.2 - 276*y*0.2 + 36*z + 48*0.2 };
          out[  4 ] = { -144*y + 204*z*0.2 + 144*0.2, -144*x - 54*z*0.2 + 144*0.2, -12*x + 6*y - 24*z + 6 };
          out[  5 ] = { 144*y - 144*z*0.2 - 144*0.2, 144*x - 126*z*0.2 - 144*0.2, 204*x*0.2 - 54*y*0.2 - 72*z + 42*0.2 };
          out[  6 ] = { -6*z, 12*z, -144*x*0.2 - 126*y*0.2 + 96*z - 42*0.2 };
          out[  7 ] = { -144*y + 54*z*0.2 + 144*0.2, -144*x - 204*z*0.2 + 144*0.2, -6*x + 12*y + 24*z - 6 };
          out[  8 ] = { 144*y + 126*z*0.2 - 144*0.2, 144*x + 144*z*0.2 - 144*0.2, 54*x*0.2 - 204*y*0.2 - 72*z + 102*0.2 };
          out[  9 ] = { -6*z, -18*z, 126*x*0.2 + 144*y*0.2 + 48*z - 102*0.2 };
          out[ 10 ] = { -18*z, -6*z, -6*x - 18*y - 36*z + 12 };
          out[ 11 ] = { 54*z, 54*z, -18*x - 6*y - 36*z + 12 };
          out[ 12 ] = { 120*z, 60*z, 54*x + 54*y + 132*z - 48 };
          out[ 13 ] = { 60*z, 120*z, 120*x + 60*y + 120*z - 60 };
          out[ 14 ] = { -300*z, -300*z, 60*x + 120*y + 120*z - 60 };
          out[ 15 ] = { -720*y + 84*z + 144, -720*x - 84*z + 144, -300*x - 300*y - 480*z + 240 };
          out[ 16 ] = { 720*y + 156*z - 144, 720*x + 144*z - 144, 84*x - 84*y - 360*z + 72 };
          out[ 17 ] = { -720*y + 144*z + 144, -720*x + 156*z + 144, 156*x + 144*y + 480*z - 132 };
          break;

        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const { return 0; }

  private:
    //std::array<R, dim+1> s_;
  };


} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_SIMPLEX_LOCALBASIS_HH
