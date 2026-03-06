// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_SIMPLEX_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_SIMPLEX_LOCALINTERPOLATION_HH

#include <algorithm>
#include <array>
#include <bitset>
#include <vector>
#include <limits>

#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules.hh>

#ifdef BDFM_USE_RT0_BASIS
#include <dune/localfunctions/raviartthomas/raviartthomas02d/raviartthomas02dlocalbasis.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas02d/raviartthomas02dlocalinterpolation.hh>
#endif

namespace Dune
{

  /**
   * \ingroup BrezziDouglasFortinMariniImpl
   * \brief Interpolation for Brezzi-Douglas-Fortin-Marini shape functions on simplices.
   *
   * \tparam D      Type of represent the field in the domain.
   * \tparam R      Type of represent the field in the domain.
   * \tparam dim    dimension of the reference element, must be >= 2.
   * \tparam order  order of the element, must be 0 or 1.
   *
   * \nosubgrouping
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMSimplexLocalInterpolation;

#ifdef BDFM_USE_RT0_BASIS
  template<class D, class R>
  class BDFMSimplexLocalInterpolation< D, R, 2, 0>
    : public RT02DLocalInterpolation< RT02DLocalBasis<D,R> >
  {
    typedef RT02DLocalInterpolation< RT02DLocalBasis<D,R> > BaseType;
  public:
    BDFMSimplexLocalInterpolation() {}
    BDFMSimplexLocalInterpolation(std::bitset<3> s) : BaseType(s){}
  };
#endif

  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMSimplexLocalInterpolation
  {
    static_assert( order < 2 , "BDFMSimplexLocalInterpolation only implement for order < 2");

    using DomainType      = FieldVector<D, dim>;
    using FaceDomainType  = FieldVector<D, dim-1>;
    using RangeType       = FieldVector<R, dim>;
    using DomainFieldType = D;
    using RangeFieldType  = R;

    static constexpr unsigned int interiorDofs = order*(3*(dim-1));
    static constexpr unsigned int faceDofs     = order == 0 ? 1 : dim;

    static constexpr std::size_t numFaces = dim+1;
    static constexpr std::size_t numDofs  = numFaces*faceDofs + interiorDofs;

  public:
    //! \brief Standard constructor
    BDFMSimplexLocalInterpolation ()
    {
      //std::fill(sign_.begin(), sign_.end(), 1.0);
    }

    /**
     * \brief Make set number s, where 0 <= s < 2^(2*dim)
     *
     * \param s  Edge orientation indicator
     */
    BDFMSimplexLocalInterpolation (const int s)
    {
      /*
      for (auto i : range(numFaces))
      {
        sign_[i] = s[i] ? -1 : 1;
      }
      */
      //std::fill(sign_.begin(), sign_.end(), 1.0);
    }

    /**
     * \brief Interpolate a given function with shape functions
     *
     * \tparam F  Function type for function which should be interpolated
     * \tparam C  Coefficient vector type
     *
     * \param f   function which should be interpolated
     * \param out  return value, vector of coefficients
     */
    template<class F, class C>
    void interpolate (const F& f, C& out) const
    {
      out.resize(numDofs);
      std::fill(out.begin(),out.end(), 0.0);

      if constexpr ( dim == 2 )
      {
        interpolate2d( f, out );
        return ;
      }

      if constexpr ( dim == 3 )
      {
        interpolate3d( f, out );
        return ;
      }

      DUNE_THROW(NotImplemented,"interpolate for this dimension not implemented");
    }

  protected:
    // 2d implementation for order 0 and 1
    template<class F, class C>
    void interpolate2d (const F& f, C& out) const
    {
      const int qOrder = 4;

      DomainType localPos;
      D weight_0 = 0.0;
      D weight_1 = 0.0;

      const Dune::QuadratureRule<D,dim-1>& rule = Dune::QuadratureRules<D,dim-1>::rule(Dune::GeometryTypes::simplex(dim-1), qOrder);
      for(const auto& qp : rule)
      {
        const D qPos   = qp.position();

        if constexpr ( order == 1 )
        {
          weight_0 = (1.0 - qPos) * qp.weight();
          weight_1 = qPos * qp.weight();
        }

        localPos[0] = qPos;
        localPos[1] = 0.0;
        auto y = f(localPos);

        // n0 = (0, 1)
        if constexpr ( order == 0 )
        {
          out[0] += y[1] * qp.weight();
        }
        if constexpr ( order == 1 )
        {
          out[0] += y[1] * weight_0;
          out[1] += y[1] * weight_1;
        }

        localPos[0] = 0.0;
        localPos[1] = qPos;
        y = f(localPos);
        // n1 = ( -1, 0 )
        if constexpr ( order == 0 )
        {
          out[1] -= y[0] * qp.weight();
        }
        if constexpr ( order == 1 )
        {
          out[2] -= y[0] * weight_0;
          out[3] -= y[0] * weight_1;
        }

        localPos[0] = 1.0 - qPos;
        localPos[1] = qPos;
        y = f(localPos);
        // n2 = ( -1, -1 )
        if constexpr ( order == 0 )
        {
          out[2] -= (y[0] + y[1]) * qp.weight();
        }
        if constexpr ( order == 1 )
        {
          const D scp = (y[0] + y[1]);
          out[4] -= scp * weight_0;
          out[5] -= scp * weight_1;
        }
      }

      if constexpr( interiorDofs > 0 )
      {
        DomainType tangent;
        const auto& rule = QuadratureRules<D, dim>::rule(GeometryTypes::simplex(dim), order+qOrder);
        for(const auto& qp : rule)
        {
          const auto& x = qp.position();
          auto y = f(x);

          // t0 = ( 1-x_1, x_0 )
          tangent[0] = 1.0 - x[1];
          tangent[1] = x[0];
          out[6] += (y * tangent) * qp.weight();

          // t1 = ( x_1, 1 - x_0 )
          tangent[0] = x[1];
          tangent[1] = 1.0 - x[0];
          out[7] += (y * tangent) * qp.weight();

          // t2 = ( -x_1, x_0 )
          tangent[0] = -x[1];
          tangent[1] =  x[0];
          out[8] += (y * tangent) * qp.weight();
        }
      }
    }

    // 3d implementation for order 0
    template<class F, class C>
    void interpolate3d (const F& f, C& out) const
    {
      assert( dim == 3 );

      const int qOrder = 4;

      DomainType localPos;
      D weight_0 = 0.0;
      D weight_1 = 0.0;
      D weight_2 = 0.0;

      const Dune::QuadratureRule<D,dim-1>& rule = Dune::QuadratureRules<D,dim-1>::rule(Dune::GeometryTypes::simplex(dim-1), qOrder);
      for(const auto& qp : rule)
      {
        const auto& qPos   = qp.position();

        if constexpr ( order == 1 )
        {
          weight_0 = (1.0 - qPos[0] - qPos[1]) * qp.weight();
          weight_1 = qPos[0] * qp.weight();
          weight_2 = qPos[1] * qp.weight();
        }

        // face0 (bottom)
        localPos = {qPos[0], qPos[1], 0.0};
        auto y = f(localPos);

        // n0 = (0, 0, 1)
        if constexpr ( order == 0 )
        {
          out[0] += y[2] * qp.weight();
        }
        if constexpr ( order == 1 )
        {
          out[0] += y[2] * weight_0;
          out[1] += y[2] * weight_1;
          out[2] += y[2] * weight_2;
        }

        // face1 (front)
        localPos = {qPos[0], 0.0, qPos[1]};
        y = f(localPos);
        // n1 = ( 0, -1, 0 )
        if constexpr ( order == 0 )
        {
          out[1] -= y[1] * qp.weight();
        }
        if constexpr ( order == 1 )
        {
          out[3] -= y[1] * weight_0;
          out[4] -= y[1] * weight_1;
          out[5] -= y[1] * weight_2;
        }

        // face2 (left)
        localPos = {0.0, qPos[0], qPos[1]};
        y = f(localPos);
        // n2 = ( 1, 0, 0 )
        if constexpr ( order == 0 )
        {
          out[2] += y[0] * qp.weight();
        }
        if constexpr ( order == 1 )
        {
          out[6] += y[0] * weight_0;
          out[7] += y[0] * weight_1;
          out[8] += y[0] * weight_2;
        }

        // face3 (rear)
        localPos = {1.0 - qPos[0] - qPos[1], qPos[0], qPos[1]};
        y = f(localPos);
        // n3 = ( 1,  1,  1 )
        if constexpr ( order == 0 )
        {
          out[3] += (y[0] + y[1] + y[2]) * qp.weight();
        }
        if constexpr ( order == 1 )
        {
          const R scp = (y[0] + y[1] + y[2]);
          out[ 9] += scp * weight_0;
          out[10] += scp * weight_1;
          out[11] += scp * weight_2;
        }
      }

      if constexpr( interiorDofs > 0 )
      {
        DomainType tangent;
        const auto& rule = QuadratureRules<D, dim>::rule(GeometryTypes::simplex(dim), order+qOrder);
        for(const auto& qp : rule)
        {
          const auto& x = qp.position();
          auto y = f(x);

          // t0 = ( 1-x_0-x_1, x_0, x_0 )
          tangent = { 1.0 - x[1] - x[2], x[0], x[0] };
          out[12] += (y * tangent) * qp.weight();

          // t1 = ( x_1, 1 - x_0 -x_2, x_1 )
          tangent = { x[1], 1 - x[0] - x[2], x[1] };
          out[13] += (y * tangent) * qp.weight();

          // t2 = ( x_2, x_2, 1 - x_0 -x_1 )
          tangent = { x[2], x[2], 1 - x[0] - x[1] };
          out[14] += (y * tangent) * qp.weight();

          // t3 = ( -x_1, x_0, 0 )
          tangent = { -x[1], x[0], 0 };
          out[15] += (y * tangent) * qp.weight();

          // t4 = ( -x_2, 0, x_0 )
          tangent = { -x[2], 0, x[0] };
          out[16] += (y * tangent) * qp.weight();

          // t5 = ( 0, -x_2, x_1 )
          tangent = { 0, -x[2], x[1] };
          out[17] += (y * tangent) * qp.weight();
        }
      }
    }
  private:
    //std::array<RangeFieldType, numFaces> sign_;
  };


  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr std::size_t BDFMSimplexLocalInterpolation<D, R, dim, order>::numFaces;

  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr unsigned int BDFMSimplexLocalInterpolation<D, R, dim, order>::interiorDofs;

  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr unsigned int BDFMSimplexLocalInterpolation<D, R, dim, order>::faceDofs;

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_SIMPLEX_LOCALINTERPOLATION_HH
