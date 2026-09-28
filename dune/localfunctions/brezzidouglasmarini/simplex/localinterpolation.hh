// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_SIMPLEX_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_SIMPLEX_LOCALINTERPOLATION_HH

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


namespace Dune
{

  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief Interpolation for Brezzi-Douglas-Marini shape functions on simplices.
   *
   * \tparam D      Type of represent the field in the domain.
   * \tparam R      Type of represent the field in the domain.
   * \tparam dim    dimension of the reference element, must be >= 2.
   * \tparam order  order of the element, must be >= 1.
   *
   * \nosubgrouping
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BDMSimplexLocalInterpolation
  {
    static_assert( AlwaysFalse<D>::value,
                  "`BDMSimplexLocalInterpolation` not implemented for chosen `dim` and `order`." );
  };


  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief First order Brezzi-Douglas-Marini shape functions on the reference triangle.
   *
   * \tparam D      Type of represent the field in the domain.
   * \tparam R      Type of represent the field in the domain.
   * \tparam dim    dimension of the reference element, must be >= 2.
   * \tparam order  order of the element, must be >= 1.
   *
   * \nosubgrouping
   */
  template<class D, class R>
  class BDMSimplexLocalInterpolation< D, R, 2, 1>
  {
    static constexpr int dim = 2;
    using DomainType      = FieldVector<D, dim>;
    using FaceDomainType  = FieldVector<D, dim-1>;
    using RangeType       = FieldVector<R, dim>;
    using DomainFieldType = D;
    using RangeFieldType  = R;

    using Scalar = RangeFieldType;
  public:
    //! \brief Standard constructor
    BDMSimplexLocalInterpolation ()
    {
      sign0 = sign1 = sign2 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    BDMSimplexLocalInterpolation (unsigned int s)
    {
      using std::sqrt;
      sign0 = sign1 = sign2 = 1.0;
      if (s & 1)
      {
        sign0 = -1.0;
      }
      if (s & 2)
      {
        sign1 = -1.0;
      }
      if (s & 4)
      {
        sign2 = -1.0;
      }
    }

    /**
     * \brief Interpolate a given function with shape functions
     *
     * \tparam F Function type for function which should be interpolated
     * \tparam C Coefficient type
     * \param f function which should be interpolated
     * \param out return value, vector of coefficients
     */
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      out.resize(6);
      std::fill(out.begin(), out.end(), 0.0);

      const int qOrder = 4;
      const Dune::QuadratureRule<Scalar,1>& rule = Dune::QuadratureRules<Scalar,1>::rule(Dune::GeometryTypes::simplex(1), qOrder);

      DomainType localPos;
      for (typename Dune::QuadratureRule<Scalar,1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        const Scalar qPos = it->position();
        const Scalar weight_0 = it->weight();
        const Scalar weight_1 = (2.0*qPos - 1.0)*it->weight();

        localPos[0] = qPos;
        localPos[1] = 0.0;
        auto y = f(localPos);
        // n0 = (0, -1)
        out[0] -= y[1]*weight_0*sign0;
        out[3] -= y[1]*weight_1;

        localPos[0] = 0.0;
        localPos[1] = qPos;
        y = f(localPos);
        // n1 = (-1, 0)
        out[1] -= y[0]*weight_0*sign1;
        out[4] += y[0]*weight_1;

        localPos[0] = 1.0 - qPos;
        localPos[1] = qPos;
        y = f(localPos);
        // n2 = (1.0, 1.0))
        const Scalar scp = y[0] + y[1];
        out[2] += scp*weight_0*sign2;
        out[5] += scp*weight_1;
      }
    }

  private:
    RangeFieldType sign0,sign1,sign2;
  };

  /**
   * \brief Second order Brezzi-Douglas-Marini shape functions on triangles.
   *
   * \tparam F Function type for function which should be interpolated
   * \tparam C Coefficient type
   * \param f function which should be interpolated
   * \param out return value, vector of coefficients
   *
   * \ingroup BrezziDouglasMariniImpl
   * \nosubgrouping
   */
  template<class D, class R>
  class BDMSimplexLocalInterpolation< D, R, 2, 2>
  {
    static constexpr int dim = 2;
    using DomainType      = FieldVector<D, dim>;
    using FaceDomainType  = FieldVector<D, dim-1>;
    using RangeType       = FieldVector<R, dim>;
    using DomainFieldType = D;
    using RangeFieldType  = R;

    using Scalar = RangeFieldType;
  public:
    //! \brief Standard constructor
    BDMSimplexLocalInterpolation()
    {
      sign0 = sign1 = sign2 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    BDMSimplexLocalInterpolation(unsigned int s)
    {
      sign0 = sign1 = sign2 = 1.0;
      if (s & 1)
      {
        sign0 = -1.0;
      }
      if (s & 2)
      {
        sign1 = -1.0;
      }
      if (s & 4)
      {
        sign2 = -1.0;
      }
    }

    /**
     * \brief Interpolate a given function with shape functions
     *
     * \tparam F Function type for function which should be interpolated
     * \tparam C Coefficient type
     * \param f function which should be interpolated
     * \param out return value, vector of coefficients
     */
    template<typename F, typename C>
    void interpolate(const F& f, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      out.resize(12);
      std::fill(out.begin(), out.end(), 0.0);

      const int qOrder = 4;
      const Dune::QuadratureRule<Scalar,1>& rule = Dune::QuadratureRules<Scalar,1>::rule(Dune::GeometryTypes::simplex(1), qOrder);

      DomainType localPos;
      for (typename Dune::QuadratureRule<Scalar,1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        const Scalar qPos = it->position();
        const Scalar weight_0 = it->weight();
        const Scalar weight_1 = (1.0 - 2.0*qPos)*it->weight();
        const Scalar weight_2 = (6.0*qPos*qPos - 6.0*qPos + 1.0)*it->weight();


        localPos[0] = qPos;
        localPos[1] = 0.0;
        auto y = f(localPos);
        // n0 = (0, -1)
        out[0] -= y[1]*weight_0*sign0;
        out[1] -= y[1]*weight_1;
        out[2] -= y[1]*weight_2*sign0;

        localPos[0] = 0.0;
        localPos[1] = qPos;
        y = f(localPos);
        // n1 = (-1, 0)
        out[3] -= y[0]*weight_0*sign1;
        out[4] += y[0]*weight_1;
        out[5] -= y[0]*weight_2*sign1;

        localPos[0] = 1.0 - qPos;
        localPos[1] = qPos;
        y = f(localPos);
        // n2 = (1, 1)
        const Scalar scp = y[0] + y[1];
        out[6] += scp*weight_0*sign2;
        out[7] += scp*weight_1;
        out[8] += scp*weight_2*sign2;
      }

      // a volume part is needed here for dofs: 9 10 11
      const QuadratureRule<D,2>& rule2 = QuadratureRules<D,2>::rule(GeometryTypes::simplex(2), qOrder);

      for (typename QuadratureRule<D,2>::const_iterator it=rule2.begin(); it!=rule2.end(); ++it)
      {
        const DomainType localPos = it->position();
        auto y = f(localPos);

        out[9] += y[0]*it->weight();
        out[10] += y[1]*it->weight();
        out[11] += (y[0]*(localPos[0]-2.0*localPos[0]*localPos[1]-localPos[0]*localPos[0])
            +y[1]*(-localPos[1]+2.0*localPos[0]*localPos[1]+localPos[1]*localPos[1]))*it->weight();
      }
    }

  private:
    RangeFieldType sign0, sign1, sign2;
  };


  /**
   * \brief First order Brezzi-Douglas-Marini shape functions on the reference
   *        tetrahedron.
   *
   * \tparam F Function type for function which should be interpolated
   * \tparam C Coefficient type
   * \param f function which should be interpolated
   * \param out return value, vector of coefficients
   *
   * \ingroup BrezziDouglasMariniImpl
   * \nosubgrouping
   */
  template<class D, class R>
  class BDMSimplexLocalInterpolation< D, R, 3, 1>
  {
    static constexpr int dim = 3;
    using DomainType      = FieldVector<D, dim>;
    using FaceDomainType  = FieldVector<D, dim-1>;
    using RangeType       = FieldVector<R, dim>;
    using DomainFieldType = D;
    using RangeFieldType  = R;

    using Scalar = RangeFieldType;

  public:
    //! \brief Standard constructor
    BDMSimplexLocalInterpolation()
    {
    }

    /**
     * \brief Make set number s, where 0 <= s < 64
     *
     * \param s Edge orientation indicator
     */
    BDMSimplexLocalInterpolation(unsigned int s)
      : BDMSimplexLocalInterpolation()
    {}

    /**
     * \brief Interpolate a given function with shape functions
     *
     * \tparam F Function type for function which should be interpolated
     * \tparam C Coefficient type
     * \param f function which should be interpolated
     * \param out return value, vector of coefficients
     */
    template<typename F, typename C>
    void interpolate(const F& f, std::vector<C>& out) const
    {
      out.resize(12);
      std::fill(out.begin(), out.end(), 0.0);

      const int qOrder = 4;
      static const int faceDim = 2;

      const Dune::QuadratureRule<Scalar,faceDim>& rule = Dune::QuadratureRules<Scalar,faceDim>::rule(Dune::GeometryTypes::simplex(faceDim), qOrder);

      DomainType localPos;

      for (typename Dune::QuadratureRule<Scalar,faceDim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        const auto qPos = it->position();

        const Scalar weight_0 = (1. - qPos[0] - qPos[1])*it->weight();
        const Scalar weight_1 = qPos[0]*it->weight();
        const Scalar weight_2 = qPos[1]*it->weight();

        // face0 (bottom)
        localPos = {qPos[0], qPos[1], 0.0};
        auto y = f(localPos);
        // n0 = (0, 0, 1)
        out[0] += y[2] * weight_0;
        out[1] += y[2] * weight_1;
        out[2] += y[2] * weight_2;

        // face1 (front)
        localPos = {qPos[0], 0.0, qPos[1]};
        y = f(localPos);
        // n1 = ( 0, -1,  0 )
        out[3] -= y[1] * weight_0;
        out[4] -= y[1] * weight_1;
        out[5] -= y[1] * weight_2;

        // face2 (left)
        localPos = {0.0, qPos[0], qPos[1]};
        y = f(localPos);
        // n2 = ( 1,  0,  0 )
        out[6] += y[0] * weight_0;
        out[7] += y[0] * weight_1;
        out[8] += y[0] * weight_2;

        // face3 (rear)
        localPos = {1.0 - qPos[0] - qPos[1], qPos[0], qPos[1]};
        y = f(localPos);
        // n3 = ( 1,  1,  1 )
        const Scalar scp = y[0] + y[1] + y[2];
        out[ 9] += scp * weight_0;
        out[10] += scp * weight_1;
        out[11] += scp * weight_2;
      }
    }
  };


#ifndef DOXYGEN
  template<class D, class R, unsigned int dim>
  class BDMSimplexLocalInterpolation<D, R, dim, 0>
  {
    static_assert(AlwaysFalse<D>::value,
                  "`BDMSimplexLocalInterpolation` not defined for order 0.");
  };
#endif //#ifndef DOXYGEN

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_SIMPLEX_LOCALINTERPOLATION_HH
