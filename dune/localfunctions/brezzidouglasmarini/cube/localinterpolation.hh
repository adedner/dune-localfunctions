// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_CUBE_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_CUBE_LOCALINTERPOLATION_HH

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
  class BDMCubeLocalInterpolation
  {
    static_assert( AlwaysFalse<D>::value,
                  "`BDMCubeLocalInterpolation` not implemented for chosen `dim` and `order`." );
  };


  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief First order Brezzi-Douglas-Marini shape functions on the reference quadrilateral.
   *
   * \tparam D      Type of represent the field in the domain.
   * \tparam R      Type of represent the field in the domain.
   * \tparam dim    dimension of the reference element, must be >= 2.
   * \tparam order  order of the element, must be >= 1.
   *
   * \nosubgrouping
   */
  template<class D, class R>
  class BDMCubeLocalInterpolation< D, R, 2, 1>
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
    BDMCubeLocalInterpolation ()
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    BDMCubeLocalInterpolation (unsigned int s)
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
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
      if (s & 8)
      {
        sign3 = -1.0;
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
      // f gives v*outer normal at a point on the edge!
      out.resize(8);
      std::fill(out.begin(), out.end(), 0.0);

      const int qOrder = 4;
      const QuadratureRule<Scalar,1>& rule = QuadratureRules<Scalar,1>::rule(GeometryTypes::cube(1), qOrder);

      DomainType localPos;
      for (typename QuadratureRule<Scalar,1>::const_iterator it = rule.begin();
           it != rule.end(); ++it)
      {
        const Scalar qPos = it->position();
        const Scalar weight_0 = it->weight();
        const Scalar weight_1 = (2.0*qPos - 1.0)*it->weight();

        localPos[0] = 0.0;
        localPos[1] = qPos;
        auto y = f(localPos);
        // n0 = (-1, 0)
        // compute y * n0 * ...
        out[0] -= y[0]*weight_0*sign0;
        out[1] -= y[0]*weight_1;

        localPos[0] = 1.0;
        localPos[1] = qPos;
        y = f(localPos);
        // n1 = (1, 0)
        // compute y * n1 * ...
        out[2] += y[0]*weight_0*sign1;
        out[3] -= y[0]*weight_1;

        localPos[0] = qPos;
        localPos[1] = 0.0;
        y = f(localPos);
        // n2 = (0, -1)
        // compute y * n2 * ...
        out[4] -= y[1]*weight_0*sign2;
        out[5] += y[1]*weight_1;

        localPos[0] = qPos;
        localPos[1] = 1.0;
        y = f(localPos);
        // n2 = (0, 1)
        // compute y * n3 * ...
        out[6] += y[1]*weight_0*sign3;
        out[7] += y[1]*weight_1;
      }
    }

  private:
    RangeFieldType sign0, sign1, sign2, sign3;
  };

  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief Second order Brezzi-Douglas-Marini shape functions on the reference quadrilateral.
   *
   * \tparam D      Type of represent the field in the domain.
   * \tparam R      Type of represent the field in the domain.
   * \tparam dim    dimension of the reference element, must be >= 2.
   * \tparam order  order of the element, must be >= 1.
   *
   * \nosubgrouping
   */
  template<class D, class R>
  class BDMCubeLocalInterpolation< D, R, 2, 2>
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
    BDMCubeLocalInterpolation()
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    BDMCubeLocalInterpolation(unsigned int s)
    {
      sign0 = sign1 = sign2 = sign3 = 1.0;
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
      if (s & 8)
      {
        sign3 = -1.0;
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
      out.resize(14);
      std::fill(out.begin(), out.end(), 0.0);

      const int qOrder = 4;
      const QuadratureRule<Scalar,1>& rule = QuadratureRules<Scalar,1>::rule(GeometryTypes::cube(1), qOrder);

      DomainType localPos;

      for (typename QuadratureRule<Scalar,1>::const_iterator it = rule.begin();
           it != rule.end(); ++it)
      {
        const Scalar qPos = it->position();
        const Scalar weight_0 = it->weight();
        const Scalar weight_1 = (2.0*qPos - 1.0)*it->weight();
        const Scalar weight_2 = (8.0*qPos*qPos - 8.0*qPos + 1.0)*it->weight();

        localPos[0] = 0.0;
        localPos[1] = qPos;
        auto y = f(localPos);
        // n0 = (-1, 0)
        // compute y * n0 ...
        out[0] -= y[0]*weight_0*sign0;
        out[1] -= y[0]*weight_1;
        out[2] -= y[0]*weight_2*sign0;

        localPos[0] = 1.0;
        localPos[1] = qPos;
        y = f(localPos);
        // n1 = (1, 0)
        // compute y * n1 ...
        out[3] += y[0]*weight_0*sign1;
        out[4] -= y[0]*weight_1;
        out[5] += y[0]*weight_2*sign1;

        localPos[0] = qPos;
        localPos[1] = 0.0;
        y = f(localPos);
        // n2 = (0, -1)
        // compute y * n2 ...
        out[6] -= y[1]*weight_0*sign2;
        out[7] += y[1]*weight_1;
        out[8] -= y[1]*weight_2*sign2;

        localPos[0] = qPos;
        localPos[1] = 1.0;
        y = f(localPos);
        // n3 = (0, 1)
        // compute y * n3 ...
        out[ 9] += y[1]*weight_0*sign3;
        out[10] += y[1]*weight_1;
        out[11] += y[1]*weight_2*sign3;
      }

      const QuadratureRule<D,2>& rule2 = QuadratureRules<D,2>::rule(GeometryTypes::cube(2), qOrder);

      for (typename QuadratureRule<D,2>::const_iterator it=rule2.begin(); it!=rule2.end(); ++it)
      {
        auto y = f(it->position());
        out[12] += y[0]*it->weight();
        out[13] += y[1]*it->weight();
      }
    }

  private:
    RangeFieldType sign0, sign1, sign2, sign3;
  };


  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief First order Brezzi-Douglas-Marini shape functions on the reference
   *        hexahedron.
   *
   * \tparam D      Type of represent the field in the domain.
   * \tparam R      Type of represent the field in the domain.
   * \tparam dim    dimension of the reference element, must be >= 2.
   * \tparam order  order of the element, must be >= 1.
   *
   * \nosubgrouping
   */
  template<class D, class R>
  class BDMCubeLocalInterpolation< D, R, 3, 1>
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
    BDMCubeLocalInterpolation()
    {
      sign0 = sign1 = sign2 = sign3 = sign4 = sign5 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 64
     *
     * \param s Edge orientation indicator
     */
    BDMCubeLocalInterpolation(unsigned int s)
    {
      sign0 = sign1 = sign2 = sign3 = sign4 = sign5 = 1.0;
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
      if (s & 8)
      {
        sign3 = -1.0;
      }
      if (s & 16)
      {
        sign4 = -1.0;
      }
      if (s & 32)
      {
        sign5 = -1.0;
      }

      n0[0] = -1.0;
      n0[1] =  0.0;
      n0[2] =  0.0;
      n1[0] =  1.0;
      n1[1] =  0.0;
      n1[2] =  0.0;
      n2[0] =  0.0;
      n2[1] = -1.0;
      n2[2] =  0.0;
      n3[0] =  0.0;
      n3[1] =  1.0;
      n3[2] =  0.0;
      n4[0] =  0.0;
      n4[1] =  0.0;
      n4[2] = -1.0;
      n5[0] =  0.0;
      n5[1] =  0.0;
      n5[2] =  1.0;
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
      DUNE_THROW( NotImplemented, "Interpolation for BDM1Cube3D finite elements is not implemented." );

      out.resize(18);
      std::fill(out.begin(), out.end(), 0.0);

      const int qOrder = 4;
      const QuadratureRule<Scalar,1>& rule = QuadratureRules<Scalar,1>::rule(GeometryTypes::cube(1), qOrder);

      for (typename QuadratureRule<Scalar,1>::const_iterator it = rule.begin();
           it != rule.end(); ++it)
      {
        // TODO: write interpolation
      }
    }

  private:
    RangeFieldType sign0, sign1, sign2, sign3, sign4, sign5;
    DomainType n0, n1, n2, n3, n4, n5;
  };

#ifndef DOXYGEN
  template<class D, class R, unsigned int dim>
  class BDMCubeLocalInterpolation<D, R, dim, 0>
  {
    static_assert(AlwaysFalse<D>::value,
                  "`BDMCubeLocalInterpolation` not defined for order 0.");
  };
#endif //#ifndef DOXYGEN

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_CUBE_LOCALINTERPOLATION_HH
