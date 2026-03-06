// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX2D_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX2D_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/geometry/quadraturerules.hh>

namespace Dune
{
  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief First order Brezzi-Douglas-Marini shape functions on the reference triangle.
   *
   * \tparam LB corresponding LocalBasis giving traits
   *
   * \nosubgrouping
   */
  template<class LB>
  class BDM1Simplex2DLocalInterpolation
  {

  public:
    //! \brief Standard constructor
    BDM1Simplex2DLocalInterpolation ()
    {
      sign0 = sign1 = sign2 = 1.0;
    }

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    BDM1Simplex2DLocalInterpolation (unsigned int s)
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
      // f gives v*outer normal at a point on the edge!
      typedef typename LB::Traits::RangeFieldType Scalar;

      out.resize(6);
      std::fill(out.begin(), out.end(), 0.0);

      const int qOrder = 4;
      const Dune::QuadratureRule<Scalar,1>& rule = Dune::QuadratureRules<Scalar,1>::rule(Dune::GeometryTypes::simplex(1), qOrder);

      for (typename Dune::QuadratureRule<Scalar,1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        const Scalar qPos = it->position();
        const Scalar weight_0 = it->weight();
        const Scalar weight_1 = (2.0*qPos - 1.0)*it->weight();

        typename LB::Traits::DomainType localPos;

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
    typename LB::Traits::RangeFieldType sign0,sign1,sign2;
  };
}

#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX2D_LOCALINTERPOLATION_HH
