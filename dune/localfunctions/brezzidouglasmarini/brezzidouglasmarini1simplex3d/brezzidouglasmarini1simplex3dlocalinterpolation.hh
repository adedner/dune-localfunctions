// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune
{

  /**
   * \brief First order Brezzi-Douglas-Marini shape functions on the reference
   *        tetrahedron.
   *
   * \tparam LB corresponding LocalBasis giving traits
   *
   * \ingroup BrezziDouglasMariniImpl
   * \nosubgrouping
   */
  template<class LB>
  class BDM1Simplex3DLocalInterpolation
  {

    // f gives v*outer normal at a point on the edge!
    typedef typename LB::Traits::RangeFieldType Scalar;

  public:
    //! \brief Standard constructor
    BDM1Simplex3DLocalInterpolation()
    {
    }

    /**
     * \brief Make set number s, where 0 <= s < 64
     *
     * \param s Edge orientation indicator
     */
    BDM1Simplex3DLocalInterpolation(unsigned int s)
      : BDM1Simplex3DLocalInterpolation()
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

      typedef typename LB::Traits::DomainType DomainType;
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
} // end namespace Dune
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALINTERPOLATION_HH
