// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALINTERPOLATION_HH

#include <type_traits>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dune
{

  /**
   * \ingroup MonomialImpl
   */
  template<class LB, unsigned int size>
  class MonomialLocalInterpolation
  {
    using D = typename LB::Traits::DomainType;
    using DF = typename LB::Traits::DomainFieldType;
    static const int dimD=LB::Traits::dimDomain;
    using R = typename LB::Traits::RangeType;
    using RF = typename LB::Traits::RangeFieldType;

    using QR = QuadratureRule<DF,dimD>;
    using QRiterator = typename QR::iterator;

  public:
    MonomialLocalInterpolation (const GeometryType &gt_,
                             const LB &lb_)
      : gt(gt_), lb(lb_), Minv(0)
        , qr(QuadratureRules<DF,dimD>::rule(gt, 2*lb.order()))
    {
      // Compute inverse of the mass matrix of the local basis, and store it in Minv
      if(size != lb.size())
        DUNE_THROW(Exception, "size template parameter does not match size of "
                   "local basis");

      const QRiterator qrend = qr.end();
      for(QRiterator qrit = qr.begin(); qrit != qrend; ++qrit) {
        std::vector<R> base;
        lb.evaluateFunction(qrit->position(),base);

        for(unsigned int i = 0; i < size; ++i)
          for(unsigned int j = 0; j < size; ++j)
            Minv[i][j] += qrit->weight() * base[i] * base[j];
      }
      Minv.invert();
    }

    /** \brief Determine coefficients interpolating a given function
     *
     * The method computes the coefficients
     * for the L^2 projection with respect to the given
     * GeometryType. Be careful: the implementation is
     * unstable for higher polynomial degrees.
     */
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      out.clear();
      out.resize(size,C(0.0));

      const QRiterator qrend = qr.end();
      for(QRiterator qrit = qr.begin(); qrit != qrend; ++qrit) {
        //TODO: mass matrix
        auto y = f(qrit->position());
        using RangeType = std::decay_t<decltype(y)>;

        std::vector<R> base;
        lb.evaluateFunction(qrit->position(),base);

        for(unsigned int i = 0; i < size; ++i)
          for(unsigned int j = 0; j < size; ++j)
            if constexpr (std::is_arithmetic<typename RangeType::field_type>::value)
              out[i] += Minv[i][j] * qrit->weight() * y * base[j];
            else
              out[i] += Minv[i][j] * qrit->weight() * y[0] * base[j];
          }
     }

  private:
    GeometryType gt;
    const LB &lb;
    FieldMatrix<RF, size, size> Minv;
    const QR &qr;
  };

}

#endif //DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALINTERPOLATION_HH
