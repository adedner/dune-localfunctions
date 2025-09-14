// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_META_PRODUCT_INTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_META_PRODUCT_INTERPOLATION_HH

#include <type_traits>
#include <vector>

#include <dune/common/dynmatrix.hh>

namespace Dune
{
  /**
   * \brief A local L2 interpolation taking a test basis and a quadrature rule.
   *
   * This class computes a local interpolation based on an L2-projection, where the coefficients
   * are of the form \f$c = M^{-1}b\f$ with
   * - \f$M = M_{ij} = \int_T \phi_i \phi_j\f$ is the mass matrix with respect to the given basis,
   * - \f$b = b_{i} = \int_T f phi_i\f$ (where \f$\phi_i\f$ are the basis functions), and
   * - \f$T\f$ the reference element domain.
   **/
  template <class B, class Q1, class Q2>
  class PrismaticProductLocalInterpolation
  {
    using LocalBasis = B;
    using Quad1 = Q1;
    using Quad2 = Q2;

    using DomainType = typename LocalBasis::Traits::DomainType;
    using RangeType = typename LocalBasis::Traits::RangeType;
    using RangeFieldType = typename LocalBasis::Traits::RangeFieldType;

  public:
    PrismaticProductLocalInterpolation (const LocalBasis& basis, const Quad1& quad1, const Quad2& quad2)
      : basis_(basis)
      , quad1_(quad1)
      , quad2_(quad2)
      , massMatrix_(basis_.size(), basis_.size(), RangeFieldType(0))
    {
      std::vector<RangeType> basisValues(basis_.size());
      DomainType pos;
      for (auto const& qp1 : quad1_) {
        for (int k = 0; k < Quad1::d; ++k)
          pos[k] = qp1.position()[k];
        for (auto const& qp2 : quad2_) {
          for (int k = 0; k < Quad2::d; ++k)
            pos[Quad1::d + k] = qp2.position()[k];
          basis_.evaluateFunction(pos, basisValues);
          for (std::size_t i = 0; i < basis_.size(); ++i)
            for (std::size_t j = 0; j < basis_.size(); ++j)
              massMatrix_[i][j] += (basisValues[i] * basisValues[j]) * qp1.weight() * qp2.weight();
        }
      }
      massMatrix_.invert();
    }

    //! Interpolate a function that implements Range operator()(Domain)
    template <class Function, class C>
    void interpolate (const Function& f, std::vector<C>& coefficients) const
    {
      using R = std::decay_t<std::invoke_result_t<Function, DomainType>>;
      std::size_t size = basis_.size();
      thread_local std::vector<R> rhs(size);
      thread_local std::vector<RangeType> basisValues(size);

      rhs.resize(size);
      for (std::size_t i = 0; i < size; ++i)
        rhs[i] = 0;

      DomainType pos;
      for (auto const& qp1 : quad1_) {
        for (int k = 0; k < Quad1::d; ++k)
          pos[k] = qp1.position()[k];
        for (auto const& qp2 : quad2_) {
          for (int k = 0; k < Quad2::d; ++k)
            pos[Quad1::d + k] = qp2.position()[k];
          basis_.evaluateFunction(pos, basisValues);
          auto factor = f(pos);
          for (std::size_t i = 0; i < size; ++i)
            rhs[i] += factor * basisValues[i] * qp1.weight() * qp2.weight();
        }
      }

      // compute the matrix vector product massMatrix^{-1} * rhs
      coefficients.resize(size);
      for (std::size_t i = 0; i < size; ++i) {
        coefficients[i] = 0;
        for (std::size_t j = 0; j < size; ++j) {
          coefficients[i] += massMatrix_[i][j] * rhs[j];
        }
      }

    }

  private:
    const LocalBasis& basis_;
    const Quad1& quad1_;
    const Quad2& quad2_;

    using MassMatrix = DynamicMatrix<RangeFieldType>;
    MassMatrix massMatrix_;
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_PRODUCT_INTERPOLATION_HH
