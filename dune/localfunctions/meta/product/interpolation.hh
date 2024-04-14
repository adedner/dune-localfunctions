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
   * @brief A local L2 interpolation taking a test basis and a quadrature
   *        rule.
   *
   * This class computes a local interpolation where the coefficients
   * are of the form:
   *     c = M^{-1}b
   * - M is the mass matrix with respect to the given basis and
   * - b = int f phi (where phi are the basis functions).
   * Thus the resulting local function u=c.varphi is defined through
   * the l2 interpolation int u phi = in f phi for all phi in the
   * base function set.
   **/
  template< class B, class Q >
  class PrismaticProductLocalInterpolation
  {
    using LocalBasis = B;
    using Quadrature = Q;

    using DomainType = typename LocalBasis::Traits::DomainType;
    using RangeType = typename LocalBasis::Traits::RangeType;
    using RangeFieldType = typename LocalBasis::Traits::RangeFieldType;

  public:
    PrismaticProductLocalInterpolation (const LocalBasis& basis, const Quadrature& quadrature)
      : basis_(basis)
      , quadrature_(quadrature)
      , massMatrix_(basis_.size(), basis_.size(), RangeFieldType(0))
    {
      std::vector<RangeType> basisValues( basis_.size() );
      for( auto& qp : quadrature_) {
        basis_.evaluateFunction( qp.position(), basisValues );
        for (unsigned int i=0; i<basis_.size(); ++i)
          for (unsigned int j=0; j<basis_.size(); ++j)
            massMatrix_[i][j] += (basisValues[i]*basisValues[j]) * qp.weight();
      }
      massMatrix_.invert();
    }

    /** \brief Interpolate a function that implements Range operator()(Domain) */
    template <class Function, class C>
    void interpolate (const Function& f, std::vector<C>& coefficients) const
    {
      using R = std::decay_t<std::invoke_result_t<Function, DomainType>>;
      const unsigned int size = basis_.size();
      thread_local std::vector<R> rhs( size );
      thread_local std::vector<RangeType> basisValues( size );

      rhs.resize( size );
      basisValues.resize( size );
      for (unsigned int i = 0; i < size; ++i)
        rhs[ i ] = C(0);

      for (auto& qp : quadrature_) {
        basis_.evaluateFunction( qp.position(), basisValues );
        auto factor = f( qp.position() );
        for (unsigned int i = 0; i < size; ++i)
          rhs[i] += factor * basisValues[i] * qp.weight();
      }

      coefficients.resize( size );
      for (unsigned int i=0; i<size; ++i) {
        coefficients[i] = 0;
        for (unsigned int j=0; j<size; ++j) {
          coefficients[i] += massMatrix_[i][j] * rhs[j];
        }
      }

    }

  private:
    const LocalBasis& basis_;
    const Quadrature& quadrature_;

    using MassMatrix = DynamicMatrix<RangeFieldType>;
    MassMatrix massMatrix_;
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_PRODUCT_INTERPOLATION_HH
