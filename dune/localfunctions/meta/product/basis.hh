// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_META_PRODUCT_BASIS_HH
#define DUNE_LOCALFUNCTIONS_META_PRODUCT_BASIS_HH

#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/localfunctions/common/localbasis.hh>

namespace Dune {

template<class LB1, class LB2, class M>
class PrismaticProductLocalBasis
{
  static_assert(LB2::Traits::dimRange == 1,
                "PrismaticProductLocalBasis works only with 1-dimensional range of the second factor");

  using Traits1 = typename LB1::Traits;
  using Traits2 = typename LB2::Traits;
  using Mapping = M;

  using DF = std::common_type_t<typename Traits1::DomainFieldType, typename Traits2::DomainFieldType>;
  static constexpr int dimDomain = LB1::Traits::dimDomain + LB2::Traits::dimDomain;
  using Domain = FieldVector<DF,dimDomain>;

  using RF = std::common_type_t<typename Traits1::RangeFieldType, typename Traits2::RangeFieldType>;
  static constexpr int dimRange = LB1::Traits::dimRange * LB2::Traits::dimRange;
  using Range = FieldVector<RF,dimRange>;
  using Jacobian = FieldMatrix<RF,dimRange,dimDomain>;

private:
  const LB1* lb1_;
  const LB2* lb2_;
  const Mapping* m_;

  // caches
  mutable std::vector<typename Traits1::RangeType> lb1Values_;
  mutable std::vector<typename Traits2::RangeType> lb2Values_;

private:
  static constexpr typename Traits1::DomainType domain1 (const Domain& x)
  {
    return Dune::unpackIntegerSequence([&](auto... i) {
      return typename Traits1::DomainType{x[i]...};
    }, std::make_index_sequence<Traits1::dimDomain>{});
  }

  static constexpr typename Traits2::DomainType domain2 (const Domain& x)
  {
    return Dune::unpackIntegerSequence([&](auto... i) {
      return typename Traits2::DomainType{x[(Traits1::dimDomain+i)]...};
    }, std::make_index_sequence<Traits2::dimDomain>{});
  }

public:
  //! Traits type to characterize the local basis
  using Traits = LocalBasisTraits<DF,dimDomain,Domain,RF,dimRange,Range,Jacobian>;

  //! Construct a PrismaticProductLocalBasis
  /**
   * \param lb1 First local basis object to construct this object from.
   * \param lb2 Second local basis object to construct this object from.
   *
   * The local basis parameters hold references to the local basis objects.
   * These reference are also copied when this object is copied.
   */
  PrismaticProductLocalBasis (const LB1& lb1, const LB2& lb2, const Mapping& mapping)
    : lb1_(&lb1)
    , lb2_(&lb2)
    , m_(&mapping)
  {}

  //! Number of shape functions
  std::size_t size () const noexcept { return m_->required_span_size(); }

  //! Polynomial order of the shape functions for quadrature
  std::size_t order () const noexcept { return std::max<std::size_t>(lb1_->order(), lb2_->order()); }

  //! Evaluate all shape functions at given position
  void evaluateFunction (const Domain& x, std::vector<Range>& out) const
  {
    lb1_->evaluateFunction(domain1(x), lb1Values_);
    lb2_->evaluateFunction(domain2(x), lb2Values_);

    out.resize(size());
    for (std::size_t i = 0; i < lb1Values_.size(); ++i)
      for (std::size_t j = 0; j < lb2Values_.size(); ++j)
        out[(*m_)(i,j)] = lb1Values_[i] * lb2Values_[j];
  }

  //! Evaluate Jacobian of all shape functions at given position
  void evaluateJacobian (const Domain& x, std::vector<Jacobian>& out) const
  {
    thread_local std::vector<typename Traits1::JacobianType> lb1Jacobians;
    thread_local std::vector<typename Traits2::JacobianType> lb2Jacobians;

    auto x1 = domain1(x);
    auto x2 = domain2(x);
    lb1_->evaluateFunction(x1, lb1Values_);
    lb1_->evaluateJacobian(x1, lb1Jacobians);
    lb2_->evaluateFunction(x2, lb2Values_);
    lb2_->evaluateJacobian(x2, lb2Jacobians);

    out.resize(size());
    for (std::size_t i = 0; i < lb1Values_.size(); ++i) {
      for (std::size_t j = 0; j < lb2Values_.size(); ++j) {
        for (int k0 = 0; k0 < Traits1::dimRange; ++k0) {
          for (int k1 = 0; k1 < Traits1::dimDomain; ++k1)
            out[(*m_)(i,j)][k0][k1] = lb1Jacobians[i][k0][k1] * lb2Values_[j][0];
          for (int k1 = 0; k1 < Traits2::dimDomain; ++k1)
            out[(*m_)(i,j)][k0][Traits1::dimDomain + k1] = lb1Values_[i][k0] * lb2Jacobians[j][0][k1];
        }
      }
    }
  }

  //! \brief Evaluate partial derivatives of all shape functions
  void partial (const std::array<unsigned int, dimDomain>& order,
                const Domain& x, std::vector<Range>& out) const
  {
    unsigned int totalOrder = 0;
    for (int i = 0; i < dimDomain; ++i)
      totalOrder += order[i];

    switch (totalOrder) {
    case 0:
      evaluateFunction(x,out);
      break;
    default: {
      DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }
  }
};

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_PRODUCT_BASIS_HH
