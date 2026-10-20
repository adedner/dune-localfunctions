// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_UTILITY_HESSIAN_HH
#define DUNE_LOCALFUNCTIONS_UTILITY_HESSIAN_HH

/**
 * \file This file provides a utility `evaluateHessian()`, which can be used as a
 * fallback implementation for basis functions that do not yet provide this as
 * a member function. It is based on the `partial()` method and computes the entries
 * of the Hessian matrix/tensor componentwise.
 *
 * For scalar basis functions, i.e., range type is a floating point type, the
 * Hessian type is a matrix, e.g., `FieldMatrix<K,dim,dim>` with `dim=dimDomain`.
 * On the other hand, if the range type is a vector type (even with just 1 component)
 * the Hessian type is a 3-tensor. One can use the `Dune::Tensor` data structure
 * for a representation. The used `HessianType` must provide a specialization of
 * `TensorTraits` to extract its rank and individual dimensions.
 */

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/tensortraits.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
namespace Impl {

// get the first component of a vector that is a scalar representation.
template <class K>
const K& getScalarValue (const FieldVector<K,1>& vector)
{
  return vector[0];
}

// return the scalar directly
template <class K>
  requires (Dune::IsNumber<K>::value)
const K& getScalarValue (const K& scalar)
{
  return scalar;
}

/**
 * \brief Evaluation of the Hessian of all basis functions, assuming the range
 * type is a scalar and thus the Hessian is a 2-tensor
 */
template <class LocalBasis, class HessianType,
          class Traits = typename LocalBasis::Traits>
void evaluateHessian2(const LocalBasis& localBasis,
                      const typename Traits::DomainType& x,
                      std::vector<HessianType>& out)
{
  using HessianTraits = TensorTraits<HessianType>;
  static_assert(HessianTraits::rank() == 2);
  static_assert(HessianTraits::rank_dynamic() == 0);
  thread_local std::vector<typename Traits::RangeType> shapePartials;

  out.resize(localBasis.size());
  for (typename HessianTraits::index_type k = 0; k < HessianTraits::static_extent(0); ++k) {
    // diagonal
    std::array<unsigned int,Traits::dimDomain> orders{};
    orders[k] = 2;
    localBasis.partial(orders,x,shapePartials);
    for (std::size_t i = 0; i < out.size(); ++i)
      out[i][k][k] = getScalarValue(shapePartials[i]);

    for (typename HessianTraits::index_type l = 0; l < k; ++l) {
      // off-diagonals
      std::array<unsigned int,Traits::dimDomain> orders{};
      orders[k] = 1; orders[l] = 1;
      localBasis.partial(orders,x,shapePartials);
      for (std::size_t i = 0; i < out.size(); ++i) {
        out[i][std::array{k,l}] = getScalarValue(shapePartials[i]);
        out[i][std::array{l,k}] = getScalarValue(shapePartials[i]);
      }
    }
  }
}

/**
 * \brief Evaluation of the Hessian of all basis functions, assuming the range
 * type is a vector and thus the Hessian is a 3-tensor
 */
template <class LocalBasis, class HessianType,
          class Traits = typename LocalBasis::Traits>
void evaluateHessian3(const LocalBasis& localBasis,
                      const typename Traits::DomainType& x,
                      std::vector<HessianType>& out)
{
  using HessianTraits = TensorTraits<HessianType>;
  using RangeTraits = TensorTraits<typename Traits::RangeType>;
  static_assert(HessianTraits::rank() == 3);
  static_assert(HessianTraits::rank_dynamic() == 0);
  static_assert(RangeTraits::rank_dynamic() == 0);
  static_assert(HessianTraits::static_extent(0) == RangeTraits::static_extent(0));
  thread_local std::vector<typename Traits::RangeType> shapePartials;

  out.resize(localBasis.size());
  for (typename HessianTraits::index_type k = 0; k < HessianTraits::static_extent(1); ++k) {
    // diagonal
    std::array<unsigned int,Traits::dimDomain> orders{};
    orders[k] = 2;
    localBasis.partial(orders,x,shapePartials);
    for (std::size_t i = 0; i < out.size(); ++i)
      for (typename HessianTraits::index_type j = 0; j < HessianTraits::static_extent(0); ++j)
        out[i][std::array{j,k,k}] = shapePartials[i][j];

    for (typename HessianTraits::index_type l = 0; l < k; ++l) {
      // off-diagonals
      std::array<unsigned int,Traits::dimDomain> orders{};
      orders[k] = 1; orders[l] = 1;
      localBasis.partial(orders,x,shapePartials);
      for (std::size_t i = 0; i < out.size(); ++i) {
        for (typename HessianTraits::index_type j = 0; j < HessianTraits::static_extent(0); ++j) {
          out[i][std::array{j,k,l}] = shapePartials[i][j];
          out[i][std::array{j,l,k}] = shapePartials[i][j];
        }
      }
    }
  }
}

} // end namespace Impl

/**
 * \brief Evaluation of the Hessian of all basis functions in a local basis.
 */
template <class LocalBasis, class HessianType,
          class Traits = typename LocalBasis::Traits>
void evaluateHessian (const LocalBasis& localBasis,
                      const typename Traits::DomainType& x,
                      std::vector<HessianType>& out)
{
  using HessianTraits = TensorTraits<HessianType>;
  if constexpr (HessianTraits::rank() == 2)
    Impl::evaluateHessian2(localBasis, x, out);
  else if constexpr (HessianTraits::rank() == 3)
    Impl::evaluateHessian3(localBasis, x, out);
  else
    DUNE_THROW(Dune::NotImplemented, "Hessian with higher rank not implemented");
}

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_UTILITY_HESSIAN_HH
