// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BASIX_LOCALBASIX_HH
#define DUNE_LOCALFUNCTIONS_BASIX_LOCALBASIX_HH

#if HAVE_BASIX

#include <utility>
#include <vector>

#include <basix/finite-element.h>
#include <basix/indexing.h>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/std/mdarray.hh>
#include <dune/common/std/mdspan.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>

namespace Dune {
namespace Impl {

  // Map the cell::type from basix to the Dune GeometryType
  GeometryType basix2dune (basix::cell::type cell_type)
  {
    switch (cell_type) {
      case basix::cell::type::point:          return GeometryTypes::vertex;
      case basix::cell::type::interval:       return GeometryTypes::line;
      case basix::cell::type::triangle:       return GeometryTypes::triangle;
      case basix::cell::type::tetrahedron:    return GeometryTypes::tetrahedron;
      case basix::cell::type::quadrilateral:  return GeometryTypes::quadrilateral;
      case basix::cell::type::hexahedron:     return GeometryTypes::hexahedron;
      case basix::cell::type::prism:          return GeometryTypes::prism;
      case basix::cell::type::pyramid:        return GeometryTypes::pyramid;
      default: return GeometryTypes::none(0);
    }
  }

  // Map the derivative-order tuple from the partial() method to the
  // derivative index used in basix
  template <std::size_t dim>
  int indexing (const std::array<unsigned int,dim>& orders)
  {
    if constexpr (dim == 1)
      return basix::indexing::idx(orders[0]);
    else if constexpr (dim == 2)
      return basix::indexing::idx(orders[0],orders[1]);
    else if constexpr (dim == 3)
      return basix::indexing::idx(orders[0],orders[1],orders[2]);
    else
      return 0;
  }

  // Map the first order derivative in direction d to the derivative index
  // used in basix
  template <std::size_t dim>
  int indexing (unsigned int d)
  {
    std::array<unsigned int,dim> orders;
    orders[d] = 1;
    return indexing(orders);
  }

} // end namespace Impl

/// \brief Type of the basis function range
enum class RangeClass
{
  scalar = 0,
  vector = 1,
  matrix = 2,
};


/**
 * \brief Definition of a local finite-element based on the basix library.
 *
 * \tparam dimDomain The dimension of the local domain the basis functions are defined on.
 * \tparam rangeClass  The class of the basis function range, i.e., RangeClass:scalar, vector or matrix.
 * \tparam F  A floating-point type used for the domain and range types of the basis functions.
 */
template <int dimDomain, RangeClass rangeClass, class F = double>
class BasixLocalFiniteElement
{
public:
  using FieldType = F;
  using Basix = basix::FiniteElement<FieldType>;

  struct LocalBasis
  {
    static constexpr std::size_t dimRange ()
    {
      switch (rangeClass) {
        case RangeClass::scalar: return 1;
        case RangeClass::vector: return dimDomain;
        case RangeClass::matrix: return dimDomain * dimDomain;
        default: return 0;
      }
    }

    using Traits = LocalBasisTraits<
      F,dimDomain,FieldVector<F,dimDomain>,   // domain
      F,dimRange(),FieldVector<F,dimRange()>, // range
      FieldMatrix<F,dimRange(),dimDomain>     // jacobian
      >;

    using Domain = typename Traits::DomainType;
    using Range = typename Traits::RangeType;
    using RangeSpan = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,dimRange()>>;
    using Jacobian = typename Traits::JacobianType;
    using JacobianSpan = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,dimRange(),dimDomain>>;

    /// \brief Return the number of basis functions
    std::size_t size () const
    {
      return basix_->dim();
    }

    /// \brief Degree of the minimal polynomial space all basis functions are embedded
    std::size_t order () const
    {
      // TODO: Is this the right degree, or to we need the embedded_subdegree()?
      return basix_->embedded_superdegree();
    }

    /// \brief Evaluate all shape functions in a point x
    void evaluateFunction (const Domain& x, RangeSpan out) const
    {
      using TabulateX = basix::element::mdspan_t<const F,2>;
      TabulateX _x{x.data(), std::array<std::size_t,2>{1,dimDomain}};

      // Add two more dimensions to the RangeSpan:
      // 1. number of derivative components (== 1, since only derivative-order 0)
      // 2. number of points (== 1)
      using TabulateOut = basix::element::mdspan_t<F,4>;
      TabulateOut _out{out.data_handle(), std::array<std::size_t,4>{1,1,out.extent(0),out.extent(1)}};
      basix_->tabulate(0, _x, _out);
    }

    void evaluateFunction (const Domain& x, std::vector<Range>& out) const
    {
      out.resize(size());

#if USE_OUT_VECTOR_FOR_MDSPAN
      RangeSpan _out{&out[0][0], size()};
      evaluateFunction(x,_out);
#else
      evaluationBuffer_.resize(size() * dimRange());
      RangeSpan _out{evaluationBuffer_.data(), size()};
      evaluateFunction(x,_out);

      for (std::size_t i = 0; i < size(); ++i)
        for (std::size_t j = 0; j < dimRange(); ++j)
          out[i][j] = _out(i,j);
#endif
    }

    /// \brief Evaluate all shape function Jacobians in a point x
    void evaluateJacobian (const Domain& x, JacobianSpan out) const
    {
      using TabulateX = basix::element::mdspan_t<const F,2>;
      TabulateX _x{x.data(), std::array<std::size_t,2>{1,dimDomain}};

      // Define an extended JacobianSpan, since we need to store also the function evaluation
      // 1. number of derivative components (== dimDomain+1, since function evaluation and all first-order derivatives are computed)
      // 2. number of points (== 1)
      using TabulateOut = basix::element::mdspan_t<F,4>;
      evaluationBuffer_.resize((dimDomain+1) * size() * dimRange());
      TabulateOut _out{evaluationBuffer_.data(), std::array<std::size_t,4>{dimDomain+1,1,size(),dimRange()}};
      basix_->tabulate(1, _x, _out);

      for (std::size_t i = 0; i < size(); ++i)
        for (std::size_t j = 0; j < dimRange(); ++j)
          for (std::size_t k = 0; k < dimDomain; ++k)
            out(i,j,k) = _out(k+1,0,i,j);
    }

    /// \brief Evaluate all shape function Jacobians in a point x
    void evaluateJacobian (const Domain& x, std::vector<Jacobian>& out) const
    {
      out.resize(size());

#if USE_OUT_VECTOR_FOR_MDSPAN
      JacobianSpan _out{&out[0][0][0], size()};
      evaluateJacobian(x, _out);
#else
      thread_local std::vector<F> jacobianEvaluationBuffer;
      jacobianEvaluationBuffer.resize(size() * dimRange() * dimDomain);
      JacobianSpan _out{jacobianEvaluationBuffer.data(), size()};
      evaluateJacobian(x, _out);

      for (std::size_t i = 0; i < size(); ++i)
        for (std::size_t j = 0; j < dimRange(); ++j)
          for (std::size_t k = 0; k < dimDomain; ++k)
            out[i][j][k] = _out(i,j,k);
#endif
    }

    /// \brief Evaluate all shape function partial derivatives with given orders in a point x
    void partial (const std::array<unsigned int,dimDomain>& order,
                  const Domain& x, RangeSpan out) const
    {
      // Add another dimension to the DomainSpan:
      // - number of points (== 1)
      using TabulateX = basix::element::mdspan_t<const F,2>;
      TabulateX _x{x.data(), std::array<std::size_t,2>{1,dimDomain}};

      int totalOrder = std::accumulate(order.begin(), order.end(), 0);
      auto shape = basix_->tabulate_shape(totalOrder, 1);

      // Add two more dimensions to the RangeSpan:
      // 1. number of derivative components (== 1, since only derivative-order 0)
      // 2. number of points (== 1)
      using TabulateOut = basix::element::mdspan_t<F,4>;
      TabulateOut _out{out.data_handle(),std::array<std::size_t,4>{shape[0],1,out.extent(0),out.extent(1)}};
      basix_->tabulate(totalOrder, _x, _out);

      for (std::size_t i = 0; i < size(); ++i)
        for (std::size_t j = 0; j < dimRange(); ++j)
          out(i,j) = _out(Impl::indexing(order),0,i,j);
    }

    /// \brief Evaluate all shape function partial derivatives with given orders in a point x
    void partial(const std::array<unsigned int,dimDomain>& order,
                  const Domain& x, std::vector<Range>& out) const
    {
      out.resize(size());

#if USE_OUT_VECTOR_FOR_MDSPAN
      RangeSpan _out{&out[0][0], size()};
      partial(order, x, _out);
#else
      thread_local std::vector<F> partialEvaluationBuffer;
      partialEvaluationBuffer.resize(size() * dimRange());
      RangeSpan _out{partialEvaluationBuffer.data(), size()};
      partial(order, x, _out);

      for (std::size_t i = 0; i < size(); ++i)
        for (std::size_t j = 0; j < dimRange(); ++j)
            out[i][j] = _out(i,j);
#endif
    }

    const Basix& basix () const { return *basix_; }

    Basix* basix_;
    mutable std::vector<F> evaluationBuffer_ = {};
  };


  struct LocalCoefficients
  {
    LocalCoefficients (Basix* basix)
      : basix_(basix)
      , localKeys_(basix_->dim())
    {
      // map from entity_dofs into LocalKeys
      auto& entity_dofs = basix_->entity_dofs();
      for (std::size_t d = 0; d < entity_dofs.size(); ++d)
        for (std::size_t s = 0; s < entity_dofs[d].size(); ++s)
          for (std::size_t i = 0; i < entity_dofs[d][s].size(); ++i)
            localKeys_[entity_dofs[d][s][i]] = LocalKey(s,Impl::basix2dune(basix_->cell_type()).dim()-d, i);
            // TODO: We might need an index mapping from the basix
            // reference element numbering to the Dune numbering
            // in the parameter `s`.
    }

    /// \brief Return the number of local keys associated to local basis functions.
    std::size_t size () const
    {
      return localKeys_.size();
    }

    /// \brief Obtain the LocalKey associated to the `i`th basis function.
    const LocalKey& localKey (std::size_t i) const
    {
      return localKeys_[i];
    }

  private:
    Basix* basix_;
    std::vector<LocalKey> localKeys_;
  };


  struct LocalInterpolation
  {
    /// \brief Determine coefficients interpolating a given function `f`
    /// and store them in the output vector `out`.
    template<class Func, class C>
    void interpolate (const Func& f, std::vector<C>& out) const
    {
      auto& p = basix_->points();
      using Points = Std::mdspan<const F, Std::extents<int, Std::dynamic_extent, dimDomain>>;
      Points points{p.first.data(), p.second[0]};

      // TODO: What is the proper type for the values-vector? It depends on the type of
      // the interpolation matrix.
      Std::mdarray<F,Std::dextents<int,1>> values{points.extent(0)};
      for (int i = 0; i < points.extent(0); ++i)
      {
        FieldVector<F,dimDomain> x;
        for (int j = 0; j < dimDomain; ++j)
          x[j] = points(i,j);

        values(i) = f(x);
        // TODO: we need to handle non-scalar return types
      }

      auto& m = basix_->interpolation_matrix();
      using InterpolationMatrix = Std::mdspan<const F, Std::dextents<int, 2>>;
      InterpolationMatrix interpolationMatrix{m.first.data(), m.second[0], m.second[1]};

      // TODO: It would be nice to have .mv function on the interpolationMatrix.
      out.resize(interpolationMatrix.extent(0));
      for (int i = 0; i < interpolationMatrix.extent(0); ++i) {
        out[i] = F(0);
        for (int j = 0; j < interpolationMatrix.extent(1); ++j)
          out[i] += interpolationMatrix(i,j) * values(j);
      }
    }

    Basix* basix_;
  };


  using Traits = LocalFiniteElementTraits<
    LocalBasis,
    LocalCoefficients,
    LocalInterpolation>;

public:
  /// \brief Construct a local finite-element from the Basix definition
  explicit BasixLocalFiniteElement (const Basix& basix)
    : basix_(basix)
    , localBasis_{&basix_}
    , localCoefficients_{&basix_}
    , localInterpolation_{&basix_}
  {}

  /// \brief Construct a local finite-element from the Basix definition
  explicit BasixLocalFiniteElement (Basix&& basix)
    : basix_(std::move(basix))
    , localBasis_{&basix_}
    , localCoefficients_{&basix_}
    , localInterpolation_{&basix_}
  {}

  /// \brief Copy constructor, needs to re-assign the internal pointers
  BasixLocalFiniteElement (const BasixLocalFiniteElement& other)
    : BasixLocalFiniteElement(other.basix_)
  {}

  /// \brief Move constructor, needs to re-assign the internal pointers
  BasixLocalFiniteElement (BasixLocalFiniteElement&& other)
    : BasixLocalFiniteElement(std::move(other.basix_))
  {}

  /// \brief Obtain a reference to the local basis
  const LocalBasis& localBasis () const { return localBasis_; }

  /// \brief Obtain a reference to the local coefficients
  const LocalCoefficients& localCoefficients () const { return localCoefficients_; }

  /// \brief Obtain a reference to the local interpolation
  const LocalInterpolation& localInterpolation () const { return  localInterpolation_; }

  /// \brief Return the dimension of the local finite-element
  std::size_t size () const
  {
    return basix_.dim();
  }

  /// \brief Return the GeometryType the local finite-element is defined on
  GeometryType type () const
  {
    return Impl::basix2dune(basix_.cell_type());
  }

  /// \brief Obtain a reference to the basix implementation
  const Basix& basix () const
  {
    return basix_;
  }

private:
  Basix basix_;

  LocalBasis localBasis_;
  LocalCoefficients localCoefficients_;
  LocalInterpolation localInterpolation_;
};

} // end namespace Dune

#endif // HAVE_BASIX
#endif // DUNE_LOCALFUNCTIONS_BASIX_LOCALBASIX_HH
