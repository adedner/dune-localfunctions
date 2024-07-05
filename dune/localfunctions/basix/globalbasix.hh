// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BASIX_GLOBALBASIX_HH
#define DUNE_LOCALFUNCTIONS_BASIX_GLOBALBASIX_HH

#if HAVE_BASIX

#include <utility>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/basix/localbasix.hh>
#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>

namespace Dune {

/**
 * \brief Definition of a global finite-element based on the basix library.
 *
 * \tparam G  The geometry of the global element the finite-element is defined on.
 * \tparam F  A floating-point type used for the domain and range types of the basis functions.
 * \tparam dimDomain The dimension of the local domain the basis functions are defined on.
 * \tparam dimRange  The dimension of the range of the basis functions. (TODO: Needs to be clarified)
 */
template <class G, class F, int dimDomain, int dimRange>
class BasixFiniteElement
{
public:
  using Geometry = G;
  using Basix = basix::FiniteElement<F>;
  using LocalFiniteElement = BasixLocalFiniteElement<F,dimDomain,dimRange>;

  struct Basis
  {
    using LocalBasis = typename LocalFiniteElement::Traits::LocalBasis;
    using Traits = LocalBasisTraits<
      F,dimDomain,FieldVector<F,dimDomain>, // domain
      F,dimRange,FieldVector<F,dimRange>,   // range
      FieldMatrix<F,dimRange,dimDomain>     // jacobian
      >;

    /// \brief Return the number of basis functions
    std::size_t size () const
    {
      return lb_->size();
    }

    /// \brief Degree of the minimal polynomial space all basis functions are embedded
    std::size_t order () const
    {
      return lb_->order();
    }

    /// \brief Evaluate all shape functions in a point x
    void evaluateFunction(const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
    {
      lb_->evaluateFunction(x, out);
      // TODO: might need a transformation
    }

    /// \brief Evaluate all shape function jacobians in a point x
    void evaluateJacobian(const typename Traits::DomainType& x,
                          std::vector<typename Traits::JacobianType>& out) const
    {
      lb_->evaluateJacobian(x,out);
      // TODO: might need a transformation
    }

    /// \brief Evaluate all shape function partial derivatives with given orders in a point x
    void partial(const std::array<unsigned int,dimDomain>& order,
                  const typename Traits::DomainType& x,
                  std::vector<typename Traits::RangeType>& out) const
    {
      lb_->partial(order,x,out);
      // TODO: might need a transformation
    }

    LocalBasis* lb_;
  };


  struct Coefficients
  {
    using LocalCoefficients = typename LocalFiniteElement::Traits::LocalCoefficients;

    /// \brief Return the number of local keys associated to local basis functions.
    std::size_t size () const
    {
      return lc_->size();
    }

    /// \brief Obtain the LocalKey associated to the `i`th basis function.
    const LocalKey& localKey (std::size_t i) const
    {
      return lc_->localKey(i);
    }

    LocalCoefficients* lc_;
  };


  struct Interpolation
  {
    using LocalInterpolation = typename LocalFiniteElement::Traits::LocalInterpolation;

    /// \brief Determine coefficients interpolating a given function `f`
    /// and store them in the output vector `out`.
    template<class Func, class C>
    void interpolate (const Func& f, std::vector<C>& out) const
    {
      li_->interpolate(f,out);
    }

    LocalInterpolation* li_;
  };


  struct Traits {
    using Basis = typename BasixFiniteElement::Basis;
    using Coefficients = typename BasixFiniteElement::Coefficients;
    using Interpolation = typename BasixFiniteElement::Interpolation;
  };

public:
  /// \brief Construct a global finite-element from an associated local finite-element
  explicit BasixFiniteElement (const LocalFiniteElement& lfe)
    : lfe_(lfe)
    , basis_{&lfe_.localBasis()}
    , coefficients_{&lfe_.localCoefficients()}
    , interpolation_{&lfe_.localInterpolation()}
  {}

  /// \brief Construct a global finite-element from an associated local finite-element
  explicit BasixFiniteElement (LocalFiniteElement&& lfe)
    : lfe_(std::move(lfe))
    , basis_{&lfe_.localBasis()}
    , coefficients_{&lfe_.localCoefficients()}
    , interpolation_{&lfe_.localInterpolation()}
  {}

  /// \brief Construct the local finite-element from the basix library
  explicit BasixFiniteElement (const Basix& basix)
    : BasixFiniteElement(LocalFiniteElement(basix))
  {}

  /// \brief Construct the local finite-element from the basix library.
  explicit BasixFiniteElement (Basix&& basix)
    : BasixFiniteElement(LocalFiniteElement(std::move(basix)))
  {}

  /// \brief Move constructor, needs to re-assign the internal pointers.
  BasixFiniteElement (const BasixFiniteElement& other)
    : BasixFiniteElement(other.lfe_)
  {}

  /// \brief Move constructor, needs to re-assign the internal pointers.
  BasixFiniteElement (BasixFiniteElement&& other)
    : BasixFiniteElement(std::move(other.lfe_))
  {}

  /// \brief Bind the global finite-element to a global element geometry and
  /// a cell permutation information.
  void bind (const Geometry& geometry, std::uint32_t cellInfo = 0)
  {
    geometry_.emplace(geometry);
    cellInfo_ = cellInfo;
  }

  /// \brief Obtain a reference to the basis.
  const Basis& basis () const { return basis_; }

  /// \brief Obtain a reference to the coefficients.
  const Coefficients& coefficients () const { return coefficients_; }

  /// \brief Obtain a reference to the interpolation.
  const Interpolation& interpolation () const { return  interpolation_; }

  /// \brief Return the dimension of the finite-element
  std::size_t size () const
  {
    return fe_.dim();
  }

  /// \brief Return the GeometryType the local finite-element is defined on
  GeometryType type () const
  {
    assert(lfe_->type() == geometry_->type());
    return lfe_->type();
  }

  /// \brief Obtain a reference to the basix implementation
  const Basix& basix () const
  {
    return lfe_->basix();
  }

private:
  const LFE* lfe_ = nullptr;
  std::optional<Geometry> geometry_ = std::nullopt;
  std::uint32_t cellInfo_ = 0;

  Basis basis_;
  Coefficients coefficients_;
  Interpolation interpolation_;
};

} // end namespace Dune

#endif // HAVE_BASIX
#endif // DUNE_LOCALFUNCTIONS_BASIX_GLOBALBASIX_HH
