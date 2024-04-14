// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_META_PRODUCT_HH
#define DUNE_LOCALFUNCTIONS_META_PRODUCT_HH

#include <cstddef>
#include <memory>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/meta/product/basis.hh>
#include <dune/localfunctions/meta/product/coefficients.hh>
#include <dune/localfunctions/meta/product/interpolation.hh>

namespace Dune {

  //! \brief Meta-finite element turning two finite elements into
  //!        a tensor-product of these two finite elements
  /**
   * \ingroup LocalFunctions
   *
   * \tparam LFE1  Type of the first local finite-element
   * \tparam LFE2  Type of the second local finite-element
   */
  template<class LFE1, class LFE2>
  class PrismaticProduct
  {
    using LB1 = typename LFE1::Traits::LocalBasisType;
    using LB2 = typename LFE2::Traits::LocalBasisType;

    using LQ1 = QuadratureRules<typename LB1::Traits::DomainFieldType, LB1::Traits::dimDomain>;
    using LQ2 = QuadratureRules<typename LB1::Traits::DomainFieldType, LB1::Traits::dimDomain>;

    using LI1 = typename LFE1::Traits::LocalInterpolationType;
    using LI2 = typename LFE2::Traits::LocalInterpolationType;

  public:
    //! types of component objects
    struct Traits {
      //! type of the Basis
      using LocalBasisType = PrismaticProductLocalBasis<LB1,LB2>;
      //! type of the Coefficients
      using LocalCoefficientsType = PrismaticProductLocalCoefficients;
      //! type of the quadrature rules to be used in L2-interpolation
      using LocalQuadratureType = TensorProductQuadratureRule<typename LocalBasisType::Traits::DomainFieldType, LocalBasisType::Traits::dimDomain>;
      //! type of the Interpolation
      using LocalInterpolationType = PrismaticProductLocalInterpolation<LocalBasisType,LocalQuadratureType>;
    };

  private:
    std::shared_ptr<const LFE1> lfe1_;
    std::shared_ptr<const LFE2> lfe2_;

    typename Traits::LocalBasisType lb_;
    typename Traits::LocalCoefficientsType lc_;
    typename Traits::LocalQuadratureType lq_;
    typename Traits::LocalInterpolationType li_;

  public:
    //! Construct a finite element
    PrismaticProduct (const LFE1& lfe1, const LFE2& lfe2)
      : PrismaticProduct(std::make_shared<const LFE1>(lfe1), std::make_shared<const LFE2>(lfe2))
    {}

    //! Construct a finite element
    PrismaticProduct (std::shared_ptr<const LFE1> lfe1, std::shared_ptr<const LFE2> lfe2)
      : lfe1_(std::move(lfe1))
      , lfe2_(std::move(lfe2))
      , lb_(lfe1_->localBasis(), lfe2_->localBasis())
      , lc_(lfe1_->localCoefficients(), lfe2_->localCoefficients(),
          referenceElement<typename LB1::Traits::DomainFieldType,LB1::Traits::dimDomain>(lfe1_->type()))
      , lq_(GeometryTypes::prismaticProduct(lfe1_->type(), lfe2_->type()),
          LQ1::rule(lfe1_->type(), lfe1_->localBasis().order()),
          LQ2::rule(lfe2_->type(), lfe2_->localBasis().order()))
      , li_(lb_, lq_)
    {}

    //! Extract basis of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return lb_;
    }

    //! Extract coefficients of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const typename Traits::LocalCoefficientsType& localCoefficients() const
    {
      return lc_;
    }

    //! Extract interpolation of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const typename Traits::LocalInterpolationType& localInterpolation() const
    {
      return li_;
    }

    //! Return the number of basis functions
    unsigned int size () const
    {
      return lfe1_->size() * lfe2_->size();
    }

    //! Return the prismatic-product of the lfe's geometry types
    const GeometryType type () const
    {
      return GeometryTypes::prismaticProduct(lfe1_->type(), lfe2_->type());
    }
  };

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_PRODUCT_HH
