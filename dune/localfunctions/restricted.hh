// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNTIONS_RESTRICTED_HH
#define DUNE_LOCALFUNTIONS_RESTRICTED_HH

namespace Dune
{
  template<class LFE>
  class RestrictedLocalFiniteElementView;

  namespace Impl
  {

    template<class LFE>
    class RestrictedLocalBasis
    {
      const RestrictedLocalFiniteElementView<LFE> & fem_;

    public:
      using Traits = typename LFE::Traits::LocalBasisType::Traits;

      RestrictedLocalBasis(const RestrictedLocalFiniteElementView<LFE> & fem) : fem_(fem) {}

      //! \brief number of shape functions
      unsigned int size () const
      {
        return fem_.size();
      }

      //! \brief Evaluate all shape functions
      inline void evaluateFunction (const typename Traits::DomainType& in,
        std::vector<typename Traits::RangeType>& out) const
      {
        if (fem_.active())
          fem_.impl().localBasis().evaluateFunction(in,out);
        else
          out.resize(0);
      }

      //! \brief Evaluate Jacobian of all shape functions
      inline void
      evaluateJacobian (const typename Traits::DomainType& in,
        std::vector<typename Traits::JacobianType>& out) const
      {
        if (fem_.active())
          fem_.impl().localBasis().evaluateJacobian(in,out);
        else
          out.resize(0);
      }

      //! \brief Evaluate partial derivatives of all shape functions
      void partial (const std::array<unsigned int, 2>& order,
        const typename Traits::DomainType& in,
        std::vector<typename Traits::RangeType>& out) const
      {
        if (fem_.active())
          fem_.impl().localBasis().partial(order,in,out);
        else
          out.resize(0);
      }

      //! \brief Polynomial order of the shape functions
      unsigned int order () const
      {
        return fem_.impl().localBasis().order();
      }

    };

    template<class LFE>
    class RestrictedLocalCoefficients
    {
      const RestrictedLocalFiniteElementView<LFE> & fem_;
    public:
      RestrictedLocalCoefficients(const RestrictedLocalFiniteElementView<LFE> & fem) : fem_(fem) {}

      //! number of coefficients
      std::size_t size () const
      {
        return fem_.size();
      }

      //! get i'th index
      const LocalKey& localKey (std::size_t i) const
      {
        return fem_.impl().localCoefficients().localKey(i);
      }
    };

    template<class LFE>
    class RestrictedLocalInterpolation
    {
      const RestrictedLocalFiniteElementView<LFE> & fem_;
    public:
      RestrictedLocalInterpolation(const RestrictedLocalFiniteElementView<LFE> & fem) : fem_(fem) {}

      template<typename F, typename C>
      void interpolate (const F& f, std::vector<C>& out) const
      {
        if (fem_.active())
          fem_.impl().localInterpolation().interpolate(f,out);
        else
          out.resize(0);
      }
    };

  }

  /**
   * \brief Restriction wrapping LocalFiniteElement classes
   *
   * This is a wrapper class to restrict the support of LocalFiniteElement
   * implementations. This requires the ability to turn-on / tunr-off
   * this implementation.
   *
   * An example use case are dsicretizations on a subdomain of the gridView.
   *
   * The implementation wraps a given LocalFiniteElement and adds
   * new methods `active()`, `activate()`, and `deactivate()`.
   *
   * If the RestrictedLocalFiniteElement is inactive, the actual
   * implementation will not be changed, but all sizes are set to 0.
   *
   * \tparam LocalFiniteElement implementation to restrict
   */
  template<class LocalFiniteElement>
  class RestrictedLocalFiniteElementView
  {

    // In each LocalFooVariant we store a std::variant<std::monostate, const FooImpl*...>, i.e. a std::variant
    // with the pointer to the Foo implementation unless LocalFiniteElementVariant stores a monostate. In this
    // case each LocalFooVariant also stores a monostate (and not a monostate*).
    using LocalBasis = Impl::RestrictedLocalBasis<LocalFiniteElement>;
    using LocalCoefficients = Impl::RestrictedLocalCoefficients<LocalFiniteElement>;
    using LocalInterpolation = Impl::RestrictedLocalInterpolation<LocalFiniteElement>;

    friend LocalBasis;
    friend LocalCoefficients;
    friend LocalInterpolation;

    RestrictedLocalFiniteElementView() = delete;

  public:

    /**
     * \brief Export LocalFiniteElementTraits
     */
    using Traits = typename Dune::LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation>;

    /**
     * \brief Construct RestrictedLocalFiniteElement
     */
    RestrictedLocalFiniteElementView(const LocalFiniteElement& impl) :
      impl_(&impl),
      localBasis_(*this),
      localCoefficients_(*this),
      localInterpolation_(*this)
    {
      std::cout << "RestrictedLocalFiniteElementView ... " << impl_ << std::endl;
    }

    /**
     * \brief Copy constructor
     */
    RestrictedLocalFiniteElementView(const RestrictedLocalFiniteElementView& other) :
      impl_(other.impl_),
      localBasis_(*this),
      localCoefficients_(*this),
      localInterpolation_(*this)
    {
      std::cout << "RestrictedLocalFiniteElementView (copy) ... " << impl_ << std::endl;
    }

    /**
     * \brief Copy assignment
     */
    RestrictedLocalFiniteElementView& operator=(const RestrictedLocalFiniteElementView& other)
    {
      impl_ = other.impl_;
      std::cout << "RestrictedLocalFiniteElementView (copy assign) ... " << impl_ << std::endl;
      return *this;
    }

    /**
     * \brief Assignment from implementation
     */
    RestrictedLocalFiniteElementView& operator=(const LocalFiniteElement& impl)
    {
      impl_ = &impl;
      std::cout << "RestrictedLocalFiniteElementView (assign) ... " << impl_ << std::endl;
      return *this;
    }


    /**
     * \brief Provide access to LocalBasis implementation of this LocalFiniteElement
     */
    const typename Traits::LocalBasisType& localBasis() const
    {
      return localBasis_;
    }

    /**
     * \brief Provide access to LocalCoefficients implementation of this LocalFiniteElement
     */
    const typename Traits::LocalCoefficientsType& localCoefficients() const
    {
      return localCoefficients_;
    }

    /**
     * \brief Provide access to LocalInterpolation implementation of this LocalFiniteElement
     */
    const typename Traits::LocalInterpolationType& localInterpolation() const
    {
      return localInterpolation_;
    }

    /**
     * \brief Number of shape functions
     */
    unsigned int size() const
    {
      return active() ? impl().size() : 0;
    }

    /**
     * \brief GeometryType for this LocalFiniteElement
     */
    constexpr GeometryType type() const
    {
      return active() ? impl().type() : type_;
    }

    /**
     * \brief get current active/inactive state
     */
    bool active() const
    {
      return active_;
    }

    /**
     * \brief set current state to active
     */
    void activate()
    {
      active_ = true;
    }

    /**
     * \brief set current state to inactive and provide alternative geometry type
     */
    void deactivate()
    {
      active_ = false;
    }

    void setDefaultGeometry(GeometryType type)
    {
      type_ = type;
    }

  private:
    const LocalFiniteElement & impl() const { return *impl_; }

    // std::reference_wrapper<const LocalFiniteElement> impl_;
    const LocalFiniteElement* impl_;
    bool active_;
    LocalBasis localBasis_;
    LocalCoefficients localCoefficients_;
    LocalInterpolation localInterpolation_;
    GeometryType type_;
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNTIONS_RESTRICTED_HH
