// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_HERMITE_HH
#define DUNE_LOCALFUNCTIONS_HERMITE_HH
#include <array>
#include <numeric>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/localfunctions/common/derivative.hh>
#include <dune/localfunctions/utility/polynomialbasiscoefficients.hh>
#include <dune/localfunctions/utility/polynomialbasis.hh>

namespace Dune
{
  namespace Impl
  {

    /**
     * \brief Implementation of hermite Polynomials using PolynomialBasisWithMatrix
     * \tparam D Type to represent the field in the domain
       \tparam R Type to represent the field in the range
       \tparam dim Dimension of the domain simplex
    */
    template <class D, class R, unsigned int dim>
    class HermiteLocalBasis:
        private PolynomialBasisWithMatrix<StandardEvaluator <VirtualMonomialBasis<dim, D>>, SparseCoeffMatrix<double, 1>, D, R>
    {
    public:
      using Eval = StandardEvaluator<VirtualMonomialBasis<dim, D>>;
      using Base = PolynomialBasisWithMatrix<Eval, SparseCoeffMatrix<double, 1>, D, R>;
      using Traits = LocalBasisTraits<D, dim, FieldVector<D, dim>, R, 1, FieldVector<R, 1>, FieldMatrix<R, 1, dim>>;

    private:
      using MBasisFactory = MonomialBasisProvider<dim, D>;

    public:
      static constexpr unsigned int coeffSize = (dim == 1) ? 4 : (dim == 2) ? 10
                                                                            : 20;
      HermiteLocalBasis()
          : Base(*MBasisFactory::template create<GeometryTypes::simplex(dim)>(3))
      {
        assert(coeffSize == this->basis().size());

        if (dim <= 3)
        {
          this->fill(PolynomialBasisCoefficients::getHermiteCoefficients<double, dim>());
        }
        else
          DUNE_THROW(Dune::NotImplemented, "only implemented for dim <= 3");
      }

      HermiteLocalBasis(HermiteLocalBasis const &) : HermiteLocalBasis() {}

      static unsigned int size() { return coeffSize; }

      unsigned int order() const { return Base::order(); }

      void evaluateFunction(const typename Traits::DomainType &in,
                            std::vector<typename Traits::RangeType> &out) const
      {
        out.resize(size());
        Base::evaluateFunction(in, out);
      }

      /** \brief Evaluate Jacobians of all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Jacobians of all shape functions at that point
       */
      void evaluateJacobian(const typename Traits::DomainType &in,
                            std::vector<typename Traits::JacobianType> &out) const
      {
        out.resize(size());
        Base::evaluateJacobian(in, out);
      }

      /** \brief Evaluate partial derivatives of all shape functions at a given point
       *
       * \param[in] order The partial derivative to be computed, as a multi-index
       * \param[in] in  The evaluation point
       * \param[out] out Jacobians of all shape functions at that point
       */
      void partial(const std::array<unsigned int, dim> &order,
                   const typename Traits::DomainType &in,
                   std::vector<typename Traits::RangeType> &out) const
      {
        out.resize(size());
        Base::partial(order, in, out);
      }
    };

    /** \brief Associations of the Hermite degrees of freedom to subentities of the
     * reference simplex
     *
     * \tparam dim Dimension of the reference simplex
     */
    template <unsigned int dim>
    class HermiteLocalCoefficients
    {
      public:
      using size_type = std::size_t;
      private:
      static constexpr size_type innerDofCodim = (dim == 2) ? 0 : 1;

    public:
      HermiteLocalCoefficients() : localKeys_(size())
      {
        if constexpr (dim > 0 && dim <= 3)
        {
          size_type numberOfVertices = dim + 1;
          size_type numberOfInnerDofs = (dim - 1) * (dim - 1); // probably incorrect for dim > 3
          for (size_type i = 0; i < numberOfVertices; ++i) // subentities: vertices
          {
            for (size_type k = 0; k < numberOfVertices; ++k) // dim + 1 dofs per subentity
              localKeys_[numberOfVertices * i + k] = LocalKey(i, dim, k);
          }
          for (size_type i = 0; i < numberOfInnerDofs; ++i)                          // subentities: element
            localKeys_[numberOfVertices * numberOfVertices + i]
                = LocalKey(i, innerDofCodim, 0); // inner dofs
        }
        else
          DUNE_THROW(NotImplemented,
                     "HermiteLocalCoefficients only implemented for dim <= 3!");
      }

      //! number of coefficients
      static constexpr size_type size()
      {
        if (dim > 3)
          DUNE_THROW(NotImplemented,
                     "HermiteLocalCoefficients only implemented for dim==1!");
        return dim == 1 ? 4 : dim == 2 ? 10
                                       : 20;
      }

      //! get i'th index
      const LocalKey &localKey(std::size_t i) const { return localKeys_[i]; }

    private:
      std::vector<LocalKey> localKeys_;
    };

    /** \brief
     * Note: The wrapper classes in c1elements/functions/functionspacebases do not use this class.
     *
     * Evaluate the degrees of freedom of a Hermite basis. This class provides a template
     * hook pattern, which allows switching between the local(reference) Interpolation and the
     * global(physical) Interpolation if the interpolation function f provides a free derivative()
     * function. If so, we have a local interpolation if derivative(f) returns the local derivative
     * and global interpolation if derivative(f) returns the global derivative.
     * If f has no free derivative() function, a finite difference approach is used to approximate
     * the local (!) derivative.
     *
     *  \tparam LocalBasis The corresponding set of shape functions
     */
    template <class LocalBasis>
    class HermiteLocalInterpolation
    {
      using size_type = std::size_t;
      using LocalCoordinate = typename LocalBasis::Traits::DomainType;
      static constexpr size_type dim = LocalBasis::Traits::dimDomain;

      static constexpr size_type numberOfVertices = dim + 1;
      static constexpr size_type innerDofCodim = (dim == 2) ? 0 : 1; // probably wrong for dim > 3
      static constexpr size_type numberOfInnerDofs
          = (dim - 1) * (dim - 1); // probably wrong for dim > 3

    public:
      /** \brief Evaluate a given function at the Lagrange nodes
       *
       * \tparam F Type of function to evaluate
       * \tparam C Type used for the values of the function
       * \param[in] ff Function to evaluate
       * \param[out] out Array of function values
       */
      template <typename F, typename C>
      void interpolate(const F &ff, std::vector<C> &out) const
      {
        auto &&f = Impl::makeFunctionWithCallOperator<LocalCoordinate>(ff);
        auto &&df = makeDerivative<LocalCoordinate>(
            ff); // makeDerivative at this position might return
                 // reference gradient for localfunctions!
                 // This is necessary for the test routine to compile, TODO change that

        out.resize(LocalBasis::size());
        auto refElement = Dune::ReferenceElements<double, dim>::simplex();
        if (dim > 0 && dim <= 3)
        {
          for (size_type i = 0; i < numberOfVertices; ++i)
          {
            LocalCoordinate x = refElement.position(i, dim);
            // std::cout << "Vertex " << i << ": " << x << std::endl;
            out[numberOfVertices * i] = f(x);

            auto const &grad = df(x);
            for (size_type d = 0; d < dim; ++d)
              out[numberOfVertices * i + d + 1] = getPartialDerivative(grad, d);
          }
          for (size_type i = 0; i < numberOfInnerDofs; ++i)
          {
            out[numberOfVertices * numberOfVertices + i] = f(refElement.position(i, innerDofCodim));
          }
        }
        else
          DUNE_THROW(Dune::NotImplemented,
                     "HermiteInterpolation only implemented for dim <= 3!");
      }

    protected:
      template <class DerivativeType, class FieldType>
      FieldType getPartialDerivative(DerivativeType const &df,
                                     std::size_t i) const
      {
        DUNE_THROW(Dune::NotImplemented,
                   "Derivative Type is neither FieldMatrix<double,1,d> nor "
                   "FieldVector<double,d>");
      }

      template <class FieldType, int d>
      FieldType getPartialDerivative(Dune::FieldVector<FieldType, d> const &df,
                                     std::size_t i) const
      {
        return df[i];
      }

      template <class FieldType, int d>
      FieldType getPartialDerivative(Dune::FieldMatrix<FieldType, 1, d> const &df,
                                     std::size_t i) const
      {
        if (df.N() == 1)
          return df[0][i];
        else if (df.M() == 1)
          return df[i][0];
        else
          DUNE_THROW(Dune::NotImplemented, "derivative of scalar function is a matrix!");
      }
    };

  } // namespace Impl

  /** \brief Hermite finite element for simplices
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   * \tparam dim dimension of the reference element
   */
  template <class D, class R, unsigned int dim>
  class HermiteLocalFiniteElement
  {

  public:
    /** \brief Export number types, dimensions, etc.
   */

    using Traits =
        LocalFiniteElementTraits<Impl::HermiteLocalBasis<D, R, dim>,
                                 Impl::HermiteLocalCoefficients<dim>,
                                 Impl::HermiteLocalInterpolation<Impl::HermiteLocalBasis<D, R, dim>>>;

    /** \brief Default constructor
     * TODO remove?
     * \deprecated This explicit implementation only exists to work around a bug
     * in clang 3.8 which disappeared in clang 6
     */
    HermiteLocalFiniteElement() {}

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType &localBasis() const { return basis_; }

    /** \brief Returns the assignment of the degrees of freedom to the element
     * subentities
     */
    const typename Traits::LocalCoefficientsType &localCoefficients() const
    {
      return coefficients_;
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType &localInterpolation() const
    {
      return interpolation_;
    }

    /** \brief The number of shape functions */
    static constexpr std::size_t size()
    {
      return dim == 1 ? 4 : dim == 2 ? 10
                                     : 20;
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type() { return GeometryTypes::simplex(dim); }

  private:
    Impl::HermiteLocalBasis<D, R, dim> basis_;
    Impl::HermiteLocalCoefficients<dim> coefficients_;
    Impl::HermiteLocalInterpolation<Impl::HermiteLocalBasis<D, R, dim>> interpolation_;
  };

} // namespace Dune

#endif //DUNE_C1ELEMENTS_HERMITE_HH
