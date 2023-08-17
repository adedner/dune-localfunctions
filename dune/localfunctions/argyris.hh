// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_ARGYRIS_HH
#define DUNE_LOCALFUNCTIONS_ARGYRIS_HH
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

#include <dune/localfunctions/tensormatvec.hh>
#include <dune/localfunctions/morley.hh>
#include <dune/localfunctions/polynomialbasiscoefficients.hh>
#include <dune/localfunctions/utility/polynomialbasis.hh>
namespace Dune
{
  namespace Impl
  {

    /**
       * \brief Implementation of Argyris shape Function using
       PolynomialBasisWithMatrix
       * \tparam D Type to represent the field in the domain
          \tparam R Type to represent the field in the range
          \tparam dim Dimension of the domain simplex
       */

    template <class D, class R>
    class ArgyrisLocalBasis:
        private PolynomialBasisWithMatrix<StandardEvaluator <VirtualMonomialBasis<2, D>>, SparseCoeffMatrix<double, 1>, D, R>
    {
    public:
      static constexpr unsigned int dim = 2;
      using Eval = StandardEvaluator<VirtualMonomialBasis<dim, D>>;
      using Base = PolynomialBasisWithMatrix<StandardEvaluator<VirtualMonomialBasis<dim, D>>,
                                             SparseCoeffMatrix<double, 1>, D, R>;
      using Traits =
          LocalBasisTraits<D, dim, FieldVector<D, dim>, R, 1, FieldVector<R, 1>, FieldMatrix<R, 1, dim>>;
      using HessianType = FieldVector<FieldMatrix<R,dim,dim>,1>;

    private:
      using MBasisFactory = MonomialBasisProvider<dim, D>;
      // Number of edges of the reference simplex(triangle)
      constexpr static std::size_t numberOfEdges = 3;
      // TODO test with int or short
      using OrientationType = double;

    public:
      static constexpr unsigned int coeffSize = 21;

      /** \brief Constructor with default edge orientation
       *
       * Default orientation means that all edges point from the vertex with lower index
       * to the vertex with higher index. The tangentials point in the direction of the edge.
       */
      ArgyrisLocalBasis()
          : Base(*MBasisFactory::template create<GeometryTypes::simplex(dim)>(
              5)) // Use monomials up to order 5!
      {
        assert(coeffSize == this->basis().size());
        this->fill(PolynomialBasisCoefficients::getArgyrisCoefficients<double>());
      }

      ArgyrisLocalBasis(std::bitset<numberOfEdges> edgeOrientation) : ArgyrisLocalBasis()
      {
        for (std::size_t i = 0; i < edgeOrientation_.size(); ++i)
          edgeOrientation_[i] = edgeOrientation[i] ? -1.0 : 1.0;
      }

      ArgyrisLocalBasis(ArgyrisLocalBasis const &other) : ArgyrisLocalBasis()
      {
        edgeOrientation_ = other.edgeOrientation_;
      }

      static unsigned int size() { return coeffSize; }

      unsigned int order() const { return Base::order(); }

      /** \brief Evaluate values of all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Jacobians of all shape functions at that point
       */
      void evaluateFunction(const typename Traits::DomainType &in,
                            std::vector<typename Traits::RangeType> &out) const
      {
        out.resize(size());
        Base::evaluateFunction(in, out);
        for (std::size_t i = 0; i < numberOfEdges; i++)
          out[18 + i] *= edgeOrientation_[i];
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
        for (std::size_t i = 0; i < numberOfEdges; i++)
          out[18 + i] *= edgeOrientation_[i];
      }

        /** \brief Evaluate Hessians of all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Hessians of all shape functions at that point
       */
      void evaluateHessian(const typename Traits::DomainType &in,
                            std::vector<HessianType> &out) const
      {
        out.resize(size());

        Base::evaluateHessian(in, out);
        for (std::size_t i = 0; i < numberOfEdges; i++)
          out[18 + i] *= edgeOrientation_[i];
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
        for (std::size_t i = 0; i < numberOfEdges; i++)
          out[18 + i] *= edgeOrientation_[i];
      }

      std::array<OrientationType, numberOfEdges> const &edgeOrientation() const
      {
        return edgeOrientation_;
      }

    private:
      // Orientations of the simplex edges
      std::array<OrientationType, numberOfEdges> edgeOrientation_;
    };

    /** \brief Associations of the Argyris degrees of freedom to
     * subentities of the reference triangle
     */
    class ArgyrisLocalCoefficients
    {
    public:
      static std::size_t const dim = 2;

      ArgyrisLocalCoefficients()
          : localKeys_(size())
      {
        for (unsigned int i = 0; i < 3; i++)             // subentities: three vertices
          for (unsigned int k = 0; k < 6; k++)           // 6 basis functions per vertex
            localKeys_[6 * i + k]
                = LocalKey(i, dim, k); //(subentity, codim, number of dof for this subentity)
        for (unsigned int i = 0; i < 3; ++i)             // the edges
          localKeys_[18 + i] = LocalKey(i, dim - 1, 0);  // one node per edge
        return;
      }

      //! number of coefficients
      static constexpr std::size_t size()
      {
        return 21;
      }

      //! get i'th index
      const LocalKey &localKey(std::size_t i) const
      {
        return localKeys_[i];
      }

    private:
      std::vector<LocalKey> localKeys_;
    };

    /** \brief Evaluate the degrees of freedom of a Argyris basis
     *
     * Note: The wrapper classes in c1elements/functions/functionspacebases do not use this class.
     *
     * This class provides a template hook pattern, which in principle allows it to be used as both
     * local(reference) Interpolation and global(physical) Interpolation, depending on
     * implemented derivatives of the interpolated function.
     *
     * The Requirements at f are a combination of the ones documented in
     * MorleyLocalInterpolation and HermiteLocalInterpolation
     * Two combinations of "legal" derivative interfaces arise:
     * 1. Local interpolation:
     *    This is given if
     *    (i) derivative(f) returns the local gradient and f has no normalDeritvative
     *      method or f.normalderivative(i) returns the reference normal derivative or
     *    (ii) neither derivative(f) nor f.normalDerivative(i) is implemented - Finite Difference!
     *      Additionally derivative(derivative(f)) has to return the local hessian, if implemented
     * 2. Global interpolation:
     *    This is given iff
     *    derivative(f) returns the global gradient
     *    and derivative(derivative(f)) is implemented and returns the global hessian.
     *    and f has a normalDerivative(i) method that returns the global normal derivative
     * All other possibilities result in a mixture of local and global dofs, in particular,
     * if derivative(f) returns a global derivative which has no free derivative() function, the
     * result are wrong!!!
     *
     *  \tparam LocalBasis The corresponding set of shape functions
     */
    template <class LocalBasis>
    class ArgyrisLocalInterpolation: MorleyLocalInterpolation<LocalBasis>
    {

      using D = typename LocalBasis::Traits::DomainFieldType;
      using LocalCoordinate = typename LocalBasis::Traits::DomainType;
      static constexpr unsigned int numberOfEdges = 3;
      using Base = MorleyLocalInterpolation<LocalBasis>;
      using Base::normalDerivative;
      using Base::s_;
      using Base::v_;

    public:
      ArgyrisLocalInterpolation(std::bitset<numberOfEdges> s = 0) : Base(s) {}

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
        // constexpr auto dim = LocalBasis::Traits::dimDomain;

        auto &&f = Impl::makeFunctionWithCallOperator<typename LocalBasis::Traits::DomainType>(ff);

        // auto df = derivative(ff);
        // auto hessf = derivative(derivative(ff));

        auto df = makeDerivative<LocalCoordinate>(ff);
        auto hessf = makeDerivative<LocalCoordinate>(df);

        out.resize(LocalBasis::size());

        // Iterate over vertices, 6 dofs per vertex
        for (unsigned int i = 0; i < 3; ++i)
        {
          auto const &derivativeValue = df(v_[i]);
          auto const &hessianValue = hessf(v_[i]);

          out[i * 6 + 0] = f(v_[i]);
          out[i * 6 + 1] = matrixToVector(derivativeValue)[0];
          out[i * 6 + 2] = matrixToVector(derivativeValue)[1];
          out[i * 6 + 3] = matrixToVector(hessianValue[0])[0];
          out[i * 6 + 4] = matrixToVector(hessianValue[0])[1];
          out[i * 6 + 5] = matrixToVector(hessianValue[1])[1];
        }

        // iterate over edges, one dof per edge
        for (unsigned int i = 0; i < numberOfEdges; ++i)
        {
          out[18 + i] = normalDerivative(ff, i) * (s_[i] ? -1 : 1);
        }
        return;
      }
    };
  } // namespace Impl
  /**
   * \brief Argyris finite element for triangles
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   */
  template <class D, class R>
  class ArgyrisLocalFiniteElement
  {
  public:
    using Traits = LocalFiniteElementTraits<Impl::ArgyrisLocalBasis<D, R>,
                                            Impl::ArgyrisLocalCoefficients,
                                            Impl::ArgyrisLocalInterpolation<Impl::ArgyrisLocalBasis<D, R>>>;

    ArgyrisLocalFiniteElement(std::bitset<3> s = 0)
        : basis_(s), interpolation_(s) {}

    const typename Traits::LocalBasisType &localBasis() const { return basis_; }

    /** \brief Returns the assignment of the degrees of freedom to the element subentities
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
      return 21;
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type()
    {
      return GeometryTypes::simplex(2);
    }

  private:
    typename Traits::LocalBasisType basis_;
    typename Traits::LocalCoefficientsType coefficients_;
    typename Traits::LocalInterpolationType interpolation_;
  };
} // Dune

#endif
