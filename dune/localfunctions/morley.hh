// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_MORLEY_HH
#define DUNE_LOCALFUNCTIONS_MORLEY_HH
#include <array>
#include <numeric>

#include <dune/common/classname.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/common/transpose.hh>
#include <dune/localfunctions/common/derivative.hh>
#include <dune/localfunctions/polynomialbasiscoefficients.hh>
#include <dune/localfunctions/utility/polynomialbasis.hh>

#include <dune/c1elements/tensormatvec.hpp>

namespace Dune
{
  namespace Impl
  {

    /**
     * \brief Implementation of Morley Polynomials using PolynomialBasisWithMatrix
     * \tparam D Type to represent the field in the domain
       \tparam R Type to represent the field in the range
    */

    template <class D, class R>
    class MorleyLocalBasis:
        private PolynomialBasisWithMatrix<StandardEvaluator<VirtualMonomialBasis<2, D>>,
                                         SparseCoeffMatrix<double, 1>, D, R>
    {
    public:
      static constexpr unsigned int dim = 2;
      using Eval = StandardEvaluator<VirtualMonomialBasis<2, D>>;
      using Base = PolynomialBasisWithMatrix<Eval, SparseCoeffMatrix<double, 1>, D, R>;
      using Traits = LocalBasisTraits<D, 2, FieldVector<D, 2>, R, 1, FieldVector<R, 1>, FieldMatrix<R, 1, 2>>;
      using HessianType = FieldVector<FieldMatrix<R,dim,dim>,1>;

    private:
      using MBasisFactory = MonomialBasisProvider<2, D>;
      // Number of edges of the reference simplex(triangle)
      constexpr static std::size_t numberOfEdges = 3;
      // TODO test with int or short
      using OrientationType = double;

    public:
      static constexpr unsigned int coeffSize = 6;

      /** \brief Constructor with default edge orientation
       *
       * Default orientation means that all edges point from the vertex with lower index
       * to the vertex with higher index. The tangentials point in the direction of the edge.
       */
      MorleyLocalBasis()
          : Base(*MBasisFactory::template create<GeometryTypes::simplex(dim)>(2)),
            edgeOrientation_({1., 1., 1.}) // order = 2
      {
        assert(coeffSize == this->basis().size());

        Base::fill(PolynomialBasisCoefficients::getMorleyCoefficients<double>());
      }

      /**
       * \brief Constructor with given edge orientation
       */
      MorleyLocalBasis(std::bitset<numberOfEdges> edgeOrientation) : MorleyLocalBasis()
      {
        for (std::size_t i = 0; i < edgeOrientation_.size(); i++)
          edgeOrientation_[i] = edgeOrientation[i] ? -1.0 : 1.0;
      }

      // Necessary in some cases, like usage via Virtual Finite Element ... creates a new
      // coefficient matrix ...
      MorleyLocalBasis(MorleyLocalBasis const &other) : MorleyLocalBasis()
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
        // First we evaluate the basis with default edge orientation, then we adapt the result
        // accordingly This works because the edge orientation only changes the "wrong" derivative
        // nodes which are multiplied by -1, and the very same happens to the nodal basis
        out.resize(size());
        Base::evaluateFunction(in, out);
        for (std::size_t i = 0; i < numberOfEdges; i++)
          out[3 + i] *= edgeOrientation_[i];
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
          out[3 + i] *= edgeOrientation_[i];
      }


      /** \brief Evaluate Jacobians of all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Jacobians of all shape functions at that point
       */
      void evaluateHessian(const typename Traits::DomainType &in,
                            std::vector<HessianType> &out) const
      {
        out.resize(size());
        Base::evaluateHessian(in, out);
        for (std::size_t i = 0; i < numberOfEdges; i++)
          out[3 + i] *= edgeOrientation_[i];
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
          out[3 + i] *= edgeOrientation_[i];
      }
      std::array<OrientationType, numberOfEdges> const &edgeOrientation() const
      {
        return edgeOrientation_;
      }

    private:
      // Orientations of the simplex edges
      std::array<OrientationType, numberOfEdges> edgeOrientation_;
    };

    /** \brief Associations of the Morley degrees of freedom to subentities of the
     * reference simplex
     *
     * \tparam dim Dimension of the reference simplex
     */
    class MorleyLocalCoefficients
    {
      using size_type = std::size_t;

    public:
      MorleyLocalCoefficients() : localKeys_(size())
      {
        for (size_type i = 0; i < 3; ++i)
        {
          localKeys_[i] = LocalKey(i, 2, 0);     // vertex dofs
          localKeys_[3 + i] = LocalKey(i, 1, 0); // edge dofs
        }
      }

      //! number of coefficients
      static constexpr size_type size()
      {
        return 6;
      }

      //! get i'th index
      const LocalKey &localKey(std::size_t i) const { return localKeys_[i]; }

    private:
      std::vector<LocalKey> localKeys_;
    };

    /** \brief Note: This class is not used by dune-functions bases!
     * Evaluates the degrees of freedom of a Morley basis
     * Note: The wrapper classes in c1elements/functions/functionspacebases do not use this class.
     *
     * This class provides a template hook pattern, which allows switching between the
     * local(reference) Interpolation and the global(physical) Interpolation.
     * If the interpolation function f provides a method normalDerivative(size_t i), this method
     * is called to evaluate the ith Dof.
     * If not we have a local interpolation if derivative(f) returns the local derivative.
     * If there is no f.normalDerivative(size_t) and derivative(f) returns a global Derivative the
     * result is invalid.
     *
     *  \tparam LocalBasis The corresponding set of shape functions
     */
    template <class LocalBasis>
    class MorleyLocalInterpolation
    {
      using size_type = std::size_t;
      using LocalCoordinate = typename LocalBasis::Traits::DomainType;
      static constexpr unsigned int numberOfEdges = 3;

    public:
      MorleyLocalInterpolation(std::bitset<numberOfEdges> s = 0) : s_(s)
      {
        auto refElement = Dune::referenceElement<double, 2>(GeometryTypes::simplex(2));

        for (std::size_t i = 0; i < numberOfEdges; i++)
        {
          v_[i] = refElement.position(i, 2); // vertices
          m_[i] = refElement.position(i, 1); // edge midpoints

          std::size_t lower = (i == 2) ? 1 : 0;
          std::size_t upper = (i == 0) ? 1 : 2;
          auto edge = refElement.position(upper, 2) - refElement.position(lower, 2);

          normal_[i] = {-edge[1], edge[0]}; // rotation by 90 degree
          // normalize normals
          if (std::abs(normal_[i].two_norm2() - 1) > 1e-12)
            normal_[i] /= normal_[i].two_norm();
        }
      }

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
        // using D = typename LocalBasis::Traits::DomainFieldType;

        auto &&f = Impl::makeFunctionWithCallOperator<LocalCoordinate>(ff);

        out.resize(LocalBasis::size());

        for (size_type i = 0; i < 3; ++i)
        {
          out[i] = f(v_[i]);
          out[3 + i] = normalDerivative(ff, i);
          if (s_[i])
            out[3 + i] *= -1.;
        }
      }

    protected:
      // Vertices of the reference simplex
      std::array<LocalCoordinate, numberOfEdges> v_;
      // Edge midpoints of the reference simplex
      std::array<LocalCoordinate, numberOfEdges> m_;

      // normals of the reference simplex
      std::array<LocalCoordinate, numberOfEdges> normal_;
      // orientation bitset
      std::bitset<numberOfEdges> s_;

      // Infrastructure for normal Derivative that allows evaluation of default oriented global
      // normal derivative, if f has this method
      template <class F>
      auto normalDerivative(F const &f, size_type i) const
      {
        return normalDerivativeImpl(f, i, PriorityTag<42>{});
      }

      template <class F,
                decltype((std::declval<F>().normalDerivative(std::declval<size_type>()), true))
                = true>
      auto normalDerivativeImpl(F const &f, size_type i, PriorityTag<4>) const
      {
        return f.normalDerivative(i);
      }

      template <class F, decltype((derivative(std::declval<F>()), true)) = true>
      auto normalDerivativeImpl(F const &f, size_type i, PriorityTag<3>) const
      {
        return matrixToVector(makeDerivative<LocalCoordinate>(f)(m_[i])).dot(normal_[i]);
      }

      template <class F>
      auto normalDerivativeImpl(F const &f, size_type i, PriorityTag<1>) const
      {
        DUNE_THROW(Dune::NotImplemented,
                   Dune::className(f)
                       + " supports neither derivative(f) nor f.normalDerivative(i)!");
        return 0;
      }
    };

  } // namespace Impl

  /** \brief Morley finite element for simplices
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   */
  template <class D, class R>
  class MorleyLocalFiniteElement
  {

  public:
    /** \brief Export number types, dimensions, etc.
     */

    using Traits
        = LocalFiniteElementTraits<Impl::MorleyLocalBasis<D, R>, Impl::MorleyLocalCoefficients,
                                   Impl::MorleyLocalInterpolation<Impl::MorleyLocalBasis<D, R>>>;

    MorleyLocalFiniteElement(std::bitset<3> s = 0)
        : basis_(s), interpolation_(s) {}

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
    static constexpr std::size_t size() { return 6; }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type() { return GeometryTypes::simplex(2); }

  private:
    typename Traits::LocalBasisType basis_;
    typename Traits::LocalCoefficientsType coefficients_;
    typename Traits::LocalInterpolationType interpolation_;
  };

} // namespace Dune

#endif // DUNE_C1ELEMENTS_MORLEY_HH
