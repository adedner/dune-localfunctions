// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_HH

/** \file
 * \brief Convenience header that includes all implementations of Lagrange finite elements
 */

#include <cassert>
#include <map>

#include <dune/geometry/type.hh>

// Headers for Lagrange elements with run-time order
#include <dune/localfunctions/utility/localfiniteelement.hh>
#include <dune/localfunctions/utility/dglocalcoefficients.hh>

#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/lagrangecoefficients.hh>
#include <dune/localfunctions/lagrange/interpolation.hh>
#include <dune/localfunctions/lagrange/lagrangebasis.hh>

// Headers for Lagrange elements with compile-time order
#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/p1.hh>
#include <dune/localfunctions/lagrange/p2.hh>
#include <dune/localfunctions/lagrange/p23d.hh>
#include <dune/localfunctions/lagrange/pk.hh>
#include <dune/localfunctions/lagrange/pk1d.hh>
#include <dune/localfunctions/lagrange/pk2d.hh>
#include <dune/localfunctions/lagrange/pk3d.hh>
#include <dune/localfunctions/lagrange/pq22d.hh>

#include <dune/localfunctions/lagrange/q1.hh>
#include <dune/localfunctions/lagrange/q2.hh>
#include <dune/localfunctions/lagrange/qk.hh>

#include <dune/localfunctions/lagrange/prismp1.hh>
#include <dune/localfunctions/lagrange/prismp2.hh>
#include <dune/localfunctions/lagrange/pyramidp1.hh>
#include <dune/localfunctions/lagrange/pyramidp2.hh>


namespace Dune
{
  /**
   * @brief Lagrange local finite elements for a given set of interpolation
   *        points.
   *
   * The class LP provides the points for the interpolation.
   * It has two template arguments, the first is the Field type to
   * use for evaluating the points the second the dimension
   * of the reference elements on which to construct the points.
   * It is instantiated with the desired order and has a template
   * method build taking a Topology to construct the points
   * (a std::vector of FieldVectors).
   * It also provides a static template method supports to indicate
   * if the point set can be build for a specified Topology.
   *
   * Examples include:
   * - EquidistantPointSet:  standard point set for Lagrange points
   * - LobattoPointSet:      an approximate Freget type point set
   *                         (provided for simplex and generalized prism
   *                         topologies (i.e. not for a 3d pyramid)
   *
   * \ingroup Lagrange
   *
   * \tparam LP a template class defining the points for the lagrange interpolation
   * \tparam dimDomain dimension of reference elements
   * \tparam D domain for basis functions
   * \tparam R range for basis functions
   * \tparam SF storage field for basis matrix
   * \tparam CF compute field for basis matrix
   **/
  template< template <class,unsigned int> class LP,
      unsigned int dimDomain, class D, class R,
      class SF=R, class CF=SF >
  class LagrangeLocalFiniteElement
    : public GenericLocalFiniteElement< LagrangeBasisFactory< LP, dimDomain, SF, CF >,
          LagrangeCoefficientsFactory<LP, dimDomain, SF >,
          LagrangeInterpolationFactory< LP, dimDomain, SF > >
  {
    typedef GenericLocalFiniteElement< LagrangeBasisFactory< LP, dimDomain, SF, CF >,
        LagrangeCoefficientsFactory<LP, dimDomain, SF >,
        LagrangeInterpolationFactory< LP, dimDomain, SF > > Base;
  public:
    typedef typename Base::Traits Traits;

    /** \todo Please doc me !
     */
    LagrangeLocalFiniteElement ( const GeometryType &gt, unsigned int order )
      : Base( gt, order )
    {}
  };


  /// \brief Cache for the generic \ref LagrangeLocalFiniteElement
  /**
   * Caches a LagrangeLocalFiniteElement instance for each GeometryType. This cache is
   * initialized with the polynomial order of the Lagrange basis functions. The concrete
   * local finite-element can be obtained by calling `cache.get(GeometryType)`.
   *
   * \tparam D          Field-type of the domain of the local basis functions
   * \tparam R          Field-type of the range of the local basis functions
   * \tparam dimDomain  Dimension of reference elements
   * \tparam LP         Lagrange point-set template <class Field, unsigned int dim>,
   *                    default: \ref EquidistantPointSet
   **/
  template <class D, class R, int dimDomain,
            template <class,unsigned int> class LP = EquidistantPointSet>
  class LagrangeLFECache
  {
  public:
    using FiniteElementType = LagrangeLocalFiniteElement<LP, dimDomain, D, R>;

    /// \brief store the polynomial order of the Lagrange basis functions. The order is
    /// required to be greater than zero.
    LagrangeLFECache (int order)
      : order_(order)
    {
      assert(order > 0);
    }

    /// \brief obtain a local finite-element for the given GeometryType and stored order.
    const FiniteElementType& get (GeometryType type)
    {
      auto it = data_.find(type);
      if (it == data_.end())
        it = data_.emplace(type,FiniteElementType(type,order_)).first;
      return it->second;
    }

  private:
    unsigned int order_;
    std::map<GeometryType, FiniteElementType> data_;
  };
}

#endif // #ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_HH
