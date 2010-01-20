// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFINED_P0_LOCALFINITEELEMENT_HH
#define DUNE_REFINED_P0_LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/p0.hh>

#include "refinedp0/refinedp0localbasis.hh"
#include "refinedp0/refinedp0localcoefficients.hh"
#include "refinedp0/refinedp0localinterpolation.hh"

/** \file
    \brief Piecewise P0 finite element
 */
namespace Dune
{

  /** \brief Local finite element that is piecewise P0 on a the once unifromly refined reference geometry
   */
  template<class D, class R, int dim>
  class RefinedP0LocalFiniteElement
  {
    RefinedP0LocalFiniteElement() {}

  public:
    // We steal the p0 traits since they exist for all dim.
    // This allows to instanciate the type and access the Traits.
    typedef typename P0LocalFiniteElement<D,R,dim>::Traits Traits;
  };

  /** \brief Local finite element that is piecewise P0 on a the once unifromly refined reference geometry
   */
  template<class D, class R>
  class RefinedP0LocalFiniteElement<D,R,2>
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<
        RefinedP0LocalBasis<D,R,2>,
        RefinedP0LocalCoefficients<2>,
        RefinedP0LocalInterpolation<RefinedP0LocalBasis<D,R,2> > > Traits;

    /** \todo Please doc me !
     */
    RefinedP0LocalFiniteElement ()
    {
      gt.makeTriangle();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

    RefinedP0LocalFiniteElement * clone () const
    {
      return new RefinedP0LocalFiniteElement(*this);
    }

  private:
    RefinedP0LocalBasis<D,R,2> basis_;
    RefinedP0LocalCoefficients<2> coefficients_;
    RefinedP0LocalInterpolation<RefinedP0LocalBasis<D,R,2> > interpolation_;
    GeometryType gt;
  };


}

#endif