// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFINED_P1_LOCALFINITEELEMENT_HH
#define DUNE_REFINED_P1_LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/p0.hh>

#include <dune/localfunctions/refined/refinedp1/refinedp1localbasis.hh>
#include <dune/localfunctions/pk2d/pk2dlocalcoefficients.hh>
#include <dune/localfunctions/pk2d/pk2dlocalinterpolation.hh>
#include <dune/localfunctions/pk2d/pk2dlocalbasis.hh>
#include <dune/localfunctions/pk3d/pk3dlocalcoefficients.hh>
#include <dune/localfunctions/pk3d/pk3dlocalinterpolation.hh>
#include <dune/localfunctions/pk3d/pk3dlocalbasis.hh>

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int dim>
  class RefinedP1LocalFiniteElement
  {
    RefinedP1LocalFiniteElement() {}

  public:
    // We steal the p0 traits since they exist for all dim.
    // This allows to instanciate the type and access the Traits.
    typedef typename P0LocalFiniteElement<D,R,dim>::Traits Traits;
  };

  /** \todo Please doc me !
   */
  template<class D, class R>
  class RefinedP1LocalFiniteElement<D,R,2>
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<RefinedP1LocalBasis<D,R,2>,
        Pk2DLocalCoefficients<2>,
        Pk2DLocalInterpolation<Pk2DLocalBasis<D,R,2> > > Traits;

    /** \todo Please doc me !
     */
    RefinedP1LocalFiniteElement ()
    {
      gt.makeTriangle();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

    RefinedP1LocalFiniteElement * clone () const
    {
      return new RefinedP1LocalFiniteElement(*this);
    }

  private:
    RefinedP1LocalBasis<D,R,2> basis;
    Pk2DLocalCoefficients<2> coefficients;
    Pk2DLocalInterpolation<Pk2DLocalBasis<D,R,2> > interpolation;
    GeometryType gt;
  };

  /** \todo Please doc me !
   */
  template<class D, class R>
  class RefinedP1LocalFiniteElement<D,R,3>
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<RefinedP1LocalBasis<D,R,3>,
        Pk3DLocalCoefficients<2>,
        Pk3DLocalInterpolation<Pk3DLocalBasis<D,R,2> > > Traits;

    /** \todo Please doc me !
     */
    RefinedP1LocalFiniteElement ()
    {
      gt.makeTetrahedron();
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

    RefinedP1LocalFiniteElement* clone () const
    {
      return new RefinedP1LocalFiniteElement(*this);
    }

  private:
    RefinedP1LocalBasis<D,R,3> basis;
    Pk3DLocalCoefficients<2> coefficients;
    Pk3DLocalInterpolation<Pk3DLocalBasis<D,R,2> > interpolation;
    GeometryType gt;
  };

}

#endif