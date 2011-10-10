// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P0LOCALFINITEELEMENT_HH
#define DUNE_P0LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "p0/p0localbasis.hh"
#include "p0/p0localcoefficients.hh"
#include "p0/p0localinterpolation.hh"

namespace Dune
{

  /** \todo The local p0 finite element on all types of reference elements
      \tparam D Domain data type
      \tparam R Range data type
      \tparam d Dimension of the reference element
   */
  template<class D, class R, int d>
  class P0LocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<P0LocalBasis<D,R,d>, P0LocalCoefficients,
        P0LocalInterpolation<P0LocalBasis<D,R,d> > > Traits;

    /** \todo Please doc me !
     */
    P0LocalFiniteElement (const GeometryType& type)
      : interpolation(type), gt(type)
    {}

    /** \todo Please doc me !
     */
    P0LocalFiniteElement (GeometryType::BasicType basicType) DUNE_DEPRECATED
      : interpolation(GeometryType(basicType,d)), gt(basicType,d)
    {}

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

    P0LocalFiniteElement* clone () const
    {
      return new P0LocalFiniteElement(*this);
    }

  private:
    P0LocalBasis<D,R,d> basis;
    P0LocalCoefficients coefficients;
    P0LocalInterpolation<P0LocalBasis<D,R,d> > interpolation;
    GeometryType gt;
  };

}

#endif
