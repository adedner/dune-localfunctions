// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "brezzidouglasmarini1simplex3d/brezzidouglasmarini1simplex3dlocalbasis.hh"
#include "brezzidouglasmarini1simplex3d/brezzidouglasmarini1simplex3dlocalcoefficients.hh"
#include "brezzidouglasmarini1simplex3d/brezzidouglasmarini1simplex3dlocalinterpolation.hh"

namespace Dune
{

  /**
   * \brief First order Brezzi-Douglas-Marini shape functions on tetrahedrons.
   *
   * \ingroup BrezziDouglasMarini
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class BDM1Simplex3DLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<
        BDM1Simplex3DLocalBasis<D,R>,
        BDM1Simplex3DLocalCoefficients,
        BDM1Simplex3DLocalInterpolation<BDM1Simplex3DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    BDM1Simplex3DLocalFiniteElement ()
    {}

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    BDM1Simplex3DLocalFiniteElement (int s) :
      basis(s),
      interpolation(s)
    {}

    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    static constexpr GeometryType type ()
    {
      return GeometryTypes::tetrahedron;
    }

  private:
    BDM1Simplex3DLocalBasis<D,R> basis;
    BDM1Simplex3DLocalCoefficients coefficients;
    BDM1Simplex3DLocalInterpolation<BDM1Simplex3DLocalBasis<D,R> > interpolation;
  };
}
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALFINITEELEMENT_HH
