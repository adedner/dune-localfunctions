// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_BDFMSIMPLEX_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_BDFMSIMPLEX_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

// if enabled use the RT0 version for BDFM0 since these are the same elements
#define BDFM_USE_RT0_BASIS

#include <dune/localfunctions/brezzidouglasfortinmarini/simplex/localbasis.hh>
#include <dune/localfunctions/brezzidouglasfortinmarini/simplex/localcoefficients.hh>
#include <dune/localfunctions/brezzidouglasfortinmarini/simplex/localinterpolation.hh>


namespace Dune
{

  /**
   * \brief Brezzi-Douglas-Fortin-Marini finite elements for cubes
   *
   *  Brezzi-Douglas-Fortin-Marini (BDFM) finite elements are reduced
   *  Brezzi-Douglas-Marini (BDM) finite elements, where the order of
   *  the normal traces is lowered by one.
   *
   *  This implementation follows the description given by
   *
   *  M. W. Scroggs, P. D. Brubeck, J. P. Dean, J. S. Dokken, I. Marsden, N. Nobre, et al.
   *  DefElement: an encyclopedia of finite element definitions, 2020-2026,
   *  https://defelement.org
   *
   *  The BDFM0 element is identical to the RT0 element.
   *
   *  For further reading see
   *  Brezzi, Fortin "Mixed and Hybrid Finite Element Methods" (1991), Chapter III Section 3
   *
   * \ingroup BrezziDouglasFortinMarini
   *
   * \tparam D      Type to represent the field in the domain.
   * \tparam R      Type to represent the field in the range.
   * \tparam dim    dimension of the reference elements, must be >= 2.
   * \tparam order  order of the element, must be >= 1.
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMSimplexLocalFiniteElement
  {
    using LocalBasis          = BDFMSimplexLocalBasis<D, R, dim, order>;
    using LocalCoefficients   = BDFMSimplexLocalCoefficients<D, R, dim, order>;
    using LocalInterpolation  = BDFMSimplexLocalInterpolation<D, R, dim, order>;

  public:
    using Traits = LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation  >;

    //! \brief Standard constructor
    BDFMSimplexLocalFiniteElement () {}

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    BDFMSimplexLocalFiniteElement (int s)
      : basis( s ), interpolation( s )
    {}

    auto localBasis () const -> const LocalBasis& { return basis; }
    auto localCoefficients () const -> const LocalCoefficients& { return coefficients; }
    auto localInterpolation () const -> const LocalInterpolation& { return interpolation; }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const { return basis.size(); }
    static constexpr auto type () -> GeometryType { return GeometryTypes::cube(dim); }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
  };

} // namespace Dune

#undef BDFM_USE_RT0_BASIS
#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_BDFMSIMPLEX_HH
