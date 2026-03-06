// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINISIMPLEX_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINISIMPLEX_HH

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

// implementation
#include <dune/localfunctions/brezzidouglasmarini/simplex/localbasis.hh>
#include <dune/localfunctions/brezzidouglasmarini/simplex/localcoefficients.hh>
#include <dune/localfunctions/brezzidouglasmarini/simplex/localinterpolation.hh>

namespace Dune
{
  /**
   * \brief Brezzi-Douglas-Marini local finite element for simplices
   *
   * \tparam D Number type to represent domain coordinates
   * \tparam R Number type to represent shape function values
   * \tparam dim Dimension of the reference elements, must be 2 or 3
   * \tparam order Polynomial order of the element
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BrezziDouglasMariniSimplexLocalFiniteElement
  {
    using LocalBasis          = BDMSimplexLocalBasis<D,R, dim, order>;
    using LocalCoefficients   = BDMSimplexLocalCoefficients<dim, order>;
    using LocalInterpolation  = BDMSimplexLocalInterpolation< D,R, dim, order>;

  public:
    using Traits = LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation  >;

    /** \brief Default constructor */
    BrezziDouglasMariniSimplexLocalFiniteElement()
    {}

    /**
     * \brief Constructor with a set of edge orientations
     *
     * \param s Bitfield of size 3 giving the orientations of the three element edges
     */
    BrezziDouglasMariniSimplexLocalFiniteElement(int s)
      : basis(s), interpolation(s)
    {}

    const LocalBasis& localBasis () const { return basis; }
    const LocalCoefficients& localCoefficients () const { return coefficients; }
    const LocalInterpolation& localInterpolation () const { return interpolation; }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const { return basis.size(); }
    static constexpr auto type () -> GeometryType { return GeometryTypes::simplex(dim); }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINISIMPLEX_HH
