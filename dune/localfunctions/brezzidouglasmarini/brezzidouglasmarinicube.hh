// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINICUBE_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINICUBE_HH

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

// implementation
#include <dune/localfunctions/brezzidouglasmarini/cube/localbasis.hh>
#include <dune/localfunctions/brezzidouglasmarini/cube/localcoefficients.hh>
#include <dune/localfunctions/brezzidouglasmarini/cube/localinterpolation.hh>

namespace Dune
{
  /**
   * \brief Brezzi-Douglas-Marini local finite element for cubes
   *
   * \tparam D Number type to represent domain coordinates
   * \tparam R Number type to represent shape function values
   * \tparam dim Dimension of the reference elements, must be 2 or 3
   * \tparam order Polynomial order of the element
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BrezziDouglasMariniCubeLocalFiniteElement
  {
    using LocalBasis          = BDMCubeLocalBasis<D,R, dim, order>;
    using LocalCoefficients   = BDMCubeLocalCoefficients<dim, order>;
    using LocalInterpolation  = BDMCubeLocalInterpolation< D,R, dim, order>;

  public:
    using Traits = LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation  >;

    /** \brief Default constructor */
    BrezziDouglasMariniCubeLocalFiniteElement()
    {}

    /**
     * \brief Constructor with a set of edge or face orientations
     *
     * \param s Bitfield of size 4 or 6 giving the orientations of the four element edges or 6 element faces
     */
    BrezziDouglasMariniCubeLocalFiniteElement(int s)
      : basis(s), interpolation(s)
    {}

    const LocalBasis& localBasis () const { return basis; }
    const LocalCoefficients& localCoefficients () const { return coefficients; }
    const LocalInterpolation& localInterpolation () const { return interpolation; }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const { return basis.size(); }
    static constexpr auto type () -> GeometryType { return GeometryTypes::cube(dim); }

  private:
    LocalBasis basis;
    LocalCoefficients coefficients;
    LocalInterpolation interpolation;
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINICUBE_HH
