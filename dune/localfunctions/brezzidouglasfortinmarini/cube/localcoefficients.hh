// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/localfunctions/common/localkey.hh>

#ifdef BDFM_USE_RT0_BASIS
#include <dune/localfunctions/raviartthomas/raviartthomas0cube3d/raviartthomas0cube3dall.hh>
#endif


namespace Dune
{

  /**
   * \ingroup BrezziDouglasFortinMariniImpl
   *
   * \brief Layout map for Brezzi-Douglas-Fortin-Marini elements on cubes.
   *
   * \tparam D      Type of represent the field in the domain.
   * \tparam R      Type of represent the field in the domain.
   * \tparam dim    dimension of the reference element, must be >= 2.
   * \tparam order  order of the element, must be >= 1.
   *
   * \nosubgrouping
   * \implements Dune::LocalCoefficientsVirtualImp
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMCubeLocalCoefficients;

#ifdef BDFM_USE_RT0_BASIS
  template<class D, class R>
  class BDFMCubeLocalCoefficients<D, R, 3, 0>
    : public RT0Cube3DLocalCoefficients
  {
  public:
    BDFMCubeLocalCoefficients() : RT0Cube3DLocalCoefficients () {}
  };
#endif

  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMCubeLocalCoefficients
  {
    static constexpr unsigned int interiorDofs = dim*binomial(dim+order-2, order-2);
    static constexpr unsigned int faceDofs     = binomial(dim+order-2, order-1);

    static constexpr std::size_t numFaces = 2*dim;
    static constexpr std::size_t numDofs  = numFaces*faceDofs + interiorDofs;

    public:
      //! \brief Standard constructor
      BDFMCubeLocalCoefficients () : li(numDofs)
      {
        for (auto j : range(numFaces))
          for (auto i : range(faceDofs))
            li[j*faceDofs + i] = LocalKey(j, 1, i);

        for (auto i : range(interiorDofs))
          li[numFaces*faceDofs + i] = LocalKey(0, 0, i);
      }

      //! \brief number of coefficients
      std::size_t size () const { return numDofs; }

      //! \brief geth i'th index
      auto localKey (std::size_t i) const -> const LocalKey& { return li[i]; }

    private:
      std::vector<LocalKey> li;
  };

  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr unsigned int BDFMCubeLocalCoefficients<D, R, dim, order>::interiorDofs;

  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr unsigned int BDFMCubeLocalCoefficients<D, R, dim, order>::faceDofs;

  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr std::size_t BDFMCubeLocalCoefficients<D, R, dim, order>::numFaces;

  // template<class D, class R, unsigned int dim, unsigned int order>
  // constexpr std::size_t BDFMCubeLocalCoefficients<D, R, dim, order>::numDofs;

#ifndef BDFM_USE_RT0_BASIS
#ifndef DOXYGEN
  template<class D, class R, unsigned int dim>
  class BDFMCubeLocalCoefficients<D, R, dim, 0>
  {
    static_assert( AlwaysFalse<D>::value,
                   "`BDFMCubeLocalCoefficients` not defined for order 0." );
  };
#endif // #ifndef DOXYGEN
#endif

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_CUBE_LOCALCOEFFICIENTS_HH
