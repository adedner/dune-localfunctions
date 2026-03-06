// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_SIMPLEX_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_SIMPLEX_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**
   * \ingroup BrezziDouglasFortinMariniImpl
   *
   * \brief Layout map for Brezzi-Douglas-Fortin-Marini elements on simplices.
   *
   * \tparam D      Type of represent the field in the domain.
   * \tparam R      Type of represent the field in the domain.
   * \tparam dim    dimension of the reference element, must be >= 2.
   * \tparam order  order of the element, must be 0 or 1.
   *
   * \nosubgrouping
   * \implements Dune::LocalCoefficientsVirtualImp
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BDFMSimplexLocalCoefficients
  {
    static constexpr unsigned int interiorDofs = order*(3*(dim-1));
    static constexpr unsigned int faceDofs     = order == 0 ? 1 : dim;

    static constexpr std::size_t numFaces = dim+1;
    static constexpr std::size_t numDofs  = numFaces*faceDofs + interiorDofs;

    public:
      //! \brief Standard constructor
      BDFMSimplexLocalCoefficients () : li(numDofs)
      {
        if constexpr ( order == 0 )
        {
          for( size_t i=0; i<numFaces; ++i )
            li[ i ] = LocalKey( i, 1, 0);
        }

        if constexpr ( order == 1 && dim == 2)
        {
          assert( li.size() == 9 );
          li[ 0] = LocalKey( 0, 1, 0);
          li[ 1] = LocalKey( 0, 1, 1);
          li[ 2] = LocalKey( 1, 1, 0);
          li[ 3] = LocalKey( 1, 1, 1);
          li[ 4] = LocalKey( 2, 1, 0);
          li[ 5] = LocalKey( 2, 1, 1);
          li[ 6] = LocalKey( 0, 0, 0);
          li[ 7] = LocalKey( 0, 0, 1);
          li[ 8] = LocalKey( 0, 0, 2);
        }

        if constexpr ( order == 1 && dim == 3 )
        {
          li[ 0] = LocalKey( 0, 1, 0);
          li[ 1] = LocalKey( 0, 1, 1);
          li[ 2] = LocalKey( 0, 1, 2);
          li[ 3] = LocalKey( 1, 1, 0);
          li[ 4] = LocalKey( 1, 1, 1);
          li[ 5] = LocalKey( 1, 1, 2);
          li[ 6] = LocalKey( 2, 1, 0);
          li[ 7] = LocalKey( 2, 1, 1);
          li[ 8] = LocalKey( 2, 1, 2);
          li[ 9] = LocalKey( 3, 1, 0);
          li[10] = LocalKey( 3, 1, 1);
          li[11] = LocalKey( 3, 1, 2);
          li[12] = LocalKey( 0, 0, 0);
          li[13] = LocalKey( 0, 0, 1);
          li[14] = LocalKey( 0, 0, 2);
          li[15] = LocalKey( 0, 0, 3);
          li[16] = LocalKey( 0, 0, 4);
          li[17] = LocalKey( 0, 0, 5);
        }

        if constexpr ( order >= 2 )
        {
          DUNE_THROW(NotImplemented, "order >= 2 not implemented yet");
        }
      }

      //! \brief number of coefficients
      std::size_t size () const { return numDofs; }

      //! \brief geth i'th index
      auto localKey (std::size_t i) const -> const LocalKey& { return li[i]; }

    private:
      std::vector<LocalKey> li;
  };

  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr unsigned int BDFMSimplexLocalCoefficients<D, R, dim, order>::interiorDofs;

  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr unsigned int BDFMSimplexLocalCoefficients<D, R, dim, order>::faceDofs;

  template<class D, class R, unsigned int dim, unsigned int order>
  constexpr std::size_t BDFMSimplexLocalCoefficients<D, R, dim, order>::numFaces;


} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASFORTINMARINI_SIMPLEX_LOCALCOEFFICIENTS_HH
