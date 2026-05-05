// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_SIMPLEX_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_SIMPLEX_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief Layout map for Brezzi-Douglas-Marini elements on simplices.
   *
   * \nosubgrouping
   * \implements Dune::LocalCoefficientsVirtualImp
   */
  template<unsigned int dim, unsigned int order>
  class BDMSimplexLocalCoefficients
  {
    static_assert( order == 1 || order == 2, "`BDMSimplexLocalCoefficients` only defined for order 1 and 2.");

    static constexpr std::size_t numFaces = dim+1;
    static constexpr std::size_t numDofs  = order == 1 ? numFaces * dim : numFaces * (dim+2);

    public:
      //! \brief Standard constructor
      BDMSimplexLocalCoefficients () : li(numDofs)
      {
        // triangle dofs
        if constexpr (dim == 2 )
        {
          if constexpr ( order == 1 )
          {
            assert( numDofs == 6 );
            for (std::size_t i=0; i<3; i++)
            {
              li[i] = LocalKey(i,1,0);
              li[3 + i] = LocalKey(i,1,1);
            }
          }

          if constexpr ( order == 2 )
          {
            assert( numDofs == 12 );
            for (std::size_t i = 0; i < 3; ++i)
            {
              li[3 * i] = LocalKey(i,1,0);
              li[3 * i + 1] = LocalKey(i,1,1);
              li[3 * i + 2] = LocalKey(i,1,2);
            }

            // last DOFs are associated with the cell (codim=0)
            li[9]  = LocalKey(0,0,0);
            li[10] = LocalKey(0,0,1);
            li[11] = LocalKey(0,0,2);
          }
        }

        // tetrahedron dofs
        if constexpr ( dim == 3 )
        {
          if constexpr ( order == 1 )
          {
            assert( numDofs == 12 );
            // LocalKey (sub-entity, codim, number on sub-entity)
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
          }
        }
      }

      //! \brief number of coefficients
      std::size_t size () const { return numDofs; }

      //! \brief geth i'th index
      const LocalKey& localKey (std::size_t i) const { return li[i]; }

    private:
      std::vector<LocalKey> li;
  };

  template<unsigned int dim, unsigned int order>
  constexpr std::size_t BDMSimplexLocalCoefficients<dim, order>::numFaces;

  template<unsigned int dim, unsigned int order>
  constexpr std::size_t BDMSimplexLocalCoefficients<dim, order>::numFaces;

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_SIMPLEX_LOCALCOEFFICIENTS_HH
