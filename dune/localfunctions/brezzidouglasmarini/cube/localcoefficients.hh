// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_CUBE_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_CUBE_LOCALCOEFFICIENTS_HH

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
  class BDMCubeLocalCoefficients
  {
    static constexpr std::size_t numFaces = 2*dim;
    static constexpr std::size_t numDofs  = order == 1 ? numFaces * dim : numFaces * (dim+1) + dim;

    public:
      //! \brief Standard constructor
      BDMCubeLocalCoefficients () : li(numDofs)
      {
        // quadrilateral dofs
        if constexpr (dim == 2 )
        {
          if constexpr ( order == 1 )
          {
            assert( numDofs == 8 );
            for (std::size_t i = 0; i < 4; ++i)
            {
              li[2*i] = LocalKey(i,1,0);
              li[2*i + 1] = LocalKey(i,1,1);
            }
          }

          if constexpr ( order == 2 )
          {
            assert( numDofs == 14 );
            for (std::size_t i = 0; i < 4; ++i)
            {
              li[3 * i] = LocalKey(i,1,0);
              li[3 * i + 1] = LocalKey(i,1,1);
              li[3 * i + 2] = LocalKey(i,1,2);
            }
            li[12] = LocalKey(0,0,0);
            li[13] = LocalKey(0,0,1);
          }
        }

        // hexahedron dofs
        if constexpr ( dim == 3 )
        {
          if constexpr ( order == 1 )
          {
            assert( numDofs == 18 );
            for (std::size_t i = 0; i < 6; ++i)
            {
              li[i] = LocalKey(i,1,0);
              li[i + 6] = LocalKey(i,1,1);
              li[i + 12] = LocalKey(i,1,2);
            }
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
  constexpr std::size_t BDMCubeLocalCoefficients<dim, order>::numFaces;

  template<unsigned int dim, unsigned int order>
  constexpr std::size_t BDMCubeLocalCoefficients<dim, order>::numDofs;

#ifndef DOXYGEN
  template<unsigned int dim>
  class BDMCubeLocalCoefficients<dim, 0>
  {
    static_assert( AlwaysFalse<double>::value,
                   "`BDMCubeLocalCoefficients` not defined for order 0." );
  };
#endif // #ifndef DOXYGEN

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_CUBE_LOCALCOEFFICIENTS_HH
