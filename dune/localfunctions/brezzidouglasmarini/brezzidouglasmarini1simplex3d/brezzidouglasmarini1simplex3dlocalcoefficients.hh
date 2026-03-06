// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include "../../common/localkey.hh"

namespace Dune
{

  /**
   * \ingroup BrezziDouglasMariniImpl
   * \brief Layout map for Brezzi-Douglas-Marini-1 elements on tetrahedra
   *
   * \nosubgrouping
   * \implements Dune::LocalCoefficientsVirtualImp
   */
  class BDM1Simplex3DLocalCoefficients
  {

  public:
    //! \brief Standard constructor
    BDM1Simplex3DLocalCoefficients() : li(12)
    {
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

    //! \brief number of coefficients
    std::size_t size() const
    {
      return 12;
    }

    //! \brief get i'th index
    const LocalKey& localKey(std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };
} // end namespace Dune
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX3D_LOCALCOEFFICIENTS_HH
