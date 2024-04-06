// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_META_PRODUCT_INTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_META_PRODUCT_INTERPOLATION_HH

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

namespace Dune {

template<class LI1, class LI2, class LB>
class PrismaticProductLocalInterpolation
{
  const LI1* li1_;
  const LI2* li2_;

public:
  //! Construct a PrismaticProductLocalInterpolation
  PrismaticProductLocalInterpolation(const LI1& li1, const LI2& li2)
    : li1_(&li1)
    , li2_(&li2)
  {}

  //! Determine coefficients interpolating a given function
  /**
   * \param f   An object supporting the expression \c y = f(x).
   * \param out Vector where to store the interpolated coefficients.
   */
  template<typename F, typename C>
  void interpolate(const F& f, std::vector<C>& out) const
  {
    // TODO
  }
};

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_PRODUCT_INTERPOLATION_HH
