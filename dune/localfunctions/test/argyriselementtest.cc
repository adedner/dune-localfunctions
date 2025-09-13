// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include "config.h"

#include <dune/localfunctions/argyris.hh>
#include <dune/localfunctions/test/c1_test-localfe.hh>

using namespace Dune;

int main(int argc, char **argv)
{
  bool success = true; // captured by Macros
  ArgyrisLocalFiniteElement<double, double> LFE_1d;
  // FE, DisableSubTests, max order for differentiabilitytest (<=2)
  TEST_FE3(LFE_1d, DisableVirtualInterface, 2);

  for (std::size_t s = 0; s < 8; ++s)
  {
    ArgyrisLocalFiniteElement<double, double> orientendArgyrisLFE(s);
    TEST_FE3(orientendArgyrisLFE, DisableVirtualInterface, 2);
  }
  return success ? 0 : 1;
}
