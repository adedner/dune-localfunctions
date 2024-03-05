// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include "config.h"

#include <dune/localfunctions/arnoldwinther.hh>
#include <dune/localfunctions/test/c1_test-localfe.hh>

using namespace Dune;

int main(int argc, char **argv)
{
  bool success = true; // captured by Macros
  ArnoldWintherLocalFiniteElement<double, double> LFE;
  // FE, DisableSubTests, max order for differentiabilitytest (<=2)
  // std::cout <<std::setprecision(2)<<std::fixed<< getDualityMatrix(LFE)<<std::endl;
  TEST_FE3(LFE, 0, 1);
  return success ? 0 : 1;
}
