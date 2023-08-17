// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <dune/localfunctions/hermite.hh>
#include <dune/localfunctions/test/c1_test-localfe.hh>

using namespace Dune;

int main(int argc, char **argv)
{
  bool success = true; // captured by Macros
  HermiteLocalFiniteElement<double, double, 1> LFE_1d;
  // FE, DisableSubTests, max order for differentiabilitytest (<=2)
  TEST_FE3(LFE_1d, DisableVirtualInterface, 2);
  HermiteLocalFiniteElement<double, double, 2> LFE_2d;
  TEST_FE3(LFE_2d, DisableVirtualInterface, 2);
  HermiteLocalFiniteElement<double, double, 3> LFE_3d;
  TEST_FE3(LFE_3d, DisableVirtualInterface, 2);

  return success ? 0 : 1;
}
