// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <iostream>

#include <dune/localfunctions/brezzidouglasfortinmarini/bdfmcube.hh>
#include <dune/localfunctions/brezzidouglasfortinmarini/bdfmsimplex.hh>

#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv)
{
  bool success = true;

  // cube tests
  Dune::BDFMCubeLocalFiniteElement<double,double, 2, 1> bdfm1cube2dlfem(1);
  TEST_FE(bdfm1cube2dlfem);

  Dune::BDFMCubeLocalFiniteElement<double,double, 2, 2> bdfm2cube2dlfem(1);
  TEST_FE(bdfm2cube2dlfem);

  Dune::BDFMCubeLocalFiniteElement<double,double, 2, 3> bdfm3cube2dlfem(1);
  TEST_FE(bdfm3cube2dlfem);

  Dune::BDFMCubeLocalFiniteElement<double,double, 3, 0> bdfm0cube3dlfem(1);
  TEST_FE(bdfm0cube3dlfem);

  // simplex tests
  Dune::BDFMSimplexLocalFiniteElement<double,double, 2, 0> bdfm0simplex2dlfem(1);
  TEST_FE(bdfm0simplex2dlfem);

  Dune::BDFMSimplexLocalFiniteElement<double,double, 2, 1> bdfm1simplex2dlfem(1);
  TEST_FE(bdfm1simplex2dlfem);

  Dune::BDFMSimplexLocalFiniteElement<double,double, 3, 0> bdfm0simplex3dlfem(1);
  TEST_FE(bdfm0simplex3dlfem);

  Dune::BDFMSimplexLocalFiniteElement<double,double, 3, 1> bdfm1simplex3dlfem(1);
  TEST_FE(bdfm1simplex3dlfem);

  return success ? 0 : 1;
}
