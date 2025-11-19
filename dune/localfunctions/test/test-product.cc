// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/meta/product.hh>

#include <dune/localfunctions/test/test-localfe.hh>

template <unsigned int k, class TimeLFE>
bool test_1d(const TimeLFE& timeLFE)
{
  bool success = true;

  std::cout << "test_1d(k=" << k << ")" << std::endl;

  using SpaceLFE = Dune::LagrangeSimplexLocalFiniteElement<double,double,1,k>;
  SpaceLFE spaceLFE;

  Dune::PrismaticProduct lfe(spaceLFE, timeLFE);
  TEST_FE(lfe);

  return success;
}

template <unsigned int k, class TimeLFE>
bool test_2d(const TimeLFE& timeLFE)
{
  bool success = true;

  std::cout << "test_2d(k=" << k << ")" << std::endl;

  { // simplex
    using SpaceLFE = Dune::LagrangeSimplexLocalFiniteElement<double,double,2,k>;
    SpaceLFE spaceLFE;

    Dune::PrismaticProduct lfe(spaceLFE, timeLFE);
    TEST_FE(lfe);
  }

  { // cube
    using SpaceLFE = Dune::LagrangeCubeLocalFiniteElement<double,double,2,k>;
    SpaceLFE spaceLFE;

    Dune::PrismaticProduct lfe(spaceLFE, timeLFE);
    TEST_FE(lfe);
  }

  return success;
}

template <unsigned int k, class TimeLFE>
bool test_3d(const TimeLFE& timeLFE)
{
  bool success = true;

  std::cout << "test_3d(k=" << k << ")" << std::endl;

  { // simplex
    using SpaceLFE = Dune::LagrangeSimplexLocalFiniteElement<double,double,3,k>;
    SpaceLFE spaceLFE;

    Dune::PrismaticProduct lfe(spaceLFE, timeLFE);
    TEST_FE(lfe);
  }

  { // cube
    using SpaceLFE = Dune::LagrangeCubeLocalFiniteElement<double,double,3,k>;
    SpaceLFE spaceLFE;

    Dune::PrismaticProduct lfe(spaceLFE, timeLFE);
    TEST_FE(lfe);
  }

  if (k < 3) { // prism
    using SpaceLFE = Dune::LagrangePrismLocalFiniteElement<double,double,k>;
    SpaceLFE spaceLFE;

    Dune::PrismaticProduct lfe(spaceLFE, timeLFE);
    TEST_FE(lfe);
  } else {
    using SpaceLFE1 = Dune::LagrangeSimplexLocalFiniteElement<double,double,2,k>;
    using SpaceLFE2 = Dune::LagrangeSimplexLocalFiniteElement<double,double,1,k>;
    SpaceLFE1 spaceLFE1;
    SpaceLFE2 spaceLFE2;

    Dune::PrismaticProduct spaceLFE(spaceLFE1, spaceLFE2);
    Dune::PrismaticProduct lfe(spaceLFE, timeLFE);
    TEST_FE(lfe);
  }

  // Pyramid shapefunctions are not differentiable on the plane where xi[0]=xi[1].
  // TODO: Implement a restricted test
  // { // pyramid
  //   using SpaceLFE = Dune::LagrangePyramidLocalFiniteElement<double,double,k>;
  //   SpaceLFE spaceLFE;

  //   Dune::PrismaticProduct lfe(spaceLFE, timeLFE);
  //   TEST_FE(lfe);
  // }

  return success;
}


template <unsigned int k, unsigned int l>
bool test_generic()
{
  bool success = true;

  std::cout << "test_generic(k=" << k << ", l=" << l << ")" << std::endl;

  { // simplex x cube
    using LFE1 = Dune::LagrangeSimplexLocalFiniteElement<double,double,2,k>;
    LFE1 lfe1;

    using LFE2 = Dune::LagrangeCubeLocalFiniteElement<double,double,2,l>;
    LFE2 lfe2;

    Dune::PrismaticProduct lfe(lfe1, lfe2);
    TEST_FE(lfe);
  }

  { // cube x cube
    using LFE1 = Dune::LagrangeCubeLocalFiniteElement<double,double,2,k>;
    LFE1 lfe1;

    using LFE2 = Dune::LagrangeCubeLocalFiniteElement<double,double,2,l>;
    LFE2 lfe2;

    Dune::PrismaticProduct lfe(lfe1, lfe2);
    TEST_FE(lfe);
  }

  return success;
}

int main(int argc, char** argv)
{
  bool success = true;

  Dune::Hybrid::forEach(Dune::StaticIntegralRange<unsigned int, 4, 1>{}, [&](auto k) {
    Dune::Hybrid::forEach(Dune::StaticIntegralRange<unsigned int, 4, 1>{}, [&](auto l) {
      using TimeLFE = Dune::LagrangeSimplexLocalFiniteElement<double,double,1,(unsigned int)(l)>;
      TimeLFE timeLFE;

      success |= test_1d<(unsigned int)(k)>(timeLFE);
      success |= test_2d<(unsigned int)(k)>(timeLFE);
      success |= test_3d<(unsigned int)(k)>(timeLFE);

      success |= test_generic<(unsigned int)(k), (unsigned int)(l)>();
    });
  });

  return success ? 0 : 1;
}
