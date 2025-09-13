// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <ios>
#include <iostream>
#include <typeinfo>
#include <fenv.h>

#include <dune/common/indices.hh>
#include <dune/common/timer.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/p0.hh>
// #include <dune/localfunctions/lagrange/pq22d.hh>
// #include <dune/localfunctions/lagrange/dynlagrangecache.hh>
// #include <dune/localfunctions/lagrange/dynlagrangesimplex.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/lagrangesimplexold.hh>

#include <dune/localfunctions/test/test-localfe.hh>

using namespace Dune;


template <unsigned int max_order>
bool benchmark(std::ostream& out, int max_iter = 10)
{
  bool success = true;
  Dune::Timer t;

  Dune::Hybrid::forEach(std::index_sequence<1,2,3>{},[&](auto dim)
  {
    out << std::endl;
    out << "dim = " << (int)(dim) << std::endl;
    out << "-----------------" << std::endl;

    { // 1. static Lagrange bases
      auto lfe = Dune::unpackIntegerSequence([&](auto... i) {
        return std::make_tuple(LagrangeSimplexLocalFiniteElement<double,double,dim,i+1>{}...);
      }, std::make_index_sequence<max_order-1>{});

      t.reset();
      for (int iter = 0; iter < max_iter; ++iter) {
        Hybrid::forEach(lfe,[&success](auto& pklfem) {
          success &= testFE(pklfem);
        });
      }
      out << "Time(static order) = " << t.elapsed()/max_iter << std::endl;
    }

    { // 4. old static Lagrange bases
      auto lfe = Dune::unpackIntegerSequence([&](auto... i) {
        return std::make_tuple(LagrangeSimplexLocalFiniteElementOld<double,double,dim,i+1>{}...);
      }, std::make_index_sequence<max_order-1>{});

      t.reset();
      for (int iter = 0; iter < max_iter; ++iter) {
        Hybrid::forEach(lfe,[&success](auto& pklfem) {
          success &= testFE(pklfem);
        });
      }
      out << "Time(old-static order) = " << t.elapsed()/max_iter << std::endl;
    }

    { // 2. dynamic version of the static Lagrange bases
      std::vector<LagrangeSimplexLocalFiniteElement<double,double,dim>> lfe;
      for (unsigned int i = 1; i < max_order; ++i) {
        lfe.emplace_back(i);
      }

      t.reset();
      for (int iter = 0; iter < max_iter; ++iter) {
        for (auto& pklfem : lfe) {
          success &= testFE(pklfem);
        }
      }
      out << "Time(dynamic order) = " << t.elapsed()/max_iter << std::endl;
    }

    { // 3. Monomial based implementation
      std::vector<LagrangeLocalFiniteElement<EquidistantPointSet,dim,double,double>> lfe;
      for (unsigned int i = 1; i < max_order; ++i) {
        lfe.emplace_back(Dune::GeometryTypes::simplex(dim),i);
      }

      t.reset();
      for (int iter = 0; iter < max_iter; ++iter) {
        for (auto& pklfem : lfe) {
          success &= testFE(pklfem);
        }
      }
      out << "Time(monomial-dynamic lagrange) = " << t.elapsed()/max_iter << std::endl;
    }
  });

  return success;
}

int main (int argc, char *argv[])
{
#if __linux__ \
  && (!defined __INTEL_COMPILER || __INTEL_COMPILER >= 1010) \
  && (!defined __clang__)
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  bool success = true;

  std::ofstream fout("benchmark.dat", std::ios_base::out);
  Dune::Timer t;

  success &= benchmark<8>(fout, 100);

  fout.close();

  return success ? 0 : 1;
}
