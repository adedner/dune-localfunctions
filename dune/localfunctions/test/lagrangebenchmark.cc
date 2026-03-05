// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <array>
#include <ios>
#include <iostream>
#include <typeinfo>
#include <vector>

#include <benchmark/benchmark.h>

#include <dune/common/indices.hh>
#include <dune/common/timer.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex1.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex2.hh>
#include <dune/localfunctions/lagrange/lagrangesimplexold.hh>

#include <dune/localfunctions/test/test-localfe.hh>

using namespace Dune;

template <class LFE, class QuadRule, class Values, class Jacobians>
static void run_benchmark(const LFE& lfe, const QuadRule& quadRule, Values& values, Jacobians& jacobians)
{
  using Traits = typename LFE::Traits;
  using LB = typename Traits::LocalBasisType;
  constexpr int dim = LB::Traits::dimDomain;
  const LB& lb = lfe.localBasis();
  for (auto&& [x,w] : quadRule) {
    lb.evaluateFunction(x,values);
    lb.evaluateJacobian(x,jacobians);
    // lb.partial(std::array<unsigned int,dim>{std::min<unsigned int>(3,lb.order())}, x,values);
  }
}


// Static order local finite elements
template <class LFE>
static void BM_static(benchmark::State& state) {
  using Traits = typename LFE::Traits;
  using LB = typename Traits::LocalBasisType;
  LFE lfe{};
  std::vector<typename LB::Traits::RangeType> values(lfe.size());
  std::vector<typename LB::Traits::JacobianType> jacobians(lfe.size());
  const auto& quadRule = Dune::QuadratureRules<double,LB::Traits::dimDomain>::rule(lfe.type(), 6);
  for (auto _ : state)
    run_benchmark(lfe,quadRule,values,jacobians);
}


// Dynamic order local finite elements
template <class LFE>
static void BM_dynamic(benchmark::State& state) {
  using Traits = typename LFE::Traits;
  using LB = typename Traits::LocalBasisType;
  LFE lfe(state.range(0));
  std::vector<typename LB::Traits::RangeType> values(lfe.size());
  std::vector<typename LB::Traits::JacobianType> jacobians(lfe.size());
  const auto& quadRule = Dune::QuadratureRules<double,LB::Traits::dimDomain>::rule(lfe.type(), 6);
  for (auto _ : state)
    run_benchmark(lfe,quadRule,values,jacobians);
}


// Dynamic order local finite elements
template <class LFE>
static void BM_dynamic_old(benchmark::State& state) {
  using Traits = typename LFE::Traits;
  using LB = typename Traits::LocalBasisType;
  LFE lfe(GeometryTypes::simplex(LB::Traits::dimDomain),state.range(0));
  std::vector<typename LB::Traits::RangeType> values(lfe.size());
  std::vector<typename LB::Traits::JacobianType> jacobians(lfe.size());
  const auto& quadRule = Dune::QuadratureRules<double,LB::Traits::dimDomain>::rule(lfe.type(), 6);
  for (auto _ : state)
    run_benchmark(lfe,quadRule,values,jacobians);
}

BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,1, 1>>)->Name("static1/dim=1/1");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,1, 2>>)->Name("static1/dim=1/2");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,1, 4>>)->Name("static1/dim=1/4");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,1, 8>>)->Name("static1/dim=1/8");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,1,16>>)->Name("static1/dim=1/16");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,1, 1>>)->Name("static2/dim=1/1");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,1, 2>>)->Name("static2/dim=1/2");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,1, 4>>)->Name("static2/dim=1/4");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,1, 8>>)->Name("static2/dim=1/8");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,1,16>>)->Name("static2/dim=1/16");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,1, 1>>)->Name("static-old/dim=1/1");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,1, 2>>)->Name("static-old/dim=1/2");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,1, 4>>)->Name("static-old/dim=1/4");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,1, 8>>)->Name("static-old/dim=1/8");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,1,16>>)->Name("static-old/dim=1/16");
BENCHMARK(BM_dynamic<LagrangeSimplexLocalFiniteElement1<double,double,1>>)->Name("dynamic1/dim=1")->RangeMultiplier(2)->Range(1, 1<<4);
BENCHMARK(BM_dynamic<LagrangeSimplexLocalFiniteElement2<double,double,1>>)->Name("dynamic2/dim=1")->RangeMultiplier(2)->Range(1, 1<<4);
BENCHMARK(BM_dynamic_old<LagrangeLocalFiniteElement<EquidistantPointSet,1,double,double>>)->Name("dynamic-old/dim=1")->RangeMultiplier(2)->Range(1, 1<<4);

BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,2, 1>>)->Name("static1/dim=2/1");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,2, 2>>)->Name("static1/dim=2/2");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,2, 4>>)->Name("static1/dim=2/4");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,2, 8>>)->Name("static1/dim=2/8");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,2, 1>>)->Name("static2/dim=2/1");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,2, 2>>)->Name("static2/dim=2/2");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,2, 4>>)->Name("static2/dim=2/4");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,2, 8>>)->Name("static2/dim=2/8");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,2, 1>>)->Name("static-old/dim=2/1");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,2, 2>>)->Name("static-old/dim=2/2");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,2, 4>>)->Name("static-old/dim=2/4");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,2, 8>>)->Name("static-old/dim=2/8");
BENCHMARK(BM_dynamic<LagrangeSimplexLocalFiniteElement1<double,double,2>>)->Name("dynamic1/dim=2")->RangeMultiplier(2)->Range(1, 1<<3);
BENCHMARK(BM_dynamic<LagrangeSimplexLocalFiniteElement2<double,double,2>>)->Name("dynamic2/dim=2")->RangeMultiplier(2)->Range(1, 1<<3);
BENCHMARK(BM_dynamic_old<LagrangeLocalFiniteElement<EquidistantPointSet,2,double,double>>)->Name("dynamic-old/dim=2")->RangeMultiplier(2)->Range(1, 1<<3);

BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,3, 1>>)->Name("static1/dim=3/1");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,3, 2>>)->Name("static1/dim=3/2");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement1<double,double,3, 4>>)->Name("static1/dim=3/4");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,3, 1>>)->Name("static2/dim=3/1");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,3, 2>>)->Name("static2/dim=3/2");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElement2<double,double,3, 4>>)->Name("static2/dim=3/4");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,3, 1>>)->Name("static-old/dim=3/1");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,3, 2>>)->Name("static-old/dim=3/2");
BENCHMARK(BM_static<LagrangeSimplexLocalFiniteElementOld<double,double,3, 4>>)->Name("static-old/dim=3/4");
BENCHMARK(BM_dynamic<LagrangeSimplexLocalFiniteElement1<double,double,3>>)->Name("dynamic1/dim=3")->RangeMultiplier(2)->Range(1, 1<<2);
BENCHMARK(BM_dynamic<LagrangeSimplexLocalFiniteElement2<double,double,3>>)->Name("dynamic2/dim=3")->RangeMultiplier(2)->Range(1, 1<<2);
BENCHMARK(BM_dynamic_old<LagrangeLocalFiniteElement<EquidistantPointSet,3,double,double>>)->Name("dynamic-old/dim=3")->RangeMultiplier(2)->Range(1, 1<<2);


BENCHMARK_MAIN();





// template <unsigned int max_order>
// bool benchmark(std::ostream& out, int max_iter = 10)
// {
//   bool success = true;
//   Dune::Timer t;

//   Dune::Hybrid::forEach(std::index_sequence<1,2,3>{},[&](auto dim)
//   {
//     out << std::endl;
//     out << "dim = " << (int)(dim) << std::endl;
//     out << "-----------------" << std::endl;

//     { // 1. static Lagrange bases
//       auto lfe = Dune::unpackIntegerSequence([&](auto... i) {
//         return std::make_tuple(LagrangeSimplexLocalFiniteElement<double,double,dim,i+1>{}...);
//       }, std::make_index_sequence<max_order-1>{});

//       t.reset();
//       for (int iter = 0; iter < max_iter; ++iter) {
//         Hybrid::forEach(lfe,[&success](auto& pklfem) {
//           success &= testFE(pklfem);
//         });
//       }
//       out << "Time(static order) = " << t.elapsed()/max_iter << std::endl;
//     }

//     { // 1. static Lagrange bases
//       auto lfe = Dune::unpackIntegerSequence([&](auto... i) {
//         return std::make_tuple(LagrangeSimplexLocalFiniteElement1<double,double,dim,i+1>{}...);
//       }, std::make_index_sequence<max_order-1>{});

//       t.reset();
//       for (int iter = 0; iter < max_iter; ++iter) {
//         Hybrid::forEach(lfe,[&success](auto& pklfem) {
//           success &= testFE(pklfem);
//         });
//       }
//       out << "Time(static order 1) = " << t.elapsed()/max_iter << std::endl;
//     }

//     { // 1. static Lagrange bases
//       auto lfe = Dune::unpackIntegerSequence([&](auto... i) {
//         return std::make_tuple(LagrangeSimplexLocalFiniteElement2<double,double,dim,i+1>{}...);
//       }, std::make_index_sequence<max_order-1>{});

//       t.reset();
//       for (int iter = 0; iter < max_iter; ++iter) {
//         Hybrid::forEach(lfe,[&success](auto& pklfem) {
//           success &= testFE(pklfem);
//         });
//       }
//       out << "Time(static order 2) = " << t.elapsed()/max_iter << std::endl;
//     }

//     { // 4. old static Lagrange bases
//       auto lfe = Dune::unpackIntegerSequence([&](auto... i) {
//         return std::make_tuple(LagrangeSimplexLocalFiniteElementOld<double,double,dim,i+1>{}...);
//       }, std::make_index_sequence<max_order-1>{});

//       t.reset();
//       for (int iter = 0; iter < max_iter; ++iter) {
//         Hybrid::forEach(lfe,[&success](auto& pklfem) {
//           success &= testFE(pklfem);
//         });
//       }
//       out << "Time(old-static order) = " << t.elapsed()/max_iter << std::endl;
//     }

//     { // 2. dynamic version of the static Lagrange bases
//       std::vector<LagrangeSimplexLocalFiniteElement<double,double,dim>> lfe;
//       for (unsigned int i = 1; i < max_order; ++i) {
//         lfe.emplace_back(i);
//       }

//       t.reset();
//       for (int iter = 0; iter < max_iter; ++iter) {
//         for (auto& pklfem : lfe) {
//           success &= testFE(pklfem);
//         }
//       }
//       out << "Time(dynamic order) = " << t.elapsed()/max_iter << std::endl;
//     }

//     { // 2. dynamic version of the static Lagrange bases
//       std::vector<LagrangeSimplexLocalFiniteElement1<double,double,dim>> lfe;
//       for (unsigned int i = 1; i < max_order; ++i) {
//         lfe.emplace_back(i);
//       }

//       t.reset();
//       for (int iter = 0; iter < max_iter; ++iter) {
//         for (auto& pklfem : lfe) {
//           success &= testFE(pklfem);
//         }
//       }
//       out << "Time(dynamic order 1) = " << t.elapsed()/max_iter << std::endl;
//     }

//     { // 2. dynamic version of the static Lagrange bases
//       std::vector<LagrangeSimplexLocalFiniteElement2<double,double,dim>> lfe;
//       for (unsigned int i = 1; i < max_order; ++i) {
//         lfe.emplace_back(i);
//       }

//       t.reset();
//       for (int iter = 0; iter < max_iter; ++iter) {
//         for (auto& pklfem : lfe) {
//           success &= testFE(pklfem);
//         }
//       }
//       out << "Time(dynamic order 2) = " << t.elapsed()/max_iter << std::endl;
//     }

//   //   { // 3. Monomial based implementation
//   //     std::vector<LagrangeLocalFiniteElement<EquidistantPointSet,dim,double,double>> lfe;
//   //     for (unsigned int i = 1; i < max_order; ++i) {
//   //       lfe.emplace_back(Dune::GeometryTypes::simplex(dim),i);
//   //     }

//   //     t.reset();
//   //     for (int iter = 0; iter < max_iter; ++iter) {
//   //       for (auto& pklfem : lfe) {
//   //         success &= testFE(pklfem);
//   //       }
//   //     }
//   //     out << "Time(monomial-dynamic lagrange) = " << t.elapsed()/max_iter << std::endl;
//   //   }

//   });

//   return success;
// }

// int main (int argc, char *argv[])
// {
// #if __linux__ \
//   && (!defined __INTEL_COMPILER || __INTEL_COMPILER >= 1010) \
//   && (!defined __clang__)
//   feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
// #endif

//   bool success = true;

//   std::ofstream fout("benchmark.dat", std::ios_base::out);

//   success &= benchmark<8>(fout, 10);

//   fout.close();

//   return success ? 0 : 1;
// }
