// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <array>
#include <limits>
#include <vector>

#include <dune/common/tensor.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/utility/hessian.hh>

using namespace Dune;

template <class Domain, class Range>
struct DerivativeTraits;

template <class Domain>
struct DerivativeTraits<Domain, double>
{
  using type = Domain;
};

template <int d, int n>
struct DerivativeTraits<FieldVector<double,d>, FieldVector<double,n>>
{
  using type = FieldMatrix<double,n,d>;
};

template <int d, int m, int n>
struct DerivativeTraits<FieldVector<double,d>, FieldMatrix<double,m,n>>
{
  using type = Tensor<double,m,n,d>;
};

template <std::size_t N1, std::size_t N2, std::size_t N3>
bool compare (std::size_t i3, const Tensor<double,N1,N2,N3>& H, const FieldMatrix<double,int(N1),int(N2)>& DJ)
{
  if (i3 >= N3)
    return false;

  using std::sqrt;
  double tol = 10*sqrt(std::numeric_limits<double>::epsilon());
  for (std::size_t i1 = 0; i1 < N1; ++i1)
    for (std::size_t i2 = 0; i2 < N2; ++i2)
      if (abs(H[std::array{i1,i2,i3}] - DJ[std::array{i1,i2}]) > tol) {
#ifndef NDEBUG
        std::cout << "error = " << abs(H[std::array{i1,i2,i3}] - DJ[std::array{i1,i2}]) << std::endl;
        std::cout << "tol = " << tol << std::endl;
#endif
        return false;
      }

  return true;
}

template <class LFE>
TestSuite testHessian (const LFE& lfe, std::string name)
{
  TestSuite subTest(name);

  using Traits = typename LFE::Traits::LocalBasisType::Traits;
  using Domain = typename Traits::DomainType;
  using Range = typename Traits::RangeType;
  using Jacobian = typename DerivativeTraits<Domain,Range>::type;
  using Hessian = typename DerivativeTraits<Domain,Jacobian>::type;

  std::vector<Hessian> shapeHessians;
  std::array<std::vector<typename Traits::JacobianType>, 2*Traits::dimDomain> shapeJacobians{};

  auto eps = sqrt(std::numeric_limits<typename Traits::DomainFieldType>::epsilon());
  for (const auto& q : QuadratureRules<double,Traits::dimDomain>::rule(lfe.type(),3))
  {
    evaluateHessian(lfe.localBasis(), q.position(), shapeHessians);

    for (int d = 0; d < Traits::dimDomain; ++d)
    {
      auto x = q.position();
      x[d] += eps;
      lfe.localBasis().evaluateJacobian(x, shapeJacobians[2*d+0]);
      x[d] -= 2*eps;
      lfe.localBasis().evaluateJacobian(x, shapeJacobians[2*d+1]);
    }

    subTest.require(shapeHessians.size() == lfe.size());
    for (std::size_t i = 0; i < lfe.size(); ++i)
    {
      for (int d = 0; d < Traits::dimDomain; ++d)
      {
        auto DJ = (shapeJacobians[2*d+0][i] - shapeJacobians[2*d+1][i])/(2*eps);
        subTest.check(compare(d,shapeHessians[i], DJ));
      }
    }
  }

  return subTest;
}

int main ()
{
  TestSuite test;

  { // k = 1
    auto p11d = LagrangeSimplexLocalFiniteElement<double,double,1,1>{};
    test.subTest(testHessian(p11d, "p11d"));

    auto p12d = LagrangeSimplexLocalFiniteElement<double,double,2,1>{};
    test.subTest(testHessian(p12d, "p12d"));

    auto p13d = LagrangeSimplexLocalFiniteElement<double,double,3,1>{};
    test.subTest(testHessian(p13d, "p13d"));
  }

  { // k = 2
    // auto p21d = LagrangeSimplexLocalFiniteElement<double,double,1,2>{};
    // test.subTest(testHessian(p21d, "p21d"));

    auto p22d = LagrangeSimplexLocalFiniteElement<double,double,2,2>{};
    test.subTest(testHessian(p22d, "p22d"));

    // auto p23d = LagrangeSimplexLocalFiniteElement<double,double,3,2>{};
    // test.subTest(testHessian(p23d, "p23d"));
  }

  { // k = 3
    // auto p21d = LagrangeSimplexLocalFiniteElement<double,double,1,2>{};
    // test.subTest(testHessian(p21d, "p21d"));

    auto p32d = LagrangeSimplexLocalFiniteElement<double,double,2,3>{};
    test.subTest(testHessian(p32d, "p32d"));

    // auto p23d = LagrangeSimplexLocalFiniteElement<double,double,3,2>{};
    // test.subTest(testHessian(p23d, "p23d"));
  }

  return test.exit();
}