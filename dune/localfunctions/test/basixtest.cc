// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <dune/localfunctions/basix/globalbasix.hh>
#include <dune/localfunctions/basix/localbasix.hh>
#include <dune/localfunctions/test/test-localfe.hh>

#include <basix/finite-element.h>
#include <basix/e-lagrange.h>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::BasixLocalFiniteElement<0, Dune::RangeClass::scalar> basixlfevertex{
    basix::element::create_lagrange<double>(basix::cell::type::point, 0,
    basix::element::lagrange_variant::equispaced, false)};
  TEST_FE(basixlfevertex);

  for (int degree = 1; degree < 2; ++degree)
  {
    std::cout << "degree " << degree << std::endl;
    Dune::BasixLocalFiniteElement<1, Dune::RangeClass::scalar> basixlfeline{
      basix::element::create_lagrange<double>(basix::cell::type::interval, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfeline);

    Dune::BasixLocalFiniteElement<2, Dune::RangeClass::scalar> basixlfetri{
      basix::element::create_lagrange<double>(basix::cell::type::triangle, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfetri);

    Dune::BasixLocalFiniteElement<3, Dune::RangeClass::scalar> basixlfetet{
      basix::element::create_lagrange<double>(basix::cell::type::tetrahedron, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfetet);

    Dune::BasixLocalFiniteElement<2, Dune::RangeClass::scalar> basixlfequad{
      basix::element::create_lagrange<double>(basix::cell::type::quadrilateral, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfequad);

    Dune::BasixLocalFiniteElement<3, Dune::RangeClass::scalar> basixlfehex{
      basix::element::create_lagrange<double>(basix::cell::type::hexahedron, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfehex);

    Dune::BasixLocalFiniteElement<3, Dune::RangeClass::scalar> basixlfeprism{
      basix::element::create_lagrange<double>(basix::cell::type::prism, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfeprism);

    Dune::BasixLocalFiniteElement<3, Dune::RangeClass::scalar> basixlfepyramid{
      basix::element::create_lagrange<double>(basix::cell::type::pyramid, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfepyramid);
  }


  auto refElem = Dune::referenceElement<double,2>(Dune::GeometryTypes::triangle);
  auto geometry = refElem.geometry<0>(0);
  using Geometry = decltype(geometry);

  Dune::BasixFiniteElement<Geometry, Dune::RangeClass::scalar> basixfetri{
    basix::element::create_lagrange<double>(basix::cell::type::triangle, 3,
    basix::element::lagrange_variant::equispaced, false)};
  basixfetri.bind(geometry);
  TEST_FE(basixfetri);

  return success ? 0 : 1;
}
