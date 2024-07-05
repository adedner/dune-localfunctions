// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <dune/localfunctions/basix/localbasix.hh>
#include <dune/localfunctions/test/test-localfe.hh>

#include <basix/finite-element.h>
#include <basix/e-lagrange.h>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::BasixLocalFiniteElement<double, 0, 1> basixlfevertex{
    basix::element::create_lagrange<double>(basix::cell::type::point, 0,
    basix::element::lagrange_variant::equispaced, false)};
  TEST_FE(basixlfevertex);

  for (int degree = 1; degree < 10; ++degree)
  {
    Dune::BasixLocalFiniteElement<double, 1, 1> basixlfeline{
      basix::element::create_lagrange<double>(basix::cell::type::interval, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfeline);

    Dune::BasixLocalFiniteElement<double, 2, 1> basixlfetri{
      basix::element::create_lagrange<double>(basix::cell::type::triangle, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfetri);

    Dune::BasixLocalFiniteElement<double, 3, 1> basixlfetet{
      basix::element::create_lagrange<double>(basix::cell::type::tetrahedron, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfetet);

    Dune::BasixLocalFiniteElement<double, 2, 1> basixlfequad{
      basix::element::create_lagrange<double>(basix::cell::type::quadrilateral, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfequad);

    Dune::BasixLocalFiniteElement<double, 3, 1> basixlfehex{
      basix::element::create_lagrange<double>(basix::cell::type::hexahedron, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfehex);

    Dune::BasixLocalFiniteElement<double, 3, 1> basixlfeprism{
      basix::element::create_lagrange<double>(basix::cell::type::prism, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfeprism);

    Dune::BasixLocalFiniteElement<double, 3, 1> basixlfepyramid{
      basix::element::create_lagrange<double>(basix::cell::type::pyramid, degree,
      basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfepyramid);
  }

  return success ? 0 : 1;
}
