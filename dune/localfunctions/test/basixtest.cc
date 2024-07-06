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

  // Dune::BasixLocalFiniteElement<0, Dune::RangeClass::scalar> basixlfevertex{
  //   basix::element::create_lagrange<double>(basix::cell::type::point, 0,
  //   basix::element::lagrange_variant::equispaced, false)};
  // TEST_FE(basixlfevertex);

  for (int degree = 1; degree < 5; ++degree)
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


  std::vector<std::pair<std::string, basix::cell::type>> cells{
    {"point", basix::cell::type::point},
    {"interval", basix::cell::type::interval},
    {"triangle", basix::cell::type::triangle},
    {"tetrahedron", basix::cell::type::tetrahedron},
    {"quadrilateral", basix::cell::type::quadrilateral},
    {"hexahedron", basix::cell::type::hexahedron},
    {"prism", basix::cell::type::prism},
    {"pyramid", basix::cell::type::pyramid}
  };

  std::cout << "Geometry:" << std::endl;
  for (auto [str,type] : cells)
  {
    std::cout << str << ":" << std::endl;

    auto [geo_data,geo_shape] = basix::cell::geometry<double>(type);
    auto geometry = basix::element::mdspan_t<double,2>{geo_data.data(),geo_shape};
    for (std::size_t i = 0; i < geo_shape[0]; ++i) {
      std::cout << "[ ";
      for (std::size_t j = 0; j < geo_shape[1]; ++j)
        std::cout << geometry(i,j) << " ";
      std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
  }


  std::cout << "Topology:" << std::endl;
  for (auto [str,type] : cells)
  {
    std::cout << str << ":" << std::endl;

    auto topology = basix::cell::topology(type);
    for (std::size_t d = 0; d < topology.size(); ++d) {
      std::cout << "  num entities(dim=" << d << ") = " << topology[d].size() << std::endl;
      for (std::size_t s = 0; s < topology[d].size(); ++s) {
        std::cout << "    entity(" << s << ") = [ ";
        for (std::size_t v = 0; v < topology[d][s].size(); ++v) {
          std::cout << topology[d][s][v] << " ";
        }
        std::cout << "]" << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  return success ? 0 : 1;
}
