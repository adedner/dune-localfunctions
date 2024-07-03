
#include <dune/localfunctions/basix/basix.hh>
#include <dune/localfunctions/test/test-localfe.hh>

#include <basix/finite-element.h>
#include <basix/e-lagrange.h>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::BasixLocalFiniteElement<double, 2, 1> basixlfe{
    basix::element::create_lagrange<double>(basix::cell::type::triangle, 3,
    basix::element::lagrange_variant::equispaced, false)};
  TEST_FE(basixlfe);

  return success ? 0 : 1;
}
