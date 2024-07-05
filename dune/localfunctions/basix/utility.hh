// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BASIX_UTILITY_HH
#define DUNE_LOCALFUNCTIONS_BASIX_UTILITY_HH

#if HAVE_BASIX

namespace Dune::Impl {

  // Map the cell::type from basix to the Dune GeometryType
  GeometryType geometryType (basix::cell::type cell_type)
  {
    switch (cell_type) {
      case basix::cell::type::point:          return GeometryTypes::vertex;
      case basix::cell::type::interval:       return GeometryTypes::line;
      case basix::cell::type::triangle:       return GeometryTypes::triangle;
      case basix::cell::type::tetrahedron:    return GeometryTypes::tetrahedron;
      case basix::cell::type::quadrilateral:  return GeometryTypes::quadrilateral;
      case basix::cell::type::hexahedron:     return GeometryTypes::hexahedron;
      case basix::cell::type::prism:          return GeometryTypes::prism;
      case basix::cell::type::pyramid:        return GeometryTypes::pyramid;
      default: return GeometryTypes::none(0);
    }
  }

  // Map the entity index of basix into the Dune numbering
  int entityIndex (basix::cell::type cell_type, int dim, int s)
  {
    // {cell_type, dimension, entity}
    static constexpr int perm[8][4][12]{
      { {0} }, // point
      { {0,1}, {0} }, // interval
      { {0,1,2}, {2,1,0}, {0} }, // triangle
      { {0,1,2,3}, {5,4,2,3,1,0}, {3,2,1,0}, {0} }, // tetrahedron
      { {0,1,2,3}, {2,0,2,3}, {0} }, // quadrilateral
      { {0,1,2,3,4,5,6,7}, {6,4,0,5,1,7,2,3,10,8,9,11}, {4,2,0,1,3,5}, {0} }, // hexahedron
      { {0,1,2,3,4,5}, {3,4,0,5,1,2,6,7,8}, {3,0,1,2,4}, {0} }, // prism
      { {0,1,2,3,4}, {2,0,4,1,5,3,6,7}, {0,3,1,2,4}, {0} } // pyramid
    };

    return perm[(int)(cell_type)][dim][s];
  }

  // Map the derivative-order tuple from the partial() method to the
  // derivative index used in basix
  template <std::size_t dim>
  int indexing (const std::array<unsigned int,dim>& orders)
  {
    if constexpr (dim == 1)
      return basix::indexing::idx(orders[0]);
    else if constexpr (dim == 2)
      return basix::indexing::idx(orders[0],orders[1]);
    else if constexpr (dim == 3)
      return basix::indexing::idx(orders[0],orders[1],orders[2]);
    else
      return 0;
  }

  // Map the first order derivative in direction d to the derivative index
  // used in basix
  template <std::size_t dim>
  int indexing (unsigned int d)
  {
    std::array<unsigned int,dim> orders;
    orders[d] = 1;
    return indexing(orders);
  }

} // end namespace Dune::Impl

#endif // HAVE_BASIX
#endif // DUNE_LOCALFUNCTIONS_BASIX_UTILITY_HH