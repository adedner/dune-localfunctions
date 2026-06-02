// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_UTILITY_TENSORWRAPPER_HH
#define DUNE_LOCALFUNCTIONS_UTILITY_TENSORWRAPPER_HH


#include <array>
#include <functional>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>
#include <dune/common/indices.hh>
#include <dune/common/referencehelper.hh>
#include <dune/common/concepts/number.hh>
#include <dune/common/std/extents.hh>
#include <dune/common/std/no_unique_address.hh>

namespace Dune {
namespace Impl {

template <class Container>
struct ElementType
{
  using field_type = typename FieldTraits<Container>::field_type;
  using type = std::conditional_t<std::is_const_v<Container>,
    std::add_const_t<field_type>, field_type>;
};


template <class Container>
struct StaticOrDynamicExtent
{
  static constexpr std::size_t value = std::dynamic_extent;
};

template <class Container>
  requires requires{ { Container::size() } -> std::convertible_to<std::size_t>; }
struct StaticOrDynamicExtent<Container>
{
  static constexpr std::size_t value = Container::size();
};

template <std::size_t E0, class ETail>
struct PushFront;

template <std::size_t E0, class I, std::size_t... EE>
struct PushFront<E0,Std::extents<I,EE...>>
{
  using type = Std::extents<I,E0,EE...>;
};

template <class Extents, class I, std::size_t... EE>
Extents push_front (std::size_t e0, Std::extents<I,EE...> const& tail)
{
  return unpackIntegerSequence([&](auto... ii) {
    return Extents{e0, tail.extent(ii)...};
  }, std::make_index_sequence<sizeof...(EE)>{});
}


template <class Container>
struct ExtentsType;

template <class Scalar>
  requires Concept::Number<Scalar>
struct ExtentsType<Scalar>
{
  using type = Std::extents<int>;
  static constexpr type make (Scalar const& c)
  {
    return type{};
  }
};

template <class Container>
  requires requires(Container const& c) { c[0]; }
struct ExtentsType<Container>
{
  using block_type = std::decay_t<decltype(std::declval<Container const&>()[0])>;
  using type = typename PushFront<StaticOrDynamicExtent<Container>::value,
    typename ExtentsType<block_type>::type>::type;

  static constexpr type make (Container const& c)
  {
    return push_front<type>(c.size(), ExtentsType<block_type>::make(c[0]));
  }
};

} // end namespace Impl


/**
 * \brief A tensor wrapper of a nested container transforms the container into a tensor like interface.
 *
 * \ingroup Tensors
 * \nosubgrouping
 *
 * \tparam Container  A vector of vector data structure
 **/
template <class Container>
class TensorWrapper
{
public:
  using extents_type = typename Impl::ExtentsType<Container>::type;
  using size_type = typename extents_type::size_type;
  using rank_type = typename extents_type::rank_type;
  using index_type = typename extents_type::rank_type;

  using element_type = typename Impl::ElementType<Container>::type;
  using value_type = std::remove_const_t<element_type>;
  using reference = std::add_lvalue_reference_t<element_type>;
  using pointer = std::add_pointer_t<std::remove_reference_t<reference>>;

public:

  /// \name TensorWrapper constructors
  /// @{

  /// \brief Constructor stores the tensor in a span and stores the view functor by value
  constexpr explicit TensorWrapper (Container& container) noexcept
    : container_(&container)
    , extents_(Impl::ExtentsType<Container>::make(container))
  {}

  /// @}


  /// \name Multi index access
  /// @{

  /// \brief Access specified element at position i0,i1,...
  template <std::convertible_to<index_type> Index,
            std::convertible_to<index_type>... Indices>
  constexpr reference operator() (Index index, Indices... indices) const
  {
    if constexpr (extents_type::rank() == 1)
      return (*container_)[index];
    else if constexpr (extents_type::rank() == 2)
      return (*container_)[index][first(indices...)];
    else if constexpr (extents_type::rank() == 3)
      return (*container_)[index][first(indices...)][second(indices...)];
    else
      return Dune::TensorWrapper{(*container_)[index]}(indices...);
  }

  /// \brief Access element at position [{i0,i1,...}]
  template <std::convertible_to<index_type> Index>
  constexpr reference operator[] (const std::array<Index,extents_type::rank()>& indices) const
  {
    if constexpr (extents_type::rank() == 1)
      return (*container_)[indices[0]];
    else if constexpr (extents_type::rank() == 2)
      return (*container_)[indices[0]][indices[1]];
    else if constexpr (extents_type::rank() == 3)
      return (*container_)[indices[0]][indices[1]][indices[2]];
    else
      return std::apply([&](auto i0, auto... ii) -> reference {
        return Dune::TensorWrapper{(*container_)[i0]}(ii...); }, indices);
  }

  /// \brief Create a slice by fixing the first dimension at `index`
  template <std::convertible_to<index_type> Index>
  constexpr decltype(auto) operator[] (Index index) const
  {
    if constexpr (extents_type::rank() == 1)
      return (*container_)[index];
    else
      return Dune::TensorWrapper{(*container_)[index]};
  }

  /// @}


  /// \name Size information
  /// @{

  /// \brief Number of dimensions of the tensor
  static constexpr rank_type rank () noexcept { return extents_type::rank(); }

  /// \brief Number of elements in the r'th dimension of the tensor
  static constexpr std::size_t static_extent (rank_type r) noexcept
  {
    return extents_type::static_extent(r);
  }

  /// \brief Number of elements in the r'th dimension of the tensor
  constexpr index_type extent (rank_type r) const noexcept { return extents_.extent(r); }

  /// \brief Number of elements in all dimensions of the tensor, \related Impl::Extents
  constexpr const extents_type& extents () const noexcept { return extents_; }

  /// @}

private:
  static constexpr auto first (auto i0, auto... ii) { return i0; }
  static constexpr auto second (auto i0, auto i1, auto... ii) { return i1; }

private:
  Container* container_;
  DUNE_NO_UNIQUE_ADDRESS extents_type extents_;
};

} // end namespace Dune

#endif