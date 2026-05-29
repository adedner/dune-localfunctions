// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_UTILITY_TENSORVIEW_HH
#define DUNE_LOCALFUNCTIONS_UTILITY_TENSORVIEW_HH


#include <array>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>
#include <dune/common/std/extents.hh>
#include <dune/common/std/no_unique_address.hh>

namespace Dune {

/**
 * \brief A tensor wrapper of arbitrary rank and individual dimension extents, static
 * or dynamic.
 *
 * \ingroup Tensors
 * \nosubgrouping
 *
 * \tparam Element  The element type stored in the tensor
 * \tparam extents  Individual static extents or std::dynamic_extent
 **/
template <class Element, std::size_t... exts>
class TensorView
{
public:
  using extents_type = Std::extents<int,exts...>;
  using size_type = typename extents_type::size_type;
  using rank_type = typename extents_type::rank_type;
  using index_type = typename extents_type::rank_type;

  using element_type = Element;
  using value_type = std::remove_const_t<element_type>;
  using reference = std::add_lvalue_reference_t<element_type>;
  using pointer = std::add_pointer_t<std::remove_reference_t<reference>>;

private:
  template <class Seq>
  struct ViewTraits;

  template <std::size_t... II>
  struct ViewTraits<std::index_sequence<II...>>
  {
    template <std::size_t> using IndexType = index_type;

    using type = std::function<reference(IndexType<II>...)>;
  };

public:
  using view_type = typename ViewTraits<std::make_index_sequence<sizeof...(exts)>>::type;

  /// \name TensorView constructors
  /// @{

  /// \brief Constructor stores the tensor in a span and stores the view functor by value
  template <class View, std::convertible_to<std::size_t>... Extents>
  constexpr TensorView (View view, Extents... extents) noexcept
    : view_(std::move(view))
    , extents_(extents...)
  {}

  /// @}

  /// \name Multi index access
  /// @{

  /// \brief Access specified element at position i0,i1,...
  template <std::convertible_to<index_type>... Indices>
  constexpr reference operator() (const Indices&... indices) const
  {
    return view_(indices...);
  }

  /// \brief Access element at position [{i0,i1,...}]
  template <std::convertible_to<index_type> Index>
  constexpr reference operator[] (const std::array<Index,extents_type::rank()>& indices) const
  {
    return std::apply([&](auto... i) -> reference { return view_(i...); }, indices);
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
  view_type view_;
  DUNE_NO_UNIQUE_ADDRESS extents_type extents_;
};

} // end namespace Dune

#endif