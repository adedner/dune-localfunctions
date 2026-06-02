// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_UTILITY_TENSORVIEW_HH
#define DUNE_LOCALFUNCTIONS_UTILITY_TENSORVIEW_HH


#include <array>
#include <dune/common/indices.hh>
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
template <class View, std::size_t ext0, std::size_t... exts>
class TensorView
{
public:
  using view_type = View;
  using extents_type = Std::extents<int,ext0,exts...>;
  using size_type = typename extents_type::size_type;
  using rank_type = typename extents_type::rank_type;
  using index_type = typename extents_type::rank_type;

  using reference = decltype(std::declval<View>()(0u,exts...));
  using element_type = std::remove_reference_t<reference>;
  using value_type = std::remove_const_t<element_type>;

public:

  /// \name TensorView constructors
  /// @{

  /// \brief Constructor stores the tensor in a span and stores the view functor by value
  template <std::convertible_to<std::size_t>... Extents>
  constexpr explicit TensorView (View view, Extents... extents) noexcept
    : view_(std::move(view))
    , extents_(extents...)
  {}

  template <class Extents>
    requires (std::is_constructible_v<extents_type, std::decay_t<Extents>>)
  constexpr TensorView (View view, Extents&& extents) noexcept
    : view_(std::move(view))
    , extents_(std::forward<Extents>(extents))
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

  /// \brief Create a slice by fixing the first dimension at `index`
  template <std::convertible_to<index_type> Index>
  constexpr decltype(auto) operator[] (const Index& index) const
  {
    if constexpr (extents_type::rank() == 1)
      return view_(index);
    else
      return unpackIntegerSequence([&](auto j0, auto... jj) {
        return Dune::TensorView{
          [i0=index, view=view_](auto... ii) -> reference { return view(i0,ii...); },
          Std::extents<index_type,exts...>(extents_.extent(jj)...)};
      }, std::make_index_sequence<1+sizeof...(exts)>{});
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

template <class View, class I, std::size_t... exts>
TensorView (View, Std::extents<I,exts...>) -> TensorView<View,exts...>;

} // end namespace Dune

#endif