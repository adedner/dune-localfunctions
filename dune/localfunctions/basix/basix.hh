// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BASIX_BASIX_HH
#define DUNE_LOCALFUNCTIONS_BASIX_BASIX_HH

#include <utility>
#include <vector>

#include <basix/finite-element.h>
#include <basix/indexing.h>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/std/mdarray.hh>
#include <dune/common/std/mdspan.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>

namespace Dune {
namespace Impl {

  GeometryType basix2dune (basix::cell::type cell_type)
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

  template <std::size_t dim>
  int indexing (const std::array<unsigned int,dim>& orders)
  {
    if constexpr (dim == 1)
      return orders[0];
    else if constexpr (dim == 2)
      return basix::indexing::idx(orders[0],orders[1]);
    else if constexpr (dim == 3)
      return basix::indexing::idx(orders[0],orders[1],orders[2]);
    else
      return 0;
  }

  template <std::size_t dim>
  int indexing (unsigned int d)
  {
    std::array<unsigned int,dim> orders;
    orders[d] = 1;
    return indexing(orders);
  }


} // end namespace Impl

template <class F, int dimDomain, int dimRange>
class BasixLocalFiniteElement
{
  using FE = basix::FiniteElement<F>;

public:
    struct LocalBasis
    {
      using Traits = LocalBasisTraits<F,dimDomain,FieldVector<F,dimDomain>,F,dimRange,FieldVector<F,dimRange>,FieldMatrix<F,dimRange,dimDomain> >;

      std::size_t size () const
      {
        return fe_->dim();
      }

      std::size_t order () const
      {
        return fe_->embedded_superdegree();
      }

      //! \brief Evaluate all shape functions
      void evaluateFunction(const typename Traits::DomainType& x,
                            std::vector<typename Traits::RangeType>& out) const
      {
        auto [tab_data,shape] = fe_->tabulate(0, basix::element::mdspan_t<const F, 2>{x.data(), 1, dimDomain});
        basix::element::mdspan_t<const F, 4> tab(tab_data.data(), shape);

        assert(shape[0] == 1);
        assert(shape[1] == 1);
        assert(shape[2] == size());
        assert(shape[3] == dimRange);

        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          for (std::size_t j = 0; j < dimRange; ++j)
            out[i][j] = tab(0,0,i,j);
      }


      //! \brief Evaluate all shape function jacobians
      void evaluateJacobian(const typename Traits::DomainType& x,
                            std::vector<typename Traits::JacobianType>& out) const
      {
        auto [tab_data,shape] = fe_->tabulate(1, basix::element::mdspan_t<const F, 2>{x.data(), 1, dimDomain});
        basix::element::mdspan_t<const F, 4> tab(tab_data.data(), shape);

        assert(shape[1] == 1);
        assert(shape[2] == size());
        assert(shape[3] == dimRange);

        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          for (std::size_t j = 0; j < dimRange; ++j)
            for (std::size_t k = 0; k < dimDomain; ++k)
              out[i][j][k] = tab(Impl::indexing<dimDomain>(k),0,i,j);
      }

      //! \brief Evaluate all shape function jacobians
      void partial(const std::array<unsigned int,dimDomain>& order,
                   const typename Traits::DomainType& x,
                   std::vector<typename Traits::RangeType>& out) const
      {
        int totalOrder = std::accumulate(order.begin(), order.end(), 0);

        auto [tab_data,shape] = fe_->tabulate(totalOrder, basix::element::mdspan_t<const F, 2>{x.data(), 1, dimDomain});
        basix::element::mdspan_t<const F, 4> tab(tab_data.data(), shape);

        assert(shape[1] == 1);
        assert(shape[2] == size());
        assert(shape[3] == dimRange);

        out.resize(size());
        for (std::size_t i = 0; i < size(); ++i)
          for (std::size_t j = 0; j < dimRange; ++j)
            out[i][j] = tab(Impl::indexing(order),0,i,j);
      }

      FE* fe_;
    };

    struct LocalCoefficients
    {
      LocalCoefficients (FE* fe)
        : fe_(fe)
        , localKeys_(fe_->dim())
      {
        auto& entity_dofs = fe_->entity_dofs();
        for (std::size_t d = 0; d < entity_dofs.size(); ++d)
          for (std::size_t s = 0; s < entity_dofs[d].size(); ++s)
            for (std::size_t i = 0; i < entity_dofs[d][s].size(); ++i)
              localKeys_[entity_dofs[d][s][i]] = LocalKey(s,Impl::basix2dune(fe_->cell_type()).dim()-d, i);
      }

      std::size_t size () const
      {
        return localKeys_.size();
      }

      //! get i'th index
      const LocalKey& localKey (std::size_t i) const
      {
        return localKeys_[i];
      }

    private:
      FE* fe_;
      std::vector<LocalKey> localKeys_;
    };


    struct LocalInterpolation
    {
      template<class Func, class C>
      void interpolate (const Func& f, std::vector<C>& out) const
      {
        auto& p = fe_->points();
        using Points = Std::mdspan<const F, Std::extents<int, Std::dynamic_extent, dimDomain>>;
        Points points{p.first.data(), p.second[0]};

        Std::mdarray<F,Std::dextents<int,1>> values{points.extent(0)};
        for (std::size_t i = 0; i < points.extent(0); ++i)
        {
          FieldVector<F,dimDomain> x;
          for (int j = 0; j < dimDomain; ++j)
            x[j] = points(i,j);

          values(i) = f(x);
        }

        auto& m = fe_->interpolation_matrix();
        using InterpolationMatrix = Std::mdspan<const F, Std::dextents<int, 2>>;
        InterpolationMatrix interpolationMatrix{m.first.data(), m.second[0], m.second[1]};

        out.resize(interpolationMatrix.extent(0));
        for (int i = 0; i < interpolationMatrix.extent(0); ++i) {
          out[i] = F(0);
          for (int j = 0; j < interpolationMatrix.extent(1); ++j)
            out[i] += interpolationMatrix(i,j) * values(j);
        }
      }

      FE* fe_;
    };


    using Traits = LocalFiniteElementTraits<
      LocalBasis,
      LocalCoefficients,
      LocalInterpolation>;

public:
  explicit BasixLocalFiniteElement (const FE& fe)
    : fe_(fe)
    , localBasis_{&fe_}
    , localCoefficients_{&fe_}
    , localInterpolation_{&fe_}
  {}

  explicit BasixLocalFiniteElement (FE&& fe)
    : fe_(std::move(fe))
    , localBasis_{&fe_}
    , localCoefficients_{&fe_}
    , localInterpolation_{&fe_}
  {}

  BasixLocalFiniteElement (const BasixLocalFiniteElement& other)
    : BasixLocalFiniteElement(other.fe_)
  {}

  BasixLocalFiniteElement (BasixLocalFiniteElement&& other)
    : BasixLocalFiniteElement(std::move(other.fe_))
  {}


  const LocalBasis& localBasis () const { return localBasis_; }
  const LocalCoefficients& localCoefficients () const { return localCoefficients_; }
  const LocalInterpolation& localInterpolation () const { return  localInterpolation_; }


  std::size_t size () const
  {
    return fe_.dim();
  }

  GeometryType type () const
  {
    return Impl::basix2dune(fe_.cell_type());
  }

private:
  FE fe_;

  LocalBasis localBasis_;
  LocalCoefficients localCoefficients_;
  LocalInterpolation localInterpolation_;
};

} // end namespace Dune






#endif // DUNE_LOCALFUNCTIONS_BASIX_BASIX_HH
