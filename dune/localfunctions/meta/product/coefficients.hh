// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_META_PRODUCT_COEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_META_PRODUCT_COEFFICIENTS_HH

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune {

class PrismaticProductLocalCoefficients
{
private:
  std::vector<LocalKey> keys_;
  std::map<GeometryType, std::size_t> layout_;

public:
  //! Construct a PrismaticProductLocalCoefficients object
  template <class LC1, class LC2, class M, class RefElem1, class RefElem2>
  PrismaticProductLocalCoefficients (const LC1& lc1, const LC2& lc2, const M& mapping,
                                     const RefElem1& refElem1, const RefElem2& refElem2)
    : keys_(mapping.required_span_size())
  {
    std::array<unsigned int, 2> maxIndex2{0u,0u};
    for (std::size_t j = 0; j < lc2.size(); ++j) {
      const auto& key2 = lc2.localKey(j);
      maxIndex2[key2.codim()] = std::max(maxIndex2[key2.codim()], key2.index());
    }

    // construct a product reference element
    using ctype = std::common_type_t<typename RefElem1::ctype, typename RefElem2::ctype>;
    constexpr int dim = RefElem1::dimension + RefElem2::dimension;
    auto refElem = referenceElement<ctype,dim>(GeometryTypes::prismaticProduct(refElem1.type(), refElem2.type()));

    for (std::size_t i = 0; i < lc1.size(); ++i) {
      for (std::size_t j = 0; j < lc2.size(); ++j) {
        const auto& key1 = lc1.localKey(i);
        const auto& key2 = lc2.localKey(j);

        unsigned int shift1 = key1.codim()+1 > refElem1.dimension ? 0 : refElem1.size(key1.codim()+1);
        unsigned int shift2 = key2.codim()+1 > 1 ? 0 : 1;

        unsigned int subEntity = shift1*shift2 + refElem1.size(key1.codim())*key2.subEntity() + key1.subEntity();
        unsigned int codim = key1.codim() + key2.codim();
        unsigned int index = key1.index() * (maxIndex2[key2.codim()]+1) + key2.index();

        keys_[mapping(i, j)] = LocalKey(subEntity, codim, index);
        layout_[refElem.type(subEntity,codim)]++;
      }
    }
  }

  //! number of coefficients
  std::size_t size () const noexcept { return keys_.size(); }

  //! get i'th index
  const LocalKey& localKey (std::size_t i) const noexcept { return keys_[i]; }

  //! A layout for the mcmgmapper
  auto layout () const
  {
    return [l=layout_](GeometryType gt, int /*dim*/) -> std::size_t
    {
      auto it = l.find(gt);
      return it != l.end() ? it->second : 0;
    };
  }
};

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_PRODUCT_COEFFICIENTS_HH
