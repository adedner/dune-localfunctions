// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PK3DLOCALCOEFFICIENTS_HH
#define DUNE_PK3DLOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
     \brief Please doc me!

     \nosubgrouping
   */
  template<unsigned int k>
  class Pk3DLocalCoefficients
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    : public LocalCoefficientsInterface
#else
    : public LocalCoefficientsInterface<Pk3DLocalCoefficients<k> >
#endif
  {
    enum {N = (k+1)*(k+2)*(k+3)/6};

  public:
    //! \brief Standard constructor
    Pk3DLocalCoefficients () : li(N)
    {
      unsigned int vertexmap[4] = {0, 1, 2, 3};
      generate_local_keys(vertexmap);
    }

    //! constructor for eight variants with order on edges flipped
    Pk3DLocalCoefficients (unsigned int vertexmap[4]) : li(N)
    {
      generate_local_keys(vertexmap);
    }

    //! number of coefficients
    int size () const
    {
      return N;
    }

    //! get i'th index
    const LocalKey& localKey (int i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;

    void generate_local_keys(unsigned int vertexmap[4])
    {
      unsigned int subindex[16];
      unsigned int codim_count[4] = {0};
      for (unsigned int m = 1; m < 16; ++m)
      {
        unsigned int codim = !(m&1) + !(m&2) + !(m&4) + !(m&8);
        subindex[m] = codim_count[codim]++;
      }

      unsigned int dof_count[16] = {0};
      unsigned int i[4];
      unsigned int n = 0;
      for (i[3] = 0; i[3] <= k; ++i[3])
        for (i[2] = 0; i[2] <= k - i[3]; ++i[2])
          for (i[1] = 0; i[1] <= k - i[2] - i[3]; ++i[1])
          {
            i[0] = k - i[1] - i[2] - i[3];
            unsigned int entity = 0;
            unsigned int codim = 0;
            for (unsigned int m = 0; m < 4; ++m)
            {
              unsigned int j = !vertexmap[i[m]];
              entity += !j << m;
              codim += j;
            }
            li[n++] = LocalKey(subindex[entity], codim, dof_count[entity]++);
          }
    }
  };

}

#endif