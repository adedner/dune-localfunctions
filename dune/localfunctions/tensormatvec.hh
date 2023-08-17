#ifndef DUNE_LOCALFUNCTIONS_TENSORMATVEC_HH
#define DUNE_LOCALFUNCTIONS_TENSORMATVEC_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune
{
  /**
   * @brief Implements the "downgrade" from a 1xN matrix to a vector. Noop if applied to a vector.
   *
   * @tparam V Datatype
   * @tparam N Dimension of Vector
   * @return FieldVector<V, N> const&
   */
  // const version
  template <class V, int N>
  FieldVector<V, N> const &matrixToVector(FieldVector<V, N> const &v)
  {
    return v;
  }

  template <class V, int N>
  FieldVector<V, N> const &matrixToVector(FieldMatrix<V, 1, N> const &m)
  {
    return m[0];
  }
  // non const version
  template <class V, int N>
  FieldVector<V, N> &matrixToVector(FieldVector<V, N> &v)
  {
    return v;
  }

  template <class V, int N>
  FieldVector<V, N> &matrixToVector(FieldMatrix<V, 1, N> &m)
  {
    return m[0];
  }

  /**
   * @brief Implements the "downgrade" from a 1xNxN Tensor to a NxN  Matrix. Noop if applied to a
   * matrix.
   *
   * @tparam V Datatype
   * @tparam N Dimension of Vector
   * @return FieldVector<V, N> const&
   */

  // const version
  template <class V, int N>
  FieldMatrix<V, N, N> const &tensorToMatrix(FieldVector<FieldMatrix<V, N, N>, 1> const &t)
  {
    return t[0];
  }

  template <class V, int N>
  FieldMatrix<V, N, N> const &tensorToMatrix(FieldMatrix<V, N, N> const &m)
  {
    return m;
  }

  // non const version
  template <class V, int N>
  FieldMatrix<V, N, N> &tensorToMatrix(FieldVector<FieldMatrix<V, N, N>, 1> &t)
  {
    return t[0];
  }

  template <class V, int N>
  FieldMatrix<V, N, N> &tensorToMatrix(FieldMatrix<V, N, N> &m)
  {
    return m;
  }
} // namespace Dune
#endif