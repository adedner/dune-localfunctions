// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_COMMON_DERIVATIVE_HH
#define DUNE_LOCALFUNCTIONS_COMMON_DERIVATIVE_HH

#include <limits>
#include <type_traits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
/** \brief This file introduces a makeDerivative(f) function which either uses derivative(f) or implements a finite difference scheme.
 */
namespace Dune
{
  namespace Impl
  {

    template <typename Signature, typename Enable = void>
    struct DerivativeTraits;

    template <typename S, typename T, int n>
    struct DerivativeTraits<S(FieldVector<T, n>),
                            std::void_t<std::common_type_t<S, T>>>
    {
      using ST = std::common_type_t<S, T>;
      using type = FieldVector<ST, n>;
    };

    template <typename S, int m, typename T, int n>
    struct DerivativeTraits<FieldVector<S, m>(FieldVector<T, n>),
                            std::void_t<std::common_type_t<S, T>>>
    {
      using ST = std::common_type_t<S, T>;
      using type = FieldMatrix<ST, m, n>;
    };

    template <typename S, typename T, int n>
    struct DerivativeTraits<FieldVector<S, 1>(FieldVector<T, n>),
                            std::void_t<std::common_type_t<S, T>>>
    {
      using ST = std::common_type_t<S, T>;
      // using type = FieldVector<ST, n>;
      using type = FieldMatrix<ST, 1, n>;
    };

    template <typename S, typename T, int n, int m>
    struct DerivativeTraits<FieldMatrix<S, 1, n>(FieldVector<T, m>),
                            std::void_t<std::common_type_t<S, T>>>
    {
      using ST = std::common_type_t<S, T>;
      using type = FieldMatrix<ST, n, m>;
    };

    // call `f.derivative(Domain)`
    template <typename Domain, typename F,
              decltype((std::declval<F>().derivative(std::declval<Domain>()), true)) = true>
    auto makeDerivativeImpl(const F &f, PriorityTag<5>)
    {
      return [f](const Domain &x)
      { return f.derivative(x); };
    }

    // call `f.jacobian(Domain, F::JacobianRangeType)`
    template <typename Domain, typename F,
              typename J = typename F::JacobianRangeType,
              decltype((std::declval<F>().jacobian(std::declval<Domain>(), std::declval<J &>()), true)) = true>
    auto makeDerivativeImpl(const F &f, PriorityTag<4>)
    {
      return [f](const Domain &x)
      {
        J jacobian;
        f.jacobian(x, jacobian);
        return jacobian;
      };
    }

    // call `f.jacobian(Domain, F::JacobianType)`
    template <typename Domain, typename F,
              typename J = typename F::JacobianType,
              decltype((std::declval<F>().jacobian(std::declval<Domain>(), std::declval<J &>()), true)) = true>
    auto makeDerivativeImpl(const F &f, PriorityTag<3>)
    {
      return [f](const Domain &x)
      {
        J jacobian;
        f.jacobian(x, jacobian);
        return jacobian;
      };
    }

    // call `derivative(f)(Domain)`
    template <typename Domain, typename F,
              decltype((derivative(std::declval<F>())(std::declval<Domain>()), true)) = true>
    auto makeDerivativeImpl(const F &f, PriorityTag<2>)
    {
      return derivative(f);
    }

    namespace
    {
      template <typename S, typename T, int n,
                std::enable_if_t<std::is_convertible_v<S, T>, int> = 0>
      void assignPartialDerivative(const S &df_i, FieldVector<T, n> &J, std::size_t i)
      {
        J[i] = df_i;
      }

      template <typename S, int m, typename T, int n,
                std::enable_if_t<std::is_convertible_v<S, T>, int> = 0>
      void assignPartialDerivative(const FieldVector<S, m> &df_i, FieldMatrix<T, m, n> &J, std::size_t i)
      {
        for (int j = 0; j < m; ++j)
          J[j][i] = df_i[j];
      }

      template <typename S, int m, typename T, int n,
                std::enable_if_t<std::is_convertible_v<S, T>, int> = 0>
      void assignPartialDerivative(const FieldMatrix<S, 1, m> &df_i, FieldMatrix<T, m, n> &J,
                                   std::size_t i)
      {
        for (int j = 0; j < m; ++j)
          J[j][i] = df_i[0][j];
      }

      template <typename S, typename T, int n,
                std::enable_if_t<std::is_convertible_v<S, T>, int> = 0>
      void assignPartialDerivative(const S &df_i, FieldMatrix<T, 1, n> &J, std::size_t i)
      {
        J[0][i] = df_i;
      }

      template <typename S, typename T, int n,
                std::enable_if_t<std::is_convertible_v<S, T>, int> = 0>
      void assignPartialDerivative(const FieldVector<S, 1> &df_i, FieldMatrix<T, 1, n> &J, std::size_t i)
      {
        J[0][i] = df_i[0];
      }
    }

    // function is not differentiable. Try simple finite-difference derivative
    template <typename Domain, typename F>
    auto makeDerivativeImpl(const F &ff, PriorityTag<0>)
    {
      auto &&f = makeFunctionWithCallOperator<Domain>(ff);
      using Range = std::decay_t<decltype(f(std::declval<Domain>()))>;
      using Jacobian = typename DerivativeTraits<Range(Domain)>::type;
      return [f](const Domain &x) -> Jacobian
      {
        using std::sqrt;
        using T = typename FieldTraits<Domain>::field_type;
        const T h = sqrt(std::numeric_limits<T>::epsilon());
        Jacobian df;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
          Domain x1{x}, x2{x};
          x1[i] -= h;
          x2[i] += h;
          assignPartialDerivative((f(x2) - f(x1)) / (2 * h), df, i);
        }
        return df;
      };
    }

    /// Construct a callable representing the derivative of a function f
    template <typename Domain, typename F>
    auto makeDerivative(const F &f)
    {
      return makeDerivativeImpl<Domain>(f, PriorityTag<42>{});
    }

  }
} // end namespace Dune::Impl

#endif // DUNE_LOCALFUNCTIONS_COMMON_DERIVATIVE_HH
