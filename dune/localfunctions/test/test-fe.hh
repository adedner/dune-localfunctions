// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

// This header is not part of the official Dune API and might be subject
// to change.  You can use this header to test external finite element
// implementations, but be warned that your tests might break with future
// Dune versions.

#ifndef DUNE_LOCALFUNCTIONS_TEST_TEST_FE_HH
#define DUNE_LOCALFUNCTIONS_TEST_TEST_FE_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <vector>

#include <dune/common/classname.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/simd/simd.hh>
#include <dune/common/simd/loop.hh>

#include <dune/geometry/quadraturerules.hh>

// This class defines a local finite element function.
// It is determined by a local finite element and
// representing the local basis and a coefficient vector.
// This provides the evaluate method needed by the interpolate()
// method.
template<class FE, class R>
class FEFunction
{
  const FE& fe;

public:
  using RangeType = R;
  using BasisRangeType = typename FE::Traits::Basis::Traits::Range;
  using DomainType = typename FE::Traits::Basis::Traits::DomainLocal;

  std::vector<R> coeff;

  FEFunction(const FE& fe_) : fe(fe_) { resetCoefficients(); }

  void resetCoefficients() {
    coeff.resize(fe.basis().size());
    for(std::size_t i=0; i<coeff.size(); ++i)
      coeff[i] = 0;
  }

  void setRandom(double max) {
    coeff.resize(fe.basis().size());
    for(std::size_t i=0; i<coeff.size(); ++i)
      for(std::size_t l = 0; l < Dune::Simd::lanes(coeff[i]); ++l)
        Dune::Simd::lane(l,coeff[i][0]) = ((1.0*std::rand()) / RAND_MAX - 0.5)*2.0*max;
  }

  RangeType operator() (const DomainType& x) const {
    RangeType y;
    std::vector<BasisRangeType> yy;
    fe.basis().evaluateFunction(x, yy);
    y = 0.0;
    for (std::size_t i=0; i<yy.size(); ++i)
      y.axpy(coeff[i][0], yy[i]);
    return y;
  }

};


// Check if interpolation is consistent with basis evaluation.
/**
 * This test generates a local coefficient vector with random values from a
 * certain range (-100..100).  It then uses the basis to wrap this coefficient
 * vector into an element-local discrete function.  This is then interpolated
 * into another coefficient vector using the interpolation of the finite
 * element.  The two coefficient vectors are then compared.
 *
 * \param FE  The finite element to check
 * \param eps Tolerance when comparing floating-point values
 * \param n   Number of times to run the check.
 */
template<class FE, class RangeFieldType = typename FE::Traits::Basis::Traits::Range::field_type>
bool testInterpolation(const FE& fe, double eps, int n=5)
{
  using std::abs;
  using std::max;

  bool success = true;
  using BasisRangeType = typename FE::Traits::Basis::Traits::Range;
  using CoeffType = Dune::FieldVector<RangeFieldType, BasisRangeType::dimension>;
  FEFunction<FE, CoeffType> f(fe);
  std::vector<CoeffType> coeff;
  for(int i=0; i<n && success; ++i) {
    // Set random coefficient vector
    f.setRandom(100);

    // Compute interpolation weights

    //////////////////////////////////////////////////////////////////////////////
    // Feed the shape functions to the 'interpolate' method in form of a callable.
    //////////////////////////////////////////////////////////////////////////////
    fe.interpolation().interpolate(f, coeff);

    // Check size of weight vector
    if (coeff.size() != fe.basis().size()) {
      std::cout << "Bug in LocalInterpolation for finite element type "
                << Dune::className<FE>() << ":" << std::endl;
      std::cout << "    Interpolation vector has size " << coeff.size()
                << std::endl;
      std::cout << "    Basis has size " << fe.basis().size() << std::endl;
      std::cout << std::endl;
      success = false;

      // skip rest of loop since that depends on matching sizes
      continue;
    }

    // Check if interpolation weights are equal to coefficients
    for(std::size_t j=0; j<coeff.size() && success; ++j) {
      if ( Dune::Simd::anyTrue(abs(coeff[j][0]-f.coeff[j][0]) >
                               (2*coeff.size()*eps)*(max(abs(f.coeff[j][0]), RangeFieldType(1.0))) ))
      {
        std::cout << std::setprecision(16);
        std::cout << "Bug in LocalInterpolation for finite element type "
                  << Dune::className<FE>() << ":" << std::endl;
        std::cout << "    Interpolation weight " << j << " differs by "
                  << abs(coeff[j][0]-f.coeff[j][0]) << " from coefficient of "
                  << "linear combination." << std::endl;
        std::cout << std::endl;
        success = false;
      }
    }
  }
  return success;
}

// check whether Jacobian agrees with FD approximation
/**
 * \param geo   The geometry the finite element is tested on.
 * \param fe    The finite element to test.
 * \param eps   Tolerance for comparing floating-point values.  When comparing
 *              numerical derivatives, this is divided by \c delta to yield an
 *              even bigger tolerance.
 * \param delta Stepsize to use when doing numerical derivatives.
 * \param order The Jacobian is checked at a number of quadrature points.
 *              This parameter determines the order of the quatrature rule
 *              used to obtain the quadrature points.
 */
template<class Geo, class FE>
bool testJacobian(const Geo &geo, const FE& fe, double eps, double delta,
                  std::size_t order = 2)
{
  typedef typename FE::Traits::Basis Basis;

  typedef typename Basis::Traits::DomainField DF;
  static const std::size_t dimDLocal = Basis::Traits::dimDomainLocal;
  typedef typename Basis::Traits::DomainLocal DomainLocal;
  static const std::size_t dimDGlobal = Basis::Traits::dimDomainGlobal;

  static const std::size_t dimR = Basis::Traits::dimRange;
  typedef typename Basis::Traits::Range Range;

  typedef typename Basis::Traits::Jacobian Jacobian;

  bool success = true;

  // ////////////////////////////////////////////////////////////
  //   Check the partial derivatives by comparing them
  //   to finite difference approximations
  // ////////////////////////////////////////////////////////////

  // A set of test points
  const Dune::QuadratureRule<DF, dimDLocal> quad =
    Dune::QuadratureRules<DF, dimDLocal>::rule(fe.type(),order);

  // Loop over all quadrature points
  for (std::size_t i=0; i < quad.size(); i++) {

    // Get a test point
    const DomainLocal& testPoint = quad[i].position();

    // Get the shape function derivatives there
    std::vector<Jacobian> jacobians;
    fe.basis().evaluateJacobian(testPoint, jacobians);
    if(jacobians.size() != fe.basis().size()) {
      std::cout << "Bug in evaluateJacobianGlobal() for finite element type "
                << Dune::className<FE>() << ":" << std::endl;
      std::cout << "    Jacobian vector has size " << jacobians.size()
                << std::endl;
      std::cout << "    Basis has size " << fe.basis().size() << std::endl;
      std::cout << std::endl;
      return false;
    }

    Dune::FieldMatrix<DF, dimDLocal, dimDGlobal> geoJT =
      geo.jacobianTransposed(testPoint);

    // Loop over all shape functions in this set
    for (std::size_t j=0; j<fe.basis().size(); ++j) {

      // basis.evaluateJacobian returns global derivatives, however we can
      // only do local derivatives, so transform the derivatives back into
      // local coordinates
      Dune::FieldMatrix<double, dimR, dimDLocal> localJacobian(0);
      for(std::size_t k = 0; k < dimR; ++k)
        for(std::size_t l = 0; l < dimDGlobal; ++l)
          for(std::size_t m = 0; m < dimDLocal; ++m)
            localJacobian[k][m] += jacobians[j][k][l] * geoJT[m][l];

      // Loop over all local directions
      for (std::size_t m = 0; m < dimDLocal; ++m) {

        // Compute an approximation to the derivative by finite differences
        DomainLocal upPos   = testPoint;
        DomainLocal downPos = testPoint;

        upPos[m]   += delta;
        downPos[m] -= delta;

        std::vector<Range> upValues, downValues;

        fe.basis().evaluateFunction(upPos,   upValues);
        fe.basis().evaluateFunction(downPos, downValues);

        //Loop over all components
        for(std::size_t k = 0; k < dimR; ++k) {

          // The current partial derivative, just for ease of notation
          double derivative = localJacobian[k][m];

          double finiteDiff = (upValues[j][k] - downValues[j][k]) / (2*delta);

          // Check
          if ( std::abs(derivative-finiteDiff) >
               eps/delta*(std::max(std::abs(finiteDiff), 1.0)) )
          {
            std::cout << std::setprecision(16);
            std::cout << "Bug in evaluateJacobian() for finite element type "
                      << Dune::className<FE>() << ":" << std::endl;
            std::cout << "    Shape function derivative does not agree with "
                      << "FD approximation" << std::endl;
            std::cout << "    Shape function " << j << " component " << k
                      << " at position " << testPoint << ": derivative in "
                      << "local direction " << m << " is "
                      << derivative << ", but " << finiteDiff << " is "
                      << "expected." << std::endl;
            std::cout << std::endl;
            success = false;
          }
        } //Loop over all components
      } // Loop over all local directions
    } // Loop over all shape functions in this set
  } // Loop over all quadrature points

  return success;
}

// call tests for given finite element
template<class Geo, class FE, bool CheckSIMD = false>
bool testFE(const Geo &geo, const FE& fe, double eps, double delta,
            unsigned order = 2)
{
  bool success = true;

  success = testInterpolation<FE>(fe, eps) and success;

  if constexpr (CheckSIMD) {
    typedef typename FE::Traits::Basis::Traits::RangeField CT;
    success = testInterpolation<FE, Dune::LoopSIMD<CT, 4>>(fe, eps) and success;
  }

  std::cout << "Interpolation test? " << (success ? "Success!" : "Failed!") << std::endl;

  success = testJacobian(geo, fe, eps, delta, order) and success;

  return success;
}

#endif // DUNE_LOCALFUNCTIONS_TEST_TEST_FE_HH
