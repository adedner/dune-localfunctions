# SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
  raviartthomas0prismlocalbasis.hh
  raviartthomas0prismlocalcoefficients.hh
  raviartthomas0prismlocalinterpolation.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas/raviartthomas0prism)
