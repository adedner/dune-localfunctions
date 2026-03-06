# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(cube)
add_subdirectory(simplex)

install(FILES
  brezzidouglasmarinicube.hh
  brezzidouglasmarinisimplex.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/brezzidouglasmarini)
