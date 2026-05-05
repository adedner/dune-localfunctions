# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(HEADERS
  bdfmcube.hh
  bdfmsimplex.hh
)

add_subdirectory(cube)
add_subdirectory(simplex)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/brezzidouglasfortinmarini)
