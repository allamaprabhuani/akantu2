/**
 * @file   scalability_test.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Feb 22 09:35:58 2011
 *
 * @brief  Test de scalability
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "static_communicator.hh"
#include "mesh_partition_scotch.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[])
{
  akantu::debug::setDebugLevel(akantu::dblWarning);
  initialize(&argc, &argv);

  /* -------------------------------------------------------------------------- */

  UInt spatial_dimension = 2;
  ElementType type = _quadrangle_4;

  UInt nex = 50, ney = 5;
  Real width = 1., height = 1.;
  if(argc == 3) {
    nex = atoi(argv[1]);
    ney = atoi(argv[2]);
  } else if (argc != 1) {
    std::cout << "Usage : " << argv[0] << " [nb_element_x nb_element_y]" << std::endl;
    exit(EXIT_FAILURE);
  }

  /* ------------------------------------------------------------------------ */
  Mesh mesh(spatial_dimension);

  std::cout << "Generating mesh..." << std::endl;
  Real height_el = height / ney;
  Real width_el = width / nex;
  UInt nnx = nex + 1, nny = ney + 1;

  Vector<Real> & nodes = const_cast<Vector<Real> &>(mesh.getNodes());
  nodes.resize(nnx * nny);

  mesh.addConnecticityType(type);
  Vector<UInt> & connectivity = const_cast<Vector<UInt> &>(mesh.getConnectivity(type));
  UInt t = 0;

  if(type == _quadrangle_4) t = nex * ney;
  if(type == _triangle_3) t = nex * ney * 2;

  connectivity.resize(t);

  for (UInt i = 0; i < nnx; ++i) {
    for (UInt j = 0; j < nny; ++j) {
      UInt n = i * nny + j;
      nodes.at(n, 0) = i * width_el;
      nodes.at(n, 1) = j * height_el;
    }
  }

  for (UInt i = 0; i < nex; ++i) {
    for (UInt j = 0; j < ney; ++j) {
      if(type == _quadrangle_4) {
	/* 
	 *    3        2
	 *     x------x
	 *     |      |
	 *     |      |
	 *     |      |
	 *     x------x
	 *    0        1
	 */
	UInt e = (i * ney + j);
	connectivity.at(e, 0) = i * nny + j;
	connectivity.at(e, 1) = (i + 1) * nny + j;
	connectivity.at(e, 2) = (i + 1) * nny + (j + 1);
	connectivity.at(e, 3) = i * nny + j + 1;
      }

      if(type == _triangle_3) {
	/*
	 *    3        2
	 *     x------x
	 *     |    / |
	 *     |  /   |
	 *     |/     |
	 *     x------x
	 *    0        1
	 */
	UInt e = (i * ney + j) * 2;
	connectivity.at(e, 0) = i * nny + j;
	connectivity.at(e, 1) = (i + 1) * nny + j;
	connectivity.at(e, 2) = (i + 1) * nny + (j + 1);

	connectivity.at(e + 1, 0) = i * nny + j;
	connectivity.at(e + 1, 1) = (i + 1) * nny + (j + 1);
	connectivity.at(e + 1, 2) = i * nny + j + 1;
      }
    }
  }

  akantu::MeshIOMSH mesh_io;
  mesh_io.write("bar.msh", mesh);

  return EXIT_SUCCESS;
}
