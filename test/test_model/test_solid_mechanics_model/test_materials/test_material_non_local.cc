/**
 * @file   test_non_local_material.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Mon Aug 22 14:07:08 2011
 *
 * @brief  test of the main part of the non local materials
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
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.hh"


static void paraviewInit(iohelper::Dumper & dumper, const SolidMechanicsModel & model);
//static void paraviewDump(iohelper::Dumper & dumper);
#endif

ByElementTypeReal quadrature_points_volumes("quadrature_points_volumes", "test");
const ElementType TYPE = _triangle_6;

int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);
  debug::setDebugLevel(akantu::dblWarning);

  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  MeshIOMSH mesh_io;
  mesh_io.read("mesh.msh", mesh);

  SolidMechanicsModel model(mesh);
  model.initModel();
  model.initVectors();
  model.readMaterials("material_non_local.dat");
  model.initMaterials();

  //  model.getFEM().getMesh().initByElementTypeVector(quadrature_points_volumes, 1, 0);
  const MaterialNonLocal<BaseWeightFunction> & mat =
    dynamic_cast<const MaterialNonLocal<BaseWeightFunction> &>(model.getMaterial(0));
  //  mat.computeQuadraturePointsNeighborhoudVolumes(quadrature_points_volumes);
  Real radius = mat.getRadius();

  UInt nb_element  = mesh.getNbElement(TYPE);
  UInt nb_tot_quad = model.getFEM().getNbQuadraturePoints(TYPE) * nb_element;

  Vector<Real> quads(0, spatial_dimension);
  quads.resize(nb_tot_quad);

  model.getFEM().interpolateOnQuadraturePoints(mesh.getNodes(),
					       quads, spatial_dimension,
					       TYPE);

  Vector<Real>::iterator<types::RVector> first_quad_1 = quads.begin(spatial_dimension);
  Vector<Real>::iterator<types::RVector> last_quad_1 = quads.end(spatial_dimension);

  std::ofstream pout;
  pout.open("bf_pairs");
  UInt q1 = 0;

  Real R = mat.getRadius();

  for(;first_quad_1 != last_quad_1; ++first_quad_1, ++q1) {
    Vector<Real>::iterator<types::RVector> first_quad_2 = quads.begin(spatial_dimension);
    //Vector<Real>::iterator<types::RVector> last_quad_2 = quads.end(spatial_dimension);
    UInt q2 = 0;
    for(;first_quad_2 != last_quad_1; ++first_quad_2, ++q2) {
      Real d = first_quad_2->distance(*first_quad_1);
      if(d <= radius) {
	Real alpha = (1 - d*d/(R*R));
	alpha = alpha*alpha;
	pout << q1 << " " << q2 << " " << alpha << std::endl;
      }
    }
  }
  pout.close();

  mat.savePairs("cl_pairs");

  ByElementTypeReal constant("constant_value", "test");
  mesh.initByElementTypeVector(constant, 1, 0);
  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator last_type = mesh.lastType(spatial_dimension);
  for(; it != last_type; ++it) {
    UInt nb_quadrature_points = model.getFEM().getNbQuadraturePoints(*it);
    UInt _nb_element = mesh.getNbElement(*it);

    Vector<Real> & constant_vect = constant(*it);
    constant_vect.resize(_nb_element * nb_quadrature_points);

    std::fill_n(constant_vect.storage(), nb_quadrature_points * _nb_element, 1.);
  }

  ByElementTypeReal constant_avg("constant_value_avg", "test");
  mesh.initByElementTypeVector(constant_avg, 1, 0);

  mat.weightedAvergageOnNeighbours(constant, constant_avg, 1);

  debug::setDebugLevel(akantu::dblTest);
  std::cout << constant(TYPE) << std::endl;
  std::cout << constant_avg(TYPE) << std::endl;
  debug::setDebugLevel(akantu::dblWarning);

#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  paraviewInit(dumper, model);
#endif

  akantu::finalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
/* iohelper::Dumper vars                                                      */
/* -------------------------------------------------------------------------- */

#ifdef AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */
template <ElementType type> static iohelper::ElemType paraviewType();
template <> static iohelper::ElemType paraviewType<_segment_2>()      { return iohelper::LINE1;     }
template <> static iohelper::ElemType paraviewType<_segment_3>()      { return iohelper::LINE2;     }
template <> static iohelper::ElemType paraviewType<_triangle_3>()     { return iohelper::TRIANGLE1; }
template <> static iohelper::ElemType paraviewType<_triangle_6>()     { return iohelper::TRIANGLE2; }
template <> static iohelper::ElemType paraviewType<_quadrangle_4>()   { return iohelper::QUAD1;     }
template <> static iohelper::ElemType paraviewType<_tetrahedron_4>()  { return iohelper::TETRA1;    }
template <> static iohelper::ElemType paraviewType<_tetrahedron_10>() { return iohelper::TETRA2;    }
template <> static iohelper::ElemType paraviewType<_hexahedron_8>()   { return iohelper::HEX1;      }

/* -------------------------------------------------------------------------- */
void paraviewInit(iohelper::Dumper & dumper, const SolidMechanicsModel & model) {
  UInt spatial_dimension = ElementClass<TYPE>::getSpatialDimension();
  UInt nb_nodes   = model.getFEM().getMesh().getNbNodes();
  UInt nb_element = model.getFEM().getMesh().getNbElement(TYPE);

  std::stringstream filename; filename << "material_non_local_" << TYPE;

  dumper.SetMode(iohelper::TEXT);
  dumper.SetParallelContext(StaticCommunicator::getStaticCommunicator()->whoAmI(),
			    StaticCommunicator::getStaticCommunicator()->getNbProc());
  dumper.SetPoints(model.getFEM().getMesh().getNodes().values,
		   spatial_dimension, nb_nodes, filename.str().c_str());
  dumper.SetConnectivity((int *)model.getFEM().getMesh().getConnectivity(TYPE).values,
			 paraviewType<TYPE>(), nb_element, iohelper::C_MODE);
  dumper.AddElemDataField(quadrature_points_volumes(TYPE).storage(),
   			  1, "volume");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
}

/* -------------------------------------------------------------------------- */
// void paraviewDump(iohelper::Dumper & dumper) {
//   dumper.Dump();
// }

/* -------------------------------------------------------------------------- */
#endif
