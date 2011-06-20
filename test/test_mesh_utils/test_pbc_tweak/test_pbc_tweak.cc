/**
 * @file   test_facet_extraction_tetra1.cc
 * @author Guillaume ANCIAUX <guillaume.anciaux@epfl.ch>
 * @date   Thu Aug 19 13:05:27 2010
 *
 * @brief  test of internal facet extraction
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "pbc_synchronizer.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 3;

  akantu::initialize(&argc, &argv);
  akantu::debug::setDebugLevel(akantu::dblInfo);

  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  mesh_io.read("cube.msh", mesh);
  mesh.computeBoundingBox();
  std::map<UInt,UInt> pbc_pair;
  MeshUtils::computePBCMap(mesh,0,pbc_pair);
  // if (flag_y)computePBCMap(mesh,1,pbc_pair);
  // if (flag_z) computePBCMap(mesh,2,pbc_pair);


  {
    std::map<UInt,UInt>::iterator it = pbc_pair.begin();
    std::map<UInt,UInt>::iterator end = pbc_pair.end();
    
    Real * coords = mesh.getNodes().values;
    UInt dim = mesh.getSpatialDimension();
    while(it != end){
      UInt i1 = (*it).first;
      UInt i2 = (*it).second;
      
      AKANTU_DEBUG_INFO("pairing " << i1 << "(" 
			<< coords[dim*i1] << "," << coords[dim*i1+1] << "," 
			<< coords[dim*i1+2]
			<< ") with"
			<< i2 << "(" 
			<< coords[dim*i2] << "," << coords[dim*i2+1] << "," 
			<< coords[dim*i2+2]
			<< ")");	
      ++it;
    }
  }

  PBCSynchronizer synch(pbc_pair);
  
  SolidMechanicsModel * model = new SolidMechanicsModel(mesh);
  model->getSynchronizerRegistry().registerSynchronizer(synch,_gst_smm_uv);
  model->getSynchronizerRegistry().registerSynchronizer(synch,_gst_smm_mass);

  model->initVectors();
  UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  
  model->getForce().clear();
  model->getVelocity().clear();
  model->getAcceleration().clear();
  model->getDisplacement().clear();
  
  model->initExplicit();
  model->initModel();
  model->readMaterials("material.dat");
  model->initMaterials();

  model->changeEquationNumberforPBC(pbc_pair);
  model->assembleMassLumped();
  model->getSynchronizerRegistry().synchronize(_gst_smm_mass);

 #ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);

  dumper.SetPoints(mesh.getNodes().values, dim, nb_nodes, "test-pbc-tweak");
  dumper.SetConnectivity((int*)mesh.getConnectivity(_tetrahedron_4).values,
			 TETRA1, mesh.getNbElement(_tetrahedron_4), C_MODE);
  dumper.AddNodeDataField(model->getMass().values,
			  3, "mass");
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();

#endif //AKANTU_USE_IOHELPER

  return EXIT_SUCCESS;
}
