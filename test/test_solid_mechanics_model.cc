/**
 * @file   test_solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 14:34:13 2010
 *
 * @brief  test of the class SolidMechanicsModel
 *
 * @section LICENSE
 *
 * <insert license here>
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"
#include "fem.hh"
/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper.h"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

int main(int argc, char *argv[])
{
  UInt max_steps = 1;
  Real epot, ekin;

  Mesh mesh(2);
  MeshIOMSH mesh_io;
  mesh_io.read("triangle.msh", mesh);

  SolidMechanicsModel * model = new SolidMechanicsModel(mesh);

  /// model initialization
  model->initVectors();
  UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  memset(model->getForce().values,        0, 2*nb_nodes*sizeof(Real));
  memset(model->getVelocity().values,     0, 2*nb_nodes*sizeof(Real));
  memset(model->getAcceleration().values, 0, 2*nb_nodes*sizeof(Real));
  memset(model->getDisplacement().values, 0, 2*nb_nodes*sizeof(Real));
  
  model->readMaterials("material.dat");
  model->initMaterials();
  model->initModel();

  Real time_step = model->getStableTimeStep();
  model->setTimeStep(time_step/10.);

  model->assembleMass();

  std::cout << *model << std::endl;

  /// boundary conditions
  // Real eps = 1e-16;
  // for (UInt i = 0; i < nb_nodes; ++i) {
  //   model->getDisplacement().values[2*i] = model->getFEM().getMesh().getNodes().values[2*i] / 100.;

  //   if(model->getFEM().getMesh().getNodes().values[2*i] <= eps) {
  //     model->getBoundary().values[2*i    ] = true;
  //     if(model->getFEM().getMesh().getNodes().values[2*i + 1] <= eps)
  // 	model->getBoundary().values[2*i + 1] = true;
  //   }
  //   if(model->getFEM().getMesh().getNodes().values[2*i + 1] <= eps) {
  //     model->getBoundary().values[2*i + 1] = true;
  //   }

  // }
  
  FEM & fem_boundary = model->getFEMBoundary();
  fem_boundary.initShapeFunctions();
  const Mesh::ConnectivityTypeList & type_list = fem_boundary.getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != fem_boundary.getElementDimension()) continue;

    //    ElementType facet_type = Mesh::getFacetElementType(*it);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
    UInt nb_quadrature_points = FEM::getNbQuadraturePoints(*it);
    UInt shape_size           = FEM::getShapeSize(*it);
    
    UInt nb_element;
    const Vector<Real> * shapes;

    nb_element   = fem_boundary.getMesh().getNbElement(*it);
    shapes       = &(fem_boundary.getShapes(*it));

    Vector<Real> * funct = new Vector<Real>(nb_element, 2*shape_size, "myfunction");

    Real * funct_val = funct->values;
    Real * shapes_val = shapes->values;

    /// compute rho * \phi_i for each nodes of each element
    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt n = 0; n < shape_size; ++n) {
	*funct_val++ = 1* *shapes_val++;
	*funct_val++ = 0;
      }
    }

    Vector<Real> * int_funct = new Vector<Real>(nb_element, 2*shape_size / nb_quadrature_points,
						    "inte_funct");
    fem_boundary.integrate(*funct, *int_funct, 2*nb_nodes_per_element, *it);
    delete funct;

    fem_boundary.assembleVector(*int_funct,model->getForce(), 2, *it);
    delete int_funct;
  }


  //  model->getDisplacement().values[1] = 0.1;


#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(TEXT);

  dumper.SetPoints(model->getFEM().getMesh().getNodes().values, 2, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(_triangle_1).values,
			 TRIANGLE1, model->getFEM().getMesh().getNbElement(_triangle_1), C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values, 2, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values, 2, "velocity");
  dumper.AddNodeDataField(model->getForce().values, 2, "force");
  dumper.AddNodeDataField(model->getResidual().values, 2, "residual");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(_triangle_1).values, 4, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(_triangle_1).values, 4, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
#endif //AKANTU_USE_IOHELPER

  model->setPotentialEnergyFlagOn();
  for(UInt s = 0; s < max_steps; ++s) {
    model->explicitPred();

    model->updateResidual();
    model->updateAcceleration();
    model->explicitCorr();

    epot = model->getPotentialEnergy();
    ekin = model->getKineticEnergy();

    std::cout << s << " " << epot << " " << ekin << " " << epot + ekin
	      << std::endl;

#ifdef AKANTU_USE_IOHELPER
    if(s % 10 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  return EXIT_SUCCESS;
}



