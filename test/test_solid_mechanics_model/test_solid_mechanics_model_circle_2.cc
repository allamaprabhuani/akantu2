/**
 * @file   test_solid_mechanics_model.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 14:34:13 2010
 *
 * @brief  test of the class SolidMechanicsModel
 *
 * @section LICENSE
 *
 * \<insert license here\>
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

static void trac(__attribute__ ((unused)) double * position,double * traction){
  memset(traction,0,sizeof(Real)*4);
  traction[0] = 1000;
  traction[3] = 1000;

  // if(fabs(position[0])< 1e-4){
  //   traction[0] = -position[1];
  // }

}

int main(int argc, char *argv[])
{
  UInt max_steps = 10000;
  Real epot, ekin;

  ElementType type = _line_2;

  Mesh mesh(2);
  MeshIOMSH mesh_io;
  mesh_io.read("circle2.msh", mesh);

  SolidMechanicsModel * model = new SolidMechanicsModel(mesh);

  /// model initialization
  model->initVectors();
  UInt nb_nodes = model->getFEM().getMesh().getNbNodes();
  memset(model->getForce().values,        0, 2*nb_nodes*sizeof(Real));
  memset(model->getVelocity().values,     0, 2*nb_nodes*sizeof(Real));
  memset(model->getAcceleration().values, 0, 2*nb_nodes*sizeof(Real));
  memset(model->getDisplacement().values, 0, 2*nb_nodes*sizeof(Real));
  memset(model->getResidual().values,     0, 2*nb_nodes*sizeof(Real));
  memset(model->getMass().values,     1, nb_nodes*sizeof(Real));

  model->readMaterials("material.dat");
  model->initMaterials();
  model->initModel();

  Real time_step = model->getStableTimeStep();
  model->setTimeStep(time_step/10.);

  model->assembleMassLumped();

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
  fem_boundary.computeNormalsOnQuadPoints();

  const Mesh::ConnectivityTypeList & type_list = fem_boundary.getMesh().getConnectivityTypeList();
  Mesh::ConnectivityTypeList::const_iterator it;
  for(it = type_list.begin(); it != type_list.end(); ++it) {
    if(Mesh::getSpatialDimension(*it) != fem_boundary.getElementDimension()) continue;

    //    ElementType facet_type = Mesh::getFacetElementType(*it);
    UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
    UInt nb_quad              = FEM::getNbQuadraturePoints(*it);

    UInt nb_element;
    const Vector<Real> * shapes;
    Vector<Real> quad_coords(0,2,"quad_coords");
    const Vector<Real> * normals_on_quad;

    nb_element   = fem_boundary.getMesh().getNbElement(*it);
    fem_boundary.interpolateOnQuadraturePoints(mesh.getNodes(), quad_coords, 2, type);
    normals_on_quad = &(fem_boundary.getNormalsOnQuadPoints(*it));

    shapes       = &(fem_boundary.getShapes(*it));

    Vector<Real> * sigma_funct = new Vector<Real>(nb_element*nb_quad, 4, "myfunction");
    Vector<Real> * funct = new Vector<Real>(nb_element*nb_quad, 2, "myfunction");

    Real * sigma_funct_val = sigma_funct->values;
    Real * shapes_val = shapes->values;

    /// compute t * \phi_i for each nodes of each element
    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt q = 0; q < nb_quad; ++q) {
	trac(quad_coords.values+el*nb_quad*2+q,sigma_funct_val);
	sigma_funct_val += 4;
      }
    }

    Math::matrix_vector(2,2,*sigma_funct,*normals_on_quad,*funct);
    funct->extendComponentsInterlaced(nb_nodes_per_element,2);

    Real * funct_val = funct->values;
    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt q = 0; q < nb_quad; ++q) {
	for (UInt n = 0; n < nb_nodes_per_element; ++n) {
	  *funct_val++ *= *shapes_val;
	  *funct_val++ *= *shapes_val++;
	}
      }
    }


    Vector<Real> * int_funct = new Vector<Real>(nb_element, 2*nb_nodes_per_element,
						    "inte_funct");
    fem_boundary.integrate(*funct, *int_funct, 2*nb_nodes_per_element, *it);
    delete funct;

    // fem_boundary.assembleVector(*int_funct,const_cast<Vector<Real> &>(model->getForce()), 2, *it);
    delete int_funct;
  }


  //  model->getDisplacement().values[1] = 0.1;


#ifdef AKANTU_USE_IOHELPER
  DumperParaview dumper;
  dumper.SetMode(BASE64);

  dumper.SetPoints(model->getFEM().getMesh().getNodes().values, 2, nb_nodes, "coordinates");
  dumper.SetConnectivity((int *)model->getFEM().getMesh().getConnectivity(_triangle_2).values,
			 TRIANGLE2, model->getFEM().getMesh().getNbElement(_triangle_2), C_MODE);
  dumper.AddNodeDataField(model->getDisplacement().values, 2, "displacements");
  dumper.AddNodeDataField(model->getVelocity().values, 2, "velocity");
  dumper.AddNodeDataField(model->getForce().values, 2, "force");
  dumper.AddNodeDataField(model->getMass().values, 1, "Mass");
  dumper.AddNodeDataField(model->getResidual().values, 2, "residual");
  dumper.AddElemDataField(model->getMaterial(0).getStrain(_triangle_2).values, 4, "strain");
  dumper.AddElemDataField(model->getMaterial(0).getStress(_triangle_2).values, 4, "stress");
  dumper.SetEmbeddedValue("displacements", 1);
  dumper.SetEmbeddedValue("force", 1);
  dumper.SetEmbeddedValue("residual", 1);
  dumper.SetEmbeddedValue("velocity", 1);
  dumper.SetPrefix("paraview/");
  dumper.Init();
  dumper.Dump();
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
    if(s % 100 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  return EXIT_SUCCESS;
}



