/**
 * @file   plasticity_implicit.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Fri Feb 24 03:42:39 2012
 *
 * @brief  This code refers to the implicit static example from the user manual
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
#include "mesh_io.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

#define bar_length 1.
#define bar_height 1.

class SurfaceLoad : public SolidMechanicsModel::SurfaceLoadFunctor {
public:
    // if a traction vector is provided:

    SurfaceLoad(UInt surface_id, Real load) : surface_id(surface_id), load(load) {
        Theta = 0.0;
    }

    void setload(Real new_load) {
        load = new_load;
    }

    void setTheta(Real d_theta) {
        Theta += d_theta;
    }

    void traction(const Vector<Real> & position,
            Vector<Real> & force,
            const Vector<Real> & normal,
            Surface surface_id) {
        if (surface_id == this->surface_id) {
            force(2) = cos(Theta) * load;
            force(1) = -sin(Theta) * load;
        }
    }

private:
    UInt surface_id;
    Real load;
    Real Theta;
};
//

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
    debug::setDebugLevel(dblWarning);
    initialize(argc, argv);

    UInt spatial_dimension = 3;
    ElementType surf_type = _triangle_3;
    //const ElementType type = _segment_2;

    Mesh mesh(spatial_dimension);
    mesh.read("cube.msh");

    std::cout << mesh << std::endl;

    SolidMechanicsModel model(mesh);

    /// model initialization
    model.initFull("material.dat", _static);

    model.setBaseName("static");
    model.addDumpField("displacement");
    model.addDumpField("force");
    model.addDumpField("residual");
    model.addDumpField("increment");
    model.addDumpField("stress");
    model.addDumpField("strain");

    // boundary conditions
    Array<Real> & position = mesh.getNodes();
    Array<bool> & boundary = model.getBoundary();
    Array<Real> & displacement = model.getDisplacement();
    Array<Real> & displacement_t = model.getDisplacement_t();
    Array<Real> & increment = model.getIncrement();
    Array<Real> & force = model.getForce();


    ///////////////////////////////////////////////////////
    // PRESSURE
    /////////////////////////////////////////////////////////
    mesh.setSurfaceIDsFromIntData("tag_0");
    Array<UInt> & surf = mesh.getSurfaceID(surf_type, _not_ghost);

    for (UInt i = 0; i < surf.getSize(); ++i) {
        surf(i)--;
    }
    //


    UInt nb_nodes = mesh.getNbNodes();

    Array<bool> * boundary_normal = new Array<bool > (nb_nodes, 3, "boundary_normal");
    boundary_normal->clear();
    Array<Real> * EulerAngles = new Array<Real > (nb_nodes, 3, "Euler_Angles");
    EulerAngles->clear();

    boundary.clear();
    displacement.clear();
    displacement_t.clear();
    increment.clear();
    Real Force_z = 2e9;

    for (UInt n = 0; n < nb_nodes; ++n) {
        if (std::abs(position(n, 2)) < Math::getTolerance()) 
            boundary(n, 2) = true;

        if (std::abs(position(n, 1)) < Math::getTolerance())
            boundary(n, 1) = true;

        if (std::abs(position(n, 0)) < Math::getTolerance()) 
            boundary(n, 0) = true;
    }

    model.updateStresses();
    model.UpdateStressesAtT();

    model.assembleStiffnessMatrix();
    model.getStiffnessMatrix().saveMatrix("stiffness.out");

    ///////////////////////////////////////////////////////
    // PRESSURE
    /////////////////////////////////////////////////////////
    UInt s_P0 = 29; //surface id where a load of 1/11 P is applied
    SurfaceLoad ss1(s_P0 - 1, Force_z);

    FEM & fem_boundary = model.getFEMBoundary();
    fem_boundary.initShapeFunctions();
    fem_boundary.computeNormalsOnControlPoints();

    force.clear();
    model.computeForcesFromFunction(ss1, _bft_traction);

    /////////////////////////////////////////////////////////

    model.updateResidual();
    model.dump();

    Real norm, norm2, norm0;
    Real tol = 1e-2;
    UInt count = 0;
    UInt tractions = 150;
    UInt Nmax = 100;
    for (UInt i = 0; i < tractions; i++) {
        increment.clear();
        displacement_t.clear();

        if (i) {

            ss1.setload((i + 1) * Force_z);

            force.clear();
            model.computeForcesFromFunction(ss1, _bft_traction);

            model.updateStresses();

            model.assembleStiffnessMatrix();
            model.getStiffnessMatrix().saveMatrix("stiffness.out");
            model.updateResidual();


        }
        std::cout << "Traction : " << i + 1 << std::endl;
        count = 0;
        model.testConvergenceResidual(tol, norm);
        model.testConvergenceIncrement(tol / 100, norm2);
        bool keep_working = true;
        norm0 = norm;
        while (keep_working && count < Nmax) {
            if (!count) {
                std::cout << " -- initial : " << " - residual norms : " << norm << " " << norm2 << std::endl;
                count++;
            } else
                std::cout << " -- Iter : " << count++ << " - residual norms : " << norm << " " << norm2 << std::endl;

            if (norm0 > 1e-9) {
                model.solveStatic(*boundary_normal, *EulerAngles);


                for (UInt n = 0; n < nb_nodes; ++n)
                    for (UInt m = 0; m < spatial_dimension; ++m) {
                        increment(n, m) = displacement(n, m) - displacement_t(n, m);
                        displacement_t(n, m) = displacement(n, m);
                    }


                model.updateStresses();
                model.assembleStiffnessMatrix();
                model.updateResidual();
                model.getStiffnessMatrix().saveMatrix("stiffness.out");

                keep_working = !model.testConvergenceResidual(tol, norm);
                keep_working = !model.testConvergenceIncrement(1e-6, norm2);

                norm /= norm0;

                if (norm < 1e-6)
                    keep_working = false;
            } else
                keep_working = false;

        }

        std::cout << " -- Iter : " << count << " - residual norms : " << norm << " " << norm2 << std::endl;

        model.UpdateStressesAtT();

        //if (!(i % 3))
        model.dump();
    }


    delete boundary_normal;
    delete EulerAngles;

    model.dump();

    finalize();

    return EXIT_SUCCESS;
}
