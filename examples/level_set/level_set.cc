/**
 * @file   level_set.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Fri Feb 24 03:42:39 2012
 *
 * @brief  Level set method example.
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
#include "level_set_model.hh"
#include "sphere.hh"
#include "plane.hh"
#include "container.hh"
#include "aka_memory.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

//#define bar_length 0.01
//Sol   #define bar_height 0.01

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
    debug::setDebugLevel(dblWarning);
    initialize(argc, argv);

    UInt spatial_dimension = 2;
    //const ElementType type = _triangle_3;

    Real * Center = new Real [spatial_dimension];

    Center[0] = 0.25;
    Center[1] = 0.5;
    //Center[2] = 0.5;
    Real R = 0.15;
    sphere Surface(spatial_dimension);
    Surface.setCenterRadius(Center, R);

    /*Center[0] = 0.25;
    Center[1] = 0.75;
    //Center[2] = 0.5;
    R = 0.1;
    sphere Surface2(spatial_dimension);
    Surface2.setCenterRadius(Center, R);

    Center[0] = 0.6;
    Center[1] = 0.25;
    //Center[2] = 0.5;
    R = 0.05;
    sphere Surface3(spatial_dimension);
    Surface3.setCenterRadius(Center, R);
    
    Center[0] = 0.7;
    Center[1] = 0.6;
    //Center[2] = 0.5;
    R = 0.05;
    sphere Surface4(spatial_dimension);
    Surface4.setCenterRadius(Center, R);

    container Surface;
    Surface.addGeometry(&Surface1);
    Surface.addGeometry(&Surface2);
    Surface.addGeometry(&Surface3);
    Surface.addGeometry(&Surface4);*/


    /*Real * Normal = new Real [spatial_dimension];
    Real * Point = new Real [spatial_dimension];
    Normal[0] = 0.0;
    Normal[1] = 1.0;
    Point[0] = 0.5;
    Point[1] = 0.5;
    plane Surface(spatial_dimension);
    Surface.setNormalPoint(Normal, Point);*/


    Mesh mesh(spatial_dimension);
    mesh.read("square.msh");
    

    LevelSetModel model(mesh);

    /// model initialization
    model.initFull();
    Real Filter = 0.05; //e=5h
    model.initPhi(Surface, true, Filter);

    model.setBaseName("level_set_SUPG");
    model.addDumpField("phi");
    model.addDumpField("boundary");
    model.addDumpField("v");

    Array<Real> & nodes = const_cast<Array<Real> &> (mesh.getNodes());
    UInt nb_nodes = nodes.getSize();


    UInt n_steps = 300;
    Array<Real> * velocity = new Array<Real > (nb_nodes, spatial_dimension, "velocity");
    velocity->clear();
    Real * v_nodes = velocity->values;
    Real time_step = 1.0 / n_steps;

    //create velocity
    //Translation
    /*for (UInt i = 0; i < nb_nodes; i++) {
        v_nodes[i * spatial_dimension] = 1.0;
        v_nodes[i * spatial_dimension + 1] = 0.0;
    }*/

    //Rotation
    Real w = 2 * 3.1415927; /// n_steps;
    for (UInt i = 0; i < nb_nodes; i++) {
        v_nodes[i * spatial_dimension] = -(nodes(i, 1) - 0.5) * w;
        v_nodes[i * spatial_dimension + 1] = (nodes(i, 0) - 0.5) * w;
    }

    //Vortex
    /*Real pi = 3.1415927; /// n_steps;
    for (UInt i = 0; i < nb_nodes; i++) {
        v_nodes[i * spatial_dimension] = -2.0 * sin(pi * nodes(i, 0)) * sin(pi * nodes(i, 0)) * sin(pi * nodes(i, 1)) * cos(pi * nodes(i, 1));
        v_nodes[i * spatial_dimension + 1] = 2.0 * sin(pi * nodes(i, 0)) * cos(pi * nodes(i, 0)) * sin(pi * nodes(i, 1)) * sin(pi * nodes(i, 1));
    }*/

    //Expansion
    /*for (UInt i = 0; i < nb_nodes; i++) {
        for (UInt j = 0; j < spatial_dimension; j++)
            v_nodes[i * spatial_dimension + j] = 10 * (nodes(i, j) - Center[j]);
    }*/
    //velocity(0, 0)=1;
    //velocity(0, 1)=1;

    /*Array<Real> & phi = model.getPhi();
    phi(0, 0) = -0.5;
    phi(1, 0) = 0.5;
    phi(2, 0) = 0.5;
    //phi(3,0)=-0.5;
     */

    model.dump();

    
    UInt freq = 5;
    Real tol = 5e-18;
    UInt Nmax = 15;
    for (UInt i = 0; i <= n_steps; i++) {
        //Oscillation velocity
        /*Real w = 2 * 3.1415927; /// n_steps;
        for (UInt m = 0; m < nb_nodes; m++) {
            Real norm = 0.0;
            for (UInt j = 0; j < spatial_dimension; j++) {
                v_nodes[m * spatial_dimension + j] = (nodes(m, j) - Center[j]);
                norm += v_nodes[m * spatial_dimension + j] * v_nodes[m * spatial_dimension + j];
            }
            norm = sqrt(norm);
            for (UInt j = 0; j < spatial_dimension; j++) {
                if(norm>Math::getTolerance())
                    v_nodes[m * spatial_dimension + j] *= 1.0/norm;
                else
                    v_nodes[m * spatial_dimension + j] =0.0;
            }
        }
        
        for (UInt m = 0; m < nb_nodes; m++) {
            Real theta=atan(fabs((nodes(m,0)-Center[0])/(nodes(m,1)-Center[1])));
            Real s=0.25*cos(8*theta)*sin(w*time);
            v_nodes[m * spatial_dimension] *= s;
            v_nodes[m * spatial_dimension + 1] *= s;
        }
        time+=time_step;*/
        model.transportLevelSet(velocity, time_step, true);
        model.dump();
        Real norm;
        model.testConvergenceResidual(1e-4, norm);
        std::cout << "Step " << i + 1 << " of " << n_steps << ", residual norm : " << norm << std::endl;
        if ((i + 1) % freq == 0)
            model.reinitializeLevelSet(time_step, tol, Nmax, 2 * Filter / (5 * 3.1415927));
    }

    delete [] Center;
    delete velocity;
    //delete Normal;
    //delete Point;

    finalize();

    return EXIT_SUCCESS;
}
