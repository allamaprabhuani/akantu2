/**
 * @file   level_set_model.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Fri Dec 14 11:45:43 2012
 *
 * @brief  Implementation of LevelSetModel class
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
#include "aka_common.hh"
#include "level_set_model.hh"
#include "aka_math.hh"
#include "static_communicator.hh"
#include "generalized_trapezoidal.hh"

#include "sparse_matrix.hh"
#include "solver.hh"
#ifdef AKANTU_USE_MUMPS
#include "solver_mumps.hh"
#include "contact/regular_grid_neighbor_structure.hh"
#endif

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
LevelSetModel::LevelSetModel(Mesh & mesh,
        UInt dim,
        const ID & id,
        const MemoryID & memory_id) :
Model(mesh, id, memory_id),
Dumpable<DumperParaview>(id),
increment_flag(false),
phi_gradient("phi_gradient", id),
phi_on_qpoints("phi_on_qpoints", id),
phi_on_qpoints_boundary("phi_on_qpoints_boundary", id),
v_on_qpoints("v_on_qpoints", id),
v_on_qpoints_boundary("v_on_qpoints_boundary", id),
v_r_on_qpoints("v_r_on_qpoints", id),
element_filter("element_filter", id),
element_filter_boundary("element_filter_boundary", id),
shapes_SUPG("shapes_SUPG", id),
supg_flag(0),
filtered_flag(false),
Filter_parameter(0.0),
spatial_dimension(dim),
mesh(mesh),
solver(NULL) {
    AKANTU_DEBUG_IN();

    createSynchronizerRegistry(this);

    if (spatial_dimension == 0) spatial_dimension = mesh.getSpatialDimension();

    //std::stringstream sstr; sstr << id << ":fem";
    registerFEMObject<MyFEMType > ("LevelSetFEM", mesh, spatial_dimension);
    registerFEMObject<MyFEMType > ("LevelSetFEM_Boundary", mesh, spatial_dimension - 1);

    this->phi = NULL;
    this->v = NULL;
    this->residual = NULL;
    this->boundary = NULL;
    this->increment = NULL;

    addDumpMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);

    //Array<UInt> & elem_filter = element_filter(_triangle_3, _not_ghost);
    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void LevelSetModel::transportLevelSet(Array<Real> * velocity, Real delta_t, bool Assembly_A) {
    AKANTU_DEBUG_IN();

    v = velocity;
    if (Assembly_A) {
        if (fabs(delta_t) < Math::getTolerance())
            AKANTU_DEBUG_ERROR("Time step is equal to zero or it has not been defined");

        time_step = delta_t;
        supg_flag = false;
        assemblePhi(_not_ghost);
        //assemblePhi(_ghost);
    }
    //getStiffnessMatrix().saveMatrix("stiffness.out");
    updateRHS();
    solveStatic();
    AKANTU_DEBUG_OUT();


}

void LevelSetModel::reinitializeLevelSet(Real delta_t, Real tol, UInt max_step, Real Epsilon) {
    AKANTU_DEBUG_IN();
    if (fabs(delta_t) < Math::getTolerance())
        AKANTU_DEBUG_ERROR("Time step is equal to zero or it has not been defined");

    AKANTU_DEBUG_ASSERT(Epsilon > 0,
            "The regularization parameter is equal to zero or has not been defined");

    //setBaseName("level_set_reinit");
    //dump();
    time_step = delta_t;
    supg_flag = false;
    Real norm = 1;
    //Real tol = 1e-17;
    UInt counter = 0;
    boundary->clear();
    //SetZeroIsoValue();
    bool reinitialize = false;
    Real ratio = 1.01;
    if (filtered_flag) {
        UInt nb_nodes = mesh.getNbNodes();
        Real Phi_min = *(phi->values);
        Real Phi_max = *(phi->values);
        for (UInt i = 1; i < nb_nodes; i++) {
            if (*(phi->values + i) > Phi_max)
                Phi_max = *(phi->values + i);

            if (*(phi->values + i) < Phi_min)
                Phi_min = *(phi->values + i);
        }
        if (Phi_max > 2 * ratio * Filter_parameter / 3.1415927 || Phi_min<-2 * ratio * Filter_parameter / 3.1415927)
            reinitialize = true;
    } else
        if (norm > tol && counter < max_step)
        reinitialize = true;

    while (reinitialize) {
        //getStiffnessMatrix().saveMatrix("stiffness.out");
        computeVReinit(Epsilon);
        assemblePhi(_not_ghost, true);
        //assemblePhi(_ghost, true);
        updateRHS(true, Epsilon);
        solveStatic();
        testConvergenceResidual(tol, norm);
        if (StaticCommunicator::getStaticCommunicator().whoAmI() == 0)
            std::cout << "Reinitialization step : " << counter << " residual norm : " << norm << std::endl;
        //dump(); 
        counter++;
        if (filtered_flag) {
            UInt nb_nodes = mesh.getNbNodes();
            Real Phi_min = *(phi->values);
            Real Phi_max = *(phi->values);
            for (UInt i = 1; i < nb_nodes; i++) {
                if (*(phi->values + i) > Phi_max)
                    Phi_max = *(phi->values + i);

                if (*(phi->values + i) < Phi_min)
                    Phi_min = *(phi->values + i);
            }
            if (Phi_max > 2 * ratio * Filter_parameter / 3.1415927 || Phi_min<-2 * ratio * Filter_parameter / 3.1415927)
                reinitialize = true;
            else
                reinitialize = false;
        } else {
            if (norm > tol)
                reinitialize = true;
        }
        if(counter >= max_step)
            reinitialize = false;



    }
    //boundary->clear();
    //setBaseName("level_set_SUPG");

    AKANTU_DEBUG_OUT();


}

/* -------------------------------------------------------------------------- */
void LevelSetModel::initModel() {
    getFEM().initShapeFunctions(_not_ghost);
    getFEM().initShapeFunctions(_ghost);
    getFEMBoundary().initShapeFunctions(_not_ghost);
    getFEMBoundary().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
void LevelSetModel::initParallel(MeshPartition * partition,
        DataAccessor * data_accessor) {
    AKANTU_DEBUG_IN();

    if (data_accessor == NULL) data_accessor = this;
    Synchronizer & synch_parallel = createParallelSynch(partition, data_accessor);

    synch_registry->registerSynchronizer(synch_parallel, _gst_htm_phi);
    synch_registry->registerSynchronizer(synch_parallel, _gst_htm_gradient_phi);

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void LevelSetModel::initPBC() {
    AKANTU_DEBUG_IN();

    Model::initPBC();
    PBCSynchronizer * synch = new PBCSynchronizer(pbc_pair);

    synch_registry->registerSynchronizer(*synch, _gst_htm_phi);
    changeLocalEquationNumberforPBC(pbc_pair, 1);

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void LevelSetModel::initArrays() {
    AKANTU_DEBUG_IN();

    UInt nb_nodes = mesh.getNbNodes();

    std::stringstream sstr_phi;
    sstr_phi << id << ":phi";
    //std::stringstream sstr_v;
    //sstr_v << id << ":v";
    std::stringstream sstr_residual;
    sstr_residual << id << ":residual";
    std::stringstream sstr_boun;
    sstr_boun << id << ":boundary";

    phi = &(alloc<Real > (sstr_phi.str(), nb_nodes, 1, REAL_INIT_VALUE));
    //v = &(alloc<Real > (sstr_v.str(), nb_nodes, spatial_dimension, REAL_INIT_VALUE));
    residual = &(alloc<Real > (sstr_residual.str(), nb_nodes, 1, REAL_INIT_VALUE));
    boundary = &(alloc<bool>(sstr_boun.str(), nb_nodes, 1, false));

    Mesh::ConnectivityTypeList::const_iterator it;

    /* -------------------------------------------------------------------------- */
    // byelementtype vectors

    mesh.initByElementTypeArray(phi_on_qpoints,
            1,
            spatial_dimension);
    mesh.initByElementTypeArray(phi_on_qpoints_boundary,
            1,
            spatial_dimension - 1);

    mesh.initByElementTypeArray(phi_gradient,
            spatial_dimension,
            spatial_dimension);

    mesh.initByElementTypeArray(v_on_qpoints,
            spatial_dimension,
            spatial_dimension);

    mesh.initByElementTypeArray(v_on_qpoints_boundary,
            spatial_dimension,
            spatial_dimension - 1);

    mesh.initByElementTypeArray(v_r_on_qpoints,
            spatial_dimension,
            spatial_dimension);
    
    
        /// for each connectivity types allocate the element filer array of the material
    mesh.initByElementTypeArray(element_filter, 1, spatial_dimension, false, _ek_regular);

    mesh.initByElementTypeArray(element_filter_boundary, 1, spatial_dimension - 1, false, _ek_regular);
    

    for (UInt g = _not_ghost; g <= _ghost; ++g) {
        GhostType gt = (GhostType) g;
        Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt, _ek_not_defined);
        Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_not_defined);

        for (; it != end; ++it) {
            UInt nb_element = mesh.getNbElement(*it, gt);
            UInt nb_quad_points = getFEM().getNbQuadraturePoints(*it, gt) * nb_element;
            //UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);

            phi_on_qpoints(*it, gt).resize(nb_quad_points);
            phi_on_qpoints(*it, gt).clear();

            phi_gradient(*it, gt).resize(nb_quad_points);
            phi_gradient(*it, gt).clear();

            v_on_qpoints(*it, gt).resize(nb_quad_points);
            v_on_qpoints(*it, gt).clear();

            v_r_on_qpoints(*it, gt).resize(nb_quad_points);
            v_r_on_qpoints(*it, gt).clear();

            for (UInt el = 0; el < nb_element; ++el)
                element_filter(*it, gt).push_back(el);

        }
    }

    for (UInt g = _not_ghost; g <= _ghost; ++g) {
        GhostType gt = (GhostType) g;
        Mesh::type_iterator it = mesh.firstType(spatial_dimension - 1, gt, _ek_not_defined);
        Mesh::type_iterator end = mesh.lastType(spatial_dimension - 1, gt, _ek_not_defined);

        for (; it != end; ++it) {
            UInt nb_element = mesh.getNbElement(*it, gt);
            UInt nb_quad_points = getFEMBoundary().getNbQuadraturePoints(*it, gt) * nb_element;
            //UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);

            v_on_qpoints_boundary(*it, gt).resize(nb_quad_points);
            v_on_qpoints_boundary(*it, gt).clear();

            phi_on_qpoints_boundary(*it, gt).resize(nb_quad_points);
            phi_on_qpoints_boundary(*it, gt).clear();

            for (UInt el = 0; el < nb_element; ++el)
                element_filter_boundary(*it, gt).push_back(el);

        }
    }

    /* -------------------------------------------------------------------------- */
    dof_synchronizer = new DOFSynchronizer(mesh, 1);
    dof_synchronizer->initLocalDOFEquationNumbers();
    dof_synchronizer->initGlobalDOFEquationNumbers();

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

LevelSetModel::~LevelSetModel() {
    AKANTU_DEBUG_IN();

    //delete [] conductivity;
    if (dof_synchronizer) delete dof_synchronizer;

    delete stiffness_matrix;
    delete jacobian_matrix;
    delete solver;

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void LevelSetModel::updateRHS(bool reinit, Real Epsilon) {
    AKANTU_DEBUG_IN();

    // start synchronization
    synch_registry->asynchronousSynchronize(_gst_htm_phi);
    // finalize communications
    synch_registry->waitEndSynchronize(_gst_htm_phi);

    //clear the array
    /// first @f$ r = q_{ext} @f$
    //  residual->clear();
    //residual->copy(*external_flux);

    /// then @f$ r -= q_{int} @f$
    // update the not ghost ones
    updateRHS(_not_ghost, reinit, Epsilon);
    // update for the received ghosts
    updateRHS(_ghost, reinit, Epsilon);


#ifndef AKANTU_NDEBUG
    getSynchronizerRegistry().synchronize(akantu::_gst_htm_gradient_phi);
#endif

    AKANTU_DEBUG_OUT();
}

void LevelSetModel::solveStatic() {
    AKANTU_DEBUG_IN();

    AKANTU_DEBUG_INFO("Solving Ku = f");
    AKANTU_DEBUG_ASSERT(stiffness_matrix != NULL,
            "You should first initialize the implicit solver and assemble the stiffness matrix");

    UInt nb_nodes = phi->getSize();
    UInt nb_degree_of_freedom = phi->getNbComponent();

    //  if(method != _static)
    jacobian_matrix->copyContent(*stiffness_matrix);

    jacobian_matrix->applyBoundary(*boundary);

    solver->setRHS(*residual);

    if (!increment) setIncrementFlagOn();

    solver->solve(*increment);

    Real * increment_val = increment->values;
    Real * phi_val = phi->values;
    bool * boundary_val = boundary->values;

    for (UInt n = 0; n < nb_nodes * nb_degree_of_freedom; ++n) {
        if (!(*boundary_val)) {
            *phi_val = *increment_val;
        }

        phi_val++;
        boundary_val++;
        increment_val++;
    }

    AKANTU_DEBUG_OUT();
}

void LevelSetModel::setIncrementFlagOn() {
    AKANTU_DEBUG_IN();

    if (!increment) {
        UInt nb_nodes = mesh.getNbNodes();
        std::stringstream sstr_inc;
        sstr_inc << id << ":increment";
        increment = &(alloc<Real > (sstr_inc.str(), nb_nodes, 1, REAL_INIT_VALUE));
    }

    increment_flag = true;

    AKANTU_DEBUG_OUT();
}

void LevelSetModel::initPhi(geometry & geo, bool filter, Real epsilon) {

    AKANTU_DEBUG_IN();
    UInt nb_nodes = mesh.getNbNodes();
    Array<Real> & node_vec = mesh.getNodes();
    Real * node = node_vec.values;
    Real * phi_node = phi->values;
    for (UInt i = 0; i < nb_nodes; i++) {
        if (filter) {
            filtered_flag = 1;
            Filter_parameter = epsilon;
            AKANTU_DEBUG_ASSERT(Filter_parameter > 0,
                    "The regularization parameter is equal to zero or has not been defined");
            Real dist = geo.distance(node);
            if (dist > Filter_parameter)
                *phi_node = 2 * Filter_parameter / 3.1415927;
            else if (dist<-Filter_parameter)
                *phi_node = -2 * Filter_parameter / 3.1415927;
            else
                *phi_node = 2 * Filter_parameter / 3.1415927 * sin(3.1415927 * dist / (2 * Filter_parameter));
        } else
            *phi_node = geo.distance(node);

        //*phi_node *= 2.0;

        node += spatial_dimension;
        phi_node++;
    }
    AKANTU_DEBUG_OUT();
    return;
}

void LevelSetModel::computeVReinit(Real Epsilon) {

    for (UInt g = _not_ghost; g <= _ghost; ++g) {
        GhostType gt = (GhostType) g;
        Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt, _ek_not_defined);
        Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_not_defined);
        for (; it != end; ++it) {

            Array<UInt> & elem_filter = element_filter(*it, gt);
            getFEM().interpolateOnQuadraturePoints(*phi, phi_on_qpoints(*it, gt), 1, *it, gt);

            //UInt nb_element = elem_filter.getSize();
            //UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
            //UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(*it, gt);

            Array<Real> & v_r_qpoints = v_r_on_qpoints(*it, gt);
            Array<Real> & grad_phi = phi_gradient(*it, gt);
            Array<Real> & phi_qpoints = phi_on_qpoints(*it, gt);

            //v_r_qpoints.resize(nb_quadrature_points * nb_element);

            //grad_phi.resize(nb_quadrature_points * nb_element);

            getFEM().gradientOnQuadraturePoints(*phi, grad_phi, 1, *it, gt, &elem_filter);



            Array<Real>::iterator< Vector<Real> > grad_phi_qpoints_it = grad_phi.begin(spatial_dimension);
            Array<Real>::iterator< Vector<Real> > grad_phi_qpoints_end = grad_phi.end(spatial_dimension);

            Array<Real>::iterator< Vector<Real> > phi_qpoints_it = phi_qpoints.begin(1);
            Array<Real>::iterator< Vector<Real> > v_r_qpoints_it = v_r_qpoints.begin(spatial_dimension);

            for (; grad_phi_qpoints_it != grad_phi_qpoints_end; ++grad_phi_qpoints_it, ++v_r_qpoints_it, ++phi_qpoints_it) {

                Vector<Real> & v_r_qpoints_element = *v_r_qpoints_it;
                Vector<Real> & grad_phi_qpoints_element = *grad_phi_qpoints_it;

                Real phi_value = (*phi_qpoints_it)[0];
                Real sign = sign_phi(phi_value, Epsilon); //phi_value / (sqrt(phi_value * phi_value + 0.01));
                /*if(phi_value>0)
                    sign=1;
                else if(phi_value<0)
                    sign=-1;
                else
                    sign=0;*/

                Real norm_grad_phi = grad_phi_qpoints_element.norm();
                //v_r_qpoints_element = grad_phi_qpoints_element;
                memcpy(v_r_qpoints_element.storage(), grad_phi_qpoints_element.storage(), spatial_dimension * sizeof (Real));

                if (norm_grad_phi > Math::getTolerance())
                    v_r_qpoints_element *= sign / norm_grad_phi;
                else
                    v_r_qpoints_element *= 0.0;
            }


        }
    }

}

/* -------------------------------------------------------------------------- */
void LevelSetModel::updateRHS(const GhostType & ghost_type, bool reinit, Real Epsilon) {
    AKANTU_DEBUG_IN();

    if (ghost_type == _not_ghost)
        residual->clear();

    if (!supg_flag)
        AKANTU_DEBUG_ERROR("You should compute the stiffness matrix at least once.");

    Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
    for (; it != last_type; ++it) {

        if (mesh.getSpatialDimension(*it) != spatial_dimension) continue;

        Array<Real> & shapes = shapes_SUPG(*it, ghost_type);

        Array<UInt> & elem_filter = element_filter(*it, ghost_type);

        UInt nb_element = elem_filter.getSize();
        UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
        UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(*it, ghost_type);

        getFEM().interpolateOnQuadraturePoints(*phi, phi_on_qpoints(*it, ghost_type), 1, *it, ghost_type);


        Array<Real> * shapes_filtered = new Array<Real > (nb_element * nb_quadrature_points,
                nb_nodes_per_element,
                "shapes filtered");

        Array<Real>::iterator< Matrix<Real> > shapes_it = shapes.begin(1, nb_nodes_per_element);

        Array<Real>::iterator< Matrix<Real> > shapes_filtered_it = shapes_filtered->begin(1, nb_nodes_per_element);

        UInt * elem_filter_val = elem_filter.storage();
        for (UInt e = 0; e < nb_element; ++e, ++elem_filter_val)
            for (UInt q = 0; q < nb_quadrature_points; ++q, ++shapes_filtered_it)
                *shapes_filtered_it = shapes_it[*elem_filter_val * nb_quadrature_points + q];



        UInt residual_size = nb_nodes_per_element;

        Array<Real> * phi_Ni = new Array<Real > (nb_element * nb_quadrature_points,
                residual_size,
                "phi_Ni");

        //Compute term phixN_i/Delta_t

        Array<Real>::iterator< Vector<Real> > Ni_it = shapes_filtered->begin(nb_nodes_per_element);

        Array<Real> & phi_qpoints = phi_on_qpoints(*it, ghost_type);

        Array<Real>::iterator< Vector<Real> > phi_qpoints_it = phi_qpoints.begin_reinterpret(1, nb_quadrature_points * nb_element);
        Array<Real>::iterator< Vector<Real> > phi_qpoints_end = phi_qpoints.end_reinterpret(1, nb_quadrature_points * nb_element);

        Array<Real>::iterator< Vector<Real> > phi_Ni_it = phi_Ni->begin(nb_nodes_per_element);

        for (; phi_qpoints_it != phi_qpoints_end; ++phi_qpoints_it, ++phi_Ni_it, ++Ni_it) {
            Real phi_q = (*phi_qpoints_it)[0];
            Real sign = 0;
            if (reinit) {
                sign = sign_phi(phi_q, Epsilon); //phi_q / (sqrt(phi_q * phi_q + 0.01));
                /*if (phi_q > 0)
                    sign = 1;
                else if (phi_q < 0)
                    sign = -1;
                else
                    sign = 0;*/
            }
            Real Filtered_value = 1.0;
            if (filtered_flag) {
                Real sin_phi = 1 - (3.1415927 * phi_q / (2 * Filter_parameter))*(3.1415927 * phi_q / (2 * Filter_parameter));
                if (sin_phi <= 0)
                    Filtered_value = 0.0;
                else
                    Filtered_value = sqrt(sin_phi);
            }



            Vector<Real> & phi_Ni_element = *phi_Ni_it;
            Vector<Real> & Ni_element = *Ni_it;
            for (UInt q = 0; q < nb_nodes_per_element; q++)
                phi_Ni_element[q] = Ni_element[q];
            phi_Ni_element *= (phi_q / time_step + sign * Filtered_value);



        }

        Array<Real> * f = new Array<Real > (nb_element,
                residual_size,
                "f");

        getFEM().integrate(*phi_Ni, *f, residual_size, *it, ghost_type, &elem_filter);

        getFEM().assembleArray(*f, *residual, dof_synchronizer->getLocalDOFEquationNumbers(), 1, *it, ghost_type, NULL, 1); //(*Phi_x_Ni, f, 1, *it, ghost_type, &elem_filter);


        //if (!reinit || (reinit && filtered_flag)) {

        phi_Ni->clear();

        Array<Real> * v_grad_phi = new Array<Real > (1, nb_element * nb_quadrature_points,
                "v_grad_phi");


        Array<Real> & grad_phi = phi_gradient(*it, ghost_type);

        Array<Real> * v_tmp = &v_on_qpoints(*it, ghost_type);
        if (reinit)
            v_tmp = &v_r_on_qpoints(*it, ghost_type);
        Array<Real> & v_qpoints = *v_tmp;


        //grad_phi.resize(nb_quadrature_points * nb_element);

        getFEM().gradientOnQuadraturePoints(*phi, grad_phi, 1, *it, ghost_type, &elem_filter);


        Ni_it = shapes_filtered->begin(nb_nodes_per_element);

        Array<Real>::iterator< Matrix<Real> > grad_phi_qpoints_it = grad_phi.begin(1, spatial_dimension);
        Array<Real>::iterator< Matrix<Real> > grad_phi_qpoints_end = grad_phi.end(1, spatial_dimension);

        Array<Real>::iterator< Matrix<Real> > v_qpoints_it = v_qpoints.begin(1, spatial_dimension);

        Array<Real>::iterator< Matrix<Real> > v_grad_phi_it = v_grad_phi->begin_reinterpret(1, 1, nb_element * nb_quadrature_points);

        phi_Ni_it = phi_Ni->begin(nb_nodes_per_element);
        Ni_it = shapes_filtered->begin(nb_nodes_per_element);

        for (; grad_phi_qpoints_it != grad_phi_qpoints_end; ++grad_phi_qpoints_it, ++v_qpoints_it, ++v_grad_phi_it, ++phi_Ni_it, ++Ni_it) {

            Matrix<Real> & v_qpoints_element = *v_qpoints_it;
            Matrix<Real> & grad_phi_qpoints_element = *grad_phi_qpoints_it;
            Matrix<Real> & v_grad_phi_element = *v_grad_phi_it;
            v_grad_phi_element.mul < false, true > (v_qpoints_element, grad_phi_qpoints_element);

            Real v_grad_phi_value = (*v_grad_phi_it)[0];
            Vector<Real> & phi_Ni_element = *phi_Ni_it;
            Vector<Real> & Ni_element = *Ni_it;
            for (UInt q = 0; q < nb_nodes_per_element; q++)
                phi_Ni_element[q] = Ni_element[q];
            phi_Ni_element *= v_grad_phi_value * (-0.5);
        }


        f->clear();

        getFEM().integrate(*phi_Ni, *f, residual_size, *it, ghost_type, &elem_filter);

        getFEM().assembleArray(*f, *residual, dof_synchronizer->getLocalDOFEquationNumbers(), 1, *it, ghost_type, NULL, 1); //(*Phi_x_Ni, f, 1, *it, ghost_type, &elem_filter);
        delete v_grad_phi;
        //}



        delete shapes_filtered;
        delete f;
        delete phi_Ni;
    }



    AKANTU_DEBUG_OUT();
}

bool LevelSetModel::testConvergenceResidual(Real tolerance, Real & norm) {
    AKANTU_DEBUG_IN();

    UInt nb_nodes = residual->getSize();

    Array<Real> * Ku = new Array<Real > (*phi, true, "Ku");
    *Ku *= *stiffness_matrix;


    norm = 0;
    Real * RHS_val = residual->values;
    Real * residual_val = Ku->values;
    bool * boundary_val = boundary->values;

    for (UInt n = 0; n < nb_nodes; ++n) {
        bool is_local_node = mesh.isLocalOrMasterNode(n);
        if (is_local_node) {

            if (!(*boundary_val)) {
                norm += (*RHS_val - *residual_val) * (*RHS_val - *residual_val);
            }
            boundary_val++;
            residual_val++;
            RHS_val++;

        } else {
            boundary_val++;
            residual_val++;
            RHS_val++;
        }
    }

    StaticCommunicator::getStaticCommunicator().allReduce(&norm, 1, _so_sum);

    norm = sqrt(norm);

    delete Ku;

    AKANTU_DEBUG_ASSERT(!Math::isnan(norm), "Something goes wrong in the solve phase");

    AKANTU_DEBUG_OUT();

    return (norm < tolerance);
}

/* -------------------------------------------------------------------------- */

void LevelSetModel::initFull() {
    //readMaterials(material_file);
    //model initialization
    initModel();
    //initialize the vectors
    initArrays();
    initFEMBoundary();
    phi->clear();
    //v->clear();
    time_step = 0.0;

    initImplicit();
}


/* -------------------------------------------------------------------------- */

/**
 * Initialize the implicit solver
 */
void LevelSetModel::initImplicit(SolverOptions & solver_options) {
    AKANTU_DEBUG_IN();

    initSolver(solver_options);

    std::stringstream sstr;
    sstr << id << ":stiffness_matrix";
    stiffness_matrix = new SparseMatrix(*jacobian_matrix, sstr.str(), memory_id);
    // } else {
    //   stiffness_matrix = jacobian_matrix;
    // }

    AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */

/**
 * Initialize the solver and create the sparse matrices needed.
 *
 */
void LevelSetModel::initSolver(SolverOptions & options) {
#if !defined(AKANTU_USE_MUMPS) // or other solver in the future \todo add AKANTU_HAS_SOLVER in CMake
    AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#else
    UInt nb_global_node = mesh.getNbGlobalNodes();

    std::stringstream sstr;
    sstr << id << ":jacobian_matrix";
    jacobian_matrix = new SparseMatrix(nb_global_node, _unsymmetric,
            1, sstr.str(), memory_id);
    //  dof_synchronizer->initGlobalDOFEquationNumbers();
    jacobian_matrix->buildProfile(mesh, *dof_synchronizer);

#ifdef AKANTU_USE_MUMPS
    std::stringstream sstr_solv;
    sstr_solv << id << ":solver";
    solver = new SolverMumps(*jacobian_matrix, sstr_solv.str());

    dof_synchronizer->initScatterGatherCommunicationScheme();
#else
    AKANTU_DEBUG_ERROR("You should at least activate one solver.");
#endif //AKANTU_USE_MUMPS*/

    solver->initialize(options);
#endif //AKANTU_HAS_SOLVER*/
}

/* -------------------------------------------------------------------------- */
void LevelSetModel::initFEMBoundary(bool create_surface) {

    if (create_surface)
        MeshUtils::buildFacets(getFEM().getMesh());

    FEM & fem_boundary = getFEMBoundary();
    fem_boundary.initShapeFunctions();
    fem_boundary.computeNormalsOnControlPoints();
}

void LevelSetModel::Shapes_SUPG(bool reinit) {

    Array<Real> & position = mesh.getNodes();

    for (UInt g = _not_ghost; g <= _ghost; ++g) {
        GhostType gt = (GhostType) g;
        Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt, _ek_not_defined);
        Mesh::type_iterator end = mesh.lastType(spatial_dimension, gt, _ek_not_defined);


        for (; it != end; ++it) {
            UInt nb_element = mesh.getNbElement(*it, gt);
            UInt nb_quad_points = getFEM().getNbQuadraturePoints(*it, gt);
            UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);

            if (!reinit)
                getFEM().interpolateOnQuadraturePoints(*v, v_on_qpoints(*it, gt), spatial_dimension, *it, gt);

            const Array<Real> & shapes = getFEM().getShapes(*it, gt);

            Array<Real> & supg_qpoints = shapes_SUPG.alloc(nb_element*nb_quad_points,
                    nb_nodes_per_element, *it, gt);

            Array<Real> * v_tmp = &v_on_qpoints(*it, gt);
            if (reinit)
                v_tmp = &v_r_on_qpoints(*it, gt);
            Array<Real> & v_qpoints = *v_tmp;

            const Array<Real> & shapes_derivatives = getFEM().getShapesDerivatives(*it, gt);

            Array<Real>::const_iterator< Matrix<Real> > shapes_derivatives_it = shapes_derivatives.begin(nb_nodes_per_element, spatial_dimension);
            Array<Real>::const_iterator< Matrix<Real> > shapes_it = shapes.begin(nb_nodes_per_element, 1);
            Array<Real>::const_iterator< Matrix<Real> > shapes_derivatives_end = shapes_derivatives.end(nb_nodes_per_element, spatial_dimension);
            Array<Real>::iterator< Matrix<Real> > v_qpoints_it = v_qpoints.begin(1, spatial_dimension);

            Array<Real>::iterator< Matrix<Real> > supg_qpoints_it = supg_qpoints.begin(nb_nodes_per_element, 1);
            
            
            Array<Real> X(0, nb_nodes_per_element * spatial_dimension);
            FEM::extractNodalToElementField(mesh, position, X, *it, _not_ghost);

            Array<Real>::iterator< Matrix<Real> > X_el = X.begin(spatial_dimension, nb_nodes_per_element);



            UInt counter_q = 0;
            Real el_size = 0;

            for (; shapes_derivatives_it != shapes_derivatives_end; ++v_qpoints_it, ++X_el, ++supg_qpoints_it, ++shapes_derivatives_it, ++counter_q, ++shapes_it) {

                const Matrix<Real> & shapes_derivatives_element = *shapes_derivatives_it;
                Matrix<Real> & v_element = *v_qpoints_it;
                Matrix<Real> & supg_element = *supg_qpoints_it;
                const Matrix<Real> & shapes_element = *shapes_it;

                if (counter_q == nb_quad_points) {
                    el_size = 0;
                    counter_q = 0;
                }


                if (!counter_q) {
                    el_size = getFEM().getElementInradius(*X_el, *it);
                }

                Real v_norm = v_element.norm();

                Real Tau_supg = 0.0;
                if (v_norm > 1e-4)
                    Tau_supg = el_size / ((spatial_dimension + 1.0) * v_norm);
                else
                    Tau_supg = 0.0;

                //Tau_supg=0.0;

                supg_element.mul < false, true > (shapes_derivatives_element, v_element, Tau_supg);

                supg_element += shapes_element;


            }
        }
    }
    supg_flag = true;

}

/*void LevelSetModel::ApplyInflowBC(bool matrix, bool reinit) {

    for (UInt g = _not_ghost; g <= _ghost; ++g) {
        GhostType ghost_type = (GhostType) g;
        Mesh::type_iterator it = mesh.firstType(spatial_dimension - 1, ghost_type, _ek_not_defined);
        Mesh::type_iterator end = mesh.lastType(spatial_dimension - 1, ghost_type, _ek_not_defined);

        //Matrix
        for (; it != end; it++) {

            if (mesh.getSpatialDimension(*it) != spatial_dimension - 1) continue;

            SparseMatrix & K = *stiffness_matrix;

            const Array<Real> & shapes = getFEMBoundary().getShapes(*it, ghost_type);

            //Compute \Grad{\phi}

            Array<UInt> & elem_filter = element_filter_boundary(*it, ghost_type);

            UInt nb_element = elem_filter.getSize();
            UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
            UInt nb_quadrature_points = getFEMBoundary().getNbQuadraturePoints(*it, ghost_type);

            Array<Real> * shapes_filtered = new Array<Real > (nb_element * nb_quadrature_points,
                    nb_nodes_per_element,
                    "shapes filtered");

            Array<Real>::const_iterator< Matrix<Real> > shapes_it = shapes.begin(1, nb_nodes_per_element);

            Array<Real>::iterator< Matrix<Real> > shapes_filtered_it = shapes_filtered->begin(1, nb_nodes_per_element);
            Array<Real>::iterator< Matrix<Real> > shapes_filtered_end = shapes_filtered->end(1, nb_nodes_per_element);

            UInt * elem_filter_val = elem_filter.storage();
            for (UInt e = 0; e < nb_element; ++e, ++elem_filter_val)
                for (UInt q = 0; q < nb_quadrature_points; ++q, ++shapes_filtered_it)
 *shapes_filtered_it = shapes_it[*elem_filter_val * nb_quadrature_points + q];


            //Compute Term v x n

            if(!reinit)
                getFEMBoundary().interpolateOnQuadraturePoints(*v, v_on_qpoints_boundary(*it, ghost_type), spatial_dimension, *it, ghost_type);

            Array<Real> * v_tmp = &v_on_qpoints(*it, ghost_type);
            if(reinit)
                v_tmp = &v_r_on_qpoints(*it, ghost_type);
            Array<Real> & v_qpoints = *v_tmp;

            Array<Real> v_n(nb_element*nb_quadrature_points, 1, "v_n");

            const Array<Real> & normals_qpoints = getFEMBoundary().getNormalsOnQuadPoints(*it, ghost_type);

            Array<Real>::iterator< Matrix<Real> > v_qpoints_it = v_qpoints.begin(1, spatial_dimension);
            Array<Real>::const_iterator< Matrix<Real> > normals_qpoints_it = normals_qpoints.begin(1, spatial_dimension);
            Array<Real>::iterator< Matrix<Real> > v_n_it = v_n.begin(1, 1);


            UInt nb_qpoints_elements = nb_element * nb_quadrature_points;
            for (UInt i = 0; i < nb_qpoints_elements; ++i, ++v_qpoints_it, ++normals_qpoints_it, ++v_n_it) {
                Matrix<Real> & v_qpoints_element = *v_qpoints_it;
                const Matrix<Real> & normals_qpoints_element = *normals_qpoints_it;
                Matrix<Real> & v_n_element = *v_n_it;
                v_n_element.mul < false, true > (v_qpoints_element, normals_qpoints_element);
            }

            Real * v_n_qpoint = v_n.values;

            if (matrix) {

                UInt Ni_Nj_size = nb_nodes_per_element;

                Array<Real> * Ni_Nj = new Array<Real > (nb_element * nb_quadrature_points,
                        Ni_Nj_size * Ni_Nj_size,
                        "Ni_Nj");

                //Compute term N_ixN_j/Delta_t

                shapes_filtered_it = shapes_filtered->begin(1, nb_nodes_per_element);

                Array<Real>::iterator< Matrix<Real> > Ni_Nj_it = Ni_Nj->begin(Ni_Nj_size, Ni_Nj_size);

                v_n_qpoint = v_n.values;

                for (; shapes_filtered_it != shapes_filtered_end; ++Ni_Nj_it, ++shapes_filtered_it, ++v_n_qpoint) {

                    Matrix<Real> & shapes_element = *shapes_filtered_it;
                    Matrix<Real> & Ni_Nj_element = *Ni_Nj_it;

                    if (*v_n_qpoint > 0.0)
                        Ni_Nj_element.mul < true, false > (shapes_element, shapes_element, *v_n_qpoint);
                    else
                        Ni_Nj_element.mul < true, false > (shapes_element, shapes_element, 0.0);

                }

                //delete tangent_stiffness_matrix;


                Array<Real> * K_e = new Array<Real > (nb_element,
                        Ni_Nj_size * Ni_Nj_size,
                        "K_e");

                getFEMBoundary().integrate(*Ni_Nj, *K_e, Ni_Nj_size * Ni_Nj_size, *it, ghost_type, &elem_filter);

                //delete Ni_Nj;

                getFEMBoundary().assembleMatrix(*K_e, K, 1, *it, ghost_type, &elem_filter);

                delete Ni_Nj;
                delete K_e;
            }


            //RHS

            getFEMBoundary().interpolateOnQuadraturePoints(*phi, phi_on_qpoints_boundary(*it, ghost_type), 1, *it, ghost_type);

            UInt residual_size = nb_nodes_per_element;

            Array<Real> * phi_Ni = new Array<Real > (nb_element * nb_quadrature_points,
                    residual_size,
                    "phi_Ni");

            //Compute term phixN_i/Delta_t

            Array<Real>::iterator< Vector<Real> > Ni_it = shapes_filtered->begin(nb_nodes_per_element);

            Array<Real> & phi_qpoints = phi_on_qpoints_boundary(*it, ghost_type);

            Array<Real>::iterator< Vector<Real> > phi_qpoints_it = phi_qpoints.begin_reinterpret(1, nb_quadrature_points * nb_element);
            Array<Real>::iterator< Vector<Real> > phi_qpoints_end = phi_qpoints.end_reinterpret(1, nb_quadrature_points * nb_element);

            Array<Real>::iterator< Vector<Real> > v_n_qpoints_it = v_n.begin_reinterpret(1, nb_quadrature_points * nb_element);

            Array<Real>::iterator< Vector<Real> > phi_Ni_it = phi_Ni->begin(nb_nodes_per_element);

            for (; phi_qpoints_it != phi_qpoints_end; ++phi_qpoints_it, ++phi_Ni_it, ++v_n_qpoints_it, ++Ni_it) {

                Real phi_q = (*phi_qpoints_it)[0];
                Real v_n_q = (*v_n_qpoints_it)[0];

                Vector<Real> & phi_Ni_element = *phi_Ni_it;
                Vector<Real> & Ni_element = *Ni_it;
                for (UInt q = 0; q < nb_nodes_per_element; q++)
                    phi_Ni_element[q] = Ni_element[q];
                if (v_n_q < 0)
                    phi_Ni_element *= -1.0 * phi_q * v_n_q;
                else
                    phi_Ni_element *= 0.0;
            }

            Array<Real> * f = new Array<Real > (nb_element,
                    residual_size,
                    "f");

            getFEMBoundary().integrate(*phi_Ni, *f, residual_size, *it, ghost_type, &elem_filter);

            getFEMBoundary().assembleArray(*f, *residual, dof_synchronizer->getLocalDOFEquationNumbers(), 1, *it, ghost_type, NULL, 1); //(*Phi_x_Ni, f, 1, *it, ghost_type, &elem_filter);

            delete shapes_filtered;
            delete f;
            delete phi_Ni;

        }


    }
    AKANTU_DEBUG_IN();
}*/

/* -------------------------------------------------------------------------- */
void LevelSetModel::assemblePhi(const GhostType & ghost_type, bool reinit) {
    AKANTU_DEBUG_IN();

    if (ghost_type == _not_ghost)
        stiffness_matrix->clear();

    if (!supg_flag) {
        Shapes_SUPG(reinit);
    }


    Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, ghost_type);
    for (; it != last_type; ++it) {

        if (mesh.getSpatialDimension(*it) != spatial_dimension) continue;

        SparseMatrix & K = *stiffness_matrix;

        const Array<Real> & shapes_derivatives = getFEM().getShapesDerivatives(*it, ghost_type);
        const Array<Real> & shapes = getFEM().getShapes(*it, ghost_type);
        Array<Real> & supg = shapes_SUPG(*it, ghost_type);


        //Compute \Grad{\phi}

        Array<UInt> & elem_filter = element_filter(*it, ghost_type);

        UInt nb_element = elem_filter.getSize();
        UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(*it);
        UInt nb_quadrature_points = getFEM().getNbQuadraturePoints(*it, ghost_type);

        Array<Real> * shapes_derivatives_filtered = new Array<Real > (nb_element * nb_quadrature_points,
                spatial_dimension * nb_nodes_per_element,
                "shapes derivatives filtered");

        Array<Real> * shapes_filtered = new Array<Real > (nb_element * nb_quadrature_points,
                nb_nodes_per_element,
                "shapes filtered");

        Array<Real> * supg_filtered = new Array<Real > (nb_element * nb_quadrature_points,
                nb_nodes_per_element,
                "supg filtered");


        Array<Real>::const_iterator< Matrix<Real> > shapes_derivatives_it = shapes_derivatives.begin(spatial_dimension,
                nb_nodes_per_element);
        Array<Real>::const_iterator< Matrix<Real> > shapes_it = shapes.begin(1, nb_nodes_per_element);

        Array<Real>::iterator< Matrix<Real> > supg_it = supg.begin(1, nb_nodes_per_element);

        Array<Real>::iterator< Matrix<Real> > shapes_derivatives_filtered_it = shapes_derivatives_filtered->begin(spatial_dimension,
                nb_nodes_per_element);
        Array<Real>::iterator< Matrix<Real> > shapes_filtered_it = shapes_filtered->begin(1, nb_nodes_per_element);

        Array<Real>::iterator< Matrix<Real> > supg_filtered_it = supg_filtered->begin(1, nb_nodes_per_element);

        UInt * elem_filter_val = elem_filter.storage();
        for (UInt e = 0; e < nb_element; ++e, ++elem_filter_val) {
            for (UInt q = 0; q < nb_quadrature_points; ++q, ++shapes_derivatives_filtered_it, ++shapes_filtered_it, ++supg_filtered_it) {
                *shapes_derivatives_filtered_it = shapes_derivatives_it[*elem_filter_val * nb_quadrature_points + q];
                *shapes_filtered_it = shapes_it[*elem_filter_val * nb_quadrature_points + q];
                *supg_filtered_it = supg_it[*elem_filter_val * nb_quadrature_points + q];
            }
        }

        UInt Ni_Nj_size = nb_nodes_per_element;

        Array<Real> * Ni_Nj = new Array<Real > (nb_element * nb_quadrature_points,
                Ni_Nj_size * Ni_Nj_size,
                "Ni_Nj");

        //Compute term N_ixN_j/Delta_t

        shapes_filtered_it = shapes_filtered->begin(1, nb_nodes_per_element);
        Array<Real>::iterator< Matrix<Real> > shapes_filtered_end = shapes_filtered->end(1, nb_nodes_per_element);

        supg_filtered_it = supg_filtered->begin(1, nb_nodes_per_element);

        Array<Real>::iterator< Matrix<Real> > Ni_Nj_it = Ni_Nj->begin(Ni_Nj_size, Ni_Nj_size);



        for (; shapes_filtered_it != shapes_filtered_end; ++Ni_Nj_it, ++shapes_filtered_it, ++supg_filtered_it) {
            Matrix<Real> & shapes_element = *shapes_filtered_it;
            Matrix<Real> & supg_element = *supg_filtered_it;
            Matrix<Real> & Ni_Nj_element = *Ni_Nj_it;

            Ni_Nj_element.mul < true, false > (supg_element, shapes_element, 1.0 / time_step);
        }

        //delete tangent_stiffness_matrix;

        Array<Real> * K_e = new Array<Real > (nb_element,
                Ni_Nj_size * Ni_Nj_size,
                "K_e");

        getFEM().integrate(*Ni_Nj, *K_e, Ni_Nj_size * Ni_Nj_size, *it, ghost_type, &elem_filter);

        //delete Ni_Nj;

        getFEM().assembleMatrix(*K_e, K, 1, *it, ghost_type, &elem_filter);






        //Compute Term v x \grad{phi} x N

        Array<Real> * v_tmp = &v_on_qpoints(*it, ghost_type);
        if (reinit)
            v_tmp = &v_r_on_qpoints(*it, ghost_type);
        Array<Real> & v_qpoints = *v_tmp;

        Array<Real>::iterator< Matrix<Real> > grad_phi_it = shapes_derivatives_filtered->begin(nb_nodes_per_element, spatial_dimension);
        Array<Real>::iterator< Matrix<Real> > grad_phi_end = shapes_derivatives_filtered->end(nb_nodes_per_element, spatial_dimension);

        Array<Real>::iterator< Matrix<Real> > v_qpoints_it = v_qpoints.begin(1, spatial_dimension);

        Array<Real> * v_x_grad_phi = new Array<Real > (nb_element * nb_quadrature_points, nb_nodes_per_element, "v_x_grad_phi");
        Array<Real>::iterator< Matrix<Real> > v_x_grad_phi_it = v_x_grad_phi->begin(nb_nodes_per_element, 1);

        for (; grad_phi_it != grad_phi_end; ++v_x_grad_phi_it, ++grad_phi_it, ++v_qpoints_it) {
            Matrix<Real> & grad_phi_element = *grad_phi_it;
            Matrix<Real> & v_element = *v_qpoints_it;
            Matrix<Real> & v_x_grad_phi_element = *v_x_grad_phi_it;

            v_x_grad_phi_element.mul < false, true > (grad_phi_element, v_element);
        }

        Ni_Nj->clear();

        //Compute term N_ixN_j/Delta_t

        supg_filtered_it = supg_filtered->begin(1, nb_nodes_per_element);
        Array<Real>::iterator< Matrix<Real> > supg_filtered_end = supg_filtered->end(1, nb_nodes_per_element);

        //supg_filtered_it = shapes_filtered->begin(1, nb_nodes_per_element); //TEST
        //Array<Real>::iterator< Matrix<Real> > supg_filtered_end = shapes_filtered->end(1, nb_nodes_per_element); //TEST

        Ni_Nj_it = Ni_Nj->begin(Ni_Nj_size, Ni_Nj_size);
        v_x_grad_phi_it = v_x_grad_phi->begin(1, nb_nodes_per_element);

        Real Value_term = 0.5;
        //if(reinit && !filtered_flag)
        //Value_term=1.0;

        for (; supg_filtered_it != supg_filtered_end; ++Ni_Nj_it, ++supg_filtered_it, ++v_x_grad_phi_it) {

            Matrix<Real> & supg_element = *supg_filtered_it;
            Matrix<Real> & v_x_grad_phi_element = *v_x_grad_phi_it;
            Matrix<Real> & Ni_Nj_element = *Ni_Nj_it;


            Ni_Nj_element.mul < true, false > (supg_element, v_x_grad_phi_element, Value_term);
            //Ni_Nj_element.mul < true, false > (v_x_grad_phi_element, supg_element, -1.0); //TEST
        }

        //delete tangent_stiffness_matrix;

        K_e->clear();

        getFEM().integrate(*Ni_Nj, *K_e, Ni_Nj_size * Ni_Nj_size, *it, ghost_type, &elem_filter);

        getFEM().assembleMatrix(*K_e, K, 1, *it, ghost_type, &elem_filter);

        //delete tangent_stiffness_matrix;

        delete shapes_derivatives_filtered;
        delete shapes_filtered;
        delete supg_filtered;
        delete Ni_Nj;
        delete K_e;
        delete v_x_grad_phi;

    }

    AKANTU_DEBUG_OUT();
}

void LevelSetModel::addDumpField(const std::string & field_id) {
#ifdef AKANTU_USE_IOHELPER
#define ADD_FIELD(field, type)						\
  addDumpFieldToDumper(BOOST_PP_STRINGIZE(field),new DumperIOHelper::NodalField<type>(*field))

    if (field_id == "phi") {
        ADD_FIELD(phi, Real);
    } else if (field_id == "velocity") {
        ADD_FIELD(v, Real);
    } else if (field_id == "residual") {
        ADD_FIELD(residual, Real);
    } else if (field_id == "boundary") {
        ADD_FIELD(boundary, bool);
    }
    //else if(field_id == "partitions"  ) {
    //  addDumpFieldToDumper(field_id,
    //			 new DumperIOHelper::ElementPartitionField(mesh,
    //								   spatial_dimension,
    //								   _not_ghost,
    //								   _ek_regular));
    //  }
#undef ADD_FIELD
#endif
}


__END_AKANTU__
