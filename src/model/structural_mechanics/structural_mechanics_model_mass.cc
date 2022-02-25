/**
 * @file   structural_mechanics_model_mass.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jul 07 2014
 * @date last modification: Thu Mar 04 2021
 *
 * @brief  function handling mass computation
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "integrator_gauss.hh"
#include "material.hh"
#include "shape_structural.hh"
#include "structural_mechanics_model.hh"

#include <aka_math.hh>
#include <cmath>

/* -------------------------------------------------------------------------- */
namespace akantu {

class ComputeRhoFunctorStruct {
public:
  explicit ComputeRhoFunctorStruct(const StructuralMechanicsModel & model)
      : model(model){};

  void operator()(Matrix<Real> & rho, const Element & element) const {
#if defined(AKANTU_STRUCTURAL_MECHANICS_OLD_MASS_COMPUTATION)
    /* Here we assume use the old mass computation.
     * We assume that `rho` has units of mass per unit length,
     * and is thus not a real density */
#   warning "Using the old mass computation for the Bernoulli beams."
    const auto& mat = model.getMaterialByElement(element);
    const Real mat_mass_per_length = mat.rho; 	//In the old mode, we assume that `rho` has the right units.

#else
    /* Here we assume that `rho` has real density units, i.e. `mass per volume`.
     * Thus to get a mass we must multiply it with the cross section. */
#   warning "Using the new mode for mass computation for the Bernulli beams."
    const auto& mat = model.getMaterialByElement(element);
    const Real mat_rho = mat.rho;		//This is the density
    const Real mat_A   = mat.A;			//This is the corss section
    const Real mat_mass_per_length = mat_rho * mat_A;	//Mass of the beam per unit length
#endif

    rho.set(mat_mass_per_length);  //The integrator _recquiers_ mass per unit length
  }

private:
  const StructuralMechanicsModel & model;
};

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleMassMatrix() {
  AKANTU_DEBUG_IN();

  if (not need_to_reassemble_mass) {
    return;
  }

  if (not getDOFManager().hasMatrix("M")) {
    getDOFManager().getNewMatrix("M", getMatrixType("M"));
  }

  this->getDOFManager().zeroMatrix("M");
  assembleMassMatrix(_not_ghost);

  //also update the lumped mass, if pressent
  if(this->mass)
  {
#   warning "Is it right to take the `_not_ghost`?"
    //It is also taken below
    this->computeShadyLumpedMass(_not_ghost);
  }

  need_to_reassemble_mass = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleMassMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  auto & fem = getFEEngineClass<MyFEEngineType>();
  ComputeRhoFunctorStruct compute_rho(*this);

  for (auto type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_structural)) {
    fem.assembleFieldMatrix(compute_rho, "M", "displacement",
                            this->getDOFManager(), type, ghost_type);
  }
  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */


namespace {

//Computes the distance between two vectors, which are supplied as pointers, in 2D
inline
Real
myDistance2(
	const Real* const __restrict__  	v1,
	const Real* const __restrict__ 		v2)
{
	using std::sqrt;

	const Real v1_0 = v1[0];
	const Real v2_0 = v2[0];
	const Real v1_1 = v1[1];
	const Real v2_1 = v2[1];

	const Real diff_0 = v1_0 - v2_0;
	const Real diff_1 = v1_1 - v2_1;

	const Real diff2_0 = diff_0 * diff_0;
	const Real diff2_1 = diff_1 * diff_1;

	const Real sum2 = diff2_0 + diff2_1;

	return sqrt(sum2);
};


//Computes the distance between two vectors, which are supplied as pointers, in 3
inline
Real
myDistance3(
	const Real* const __restrict__  	v1,
	const Real* const __restrict__ 		v2)
{
	using std::sqrt;

	const Real v1_0 = v1[0];
	const Real v2_0 = v2[0];
	const Real v1_1 = v1[1];
	const Real v2_1 = v2[1];
	const Real v1_2 = v1[2];
	const Real v2_2 = v2[2];

	const Real diff_0 = v1_0 - v2_0;
	const Real diff_1 = v1_1 - v2_1;
	const Real diff_2 = v1_2 - v2_2;

	const Real diff2_0 = diff_0 * diff_0;
	const Real diff2_1 = diff_1 * diff_1;
	const Real diff2_2 = diff_2 * diff_2;

	const Real sum2 = diff2_0 + diff2_1;
	const Real sum3 = sum2 + diff2_2;

	return sqrt(sum3);
};


template<class Array_t>
inline
Real
myDistance(
	const Array_t& 		v1,
	const Array_t& 		v2)
{
	AKANTU_DEBUG_ASSERT(v1.size() == v2.size(), "The two vectors did not have the same size; size(v1) = " << v1.size() << ", size(v2) = " << v2.size());
	switch(v1.size() )
	{
	  case 2:  return myDistance2(v1.storage(), v2.storage());
	  case 3:  return myDistance3(v1.storage(), v2.storage());

	  default:
	  	AKANTU_DEBUG_ASSERT(false, "Unkown dimension " << v1.size());
		return NAN;
	}; //end switch
};
}; //End anonymous namespace


void
StructuralMechanicsModel::computeShadyLumpedMass(
	GhostType 	ghost_type)
{
	using std::pow;

	//Test if it can be stored.
	if(not this->mass)
	   { AKANTU_EXCEPTION("Can not compute the shady lumped mass, the `mass` member is not allocated."); };

	//Preparing
	Array<Real>&	LumpedMass = *this->mass;
	const Mesh&	MyMesh     = this->mesh;
	const UInt 	nbNodes    = MyMesh.getNbNodes();
	const UInt 	nbDim      = MyMesh.getSpatialDimension();
	const auto& 	Nodes      = MyMesh.getNodes();
	const auto&     NodeIT     = Nodes.begin(nbDim);	//Iterator pointing to the beginning of the node range

	AKANTU_DEBUG_ASSERT(nbNodes == LumpedMass.size(), "There is something wrong in the allocation.");		//make some tests
	
	LumpedMass.set(0);		//Set the lumped mass to zero


	/* We now compute the mass and the volume, but not the inertia
	 */
	for(auto elType : MyMesh.elementTypes(nbDim, ghost_type, _ek_structural))
	{
		const Array<UInt>&  elMatMap 	= this->element_material(elType, ghost_type);		//get the material map
		const Array<UInt>&  Conns     	= MyMesh.getConnectivity(elType, ghost_type);		//The connectivity
		const UInt    	    nbConns   	= Conns.size();
		  AKANTU_DEBUG_ASSERT(elMatMap.size() == nbConns, "Error in the element map detected.");

		//Now iterate through all the connectivity
		for(UInt it = 0; it != nbConns; ++it)
		{
			const UInt 	node1    = Conns(it, 0);		//Nodes forming the connection
			const UInt 	node2    = Conns(it, 1);
			const UInt      matID    = elMatMap(it);		//Which material it is
			const auto&     material = this->materials.at(matID);	//the material object.
			const Real      bRho     = material.rho;
			const Real      bXSurf   = material.A;
			  AKANTU_DEBUG_ASSERT(node1 < nbNodes, "Node 1 (" << node1 << ") is out of range (" << nbNodes << ")");
			  AKANTU_DEBUG_ASSERT(node1 < nbNodes, "Node 2 (" << node2 << ") is out of range (" << nbNodes << ")");

			const auto node1_it = NodeIT.operator+(node1);	//Get iterators that points to the correct position
			const auto node2_it = NodeIT.operator+(node2);

			const Vector<Real>&  nCoord1 = *node1_it;	//Get a representation of the node's coordinates
			const Vector<Real>&  nCoord2 = *node2_it;

			const Real bLength = myDistance(nCoord1, nCoord2);	//get the length of the beam

			const Real BVol   = bLength * bXSurf;		//Volume of the beam
	
#			if defined(AKANTU_STRUCTURAL_MECHANICS_CONSISTENT_LUMPED_MASS)
			/* In the original calculation of the consistent mass matrix `bRho` is assumed to be
			 * a "mass per unit length", which means [kg / m^2].
			 * Which results in shitty units, i.e. the entries inside the consisent mass 
			 * matrix have different units.
			 * Howeverm they become independent of the corssectional area.
			 */
#			warning "StructuralMechanicsModel: Ignoring cross sectional area for lumped mass."
			const Real BMass  = bLength * bRho;		//Consistent mass of the beam

#			else
			/* A bit more logical, at least in my view, `bRho` should be a density,
			 * this also results in proper units, at least from my perspective.
			 * However, contrary to the consistent mass matrix, the mass now depends on 
			 * the corss sectional area.
			 */
#			warning "StructuralMechanicsModel: Considering cross sectional area for lumped mass."
			const Real BMass  = BVol * bRho;		//Mass of the beam

#			endif

			//Now distribute the mass at the rigth places
			for(UInt d = 0; d != nbDim; ++d)
			{
				LumpedMass(node1, d) += 0.5 * BMass;
				LumpedMass(node2, d) += 0.5 * BMass;
			}; //end for(d):
			
			//Now we sum up the volume
			// We will only do it at one particular place
			// and handle the 3D case later
			LumpedMass(node1, nbDim) += 0.5 * BVol;	//this entry never point to mass
			LumpedMass(node2, nbDim) += 0.5 * BVol;
		}; //End for(it): going through all the connections
	
	}; //end for(elType):
	
	/*
	 * We now compute the inertia.
	 */
	if(nbDim == 2)
	{
		/* This is the 2D case, so we are assuming that the mass is inside a disc.
		 * Which is given as
		 * \begin{align}
		 * 	I := m \cdot r^2
		 * \end{align}
		 *
		 * The radius is obtained by assuming that the volume, that is associated to a beam,
		 * forms a uniformly disc.
		 * From this volume we can compute the radius.
		 */
		for(UInt it = 0; it != nbNodes; ++it)
		{
			const Real  nVol  = LumpedMass(it, nbDim);	//By definition the volume
			const Real  nMass = LumpedMass(it, 0);		//The mass of the node
			const Real  r2    = nVol / M_PI;		//The square of the volume

			const Real Inertia = nMass * r2;
			LumpedMass(it, nbDim) = Inertia;
		}; //End for(it): going through all the nodes
	}
	else if(nbDim == 3)
	{
		/* This is essentially the same as for 2D only that we assume here,
		 * that the mass is uniformly distributed iniside a sphere.
		 * And thus we have
		 * \begin{align}
		 * 	I := \frac{2}{5} m \cdot r^2
		 * \end{align}
		 *
		 * We also have to set three values, for the three axis.
		 * But since we assume a sphere, they are all the same.
		 */

		const UInt nbMassComponents = LumpedMass.getNbComponent();

		for(UInt it = 0; it != nbNodes; ++it)
		{
			const Real  nVol  = LumpedMass(it, nbDim);	//Associated volume
			const Real  nMass = LumpedMass(it, 0);		//The volume; is the same for all elements
			const Real  r2    = pow((nVol * 3.) / (4 * M_PI), 2. / 3.);    

			const Real  Inertia = (2. / 5) * nMass * r2;

			//now setting the remaining components
			for(UInt d = nbDim; d != nbMassComponents; ++d)
			{
				LumpedMass(it, d) = Inertia;
			}; //end for(d):
		}; //End for(it):
	}
	else
	{
		AKANTU_EXCEPTION("The dimensionsion count " << nbDim << " is not known.");
	};

	return;
}; //End: computeShadyLumpedMass

} // namespace akantu
