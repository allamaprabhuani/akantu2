Contact Mechanics Internodes Model
==================================

The contact mechanics internodes model is a implementation of
:cpp:class:`Model <akantu::Model>` interface to handle contact between
two solid bodies using the INTERNODES (reference) method.

Overview INTERNODES method
--------------------

The INTERNODES method in the context of cantact mechanics is algorithm to solve the
contact problem using nodes belonging to the interface. Information such as
displacement is transfered from one body to the other with a interpolation
using rescaled radial basis function. 

The initial configuration should already exhibit some interprenation to start
the algorithm. A set of predefined possible contact nodes of both bodies will
be taken in consideration to solve the problem. The algorithm will solve the
intenodes method.


Using the contact mechanics internodes model
--------------------------------------------

The :cpp:class:`ContactMechanicsInternodesModel
<akantu::ContactMechanicsInternodesModel>` works simalarly to the
:cpp:class:`ContactMechanicssModel <akantu::ContactMechanicsInternodesModel>`.
The class can be instanciated with an existing (see \ref{sect:common:mesh}) as follows::

   ContactMechanicsInternodesModel contact_model(mesh, spatial_dimension);

In contrast to the the standard :cpp:class:`ContactMechanicssModel
<akantu::ContactMechanicsInternodesModel>` the INTERNODES model handles the
contact resolution with the :cpp:class:`ContactDetectorInternodes
<akantu::ContactDetectorInternodes>` and the coupling with a
solid model. Currently this method is fixed to use
:cpp:class:`SolidMechanicsModel <akantu::SolidMechanicsModel>`.


The initilization can be done with::

  contact_model.initFull(_analysis_method =_static);

The code belows gives the associated detector and solid model::

  detector = contact_model.getContactDetector()
  solid = contact_model.getSolidMechanicsModel()

The boundary conditions can then be set directly through the solid model.


To initizial the contact detection of the model a configuration section can be
added to the 'material.dat' file. The possible nodes for the two interfaces can
be specified as follows.

.. code-block::

   contact_detector [
     master = contact_bottom
     slave = contact_top
   ]


Currently, one iteration of the INTERNODES method is possible::

.. code-block:: c++

   contact_model.solveStep();
