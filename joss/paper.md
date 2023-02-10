---
title: 'Akantu: an HPC finite-element library for contact and fracture simulations'
tags:
  - C++
  - cohesive element
  - contact
  - fracture
  - python

authors:
  - name: Nicolas Richart
    orcid: 0000-0002-1463-4405
    affiliation: 1
  - name: Guillaume Anciaux
    orcid: 0000-0002-9624-5621
    affiliation: 1
  - name: Emil Gallyamov
    affiliation: 1
  - name: Lucas Frérot
    orcid: 0000-0002-4138-1052
    affiliation: "1, 2"
  - name: David Kammer
    orcid: 0000-0003-3782-9368
    affiliation: "1, 3"
  - name: Mohit Pundir
    affiliation: "1, 3"
  - name: Marco Vocialta
    affiliation: 1
  - name: Aurelia Cuba Ramos
    affiliation: 1
  - name: Mauro Corrado
    affiliation: "1, 4"
  - name: Fabian Barras
    orcid: 0000-0003-1109-0200
    affiliation: "1, 5"
  - name: Jean-François Molinari
    orcid: 0000-0002-1728-1844
    affiliation: 1

affiliations:
  - name: Civil Engineering Institute, École Polytechnique Fédérale de Lausanne, Switzerland
    index: 1
  - name: Department of Microsystems Engineering, Univeristy of Freiburg, Germany
    index: 2
  - name: Institute for Building Materials, ETH Zurich, Switzerland
    index: 3
  - name: Department of Structural, Geotechnical and Building Engineering, Politecnico di Torino, Italy
    index: 4
  - name: The Njord Centre Department of Physics, Department of Geosciences, University of Oslo, Norway
    index: 5

date: February 2023
bibliography: paper.bib
---

# Summary
Complex, nonlinear and transient phenomena are at the heart of modern research
in mechanics of materials. For example, the buildup and release of elastic
energy at geological fault is what causes earthquakes, and the intricate details
of the slip zone, the propagation of slip fronts and waves radiated through the
various geological media are still active areas of research
[@kammer_propagation_2012;@kammer_existence_2014]. Similarly, understanding
fracture in heterogeneous materials such as concrete, masonry or ceramics
necessitates the modeling of interaction of crack fronts with complex materials
[@taheri_mousavi_dynamic_2015;@yilmaz_damage_2017;@cuba_ramos_hpc_2018], the
representation of residual shear stresses in the contact of newly-formed crack
surfaces [@zhang_micro-mechanical_2017;@pundir_coupling_2021], and the accurate
characterization of transient dynamics
[@vocialta_numerical_2018;@corrado_effects_2016] and material structure
evolution [@cuba_ramos_hpc_2018;@gallyamov_multi-scale_2020].

The finite-element method is now ubiquitous in virtually all areas of solid
mechanics. With meticulous care on code architecture and performance, we show
with our finite-element library Akantu that it can handle the requirements
mentioned above for state-of-the-art research in mechanics of materials. Akantu
is designed from the ground up for high-performance, highly distributed
computations, while retaining the necessary flexibility to handle:

- crack propagation with cohesive elements
- non-local damage models
- plastic and visco-plastic constitutive laws
- large deformations
- contact constraints (including friction)
- structural elements (beams and shells)
- one-dimensional elements embedded in a three-dimensional mesh (e.g.
  reinforcements in concrete)
- interaction between contact and cohesive elements (residual crack shear
  strength)

# Statement of need
Understanding the interplay between material constitutive behavior and interface
processes such as crack propagation, contact and friction is fundamental to the
study of, among others, earthquakes, concrete structures, ceramics, and
poly-crystalline failure. Thanks to its versatility, the finite-element method
(FEM) has become a staple in these areas. However, codes that can handle cutting
edge simulations with interaction of material behavior and interface processes
in a high-performance computing (HPC) setting are rare, particularity in the
open-source space. Driving these state-of-the-art research simulations to the
exascale era is the primary *raison d'être* of Akantu.

At its heart, Akantu leverages a SOA (structure of array) architecture in order
to take advantage of an object-oriented high-level abstraction, and maintain
performance in the critical areas of the code. In addition Akantu, benefits from
distributed memory parallelization, and on the contrary to many finite-element
codes, Akantu has a element-centric parallelization that makes it easier to
implement algorithms like the dynamic insertion of extrinsic cohesive elements.


# Scaling analysis
![Time to solution with and without cohesive
insertion.\label{fig:tts}](results/TTS.svg) High performance and
scalability is a necessity for the resolution of fracture and contact
simulations. To illustrate the possibilities offered by Akantu, a 3D simulation is presented
where a cube composed of 4'392'180 tetrahedra and 734'594 nodes is being compressed and sheared.
This simulation only serves to demonstrate how Akantu behaves in a situation
where massive fragmentation takes place: about 460'000 cohesive elements are inserted during the run.
This simulation was run on 1 up to 720 cores, on a cluster composed of Intel Xeon nodes with 2 sockets of 36
cores, 512Gb of RAM and dual 25Gb Ethernet links. The time to solution (TTS)
averaged over six different runs is computed for each core count on the $x$-axis in Figure \autoref{fig:tts}.
The overhead due to cohesive element insertion is also highlighted by providing the timings when cohesive element insertions
are precluded.

# Publications
The following publications have been made possible with Akantu:

- @kammer_propagation_2012
- @kammer_existence_2014
- @wolff_non-local_2014
- @richart_implementation_2015
- @taheri_mousavi_dynamic_2015
- @cuba_ramos_new_2015
- @radiguet_role_2015
- @vocialta_influence_2015
- @corrado_effects_2016
- @kammer_length_2016
- @svetlizky_properties_2016
- @vocialta_3d_2016
- @yilmaz_mesoscale_2017
- @yilmaz_damage_2017
- @zhang_micro-mechanical_2017
- @cuba_ramos_hpc_2018
- @vocialta_numerical_2018
- @yilmaz_influence_2018
- @zhang_numerical_2018
- @zhang_numerical_2019
- @frerot_fourier_2019
- @gallyamov_multi-scale_2020
- @milanese_mechanistic_2020
- @albertini_three-dimensional_2021
- @brun_hybrid_2021
- @rezakhani_meso-scale_2021
- @pundir_coupling_2021

# Acknowledgement

The development of Akantu would not have been possible without the support of
the European Research Council ERCstg UFO-240332, the Swiss Federal Office of
Energy contract No. SI/500852-01, Swiss National Science Foundation grant number
CRSII5_17108

# References
