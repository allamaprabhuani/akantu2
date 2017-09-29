/* 
Copyright 2008 Guillaume ANCIAUX (guillaume.anciaux@epfl.ch)

This file is part of ParaViewHelper.

ParaViewHelper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ParaViewHelper is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ParaViewHelper.  If not, see <http://www.gnu.org/licenses/>.

*/

extern "C" {
  extern void BuildGridMesh(double ** pos,int ** conn, 
		   int * nb_nodes,int * nb_elements,
		   double L, double D, double W,
		   double NL,double ND, double NW);
}

#include "../src/dumper_paraview.h"

int main(){

  double * positions;
  int * connectivity;

  int nb_nodes;
  int nb_elements;
  int D=1,W=1,L=1,ND=1,NL=1,NW=1;

  BuildGridMesh(&positions,&connectivity,&nb_nodes,&nb_elements,
		L,D,W,NL,ND,NW);
  

  DumperParaview dumper;
  dumper.SetPoints(positions,3,nb_nodes,"cube-mesh");
  dumper.SetConnectivity(connectivity,TETRA2,nb_elements,FORTRAN_MODE);
  dumper.AddNodeDataField(positions,3,"positions");

  dumper.SetMode(TEXT);
  dumper.Init();

  dumper.Dump();
}
