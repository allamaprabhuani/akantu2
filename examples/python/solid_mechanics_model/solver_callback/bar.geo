// Mesh size
h  = 0.05;

h1 = h;
h2 = h;

// Dimensions of the bar
Lx = 10;
Ly = 1;

// ------------------------------------------
// Geometry
// ------------------------------------------

Point(101) = { 0.0, -Ly/2, 0.0, h1};
Point(102) = { Lx,  -Ly/2, 0.0, h2};

Point(103) = { Lx,  Ly/2., 0.0,  h2};
Point(104) = { 0.0, Ly/2., 0.0,  h1};

Line(201) = {101, 102};
Line(202) = {102, 103};
Line(203) = {103, 104};
Line(204) = {104, 101};

Line Loop(301) = {201, 202, 203, 204};
Plane Surface(301) = {301};

Transfinite Surface "*";
Recombine Surface "*";
Physical Surface("Body") = {301};

Physical Line("YBlocked") = {201, 203};
Physical Line("Right") = {202};
Physical Line("Left") = {204};
