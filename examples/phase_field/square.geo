// Mesh size
h  = 1.0;

// Dimensions of the square
Lx = 1.0;
Ly = 1.0;

// ------------------------------------------
// Geometry
// ------------------------------------------
Point(1) = { 0.0, 0.0, 0.0, h};
Point(2) = { Lx,  0.0, 0.0, h};
Point(3) = { Lx,  Ly, 0.0,  h};
Point(4) = { 0.0, Ly, 0.0,  h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1:4};

Plane Surface(1) = {1};

Physical Surface(1) = {1};
Physical Line("bottom")  = {1};
Physical Line("left")    = {4};
Physical Line("right")   = {2};
Physical Line("top")     = {3};
