// Element size
h = 5;

// Dimension of square
L = 1;

// ------------------------------------------
// Geometry
// ------------------------------------------

// Points
Point(1) = {0, 0, 0, h};
Point(2) = {-2*L, 0.5*L, 0, h};
Point(3) = {-2*L, -0.5*L, 0, h};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};

// Geometric and Physical Surface
Line Loop(1) = {1, 2, 3};

Plane Surface(1) = {1};

