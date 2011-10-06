// Element size
h = 0.2;
h1 = 0.1;  

// Dimension of square
L = 1;
L1 = L/4.;

// ------------------------------------------
// Geometry
// ------------------------------------------

// Points
Point(1) = {0, 0, 0, h};
Point(2) = {L/2, 0, 0, h};
Point(3) = {L/2, L, 0, h};
Point(4) = {0, L, 0, h};

Point(5) = {1.01*L/2,    L/2-L1/2, 0, h1};
Point(6) = {1.01*L/2+L1, L/2-L1/2, 0, h1};
Point(7) = {1.01*L/2+L1, L/2+L1/2, 0, h1};
Point(8) = {1.01*L/2,    L/2+L1/2, 0, h1};

// Lines
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

Line(5) = {8, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};

// Geometric and Physical Surface
Line Loop(9) = {1, 2, 3, 4};
Plane Surface(9) = {9};

//Transfinite Surface {9};
//Recombine Surface {9};

Line Loop(10) = {5, 6, 7, 8};
Plane Surface(10) = {10};

//Transfinite Surface {10};
//Recombine Surface {10};
