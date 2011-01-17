// MESH WITH TWO SQUARES

// Dimension of squares
L = 1;

// Distance between the two bodies
g = 0.01;

// Element size
h = 0.25;   // generic
h1 = 0.1; // contact surface of bottom square
h2 = 0.07; // contact surface of top square

// Points
Point(1) = {0, -L, 0, h};
Point(2) = {L, -L, 0, h};
Point(3) = {L, 0, 0, h1};
Point(4) = {0, 0, 0, h1};

Point(5) = {0.02*L, 0.+g, 0, h2};
Point(6) = {0.98*L, 0.+g, 0, h2};
Point(7) = {0.98*L, L, 0, h};
Point(8) = {0.02*L, L, 0, h};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Geometric and Physical Surface
Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
