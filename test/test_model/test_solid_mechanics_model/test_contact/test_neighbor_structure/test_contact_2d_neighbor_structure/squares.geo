// 2D MESH WITH 3 SQUARES

// element size for first square
h1 = 0.2;

// element size for second square
h2 = 0.166;

// element size for second square
h3 = 0.125;

// -------- //
// GEOMETRY //
// -------- //

// square 1
Point(1) = {0, -1, 0, h1};
Point(2) = {2, -1, 0, h1};
Point(3) = {2, 0, 0, h1};
Point(4) = {0, 0, 0, h1};

// square 2
Point(5) = {0.75, 0.05, 0, h2};
Point(6) = {1.25, 0.05, 0, h2};
Point(7) = {1.25, 0.55, 0, h2};
Point(8) = {0.75, 0.55, 0, h2};

// square 3
Point(9) = {1.75, 0.3, 0, h3};
Point(10) = {2.25, -0.2, 0, h3};
Point(11) = {2.75, 0.3, 0, h3};
Point(12) = {2.25, 0.8, 0, h3};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};
Line Loop(3) = {9, 10, 11, 12};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
