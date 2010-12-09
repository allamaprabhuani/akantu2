// 2D MESH WITH TWO SQUARES

// element size for first square
h1 = 0.2;

// element size for second square
h2 = 0.11;

// -------- //
// GEOMETRY //
// -------. //

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

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
