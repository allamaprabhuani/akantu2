// Element size
h = 0.2;
h1 = 0.06;

// Dimension of square
L = 1;
L1 = L/4.;

// ------------------------------------------
// Geometry
// ------------------------------------------

// Points
Point(1) = {0, 0, L/2, h};
Point(2) = {L/2, 0, L/2, h};
Point(3) = {L/2, L, L/2, h};
Point(4) = {0, L, L/2, h};

Point(5) = {0, 0, -L/2, h};
Point(6) = {L/2, 0, -L/2, h};
Point(7) = {L/2, L, -L/2, h};
Point(8) = {0, L, -L/2, h};

Point(9) = {1.01*L/2,    L/2-L1/2, L1/2, h1};
Point(10) = {1.01*L/2+L1, L/2-L1/2, L1/2, h1};
Point(11) = {1.01*L/2+L1, L/2+L1/2, L1/2, h1};
Point(12) = {1.01*L/2,    L/2+L1/2, L1/2, h1};

Point(13) = {1.01*L/2,    L/2-L1/2, -L1/2, h1};
Point(14) = {1.01*L/2+L1, L/2-L1/2, -L1/2, h1};
Point(15) = {1.01*L/2+L1, L/2+L1/2, -L1/2, h1};
Point(16) = {1.01*L/2,    L/2+L1/2, -L1/2, h1};

// Lines
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

Line(5) = {8, 7};
Line(6) = {7, 6};
Line(7) = {6, 5};
Line(8) = {5, 8};

Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line(13) = {12, 9};
Line(14) = {9, 10};
Line(15) = {10, 11};
Line(16) = {11, 12};

Line(17) = {16, 15};
Line(18) = {15, 14};
Line(19) = {14, 13};
Line(20) = {13, 16};

Line(21) = {9, 13};
Line(22) = {10, 14};
Line(23) = {11, 15};
Line(24) = {12, 16};

// Geometric and Physical Surface
Line Loop(25) = {1, 2, 3, 4};
Plane Surface(25) = {25};
Line Loop(26) = {5, 6, 7, 8};
Plane Surface(26) = {26};
Line Loop(27) = {-5, -12, -4, 11};
Plane Surface(27) = {27};
Line Loop(28) = {-6, -11, -3, 10};
Plane Surface(28) = {28};
Line Loop(29) = {9, -7, -10, -2};
Plane Surface(29) = {29};
Line Loop(30) = {12, -8, -9, -1};
Plane Surface(30) = {30};

Line Loop(31) = {13, 14, 15, 16};
Plane Surface(31) = {31};
Line Loop(32) = {17, 18, 19, 20};
Plane Surface(32) = {32};
Line Loop(33) = {-17, -24, -16, 23};
Plane Surface(33) = {33};
Line Loop(34) = {-18, -23, -15, 22};
Plane Surface(34) = {34};
Line Loop(35) = {21, -19, -22, -14};
Plane Surface(35) = {35};
Line Loop(36) = {24, -20, -21, -13};
Plane Surface(36) = {36};

Surface Loop(37) = {25, 26, 27, 28, 29, 30};
Volume(37) = {37};

Surface Loop(38) = {31, 32, 33, 34, 35, 36};
Volume(38) = {38};

Transfinite Surface "*";
Recombine Surface "*";

Transfinite Volume "*";
Recombine Volume "*";
