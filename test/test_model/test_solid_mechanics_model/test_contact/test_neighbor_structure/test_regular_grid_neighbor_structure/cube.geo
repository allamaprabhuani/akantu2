h  = 0.05;

Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};
Point(5) = {0, 0, 0.2, h};
Point(6) = {1, 0, 0.2, h};
Point(7) = {1, 1, 0.2, h};
Point(8) = {0, 1, 0.2, h};

Line(1) = {5, 6};
Line(2) = {6, 7};
Line(3) = {7, 8};
Line(4) = {8, 5};
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};
Line(9) = {5, 1};
Line(10) = {6, 2};
Line(11) = {7, 3};
Line(12) = {8, 4};

Line Loop(25) = {4, 1, 2, 3};
Plane Surface(26) = {25};
Line Loop(27) = {7, 8, 5, 6};
Plane Surface(28) = {27};
Line Loop(29) = {12, -7, -11, 3};
Plane Surface(30) = {29};
Line Loop(31) = {12, 8, -9, -4};
Plane Surface(32) = {31};
Line Loop(33) = {1, 10, -5, -9};
Plane Surface(34) = {33};
Line Loop(35) = {2, 11, -6, -10};
Plane Surface(36) = {35};

Surface Loop(49) = {32, 30, 28, 34, 26, 36};
Volume(50) = {49};
