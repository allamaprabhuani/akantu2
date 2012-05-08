h = 1.0;

Point(1) = {0, 0, 0, h};
Point(2) = {0, 1, 0, h};
Point(3) = {-1, 1, 0, h};
Point(4) = {-1, 0, 0, h};
Point(5) = {0, 0, 0, h};
Point(6) = {1, 0, 0, h};
Point(7) = {1, 1, 0, h};
Point(8) = {0, 1, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line Loop(9) = {5, 6, 7, 8};
Plane Surface(10) = {9};
Line Loop(11) = {1, 2, 3, 4};
Plane Surface(12) = {11};
Physical Surface(0) = {12};
Physical Surface(1) = {10};
