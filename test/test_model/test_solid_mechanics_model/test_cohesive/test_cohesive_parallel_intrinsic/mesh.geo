h = 1.;

Point(1) = {0, 0, 0, h};
Point(2) = {0, 1, 0, h};
Point(3) = {-1, 1, 0, h};
Point(4) = {-1, 0, 0, h};
Point(6) = {1, 0, 0, h};
Point(7) = {1, 1, 0, h};
Point(8) = {1, -1, 0, h};
Point(9) = {0, -1, 0, h};
Point(10) = {-1, -1, 0, h};
Line(1) = {10, 9};
Line(2) = {9, 8};
Line(3) = {8, 6};
Line(4) = {6, 7};
Line(5) = {7, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 10};
Line Loop(9) = {3, 4, 5, 6, 7, 8, 1, 2};
Plane Surface(10) = {9};
Physical Surface(11) = {10};

Transfinite Line {1, 2, 3, 4, 5, 6, 7, 8} = 2;