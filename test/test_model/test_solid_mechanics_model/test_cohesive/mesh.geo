h =0.2;

Point(1) = {0, 0, 0, h};
Point(2) = {0, 1, 0, h};
Point(3) = {-1, 1, 0, h};
Point(4) = {-1, 0, 0, h};
Point(6) = {1, 0, 0, h};
Point(7) = {1, 1, 0, h};
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 6};
Line(5) = {6, 7};
Line(6) = {7, 2};

Line Loop(7) = {1, 2, 3, 4, 5, 6};
Plane Surface(8) = {7};
Physical Surface(0) = {8};
