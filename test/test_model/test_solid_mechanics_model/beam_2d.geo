h = 0.5;

Point(1) = {0,  0, 0, h};
Point(2) = {10, 0, 0, h};
Point(3) = {10, 1, 0, h};
Point(4) = {5,  1, 0, h};
Point(5) = {0,  1, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};
Line Loop(6) = {1, 2, 3, 4, 5};
Plane Surface(7) = {6};
Physical Surface(1) = {7};
