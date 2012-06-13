h = 0.0015;
L = 0.015;

Point(1) = {-L, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, h, 0, h};
Point(4) = {-L, h, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Physical Surface(7) = {6};

Transfinite Surface "*";
