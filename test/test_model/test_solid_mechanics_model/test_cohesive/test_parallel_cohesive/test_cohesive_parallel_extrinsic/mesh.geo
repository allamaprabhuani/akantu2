h = 1.;

Point(1) = { 1, 1, 0, h};
Point(2) = {-1, 1, 0, h};
Point(3) = {-1,-1, 0, h};
Point(4) = { 1,-1, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Physical Surface(7) = {6};

Transfinite Line {1, 3} = 3;
Transfinite Line {2, 4} = 3;
// Transfinite Surface "*";
