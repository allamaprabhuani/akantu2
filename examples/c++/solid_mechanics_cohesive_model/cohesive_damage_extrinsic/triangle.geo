h = .1;
L = 0.4;

Point(1) = { L, L, 0, h};
Point(2) = {0., L, 0, h};
Point(3) = {0.,0., 0, h};
Point(4) = {L,0., 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Physical Surface(7) = {6};

Transfinite Line {2, 4} = 2;
Transfinite Line {1, 3} = 3;
Transfinite Surface {6} = {1, 2, 3, 4};


Physical Curve("left", 8) = {2};
Physical Curve("right", 9) = {4};
Physical Point("point", 10) = {3};
