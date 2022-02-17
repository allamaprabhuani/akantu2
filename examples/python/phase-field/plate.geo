h1 = 1e-4;
h2 = 1e-3;
L = 32e-3;
H = 16e-3;
l = 4e-3;
Point(1) = {0, 0, 0, h1};
Point(2) = {L, 0, 0, h1};
Point(3) = {L, H/2, 0, h2};
Point(4) = {0, H/2, 0, h2};
Point(5) = {l, 0, 0, h1};

Point(6) =  {0, 0, 0, h1};
Point(7) =  {L, -H/2, 0, h2};
Point(8) =  {0, -H/2, 0, h2};


Line(1) = {1, 5};
Line(2) = {4, 1};
Line(3) = {3, 4};
Line(4) = {2, 3};
Line(5) = {5, 2};

Line Loop(1) = {2, 3, 4, 5, 1};
Plane Surface(1) = {1};

Line(6) =  {5, 6};
Line(7) =  {6, 8};
Line(8) =  {8, 7};
Line(9) =  {7, 2};
Line Loop(2) = {6, 7, 8, 9, -5};
Plane Surface(2) = {2};


Physical Surface(8) = {1,2};
Physical Line("left") = {2,7};
Physical Line("bottom") = {8};
Physical Line("top") = {3};
Physical Line("right") = {4,9};
