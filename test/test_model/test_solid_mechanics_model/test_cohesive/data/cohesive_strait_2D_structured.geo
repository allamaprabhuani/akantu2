h = 0.4;

y = 1;
x = 2;

Point(1) = { x/2., y/2, 0, h};
Point(2) = {-x/2., y/2, 0, h};
Point(3) = {-x/2.,-y/2, 0, h};
Point(4) = { x/2.,-y/2, 0, h};
Point(5) = {-x/2.,   0, 0, h};
Point(6) = { x/2.,   0, 0, h};

Point(7) = {    0,   0, 0, h};
Point(8) = {    0, y/2, 0, h};
Point(9) = {    0,-y/2, 0, h};


Line(11) = {1, 8};
Line(12) = {8, 2};
Line(2) = {2, 5};
Line(31) = {5, 7};
Line(32) = {7, 6};
Line(4) = {6, 1};

Line(5) = {5, 3};
Line(61) = {3, 9};
Line(62) = {9, 4};
Line(7) = {4, 6};

Line(8) = {7, 8};
Line(9) = {7, 9};

Line Loop(11) = {11, -8, 32, 4};
Line Loop(12) = {12, 2, 31, 8};
Line Loop(21) = {-32, 9, 62, 7};
Line Loop(22) = {-31, 5, 61, -9};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(21) = {21};
Plane Surface(22) = {22};

Physical Line("fixed") = {61, 62};
Physical Line("loading") = {11, 12};
Physical Line("insertion") = {31, 32};
Physical Line("sides") = {2, 5, 7, 4};
Physical Line("middle") = {8, 9};

Physical Surface("body") = {11, 12, 21, 22};

Recombine Surface "*";
Transfinite Surface "*";
Mesh.SecondOrderIncomplete = 1;
