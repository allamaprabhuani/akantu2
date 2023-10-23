h = 0.4;

z = 1;
y = 2;
x = 2;

Point(1) = { x/2, y/2, -z/2, h};
Point(2) = {-x/2, y/2, -z/2, h};
Point(3) = {-x/2,-y/2, -z/2, h};
Point(4) = { x/2,-y/2, -z/2, h};

Point(5) = {-x/2, 0,  -z/2, h};
Point(6) = { x/2, 0,  -z/2, h};

Point(7) = {    0,   0, -z/2, h};
Point(8) = {    0, y/2, -z/2, h};
Point(9) = {    0,-y/2, -z/2, h};


Line(11) = {1, 8};
Line(12) = {8, 2};
Line(2)  = {2, 5};
Line(31) = {5, 7};
Line(32) = {7, 6};
Line(4)  = {6, 1};

Line(5)  = {5, 3};
Line(61) = {3, 9};
Line(62) = {9, 4};
Line(7)  = {4, 6};

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

Extrude {0, 0, z} {
  Surface{11}; Surface{12}; Surface{21}; Surface{22};
}

Physical Line("fixed") = {145, 123};
Physical Line("loading") = {93, 71};
Physical Line("insertion") = {101, 79};
Physical Line("sides") = {97, 141, 83, 127};
Physical Line("middle") = {75, 119};

Physical Surface("body") = {1, 2, 3, 4};

Transfinite Surface "*";
Transfinite Volume "*";

Recombine Surface "*";

Mesh.SecondOrderIncomplete = 1;
