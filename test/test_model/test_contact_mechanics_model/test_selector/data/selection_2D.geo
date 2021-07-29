h = 0.4;

y = 1;
x = 2;

Point(1) = { x/2., y/2, 0, h};
Point(2) = {-x/2., y/2, 0, h};
Point(3) = {-x/2.,-y/2, 0, h};
Point(4) = { x/2.,-y/2, 0, h};
Point(5) = {-x/2.,   0, 0, h};
Point(6) = { x/2.,   0, 0, h};
Point(7) = {-x/2.,   0, 0, h};
Point(8) = { x/2.,   0, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 6};
Line(4) = {6, 1};

Line(5) = {7, 3};
Line(6) = {3, 4};
Line(7) = {4, 8};
Line(8) = {8, 7};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Line("fixed") = {6};
Physical Line("loading") = {1};
Physical Line("master") = {3};
Physical Line("slave") = {8};

Physical Line("sides") = {2, 5, 7, 4};

Physical Surface("master_body") = {1};
Physical Surface("slave_body") = {2};
