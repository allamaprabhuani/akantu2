h_contact = 0.00625;
h_topcube = 0.00625;
h_plate = 0.8;	 

Point(1) = {0, 0, 0, h_topcube};
Point(2) = {0, -1, 0, h_contact};
Point(3) = {1, -1, 0, h_contact};
Point(4) = {1, 0, 0, h_topcube};

Point(5) = {-0.25, -1, 0, h_plate};
Point(6) = {-0.25, -1.5, 0, h_plate};
Point(7) = { 2, -1.5, 0, h_plate};
Point(8) = { 2, -1, 0, h_plate};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line Loop(9) = {4, 1, 2, 3};
Plane Surface(9) = {9};

Line Loop(11) = {8, 5, 6, 7};
Plane Surface(11) = {11};

Transfinite Surface "*";
Recombine Surface "*";
