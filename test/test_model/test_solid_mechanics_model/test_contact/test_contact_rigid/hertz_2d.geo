h_contact = 0.005;
h_sphere = 0.1;
h_plate = 0.4;	 

Point(1) = {0,  0, 0, h_sphere};
Point(2) = {1,  0, 0, h_sphere};
Point(3) = {0, -1, 0, h_contact};

Point(4) = {-0.25, -1,   0, h_plate};
Point(5) = {-0.25, -1.2, 0, h_plate};
Point(6) = { 1.25, -1.2, 0, h_plate};
Point(7) = { 1.25, -1,   0, h_plate};

Circle(1) = {3, 1, 2};
Line(2) = {2, 1};
Line(3) = {1, 3};

Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 4};

Line Loop(9) = {1, 2, 3};
Plane Surface(9) = {9};

Line Loop(11) = {7, 4, 5, 6};
Plane Surface(11) = {11};
