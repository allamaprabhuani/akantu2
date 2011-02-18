h_contact = 0.05;
h_sphere = 0.2;
h_plate = 0.05;	 

Point(1) = {0, 0, 0, h_sphere};
Point(2) = {0, -1, 0, h_contact};
Point(3) = {1, -1, 0, h_contact};
Point(4) = {1, 0, 0, h_sphere};

Point(5) = {-0.25, -1, 0, h_plate};
Point(6) = {-0.25, -1.2, 0, h_plate};
Point(7) = { 1.25, -1.2, 0, h_plate};
Point(8) = { 1.25, -1, 0, h_plate};

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
