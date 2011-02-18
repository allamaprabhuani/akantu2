h_contact = 0.05; 
h_sphere = 0.2;
h_plate = 0.4;	 

Point(1) = {0,  0, 0, h_sphere};
Point(2) = {1,  0, 0, h_sphere};
Point(3) = {0, -1, 0, h_contact};
Point(4) = {0,  0, 1, h_sphere};

Point(5) = {-0.25, -1,   -0.2, h_plate};
Point(6) = {-0.25, -1.2, -0.2, h_plate};
Point(7) = { 1.25, -1.2, -0.2, h_plate};
Point(8) = { 1.25, -1,   -0.2, h_plate};

Point(9) = { -0.25, -1,   1.2, h_plate};
Point(10) = {-0.25, -1.2, 1.2, h_plate};
Point(11) = { 1.25, -1.2, 1.2, h_plate};
Point(12) = { 1.25, -1,   1.2, h_plate};

Line(1) = {2, 1};
Line(2) = {1, 4};
Circle(3) = {4, 1, 2};
Line(4) = {1, 3};
Circle(5) = {3, 1, 4};
Circle(6) = {3, 1, 2};

Line(7) = {5, 9};
Line(8) = {9, 12};
Line(9) = {12, 8};
Line(10) = {8, 5};

Line(11) = {10, 6};
Line(12) = {6, 7};
Line(13) = {7, 11};
Line(14) = {11, 10};
Line(15) = {9, 10};
Line(16) = {5, 6};
Line(17) = {8, 7};
Line(18) = {12, 11};

Line Loop(20) = {6, -3, -5};
Ruled Surface(20) = {20};

Line Loop(22) = {10, 7, 8, 9};
Plane Surface(22) = {22};
Line Loop(24) = {15, -14, -18, -8};
Plane Surface(24) = {24};
Line Loop(26) = {-13, -17, -9, 18};
Plane Surface(26) = {26};
Line Loop(28) = {-16, -10, 17, -12};
Plane Surface(28) = {28};
Line Loop(30) = {16, -11, -15, -7};
Plane Surface(30) = {30};

Line Loop(32) = {5, -2, 4};
Plane Surface(32) = {32};
Line Loop(34) = {-4, -1, -6};
Plane Surface(34) = {34};
Line Loop(36) = {1, 2, 3};
Plane Surface(36) = {36};

Line Loop(40) = {14, 11, 12, 13};
Plane Surface(40) = {40};

Surface Loop(38) = {36, 34, 32, 20};
Volume(38) = {38};

Surface Loop(42) = {22, 28, 30, 24, 26, 40};
Volume(42) = {42};
