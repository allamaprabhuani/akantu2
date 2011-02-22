h_contact = 0.005; 
h_sphere = 0.3; 
h_plate = 0.4;	 

Point(1) = {0,  0, 0, h_sphere};
Point(2) = {1,  0, 0, h_sphere};
Point(3) = {0, -1, 0, h_contact};
Point(4) = {0,  0, 1, h_sphere};

Point(5) = {-0.5, -1, -0.5, h_plate};
Point(6) = {-0.5, -1.2, -0.5, h_plate};
Point(7) = { 0.5, -1.2, -0.5, h_plate};
Point(8) = { 0.5, -1, -0.5, h_plate};

Point(9) = { -0.5, -1, 0.5, h_plate};
Point(10) = {-0.5, -1.2, 0.5, h_plate};
Point(11) = { 0.5, -1.2, 0.5, h_plate};
Point(12) = { 0.5, -1, 0.5, h_plate};

Point(13) = {-1, 0, 0, h_sphere};
Point(14) = {0, 0, -1, h_sphere};



Circle(1) = {2, 1, 14};
Circle(2) = {14, 1, 13};
Circle(3) = {4, 1, 2};
Circle(4) = {3, 1, 14};
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

Circle(19) = {13, 1, 4};
Circle(20) = {3, 1, 13};

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
Line Loop(32) = {14, 11, 12, 13};
Plane Surface(32) = {32};

Line Loop(34) = {-5, 6, -3};
Ruled Surface(34) = {34};
Line Loop(36) = {4, -1, -6};
Ruled Surface(36) = {36};
Line Loop(38) = {20, -2, -4};
Ruled Surface(38) = {38};
Line Loop(40) = {5, -19, -20};
Ruled Surface(40) = {40};
Line Loop(42) = {3, 1, 2, 19};
Plane Surface(42) = {42};


Surface Loop(50) = {34, 36, 38, 40, 42};
Volume(50) = {50};

Surface Loop(52) = {22, 28, 30, 24, 26, 32};
Volume(52) = {52};
