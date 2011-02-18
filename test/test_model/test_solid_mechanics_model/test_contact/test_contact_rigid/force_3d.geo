h_contact = 0.2; 
h_sphere = 0.7;
h_plate = 0.4;	 

Point(1) = {0,  0, 0, h_sphere};
Point(2) = {0, -1, 0, h_contact};
Point(3) = {1, -1, 0, h_contact};
Point(4) = {1,  0, 0, h_sphere};

Point(5) = {0,  0, 1, h_sphere};
Point(6) = {0, -1, 1, h_contact};
Point(7) = {1, -1, 1, h_contact};
Point(8) = {1,  0, 1, h_sphere};

Point(9)  = {-0.25, -1,   -0.25, h_plate};
Point(10) = {-0.25, -1.2, -0.25, h_plate};
Point(11) = { 1.25, -1.2, -0.25, h_plate};
Point(12) = { 1.25, -1,   -0.25, h_plate};

Point(13) = {-0.25, -1,   1.25, h_plate};
Point(14) = {-0.25, -1.2, 1.25, h_plate};
Point(15) = { 1.25, -1.2, 1.25, h_plate};
Point(16) = { 1.25, -1,   1.25, h_plate};

Line(1) = {1, 5};
Line(2) = {5, 8};
Line(3) = {8, 4};
Line(4) = {4, 1};
Line(5) = {6, 2};
Line(6) = {2, 3};
Line(7) = {3, 7};
Line(8) = {7, 6};
Line(9) = {5, 6};
Line(10) = {1, 2};
Line(11) = {4, 3};
Line(12) = {8, 7};

Line(13) = {9, 13};
Line(14) = {13, 16};
Line(15) = {16, 12};
Line(16) = {12, 9};
Line(17) = {14, 10};
Line(18) = {10, 11};
Line(19) = {11, 15};
Line(20) = {15, 14};
Line(21) = {13, 14};
Line(22) = {9, 10};
Line(23) = {12, 11};
Line(24) = {16, 15};

Line Loop(26) = {4, 1, 2, 3};
Plane Surface(26) = {26};
Line Loop(28) = {9, -8, -12, -2};
Plane Surface(28) = {28};
Line Loop(30) = {-7, -11, -3, 12};
Plane Surface(30) = {30};
Line Loop(32) = {-10, -4, 11, -6};
Plane Surface(32) = {32};
Line Loop(34) = {10, -5, -9, -1};
Plane Surface(34) = {34};
Line Loop(36) = {8, 5, 6, 7};
Plane Surface(36) = {36};

Line Loop(38) = {16, 13, 14, 15};
Plane Surface(38) = {38};
Line Loop(40) = {21, -20, -24, -14};
Plane Surface(40) = {40};
Line Loop(42) = {-19, -23, -15, 24};
Plane Surface(42) = {42};
Line Loop(44) = {-22, -16, 23, -18};
Plane Surface(44) = {44};
Line Loop(46) = {22, -17, -21, -13};
Plane Surface(46) = {46};
Line Loop(48) = {20, 17, 18, 19};
Plane Surface(48) = {48};

Surface Loop(50) = {26, 28, 30, 32, 34, 36};
Volume(50) = {50};

Surface Loop(52) = {38, 40, 42, 44, 46, 48};
Volume(52) = {52};
