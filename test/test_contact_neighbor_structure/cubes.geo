h  = 0.05;
h1 = 0.2;

Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};
Point(5) = {0, 0, 1, h};
Point(6) = {1, 0, 1, h};
Point(7) = {1, 1, 1, h};
Point(8) = {0, 1, 1, h};

Point(9)  = {0.45, 0.45, 1.001, h1};
Point(10) = {0.55, 0.45, 1.001, h1};
Point(11) = {0.55, 0.55, 1.001, h1};
Point(12) = {0.45, 0.55, 1.001, h1};
Point(13) = {0.45, 0.45, 1.101, h1};
Point(14) = {0.55, 0.45, 1.101, h1};
Point(15) = {0.55, 0.55, 1.101, h1};
Point(16) = {0.45, 0.55, 1.101, h1};

Point(17) = {0.90, 0.45, 1.001, h1};
Point(18) = {1.00, 0.45, 1.001, h1};
Point(19) = {1.00, 0.55, 1.001, h1};
Point(20) = {0.90, 0.55, 1.001, h1};
Point(21) = {0.90, 0.45, 1.101, h1};
Point(22) = {1.00, 0.45, 1.101, h1};
Point(23) = {1.00, 0.55, 1.101, h1};
Point(24) = {0.90, 0.55, 1.101, h1};


Line(1) = {5, 6};
Line(2) = {6, 7};
Line(3) = {7, 8};
Line(4) = {8, 5};
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};
Line(9) = {5, 1};
Line(10) = {6, 2};
Line(11) = {7, 3};
Line(12) = {8, 4};

Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 13};
Line(17) = {9, 10};
Line(18) = {10, 11};
Line(19) = {11, 12};
Line(20) = {12, 9};
Line(21) = {13, 9};
Line(22) = {14, 10};
Line(23) = {15, 11};
Line(24) = {16, 12};

Line(25) = {21, 22};
Line(26) = {22, 23};
Line(27) = {23, 24};
Line(28) = {24, 21};
Line(29) = {17, 18};
Line(30) = {18, 19};
Line(31) = {19, 20};
Line(32) = {20, 17};
Line(33) = {21, 17};
Line(34) = {22, 18};
Line(35) = {23, 19};
Line(36) = {24, 20};


Line Loop(25) = {4, 1, 2, 3};
Plane Surface(26) = {25};
Line Loop(27) = {7, 8, 5, 6};
Plane Surface(28) = {27};
Line Loop(29) = {12, -7, -11, 3};
Plane Surface(30) = {29};
Line Loop(31) = {12, 8, -9, -4};
Plane Surface(32) = {31};
Line Loop(33) = {1, 10, -5, -9};
Plane Surface(34) = {33};
Line Loop(35) = {2, 11, -6, -10};
Plane Surface(36) = {35};

Line Loop(37) = {16, 13, 14, 15};
Plane Surface(38) = {37};
Line Loop(39) = {20, 17, 18, 19};
Plane Surface(40) = {39};
Line Loop(41) = {24, 20, -21, -16};
Plane Surface(42) = {41};
Line Loop(43) = {13, 22, -17, -21};
Plane Surface(44) = {43};
Line Loop(45) = {14, 23, -18, -22};
Plane Surface(46) = {45};
Line Loop(47) = {15, 24, -19, -23};
Plane Surface(48) = {47};

Line Loop(49) = {28, 25, 26, 27};
Plane Surface(50) = {49};
Line Loop(51) = {32, 29, 30, 31};
Plane Surface(52) = {51};
Line Loop(53) = {36, 32, -33, -28};
Plane Surface(54) = {53};
Line Loop(55) = {25, 34, -29, -33};
Plane Surface(56) = {55};
Line Loop(57) = {26, 35, -30, -34};
Plane Surface(58) = {57};
Line Loop(59) = {27, 36, -31, -35};
Plane Surface(60) = {59};


Surface Loop(49) = {32, 30, 28, 34, 26, 36};
Volume(50) = {49};

Surface Loop(51) = {42, 48, 38, 44, 46, 40};
Volume(52) = {51};

Surface Loop(53) = {54, 60, 50, 56, 58, 52};
Volume(54) = {53};
