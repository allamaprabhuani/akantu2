h  = 0.05;
h1 = 0.4;

Point(1) = {0, 0, 0.7, h};
Point(2) = {1, 0, 0.7, h};
Point(3) = {1, 1, 0.7, h};
Point(4) = {0, 1, 0.7, h};
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

Surface Loop(49) = {32, 30, 28, 34, 26, 36};
Volume(50) = {49};

Surface Loop(51) = {42, 48, 38, 44, 46, 40};
Volume(52) = {51};

Physical Volume(0) = {50};
Physical Volume(1) = {52};
