h = 1;

Point(1) = {0.0, 0.0, -1.0, h};
Point(2) = {1.0, 0.0, -1.0, h};
Point(3) = {0.0, 1.0, -1.0, h};
Point(4) = {1.0, 1.0, -1.0, h};
Point(5) = {-1.0,0.0, -1.0, h};
Point(6) = {0.0,-1.0, -1.0, h};
Point(7) = {-1.0,-1.0, -1.0, h};
Point(8) = {1.0, -1.0, -1.0, h};
Point(9) = {-1.0, 1.0, -1.0, h};

Point(11) = {0.0, 0.0, 0.0, h};
Point(12) = {1.0, 0.0, 0.0, h};
Point(13) = {0.0, 1.0, 0.0, h};
Point(14) = {1.0, 1.0, 0.0, h};
Point(15) = {-1.0,0.0, 0.0, h};
Point(16) = {0.0,-1.0, 0.0, h};
Point(17) = {-1.0,-1.0, 0.0, h};
Point(18) = {1.0, -1.0, 0.0, h};
Point(19) = {-1.0, 1.0, 0.0, h};

Point(21) = {0.0, 0.0, 1.0, h};
Point(22) = {1.0, 0.0, 1.0, h};
Point(23) = {0.0, 1.0, 1.0, h};
Point(24) = {1.0, 1.0, 1.0, h};
Point(25) = {-1.0,0.0, 1.0, h};
Point(26) = {0.0,-1.0, 1.0, h};
Point(27) = {-1.0,-1.0, 1.0, h};
Point(28) = {1.0, -1.0, 1.0, h};
Point(29) = {-1.0, 1.0, 1.0, h};

Line(1) = {13, 23};
Line(2) = {23, 24};
Line(3) = {24, 14};
Line(4) = {14, 13};
Line(5) = {23, 29};
Line(6) = {29, 19};
Line(7) = {19, 13};
Line(8) = {14, 13};
Line(9) = {13, 3};
Line(10) = {3, 4};
Line(11) = {4, 14};
Line(12) = {19, 9};
Line(13) = {9, 3};

Line(14) = {29, 25};
Line(15) = {25, 21};
Line(16) = {21, 23};
Line(17) = {21, 22};
Line(18) = {22, 24};
Line(19) = {25, 27};
Line(20) = {27, 26};
Line(21) = {26, 21};
Line(22) = {26, 28};
Line(23) = {28, 22};
Line(24) = {22, 12};
Line(25) = {12, 14};
Line(26) = {12, 2};
Line(27) = {2, 4};
Line(28) = {28, 18};
Line(29) = {18, 12};
Line(30) = {18, 8};
Line(31) = {8, 2};
Line(32) = {2, 1};
Line(33) = {1, 3};
Line(34) = {8, 6};
Line(35) = {6, 1};
Line(36) = {1, 5};
Line(37) = {5, 9};
Line(38) = {5, 7};
Line(39) = {7, 6};
Line(40) = {5, 15};
Line(41) = {15, 19};
Line(42) = {15, 25};
Line(43) = {7, 17};
Line(44) = {17, 15};
Line(45) = {17, 27};
Line(46) = {6, 16};
Line(47) = {16, 17};
Line(48) = {16, 26};
Line(49) = {16, 18};
Line Loop(50) = {5, 14, 15, 16};
Plane Surface(51) = {50};
Line Loop(52) = {2, -18, -17, 16};
Plane Surface(53) = {52};
Line Loop(54) = {15, -21, -20, -19};
Plane Surface(55) = {54};
Line Loop(56) = {17, -23, -22, 21};
Plane Surface(57) = {56};
Line Loop(58) = {3, -25, -24, 18};
Plane Surface(59) = {58};
Line Loop(60) = {11, -25, 26, 27};
Plane Surface(61) = {60};
Line Loop(62) = {24, -29, -28, 23};
Plane Surface(63) = {62};
Line Loop(64) = {26, -31, -30, 29};
Plane Surface(65) = {64};
Line Loop(66) = {33, 10, -27, 32};
Plane Surface(67) = {66};
Line Loop(68) = {32, -35, -34, 31};
Plane Surface(69) = {68};
Line Loop(70) = {36, 38, 39, 35};
Plane Surface(71) = {70};
Line Loop(72) = {37, 13, -33, 36};
Plane Surface(73) = {72};
Line Loop(74) = {14, -42, 41, -6};
Plane Surface(75) = {74};
Line Loop(76) = {19, -45, 44, 42};
Plane Surface(77) = {76};
Line Loop(78) = {41, 12, -37, 40};
Plane Surface(79) = {78};
Line Loop(80) = {44, -40, 38, 43};
Plane Surface(81) = {80};
Line Loop(82) = {2, 3, 4, 1};
Plane Surface(83) = {82};
Line Loop(84) = {4, 9, 10, 11};
Plane Surface(85) = {84};
Line Loop(86) = {5, 6, 7, 1};
Plane Surface(87) = {86};
Line Loop(88) = {7, 9, -13, -12};
Plane Surface(89) = {88};
Line Loop(90) = {20, -48, 47, 45};
Plane Surface(91) = {90};
Line Loop(92) = {47, -43, 39, 46};
Plane Surface(93) = {92};
Line Loop(94) = {22, 28, -49, 48};
Plane Surface(95) = {94};
Line Loop(96) = {49, 30, 34, 46};
Plane Surface(97) = {96};
Surface Loop(98) = {77, 55, 51, 87, 75, 79, 89, 85, 83, 53, 59, 61, 65, 69, 67, 73, 71, 81, 93, 91, 95, 57, 63, 97};
Volume(99) = {98};
