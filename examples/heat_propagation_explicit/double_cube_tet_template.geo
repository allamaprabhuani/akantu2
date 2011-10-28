h = 6.03863042319*3;
L = 24*4.047971240341177*2;
LZ = 50*4.047971240341177;
contact = SSSSSS;
L2 = contact/512*L;

p = 15;

Point(1) = {0.0, 0.0, LZ, p};
Point(2) = {L, 0.0, LZ, p};
Point(3) = {L, L, LZ, p};
Point(4) = {0.0, L, LZ, p};

Point(5) = {0.0, 0.0, 0.0, p};
Point(6) = {L, 0.0, 0.0, p};
Point(7) = {L, L, 0.0, p};
Point(8) = {0.0, L, 0.0, p};


Point(11) = {0.0, 0.0, -h, p};
Point(12) = {L, 0.0, -h, p};
Point(13) = {L, L, -h, p};
Point(14) = {0.0, L, -h, p};

Point(15) = {0.0, 0.0, -LZ, p};
Point(16) = {L, 0.0, -LZ, p};
Point(17) = {L, L, -LZ, p};
Point(18) = {0.0, L, -LZ, p};


Point(25) = {L/2-L2/2, L/2-L2/2, 0.0, p};
Point(26) = {L/2+L2/2, L/2-L2/2, 0.0, p};
Point(27) = {L/2+L2/2, L/2+L2/2, 0.0, p};
Point(28) = {L/2-L2/2, L/2+L2/2, 0.0, p};

Point(31) = {L/2-L2/2, L/2-L2/2, -h, p};
Point(32) = {L/2+L2/2, L/2-L2/2, -h, p};
Point(33) = {L/2+L2/2, L/2+L2/2, -h, p};
Point(34) = {L/2-L2/2, L/2+L2/2, -h, p};





Line(1) = {15, 18};
Line(2) = {18, 17};
Line(3) = {17, 16};
Line(4) = {16, 15};
Line(5) = {11, 15};
Line(6) = {12, 11};
Line(7) = {12, 16};
Line(8) = {12, 13};
Line(9) = {13, 17};
Line(10) = {13, 14};
Line(11) = {14, 18};
Line(12) = {14, 11};
Line(13) = {4, 3};
Line(14) = {3, 2};
Line(15) = {2, 1};
Line(16) = {1, 4};
Line(17) = {1, 5};
Line(18) = {5, 8};
Line(19) = {8, 7};
Line(20) = {4, 8};
Line(21) = {3, 7};
Line(22) = {7, 6};
Line(23) = {6, 2};
Line(24) = {6, 5};
Line(25) = {25, 26};
Line(26) = {26, 27};
Line(27) = {27, 28};
Line(28) = {28, 25};
Line(29) = {33, 34};
Line(30) = {34, 31};
Line(31) = {31, 32};
Line(32) = {32, 33};
Line(46) = {25, 31};
Line(47) = {28, 34};
Line(48) = {27, 33};
Line(49) = {26, 32};


//topcube
Line Loop(50) = {16, 20, -18, -17};
Plane Surface(51) = {50};
Line Loop(52) = {15, 17, -24, 23};
Plane Surface(53) = {52};
Line Loop(54) = {14, -23, -22, -21};
Plane Surface(55) = {54};
Line Loop(56) = {13, 21, -19, -20};
Plane Surface(57) = {56};
Line Loop(58) = {15, 16, 13, 14};
Plane Surface(59) = {58};

//Transfinite Surface {51};
//Transfinite Surface {53};
//Transfinite Surface {55};
//Transfinite Surface {57};
//Transfinite Surface {59};

//bottomcube
Line Loop(60) = {8, 9, 3, -7};
Plane Surface(61) = {60};
Line Loop(62) = {10, 11, 2, -9};
Plane Surface(63) = {62};
Line Loop(64) = {12, 5, 1, -11};
Plane Surface(65) = {64};
Line Loop(66) = {7, 4, -5, -6};
Plane Surface(67) = {66};
Line Loop(68) = {2, 3, 4, 1};
Plane Surface(69) = {68};

//Transfinite Surface {61};
//Transfinite Surface {63};
//Transfinite Surface {65};
//Transfinite Surface {67};
//Transfinite Surface {69};

//middle cube
Line Loop(70) = {26, 48, -32, -49};
Plane Surface(71) = {70};
Line Loop(72) = {25, 49, -31, -46};
Plane Surface(73) = {72};
Line Loop(74) = {28, 46, -30, -47};
Plane Surface(75) = {74};
Line Loop(76) = {27, 47, -29, -48};
Plane Surface(77) = {76};

//junction surfaces
Line Loop(78) = {24, 18, 19, 22};
Line Loop(79) = {25, 26, 27, 28};
Plane Surface(80) = {78, 79};
Line Loop(81) = {6, -12, -10, -8};
Line Loop(82) = {31, 32, 29, 30};
Plane Surface(83) = {81, 82};


Surface Loop(84) = {63, 83, 67, 61, 69, 65, 75, 80, 53, 59, 51, 57, 55, 77, 71, 73};
Volume(85) = {84};

Physical Volume(500) = {85};
//Coherence;
