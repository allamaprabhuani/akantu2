cl1 = 0.02;
cl2 = 0.02;
cl3 = 0.05;
cl4 = 0.02;
radius = 0.1;
depth = -0.1;
y = 0.0;
Dz = 1;

Point(1) = {0, 0, 0, cl1};
Point(2) = {radius/2, radius/2, 0, cl2};
Point(3) = {-radius/2, radius/2, 0, cl2};
Point(4) = {0, radius/2, 0, cl2};
Point(5) = {0, radius/2, -radius/2, cl2};
Point(6) = {0, radius/2, radius/2, cl2};

Point(12) = {0, y, 0, cl4};
Point(13) = {radius, y, radius, cl4};
Point(14) = {-radius, y, radius, cl4};
Point(15) = {-radius, y, -radius, cl4};
Point(16) = {radius, y, -radius, cl4};

Point(17) = {radius, y + depth, radius, cl2};
Point(18) = {-radius, y + depth, radius, cl2};
Point(19) = {-radius, y + depth, -radius, cl2};
Point(20) = {radius, y + depth, -radius, cl2};


Circle(1) = {3, 4, 1};
Circle(2) = {1, 4, 2};
Circle(3) = {6, 4, 1};
Circle(4) = {1, 4, 5};
Circle(5) = {3, 4, 6};
Circle(6) = {6, 4, 2};
Circle(7) = {2, 4, 5};
Circle(8) = {5, 4, 3};

Line Loop(1) = {3, 2, -6};
Ruled Surface(1) = {1};
Line Loop(2) = {-2, -7, 4};
Ruled Surface(2) = {2};
Line Loop(3) = {-4, -8, -1};
Ruled Surface(3) = {3};
Line Loop(4) = {1, -3, -5};
Ruled Surface(4) = {4};
Line Loop(5) = {6, 7, 8, 5};
Plane Surface(5) = {5};
Surface Loop(1) = {3, 2, 1, 4, 5};
Volume(1) = {1};

Line(95) = {17, 18};
Line(96) = {18, 19};
Line(97) = {19, 20};
Line(98) = {20, 17};
Line(100) = {13, 14};
Line(101) = {13, 17};
Line(102) = {18, 14};
Line(103) = {14, 15};
Line(104) = {19, 15};
Line(105) = {15, 16};
Line(106) = {20, 16};
Line(107) = {16, 13};

Line Loop(7) = {95, 96, 97, 98};
Plane Surface(8) = {7};
Line Loop(8) = {100, 103, 105, 107};
Plane Surface(9) = {8};
Line Loop(9) = {-106, 98, -101, -107};
Plane Surface(10) = {9};
Line Loop(10) = {101, 95, 102, -100};
Plane Surface(11) = {10};
Line Loop(11) = {-102, 96, 104, -103};
Plane Surface(12) = {11};
Line Loop(12) = {-104, 97, 106, -105};
Plane Surface(13) = {12};

Surface Loop(101) = {8, 9, 10, 11, 12, 13};
Volume(3) = {101};

Physical Surface("top_surface") = {5};
Physical Surface("flat") = {9};
Physical Surface("bottom_surface") = {8};
Physical Volume("top_body") = {1};
Physical Volume("bot_body") = {3};
Physical Surface("curved") = {1, 2, 4, 3};
