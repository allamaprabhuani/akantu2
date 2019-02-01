cl1 = 0.01;
cl2 = 0.1;
cl3 = 0.1;
cl4 = 0.005;
radius = 0.25;
depth = -0.25;
y = 0.1;
Dz = 1;
Point(1) = {0, y, 0, cl1};
Point(2) = {radius, radius + y, 0, cl2};
Point(3) = {-radius, radius + y, 0, cl2};
Point(4) = {0, radius + y, 0, cl2};
Point(5) = {0, radius + y, -radius, cl2};
Point(6) = {0, radius + y, radius, cl2};

/*Point(7) = {0, y, 0, cl4};
Point(8) = {0.25, y, 0, cl4};
Point(9) = {-0.25, y, 0, cl4};
Point(10) = {0, y, -0.25, cl4};
Point(11) = {0, y, 0.25, cl4};*/

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

/*Circle(9) = {9, 7, 11};
Circle(10) = {11, 7, 8};
Circle(11) = {8, 7, 10};
Circle(12) = {10, 7, 9};
Line Loop(6) = {10, 11, 12, 9};
Plane Surface(6) = {6};
Extrude {0, -cl3, 0} {
  Surface{6};
}*/

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
Physical Surface("contact_surface") = {9};
Physical Surface("bottom_surface") = {8};
Physical Volume("top_body") = {1};
Physical Volume("bot_body") = {3};
Physical Surface("rigid") = {1, 2, 4, 3};
