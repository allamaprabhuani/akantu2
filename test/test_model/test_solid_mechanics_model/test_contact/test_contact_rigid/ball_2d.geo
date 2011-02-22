h_ball = 0.05;
h_wall = 0.2;

a = 0.5;
da = 0.1;
b = 0.7;
db =0.1;

c = 0.2;
d = 0.1;
r = 0.05;

Point(1) = {0, 0, 0, h_wall};
Point(2) = {a, 0, 0, h_wall};
Point(3) = {a, b, 0, h_wall};
Point(4) = {0, b, 0, h_wall};

Point(5)  = {0,     -db, 0, h_wall};
Point(6)  = {a,     -db, 0, h_wall};
Point(7)  = {a+da,    0, 0, h_wall};
Point(8)  = {a+da,    b, 0, h_wall};
Point(9)  = {a,    b+db, 0, h_wall};
Point(10) = {0,    b+db, 0, h_wall};
Point(11) = {-da,     b, 0, h_wall};
Point(12) = {-da,     0, 0, h_wall};

Point(13) = {c, d, 0, h_ball};
Point(14) = {c+r, d, 0, h_ball};
Point(15) = {c-r, d, 0, h_ball};

Point(16) = {0, 0, 0, h_wall};
Point(17) = {a, 0, 0, h_wall};
Point(18) = {a, b, 0, h_wall};
Point(19) = {0, b, 0, h_wall};


Circle(1) = {14, 13, 15};
Circle(2) = {15, 13, 14};

Line(3) = {2, 1};
Line(4) = {1, 5};
Line(5) = {5, 6};
Line(6) = {6, 2};

Line(7)  = {18, 17};
Line(8)  = {17, 7};
Line(9)  = {7, 8};
Line(10) = {8, 18};

Line(11) = {4, 3};
Line(12) = {3, 9};
Line(13) = {9, 10};
Line(14) = {10, 4};

Line(15) = {16, 19};
Line(16) = {19, 11};
Line(17) = {11, 12};
Line(18) = {12, 16};

Line Loop(20) = {3, 4, 5, 6};
Plane Surface(20) = {20};

Line Loop(22) = {7, 8, 9, 10};
Plane Surface(22) = {22};

Line Loop(24) = {11, 12, 13, 14};
Plane Surface(24) = {24};

Line Loop(26) = {15, 16, 17, 18};
Plane Surface(26) = {26};

Line Loop(28) = {1, 2};
Plane Surface(28) = {28};
