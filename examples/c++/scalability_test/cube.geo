h = .002;

Point(1) = {1, 1, -1, h};
Point(2) = {-1, 1, -1, h};
Point(3) = {-1, -1, -1, h};
Point(4) = {1, -1, -1, h};

Point(5) = {0, -1, -1, h};
Point(6) = {0, 1, -1, h};

Line(1) = {1, 6};
Line(2) = {6, 2};
Line(3) = {2, 3};
Line(4) = {3, 5};
Line(5) = {5, 4};
Line(6) = {4, 1};
Line(7) = {5, 6};
Line Loop(5) = {1, -7, 5, 6};
Line Loop(6) = {2, 3, 4, 7};

Plane Surface(5) = {-5};
Plane Surface(6) = {-6};

Point(11) = {1, 1, 1, h};
Point(12) = {-1, 1, 1, h};
Point(13) = {-1, -1, 1, h};
Point(14) = {1, -1, 1, h};

Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 11};
Line Loop(15) = {11, 12, 13, 14};
Plane Surface(16) = {15};

Line(15) = {2, 12};
Line(16) = {11, 1};
Line(17) = {3, 13};
Line(18) = {14, 4};

Curve Loop(16) = {17, -12, -15, 3};
Plane Surface(17) = {16};

Curve Loop(17) = {4, 5, -18, -13, -17};
Plane Surface(18) = {17};

Curve Loop(18) = {18, 6, -16, -14};
Plane Surface(19) = {18};

Curve Loop(19) = {1, 2, 15, -11, 16};
Plane Surface(20) = {19};

Surface Loop(1) = {5, 6, 20, 17, 18, 19, 16};
Volume(1) = {1};

Physical Volume("Body") = {1};
Physical Surface("top") = {16};
Physical Surface("bottom") = {5, 6};
Physical Curve("mid line") = {7};
Physical Curve("right line") = {6};
//+
Physical Surface("side left") = {17};
//+
Physical Surface("side right") = {19};
