h = 4.0;
height = 1.0;
length = 1.0;
epsilon=1e-3;

Point(1) = {0,  0, 0, h*0.1};
Point(2) = {length, 0, 0, h*0.1};
Point(3) = {length, height/2-epsilon, 0, h*0.1};
Point(4) = {length, height, 0, h*0.2};
Point(5) = {0, height, 0, h*0.2};
Point(6) = {0, height/2+epsilon, 0, h*0.2};
Point(7) = {length, height/2+epsilon, 0, h*0.2};
Point(8) = {0, height/2-epsilon, 0, h*0.1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {7, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {8, 1};
Line(8) = {3, 8};

Line Loop(1) = {1, 2, 8, 7};
Line Loop(2) = {3, 4, 5, 6};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Line ("top") = {4};
Physical Line ("bottom") = {1};

Physical Line ("contact_top") = {6};
Physical Line ("contact_bottom") = {8};

Physical Surface("bot_body") = {1};
Physical Surface("top_body") = {2};

Transfinite Surface{1, 2};
Recombine Surface{1, 2};