h = 0.5;

Point(1) = { 1., 0., 0., h};
Point(3) = { 0., 0., 0., h};
Point(4) = { 0., 0., 0., h};
Point(2) = {-1., 0., 0., h};

Line(1) = {2, 3};
Line(2) = {4, 1};

Physical Point ("loading") = {1};
Physical Point ("fixed") = {2};

Physical Point ("master") = {3};
Physical Point ("slave") = {4};

Physical Line ("master_body") = {1};
Physical Line ("slave_body") = {2};