h = 1.0;
upper_block_factor = 1.0;
lower_block_factor = 1.0;

height = 1;
upper_length = 1;
lower_height = 5*upper_length;
lower_length = 7.5*upper_length;
epsilon= 0;

nb_cells = 2;

//Point(1) = {-lower_length,  0, 0, lower_block_factor};
Point(2) = {lower_length, lower_height/2, 0, lower_block_factor};
Point(3) = {lower_length, lower_height-epsilon, 0, lower_block_factor};
Point(4) = {upper_length, lower_height+height, 0, upper_block_factor};
Point(5) = {0, lower_height+height, 0, upper_block_factor};
Point(6) = {0, lower_height+epsilon, 0, upper_block_factor};
Point(7) = {upper_length, lower_height+epsilon, 0, upper_block_factor};
//Point(8) = {-lower_length, lower_height-epsilon, 0, lower_block_factor};
//Point(9) = {4*upper_length, lower_height/2, 0, lower_block_factor};
//Point(10) = {4*upper_length, lower_height, 0, lower_block_factor};
Point(11) = {-lower_length/4, lower_height/2, 0, lower_block_factor};
Point(12) = {-lower_length/4, lower_height, 0, lower_block_factor};



//Line(1) = {1, 11};
Line(2) = {2, 3};
Line(3) = {7, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {12, 11};
//Line(8) = {3, 10};
//Line(10) = {9, 2};
Line(11) = {3, 12};
Line(12) = {11, 2};
//Line(13) = {12, 8};

Line Loop(1) = {12,  2,  11, 7};
Line Loop(2) = {3, 4, 5, 6};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Line ("loading") = {4};
Physical Line ("fixed") = {12};

Physical Line ("contact_top") = {6};
Physical Line ("contact_bottom") = {11};


Physical Surface("lower") = {1};
Physical Surface("upper") = {2};

Transfinite Line {12, -11} = 20 Using Progression 1;

Transfinite Line {7, -2} = 5 Using Progression 1;

Transfinite Line {6, -4} = nb_cells Using Progression 1;
Transfinite Line {-3, 5} = nb_cells Using Progression 1;

Transfinite Surface{1} = {11, 2, 3, 12};
Transfinite Surface{2};

Recombine Surface{1, 2};