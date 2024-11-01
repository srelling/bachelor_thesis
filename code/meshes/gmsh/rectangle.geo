Point(1) = {-6,  2, 0, 0.5};
Point(2) = {-6, -2, 0, 0.5};
Point(3) = { 6, -2, 0, 0.5};
Point(4) = { 6,  2, 0, 0.5};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Curve Loop( 9) = {1, 2, 3, 4};
Plane Surface(1) = {9};
Physical Curve("HorEdges", 11) = {1, 3};
Physical Curve("LeftEdge", 12) = {4};
Physical Curve("RightEdge", 13) = {2};
Physical Surface("Rectangle", 4) = {1};

