// Gmsh project created on Tue Dec 10 01:15:56 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, -0, 1.0};
//+
Point(2) = {0, 1, -0, 1.0};
//+
Point(3) = {2, 1, -0, 1.0};
//+
Point(4) = {2, 0, -0, 1.0};
//+
Point(5) = {0, 0.4, -0, 1.0};
//+
Point(6) = {0.5, 0.4, -0, 1.0};
//+
Point(7) = {1.5, 0.4, -0, 1.0};
//+
Point(8) = {2, 0.4, -0, 1.0};
//+
Point(9) = {0, 0.6, -0, 1.0};
//+
Point(10) = {0.6, 0.6, -0, 1.0};
//+
Point(11) = {1.4, 0.6, -0, 1.0};
//+
Point(12) = {2, 0.6, -0, 1.0};
//+
Point(13) = {0.6, 1, -0, 1.0};
//+
Point(14) = {1.4, 1, -0, 1.0};
//+
Point(15) = {0.5, 0.6, -0, 1.0};
//+
Point(16) = {1.5, 0.6, 0, 1.0};
//+
Line(1) = {1, 5};
//+
Line(2) = {5, 9};
//+
Line(3) = {9, 2};
//+
Line(4) = {2, 13};
//+
Line(5) = {13, 14};
//+
Line(6) = {14, 3};
//+
Line(7) = {3, 12};
//+
Line(8) = {12, 8};
//+
Line(9) = {8, 4};
//+
Line(10) = {4, 1};
//+
Line(11) = {5, 6};
//+
Line(12) = {6, 7};
//+
Line(13) = {7, 8};
//+
Line(14) = {9, 15};
//+
Line(15) = {15, 10};
//+
Line(16) = {10, 11};
//+
Line(17) = {11, 16};
//+
Line(18) = {16, 12};
//+
Line(19) = {6, 15};
//+
Line(20) = {7, 16};
//+
Line(21) = {10, 13};
//+
Line(22) = {11, 14};
//+
Point(17) = {0.95, 0.1, 0, 1.0};
//+
Point(18) = {0.95, 0.2, 0, 1.0};
//+
Point(19) = {1.05, 0.2, 0, 1.0};
//+
Point(20) = {1.05, 0.1, 0, 1.0};
//+
Line(23) = {17, 18};
//+
Line(24) = {18, 19};
//+
Line(25) = {19, 20};
//+
Line(26) = {20, 17};
//+
Curve Loop(1) = {1, 11, 12, 13, 9, 10};
//+
Curve Loop(2) = {24, 25, 26, 23};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {2, 14, -19, -11};
//+
Plane Surface(2) = {3};
//+
Curve Loop(4) = {12, 20, -17, -16, -15, -19};
//+
Plane Surface(3) = {4};
//+
Curve Loop(5) = {20, 18, 8, -13};
//+
Plane Surface(4) = {5};
//+
Curve Loop(6) = {3, 4, -21, -15, -14};
//+
Plane Surface(5) = {6};
//+
Curve Loop(7) = {21, 5, -22, -16};
//+
Plane Surface(6) = {7};
//+
Curve Loop(8) = {22, 6, 7, -18, -17};
//+
Plane Surface(7) = {8};
//+
Curve Loop(9) = {23, 24, 25, 26};
//+
Plane Surface(8) = {9};
//+
Physical Surface("Bot", 1) = {1};
//+
Physical Surface("Mid Left", 2) = {2};
//+
Physical Surface("Mid Mid", 3) = {3};
//+
Physical Surface("Mid Right", 4) = {4};
//+
Physical Surface("Top Left", 5) = {5};
//+
Physical Surface("Top Mid", 6) = {6};
//+
Physical Surface("Top Right", 7) = {7};
//+
Physical Surface("Waste", 8) = {8};
