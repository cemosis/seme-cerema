myh=1;//0.1;//0.1;
xsize=12;//5;
ysize=12;
zsize=12;
Point(1) = {0, 0, 0, myh};
Point(2) = {xsize, 0, 0, myh};
Point(3) = {xsize, ysize, 0, myh};
Point(4) = {0, ysize, 0, myh};
Point(5) = {0, 0, zsize, myh};
Point(6) = {xsize, 0, zsize, myh};
Point(7) = {xsize, ysize, zsize, myh};
Point(8) = {0, ysize, zsize, myh};
Line(1) = {4, 1};
Line(2) = {2, 1};
Line(3) = {2, 3};
Line(4) = {3, 4};

Line(7) = {8, 7};
Line(8) = {7, 6};
Line(9) = {6, 5};
Line(10) = {5, 8};

Line(11) = {4, 8};
Line(12) = {1, 5};
Line(13) = {3, 7};
Line(14) = {2, 6};


Line Loop(5) = {4, 1, -2, 3};
Plane Surface(6) = {5};
Line Loop(15) = {7, 8, 9, 10};
Plane Surface(16) = {15};
Line Loop(17) = {10, -11, 1, 12};
Plane Surface(18) = {17};
Line Loop(19) = {3, 13, 8, -14};
Plane Surface(20) = {19};
Line Loop(21) = {4, 11, 7, -13};
Plane Surface(22) = {21};
Line Loop(23) = {9, -12, -2, 14};
Plane Surface(24) = {23};
Surface Loop(25) = {16, 22, 6, 18, 24, 20};
Volume(26) = {25};


/* // uncomment this to get hypercube mesh

Transfinite Line {1,3,8,10} = ysize/myh+1;// Using Progression 1.0;
Transfinite Line {2,4,7,9} = xsize/myh+1;// Using Progression 1.0;
Transfinite Line {11,12,13,14} = zsize/myh+1;// Using Progression 1.0;

Transfinite Surface {6} = {1,2,3,4};
Recombine Surface {6};
Transfinite Surface {16} = {5,6,7,8};
Recombine Surface {16};
Transfinite Surface {18} = {5,8,4,1};
Recombine Surface {18};
Transfinite Surface {20} = {3,2,6,7};
Recombine Surface {20};
Transfinite Surface {22} = {3,7,4,8};
Recombine Surface {22};
Transfinite Surface {24} = {1,5,6,2};
Recombine Surface {24};

Transfinite Volume {26} = {1,2,3,4,5,6,7,8};
Recombine Volume {26};
*/

/*
Transfinite Line {1,3,8,10} = ysize/myh+1;// Using Progression 1.0;
Transfinite Line {2,4,7,9} = xsize/myh+1;// Using Progression 1.0;
Transfinite Line {11,12,13,14} = zsize/myh+1;// Using Progression 1.0;

Transfinite Surface {6} = {1,2,3,4};
Transfinite Surface {16} = {5,6,7,8};
Transfinite Surface {18} = {5,8,4,1};
Transfinite Surface {20} = {3,2,6,7};
Transfinite Surface {22} = {3,7,4,8};
Transfinite Surface {24} = {1,5,6,2};

Transfinite Volume {26} = {1,2,3,4,5,6,7,8};


useHexaedralElt=0;

If ( useHexaedralElt )
Recombine Surface {6};
Recombine Surface {16};
Recombine Surface {18};
Recombine Surface {20};
Recombine Surface {22};
Recombine Surface {24};
Recombine Volume {26};
EndIf
*/

Physical Line("bordX") = {2};
Physical Line("bordY") = {1};
Physical Line("bordZ") = {12};
Physical Line("autrebords") = {3,4,5,6,7,8,9,10,11,13,14};

Physical Surface("PartialOmega") = {6,16,18,20,22,24};
Physical Volume("Omega") = {26};
