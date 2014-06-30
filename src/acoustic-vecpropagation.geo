myhsTheta=2*Pi/(4);//0.1;
myhsPhi=Pi/(4*2);//2+1
Point(1) = {0, 0, 0, 1};
Point(2) = {Pi, 0, 0, 1};
Point(3) = {Pi, 2*Pi, 0, 1};
Point(4) = {0, 2*Pi, 0, 1};
Line(1) = {4, 1};
Line(2) = {2, 1};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line Loop(5) = {4, 1, -2, 3};
Plane Surface(6) = {5};
Transfinite Line {1,3} = Pi/(myhsPhi)+1;// Using Progression 1.0;
Transfinite Line {2,4} = 2*Pi/myhsTheta+1;// Using Progression 1.0;
Transfinite Surface {6} = {1,2,3,4};
Recombine Surface {6};

Physical Line("bordX") = {2};
Physical Line("bordY") = {1};
Physical Line("bordXBis") = {4};
Physical Line("bordYBis") = {3};
Physical Surface("Omega") = {6};
