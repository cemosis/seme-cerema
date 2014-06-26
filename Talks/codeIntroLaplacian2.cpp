    // forme lineaire : $l(v)=\int_\Omega v$
    auto l = form1(_test=Vh);
    l = integrate(_range=elements(mesh),
		  _expr=id(v));

    // forme bilineaire : $a(u,v)=\int_\Omega \nabla u \cdot \nabla v$
    auto a = form2(_trial=Vh,_test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=inner(gradt(u),grad(v)));

    // conditions aux limites : $u = 0 \ \text{sur} \ \partial\Omega$
    a+=on(_range=boundaryfaces(mesh),
	  _rhs=l,_element=u,
	  _expr=cst(0.) );

    // resolution du systeme algebrique : $A U = F$
    a.solve(_rhs=l,_solution=u);
