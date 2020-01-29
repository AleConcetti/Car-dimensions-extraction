function [K, iac]=Calibration(vh,vv,vl,v1,v2)
    syms fx fy ux uy;
    K=[fx 0 ux; 0 fy uy; 0 0 1];
    wstar=K*(K.');
    w=inv(wstar);

    eq1 = vh'*w*vv;
    eq2 = vh'*w*vl;
    eq3 = vl'*w*vv;
    eq4 = v1'*w*v2;
    
    eqs=[eq1;eq2;eq3;eq4];

    k_par=solve(eqs, [fx fy ux uy]);
    
    fx=double(k_par.fx);
    fx=abs(fx(1));
    fy=double(k_par.fy);
    fy=abs(fy(1));
    ux=double(k_par.ux);
    ux=abs(ux(1));
    uy=double(k_par.uy);
    uy=abs(uy(1));
    K=[fx  0   ux;...
        0  fy  uy;...
        0  0   1 ];
    
    diac=K*(K.');
    diac=diac/diac(3,3);
    
    iac=inv(diac); 
    iac=iac/iac(3,3);

end


