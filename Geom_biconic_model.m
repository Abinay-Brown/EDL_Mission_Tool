function [CL, CD, A] = Geom_biconic_model(params)
    Cpmax = params.consts.Cpmax;
    alpha = params.consts.alpha;
    rn = params.consts.rn;
    rc1 = params.consts.rc1;
    rc2 = params.consts.rc2;
    d1 = params.consts.dc1;
    d2 = params.consts.dc2;
        
    CN1 = (Cpmax*((rc1^2)-(rn*cos(d1))^2)*sin(2*alpha)*cos(d1)^2)/(2*rc1^2);
    CN2 = (Cpmax*((rc2^2)-(rc1^2))*sin(2*alpha)*cos(d2)^2)/(2*rc2^2);
    CN = ((rc1/rc2)^2)*CN1 + CN2;
     
    CAnose = 0.5*Cpmax*(1-(sin(d1)^4));
    
    CA1 = Cpmax*((rc1^2)-(rn*cos(d1))^2)*(((cos(d1)*sin(alpha))^2)+(2*(cos(alpha)*sin(d1))^2))/(2*rc1^2);
    CA2 = Cpmax*((rc2^2)-(rc1^2))*(((cos(d2)*sin(alpha))^2)+(2*(cos(alpha)*sin(d2))^2))/(2*rc2^2);
    CA = ((rn/rc2)^2)*CAnose + ((rc1/rc2)^2)*CA1 + CA2;
        
    CL = CN*cos(alpha) - CA*sin(alpha);
    CD = CN*sin(alpha) + CA*cos(alpha);
    A = pi*rc2^2;
end