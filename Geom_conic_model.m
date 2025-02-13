function [CL, CD, A] = Geom_conic_model(params)
    Cpmax = params.consts.Cpmax;
    alpha = params.consts.alpha;
    rn = params.consts.rn;
    rc = params.consts.rc;
    dc = params.consts.dc;

    CN = Cpmax*(1 - ((rn/rc)^2)*(cos(dc)^2))*((cos(dc)^2)*sin(alpha)*cos(alpha));
    CA = Cpmax*((0.5*(1-(sin(dc)^4))*(rn/rc)^2) + ((sin(dc)^2)*(cos(alpha)^2) + 0.5*(sin(alpha)^2)*(cos(dc)^2))*(1-((rn/rc)^2)*(cos(dc)^2)));
    CL = CN*cos(alpha) - CA*sin(alpha);
    CD = CN*sin(alpha) + CA*cos(alpha);
    A = pi*rc^2;
end