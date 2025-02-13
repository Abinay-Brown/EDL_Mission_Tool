function g = Grav_J2_model(h, phi, params)
%% h: altitude (m)
%% params: [mu (gravity param (km^3/s^2)), J2 (coeff), w (rot rate rad/s), Req (Radius eq (km), K (ellipticity)]
%% g: gravitational acceleration (m/s^2)
mu = params.consts.mu;
J2 = params.consts.J2;
w = params.consts.w;
Req = params.consts.Req;
K = params.consts.k;

Re = Req * (1 - (K*(sin(phi)^2))) * 1000;
D = K * (1 - (h/Re))*sin(2*phi); 
R = Re + h;
g = (10^9)*(mu/(R^2))*(1 - (3*J2/4)*(1 - (3*cos(2*phi)))) - R*w*w*cos(phi)*cos(phi - D);

end