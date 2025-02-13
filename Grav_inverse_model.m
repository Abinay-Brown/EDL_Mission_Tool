function g = Grav_inverse_model(h, params)
%% h: altitude (m)
%% params: [g0 (surface gravity constant), Re (radius (km))]
%% g: gravitational acceleration (g)
g0 = params.consts.g0;
Re = params.consts.Re;

r = Re + h;
g = g0 .* (Re./r).^2;
    
end