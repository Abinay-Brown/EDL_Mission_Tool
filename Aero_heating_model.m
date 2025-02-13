function [qdot, qint] = Aero_heating_model(t, rho, V, params)
%% t: time history (sec)
%% rho: density history (kg/m^3)
%% V: velocity (m/s)
%% qdot: W/cm^2
%% qint: J/cm^2

K = params.consts.Ksg;  % Sutton Graves constant
rn = params.consts.rn; % m

qdot = K .* sqrt(rho./rn) .* V.^3 / (10^4);
qint = cumtrapz(t, qdot);

    
end