function rho = Atm_exponential_model(h, params)
%% h: Geometric Height (m)
%% params: [H (scale height, km), rho0 (surface air density, kg/m^3)]
H = params.consts.H;
rho0 = params.consts.rho0;
rho = rho0 .* exp(-h./H./1000);

end