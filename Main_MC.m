%% Planetary Entry Descent Computer Project
%% Name: Abinay Brown
%% Email: abrown472@gatech.edu
clear;
clc;
% close all;
format longg;
Planetary_Constants;
Constants;
%% USER DEFINED SETTINGS
Planet = "Earth";
nonplanar = true;
geometry = "Conic";
Cpmax = 2.0;
mass = 46;

%%
if strcmpi(Planet, "Earth")
    params.consts = Earth;
    params.atm_model = @Atm_exponential_model;
    params.atm_model = @Atm_1962_1976_model;
    
elseif strcmpi(Planet, "Mars")
    params.consts = Mars;
    params.atm_model = @Atm_exponential_model; 
elseif strcmpi(Planet, "Venus")
    params.consts = Venus;
    params.atm_model = @Atm_exponential_model;
end

params.nonplanar = nonplanar;
if nonplanar == true
    params.gravity_model = @Grav_J2_model;
    params.EOM = @EOM_3DOF_nonplanar;
else
    params.gravity_model = @Grav_inverse_model;
    params.EOM = @EOM_3DOF_planar;
end

params.consts.Cpmax = Cpmax;
if strcmpi(geometry, "Conic")
    % Using Conic Model
    params.geom_model = @Geom_conic_model;
    params.consts.alpha = 0  * (pi/180);
    params.consts.rn = 0.2202;
    params.consts.rc = 0.8128/2;
    params.consts.dc = 60 * (pi/180);

elseif strcmpi(geometry, "Biconic")
    % Using Biconic Model
    params.geom_model = @Geom_conic_model;
    params.consts.alpha = 11.2 * (pi/180);
    params.consts.rn = 0.2202;
    params.consts.rc1 = 0.8128/2;
    params.consts.rc2 = 0.8128/2;
    params.consts.dc1 = 60 * (pi/180);
    params.consts.dc2 = 60 * (pi/180);
end

[CL, CD, A] = params.geom_model(params);
params.consts.m = mass;
params.consts.LD = CL/CD;
params.consts.beta = params.consts.m/CD/A;
% params.consts.beta = 63.7;
params.consts.A = A;
params.consts.sig = 0;

params.init_cond.V = 12.791*1000;
params.init_cond.y = -8.00 * (pi/180);
params.init_cond.psi = 96  * (pi/180);
params.init_cond.r = 125*1000 + params.consts.Re;
params.init_cond.phi = 41.59 * (pi/180);
params.init_cond.theta = 236.65 * (pi/180);
params.mc.flag = false;

%% Monte Carlo Analysis
params.mc.flag = true;
params.mc.V_3sigma = 300;
params.mc.y_3sigma = 0.01*(pi/180);
params.mc.rho_3sigma = 0.02;
params.mc.rho_uncertainty_model = true;
params.mc.lat_3sigma = 0.09 * (pi/180);
params.mc.lon_3sigma = 0.09 * (pi/180);
params.mc.samples = 300;
[lats, lons] = Solver_Monte_Carlo(params);
% scatter(lons, lats);
%% Variational Analysis
% Re = consts.Req .* (1 - (consts.k.*(sin(sol(:, 5)).^2))) * 1000;
% params.vary_vel = [-1.5, -0.5, 0, 0.5, 1.5]*10;
% params.vary_gamma = [-1.5, -0.5, 0, 0.5, 1.5]*10;
% params.vary_beta = [-1.5, -0.5, 0, 0.5, 1.5]*30;
% 
% params.vary_sigma = [10, 45] *(pi/180);
% params.vary_LD = [0, 75];
% 
% 
% Solver_Variational_Analysis(params);