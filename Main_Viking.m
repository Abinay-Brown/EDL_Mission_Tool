%% Planetary Entry Descent Computer Project
%% Name: Abinay Brown
%% Email: abrown472@gatech.edu
clear;
clc;
close all;
format longg;
Planetary_Constants;
Constants;
%% USER DEFINED SETTINGS
Planet = "Mars";
nonplanar = false;
geometry = "Conic";
Cpmax = 2.0;
mass = 930;

%%
if strcmpi(Planet, "Earth")
    params.consts = Earth;
    params.atm_model = @Atm_exponential_model;
    params.atm_params = [params.consts.H, params.consts.rho0];
    % params.atm_model = @Atm_1962_1976_model;
    % params.atm_params = [consts.Rbar, consts.r0, consts.g0];
    
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
    params.consts.alpha = -11.2  * (pi/180);
    params.consts.rn = 0.8763;
    params.consts.rc = 3.505/2;
    params.consts.dc = 70 * (pi/180);

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

params.init_cond.V = 4.7*1000;
params.init_cond.y = 17.6 * (pi/180);
params.init_cond.h = 120*1000;

[t1, sol1, flag] = Solver_EOM_3DOF(params);
[qdot1, qint1] = Aero_heating_model(t1, sol1.rho, sol1.V, params);

save("Mars.mat", "t1", "sol1");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
format longg;
Planetary_Constants;
Constants;
%% USER DEFINED SETTINGS
Planet = "Mars";
nonplanar = true;
geometry = "Conic";
Cpmax = 2.0;
mass = 930;

%%
if strcmpi(Planet, "Earth")
    params.consts = Earth;
    params.atm_model = @Atm_exponential_model;
    params.atm_params = [params.consts.H, params.consts.rho0];
    % params.atm_model = @Atm_1962_1976_model;
    % params.atm_params = [consts.Rbar, consts.r0, consts.g0];
    
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
    params.consts.alpha = -11.2  * (pi/180);
    params.consts.rn = 0.8763;
    params.consts.rc = 3.505/2;
    params.consts.dc = 70 * (pi/180);

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

params.init_cond.V = 4.7*1000;
params.init_cond.y = -17.6 * (pi/180);
params.init_cond.psi = 0  * (pi/180);
params.init_cond.r = 120*1000 + params.consts.Re;
params.init_cond.phi = 0;
params.init_cond.theta = 0;
params.mc.flag = false;

[t2, sol2, flag] = Solver_EOM_3DOF(params);
[qdot2, qint2] = Aero_heating_model(t2, sol2.rho, sol2.V, params);
% plot(sol2.V, sol2.h);
% hold on;
save("Mars.mat", "t1", "t2", "sol1", "sol2", "qdot1", "qdot2", "qint1", "qint2");
Plots_Viking;