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

%%
params.consts = Earth;
params.atm_model = @Atm_1962_1976_model;
    
params.nonplanar = false;
params.gravity_model = @Grav_inverse_model;
params.EOM = @EOM_3DOF_planar;
params.consts.Cpmax = 2.0;
params.geom_model = @Geom_conic_model;
params.consts.alpha = 0  * (pi/180);
params.consts.rn = 0.2202;
params.consts.rc = 0.8128/2;
params.consts.dc = 60 * (pi/180);
[CL, CD, A] = params.geom_model(params);
params.consts.m = 46;
params.consts.LD = CL/CD;
params.consts.beta = params.consts.m/CD/A;
params.consts.A = A;
params.consts.sig = 0;
arr = [CL, CD, A, CL/CD, params.consts.m/CD/A]';
T = array2table(arr, "VariableNames", {'Results'}, "RowNames", {'CL', 'CD', 'A (m^2)', 'LD', 'B (kg/m^2)'});
disp(T);
params.init_cond.V = 12.791*1000;
params.init_cond.y = 8.00 * (pi/180);
params.init_cond.h = 125*1000;

params.deorbit.hp = 300;
params.deorbit.e = 0.8;
params.deorbit.thetad = 180 * (pi/180);
params.deorbit.nmax_constraint = 32.71;


params.vary_vel = [-1.5, -0.5, 0.5]*10;
params.vary_gamma = [-1.5, -0.5, 0.5]*10;
params.vary_beta = [-1.5, -0.5, 0.5]*17;
params.vary_LD = [0.1, 0.15, 0.35];
params.vary_sigma = [30, 45, 60] * (pi/180); 
disp("Solving for Deorbit Parameters please wait....")
deorbit_results = Solver_Deorbit(params);
arr = [deorbit_results.Ve, deorbit_results.ye, deorbit_results.deltaV, deorbit_results.ymin, deorbit_results.ymax]';
T = array2table(arr, "VariableNames",{'Results'}, "RowNames", {'Ve (km/s)', 'ye (deg)', 'Delta-V (km/s)', 'ymin (deg)', 'ymax (deg)'});
disp(T);

[t1, sol1, flag] = Solver_EOM_3DOF(params);
[qdot, qint] = Aero_heating_model(t1, sol1.rho, sol1.V, params);
acc = gradient(sol1.V, 0.5);
[nmax, idx] = max(abs(acc./params.consts.g0));
hnmax = sol1.h(idx);
Vnmax = sol1.V(idx);
qdotmax = max(qdot);
qintmax = max(qint);
arr = [nmax, hnmax/1000, Vnmax/1000, qdotmax, qintmax]';
T = array2table(arr, "VariableNames", {'Results'}, 'RowNames', {'nmax (gs)', 'hnmax (km)', 'Vnmax (km/s)', 'qdotmax (W/cm^2)', 'qintmax (J/cm^2)'});
disp(T);
[results] = Solver_Variational_Analysis(params);
arr = [(params.init_cond.V + (params.init_cond.V *params.vary_vel(1)/100))/10^3, results.vary_vel.qdot1, results.vary_vel.qint1;
    (params.init_cond.V + (params.init_cond.V *params.vary_vel(2)/100))/10^3, results.vary_vel.qdot2, results.vary_vel.qint2;
    (params.init_cond.V + (params.init_cond.V *params.vary_vel(3)/100))/10^3, results.vary_vel.qdot3, results.vary_vel.qint3;
    (params.init_cond.y + (params.init_cond.y *params.vary_gamma(1)/100))*(180/pi), results.vary_gamma.qdot1, results.vary_gamma.qint1;
    (params.init_cond.y + (params.init_cond.y *params.vary_gamma(2)/100))*(180/pi), results.vary_gamma.qdot2, results.vary_gamma.qint2;
    (params.init_cond.y + (params.init_cond.y *params.vary_gamma(3)/100))*(180/pi), results.vary_gamma.qdot3, results.vary_gamma.qint3;
    params.consts.beta + (params.consts.beta *params.vary_beta(1)/100), results.vary_beta.qdot1, results.vary_beta.qint1;
    params.consts.beta + (params.consts.beta *params.vary_beta(2)/100), results.vary_beta.qdot2, results.vary_beta.qint2;
    params.consts.beta + (params.consts.beta *params.vary_beta(3)/100), results.vary_beta.qdot3, results.vary_beta.qint3;
    params.vary_LD(1), results.vary_LDsig.qdot1, results.vary_LDsig.qint1;
    params.vary_LD(2), results.vary_LDsig.qdot2, results.vary_LDsig.qint2;
    params.vary_LD(3), results.vary_LDsig.qdot3, results.vary_LDsig.qint3];
T = array2table(arr, 'VariableNames', {'Parameter Value', 'Heating Rate (W/cm^2)','Heating Load (J/cm^2)'},...
    'RowName',{'Ve1 (km/s)','Ve2 (km/s)', 'Ve3 (km/s)', 'ye1 (deg)','ye2 (deg)', 'ye3 (deg)','B1 (kg/m^2)','B2 (kg/m^2)', 'B3 (kg/m^2)', 'LD1/Bank1','LD2/Bank2', 'LD3/Bank3'});
disp(T);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NONPLANAR EQN AND MONTE CARLO
%% USER DEFINED SETTINGS
nonplanar = true;
 
%%
params.consts = Earth;
params.atm_model = @Atm_1962_1976_model;
params.nonplanar = nonplanar;
params.gravity_model = @Grav_J2_model;
params.EOM = @EOM_3DOF_nonplanar;

params.consts.Cpmax = 2;
params.geom_model = @Geom_conic_model;
params.consts.alpha = 0  * (pi/180);
params.consts.rn = 0.2202;
params.consts.rc = 0.8128/2;
params.consts.dc = 60 * (pi/180);

[CL, CD, A] = params.geom_model(params);
params.consts.m = 46;
params.consts.LD = CL/CD;
params.consts.beta = params.consts.m/CD/A;
params.consts.A = A;
params.consts.sig = 0;

params.init_cond.V = 12.791*1000;
params.init_cond.y = -8.00 * (pi/180);
params.init_cond.psi = 96  * (pi/180);
params.init_cond.r = 125*1000 + params.consts.Re;
params.init_cond.phi = 41.59 * (pi/180);
params.init_cond.theta = 236.65 * (pi/180);
params.mc.flag = false;

[t2, sol2, flag] = Solver_EOM_3DOF(params);
params.LB_mod.alpha1 = 3 * (pi/180);
params.LB_mod.sigma1 = 10 * (pi/180);
params.LB_mod.alpha2 = 5 * (pi/180);
params.LB_mod.sigma2 = 20 * (pi/180);
params.LB_mod.alpha3 = 10 * (pi/180);
params.LB_mod.sigma3 = 30 * (pi/180);
[results] = Solver_Lift_Bank_Mod(params);
arr = [params.LB_mod.alpha1*(180/pi), params.LB_mod.sigma1*(180/pi), results.case1.S, results.case1.l, results.case1.nmax, results.case1.hnmax, results.case1.Vnmax;
    params.LB_mod.alpha2*(180/pi), params.LB_mod.sigma2*(180/pi), results.case2.S, results.case2.l, results.case2.nmax, results.case2.hnmax, results.case2.Vnmax;
    params.LB_mod.alpha3*(180/pi), params.LB_mod.sigma3*(180/pi), results.case3.S, results.case3.l, results.case3.nmax, results.case3.hnmax, results.case3.Vnmax];
T = array2table(arr, "VariableNames", {'alpha (deg)', 'sigma (deg)', 'Downrange (km)', 'Crossrange (km)', 'nmax (gs)', 'hnmax (km)', 'Vnmax (km/s)'});
disp(T);

disp("Running Monte Carlo Please Wait");
params.mc.flag = true;
params.mc.V_3sigma = 300;
params.mc.y_3sigma = 0.01*(pi/180);
params.mc.rho_3sigma = 0.02;
params.mc.rho_uncertainty_model = true;
params.mc.lat_3sigma = 0.09 * (pi/180);
params.mc.lon_3sigma = 0.09 * (pi/180);
params.mc.samples = 300;

[lats, lons] = Solver_Monte_Carlo(params);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataStardust;
figure;
plot(sol1.V/1000, sol1.h/1000, '-k', 'LineWidth', 1.5); % Line 2
hold on;
plot(sol2.V/1000, sol2.h/1000, ':k', 'LineWidth', 1.5);
hold on;
plot(Stardust.V, Stardust.h, '--k', 'LineWidth', 1.5);

ax = gca;
ax.FontSize = 12; % Set font size for axis labels and ticks
ax.LineWidth = 1.5; % Set line width for the axis lines

xlabel('Velocity (km/s)', 'FontSize', 14, 'FontWeight', 'bold'); % X-axis label
ylabel('Altitude (km)', 'FontSize', 14, 'FontWeight', 'bold'); % Y-axis label

legend({'Planar EOM', 'Non-Planar EOM', 'Nock et al'}, 'FontSize', 12, 'Location', 'best');

set(gcf, 'Position', [100, 100, 600, 400]); % [left, bottom, width, height]
saveas(gcf, 'Stardust_heating.png');
