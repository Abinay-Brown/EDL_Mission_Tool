clear; close all; clc;
Planetary_Constants;
planet = inputdlg("Enter Entry Planet (1: Earth, 2:Mars, 3: Venus, 4: Stardust)");
if planet{1} ~= '4'
    nonplanar = inputdlg("Planar or Nonplanar Equations of Motion? (0: Planar,  1: Non-Planar)");
    if nonplanar{1} == '1'
        params.nonplanar = true;
    else
        params.nonplanar = false;
    end
    if planet{1} == '1'
        params.consts = Earth;
        atms = inputdlg("Atmosphere model (1: Exponential, 2: 1976-1962 model)");
        if atms{1} == '1'
            params.atm_model = @Atm_exponential_model;
        elseif atms{1} == '2'
            params.atm_model = @Atm_1962_1976_model;
        end
    elseif planet{1} == '2'
        params.consts = Mars;
        params.atm_model = @Atm_exponential_model; 
    
    elseif planet{1} == '3'
        params.consts = Venus;
        params.atm_model = @Atm_exponential_model;
    end
    if nonplanar{1} == '1'
        params.gravity_model = @Grav_J2_model;
        params.EOM = @EOM_3DOF_nonplanar;
    elseif nonplanar{1} == '0'
        params.gravity_model = @Grav_inverse_model;
        params.EOM = @EOM_3DOF_planar;
    end
    geometry = inputdlg("Enter Geometry (1: Conic, 2: Biconic)");
    if geometry{1} == '1'
        params.geom_model = @Geom_conic_model;
        geom = inputdlg({'mass (kg)', 'alpha (deg)','Cpmax','Nose Radius (m)', 'Cone Radius (m)', 'Cone Angle (deg)' });
        params.geom_model = @Geom_conic_model;
        params.consts.m = str2double(geom{1});
        params.consts.alpha = str2double(geom{2}) * (pi/180);
        params.consts.Cpmax = str2double(geom{3}); 
        params.consts.rn = str2double(geom{4});
        params.consts.rc = str2double(geom{5});
        params.consts.dc = str2double(geom{6})* (pi/180);
        
    elseif geometry{1} == '2'
        params.geom_model = @Geom_biconic_model;
        geom = inputdlg({'mass (kg)', 'alpha (deg)','Cpmax','Nose Radius (m)', 'Cone Radius 1 (m)', 'Cone Radius 2 (m)', 'Cone Angle 1 (deg)', 'Cone Angle 2 (deg)' });
        params.geom_model = @Geom_biconic_model;
        params.consts.m = str2double(geom{1});
        params.consts.alpha = str2double(geom{2}) * (pi/180);
        params.consts.Cpmax = str2double(geom{3}); 
        params.consts.rn = str2double(geom{4});
        params.consts.rc1 = str2double(geom{5});
        params.consts.rc2 = str2double(geom{6});
        params.consts.dc1 = str2double(geom{7})* (pi/180);
        params.consts.dc2 = str2double(geom{8}) * (pi/180);
    end
    
    [CL, CD, A] = params.geom_model(params);
    params.consts.LD = CL/CD;
    params.consts.beta = params.consts.m/CD/A;
    params.consts.A = A;
    params.consts.sig = 0;
    arr = [CL, CD, A, CL/CD, params.consts.m/CD/A]';
    T = array2table(arr, "VariableNames", {'Results'}, "RowNames", {'CL', 'CD', 'A (m^2)', 'LD', 'B (kg/m^2)'});
    disp(T);
    
    params.mc.flag = false;
    deorbit = inputdlg('1: Determine Deorbit Conditions or 2: Bypass to entry');
    if deorbit{1} == '1'
        deorbit_cond = inputdlg({'Periapsis height (km)', 'Eccentricity', 'Nmax Limit (gs)'});
        params_copy = params;
        params_copy.nonplanar = false;
        params_copy.deorbit.hp = str2double(deorbit_cond{1});
        params_copy.deorbit.e = str2double(deorbit_cond{2});
        params_copy.deorbit.thetad = 180 * (pi/180);
        params_copy.deorbit.nmax_constraint = str2double(deorbit_cond{3});
        params_copy.gravity_model = @Grav_inverse_model;
        params_copy.atm_model = @Atm_exponential_model;
        params_copy.EOM = @EOM_3DOF_planar;
        disp("Solving for Deorbit Parameters please wait....")
        [deorbit_results] = Solver_Deorbit(params_copy);
        arr = [deorbit_results.Ve, deorbit_results.ye, deorbit_results.deltaV, deorbit_results.ymin, deorbit_results.ymax]';
        T = array2table(arr, "VariableNames",{'Results'}, "RowNames", {'Ve (km/s)', 'ye (deg)', 'Delta-V (km/s)', 'ymin (deg)', 'ymax (deg)'});
        disp(T);
        use_entry = inputdlg('(1) Use Deorbit Entry Results  or (0) Bypass to Entry');
        
    end
    if nonplanar{1} == '1'
        
        if use_entry{1} == '1'
            entry_cond = inputdlg({'azimuth psi (deg)', 'entry h (km)', 'entry lat phi (deg)', 'entry lon theta (deg)'}); 
            params.init_cond.V = deorbit_results.Ve * 1000;
            params.init_cond.y = - deorbit_results.ye * (pi/180);
            params.init_cond.psi = str2double(entry_cond{1})  * (pi/180);
            params.init_cond.r = str2double(entry_cond{2})*1000 + params.consts.Re;
            params.init_cond.phi = str2double(entry_cond{3}) * (pi/180);
            params.init_cond.theta = str2double(entry_cond{4}) * (pi/180);
        else
            entry_cond = inputdlg({'V (km/s)', '+ve y (deg)', 'azimuth psi (deg)', 'entry h (km)', 'entry lat phi (deg)', 'entry lon theta (deg)'}); 
            params.init_cond.V = str2double(entry_cond{1}) * 1000;
            params.init_cond.y = -str2double(entry_cond{2}) * (pi/180);
            params.init_cond.psi = str2double(entry_cond{3})  * (pi/180);
            params.init_cond.r = str2double(entry_cond{4})*1000 + params.consts.Re;
            params.init_cond.phi = str2double(entry_cond{5}) * (pi/180);
            params.init_cond.theta = str2double(entry_cond{6}) * (pi/180);
        end
    elseif nonplanar{1} == '0'
        if use_entry{1} == '1'
            entry_alt = inputdlg("Enter entry altitude (km)"); 
            params.init_cond.V = deorbit_results.Ve * 1000;
            params.init_cond.y = deorbit_results.ye * (pi/180);
            params.init_cond.h = str2double(entry_alt{1})*1000;
        else
            entry_cond = inputdlg({'Ve (km/s)', '+ve ye (deg)', 'h (km)'});
            params.init_cond.V = str2double(entry_cond{1}) * 1000;
            params.init_cond.y = str2double(entry_cond{2}) * (pi/180);
            params.init_cond.h = str2double(entry_cond{3}) * 1000;
        end
    end
    [t, sol, flag] = Solver_EOM_3DOF(params);
    [qdot, qint] = Aero_heating_model(t, sol.rho, sol.V, params);

    figure(1);
    plot(sol.V/1000, sol.h/1000, '-k', 'LineWidth', 1.5); % Line 1
    
    ax = gca;
    ax.FontSize = 12; % Set font size for axis labels and ticks
    ax.LineWidth = 1.5; % Set line width for the axis lines
    
    xlabel('Velocity (km/s)', 'FontSize', 14, 'FontWeight', 'bold'); % X-axis label
    ylabel('Altitude (km)', 'FontSize', 14, 'FontWeight', 'bold'); % Y-axis label
    
    acc = gradient(sol.V, 0.5);
    [nmax, idx] = max(abs(acc./params.consts.g0));
    hnmax = sol.h(idx);
    Vnmax = sol.V(idx);
    qdotmax = max(qdot);
    qintmax = max(qint);
    arr = [nmax, hnmax/1000, Vnmax/1000, qdotmax, qintmax]';
    T = array2table(arr, "VariableNames", {'Results'}, 'RowNames', {'nmax (gs)', 'hnmax (km)', 'Vnmax (km/s)', 'qdotmax (W/cm^2)', 'qintmax (J/cm^2)'});
    disp(T);

    use_vary = inputdlg('Run Parameter Varying Analysis for Heating? (1: Yes, 0: No)');
    if use_vary{1} == '1'
        vary_params = inputdlg({'V1 (%change)', 'V2 (%change)', 'V3 (%change)', 'y1 (%change)', 'y2 (%change)', 'y3 (%change)', 'B1 (%change)', 'B2 (%change)',...
            'B3 (%change)', 'L/D1', 'L/D2', 'L/D3', 'Bank1 (deg)', 'Bank2 (deg)', 'Bank3 (deg)'});

        params.vary_vel = [str2double(vary_params{1}), str2double(vary_params{2}), str2double(vary_params{3})];
        params.vary_gamma = [str2double(vary_params{4}), str2double(vary_params{5}), str2double(vary_params{6})];
        params.vary_beta = [str2double(vary_params{7}), str2double(vary_params{8}), str2double(vary_params{9})];
        params.vary_LD = [str2double(vary_params{10}), str2double(vary_params{11}), str2double(vary_params{12})];
        params.vary_sigma = [str2double(vary_params{13}), str2double(vary_params{14}), str2double(vary_params{15})] * (pi/180);
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
    end
    
    if nonplanar{1} == '1'
        use_mod = inputdlg("Perform Lift Bank Modulation (1: yes, 0: no)");
        if use_mod{1} == '1'
            LB = inputdlg({'Case 1 Alpha (deg)', 'Case 1 sigma (deg)', 'Case 2 Alpha (deg)', 'Case 2 sigma (deg)', 'Case 3 Alpha (deg)', 'Case 3 sigma (deg)'});
            params.LB_mod.alpha1 = str2double(LB{1}) * (pi/180);
            params.LB_mod.sigma1 = str2double(LB{2}) * (pi/180);
            params.LB_mod.alpha2 = str2double(LB{3}) * (pi/180);
            params.LB_mod.sigma2 = str2double(LB{4}) * (pi/180);
            params.LB_mod.alpha3 = str2double(LB{5}) * (pi/180);
            params.LB_mod.sigma3 = str2double(LB{6}) * (pi/180);
            [results] = Solver_Lift_Bank_Mod(params);
            arr = [params.LB_mod.alpha1*(180/pi), params.LB_mod.sigma1*(180/pi), results.case1.S, results.case1.l, results.case1.nmax, results.case1.hnmax, results.case1.Vnmax;
                   params.LB_mod.alpha2*(180/pi), params.LB_mod.sigma2*(180/pi), results.case2.S, results.case2.l, results.case2.nmax, results.case2.hnmax, results.case2.Vnmax;
                   params.LB_mod.alpha3*(180/pi), params.LB_mod.sigma3*(180/pi), results.case3.S, results.case3.l, results.case3.nmax, results.case3.hnmax, results.case3.Vnmax];
            T = array2table(arr, "VariableNames", {'alpha (deg)', 'sigma (deg)', 'Downrange (km)', 'Crossrange (km)', 'nmax (gs)', 'hnmax (km)', 'Vnmax (km/s)'});
            disp(T);
        end
    end
    if nonplanar{1} == '1'
        use_monte = inputdlg("Perform Monte-Carlo Runs (1: yes, 0: no)");
        if use_monte{1} == '1'
            MC = inputdlg({'Vel (3-sig) error bounds (m/s)', 'y (3-sig) error bounds (deg)', 'rho (3-sig) error bounds (kg/m^3)',...
                'Lat (3-sig) error bounds (deg)', 'Lon (3-sig) error bounds (deg)', 'Total Samples'});
            params.mc.flag = true;
            params.mc.V_3sigma = str2double(MC{1});
            params.mc.y_3sigma = str2double(MC{2})*(pi/180);
            params.mc.rho_3sigma = str2double(MC{3});
            params.mc.rho_uncertainty_model = true;
            params.mc.lat_3sigma =str2double(MC{4}) * (pi/180);
            params.mc.lon_3sigma = str2double(MC{5}) * (pi/180);
            params.mc.samples = str2double(MC{6});
            [lats, lons] = Solver_Monte_Carlo(params);
        end
    end
elseif planet{1} == '4'
    Main_Stardust;
end

