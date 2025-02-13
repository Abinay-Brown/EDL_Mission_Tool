function [results] = Solver_Variational_Analysis(params)
    dt = 0.1;
    tspan = 0:dt:1000;
    
    %% Varying Entry Velocity
    fig1 = figure();
    legend_sent = [];
 
    params_copy = params;
    params_copy.init_cond.V = params_copy.init_cond.V + (params_copy.init_cond.V * params.vary_vel(1)/100);
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot1, qint1] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot1, '--k', 'LineWidth', 1.5);
    sent1 = sprintf("Ve = %0.3f km/s", params_copy.init_cond.V/1000);
    hold on;
    results.vary_vel.qdot1 = max(qdot1);
    results.vary_vel.qint1 = max(qint1);

    params_copy = params;
    params_copy.init_cond.V = params_copy.init_cond.V + (params_copy.init_cond.V * params.vary_vel(2)/100);
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot2, qint2] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot2, ':k', 'LineWidth', 1.5);
    sent2 = sprintf("Ve = %0.3f km/s", params_copy.init_cond.V/1000);
    hold on;
    results.vary_vel.qdot2 = max(qdot2);
    results.vary_vel.qint2 = max(qint2);

    params_copy = params;
    params_copy.init_cond.V = params_copy.init_cond.V + (params_copy.init_cond.V * params.vary_vel(3)/100);
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot3, qint3] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot3, '-k', 'LineWidth', 1.5);
    sent3 = sprintf("Ve = %0.3f km/s", params_copy.init_cond.V/1000);
    hold on;
    results.vary_vel.qdot3 = max(qdot3);
    results.vary_vel.qint3 = max(qint3);
    ax = gca;
    ax.FontSize = 12; % Set font size for axis labels and ticks
    ax.LineWidth = 1.5; % Set line width for the axis lines

    legend([sent1, sent2, sent3]);
    xlabel('Time (sec)')
    ylabel('Stagnation Point Heating Rate W/cm^2');
    xlim([0 200]);

    %% Varying Entry Flight Path Angle
    fig2 = figure();
    legend_sent = [];
 
    params_copy = params;
    params_copy.init_cond.y = params_copy.init_cond.y + (params_copy.init_cond.y * params.vary_gamma(1)/100);
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot1, qint1] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot1, '--k', 'LineWidth', 1.5);
    sent1 = sprintf("ye = %0.3f deg", params_copy.init_cond.y*(180/pi));
    hold on;
    results.vary_gamma.qdot1 = max(qdot1);
    results.vary_gamma.qint1 = max(qint1);
    
    params_copy = params;
    params_copy.init_cond.y = params_copy.init_cond.y + (params_copy.init_cond.y * params.vary_gamma(2)/100);
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot2, qint2] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot2, ':k', 'LineWidth', 1.5);
    sent2 = sprintf("ye = %0.3f deg", params_copy.init_cond.y*(180/pi));
    hold on;
    results.vary_gamma.qdot2 = max(qdot2);
    results.vary_gamma.qint2 = max(qint2);

    params_copy = params;
    params_copy.init_cond.y = params_copy.init_cond.y + (params_copy.init_cond.y * params.vary_gamma(3)/100);
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot3, qint3] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot3, '-k', 'LineWidth', 1.5);
    sent3 = sprintf("ye = %0.3f deg", params_copy.init_cond.y*(180/pi));
    hold on;
    results.vary_gamma.qdot3 = max(qdot3);
    results.vary_gamma.qint3 = max(qint3);
    ax = gca;
    ax.FontSize = 12; % Set font size for axis labels and ticks
    ax.LineWidth = 1.5; % Set line width for the axis lines

    legend([sent1, sent2, sent3]);
    xlabel('Time (sec)')
    ylabel('Stagnation Point Heating Rate W/cm^2');
    xlim([0 200]);
    
    %% Varying Entry Beta
    fig2 = figure();
    legend_sent = [];
 
    params_copy = params;
    params_copy.consts.beta = params_copy.consts.beta + (params_copy.consts.beta * params.vary_beta(1)/100);
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot1, qint1] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot1, '--k', 'LineWidth', 1.5);
    sent1 = sprintf("B = %0.3f kg/m^2", params_copy.consts.beta);
    hold on;
    results.vary_beta.qdot1 = max(qdot1);
    results.vary_beta.qint1 = max(qint1);
    
    params_copy = params;
    params_copy.consts.beta = params_copy.consts.beta + (params_copy.consts.beta * params.vary_beta(2)/100);
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot2, qint2] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot2, ':k', 'LineWidth', 1.5);
    sent2 = sprintf("B = %0.3f kg/m^2", params_copy.consts.beta);
    hold on;
    results.vary_beta.qdot2 = max(qdot2);
    results.vary_beta.qint2 = max(qint2);

    params_copy = params;
    params_copy.consts.beta = params_copy.consts.beta + (params_copy.consts.beta * params.vary_beta(3)/100);
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot3, qint3] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot3, '-k', 'LineWidth', 1.5);
    sent3 = sprintf("B = %0.3f kg/m^2", params_copy.consts.beta);
    hold on;
    results.vary_beta.qdot3 = max(qdot3);
    results.vary_beta.qint3 = max(qint3);
    ax = gca;
    ax.FontSize = 12; % Set font size for axis labels and ticks
    ax.LineWidth = 1.5; % Set line width for the axis lines

    legend([sent1, sent2, sent3]);
    xlabel('Time (sec)')
    ylabel('Stagnation Point Heating Rate W/cm^2');
    xlim([0 200]);

    %% Varying LD and sigma
    fig2 = figure();
    legend_sent = [];
 
    params_copy = params;
    params_copy.consts.LD = params_copy.vary_LD(1); 
    params_copy.consts.sig = params_copy.vary_sigma(1); 
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot1, qint1] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot1, '--k', 'LineWidth', 1.5);
    sent1 = sprintf("L/D = %0.3f, Bank = %0.2fdeg", params_copy.consts.LD, params_copy.consts.sig * (180/pi));
    hold on;
    results.vary_LDsig.qdot1 = max(qdot1);
    results.vary_LDsig.qint1 = max(qint1);
    
    params_copy = params;
    params_copy.consts.LD = params_copy.vary_LD(2); 
    params_copy.consts.sig = params_copy.vary_sigma(2); 
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot2, qint2] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot2, ':k', 'LineWidth', 1.5);
    sent2 = sprintf("L/D = %0.3f, Bank = %0.2fdeg", params_copy.consts.LD, params_copy.consts.sig * (180/pi));
    hold on;
    results.vary_LDsig.qdot2 = max(qdot2);
    results.vary_LDsig.qint2 = max(qint2);
    

    params_copy = params;
    params_copy.consts.LD = params_copy.vary_LD(3); 
    params_copy.consts.sig = params_copy.vary_sigma(3);
    [t, sol] = Solver_EOM_3DOF(params_copy);
    [qdot3, qint3] = Aero_heating_model(t, sol.rho, sol.V, params);
    plot(t, qdot3, '-k', 'LineWidth', 1.5);
    sent3 = sprintf("L/D = %0.3f, Bank = %0.2fdeg", params_copy.consts.LD, params_copy.consts.sig * (180/pi));
    hold on;
    results.vary_LDsig.qdot3 = max(qdot3);
    results.vary_LDsig.qint3 = max(qint3);
    
    ax = gca;
    ax.FontSize = 12; % Set font size for axis labels and ticks
    ax.LineWidth = 1.5; % Set line width for the axis lines

    legend([sent1, sent2, sent3]);
    xlabel('Time (sec)')
    ylabel('Stagnation Point Heating Rate W/cm^2');
    xlim([0 200]);
    
end