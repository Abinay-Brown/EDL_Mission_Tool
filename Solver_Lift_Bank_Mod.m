function [results] = Solver_Lift_Bank_Mod(params)
    
    %% Case1
    params_copy = params;
    params_copy.consts.alpha = params.LB_mod.alpha1;
    params_copy.consts.sig = params.LB_mod.sigma1;

    [CL, CD, A] = params.geom_model(params_copy);
    params_copy.consts.LD = CL/CD;
    params_copy.consts.beta = params_copy.consts.m/CD/A;
    params_copy.consts.A = A;
    
    [t, sol] = Solver_EOM_3DOF(params_copy);
    acc = gradient(sol.V, t(2) - t(1));
    [nmax, idx] = max(abs(acc./params.consts.g0));
    hnmax = sol.h(idx);
    Vnmax = sol.V(idx);
    ds_dh = 1./sin(sol.y);
    s = cumtrapz(sol.h, ds_dh);
    dS_ds = cos(sol.y);
    S = cumtrapz(s, dS_ds); % Downrange (m)
    dl_dS = sin(sol.psi);
    l = cumtrapz(S, dl_dS); % Crossrange (m)
    results.case1.S = S(end)/1000;
    results.case1.l = l(end)/1000;
    results.case1.nmax = nmax;
    results.case1.hnmax = hnmax/1000;
    results.case1.Vnmax = Vnmax/1000;

    %% Case 2
    params_copy = params;
    params_copy.consts.alpha = params.LB_mod.alpha2;
    params_copy.consts.sig = params.LB_mod.sigma2;

    [CL, CD, A] = params.geom_model(params_copy);
    params_copy.consts.LD = CL/CD;
    params_copy.consts.beta = params_copy.consts.m/CD/A;
    params_copy.consts.A = A;

    [t, sol] = Solver_EOM_3DOF(params_copy);
    acc = gradient(sol.V, t(2) - t(1));
    [nmax, idx] = max(abs(acc./params.consts.g0));
    hnmax = sol.h(idx);
    Vnmax = sol.V(idx);
    ds_dh = 1./sin(sol.y);
    s = cumtrapz(sol.h, ds_dh);
    dS_ds = cos(sol.y);
    S = cumtrapz(s, dS_ds); % Downrange (m)
    dl_dS = sin(sol.psi);
    l = cumtrapz(S, dl_dS); % Crossrange (m)
    results.case2.S = S(end)/1000;
    results.case2.l = l(end)/1000;
    results.case2.nmax = nmax;
    results.case2.hnmax = hnmax/1000;
    results.case2.Vnmax = Vnmax/1000;

    %% Case 3
    params_copy = params;
    params_copy.consts.alpha = params.LB_mod.alpha3;
    params_copy.consts.sig = params.LB_mod.sigma3;

    [CL, CD, A] = params.geom_model(params_copy);
    params_copy.consts.LD = CL/CD;
    params_copy.consts.beta = params_copy.consts.m/CD/A;
    params_copy.consts.A = A;

    [t, sol] = Solver_EOM_3DOF(params_copy);
    acc = gradient(sol.V, t(2) - t(1));
    [nmax, idx] = max(abs(acc./params.consts.g0));
    hnmax = sol.h(idx);
    Vnmax = sol.V(idx);
    ds_dh = 1./sin(sol.y);
    s = cumtrapz(sol.h, ds_dh);
    dS_ds = cos(sol.y);
    S = cumtrapz(s, dS_ds); % Downrange (m)
    dl_dS = sin(sol.psi);
    l = cumtrapz(S, dl_dS); % Crossrange (m)
    results.case3.S = S(end)/1000;
    results.case3.l = l(end)/1000;
    results.case3.nmax = nmax;
    results.case3.hnmax = hnmax/1000;
    results.case3.Vnmax = Vnmax/1000;

end