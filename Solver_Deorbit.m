function [sol] = Solver_Deorbit(params)
    hp = params.deorbit.hp; % km
    e = params.deorbit.e;
    thetad = params.deorbit.thetad; % radians
    nmax_constraint = params.deorbit.nmax_constraint;
    g0 = params.consts.g0;
    Re = params.consts.Re/1000;
    Ratm = Re + params.consts.hinterface; % km
    H = params.consts.H;
    
    mu = params.consts.mu; 
    rp = Re + hp;
    a1 = rp / (1 - e);
    ra = a1 * (1 + e);
    
    eps = -mu/(2*a1);
    p = a1*(1-e^2);
    h = sqrt(mu*p);

    rd = p / (1 + e*cos(thetad));
    lam = rd/Ratm;
    alpha1 = a1/Ratm;
    u1 = sqrt((2*alpha1 - lam)/alpha1);
    v1 = u1 * sqrt(mu/rd);
    y1 = acos(alpha1*sqrt((1-e^2)/(lam*(2*alpha1-lam))));
    ue_constraint = (2*lam*(lam-1)*(1-e^2)*alpha1^2)/((lam*(1-e^2)*alpha1^2)-2*alpha1 + lam);
    ue_constraint = sqrt(ue_constraint);
    Ve = ue_constraint * sqrt(mu/rd);
    Ve = 0.9999998*Ve;
    ue = Ve / sqrt(mu/rd);
    Vatm = Ve;
    
    y = asin(H*1000 * (2*g0*exp(1)) * nmax_constraint / ((Vatm*1000)^2));
    

    %% solving for y upper limit (bisection method)
    
    ya = y;
    yb = pi/2;
    yc = (ya +yb)/2;
    tol = 0.001;
    nmaxc = 2000;
    while true
        if nmaxc > nmax_constraint
            yb = yc;
            yc = (ya + yb)/2;
        end
        if nmaxc < nmax_constraint
            ya = yc;
            yb = pi/2;
            yc = (ya + yb)/2;
        end

        params_copy = params;
        params_copy.init_cond.V = Vatm*1000;
        params_copy.init_cond.y =  ya;
        params_copy.init_cond.h = params.consts.hinterface*1000;
        [t, sol, flaga] = Solver_EOM_3DOF(params_copy);
        dvdt = -(sol.rho.*sol.V.^2)./(2*params_copy.consts.beta) + g0*sin(sol.y);
        nmaxa = max(abs(dvdt/g0));

        params_copy = params;
        params_copy.init_cond.V = Vatm*1000;
        params_copy.init_cond.y =  yb;
        params_copy.init_cond.h = params.consts.hinterface*1000;
        [t, sol, flagb] = Solver_EOM_3DOF(params_copy);
        dvdt = -(sol.rho.*sol.V.^2)./(2*params_copy.consts.beta) + g0*sin(sol.y);
        nmaxb = max(abs(dvdt/g0));

        params_copy = params;
        params_copy.init_cond.V = Vatm*1000;
        params_copy.init_cond.y =  yc;
        params_copy.init_cond.h = params.consts.hinterface*1000;
        [t, sol, flagc] = Solver_EOM_3DOF(params_copy);
        dvdt = -(sol.rho.*sol.V.^2)./(2*params_copy.consts.beta) + g0*sin(sol.y);
        nmaxc = max(abs(dvdt/g0));

        if abs(nmaxc- nmax_constraint) < tol
            
            ymax = yc;
            break;
        end
        
    end

    %% solving for y lower limit 
    ya = 0;
    yb = ymax;
    yc = (ya + yb)/2;
    tol = 0.001;
    nmaxc = 2000;
    flaga = false;
    flagb = false;
    flagc = false;
    nmaxc_prev = 0;
    while true

        params_copy = params;
        params_copy.init_cond.V = Vatm*1000;
        params_copy.init_cond.y =  ya;
        params_copy.init_cond.h = params.consts.hinterface*1000;
        [t, sol, flaga] = Solver_EOM_3DOF(params_copy);
        if flaga == false
            dvdt = -(sol.rho.*sol.V.^2)./(2*params_copy.consts.beta) + g0*sin(sol.y);
            nmaxa = max(abs(dvdt/g0));
        end
        params_copy = params;
        params_copy.init_cond.V = Vatm*1000;
        params_copy.init_cond.y =  yb;
        params_copy.init_cond.h = params.consts.hinterface*1000;
        [t, sol, flagb] = Solver_EOM_3DOF(params_copy);
        if flagb == false
            dvdt = -(sol.rho.*sol.V.^2)./(2*params_copy.consts.beta) + g0*sin(sol.y);
            nmaxb = max(abs(dvdt/g0));
        end
        params_copy = params;
        params_copy.init_cond.V = Vatm*1000;
        params_copy.init_cond.y =  yc;
        params_copy.init_cond.h = params.consts.hinterface*1000;
        [t, sol, flagc] = Solver_EOM_3DOF(params_copy);
        if flagc == false
            dvdt = -(sol.rho.*sol.V.^2)./(2*params_copy.consts.beta) + g0*sin(sol.y);
            nmaxc = max(abs(dvdt/g0));
        end
        if flaga == true && flagc == false
            yb = yc;
            yc = (ya + yb)/2;
        end

        if flaga == true && flagc == true
            ya = yc;
            yc = (ya + yb)/2;
        end

        if flagc == false
            if abs(nmaxc - nmaxc_prev) < tol
                ymin = yc;
                break;
            end
            nmaxc_prev = nmaxc;
            
        end    

    end

    cos_ye = (alpha1/ue)*sqrt(lam*(1-e^2)*(ue*ue + 2*(1-lam))/(2*alpha1 - lam));
    ye = acos(cos_ye);

    % disp(ye);

    %% Bisection for entry velocity
    ymax = 9.26 * (pi/180);
    ymin = 4.864 * (pi/180);
    
    Ve = ue_constraint * sqrt(mu/rd);
    Veb = 0.999999998*Ve;
    
    ueb = Veb /sqrt(mu/rd);
    cos_yeb = (alpha1/ueb)*sqrt(lam*(1-e^2)*(ueb*ueb + 2*(1-lam))/(2*alpha1 - lam));
    yeb = acos(cos_yeb);
    
    if ((yeb > ymin) && (yeb < ymax))
        
    else
        while true
        
            ueb = Veb /sqrt(mu/rd);
            cos_yeb = (alpha1/ueb)*sqrt(lam*(1-e^2)*(ueb*ueb + 2*(1-lam))/(2*alpha1 - lam));
            yeb = acos(cos_yeb);
            
            
            if ((yeb > (ymin+ymax)/2) && (yeb < ymax))
                
                break;
            else
                Veb = Veb - 10^-7;
            
            end

        
        end
    u2 = sqrt((ueb^2) + 2*(1-lam));
    deltau = u1 - u2;
    deltaV = deltau*sqrt(mu/rd);
    sol.deltaV = deltaV;
    sol.Ve = ueb * sqrt(mu/rd);
    sol.ymin = ymin * (180/pi);
    sol.ye = yeb * (180/pi);
    sol.ymax = ymax * (180/pi);

end