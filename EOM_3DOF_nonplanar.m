function state_dot = EOM_3DOF_nonplanar(t, state, params)
    V = state(1);
    y = state(2);
    psi = state(3);
    r = state(4);
    phi = state(5);
    theta = state(6);
    

    A = params.consts.A;
    K = params.consts.k;
    m = params.consts.m;
    Req = params.consts.Req;
    Re = Req * (1 - (K*(sin(phi)^2))) * 1000;
    h = r - Re;
    
    % Auxiliary Variables
    w = params.consts.w;
    if isreal(h)
        rho = params.atm_model(h, params);
        g = params.gravity_model(h, phi, params);
    else
        state_dot = [0, 0, 0, 0, 0, 0];
        return;
    end
    % disp(Re);
    if params.mc.flag == true
        if params.mc.rho_uncertainty_model
            if h > 70000
                rho = rho + (params.mc.rho_3sigma * randn);
            elseif h < 70000 && h > 30000
                rho = rho + (params.mc.rho_3sigma * randn * (h - 27500)/(70000-30000));
            else
                rho = rho + (params.mc.rho_3sigma * randn * 0.05);
            end
        else
            rho = rho + (params.mc.rho_3sigma * randn);
        end
    end
    
    sig = params.consts.sig;
    
    CD = params.consts.m/params.consts.beta/params.consts.A;
    CL = params.consts.LD * CD;
    D = 0.5*rho*V*V*A*CD;
    
    if params.consts.LD == 0
        L = 0;
    else
        L = 0.5*rho*V*V*A*CL;
    end
    
    Vdot = ((-D)/m) - (g*sin(y)) + (w*w*r*cos(theta)*(sin(y)*cos(theta)-cos(y)*sin(theta)*sin(psi)));
    ydot = (((L)/m)*cos(sig) - (g - (V*V)/r)*cos(y) + (2*w*V*cos(theta)*cos(psi)) + (w^2) + cos(theta)*(cos(y)*cos(theta) + sin(y)*sin(theta)*sin(psi)))/V;
    psidot = (((L)/m)*(sin(sig)/cos(y)) - (V*V/r*cos(y)*cos(psi)*tan(theta)) + 2*w*V*(tan(y)*cos(theta)*sin(psi)-sin(theta)) - (w*w*r/cos(y)*sin(theta)*cos(theta)*cos(psi)))/V;
    rdot = V*sin(y);
    phidot = (V*cos(y)*cos(psi)/r/cos(theta));
    thetadot = (V*cos(y)*sin(psi)/r);

    state_dot = [Vdot, ydot, psidot, rdot, phidot, thetadot]';
end