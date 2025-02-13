function state_dot = EOM_3DOF_planar(t, state, params)
    % Load State Variables 
    V = state(1);   % Reentry Velocity (meters/second)
    y = state(2);   % Flight Path Angle (radians)
    h = state(3);   % Altitude (meters)
    
    % Load Parameters
    beta   = params.consts.beta;  % Ballistic Coefficient
    LD = params.consts.LD;
    sig = params.consts.sig;
    
    % Determine Auxiliary Variables
    if isreal(h)    
        rho = params.atm_model(h, params);
        g = params.gravity_model(h, params);
        r = params.consts.Re + h;
        else
        
        state_dot = [0, 0, 0, 0, 0, 0];
        return;
    end
    % Equations of Motion
    Vdot = -((rho.*V.^2)/(2.*beta)) + g.*sin(y);
    ydot = ((-(V.^2).*cos(y)./r) - (rho*V*V/2/beta*LD*cos(sig))+ g.*cos(y))./V;
    hdot =  - V.*sin(y);
    
    state_dot = [Vdot; ydot; hdot];
end
