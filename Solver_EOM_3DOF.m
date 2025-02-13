function [t, res, flag] = Solver_EOM_3DOF(params)
    dt = 0.5;
    tspan = 0:dt:1000;
    
    % options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
    if params.nonplanar == false
        x0 = [params.init_cond.V, params.init_cond.y, params.init_cond.h];
    else
        x0 = [params.init_cond.V, params.init_cond.y, params.init_cond.psi, params.init_cond.r, params.init_cond.theta, params.init_cond.phi];
    end
    func = @(t, y) params.EOM(t, y, params);
    [t, sol, flag] = ode4(func, tspan, x0, dt);
    if ~any(isreal(sol(:,:)))
        flag = true;
    end
    
    if flag == false 
        if params.nonplanar == false
            sol = sol(sol(:,3)>0, :);
            t = t(sol(:, 3)> 0);
            res.V = sol(:, 1);
            res.y = sol(:, 2);
            res.h = sol(:, 3);
            if any(res.h > params.init_cond.h)
                flag = true;
                return;
            end
        else
            Re = params.consts.Req .* (1 - (params.consts.k.*(sin(sol(:, 5)).^2))) * 1000;
            
            res.h = sol(:, 4) - Re;
            
            res.V = sol(:, 1);
            res.y = sol(:, 2);
            res.psi = sol(:, 3);
            res.r = sol(:, 4);
            res.theta = sol(:, 5);
            res.phi = sol(:, 6);
            % if any(res.h > params.init_cond.h)
            %     flag = true;
            %     return;
            % end
        end
    
        rhos = zeros(length(res.h), 1);
      
        for i = 1:length(res.h)
           
            rho_val = params.atm_model(res.h(i), params);
            rhos(i) = rho_val;
        end
        res.rho = rhos;
    else
        if params.nonplanar == false
            res.V = sol(:, 1);
            res.y = sol(:, 2);
            res.h = sol(:, 3);
    
        else
            Re = params.consts.Req .* (1 - (params.consts.k.*(sin(sol(:, 5)).^2))) * 1000;
            
            res.V = sol(:, 1);
            res.y = sol(:, 2);
            res.psi = sol(:, 3);
            res.r = sol(:, 4);
            res.h = sol(:, 4) - Re;
            res.theta = sol(:, 5);
            res.phi = sol(:, 6);
        end
    
    end
    
end