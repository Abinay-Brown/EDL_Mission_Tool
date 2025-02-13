function rho = Atm_1962_1976_model(h, params)
%% Function
%% Inputs:
%% h: Geometric Altitude in (m)
%% Outputs:
%% rho: density(kg/m^3)
    Rbar = params.consts.Rbar;
    r0 = params.consts.r0;
    g0 = params.consts.g0;
    
    Z = [86, 100, 110, 120, 150] * 1000; % Geometric Altitude Ranges (m)

    % Use 1962 atmosphere model if altitude is above 86km    
    if h > Z(1) || h == Z(1)
        % 1962 Standard Atmosphere Data (Geometric Altitude)
        MW = [28.9644, 28.88, 28.56, 28.08] / 1000; % kg/mol
        Pi = [0.34313, 0.030075, 0.0073544, 0.0025217]; % N/m^2
        Ti = [186.945, 210.65, 260.65, 360.65]; % Molecular Scale Temp (K)
        T = [186.946, 210.02, 257.0, 349.49]; % Kinetic Temp (K)
        Li = [1.6481, 5.0, 10.0, 20.0] / 1000; % Molecular Lapse Rate (K/m)
    
        for i = 1:length(Z) - 1
            if (h > Z(i) && h < Z(i + 1)) || h == Z(i)
                break;
            end
        end
        
        b = 3.31 * 10^-7;

        MWval = interp1(Z(1:end-1), MW, h, 'linear', 'extrap');
        Tmval = interp1(Z(1:end-1), T, h, 'pchip', 'extrap');
        
        R = Rbar./MWval;
        
        Tm = (MWval./MW(1)) .* Tmval;
        
        powp = -((g0./R./Li(i)).*(1 + b.*((Ti(i)./Li(i)) - Z(i))));
        P = (Pi(i) .* (Tm ./ Ti(i)).^powp) .* exp((g0.*b./R./Li(i)) .* (h - Z(i)));
        
        
        rhoi = Pi(i) ./ Ti(i) ./ R; 

        powr = -((g0./R./Li(i)).*((R.*Li(i)./g0) + 1 + b.*((Ti(i)./Li(i)) - Z(i))));
        rho = (rhoi .* (Tm ./ Ti(i)).^powr) .* exp((g0.*b./R./Li(i)) .* (h - Z(i)));
        
        T = Tm;
    elseif h < Z(1) 
        % Use 1976 atmosphere model if altitude is below 86km
        Z = [0, 11.0102, 20.0631, 32.1619, 47.3501, 51.4125, 71.802, 86] * 1000;
        hi = [0, 11, 20, 32, 47, 51, 71, 84.852] * 1000;
        Li = [-6.5, 0.0, 1.0, 2.8, 0, -2.8, -2.0] / 1000;
        Ti = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65];
        Pi = [101325, 22631.95, 5474.79, 868.01, 110.9, 66.94, 3.956];

        for i = 1:length(Z) - 1
            if (h > Z(i) && h < Z(i + 1)) || h == Z(i)
                break;
            end
        end
        hg = (r0 .* h)./ (r0 + h);
        R = Rbar ./ (28.9644 ./1000);
        
        T = Ti(i) + Li(i) .* (hg - hi(i));
        if Li(i) == 0
            P = Pi(i) .* exp(-g0 .* (hg - hi(i))./R./Ti(i));
        else
            P = Pi(i) .* (Ti(i)./(Ti(i) + Li(i).*(hg - hi(i)))).^(g0./R./Li(i));
        end
        
        rho = P ./ R ./ T;
    else
        rho = 0;
        
    end
    
end