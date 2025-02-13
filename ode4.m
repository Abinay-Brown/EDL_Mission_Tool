function [tspan, sol, flag] = ode4(odefun, tspan, y0, h)
    y = y0';

    sol = zeros(length(tspan), length(y0));
    sol(1, :) = y0;
    for i = 1:length(tspan)-1
        k1 = odefun(tspan(i), y);
        k2 = odefun(tspan(i) + h/2, y + h*k1/2);
        k3 = odefun(tspan(i) + h/2, y + h*k2/2);
        k4 = odefun(tspan(i) + h, y + h*k3);
        y = y + h/6 * (k1 + 2*k2 + 2*k3 + k4);
        if any([all(k1 ==0), all(k2 ==0), all(k3==0), all(k4==0)] == 1)
            tspan = tspan(1:i);
            flag = true;
            break;
        else
            flag = false;
        end

        sol(i+1, :) = y;
        % disp(y');
    end
end
