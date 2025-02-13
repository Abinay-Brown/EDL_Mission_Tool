function [lat_samp, lon_samp] = Solver_Monte_Carlo(params)
lat_samp = [];
lon_samp = [];
disp("Running Monte-Carlo Simulation, Please Wait...");
for i = 1: params.mc.samples
    params_copy = params;
    params_copy.init_cond.V = params_copy.init_cond.V + (params_copy.mc.V_3sigma*randn);
    params_copy.init_cond.y = params_copy.init_cond.y + (params_copy.mc.y_3sigma*randn);
    params_copy.init_cond.phi = params_copy.init_cond.phi + (params_copy.mc.lat_3sigma*randn);
    params_copy.init_cond.theta = params_copy.init_cond.theta + (params_copy.mc.lon_3sigma*randn);
   
    [t, res, flag] = Solver_EOM_3DOF(params_copy);
    if flag == false && ~isnan(res.phi(end)) && ~isnan(res.theta(end))
    lat_samp = [lat_samp, res.phi(end) * (180/pi)];
    lon_samp = [lon_samp, res.theta(end) * (180/pi)];
    end
    if flag == true
        flag == false;
    end
    
end

meanLat = mean(lat_samp);
meanLon = mean(lon_samp);
latCentered = lat_samp - meanLat;
lonCentered = lon_samp - meanLon;

% Covariance matrix
covMatrix = cov(latCentered, lonCentered);
[eigVec, eigVal] = eig(covMatrix);

eigValues = diag(eigVal);
[majorAxisLength, idx] = max(eigValues);
majorAxisLength = 3*majorAxisLength;
minorAxisLength = min(eigValues);
minorAxisLength = 3* minorAxisLength;
majorAxisVector = eigVec(:, idx);

% Calculate the angle of the major axis
theta = atan2(majorAxisVector(2), majorAxisVector(1));

% Generate ellipse points
t = linspace(0, 2*pi, 100);
ellipse_x = sqrt(majorAxisLength) * cos(t);
ellipse_y = sqrt(minorAxisLength) * sin(t);
% Rotate the ellipse
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
ellipse_points = R * [ellipse_x; ellipse_y];

% Shift ellipse to the mean location
ellipseLat = ellipse_points(1, :) + meanLat;
ellipseLon = ellipse_points(2, :) + meanLon;
figure;
scatter(ellipseLon, ellipseLat, '.k');
hold on;
scatter(lon_samp, lat_samp, 'ok');

ax = gca;
ax.FontSize = 12; % Set font size for axis labels and ticks
ax.LineWidth = 1.5; % Set line width for the axis lines
set(gcf, 'Position', [100, 100, 600, 400]); % [left, bottom, width, height]
% saveas(gcf, 'Stardust_Vel_Time.png');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
legend('3-sigma Landing Ellipse', 'Landing Sites');
end