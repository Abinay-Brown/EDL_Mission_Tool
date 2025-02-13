% Example MATLAB code for plotting three lines in AIAA journal paper format
load Mars.mat
DataViking;
% Create a new figure
figure(1);
plot(sol1.V/1000, sol1.h/1000, '-k', 'LineWidth', 1.5); % Line 1
hold on;
plot(sol2.V/1000, sol2.h/1000, '--k', 'LineWidth', 1.5); % Line 2
ylim([0, 120]);
plot(Viking_sol(:, 1), Viking_sol(:, 2), ':k', 'LineWidth', 1.5); % Line 3

ax = gca;
ax.FontSize = 12; % Set font size for axis labels and ticks
ax.LineWidth = 1.5; % Set line width for the axis lines

xlabel('Velocity (km/s)', 'FontSize', 14, 'FontWeight', 'bold'); % X-axis label
ylabel('Altitude (km)', 'FontSize', 14, 'FontWeight', 'bold'); % Y-axis label

legend({'Planar EOM', 'Non-Planar EOM', 'Rasky & Subrahmanyam'}, 'FontSize', 12, 'Location', 'best');

set(gcf, 'Position', [100, 100, 600, 400]); % [left, bottom, width, height]
saveas(gcf, 'Viking_trajectory.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);
plot(t1, qdot1, '-k', 'LineWidth', 1.5); % Line 1
hold on;
plot(t2, qdot2, '--k', 'LineWidth', 1.5); % Line 2
ylim([0, 120]);
plot(Viking_qdot(:, 1), Viking_qdot(:, 2), ':k', 'LineWidth', 1.5); % Line 3
xlim([0 300])
ylim([0, 30])
ax = gca;
ax.FontSize = 12; % Set font size for axis labels and ticks
ax.LineWidth = 1.5; % Set line width for the axis lines

xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold'); % X-axis label
ylabel('Stagnation Point Heat Rate (W/cm^2)', 'FontSize', 14, 'FontWeight', 'bold'); % Y-axis label

legend({'Planar EOM', 'Non-Planar EOM', 'Rasky & Subrahmanyam'}, 'FontSize', 12, 'Location', 'best');

set(gcf, 'Position', [100, 100, 600, 400]); % [left, bottom, width, height]
saveas(gcf, 'Viking_heating.png');
