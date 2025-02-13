% Example MATLAB code for plotting three lines in AIAA journal paper format
load Stardust.mat
DataStardust;
% Create a new figure
figure(1);
plot(t1, sol1.h/1000, '-k', 'LineWidth', 1.5); % Line 1
hold on;
plot(t2, sol2.h/1000, '--k', 'LineWidth', 1.5); % Line 2
ylim([0, 120]);
plot(SD_Alt_Time(:, 1), SD_Alt_Time(:, 2), ':k', 'LineWidth', 1.5); % Line 3

ax = gca;
ax.FontSize = 12; % Set font size for axis labels and ticks
ax.LineWidth = 1.5; % Set line width for the axis lines

xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold'); % X-axis label
ylabel('Altitude (km)', 'FontSize', 14, 'FontWeight', 'bold'); % Y-axis label

legend({'Planar EOM', 'Non-Planar EOM', 'Nock et al'}, 'FontSize', 12, 'Location', 'best');

set(gcf, 'Position', [100, 100, 600, 400]); % [left, bottom, width, height]
saveas(gcf, 'Stardust_Alt_Time.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);
plot(t1, sol1.V/1000, '-k', 'LineWidth', 1.5); % Line 1
hold on;
plot(t2, sol2.V/1000, '--k', 'LineWidth', 1.5); % Line 2
plot(SD_Vel_Time(:, 1), SD_Vel_Time(:, 2), ':k', 'LineWidth', 1.5); % Line 3

ax = gca;
ax.FontSize = 12; % Set font size for axis labels and ticks
ax.LineWidth = 1.5; % Set line width for the axis lines

xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold'); % X-axis label
ylabel('Velocity (km/s)', 'FontSize', 14, 'FontWeight', 'bold'); % Y-axis label

legend({'Planar EOM', 'Non-Planar EOM', 'Nock et al'}, 'FontSize', 12, 'Location', 'best');

set(gcf, 'Position', [100, 100, 600, 400]); % [left, bottom, width, height]
saveas(gcf, 'Stardust_Vel_Time.png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
plot(t1, qdot1, '-k', 'LineWidth', 1.5); % Line 1
hold on;
plot(t2, qdot2, '--k', 'LineWidth', 1.5); % Line 2

plot(SD_HR_Time(:, 1), SD_HR_Time(:, 2), ':k', 'LineWidth', 1.5); % Line 3

ax = gca;
ax.FontSize = 12; % Set font size for axis labels and ticks
ax.LineWidth = 1.5; % Set line width for the axis lines

xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold'); % X-axis label
ylabel('Heat Rate (W/cm^2)', 'FontSize', 14, 'FontWeight', 'bold'); % Y-axis label

legend({'Planar EOM', 'Non-Planar EOM', 'Nock et al'}, 'FontSize', 12, 'Location', 'best');

set(gcf, 'Position', [100, 100, 600, 400]); % [left, bottom, width, height]
saveas(gcf, 'Stardust_heating.png');
