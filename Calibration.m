clear
clc
% Number of Monte Carlo repetitions (MAKE SURE TO USE LARGE STEPS)
N = 1000:500:10000; % Maximum of 10,000 trials

% Initialize variables
results = zeros(1, 3); % To store temporary Monte Carlo results
avgvar = zeros(1, length(N)); % To store the average variance for each N
calibrating_trial = NaN; % To store the number of trials at the stopping criterion

% Calibration loop
z = 0; % Counter for indexing
for j = N % Iterate through the range of Monte Carlo repetitions
    z = z + 1;

    % Compute 3 estimates of the variance for each N
    for i = 1:3
        z
        i
        results(i) = TwoGeneratorsIndirect(j); % Call the Monte Carlo function.
    end

    % Compute the mean variance of the 3 estimates
    avgvar(z) = mean(results);

    % Check the stopping criterion
    if avgvar(z) < 5e-6
        calibrating_trial = N(z); % Store the value of N that meets the criterion
        break; % Exit the calibration loop
    end
end

% Plotting
figure;

% Plot the average variance against N
plot(N(1:z), avgvar(1:z), 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;

% Plot the stopping criterion as a horizontal red dashed line
plot(N, 5e-6 * ones(size(N)), 'r--', 'LineWidth', 1.5);

% Highlight the point where the stopping criterion is met
if ~isnan(calibrating_trial)
    scatter(calibrating_trial, avgvar(z), 100, 'r', 'filled', ...
        'DisplayName', 'Stopping Criterion Met');
end

% Add labels, legend, and title
xlabel('Number of Monte Carlo Repetitions (N)', 'FontSize', 12);
ylabel('Average Variance', 'FontSize', 12);
title('Monte Carlo Model Calibration', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Average Variance', 'Stopping Criterion'}, 'Location', 'northeast');
grid on;

% Annotate the stopping criterion
if ~isnan(calibrating_trial)
    text(calibrating_trial, avgvar(z) + 0.5e-6, ...
        sprintf(' N = %d', calibrating_trial), ...
        'FontSize', 10, 'Color', 'r', 'HorizontalAlignment', 'center');
end
