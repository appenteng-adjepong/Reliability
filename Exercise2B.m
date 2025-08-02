clear
clc

Trials = 10:10:1000;
avg_numbers = zeros(length(Trials),4); % array will look like [avgfails(T1) avgfails(T2) avgrepair(T1) avgrepair(T2)]
for i = 1:length(Trials)
    i
    avg_numbers(i,:) = TwoGeneratorsDirect(Trials(i));

end

avg_failure_difference = abs(avg_numbers(:,1)-avg_numbers(:,2));
avg_repair_difference = abs(avg_numbers(:,3)-avg_numbers(:,4));
 
% Plotting
figure;

% Plot the average failure difference for each trial
plot(Trials, avg_failure_difference, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;
% Plot the average repair difference for each trial
plot(Trials, avg_repair_difference, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

%set(gca, 'XScale', 'log', 'YScale', 'log');

% Add labels, legend, and title
xlabel('Number of trials', 'FontSize', 12);
ylabel('Average difference(failure/repair)', 'FontSize', 12);
title('Difference In Failure/Repair Numbers Against Number Of Trials', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Average failure difference', 'Average repair difference'}, 'Location', 'northeast');
grid on;
%% 2B explanation
% Although the failure and repair rates of both generators are identical,
% the randomness of the exponential distribution used to model these processes
% causes one generator to fail and be repaired more frequently in finite simulations. 
% This occurs because the intervals between failures and repairs are sampled independently
% for each generator, leading to statistical variations over a small number of trials or 
% limited simulation time. While the law of large numbers ensures that failure and 
% repair counts converge for both generators over infinite trials, in a small 
% set of runs (e.g., 10 trials), these random fluctuations result in one generator 
% appearing to fail or recover more often. The generator that fails more frequently
% will also be repaired more often, as repairs follow failures. Increasing the number
% of trials would reduce this variation and bring the results closer to equality.
