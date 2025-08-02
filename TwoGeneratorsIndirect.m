%This code develops an INDIRECT Monte Carlo solution to Exercise 1 of the 
% second midterm exam.

% Comment these to calibrate the model
% clear
% clc
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = TwoGeneratorsIndirect (nt) % Uncomment the function to calibrate the model

%Reliability parameters
% All repair rates are per hour
    m_b =  1/60; % Repair rate when two generators have failed
    m   =  1/10; % Repair rate for a single generator
    l_h =  1e-3; % Failure rate of 1 generator under half load 
    l_f =  5*l_h; % Failure rate of 1 generator under full load
    l_c =  3e-4; % Common cause failure rate

    
    %Initialize Monte Carlo simulation parameters
    tm = 4*365*24;       % mission time
    % nt = 8e3;      % number of Monte Carlo trials % comment out to calibrate the model
    dt = 10;    %time resolution
    time_axis = 0:dt:tm;                    % definition of the time axis for the simulation
    counter = zeros(nt, length(time_axis));  % one counter for each simulation where we store the system state
    first_failure_time = inf*ones(1, nt);    % variable to store the time of first system failure for reliability computation
    % counter_tansf_repair = zeros(nt, 1);     % one counter for each simulation where we store the number of transformer repairs 
        % Definition of rate transition matrix 
    A = [0,          0,      m_b;...
           l_f + l_c, 0,      m;...
           l_c,       2*l_h, 0];

    % System transition rate
    LAMBDA =  sum(A, 2);
 
    for i = 1:nt  %the Monte Carlo loop
        % fprintf("Trial %d\n", i)
        % parameter initialization for each trial
        t = 0;
        failure_time = 0;
        repair_time = 0;
        state = 3; % For the sake of the chosen method and MATLAB indexing, this was changed to 
        % represent two generators working at the initial time. States 2 and 1 represent one failed and 
        % both failed transformers respectively.

        flag_first_failure = 0;  %0 if first failure has not yet occurred

        while t < tm            %check if mission time is reached and stop the current trial
        %**********************************************************************    

        % WRITE YOUR INDIRECT MONTE-CARLO ALGORITHM HERE.

        % YOU CAN USE THE LAMBDA MATRIX AS IN THE EXERCISE WITH THE ALARM
        % CLOCK, OR YOU CAN GO FROM STATE TO STATE AND CALCULATE THE
        % TRANSITION RATES TO THE OTHER STATES.

        % BE MINDFUL THAT THE RATE FOR OPERATING UNDER FULL LOAD CHANGES!
            
             % Sampling of the transition times
            t = t - log(1 - rand())/LAMBDA(state);
            
            % Keep track of previous states in order to know which specific transition back into 
            % state 3 (both operational). Useful in calculating transformer repairs
            % prev_state = state; % UNCOMMENT FOR COMPUTING TRANSFORMER REPAIRS

            % Using cumulative probabilities to determine the next state
            state = find(cumsum(A(state, :)/LAMBDA(state)) > rand(), 1, "first");
            % Update the counters
            % Both transformers are repaired if they are in a failed state simultaneously. For this parallel setup, 
            % a system failure is only recorded when both transformers are in a failed state (1).
            if state == 1
                failure_time = t;
                counter(i, time_axis >= failure_time) = 1;
                if flag_first_failure == 0
                    first_failure_time(i) = t;
                    flag_first_failure = 1;
                end
            elseif state == 2
                failure_time = t; % There is a component failure, but the system has not failed
                counter(i, time_axis >= failure_time) = 0;
            else % state == 3
                repair_time = t;
                counter(i, time_axis >= repair_time) = 0;
                %% TRANSFORMER REPAIRS (UNCOMMENT TO COMPUTE)
                % Count two transformer repairs if state transition occurs from
                % both failed (1) to both operational (3), and one transformer
                % repair if state transition occurs from one failed (2) to
                % both operational (3).
                % if prev_state == 1
                %     counter_tansf_repair(i, 1) = counter_tansf_repair(i, 1) + 2;
                % elseif prev_state == 2
                %     counter_tansf_repair(i, 1) = counter_tansf_repair(i, 1) + 1;
                % end
            end
        %*********************************************************************       
        end
      %*********************************************************************  
    end   %end of the Monte Carlo loop

    unrel_counter = zeros(nt, length(time_axis)); %unreliability counter: the system should work continuously so we register the first system failure
        for j = 1:nt
            unrel_counter(j, time_axis >= first_failure_time(1, j)) = 1;
        end

    %  For the sake of direct comparison with the solution obtained from the
    % steady state probabilities, the time-dependent availability and
    % reliability were plotted
    unav = sum(counter, 1)/nt;   % availability computation
    unrel = sum(unrel_counter, 1)/nt;   % reliability computation
    
    av = 1 - unav;   % availability computation
    rel = 1 - unrel;   % reliability computation
    
    % Find the first time index at which reliability reaches zero
    zero_rel_time = find(rel == 0, 1);
    fprintf("\nThe reliability reaches zero after %d hours. Monte Carlo trials: %d\n", zero_rel_time, nt)

    % CALIBRATION CRITERION
    % repairs = sum(counter_tansf_repair)/ nt; % Mean number of transformer repairs
    
    % Calculate the variance of the unavailability at each estimate.
    F = var(unav);
    %%  PLOTTING THE RESULTS (COMMENT OUT THIS SECTION DURING CALIBRATION)
    % figure;
    % plot(time_axis, counter(1, :),'k.', DisplayName="Single Trial")
    % hold on
    % plot(time_axis,av,'b', DisplayName="Monte Carlo Availability")
    % 
    % %%%%%%%%%%%%%%%%%%%%%%
    % % PLOT THE STEADY STATE UNAVAILABILITY FROM QUESTION 1A)
    % %%%%%%%%%%%%%%%%%%%%%%
    % 
    % % The limit availability is obtained from the steady-state probabilities
    % % and plotted on the same graph as the Monte Carlo availability
    % 
    % limit_av = 0.9769*ones(length(time_axis), 1);  
    % 
    % plot(time_axis, limit_av, "r--", LineWidth=0.7, DisplayName="Limit Availability")
    % plot(time_axis,rel,'k', DisplayName="Reliability")
    % xlabel('Time')
    % ylabel('Availability/ Reliability')
    % grid on; %axis([0 max(time_axis) 0 lim_unavail*5]) 
    % legend(Location="best")

% END OF THE CALIBRATION FUNCTION
end



