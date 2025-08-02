%This code develops a DIRECT Monte Carlo solution to Exercise 2 of the 
% second midterm exam.

% Comment these to calibrate the model
% clear
% clc
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = TwoGeneratorsDirect(nt) % Uncomment the function to calibrate the model

%% Reliability parameters (all rates are per hour)
    m_b =  1/60; % Repair rate when two generators have failed
    m   =  1/10; % Repair rate for a single generator
    l_h =  1e-3; % Failure rate of 1 generator under half load 
    l_f =  5*l_h; % Failure rate of 1 generator under full load
    l_c =  3e-4; % Common cause failure rate

    %% Initialize Monte Carlo simulation parameters
   tm = 24 * 365 * 4;      % mission time
    % nt = 8e3;              % number of Monte Carlo trials % comment if used as a function
    dt = 10;                % time resolution
    time_axis = 0:dt:tm;    % definition of the time axis for the simulation
    counter = zeros(nt, length(time_axis)); % counter for system state
    first_failure_time = inf * ones(1, nt); % time of first system failure
    unrel_counter = zeros(nt, length(time_axis));
    generator_failure_count = zeros(nt,2); %First index for transformer 1, second index for transformer 2 
    generator_repair_count = zeros(nt,2);
 %% 
 
for i = 1:nt  %the Monte Carlo loop
    % fprintf("Trial %d\n", i)
    %parameter initialization for each trial
    t = 0;
    failure_time = 0;
    repair_time = 0;
    state = [1, 1];     %two generators work at time 
    old_state = [1, 1]; % MAKE SURE TO ALSO STORE THE OLD STATES OF THE GENERATORS

    flag_first_failure = 0;  %0 if first failure has not yet occurred

    while t < tm            %check if mission time is reached and stop the current trial
    %**********************************************************************    
        %% Time sampling
        % DEVELOP THE DIRECT ALGORITHM TO SOLVE THE PROBLEM.

        % MAKE SURE THAT YOU SAMPLE INDIVIDUALLY FOR EACH GENERATOR, AND
        % THEN UNDERSTAND THE STATE CHANGE OF THE SYSTEM.

        % WHEN ONE GENERATOR FAILS AND THE OTHER OPERATES UNDER FULL LOAD,
        % YOU HAVE TO RE-ADJUST THE TIME TO FAILURE WITH THE FOLLOWING
        % EQUATION: timeNew = l_new / l_old * (time_difference);

        % Where l_new is l_f if the other generator failed, or l_h if the
        % other generator was repaired.

        % WHEN A GENERATOR IS REPAIRED, YOU HAVE TO RESAMPLE THE FAILURE
        % TIME AGAIN!

            if sum(state) == 0 % If both transformers have failed
                generator_failure_count(i,:) = generator_failure_count(i,:) +  ones(1,2);
                time_samples=[inf,inf,exprnd(1/(m_b))]; % Sampling when the repair occurs
                failure_time=t;
                first_index= find(time_axis>=failure_time,1);  %identify the FIRST counter cell to be filled with ones
                counter(i,first_index)=1;          %1 is stored in the counter at the proper time bin, if the system is failed
                    if flag_first_failure==0       %first failure time for 
                        first_failure_time(i)=t;   %system unreliability
                        flag_first_failure=1;
                    end
            elseif state(1)==0 % If only transformer 1 has failed
                generator_failure_count(i,:) = generator_failure_count(i,:) +  [1 0];
                time_samples(2) = l_h/l_f * (time_samples(2) - time_samples(1)); % Sampling the new time until transformer 2 fails due to the increased rate l_f
                time_samples(3)=time_samples(3)-time_samples(1); % New time to state change
                time_samples(1)=exprnd(1/(m)); % Sampling the repair time for transformer 1
             elseif state(2)==0 % If only transformer 2 has failed
                generator_failure_count(i,:) = generator_failure_count(i,:) +  [0 1];
                time_samples(1)= l_h/l_f * (time_samples(1) - time_samples(2)); % Sampling the new time until transformer 1 fails due to the increased rate l_f
                time_samples(3)=time_samples(3)-time_samples(2); % New time to state change
                time_samples(2)=exprnd(1/(m)); % Sampling the repair time for transformer 2
            else % Both transformers are operating
                if old_state(1)==0 && old_state(2)==1 % If transformer 1 was repaired in the previous state change
                    generator_repair_count(i,:) = generator_repair_count(i,:) + [1 0];
                    time_samples(2)=l_f/l_h * (time_samples(2) - time_samples(1));  % Sampling the new time until transformer 2 fails due to the decreased rate l_h
                    time_samples(3)=time_samples(3)-time_samples(1); % New time to state change
                    time_samples(1)=exprnd(1/(l_h)); % Sampling the new failure rate for transformer 1
                elseif old_state(1)==1 && old_state(2)==0 % If transformer 2 was repaired in the previous state change
                    generator_repair_count(i,:) = generator_repair_count(i,:) + [0 1];
                    time_samples(1)=l_f/l_h * (time_samples(1) - time_samples(2)); % Sampling the new time until transformer 1 fails due to the decreased rate l_h
                    time_samples(3)=time_samples(3)-time_samples(2); % New time to state change
                    time_samples(2)=exprnd(1/(l_h)); % Sampling the new failure rate for
                else % If both transformers were recovered in the previous state change, or it is the first iteration of each trial
                     time_samples=[exprnd(1/(l_h)),exprnd(1/(l_h)),exprnd(1/(l_c))]; % Sampling the time for state change
                    if old_state(1)==0 && old_state(2)==0 % If both transformers were recovered in the previous state change
                        generator_repair_count(i,:) = generator_repair_count(i,:) + [1 1];
                        repair_time=t; % Storing the repair time for the system
                        last_index =find(time_axis<repair_time,1,'last');  % identify the LAST counter cell to be filled with ones
                        counter(i,first_index+1:last_index)=1;    % 1 is stored in the counter cells in the proper time bins, if the system is failed                
                    end
                end
            end
    
    % Updating the time and finding the state 
       old_state = state;
       
       [time_update, newstate] = min(time_samples); % Finding the first state change: [change1, change2, change3] = [G1change, G2change, BothGeneratorsChange]
                              % HINT: Use the min() function as it provides both the first time of change and the new state.    
                              
    t = t + time_update; % Progressing the time
        
    if newstate==1 % If transformer 1 has a state change
        if state(1)==0
            state(1)=1;
        else
            state(1)=0;
        end
    elseif newstate==2 % If transformer 2 has a state change
        if state(2)==0
            state(2)=1;
        else
            state(2)=0;
        end
    else % If both transformers have a state change
        if state(1)==0 && state(2)==0
            state(1)=1;
            state(2)=1;
        else
            state(1)=0;
            state(2)=0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UPDATE THE STATE CHANGE FOR THE GENERATORS AND THE SYSTEM
    % MAKE SURE TO ALSO UPDATE THE PREVIOUS STATES IN ORDER TO KNOW WHICH
    % GENERATOR WAS REPAIRED.
   
    end  %while loop wrt time; end of each MC trial
    unrel_counter(i, time_axis >= first_failure_time(i)) = 1;  %unreliability counter: the system should work continuously so we register the first system failure

end   %end of the Monte Carlo loop                                                      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjusting the counters and plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unav = sum(counter, 1)/nt;          % unavailability computation
unrel = sum(unrel_counter, 1)/nt;   % unreliability computation
av = 1 - unav;
rel = 1 - unrel;
total_generator_repair = sum(generator_repair_count,1);
total_generator_failure = sum(generator_failure_count,1);

% CALIBRATION CRITERION
% F = var(av); % Calculate the variance of the unavailability at each estimate.
% this should be commented out to run Exercise2B

% Exercise 2B
B = total_generator_failure/nt; % Uncomment for Exercise 2B
B_2 = total_generator_repair/nt;
F = [B B_2]; %Uncomment for Exercise2B

zero_rel_time = find(rel == 0, 1);
fprintf("\nThe reliability reaches zero after %d hours. Monte Carlo trials: %d\n", zero_rel_time, nt)

% PLOTTING THE RESULTS -  Comment this out when calibrating the model!!!
% figure;
% plot(time_axis, counter(1, :),'k.', DisplayName="Single Trial")
% hold on
% plot(time_axis, av,'b', DisplayName="Monte Carlo Availability")

    %%%%%%%%%%%%%%%%%%%%%%
    % PLOT THE STEADY STATE UNAVAILABILITY FROM QUESTION 1A) AND FROM \
    % THE DIRECT SIMULATION!
    %%%%%%%%%%%%%%%%%%%%%%

% limit_av = 0.9769*ones(length(time_axis), 1);  
% plot(time_axis, limit_av, "r--", LineWidth=0.7, DisplayName="Limit Availability")
% plot(time_axis, rel,'k', DisplayName="Reliability")
% xlabel('Time')
% ylabel('Availability/ Reliability')
% grid on; %axis([0 max(time_axis) 0 lim_unavail*5]) 
% legend(Location="best")
% FUNCTION END
end

