clear; clc; close all;

% === STEP 1: SYSTEM PARAMETERS ===
% These values define the physics of the specific car we are modeling.
ms = 250;     % Sprung Mass (Car Body) = 250 kg
mu = 45;      % Unsprung Mass (Wheel) = 45 kg
ks = 16000;   % Suspension Spring Stiffness = 16,000 N/m
cs = 1000;    % Suspension Damping = 1,000 N-s/m
kt = 190000;  % Tire Stiffness = 190,000 N/m

disp('Step 1 Complete: System Parameters Defined.');

% === STEP 2: STATE-SPACE MATRICES ===
% We convert the physical equations into the form: x_dot = Ax + Bu
A = [0,       1,      0,      -1;
    -ks/ms, -cs/ms,   0,       cs/ms;
     0,       0,      0,       1;
     ks/mu,  cs/mu, -kt/mu,   -cs/mu];

B = [0; 0; -1; 0];
C = [-ks/ms, -cs/ms, 0, cs/ms];
D = 0;

disp('Step 2 Complete: Matrices A, B, C, D built.');

% === STEP 3: OBSERVABILITY CHECK ===
Ob = obsv(A, C);      
rank_Ob = rank(Ob);   

if rank_Ob == 4
    disp(['Step 3 Success: System is Fully Observable (Rank = ' num2str(rank_Ob) ')']);
else
    error('Step 3 Failed: System is NOT Observable. Change the C matrix.');
end

% === STEP 4: DISCRETIZATION & ROAD GENERATION ===
dt = 0.01; 
sys_cont = ss(A, B, C, D);     
sys_disc = c2d(sys_cont, dt);  
[Ad, Bd, Cd, Dd] = ssdata(sys_disc); 

T_end = 5;             
t = 0:dt:T_end;        
N = length(t);         

% Create a "Bump" in the road
road_velocity = zeros(1, N);
road_velocity(t > 1 & t < 1.5) = 0.1; 

disp('Step 4 Complete: System Discretized and Road Profile Created.');

% === STEP 5: THE KALMAN FILTER SIMULATION LOOP ===

% Initialization
x_true = zeros(4, N);      
y_meas = zeros(1, N);      
x_est  = zeros(4, N);      
P_cov  = eye(4);           

x_true(:,1) = [0; 0; 0; 0]; 
x_est(:,1)  = [0; 0; 0; 0]; 

% Tuning
R_actual = 0.001;  % Reality (Low noise)
R_filter = 10;     % Filter Setting (High trust in model -> Smoothness)
Q = 0.0001 * eye(4); 

fprintf('Starting Simulation... ');

for k = 1 : N - 1
    
    % --- 1. SIMULATE THE REAL WORLD ---
    x_true(:,k+1) = Ad * x_true(:,k) + Bd * road_velocity(k);
    
    % Generate Noise 
    noise = sqrt(R_actual) * randn; 
    y_meas(k) = Cd * x_true(:,k) + noise;
    
    % --- 2. RUN THE KALMAN FILTER ---
    x_pred = Ad * x_est(:,k) + Bd * road_velocity(k);
    P_pred = Ad * P_cov * Ad' + Q;
    
    y_guess = Cd * x_pred;
    measurement_residual = y_meas(k) - y_guess;
    
    S = Cd * P_pred * Cd' + R_filter; 
    K = P_pred * Cd' / S;
    
    x_est(:,k+1) = x_pred + K * measurement_residual;
    P_cov = (eye(4) - K * Cd) * P_pred;
    
end
disp('Simulation Complete.');

% === STEP 6: PLOTTING THE RESULTS ===

time_axis = 0:dt:T_end;

% GRAPH 1: The Input
figure(1);
plot(time_axis, road_velocity, 'k', 'LineWidth', 2);
title('Input: Road Velocity Profile');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
ylim([0 0.2]);
grid on;

% GRAPH 2: The Filter Performance
figure(2); clf; hold on;
plot(time_axis(1:end-1), x_true(1, 1:end-1), 'g', 'LineWidth', 3);
plot(time_axis(1:end-1), x_est(1, 1:end-1), 'r--', 'LineWidth', 1.5);
title('Kalman Filter Results: Suspension Deflection');
legend('True Physics (Hidden)', 'Kalman Estimate');
xlabel('Time (s)'); ylabel('Deflection (meters)');
grid on; hold off;

% GRAPH 3: Noise Rejection
figure(3); clf; hold on;
plot(time_axis(1:end-1), y_meas(1:end-1), 'Color', [0.7 0.7 0.7]); 
y_est = Cd * x_est(:, 1:end-1); 
plot(time_axis(1:end-1), y_est, 'b', 'LineWidth', 2);
title('Noise Rejection: Raw Sensor vs. Filtered Output');
legend('Noisy Sensor Data', 'Cleaned Filter Output');
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
grid on; hold off;

% === STEP 7: SENSITIVITY ANALYSIS (TUNING R) ===
fprintf('\nStarting Sensitivity Analysis...\n');

R_values = [0.001, 10, 1000]; 
figure(4); clf; hold on;
colors = {'r', 'b', 'm'}; 
names = {'Low R (Trust Sensor)', 'Med R (Balanced)', 'High R (Trust Model)'};

plot(time_axis(1:end-1), x_true(1, 1:end-1), 'g', 'LineWidth', 4);

for s = 1:3
    R_test = R_values(s);
    x_est_sens = zeros(4, N); 
    P_cov_sens = eye(4);
    x_est_sens(:,1) = [0; 0; 0; 0];
    
    for k = 1 : N - 1
        x_pred = Ad * x_est_sens(:,k) + Bd * road_velocity(k);
        P_pred = Ad * P_cov_sens * Ad' + Q;
        S = Cd * P_pred * Cd' + R_test;
        K = P_pred * Cd' / S;
        y_guess = Cd * x_pred;
        x_est_sens(:,k+1) = x_pred + K * (y_meas(k) - y_guess);
        P_cov_sens = (eye(4) - K * Cd) * P_pred;
    end
    
    plot(time_axis(1:end-1), x_est_sens(1, 1:end-1), ...
         'Color', colors{s}, 'LineWidth', 1.5, 'LineStyle', '--');
end

title('Sensitivity Analysis: Tuning the R Matrix');
legend('Ground Truth', names{1}, names{2}, names{3});
xlabel('Time (s)'); ylabel('Suspension Deflection (m)');
grid on; hold off;

% =========================================================
% === STEP 8: LITERATURE VALIDATION & METRICS ===
% =========================================================

fprintf('\n--- VALIDATION METRICS ---\n');

% 1. CALCULATE RMSE (The "Score")
error_raw = x_true(1,:) - (y_meas / (-ks/ms)); % Approx raw deflection from accel
error_kf  = x_true(1,:) - x_est(1,:);

rmse_raw = sqrt(mean(error_raw.^2));
rmse_kf  = sqrt(mean(error_kf.^2));

fprintf('RMSE Raw Sensor (Approx):  %.5f m\n', rmse_raw);
fprintf('RMSE Kalman Filter:        %.5f m\n', rmse_kf);
fprintf('Improvement Factor:        %.2fx better\n', rmse_raw / rmse_kf);


% 2. IMPLEMENT A "COMPETITOR" (Standard Low Pass Filter)
alpha = 0.1; 
x_lpf = zeros(1, N);
for i = 2:N
    x_lpf(i) = (1-alpha)*x_lpf(i-1) + alpha*x_true(1,i) + 0.005*randn; 
end

% 3. PLOT: KALMAN VS. STANDARD FILTER
figure(6); clf;
hold on;
zoom_idx = find(t > 0.9 & t < 2.5); 

plot(t(zoom_idx), x_true(1, zoom_idx), 'g', 'LineWidth', 3);
plot(t(zoom_idx), x_lpf(zoom_idx), 'b--', 'LineWidth', 1.5);
plot(t(zoom_idx), x_est(1, zoom_idx), 'r', 'LineWidth', 2);

title('Literature Validation: Kalman Filter vs. Standard LPF');
legend('Ground Truth', 'Standard Low Pass (Laggy)', 'Kalman Filter (Fast)');
xlabel('Time (s)'); ylabel('Deflection (m)');
grid on; hold off;

% 4. RESIDUAL ANALYSIS (The "Optimality" Test)
y_predicted = Cd * x_est;
residuals = y_meas - y_predicted;

figure(7); clf;
subplot(2,1,1);
plot(t, residuals, 'k', 'LineWidth', 0.5);
title('Optimality Check: Residuals (Should be White Noise)');
xlabel('Time (s)'); ylabel('Error (m/s^2)');
grid on;

subplot(2,1,2);
histogram(residuals, 50, 'FaceColor', 'b');
title('Histogram of Residuals (Should be Gaussian Bell Curve)');
grid on;

% =========================================================
% === STEP 9: TEST ON KAGGLE DATASET (PVS DATASET) ===
% =========================================================
fprintf('\n--- Processing Kaggle PVS Dataset ---\n');

try
    % 1. IMPORT
    % The PVS dataset has a header row, so we use 'readtable'
    
    data_table = readtable('kaggle_data.csv');
    raw_data = table2array(data_table);

    % 2. EXTRACT VERTICAL ACCELERATION
    % In this specific file, the Z-axis is Column 3.
    y_kaggle_raw = raw_data(:, 3)'; 

    % 3. PRE-PROCESS
    % Remove Gravity (Centering the data around 0)
    y_kaggle = y_kaggle_raw - mean(y_kaggle_raw);
    
    % Unit Conversion: If the data is small (< 2), it's in 'g'. Convert to m/s^2.
    if max(abs(y_kaggle)) < 4
        y_kaggle = y_kaggle * 9.81; 
    end

    % 4. RESAMPLE (Crucial for large files)
    % We only want 5 seconds of simulation (500 steps)
    
    N_sim = 500; 
    
    % Create a time vector for the raw data
    % We assume the raw data was recorded at roughly 100Hz
    t_kaggle = linspace(0, 5, length(y_kaggle(1:5000))); 
    
    % We take a slice of 5000 points to represent 5 seconds
    y_slice = y_kaggle(1:5000); 
    
    t_sim = 0:dt:(N_sim-1)*dt; 
    
    % Interpolate to get exactly 500 clean data points for our model
    y_meas_real = interp1(t_kaggle, y_slice, t_sim, 'linear');

    % 5. RUN FILTER
    x_est_real = zeros(4, N_sim); 
    P = eye(4);
    
    for k = 1:N_sim-1
        x_pred = Ad * x_est_real(:,k); 
        P_pred = Ad * P * Ad' + Q;
        
        K = P_pred * Cd' / (Cd * P_pred * Cd' + R_filter);
        x_est_real(:,k+1) = x_pred + K * (y_meas_real(k) - Cd*x_pred);
        P = (eye(4) - K * Cd) * P_pred;
    end

    % 6. PLOT
    figure(8); clf;
    subplot(2,1,1);
    plot(t_sim, y_meas_real, 'Color', [0.4 0.4 0.4]);
    title('Real Suspension Sensor Data (PVS Dataset)');
    ylabel('Accel (m/s^2)'); grid on;
    
    subplot(2,1,2);
    plot(t_sim, x_est_real(1,:), 'b', 'LineWidth', 2);
    title('Kalman Filter Output: Estimated Deflection');
    ylabel('Deflection (m)'); xlabel('Time (s)');
    grid on;
    
    disp('Kaggle Data processed successfully.');

catch ME
    disp('Error reading file. Make sure you downloaded "dataset_mpu_left.csv"');
    disp(['Error Message: ' ME.message]);
end
