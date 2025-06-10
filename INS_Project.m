%% AE 484 Inertial Navigation Systems - Complete Project Implementation
% This script implements integrated navigation with IMU/GPS/Barometer fusion
% Including both EKF/UKF and loosely/tightly coupled approaches
%
% State vector (16 states): x = [p; v; q; b_a; b_g]
%   p: position (3×1) in navigation frame
%   v: velocity (3×1) in navigation frame  
%   q: quaternion (4×1) [qw qx qy qz] from body to nav
%   b_a: accelerometer bias (3×1)    
%   b_g: gyroscope bias (3×1)

clear; close all; clc;
tic;  % Start timer

%% 1. Setup & Configuration
dataFolder = 'AGZ';

% System parameters
g0 = 9.80665;           % standard gravity at sea level (m/s^2)
g_nav = [0; 0; -g0];    % gravity in navigation frame (ENU)
R_earth = 6378137;      % Earth radius (m)

% IMU noise parameters (tune based on sensor specs)
sigma_acc = 0.1;        % accelerometer noise (m/s^2)
sigma_gyro = 0.01;      % gyroscope noise (rad/s)
sigma_ba = 1e-3;        % accel bias random walk (m/s^2/sqrt(s))
sigma_bg = 1e-4;        % gyro bias random walk (rad/s/sqrt(s))

% GPS/Baro measurement noise
R_gps = 5.0;            % GPS position noise (m)
R_baro = 2.0;           % Barometer altitude noise (m)

%% 2. Load & Preprocess Data
fprintf('=== Loading sensor data ===\n');

% Load IMU accelerometer data
accelDataRaw = readmatrix(fullfile(dataFolder,'Log Files','RawAccel.csv'));
% Based on CSV info: Timestamp(1), Error_count(2), x(3), y(4), z(5)
timeAccel_us = accelDataRaw(:,1);
accelData = accelDataRaw(:,3:5);  % N×3 accelerometer data

% Load IMU gyroscope data  
gyroDataRaw = readmatrix(fullfile(dataFolder,'Log Files','RawGyro.csv'));
timeGyro_us = gyroDataRaw(:,1);
gyroData = gyroDataRaw(:,3:5);    % N×3 gyroscope data

% Load GPS data
gpsDataRaw = readmatrix(fullfile(dataFolder,'Log Files','OnboardGPS.csv'));
if size(gpsDataRaw,2) >= 5
    timeGPS_us = gpsDataRaw(:,1);
    % Assuming columns: time, img_id, lon, lat, alt (verify with your data)
    lon_gps = gpsDataRaw(:,3);
    lat_gps = gpsDataRaw(:,4);
    alt_gps = gpsDataRaw(:,5);
    
    % Convert GPS to local ENU coordinates
    lat0 = lat_gps(1);
    lon0 = lon_gps(1);
    alt0 = alt_gps(1);
    
    [gpsENU_E, gpsENU_N, gpsENU_U] = geodetic2enu(lat_gps, lon_gps, alt_gps, ...
                                                   lat0, lon0, alt0, wgs84Ellipsoid);
    gpsENU = [gpsENU_E, gpsENU_N, gpsENU_U];
else
    error('GPS data format not as expected. Check OnboardGPS.csv structure');
end

% Load barometer data
baroDataRaw = readmatrix(fullfile(dataFolder,'Log Files','BarometricPressure.csv'));
timeBaro_us = baroDataRaw(:,1);
baroPressure = baroDataRaw(:,2);

% Convert pressure to altitude
p0 = mean(baroPressure(1:min(50,length(baroPressure))));
baroAlt = 44330 * (1 - (baroPressure/p0).^(1/5.255));
baroAlt = baroAlt - baroAlt(1);  % relative altitude

% Load ground truth
gtRaw = readmatrix(fullfile(dataFolder,'Log Files','GroundTruthAGL.csv'));
gtImgIDs = gtRaw(:,1);
gtPos = gtRaw(:,2:4);  % absolute positions
gtPos = gtPos - gtPos(1,:);  % relative positions

% Convert all timestamps to seconds
timeAccel_s = timeAccel_us / 1e6;
timeGyro_s = timeGyro_us / 1e6;
timeGPS_s = timeGPS_us / 1e6;
timeBaro_s = timeBaro_us / 1e6;

% Create common timeline (use IMU time as master)
timeMaster_s = timeAccel_s;
dt_mean = mean(diff(timeMaster_s));
fs = 1/dt_mean;
N = length(timeMaster_s);

fprintf('Sampling frequency: %.1f Hz\n', fs);
fprintf('Total duration: %.1f seconds\n', timeMaster_s(end)-timeMaster_s(1));

% Interpolate all sensors to common timeline
gpsInterp = interp1(timeGPS_s, gpsENU, timeMaster_s, 'linear', 'extrap');
baroInterp = interp1(timeBaro_s, baroAlt, timeMaster_s, 'linear', 'extrap');

%% 2.5 Initial Attitude Alignment and Bias Estimation
fprintf('\n=== Performing initial alignment ===\n');

% Use first 2 seconds of data (assumed stationary)
n_align = min(round(2 * fs), 200);
accel_static = accelData(1:n_align, :);
gyro_static = gyroData(1:n_align, :);

% Estimate initial biases
gyro_bias0 = mean(gyro_static, 1)';  % 3×1

% For accelerometer bias: in stationary condition, we measure gravity
% The measurement should be [0, 0, g] in body frame if level
% Average measurement
f_static_mean = mean(accel_static, 1)';

% Compute initial attitude from averaged accelerometer
g_body = f_static_mean / norm(f_static_mean);  % normalize

% Calculate roll and pitch from gravity vector
% For z-down body frame: gravity points up in body frame
roll0 = atan2(g_body(2), g_body(3));
pitch0 = atan2(-g_body(1), sqrt(g_body(2)^2 + g_body(3)^2));
yaw0 = 0;  % Cannot determine from accelerometer

% Convert to quaternion
quat0 = euler2quat_custom([roll0, pitch0, yaw0]);

% Now compute expected gravity in body frame given the attitude
C_nb = quat2dcm_custom(quat0)';  % nav to body
g_expected = C_nb * [0; 0; -g0];  % gravity in body frame

% Accelerometer bias is the difference
accel_bias0 = f_static_mean - g_expected;

% Initial state
pos0 = [0; 0; 0];
vel0 = [0; 0; 0];

fprintf('Initial alignment results:\n');
fprintf('  Attitude (RPY): [%.2f, %.2f, %.2f] deg\n', rad2deg([roll0, pitch0, yaw0]));
fprintf('  Accel bias: [%.3f, %.3f, %.3f] m/s^2\n', accel_bias0);
fprintf('  Gyro bias: [%.4f, %.4f, %.4f] rad/s\n', gyro_bias0);

%% 3. Pure INS Mechanization
fprintf('\n=== Running pure INS mechanization ===\n');

pos_INS = zeros(N, 3);
vel_INS = zeros(N, 3);
quat_INS = zeros(N, 4);
att_INS = zeros(N, 3);  % Euler angles for visualization

pos_INS(1,:) = pos0';
vel_INS(1,:) = vel0';
quat_INS(1,:) = quat0;
att_INS(1,:) = [roll0, pitch0, yaw0];

for k = 1:N-1
    dt = timeMaster_s(k+1) - timeMaster_s(k);
    
    % Remove biases from measurements
    f_b = accelData(k,:)' - accel_bias0;
    omega_b = gyroData(k,:)' - gyro_bias0;
    
    % Current quaternion
    q = quat_INS(k,:);
    
    % Attitude update
    Omega = [0, -omega_b(1), -omega_b(2), -omega_b(3);
             omega_b(1), 0, omega_b(3), -omega_b(2);
             omega_b(2), -omega_b(3), 0, omega_b(1);
             omega_b(3), omega_b(2), -omega_b(1), 0];
    
    q_new = (eye(4) + 0.5*Omega*dt) * q';
    q_new = q_new / norm(q_new);
    quat_INS(k+1,:) = q_new';
    
    % Transform to navigation frame
    C_bn = quat2dcm_custom(q_new');
    f_n = C_bn * f_b;
    
    % Velocity and position update
    a_n = f_n + g_nav;
    vel_INS(k+1,:) = vel_INS(k,:) + a_n' * dt;
    pos_INS(k+1,:) = pos_INS(k,:) + vel_INS(k,:) * dt + 0.5 * a_n' * dt^2;
    
    % Store Euler angles
    att_INS(k+1,:) = quat2euler_custom(q_new');
end

%% 4. Extended Kalman Filter (Loosely Coupled)
fprintf('\n=== Running EKF (loosely coupled) ===\n');

n = 16;  % state dimension
x_ekf = zeros(n, N);
P_ekf = zeros(n, n, N);

% Initial state and covariance
x_ekf(:,1) = [pos0; vel0; quat0'; accel_bias0; gyro_bias0];
P_ekf(:,:,1) = diag([ones(1,3)*10, ones(1,3)*1, ones(1,4)*0.1, ...
                     ones(1,3)*0.1, ones(1,3)*0.01]);

% Process noise
Q = zeros(n);
Q(4:6,4:6) = (sigma_acc^2) * eye(3);      % velocity noise
Q(7:10,7:10) = (sigma_gyro^2) * eye(4);   % attitude noise
Q(11:13,11:13) = (sigma_ba^2) * eye(3);   % accel bias noise
Q(14:16,14:16) = (sigma_bg^2) * eye(3);   % gyro bias noise

% Measurement noise
R = diag([R_gps^2, R_gps^2, R_baro^2]);

% Initialize NIS for consistency check
NIS_ekf = zeros(N-1,1);
NIS_ukf = zeros(N-1,1);

for k = 1:N-1
    dt = timeMaster_s(k+1) - timeMaster_s(k);
    
    % Extract state
    p = x_ekf(1:3,k);
    v = x_ekf(4:6,k);
    q = x_ekf(7:10,k);
    ba = x_ekf(11:13,k);
    bg = x_ekf(14:16,k);
    
    % IMU measurements (corrected)
    f_b = accelData(k,:)' - ba;
    omega_b = gyroData(k,:)' - bg;
    
    % === Prediction Step ===
    
    % Attitude propagation
    Omega = [0, -omega_b(1), -omega_b(2), -omega_b(3);
             omega_b(1), 0, omega_b(3), -omega_b(2);
             omega_b(2), -omega_b(3), 0, omega_b(1);
             omega_b(3), omega_b(2), -omega_b(1), 0];
    
    q_pred = (eye(4) + 0.5*Omega*dt) * q;
    q_pred = q_pred / norm(q_pred);
    
    % Navigation frame acceleration
    C_bn = quat2dcm_custom(q);
    f_n = C_bn * f_b;
    a_n = f_n + g_nav;
    
    % Position and velocity propagation
    p_pred = p + v*dt + 0.5*a_n*dt^2;
    v_pred = v + a_n*dt;
    
    % Bias propagation (random walk)
    ba_pred = ba;
    bg_pred = bg;
    
    % Predicted state
    x_pred = [p_pred; v_pred; q_pred; ba_pred; bg_pred];
    
    % State transition matrix
    F = eye(n);
    F(1:3,4:6) = dt * eye(3);
    F(4:6,7:10) = dt * dfn_dq(f_b, q);
    F(4:6,11:13) = -dt * C_bn;
    F(7:10,7:10) = eye(4) + 0.5*dt*Omega;
    F(7:10,14:16) = -0.5*dt*dOmega_domega(q, omega_b);
    
    % Covariance propagation
    P_pred = F * P_ekf(:,:,k) * F' + Q*dt;
    
    % === Update Step ===
    
    % Measurement: GPS East/North + Baro Up
    z = [gpsInterp(k,1); gpsInterp(k,2); baroInterp(k)];
    
    % Measurement model
    H = zeros(3,n);
    H(1:3,1:3) = eye(3);
    
    % Innovation
    y = z - H*x_pred;
    S = H*P_pred*H' + R;
    K = P_pred*H' / S;
    
    % Store NIS for consistency check
    NIS_ekf(k) = y' / S * y;
    
    % State update
    x_ekf(:,k+1) = x_pred + K*y;
    
    % Normalize quaternion
    x_ekf(7:10,k+1) = x_ekf(7:10,k+1) / norm(x_ekf(7:10,k+1));
    
    % Covariance update
    P_ekf(:,:,k+1) = (eye(n) - K*H) * P_pred;
end

%% 5. Unscented Kalman Filter (Loosely Coupled)
fprintf('\n=== Running UKF (loosely coupled) ===\n');

% UKF parameters
alpha = 1e-3;
beta = 2;
kappa = 0;
lambda = alpha^2 * (n + kappa) - n;

% Weights
Wm = zeros(2*n+1,1);
Wc = zeros(2*n+1,1);
Wm(1) = lambda/(n+lambda);
Wc(1) = Wm(1) + (1 - alpha^2 + beta);
for i = 2:2*n+1
    Wm(i) = 1/(2*(n+lambda));
    Wc(i) = Wm(i);
end

x_ukf = zeros(n, N);
P_ukf = zeros(n, n, N);
x_ukf(:,1) = [pos0; vel0; quat0'; accel_bias0; gyro_bias0];
P_ukf(:,:,1) = P_ekf(:,:,1);
NIS_ukf = zeros(N-1,1);

for k = 1:N-1
    dt = timeMaster_s(k+1) - timeMaster_s(k);
    
    % Generate sigma points
    x_k = x_ukf(:,k);
    P_k = P_ukf(:,:,k);
    
    sqrtP = real(sqrtm((n+lambda)*P_k));
    chi = zeros(n, 2*n+1);
    chi(:,1) = x_k;
    for i = 1:n
        chi(:,i+1) = x_k + sqrtP(:,i);
        chi(:,i+n+1) = x_k - sqrtP(:,i);
    end
    
    % Propagate sigma points
    chi_pred = zeros(n, 2*n+1);
    for i = 1:2*n+1
        % Extract states
        p_i = chi(1:3,i);
        v_i = chi(4:6,i);
        q_i = chi(7:10,i);
        ba_i = chi(11:13,i);
        bg_i = chi(14:16,i);
        
        % IMU measurements
        f_b = accelData(k,:)' - ba_i;
        omega_b = gyroData(k,:)' - bg_i;
        
        % Propagate
        Omega = [0, -omega_b(1), -omega_b(2), -omega_b(3);
                 omega_b(1), 0, omega_b(3), -omega_b(2);
                 omega_b(2), -omega_b(3), 0, omega_b(1);
                 omega_b(3), omega_b(2), -omega_b(1), 0];
        
        q_i_pred = (eye(4) + 0.5*Omega*dt) * q_i;
        q_i_pred = q_i_pred / norm(q_i_pred);
        
        C_bn_i = quat2dcm_custom(q_i);
        f_n_i = C_bn_i * f_b;
        a_n_i = f_n_i + g_nav;
        
        p_i_pred = p_i + v_i*dt + 0.5*a_n_i*dt^2;
        v_i_pred = v_i + a_n_i*dt;
        
        chi_pred(:,i) = [p_i_pred; v_i_pred; q_i_pred; ba_i; bg_i];
    end
    
    % Predicted mean and covariance
    x_pred = zeros(n,1);
    for i = 1:2*n+1
        x_pred = x_pred + Wm(i)*chi_pred(:,i);
    end
    x_pred(7:10) = x_pred(7:10) / norm(x_pred(7:10));
    
    P_pred = Q*dt;
    for i = 1:2*n+1
        dx = chi_pred(:,i) - x_pred;
        P_pred = P_pred + Wc(i)*(dx*dx');
    end
    
    % Measurement update
    z = [gpsInterp(k,1); gpsInterp(k,2); baroInterp(k)];
    
    % Predict measurements
    Z = zeros(3, 2*n+1);
    for i = 1:2*n+1
        Z(:,i) = chi_pred(1:3,i);
    end
    
    z_pred = zeros(3,1);
    for i = 1:2*n+1
        z_pred = z_pred + Wm(i)*Z(:,i);
    end
    
    Pzz = R;
    Pxz = zeros(n,3);
    for i = 1:2*n+1
        dz = Z(:,i) - z_pred;
        dx = chi_pred(:,i) - x_pred;
        Pzz = Pzz + Wc(i)*(dz*dz');
        Pxz = Pxz + Wc(i)*(dx*dz');
    end
    
    K = Pxz / Pzz;
    innovation = z - z_pred;
    NIS_ukf(k) = innovation' / Pzz * innovation;
    x_ukf(:,k+1) = x_pred + K*innovation;
    x_ukf(7:10,k+1) = x_ukf(7:10,k+1) / norm(x_ukf(7:10,k+1));
    P_ukf(:,:,k+1) = P_pred - K*Pzz*K';
end

%% 6. Tightly Coupled EKF (Bonus Implementation)
fprintf('\n=== Running EKF (tightly coupled) ===\n');

% For tightly coupled, we need pseudorange measurements
% Since dataset only has GPS positions, we'll simulate pseudoranges
% In real implementation, you'd use RINEX data

% Simulate satellite positions (simplified)
n_sats = 8;  % number of visible satellites
sat_pos = zeros(3, n_sats, N);
for i = 1:n_sats
    % Simplified satellite motion (circular orbits)
    angle = 2*pi*i/n_sats + timeMaster_s*2*pi/(12*3600);  % 12hr orbit
    radius = 26600e3;  % GPS orbit radius
    sat_pos(1,i,:) = radius * cos(angle);
    sat_pos(2,i,:) = radius * sin(angle) * cos(pi/6);
    sat_pos(3,i,:) = radius * sin(angle) * sin(pi/6);
end

% Augment state with clock bias and drift
n_tc = 18;  % [pos(3), vel(3), quat(4), ba(3), bg(3), clk_bias(1), clk_drift(1)]
x_ekf_tc = zeros(n_tc, N);
P_ekf_tc = zeros(n_tc, n_tc, N);

% Initial state
x_ekf_tc(1:16,1) = x_ekf(:,1);
x_ekf_tc(17:18,1) = [0; 0];  % clock bias and drift
P_ekf_tc(:,:,1) = blkdiag(P_ekf(:,:,1), diag([100^2, 1^2]));

% Process noise for clock
Q_tc = blkdiag(Q, diag([0.1^2, 0.01^2]));

% Measurement noise for pseudoranges
R_pr = 5^2 * eye(n_sats);  % 5m pseudorange noise

for k = 1:N-1
    dt = timeMaster_s(k+1) - timeMaster_s(k);
    
    % Prediction (similar to loosely coupled for INS states)
    x_pred_tc = zeros(n_tc,1);
    x_pred_tc(1:16) = ekf_predict_state(x_ekf_tc(1:16,k), accelData(k,:)', gyroData(k,:)', dt);
    x_pred_tc(17) = x_ekf_tc(17,k) + x_ekf_tc(18,k)*dt;  % clock bias
    x_pred_tc(18) = x_ekf_tc(18,k);  % clock drift
    
    % State transition matrix (recompute for tightly coupled)
    F_ins = eye(16);
    F_ins(1:3,4:6) = dt * eye(3);
    % Note: For simplicity, using linearized F matrix
    
    F_tc = eye(n_tc);
    F_tc(1:16,1:16) = F_ins;  % INS portion
    F_tc(17,18) = dt;         % clock bias/drift
    
    P_pred_tc = F_tc * P_ekf_tc(:,:,k) * F_tc' + Q_tc*dt;
    
    % Measurement update with pseudoranges
    H_pr = zeros(n_sats, n_tc);
    z_pr = zeros(n_sats, 1);
    z_pred_pr = zeros(n_sats, 1);
    
    p_pred = x_pred_tc(1:3);
    clk_bias = x_pred_tc(17);
    
    for i = 1:n_sats
        % True range (simulated)
        true_range = norm(sat_pos(:,i,k) - gpsInterp(k,:)');
        z_pr(i) = true_range + randn*sqrt(R_pr(1,1));  % add noise
        
        % Predicted range
        pred_range = norm(sat_pos(:,i,k) - p_pred);
        z_pred_pr(i) = pred_range + clk_bias;
        
        % Measurement Jacobian
        los = (p_pred - sat_pos(:,i,k)) / pred_range;  % line of sight
        H_pr(i,1:3) = los';
        H_pr(i,17) = 1;  % clock bias
    end
    
    % Kalman update
    y_pr = z_pr - z_pred_pr;
    S_pr = H_pr*P_pred_tc*H_pr' + R_pr;
    K_pr = P_pred_tc*H_pr' / S_pr;
    
    x_ekf_tc(:,k+1) = x_pred_tc + K_pr*y_pr;
    x_ekf_tc(7:10,k+1) = x_ekf_tc(7:10,k+1) / norm(x_ekf_tc(7:10,k+1));
    P_ekf_tc(:,:,k+1) = (eye(n_tc) - K_pr*H_pr) * P_pred_tc;
end

%% 7. Results Processing and Visualization

% Find valid ground truth indices
% Map ground truth to time using image IDs if available
if size(gpsDataRaw,2) >= 2 && any(gpsDataRaw(:,2) > 0)
    % Map GT image IDs to GPS timestamps
    imgIDs_GPS = gpsDataRaw(:,2);
    gtTime_s = zeros(size(gtImgIDs));
    for i = 1:length(gtImgIDs)
        idx = find(imgIDs_GPS == gtImgIDs(i), 1);
        if ~isempty(idx)
            gtTime_s(i) = timeGPS_s(idx);
        else
            gtTime_s(i) = NaN;
        end
    end
    % Remove invalid mappings
    validMap = ~isnan(gtTime_s);
    gtTime_s = gtTime_s(validMap);
    gtPos = gtPos(validMap,:);
    
    % Interpolate to master timeline
    idxValidGT = (timeMaster_s >= min(gtTime_s)) & (timeMaster_s <= max(gtTime_s));
    gtPos_interp = interp1(gtTime_s, gtPos, timeMaster_s(idxValidGT), 'linear');
else
    % If no image ID mapping, assume GT is sampled at same rate as GPS
    if length(gtPos) == length(timeGPS_s)
        gtTime_s = timeGPS_s;
        idxValidGT = (timeMaster_s >= min(gtTime_s)) & (timeMaster_s <= max(gtTime_s));
        gtPos_interp = interp1(gtTime_s, gtPos, timeMaster_s(idxValidGT), 'linear');
    else
        % Last resort: assume GT aligns with master timeline
        fprintf('Warning: Cannot align ground truth with timestamps. Using direct alignment.\n');
        nGT = min(size(gtPos,1), N);
        idxValidGT = false(N,1);
        idxValidGT(1:nGT) = true;
        gtPos_interp = gtPos(1:nGT,:);
    end
end

% Extract estimates at valid times
pos_INS_valid = pos_INS(idxValidGT,:);
pos_ekf_valid = x_ekf(1:3,idxValidGT)';
vel_ekf_valid = x_ekf(4:6,idxValidGT)';
pos_ukf_valid = x_ukf(1:3,idxValidGT)';
vel_ukf_valid = x_ukf(4:6,idxValidGT)';
pos_ekf_tc_valid = x_ekf_tc(1:3,idxValidGT)';

% Convert quaternions to Euler angles for visualization
euler_INS = zeros(sum(idxValidGT),3);
euler_ekf = zeros(sum(idxValidGT),3);
euler_ukf = zeros(sum(idxValidGT),3);

idx_valid = find(idxValidGT);
for i = 1:length(idx_valid)
    euler_INS(i,:) = quat2euler_custom(quat_INS(idx_valid(i),:));
    euler_ekf(i,:) = quat2euler_custom(x_ekf(7:10,idx_valid(i))');
    euler_ukf(i,:) = quat2euler_custom(x_ukf(7:10,idx_valid(i))');
end

%% 8. Generate All Required Plots

% === 3D Trajectory Comparison ===
figure('Name','3D Trajectory Comparison','Position',[100,100,800,600]);
plot3(gtPos_interp(:,1), gtPos_interp(:,2), gtPos_interp(:,3), 'g-', 'LineWidth', 2);
hold on;
plot3(pos_INS_valid(:,1), pos_INS_valid(:,2), pos_INS_valid(:,3), 'r--', 'LineWidth', 1.5);
plot3(pos_ekf_valid(:,1), pos_ekf_valid(:,2), pos_ekf_valid(:,3), 'b-', 'LineWidth', 1.5);
plot3(pos_ukf_valid(:,1), pos_ukf_valid(:,2), pos_ukf_valid(:,3), 'm-.', 'LineWidth', 1.5);
plot3(pos_ekf_tc_valid(:,1), pos_ekf_tc_valid(:,2), pos_ekf_tc_valid(:,3), 'c:', 'LineWidth', 1.5);
legend('Ground Truth','Pure INS','EKF (LC)','UKF (LC)','EKF (TC)','Location','best');
xlabel('East (m)'); ylabel('North (m)'); zlabel('Up (m)');
title('3D Trajectory Comparison');
grid on; axis equal;
view(45,30);

% === 2D Horizontal Trajectory ===
figure('Name','Horizontal Trajectory Comparison','Position',[100,100,800,600]);
plot(gtPos_interp(:,1), gtPos_interp(:,2), 'g-', 'LineWidth', 2);
hold on;
plot(pos_INS_valid(:,1), pos_INS_valid(:,2), 'r--', 'LineWidth', 1.5);
plot(pos_ekf_valid(:,1), pos_ekf_valid(:,2), 'b-', 'LineWidth', 1.5);
plot(pos_ukf_valid(:,1), pos_ukf_valid(:,2), 'm-.', 'LineWidth', 1.5);
plot(pos_ekf_tc_valid(:,1), pos_ekf_tc_valid(:,2), 'c:', 'LineWidth', 1.5);
legend('Ground Truth','Pure INS','EKF (LC)','UKF (LC)','EKF (TC)','Location','best');
xlabel('East (m)'); ylabel('North (m)');
title('Horizontal (XY) Trajectory Comparison');
grid on; axis equal;

% === Position Errors with 3-sigma bounds (All 3 axes) ===
t_valid = timeMaster_s(idxValidGT) - timeMaster_s(find(idxValidGT,1));

% Position errors
err_pos_ekf = pos_ekf_valid - gtPos_interp;
err_pos_ukf = pos_ukf_valid - gtPos_interp;
err_pos_ins = pos_INS_valid - gtPos_interp;

% 3-sigma bounds
sigma_pos_ekf = zeros(3, 3, sum(idxValidGT));
sigma_pos_ukf = zeros(3, 3, sum(idxValidGT));
idx_valid = find(idxValidGT);
for i = 1:length(idx_valid)
    sigma_pos_ekf(:,:,i) = 3*sqrt(P_ekf(1:3,1:3,idx_valid(i)));
    sigma_pos_ukf(:,:,i) = 3*sqrt(P_ukf(1:3,1:3,idx_valid(i)));
end

figure('Name','Position Errors with 3σ Bounds','Position',[100,100,1200,800]);

% East error
subplot(3,2,1);
plot(t_valid, err_pos_ekf(:,1), 'b-', 'LineWidth', 1);
hold on;
plot(t_valid, squeeze(sigma_pos_ekf(1,1,:)), 'b--');
plot(t_valid, -squeeze(sigma_pos_ekf(1,1,:)), 'b--');
xlabel('Time (s)'); ylabel('East Error (m)');
title('EKF East Position Error');
legend('Error','±3σ','Location','best');
grid on;

subplot(3,2,2);
plot(t_valid, err_pos_ukf(:,1), 'm-', 'LineWidth', 1);
hold on;
plot(t_valid, squeeze(sigma_pos_ukf(1,1,:)), 'm--');
plot(t_valid, -squeeze(sigma_pos_ukf(1,1,:)), 'm--');
xlabel('Time (s)'); ylabel('East Error (m)');
title('UKF East Position Error');
legend('Error','±3σ','Location','best');
grid on;

% North error
subplot(3,2,3);
plot(t_valid, err_pos_ekf(:,2), 'b-', 'LineWidth', 1);
hold on;
plot(t_valid, squeeze(sigma_pos_ekf(2,2,:)), 'b--');
plot(t_valid, -squeeze(sigma_pos_ekf(2,2,:)), 'b--');
xlabel('Time (s)'); ylabel('North Error (m)');
title('EKF North Position Error');
grid on;

subplot(3,2,4);
plot(t_valid, err_pos_ukf(:,2), 'm-', 'LineWidth', 1);
hold on;
plot(t_valid, squeeze(sigma_pos_ukf(2,2,:)), 'm--');
plot(t_valid, -squeeze(sigma_pos_ukf(2,2,:)), 'm--');
xlabel('Time (s)'); ylabel('North Error (m)');
title('UKF North Position Error');
grid on;

% Up error
subplot(3,2,5);
plot(t_valid, err_pos_ekf(:,3), 'b-', 'LineWidth', 1);
hold on;
plot(t_valid, squeeze(sigma_pos_ekf(3,3,:)), 'b--');
plot(t_valid, -squeeze(sigma_pos_ekf(3,3,:)), 'b--');
xlabel('Time (s)'); ylabel('Up Error (m)');
title('EKF Up Position Error');
grid on;

subplot(3,2,6);
plot(t_valid, err_pos_ukf(:,3), 'm-', 'LineWidth', 1);
hold on;
plot(t_valid, squeeze(sigma_pos_ukf(3,3,:)), 'm--');
plot(t_valid, -squeeze(sigma_pos_ukf(3,3,:)), 'm--');
xlabel('Time (s)'); ylabel('Up Error (m)');
title('UKF Up Position Error');
grid on;

% === Velocity Errors with 3-sigma bounds ===
% Need ground truth velocity (differentiate position)
vel_gt = zeros(size(gtPos_interp));
for i = 2:size(gtPos_interp,1)-1
    vel_gt(i,:) = (gtPos_interp(i+1,:) - gtPos_interp(i-1,:)) / (2*dt_mean);
end

err_vel_ekf = vel_ekf_valid - vel_gt;
err_vel_ukf = vel_ukf_valid - vel_gt;

sigma_vel_ekf = zeros(3, 3, sum(idxValidGT));
sigma_vel_ukf = zeros(3, 3, sum(idxValidGT));
for i = 1:length(idx_valid)
    sigma_vel_ekf(:,:,i) = 3*sqrt(P_ekf(4:6,4:6,idx_valid(i)));
    sigma_vel_ukf(:,:,i) = 3*sqrt(P_ukf(4:6,4:6,idx_valid(i)));
end

figure('Name','Velocity Errors with 3σ Bounds','Position',[100,100,1200,800]);

% Similar subplot structure for velocity...
vel_labels = {'East', 'North', 'Up'};
for i = 1:3
    subplot(3,2,2*i-1);
    plot(t_valid, err_vel_ekf(:,i), 'b-', 'LineWidth', 1);
    hold on;
    plot(t_valid, squeeze(sigma_vel_ekf(i,i,:)), 'b--');
    plot(t_valid, -squeeze(sigma_vel_ekf(i,i,:)), 'b--');
    xlabel('Time (s)'); ylabel(sprintf('%s Vel Error (m/s)', vel_labels{i}));
    title(sprintf('EKF %s Velocity Error', vel_labels{i}));
    legend('Error','±3σ','Location','best');
    grid on;
    
    subplot(3,2,2*i);
    plot(t_valid, err_vel_ukf(:,i), 'm-', 'LineWidth', 1);
    hold on;
    plot(t_valid, squeeze(sigma_vel_ukf(i,i,:)), 'm--');
    plot(t_valid, -squeeze(sigma_vel_ukf(i,i,:)), 'm--');
    xlabel('Time (s)'); ylabel(sprintf('%s Vel Error (m/s)', vel_labels{i}));
    title(sprintf('UKF %s Velocity Error', vel_labels{i}));
    legend('Error','±3σ','Location','best');
    grid on;
end

% === Attitude Errors ===
% If ground truth attitude is not available, compute from trajectory
% Compute ground truth attitude from velocity
vel_gt_smooth = zeros(size(vel_gt));
for i = 3:size(vel_gt,1)-2
    vel_gt_smooth(i,:) = mean(vel_gt(i-2:i+2,:), 1);
end

% Compute heading from velocity
yaw_gt = atan2(vel_gt_smooth(:,2), vel_gt_smooth(:,1));

% For roll and pitch, we need IMU data or assume zero
roll_gt = zeros(size(yaw_gt));
pitch_gt = zeros(size(yaw_gt));
euler_gt = [roll_gt, pitch_gt, yaw_gt];

% Compute attitude errors
err_att_ekf = euler_ekf - euler_gt;
err_att_ukf = euler_ukf - euler_gt;
err_att_ins = att_INS(idx_valid,:) - euler_gt;

% Wrap angles to [-pi, pi]
err_att_ekf = wrapToPi(err_att_ekf);
err_att_ukf = wrapToPi(err_att_ukf);
err_att_ins = wrapToPi(err_att_ins);

figure('Name','Attitude Errors vs Ground Truth','Position',[100,100,1200,800]);

att_labels = {'Roll', 'Pitch', 'Yaw'};
for i = 1:3
    subplot(3,1,i);
    plot(t_valid, rad2deg(err_att_ins(:,i)), 'r--', 'LineWidth', 1);
    hold on;
    plot(t_valid, rad2deg(err_att_ekf(:,i)), 'b-', 'LineWidth', 1);
    plot(t_valid, rad2deg(err_att_ukf(:,i)), 'm-.', 'LineWidth', 1);
    xlabel('Time (s)'); ylabel(sprintf('%s Error (deg)', att_labels{i}));
    title(sprintf('%s Error Comparison', att_labels{i}));
    legend('Pure INS', 'EKF', 'UKF', 'Location','best');
    grid on;
    if i < 3  % Roll and pitch ground truth is approximate
        ylim([-180, 180]);
    end
end

% === Pure INS Error Analysis ===
figure('Name','Pure INS Position Errors','Position',[100,100,800,600]);

subplot(3,1,1);
plot(t_valid, err_pos_ins(:,1), 'r-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('East Error (m)');
title('Pure INS East Position Error (No Corrections)');
grid on;

subplot(3,1,2);
plot(t_valid, err_pos_ins(:,2), 'r-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('North Error (m)');
title('Pure INS North Position Error (No Corrections)');
grid on;

subplot(3,1,3);
plot(t_valid, err_pos_ins(:,3), 'r-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Up Error (m)');
title('Pure INS Up Position Error (No Corrections)');
grid on;

% === Loosely vs Tightly Coupled Comparison ===
figure('Name','Loosely vs Tightly Coupled Integration','Position',[100,100,1200,600]);

% Position error comparison
subplot(2,2,1);
err_ekf_lc = sqrt(sum(err_pos_ekf.^2, 2));
err_ekf_tc = sqrt(sum((pos_ekf_tc_valid - gtPos_interp).^2, 2));
plot(t_valid, err_ekf_lc, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_valid, err_ekf_tc, 'c--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Position Error (m)');
title('Position Error: Loosely vs Tightly Coupled');
legend('EKF (Loosely Coupled)', 'EKF (Tightly Coupled)', 'Location','best');
grid on;

% Horizontal trajectory
subplot(2,2,2);
plot(pos_ekf_valid(:,1), pos_ekf_valid(:,2), 'b-', 'LineWidth', 1.5);
hold on;
plot(pos_ekf_tc_valid(:,1), pos_ekf_tc_valid(:,2), 'c--', 'LineWidth', 1.5);
plot(gtPos_interp(:,1), gtPos_interp(:,2), 'g-', 'LineWidth', 2);
xlabel('East (m)'); ylabel('North (m)');
title('Horizontal Trajectory Comparison');
legend('EKF (LC)', 'EKF (TC)', 'Ground Truth', 'Location','best');
grid on; axis equal;

% Clock bias estimation (tightly coupled only)
subplot(2,2,3);
clock_bias = x_ekf_tc(17,idxValidGT);
plot(t_valid, clock_bias, 'c-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Clock Bias (m)');
title('GPS Clock Bias Estimation (Tightly Coupled)');
grid on;

% RMSE comparison over time
subplot(2,2,4);
window = 100;  % samples
rmse_lc = zeros(length(t_valid)-window+1, 1);
rmse_tc = zeros(length(t_valid)-window+1, 1);
for i = 1:length(rmse_lc)
    rmse_lc(i) = sqrt(mean(err_ekf_lc(i:i+window-1).^2));
    rmse_tc(i) = sqrt(mean(err_ekf_tc(i:i+window-1).^2));
end
plot(t_valid(1:end-window+1), rmse_lc, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_valid(1:end-window+1), rmse_tc, 'c--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('RMSE (m)');
title('Windowed RMSE Comparison');
legend('Loosely Coupled', 'Tightly Coupled', 'Location','best');
grid on;

% === RMSE Summary ===
RMSE_ins = sqrt(mean(err_pos_ins.^2, 1));
RMSE_ekf = sqrt(mean(err_pos_ekf.^2, 1));
RMSE_ukf = sqrt(mean(err_pos_ukf.^2, 1));
RMSE_ekf_tc = sqrt(mean((pos_ekf_tc_valid - gtPos_interp).^2, 1));

fprintf('\n=== RMSE Results ===\n');
fprintf('Pure INS:  East=%.2fm, North=%.2fm, Up=%.2fm, Total=%.2fm\n', ...
    RMSE_ins(1), RMSE_ins(2), RMSE_ins(3), norm(RMSE_ins));
fprintf('EKF (LC):  East=%.2fm, North=%.2fm, Up=%.2fm, Total=%.2fm\n', ...
    RMSE_ekf(1), RMSE_ekf(2), RMSE_ekf(3), norm(RMSE_ekf));
fprintf('UKF (LC):  East=%.2fm, North=%.2fm, Up=%.2fm, Total=%.2fm\n', ...
    RMSE_ukf(1), RMSE_ukf(2), RMSE_ukf(3), norm(RMSE_ukf));
fprintf('EKF (TC):  East=%.2fm, North=%.2fm, Up=%.2fm, Total=%.2fm\n', ...
    RMSE_ekf_tc(1), RMSE_ekf_tc(2), RMSE_ekf_tc(3), norm(RMSE_ekf_tc));

% === Bias Estimates Plot ===
figure('Name','IMU Bias Estimates','Position',[100,100,1200,600]);

subplot(2,1,1);
plot(t_valid, x_ekf(11:13,idxValidGT)', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Accel Bias (m/s²)');
title('EKF Accelerometer Bias Estimates');
legend('b_{ax}','b_{ay}','b_{az}','Location','best');
grid on;

subplot(2,1,2);
plot(t_valid, rad2deg(x_ekf(14:16,idxValidGT)'), 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Gyro Bias (deg/s)');
title('EKF Gyroscope Bias Estimates');
legend('b_{gx}','b_{gy}','b_{gz}','Location','best');
grid on;

% === Filter Consistency Check (NIS) ===
figure('Name','Filter Consistency Check','Position',[100,100,800,600]);

% Chi-squared 95% confidence bounds for 3 DOF
chi2_lower = chi2inv(0.025, 3);
chi2_upper = chi2inv(0.975, 3);

subplot(2,1,1);
plot(timeMaster_s(2:end), NIS_ekf, 'b-', 'LineWidth', 0.5);
hold on;
plot(timeMaster_s([2,end]), [chi2_lower, chi2_lower], 'r--', 'LineWidth', 1.5);
plot(timeMaster_s([2,end]), [chi2_upper, chi2_upper], 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('NIS');
title('EKF Normalized Innovation Squared (NIS)');
legend('NIS', '95% Bounds', 'Location','best');
grid on;
ylim([0, max(20, max(NIS_ekf))]);

subplot(2,1,2);
plot(timeMaster_s(2:end), NIS_ukf, 'm-', 'LineWidth', 0.5);
hold on;
plot(timeMaster_s([2,end]), [chi2_lower, chi2_lower], 'r--', 'LineWidth', 1.5);
plot(timeMaster_s([2,end]), [chi2_upper, chi2_upper], 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('NIS');
title('UKF Normalized Innovation Squared (NIS)');
legend('NIS', '95% Bounds', 'Location','best');
grid on;
ylim([0, max(20, max(NIS_ukf))]);

% Print consistency statistics
fprintf('\n=== Filter Consistency (NIS) ===\n');
fprintf('EKF: %.1f%% of NIS values within 95%% bounds\n', ...
    100*sum(NIS_ekf > chi2_lower & NIS_ekf < chi2_upper)/length(NIS_ekf));
fprintf('UKF: %.1f%% of NIS values within 95%% bounds\n', ...
    100*sum(NIS_ukf > chi2_lower & NIS_ukf < chi2_upper)/length(NIS_ukf));

%% Helper Functions

function q = euler2quat_custom(euler)
    % Convert Euler angles to quaternion
    r = euler(1); p = euler(2); y = euler(3);
    
    cr = cos(r/2); sr = sin(r/2);
    cp = cos(p/2); sp = sin(p/2);
    cy = cos(y/2); sy = sin(y/2);
    
    qw = cr*cp*cy + sr*sp*sy;
    qx = sr*cp*cy - cr*sp*sy;
    qy = cr*sp*cy + sr*cp*sy;
    qz = cr*cp*sy - sr*sp*cy;
    
    q = [qw, qx, qy, qz];
    q = q / norm(q);
end

function euler = quat2euler_custom(q)
    % Convert quaternion to Euler angles
    qw = q(1); qx = q(2); qy = q(3); qz = q(4);
    
    % Roll
    sinr_cosp = 2 * (qw*qx + qy*qz);
    cosr_cosp = 1 - 2 * (qx^2 + qy^2);
    roll = atan2(sinr_cosp, cosr_cosp);
    
    % Pitch
    sinp = 2 * (qw*qy - qz*qx);
    if abs(sinp) >= 1
        pitch = sign(sinp) * pi/2;
    else
        pitch = asin(sinp);
    end
    
    % Yaw
    siny_cosp = 2 * (qw*qz + qx*qy);
    cosy_cosp = 1 - 2 * (qy^2 + qz^2);
    yaw = atan2(siny_cosp, cosy_cosp);
    
    euler = [roll, pitch, yaw];
end

function C = quat2dcm_custom(q)
    % Convert quaternion to DCM
    qw = q(1); qx = q(2); qy = q(3); qz = q(4);
    
    C = [ 1-2*(qy^2+qz^2),   2*(qx*qy-qz*qw),   2*(qx*qz+qy*qw);
          2*(qx*qy+qz*qw),   1-2*(qx^2+qz^2),   2*(qy*qz-qx*qw);
          2*(qx*qz-qy*qw),   2*(qy*qz+qx*qw),   1-2*(qx^2+qy^2) ];
end

function S = skew(v)
    % Skew-symmetric matrix
    S = [  0   -v(3)  v(2);
          v(3)   0   -v(1);
         -v(2)  v(1)   0  ];
end

function J = dfn_dq(f_b, q)
    % Jacobian of f_n with respect to quaternion
    qw = q(1); qx = q(2); qy = q(3); qz = q(4);
    fx = f_b(1); fy = f_b(2); fz = f_b(3);
    
    J = 2 * [ qw*fx - qz*fy + qy*fz,  qx*fx + qy*fy + qz*fz, -qy*fx + qx*fy + qw*fz, -qz*fx - qw*fy + qx*fz;
              qz*fx + qw*fy - qx*fz,  qy*fx - qx*fy - qw*fz,  qx*fx + qy*fy + qz*fz,  qw*fx - qz*fy + qy*fz;
             -qy*fx + qx*fy + qw*fz,  qz*fx + qw*fy - qx*fz, -qw*fx + qz*fy - qy*fz,  qx*fx + qy*fy + qz*fz ];
end

function W = dOmega_domega(q, omega_b)
    % Jacobian of Omega matrix with respect to omega
    % Returns the derivative of Omega*q with respect to omega_b
    qw = q(1); qx = q(2); qy = q(3); qz = q(4);
    
    W = [ qx,  qy,  qz;
         -qw,  qz, -qy;
         -qz, -qw,  qx;
          qy, -qx, -qw];
end

function x_pred = ekf_predict_state(x, f_b, omega_b, dt)
    % State prediction for EKF/UKF
    p = x(1:3);
    v = x(4:6);
    q = x(7:10);
    ba = x(11:13);
    bg = x(14:16);
    
    % Corrected measurements
    f_b_corr = f_b - ba;
    omega_b_corr = omega_b - bg;
    
    % Attitude update
    Omega = [0, -omega_b_corr(1), -omega_b_corr(2), -omega_b_corr(3);
             omega_b_corr(1), 0, omega_b_corr(3), -omega_b_corr(2);
             omega_b_corr(2), -omega_b_corr(3), 0, omega_b_corr(1);
             omega_b_corr(3), omega_b_corr(2), -omega_b_corr(1), 0];
    
    q_pred = (eye(4) + 0.5*Omega*dt) * q;
    q_pred = q_pred / norm(q_pred);
    
    % Navigation update
    C_bn = quat2dcm_custom(q);
    f_n = C_bn * f_b_corr;
    g0 = 9.81;
    a_n = f_n + [0; 0; -g0];
    
    p_pred = p + v*dt + 0.5*a_n*dt^2;
    v_pred = v + a_n*dt;
    
    x_pred = [p_pred; v_pred; q_pred; ba; bg];
end

function [E, N, U] = geodetic2enu(lat, lon, h, lat0, lon0, h0, ellipsoid)
    % Simple geodetic to ENU conversion
    % For full implementation, use MATLAB Mapping Toolbox
    R = 6378137;  % Earth radius
    E = R * deg2rad(lon - lon0) .* cos(deg2rad(lat0));
    N = R * deg2rad(lat - lat0);
    U = h - h0;
end

function angles = wrapToPi(angles)
    % Wrap angles to [-pi, pi]
    angles = mod(angles + pi, 2*pi) - pi;
end

% WGS84 ellipsoid (simplified)
wgs84Ellipsoid = struct('a', 6378137, 'f', 1/298.257223563);