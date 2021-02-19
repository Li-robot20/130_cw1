function [Xk_new, Pk_new] = GNSS_KF_solver(X0, P0, time, sate_index, pseudo_range, pseudo_range_rate) 
% Inputs:
%   X0                  State vector
%   P0                  Error covariance
%   time                Current simulation time(s)
%   sate_index          Satellite number
%   Pseudoranges        Corresponding Pseudoranges data
%   Pseudorangesrates   Corresponding Pseudorangesrates data
% Outputs:
%   Xk_new              Updated state vector
%   Pk_new              Updated error covariance


% trasition matrix
tao_s = 0.5;
I3 = eye(3);
I0 = [0 0 0;
    0 0 0;
    0 0 0];
Phi = [I3, tao_s*I3, I0(:,1:2);
       I0, I3, I0(:,1:2);
       I0(1,:), I0(1,:), 1, tao_s;
       I0(1,:), I0(1,:), 0, 1];

% noise covariance matrix
Sa = 5;
Scp = 0.01;
Scf = 0.04;
Q = [1/3*Sa*tao_s^3*I3, 1/2*Sa*tao_s^2*I3, I0(:,1:2);
     1/2*Sa*tao_s^2*I3, Sa*tao_s*I3, I0(:,1:2);
     I0(1,:), I0(1,:), Scp*tao_s + 1/3*Scf*tao_s^3, 1/2*Scf*tao_s^2;
     I0(1,:), I0(1,:), 1/2*Scf*tao_s^2, Scf*tao_s];
 
 % propagate the state
 Xk = Phi * X0;
 
 % propagate the error covariance
 Pk = Phi*P0*Phi.' + Q;
 

index = sate_index;
ra = [];
ra_dot = [];
ua = [];
r_ea = Xk(1:3).';
v_ea = Xk(4:6).';
r_aj = 0;
for i = 1:length(index)
    %ECEF satellite position and velocity
    [r_ej, v_ej] = Satellite_position_and_velocity(0.5*(time-1),index(i));
    
    % Predict the ranges from the position until convergence
    r_aj = predict_raj(r_aj, r_ej, r_ea.');
    r_aj = predict_raj(r_aj, r_ej, r_ea.');
    r_aj = predict_raj(r_aj, r_ej, r_ea.');
    ra = [ra; r_aj];
    
    % Predict the range rates
    u_aj = unit_vector(r_aj, r_ej, r_ea.');
    ua = [ua; u_aj.'];
    raj_dot = predict_raj_dot(u_aj, r_aj, v_ej, r_ej, r_ea.', v_ea.');
    ra_dot = [ra_dot; raj_dot];
end

% following measurement calculation is based on the number of states
num = length(index);

% measurement matrix
Hk = [-ua, zeros(num,3), ones(num,1), zeros(num,1);
      zeros(num,3), -ua, zeros(num,1), ones(num,1)];
  
% measurement noise covariance
sigma_p = 10;
sigma_r = 0.05;
Rk = zeros(2*num, 2*num);
for i = 1:num
    Rk(i,i) = sigma_p^2;
    Rk(i+num, i+num) = sigma_r^2;
end

% Kalman gain matrix
Kk = Pk * Hk.' * inv(Hk * Pk * Hk.' + Rk);

% measurement innovation vector
z = zeros(2*num,1);
for i = 1:num
    z(i) = pseudo_range(i+1) - ra(i) - Xk(7);
    z(i+num) = pseudo_range_rate(i+1) - ra_dot(i) - Xk(8);
end

% update the state
Xk_new = Xk + Kk * z;

% update the error covariance
Pk_new = (eye(num) - Kk * Hk) * Pk;

end