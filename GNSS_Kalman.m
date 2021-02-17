function [Xk_new, Pk_new] = GNSS_Kalman(X0, P0, time, pseudo_range, pseudo_range_rate) 
%b trasition matrix
tao_s = 0.5;
I3 = eye(3);
I0 = [0 0 0;
    0 0 0;
    0 0 0];
Phi = [I3, tao_s*I3, I0(:,1:2);
       I0, I3, I0(:,1:2);
       I0(1,:), I0(1,:), 1, tao_s;
       I0(1,:), I0(1,:), 0, 1];

%c noise covariance matrix
Sa = 5;
Scp = 0.01;
Scf = 0.04;
Q = [1/3*Sa*tao_s^3*I3, 1/2*Sa*tao_s^2*I3, I0(:,1:2);
     1/2*Sa*tao_s^2*I3, Sa*tao_s*I3, I0(:,1:2);
     I0(1,:), I0(1,:), Scp*tao_s + 1/3*Scf*tao_s^3, 1/2*Scf*tao_s^2;
     I0(1,:), I0(1,:), 1/2*Scf*tao_s^2, Scf*tao_s];
 
 %d propagate the state
 Xk = Phi * X0;
 
 %e propagate the error covariance
 Pk = Phi*P0*Phi.' + Q;
 
 %f, g 
index = [5, 6, 7, 9, 10, 11, 15, 30];
re = [];
ve = [];
ra = [];
ra_dot = [];
ua = [];
r_ea = Xk(1:3).';
v_ea = Xk(4:6).';
r_aj = 0;
for i = 1:length(index)
    [r_ej, v_ej] = Satellite_position_and_velocity(0.5*(time-1),index(i));
    re = [re; r_ej];
    ve = [ve; v_ej];
    r_aj = calculate_raj(r_aj, r_ej, r_ea.');
    r_aj = calculate_raj(r_aj, r_ej, r_ea.');
    r_aj = calculate_raj(r_aj, r_ej, r_ea.');

    ra = [ra; r_aj];
    u_aj = calculate_uaj(r_aj, r_ej, r_ea.');
    ua = [ua; u_aj.'];
    raj_dot = calculate_raj_dot(r_aj, r_ej.', r_ea.', v_ej.', v_ea.', u_aj);
    ra_dot = [ra_dot; raj_dot];
end

%h measurement matrix
Hk = [-ua, zeros(8,3), ones(8,1), zeros(8,1);
      zeros(8,3), -ua, zeros(8,1), ones(8,1)];
  
%j measurement noise covariance
sigma_p = 10;
sigma_r = 0.05;
Rk = zeros(16, 16);
for i = 1:8
    Rk(i,i) = sigma_p^2;
    Rk(i+8, i+8) = sigma_r^2;
end

%k Kalman gain matrix
Kk = Pk * Hk.' * inv(Hk * Pk * Hk.' + Rk);

z = zeros(16,1);
for i = 1:8
    z(i) = pseudo_range(i+1) - ra(i) - Xk(7);
    z(i+8) = pseudo_range_rate(i+1) - ra_dot(i) - Xk(8);
end

%m update the state
Xk_new = Xk + Kk * z;

%n update error covariance
Pk_new = (eye(8) - Kk * Hk) * Pk;

end