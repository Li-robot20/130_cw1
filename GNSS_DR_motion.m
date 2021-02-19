% heading and forward speed solution
phi = Heading * deg_to_rad;
Vk = (DR(:, 2) + DR(:, 3) + DR(:, 4) + DR(:, 5))/4;

result_1 = [0, Motion_1(1,2:5), Heading(1)];
Motion_4 = [result_1];

% initialise position and velocity
L_k = Motion_1(1,2) * deg_to_rad;
Lambda_k = Motion_1(1,3) * deg_to_rad;
V_Nk = Motion_1(1,4);
V_Ek = Motion_1(1,5);
h_k = Heights(1);

% initialise the state vector
X0 = [0; 0; 0; 0];

% initialise the state estimation error covariance matrix
tao_v = 0.1;
tao_r = 10;
[R_N,R_E] = Radii_of_curvature(L_k);
element = [tao_v^2, tao_v^2, tao_r^2/(R_N + h_k)^2, tao_r^2/(cos(L_k)*cos(L_k)*(R_E + h_k)^2)];
P0 = diag(element);

for i = 2:length(DR)
    %Dead reckoning
    h_k = Heights(i);
    V = 0.5 * [cos(phi(i)) + cos(phi(i-1)); sin(phi(i)) + sin(phi(i-1))] * Vk(i);    
    [R_N,R_E] = Radii_of_curvature(L_k); 
    
    % Dead reckoning solution for position
    Lambda_k = Lambda_k + V(2) * 0.5 / ((R_E + h_k) * cos(L_k));
    L_k = L_k + V(1) * 0.5 / (R_N + h_k);           
    
    % Dead reckoning solution for velocity
    V_Nk = 1.7 * V(1) - 0.7 * V_Nk;
    V_Ek = 1.7 * V(2) - 0.7 * V_Ek;
    
    % Compute the transition matrix
    Phi = eye(4);
    tao_s = 0.5;
    Phi(3,1) = tao_s/(R_N + h_k);
    Phi(4,2) = tao_s/((R_E + h_k) * cos(L_k));
    
    % Compute the system noise covariance matrix
    Q = zeros(4,4);
    S_DR = 0.2;
    Q(1,1) = S_DR * tao_s;
    Q(1,3) = 0.5 * S_DR * tao_s^2 / (R_N + h_k);
    Q(2,2) = S_DR * tao_s;
    Q(2,4) = 0.5 * S_DR * tao_s^2 / (cos(L_k)*(R_E + h_k));
    Q(3,1) = 0.5 * S_DR * tao_s^2 / (R_N + h_k);
    Q(3,3) = 1/3 * S_DR * tao_s^3 / (R_N + h_k)^2;
    Q(4,2) = 0.5 * S_DR * tao_s^2 / (cos(L_k)*(R_E + h_k));
    Q(4,4) = 1/3 * S_DR * tao_s^3 / (cos(L_k)*cos(L_k)*(R_E + h_k)^2);
    
    % Propagate the state estimates and error covariance matrix
    Xk = Phi * X0;
    Pk = Phi * P0 * Phi.'+ Q;

    % Compute the measurement matrix
    Hk = zeros(4,4);
    Hk(1,3) = -1;
    Hk(2,4) = -1;
    Hk(3,1) = -1;
    Hk(4,2) = -1;

    % Compute the measurement noise covariance matrix    
    tao_Gr = 15;
    tao_Gv = 0.02;
    L_GNSS = Motion_1(i,2) * deg_to_rad;
    [R_N,R_E] = Radii_of_curvature(L_k);
    Rk = zeros(4,4);    
    Rk(1,1) = tao_Gr^2/(R_N + h_k)^2;
    Rk(2,2) = tao_Gr^2/(cos(L_GNSS)*cos(L_GNSS)*(R_E + h_k)^2);
    Rk(3,3) = tao_Gv^2;
    Rk(4,4) = tao_Gv^2;

    % Compute the Kalman gain matrix
    Kk = Pk * Hk.' * inv(Hk * Pk * Hk.' + Rk);

    % the measurement innovation vector
    z = zeros(4,1);
    z(1) = Motion_1(i,2) * deg_to_rad - L_k + Xk(3);
    z(2) = Motion_1(i,3) * deg_to_rad - Lambda_k + Xk(4);
    z(3) = Motion_1(i,4) - V_Nk + Xk(1);
    z(4) = Motion_1(i,5) - V_Ek + Xk(2);

    % Update the state estimates
    Xk_new = Xk + Kk * z;

    % Update the error covariance matrix
    Pk_new = (eye(4) - Kk * Hk) * Pk;
    P0 = Pk_new;

    % Correct the DR solution
    L_k = L_k - Xk_new(3);
    Lambda_k = Lambda_k - Xk_new(4);
    V_Nk = V_Nk - Xk_new(1);
    V_Ek = V_Ek - Xk_new(2);

    % store results at each epoch into motion file
    result_single = [0.5*(i-1), L_k * rad_to_deg, Lambda_k * rad_to_deg, V_Nk, V_Ek, Heading(i)];
    Motion_4 = [Motion_4; result_single];
    
end

