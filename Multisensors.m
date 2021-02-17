phi = Motion(:, 6) * deg_to_rad;
Vk = (DR(:, 2) + DR(:, 3) + DR(:, 4) + DR(:, 5))/4;
V_N0 = Motion(1,4);
V_E0 = Motion(1,5);

result_1 = [0, Motion(1,2:5)];
Motion_3 = [result_1];

X0 = [0; 0; 0; 0];
tao_v = 0.1;
tao_r = 10;

L0 = Motion(1,2) * deg_to_rad;
Lambda0 = Motion(1,3) * deg_to_rad;

h0 = Motion(1,7);
V = 0.5 * [cos(phi(2)) + cos(phi(1)); sin(phi(2)) + sin(phi(1))] * Vk(1);    
[R_N,R_E] = Radii_of_curvature(L0);
element = [tao_v^2, tao_v^2, tao_r^2/(R_N + h0)^2, tao_r^2/(cos(L0)*cos(L0)*(R_E + h0)^2)];
P0 = diag(element);

for i = 2:length(DR)
    %DR
    h0 = Motion(i,7);
    V = 0.5 * [cos(phi(i)) + cos(phi(i-1)); sin(phi(i)) + sin(phi(i-1))] * Vk(i);    
    [R_N,R_E] = Radii_of_curvature(L0); 
    Lambda0 = Lambda0 + V(2) * 0.5 / ((R_E + h0) * cos(L0));
    L0 = L0 + V(1) * 0.5 / (R_N + h0);           
    
    
    V_Nk = 1.7 * V(1) - 0.7 * V_N0;
    V_Ek = 1.7 * V(2) - 0.7 * V_E0;
    
    %phi
    Phi = eye(4);
    tao_s = 0.5;
    Phi(3,1) = tao_s/(R_N + h0);
    Phi(4,2) = tao_s/((R_E + h0) * cos(L0));
    
    %Q
    Q = zeros(4,4);
    S_DR = 0.2;
    Q(1,1) = S_DR * tao_s;
    Q(1,3) = 0.5 * S_DR * tao_s^2 / (R_N + h0);
    Q(2,2) = S_DR * tao_s;
    Q(2,4) = 0.5 * S_DR * tao_s^2 / (cos(L0)*(R_E + h0));
    Q(3,1) = 0.5 * S_DR * tao_s^2 / (R_N + h0);
    Q(3,3) = 1/3 * S_DR * tao_s^3 / (R_N + h0)^2;
    Q(4,2) = 0.5 * S_DR * tao_s^2 / (cos(L0)*(R_E + h0));
    Q(4,4) = 1/3 * S_DR * tao_s^3 / (cos(L0)*cos(L0)*(R_E + h0)^2);

    Xk = Phi * X0;
    Pk = Phi * P0 * Phi.'+ Q;

    %5 measurement matrix
    Hk = zeros(4,4);
    Hk(1,3) = -1;
    Hk(2,4) = -1;
    Hk(3,1) = -1;
    Hk(4,2) = -1;

    %6 measurement noise covariance    
    tao_Gr = 5;
    tao_Gv = 0.02;
    L00 = Motion(i,2) * deg_to_rad;
    [R_N,R_E] = Radii_of_curvature(L0);
    Rk = zeros(4,4);    
    Rk(1,1) = tao_Gr^2/(R_N + h0)^2;
    Rk(2,2) = tao_Gr^2/(cos(L00)*cos(L00)*(R_E + h0)^2);
    Rk(3,3) = tao_Gv^2;
    Rk(4,4) = tao_Gv^2;

    %7 Kalman gain matrix
    Kk = Pk * Hk.' * inv(Hk * Pk * Hk.' + Rk);

    %8 measurement innovation vector
    z = zeros(4,1);

    z(1) = Motion(i,2) * deg_to_rad - L0 + Xk(3);
    z(2) = Motion(i,3) * deg_to_rad - Lambda0 + Xk(4);
    z(3) = Motion(i,4) - V_Nk + Xk(1);
    z(4) = Motion(i,5) - V_Ek + Xk(2);

    %9 update the state
    Xk_new = Xk + Kk * z;
    %X0 = Xk_new;

    %10 update the error covariacne
    Pk_new = (eye(4) - Kk * Hk) * Pk;
    P0 = Pk_new;

    %11 correct the DR solution
    L0 = L0 - Xk_new(3);
    Lambda0 = Lambda0 - Xk_new(4);
    V_Nk = V_Nk - Xk_new(1);
    V_Ek = V_Ek - Xk_new(2);
    
    V_N0 = V_Nk;
    V_E0 = V_Ek;

    result_single = [0.5*(i-1), L0 * rad_to_deg, Lambda0 * rad_to_deg, V_Nk, V_Ek];
    Motion_3 = [Motion_3; result_single];
    
end

