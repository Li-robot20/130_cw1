% heading and forward speed solution
phi = Heading * deg_to_rad;
Vk = (DR(:, 2) + DR(:, 3) + DR(:, 4) + DR(:, 5))/4;

result_1 = [0, Motion_1(1,2:5), Heading(1)];
Motion_3 = [result_1];

% initialise position and velocity
L_k = Motion_1(1,2) * deg_to_rad;
Lambda_k = Motion_1(1,3) * deg_to_rad;
V_Nk = Motion_1(1,4);
V_Ek = Motion_1(1,5);
h_k = Heights(1);

%Dead reckoning
for i = 2:length(DR)
    h_k = Heights(i);
    V = 0.5 * [cos(phi(i)) + cos(phi(i-1)); sin(phi(i)) + sin(phi(i-1))] * Vk(i);    
    [R_N,R_E] = Radii_of_curvature(L_k); 
    
    % Dead reckoning solution for position
    Lambda_k = Lambda_k + V(2) * 0.5 / ((R_E + h_k) * cos(L_k));
    L_k = L_k + V(1) * 0.5 / (R_N + h_k);           
    
    % Dead reckoning solution for velocity
    V_Nk = 1.7 * V(1) - 0.7 * V_Nk;
    V_Ek = 1.7 * V(2) - 0.7 * V_Ek;

    % store results at each epoch into motion file
    result_single = [0.5*(i-1), L_k * rad_to_deg, Lambda_k * rad_to_deg, V_Nk, V_Ek, Heading(i)];
    Motion_3 = [Motion_3; result_single];
    
end