% read data from Dead reckoning file
gyro = DR( : ,6);
magne = DR( : ,7);

% initialise state vector and error covariance
x_kmins1 = zeros(2,1);
P_kmins1 = [0.01^2, 0;
            0, 0.001^2];

% initialise gyro solution
gyro_init = magne(1);
Heading = [gyro_init];

for i = 1:length(DR)-1
    % Integrate Gyro and magnetometer to obtain the heading solution
    [gyro_sol, x_k_plus, P_k_plus] = Gyro_magne(gyro, magne, i, gyro_init, x_kmins1, P_kmins1);
    Heading = [Heading; gyro_sol];    
    gyro_int = gyro_sol;
    x_kmins1 = x_k_plus;
    P_kmins1 = P_k_plus;
end

