num_total = length(Pseudoranges);
Motion = zeros(length(Pseudoranges)-1, 7);
n = 0;
outlier =[];
outlier_second = [];
while n<=length(Pseudoranges)-2
    Motion(n+1,1) = Pseudoranges(n+2, 1);
    [Motion(n+1,2), Motion(n+1,3), Motion(n+1,7), Motion(n+1,4), Motion(n+1,5), v_down, outlier_temp] = solver(0.5*n, Pseudoranges, Pseudorangerates, omega_ie, c);
    outlier = [outlier; outlier_temp];
    n = n+1;
end

n = 0;
while n<=length(Pseudoranges)-2
    Motion(n+1,1) = Pseudoranges(n+2, 1);
    [Motion(n+1,2), Motion(n+1,3), Motion(n+1,7), Motion(n+1,4), Motion(n+1,5), v_down, outlier_temp] = solver_no_outlier(0.5*n, Pseudoranges, Pseudorangerates, omega_ie, c, outlier);
    outlier_second = [outlier_second; outlier_temp];
    n = n+1;
end

%heading 
Motion(1,6) = DR(1,7);
gyro = DR( : ,6);
magne = DR( : ,7);
x_kmins1 = zeros(2,1);
P_kmins1 = [0.01^2, 0;
            0, 0.001^2];
gyro_init = magne(1);
for i = 1:length(DR)-1
    [gyro_sol, x_k_plus, P_k_plus] = Gyro_magne(gyro, magne, i, gyro_init, x_kmins1, P_kmins1);
    Motion(i+1,6) = gyro_sol;
    gyro_int = gyro_sol;
    x_kmins1 = x_k_plus;
    P_kmins1 = P_k_plus;
end


