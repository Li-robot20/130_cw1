function [latitude, longitude, height, v_north, v_east, v_down, outlier] = solver(time, Pseudoranges, Pseudorangesrates, omega_ie, c)
r_ea = [0;0;0];
v_ea = [0;0;0];

re = [];
ve = [];
ra = [];

outlier =[];

%Satellite number
index = [5, 6, 7, 9, 10, 11, 15, 30];
for i = 1:length(index)
    [r_ej, v_ej] = Satellite_position_and_velocity(time, index(i));
    r_aj = 0;
    re = [re; r_ej];
    ve = [ve; v_ej];
    ra = [ra; r_aj];
end

r_e = Pseudoranges(time/0.5+2, 2:length(index)+1);
r_e_dot = Pseudorangesrates(time/0.5+2, 2:length(index)+1);

rho_c = 100000;
x = [r_ea; rho_c];

rho_c_dot = 0;
x_dot = [v_ea; rho_c_dot];

distance = Inf;
while distance>0.1
    ua = [];
    ra_dot = [];
    z = [];
    z_dot = [];
    H = [];
    for i = 1:length(index)
        ra(i) = C2range(ra(i), re(i, :), r_ea, omega_ie, c);
        
        u_aj = unit_vector(ra(i), re(i, :), r_ea);
        ua = [ua; u_aj.'];
        
        r_aj_dot = predicted_range_rate(u_aj, ra(i), ve(i, :), re(i, :), r_ea, v_ea);
        ra_dot = [ra_dot; r_aj_dot];
        
        
        z_j = r_e(i) - ra(i) - x(4);
        z = [z; z_j];
        z_j_dot = r_e_dot(i) - r_aj_dot - x_dot(4);
        z_dot = [z_dot; z_j_dot];
        
        H_j = [-u_aj.', 1];
        H = [H; H_j];
    end
        
    x_new = x + inv(H.'*H)*H.'*z;
    distance = sqrt(sum((x_new - x).^ 2));
    x = x_new;
    r_ea = x(1:3);
    
    x_dot_new = x_dot + inv(H.'*H)*H.'*z_dot;
    x_dot = x_dot_new;
    v_ea = x_dot(1:3);
end

%outlier
v = (H*inv(H.'*H)*H.'-eye(8))*z;
Cv = (eye(8)-H*inv(H.'*H)*H.')*25;
for i=1:8
    if abs(v(i))>sqrt(Cv(i,i))*6
        outlier_temp = [time, i, v(i)];
        outlier = [outlier; outlier_temp];
    end
end
    
[L_b,lambda_b,h_b] = pv_ECEF_to_NED(x_new(1:3), [0;0;0]);
latitude = L_b*180/pi;
longitude = lambda_b*180/pi;
height = h_b;

[~,~,~,v_eb_n] = pv_ECEF_to_NED(r_ea, [x_dot_new(1);x_dot_new(2);x_dot_new(3)]);
v_north = v_eb_n(1);
v_east = v_eb_n(2);
v_down = v_eb_n(3);