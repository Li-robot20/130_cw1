function [latitude, longitude, height, v_north, v_east, v_down, outlier] = GNSS_LS_solver_remove_outlier(time, sate_index, Pseudoranges, Pseudorangesrates, detected_outlier)
% Inputs:
%   time                 Current simulation time(s)
%   sate_index           Satellite number
%   Pseudoranges         Corresponding Pseudoranges data
%   Pseudorangesrates    Corresponding Pseudorangesrates data
%   detected_outlier     detected outlier information, including time and Satellite number
% Outputs:
%   L_b                 latitude (rad)
%   lambda_b            longitude (rad)
%   h_b                 height (m)
%   v_north             north velocity of body frame w.r.t. ECEF frame
%   v_east              east velocity of body frame w.r.t. ECEF frame
%   v_down              down velocity of body frame w.r.t. ECEF frame    
%   outlier             detected outlier information, including time and Satellite number 

% init position and velocity
r_ea = [0;0;0];
v_ea = [0;0;0];

%ECEF satellite position and velocity
index = sate_index;
re = [];
ve = [];
ra = [];
for i = 1:length(index)
    [r_ej, v_ej] = Satellite_position_and_velocity(time, index(i));
    r_aj = 0;
    re = [re; r_ej];
    ve = [ve; v_ej];
    ra = [ra; r_aj];
end

% ECEF satellite position and velocity from Pseudo data
r_e = Pseudoranges(time/0.5+2, 2:length(index)+1);
r_e_dot = Pseudorangesrates(time/0.5+2, 2:length(index)+1);

% predicted state vector
rho_c = 100000;
x = [r_ea; rho_c];
rho_c_dot = 200;
x_dot = [v_ea; rho_c_dot];

% Iterate until difference from previous iteration is small
distance = Inf;
while distance>0.01
    ua = [];
    ra_dot = [];
    z = [];
    z_dot = [];
    H = [];
    for i = 1:length(index)
        % Predict the ranges from the position       
        ra(i) = predict_raj(ra(i), re(i, :), r_ea);
        
        % Compute the line-of-sight unit vector
        u_aj = unit_vector(ra(i), re(i, :), r_ea);
        ua = [ua; u_aj.'];
        
        % Predict the range rates
        r_aj_dot = predict_raj_dot(u_aj, ra(i), ve(i, :), re(i, :), r_ea, v_ea);
        ra_dot = [ra_dot; r_aj_dot];
        
        % measurement innovation vector 
        if ismember(time, detected_outlier(:,1)) & ismember(i, detected_outlier(:,2))
            % calculate without the measurement which is detected as outlier
            z_j = 0;
        else
            z_j = r_e(i) - ra(i) - x(4);
        end
        z = [z; z_j];
        z_j_dot = r_e_dot(i) - r_aj_dot - x_dot(4);
        z_dot = [z_dot; z_j_dot];
        
        % measurement matrix
        H_j = [-u_aj.', 1];
        H = [H; H_j];
    end
    
    % position and receiver clock drift solution
    x_new = x + inv(H.'*H)*H.'*z;
    distance = sqrt(sum((x_new - x).^ 2));
    x = x_new;
    r_ea = x(1:3);
    
    % velocity and receiver clock drift solution
    x_dot_new = x_dot + inv(H.'*H)*H.'*z_dot;
    x_dot = x_dot_new;
    v_ea = x_dot(1:3);
end

% Detect outlier
outlier =[];

% Compute the residuals vector
v = (H*inv(H.'*H)*H.'-eye(8))*z;

% Compute the residuals covariance matrix
Cv = (eye(8)-H*inv(H.'*H)*H.')*25;

% Compare normalised residuals with a threshold
for i=1:length(index)
    if abs(v(i))>sqrt(Cv(i,i))*6
        outlier_temp = [time, i];
        outlier = [outlier; outlier_temp];
    end
end
    

% Converts Cartesian  to curvilinear position and velocity resolving axes from ECEF to NED
[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_new(1:3), x_dot_new(1:3));
latitude = L_b*180/pi;
longitude = lambda_b*180/pi;
height = h_b;
v_north = v_eb_n(1);
v_east = v_eb_n(2);
v_down = v_eb_n(3);