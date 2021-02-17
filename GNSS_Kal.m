L_init = Motion(1,2) * deg_to_rad;
lambda_init = Motion(1,3) * deg_to_rad;
h_init = 0;
v_init =[Motion(1,4:5), 0].';
[r_init_e,v_init_e] = pv_NED_to_ECEF(L_init,lambda_init,h_init,v_init);

X0 = [r_init_e; v_init_e; 9901.1; 99.9];
P0 =  zeros(8);
P0(1,1) = 100;
P0(2,2) = 100;
P0(3,3) = 100;
P0(4,4) = 0.01;
P0(5,5) = 0.01;
P0(6,6) = 0.01;
P0(7,7) = 100;
P0(8,8) = 0.01;

i = 1;
Motion_2 = [];
while i<length(Pseudoranges)
    [Xk_new, Pk_new] = GNSS_Kalman(X0, P0, i, Pseudoranges(i+1,:), Pseudorangerates(i+1,:));
    X0 = Xk_new;
    P0 = Pk_new;
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(Xk_new(1:3), Xk_new(4:6));
    Lati = L_b * rad_to_deg;
    Longti = lambda_b * rad_to_deg;
    epoch = [0.5*(i-1),Lati,Longti,v_eb_n(1:2).'];
    Motion_2 = [Motion_2; epoch];
    i = i + 1;
end
    
