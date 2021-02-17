function r_aj_dot = predicted_range_rate(u_aj,r_aj,v_ej,r_ej,r_ea,v_ea)
omega_ie = 7.292115*10^-5;
A = [0, -omega_ie, 0;
    omega_ie, 0, 0;
    0, 0, 0];
c = 299792458;
C = [1, omega_ie*r_aj/c, 0;
    -omega_ie*r_aj/c, 1, 0;
    0, 0, 1];
r_aj_dot = u_aj.'*(C*(v_ej.' + A*r_ej.') - (v_ea + A*r_ea));
