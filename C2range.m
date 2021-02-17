function r = C2range(r_aj, r_ej, r_ea, omega_ie, c)
%omega_ie = 7.292115*10^-5;
%c = 299792458;
C = [1, omega_ie*r_aj/c, 0;
    -omega_ie*r_aj/c, 1, 0;
    0, 0, 1];
r = sqrt((C*r_ej.' - r_ea).' * (C*r_ej.' - r_ea));
end