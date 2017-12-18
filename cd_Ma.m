function cd = cd_Ma(Ma)

load('c_D_vs_Ma')
c_D_data = c_D_vs_Ma(:,2);
Ma_data = c_D_vs_Ma(:,1);
cd = interp1(Ma_data,c_D_data,Ma,'pchip');

return