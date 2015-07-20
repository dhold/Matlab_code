function [out_x,out_y,out_z] = full_cont_3rd_order(xxxx_av,yyxx_av,...
                xyyx_av,yxyx_av ,xxxyz_av,pol,k_1,k_2,k_3,k_4)
%calculates the contributions to each output polarization direction given
%appropriate dipole averages (fed to function) and 3 beam wavevectors and
%output wavevector k_4

out_x = ...
(xxxx_av*pol{1}(1)*pol{2}(1)*pol{3}(1)+... 
yyxx_av*(pol{1}(2)*pol{2}(2)*pol{3}(1) +pol{1}(3)*pol{2}(3)*pol{3}(1))+...
xyyx_av*(pol{1}(1)*pol{2}(2)*pol{3}(2) +pol{1}(1)*pol{2}(3)*pol{3}(3))+...
yxyx_av*(pol{1}(2)*pol{2}(1)*pol{3}(2) +pol{1}(3)*pol{2}(1)*pol{3}(3))+...
(xxxyz_av(1)*k_1(2)+xxxyz_av(2)*k_2(2)+xxxyz_av(3)*k_3(2)+xxxyz_av(4)*k_4(2))...
  *pol{1}(3)*pol{2}(3)*pol{3}(3));

out_y = ...
(xxxx_av*pol{1}(2)*pol{2}(2)*pol{3}(2)+...
yyxx_av*(pol{1}(1)*pol{2}(1)*pol{3}(2) +pol{1}(3)*pol{2}(3)*pol{3}(2))+...
xyyx_av*(pol{1}(2)*pol{2}(1)*pol{3}(1) +pol{1}(2)*pol{2}(3)*pol{3}(3))+...
yxyx_av*(pol{1}(1)*pol{2}(2)*pol{3}(1) +pol{1}(3)*pol{2}(2)*pol{3}(3))+...
(xxxyz_av(1)*k_1(3)+xxxyz_av(2)*k_2(3)+xxxyz_av(3)*k_3(3)+xxxyz_av(4)*k_4(3))...
  *pol{1}(1)*pol{2}(1)*pol{3}(1));

out_z = ...
(xxxx_av*pol{1}(3)*pol{2}(3)*pol{3}(3)+...
yyxx_av*(pol{1}(2)*pol{2}(2)*pol{3}(3) +pol{1}(1)*pol{2}(1)*pol{3}(3))+...
xyyx_av*(pol{1}(3)*pol{2}(2)*pol{3}(2) +pol{1}(3)*pol{2}(1)*pol{3}(1))+...
yxyx_av*(pol{1}(2)*pol{2}(3)*pol{3}(2) +pol{1}(1)*pol{2}(3)*pol{3}(1))+...
(xxxyz_av(1)*k_1(1)+xxxyz_av(2)*k_2(1)+xxxyz_av(3)*k_3(1)+xxxyz_av(4)*k_4(1))...
  *pol{1}(2)*pol{2}(2)*pol{3}(2));