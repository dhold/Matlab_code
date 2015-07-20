function [out_x,out_y,out_z] = full_cont_3rd_order2(mu,R,jj,pol,kk)
%calculates the contributions to each output polarization direction given
%appropriate dipole function, positions R, and which site each interaction
%is taking place on in "jj" and 3 beam wavevectors in $kk$
%output wavevector k_4 is sum(kk)
k_1 = kk(1,:); k_2 = kk(2,:); k_3 = kk(3,:);  %kk is k_j along rows
k_4 = sum(kk); %auto sums along columns as desired

%principle average
   xxxx_av = (1/15)*(dot(mu(jj(1),:),mu(jj(2),:))*dot(mu(jj(3),:),mu(jj(4),:))+...
                dot(mu(jj(1),:),mu(jj(3),:))*dot(mu(jj(2),:),mu(jj(4),:))+...
                dot(mu(jj(2),:),mu(jj(3),:))*dot(mu(jj(1),:),mu(jj(4),:))); %also yyyy etc
        %only non zero 4 comp averages with last index x
                  
        yyxx_av = (1/30)*(4*dot(mu(jj(1),:),mu(jj(2),:))*dot(mu(jj(3),:),mu(jj(4),:))-...
                dot(mu(jj(1),:),mu(jj(3),:))*dot(mu(jj(2),:),mu(jj(4),:))-...
                dot(mu(jj(2),:),mu(jj(3),:))*dot(mu(jj(1),:),mu(jj(4),:)));    %also xxzz etc
        xyyx_av = (1/30)*(-dot(mu(jj(1),:),mu(jj(2),:))*dot(mu(jj(3),:),mu(jj(4),:))+4*...
                dot(mu(jj(1),:),mu(jj(3),:))*dot(mu(jj(2),:),mu(jj(4),:))-...
                dot(mu(jj(2),:),mu(jj(3),:))*dot(mu(jj(1),:),mu(jj(4),:)));               
        yxyx_av = (1/30)*(-dot(mu(jj(1),:),mu(jj(2),:))*dot(mu(jj(3),:),mu(jj(4),:))-...
                dot(mu(jj(1),:),mu(jj(3),:))*dot(mu(jj(2),:),mu(jj(4),:))+4*...
                dot(mu(jj(2),:),mu(jj(3),:))*dot(mu(jj(1),:),mu(jj(4),:)));           
            
            %now xxxyz averages, also yyyzx and zzzxy, this is very general,
        %possibly overkill but takes little computational time so is
        %included, works 9or should at least) for any input beam 
        %configuration and polarizations
        % diff one for each
        
        xxxyz_av(1) = (1i/30)*...
            (dot(mu(jj(1),:),mu(jj(2),:))*dot(mu(jj(3),:),cross(mu(jj(4),:),R(jj(1),:)))+...
             dot(mu(jj(1),:),mu(jj(3),:))*dot(mu(jj(2),:),cross(mu(jj(4),:),R(jj(1),:)))+...
             dot(mu(jj(2),:),mu(jj(3),:))*dot(mu(jj(1),:),cross(mu(jj(4),:),R(jj(1),:))));
        xxxyz_av(2) = (1i/30)*...
            (dot(mu(jj(1),:),mu(jj(2),:))*dot(mu(jj(3),:),cross(mu(jj(4),:),R(jj(2),:)))+...
             dot(mu(jj(1),:),mu(jj(3),:))*dot(mu(jj(2),:),cross(mu(jj(4),:),R(jj(2),:)))+...
             dot(mu(jj(2),:),mu(jj(3),:))*dot(mu(jj(1),:),cross(mu(jj(4),:),R(jj(2),:))));        
        xxxyz_av(3) = (1i/30)*...
            (dot(mu(jj(1),:),mu(jj(2),:))*dot(mu(jj(3),:),cross(mu(jj(4),:),R(jj(3),:)))+...
             dot(mu(jj(1),:),mu(jj(3),:))*dot(mu(jj(2),:),cross(mu(jj(4),:),R(jj(3),:)))+...
             dot(mu(jj(2),:),mu(jj(3),:))*dot(mu(jj(1),:),cross(mu(jj(4),:),R(jj(3),:))));        
        xxxyz_av(4) = -(1i/30)*...
            (dot(mu(jj(1),:),mu(jj(2),:))*dot(mu(jj(3),:),cross(mu(jj(4),:),R(jj(4),:)))+...
             dot(mu(jj(1),:),mu(jj(3),:))*dot(mu(jj(2),:),cross(mu(jj(4),:),R(jj(4),:)))+...
             dot(mu(jj(2),:),mu(jj(3),:))*dot(mu(jj(1),:),cross(mu(jj(4),:),R(jj(4),:))));
%-sign for the last interaction just from the derivation

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