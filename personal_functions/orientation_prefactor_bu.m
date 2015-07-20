function [Yori,Y2] = orientation_prefactor(dip_ang,beam_ang,mdim_ang,bmag_ang)

%calculates the oritentation prefactor given a set of angles between the
%electric dipole, and between the beams.  Also if given 4 arguments
%will calculate the extra terms related to the magnetic dipole moments and
%beam prefactors
%
%This code works with sets angles between different dipoles and beams
% dip_ang should be of the form [theta_12,theta_13,theta_14,theta_23...
% theta_24, theta_34; .... ; ...] and similar for the beam angles.
% the mdim_ang and bmag_ang are the angles between the magnetic dipole of
% one site and the electric dipole of another (and similar for the beam)
% mdim_ang = [theta'_12,theta'_13,theta'_14;...] and
% other terms (theta'_21,theta'_23,theta'_24 etc) in the 3rd dimension


       tmp = cos(dip_ang);
       g1 = tmp(:,6).*tmp(:,1); g2 = tmp(:,5).*tmp(:,2); g3 = tmp(:,3).*tmp(:,4);
       
       tmp2 = cos(beam_ang);

    sigma = tmp2(:,6)*tmp2(:,1) + tmp2(:,5)*tmp2(:,2) + tmp2(:,3)*tmp2(:,4);
    %5*cos(theta_43)*cos(theta_21) - sigma
    f1 = 5*tmp2(:,6)*tmp2(:,1) -sigma; f2 = 5*tmp2(:,5)*tmp2(:,2) - sigma;
    f3 = 5*tmp2(:,3)*tmp2(:,4) - sigma;

    Yori = f1*g1.' + f2*g2.' + f3*g3.';
    
if nargin > 3 %terms where one electric dipole interaction is replaced with
    %a magnetic dipole interaction.  Come to think of it I can just call
    %the function with the angles between the jth interaction term and all
    %the others replaced with 
    
    tmp3 = cos(mdim_ang); tmp4 = cos(bmag_ang); %n1 by n2 by 4
    g11 = tmp(:,6).*tmp3(:,1,1); g21 = tmp(:,5).*tmp3(:,2,1); g31 = tmp3(:,3,1).*tmp(:,4);
    g12 = tmp(:,6).*tmp3(:,1,2); g22 = tmp(:,5).*tmp3(:,3,2); g31 = tmp3(:,3,1).*tmp(:,4);
    g13 = tmp(:,1).*tmp3(:,3,3); g23 = tmp(:,5).*tmp3(:,2,1); g31 = tmp3(:,3,1).*tmp(:,4);
    g14 = tmp(:,1).*tmp3(:,3,4); g24 = tmp(:,5).*tmp3(:,2,1); g31 = tmp3(:,3,1).*tmp(:,4);    


    Y2 = f11*g11.' + f21*g21.' + f31*g31.'+...
        f12*g12.' + f22*g22.' + f32*g32.'+...
        f13*g13.' + f23*g23.' + f33*g33.'+...
        f14*g14.' + f24*g24.' + f34*g34.';
else
    Y2 = [];
end
    
    
end

