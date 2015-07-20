function [Yori] = orientation_prefactor(dip_ang,beam_ang)

%calculates the oritentation prefactor given a set of angles between the
%electric dipole, and between the beams.  Also if given 4 arguments
%will calculate the extra terms related to the magnetic dipole moments and
%beam prefactors
%
%This code works with sets angles between different dipoles and beams
% dip_ang should be of the form [theta_12,theta_13,theta_14,theta_23...
% theta_24, theta_34; .... ; ...] and similar for the beam angles.
%
% To include the magnetic dipole effects just pass four different terms
% with the first el dipole replaced with mag dipole angles (relative to the
% other electric dipoles) and mag field of the first beam relative to the
% other beams then second etc and sum over these four terms to reduce the
% output


       tmp = cos(dip_ang);
       g1 = tmp(:,6).*tmp(:,1); g2 = tmp(:,5).*tmp(:,2); g3 = tmp(:,3).*tmp(:,4);
       
       tmp2 = cos(beam_ang);
    sigma = tmp2(:,6)*tmp2(:,1) + tmp2(:,5)*tmp2(:,2) + tmp2(:,3)*tmp2(:,4);
    %5*cos(theta_43)*cos(theta_21) - sigma
    f1 = 5*tmp2(:,6)*tmp2(:,1) -sigma; f2 = 5*tmp2(:,5)*tmp2(:,2) - sigma;
    f3 = 5*tmp2(:,3)*tmp2(:,4) - sigma;

    Yori = f1*g1.' + f2*g2.' + f3*g3.';
    
    
end

