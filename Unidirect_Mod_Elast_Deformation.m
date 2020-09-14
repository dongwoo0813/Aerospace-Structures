clc
clear all
close all

% This is the MATLAB file which computes the linear modulus of elasticity
% and ultimate strength of unidirectional fiber lamina in the longitudinal and 
% transverse direction WITH DEFORMATION. This file assumes circular fibers.

% This MATLAB file is done as the assignment in AE418 course at
% Embry-Riddle Aeronautical University by Edward Jung (Dongwoo Han).
% Please do not copy and use this file without any notice.

% This file takes the values such as Elastic Modulus and Shear Modulus as
% inputs to compute linear modulus of elasticity and ultimate strength.

answer = input('In a quotation (''), type in the type of the composite such as carbon epoxy : ' );
t = input('\nType in the value of thickness, t (in inch): ');
E_F = input('\nType in the value of Fiber Modulus, E_F (in ksi): ');  
sigma_FU = input('\nType in the value of Ultimate Fiber Strength, sigma_FU (in ksi): ');
E_M = input('\nType in the value of Matrix Modulus, E_M (in ksi): ');
sigma_M_e = input('\nType in the value of Matrix Strengh at Fiber Breakage sigma_M_epsilon (in ksi): ');
sigma_MU = input('\nType in the value of Ultimate Matrix Strength, sigma_MU (in ksi): ');
V_F = input('\nType in the value of Fiber Volume Fraction, V_F: ');
v_LT = input('\nType in the value of Major Poissons ratio for Transverse Strain caused by Longitudinal Stress, v_LT: ');
v_TL = input('\nType in the value of Minor Poissons ratio for Longitudinal Strain caused by Transverse Stress, v_TL: ');
G_LT = input('\nType in the value of Shear Modulus with respect to axes of symmetry, G_LT (in ksi): ');

L = input('\nType in the value of Length (in): ');
w = input('\nType in the value of Width (in): ');
P_T = input('\nType in the value of Transverse Force (kip): ');
P_L = input('\nType in the value of Longitudinal Force (kip): ');
SF_LT = input('\nType in the value of Shear Force (kip): ');

E_L = E_F*V_F + E_M*(1-V_F);                    % Calculates the Longitudinal Modulus.
xi = 2;                                         % Filler Geometry Parameter, 2 for circular fiber, 2*l/d for calculation of the longitudinal modulus
eta = ((E_F/E_M)-1)/((E_F/E_M) + xi);           % Calculates the Fiber Shape Factor.
E_T = E_M*(1 + xi*eta*V_F)/(1 - eta*V_F);       % Calculates the Transverse Modulus.
sigma_LU = sigma_FU*V_F + sigma_M_e*(1-V_F);    % Calculates the Longitudinal Strength.

% For Transverse Strength, it is assumed that it is directly proportional.
% However, it uses stress concentration factor to account for holes in
% matrix due to fibers.
S = (1- (V_F*(1-(E_M/E_F))))/(1- (((4*V_F/pi)^0.5)*(1-(E_M/E_F)))); % Stress Concentration Factor
sigma_TU = sigma_MU/S;                          % Calculates the Transverse Strength.

% Compliance Matrix
S11 = 1/E_L;
S12 = -v_LT/E_L;
% S21 = -v_TL/E_T;
S22 = 1/E_T;
S66 = 1/G_LT;
S_Matrix = [S11 S12 0; S12 S22 0; 0 0 S66];

sigma_L = P_L/(t*w);    % Longitudinal Stress
sigma_T = P_T/(t*L);    % Tranverse Stress
tau_LT = SF_LT/(w*L);   % Shear Stress

ans = S_Matrix*[sigma_L; sigma_T; tau_LT];
epsilon_L = ans(1,1); 
epsilon_T = ans(2,1);
gamma_LT = ans(3,1); 

delta_L = epsilon_L*L; % Longitudinal Deformation
delta_T = epsilon_T*w; % Transverse Deformation

fid = fopen('Lamina Output Files with Deformation.txt', 'w');

fprintf(fid, 'Composite Type: %-20s \r\n \r\n', answer);
fprintf(fid, 'INPUT DATA:\r\n');
fprintf(fid, 'Fiber Modulus (ksi) = %5.2f \r\n', E_F);
fprintf(fid, 'Matrix Modulus (ksi) = %5.2f \r\n', E_M);
fprintf(fid, 'Fiber Volume Fraction = %5.2f \r\n', V_F);
fprintf(fid, 'Ultimte Fiber Strength (ksi) = %5.2f \r\n', sigma_FU);
fprintf(fid, 'Matrix Strength at fiber breakage(ksi) = %5.2f \r\n', sigma_M_e);
fprintf(fid, 'Ultimate Matrix Strength (ksi) = %5.2f \r\n', sigma_MU);
fprintf(fid, 'Poissons LT = %5.2f \r\n', v_LT);
fprintf(fid, 'Poissons TL = %5.2f \r\n', v_TL);
fprintf(fid, 'Shear Modulus LT (ksi) = %5.2f \r\n', G_LT);
fprintf(fid, 'Lamina thickness (in.) = %6.4f \r\n', t);
fprintf(fid, 'Longitudinal force (kips) = %5.2f \r\n',P_L);
fprintf(fid, 'Transverse force (kips) = %5.2f \r\n',P_T);
fprintf(fid, 'Lamina Length (in.) = %5.2f \r\n', L);
fprintf(fid, 'Lamina Width (in.) = %5.2f \r\n\r\n',w);


fprintf(fid, 'OUTPUT DATA: \r\n\r\n');
fprintf(fid, 'S Matrix = \r\n');
fprintf(fid, '(ksi^-1) \r\n');
fprintf(fid, '%+10f  ',S11);
fprintf(fid, '%+10f  ', S12);
fprintf(fid, '%+10f\r\n', 0);
fprintf(fid, '%+10f  ',S21);
fprintf(fid, '%+10f  ', S22);
fprintf(fid, '%+10f\r\n', 0);
fprintf(fid, '%+10f  ',0);
fprintf(fid, '%+10f    ', 0);
fprintf(fid, '%+10f   \r\n\r\n', S66);



fprintf(fid, 'sigma_L (ksi) = %5.2f \r\n', sigma_L); 
fprintf(fid, 'sigma_T (ksi) = %5.2f \r\n', sigma_T);
fprintf(fid, 'tau_LT (ksi) = %5.2f \r\n\n', tau_LT);
fprintf(fid, 'epsilon_L = %6.5f \r\n', epsilon_L);
fprintf(fid, 'epsilon_T = %6.5f \r\n', epsilon_T);
fprintf(fid, 'gamma_LT = %6.5f \r\n\n', gamma_LT);
fprintf(fid, 'Longitudinal Deformation (in.) = %6.5f \r\n', delta_L);
fprintf(fid, 'Transverse Deformation (in.)= %6.5f \r', delta_T);
fclose(fid);
