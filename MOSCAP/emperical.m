clear;close all;
% Defining constants
epsilon0 = 8.854 * 10^-12;epsilon_si = 11.68;   % S.I. units
epsilon_sio2 = 3.9;
nm = 1e-9;                                      % nano meter -> m
k = 1.38e-23;                                   % S.I. units
q = 1.6e-19;                                    % S.I. units
T = 300;                                        % K

% p-type
doping_na = 1e17;           %(cm^-3)
doping_na = doping_na*1e6;  %(m^-3)
n_i = 1.5e10;               %(cm^-3)
n_i = n_i*1e6;              %(m^-3)


np0 = n_i^2/doping_na;
pp0 = doping_na;
epsilon = epsilon0*epsilon_si;
epsilon_ox=epsilon_sio2*epsilon0;

Vt = k*T/q;
phi_f = Vt*log(pp0/n_i);
Vs = -0.25:0.01:1;

expr = cosh((Vs-phi_f)/Vt) + (Vs/Vt)*sinh(phi_f/Vt) - cosh(phi_f/Vt);
Q = 2*sign(Vs).*sqrt(q*epsilon*n_i*Vt*expr);

semilogy(Vs,abs(Q));