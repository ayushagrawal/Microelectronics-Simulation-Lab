close all;clear;

%% Defining constants
epsilon0 = 8.854 * 10^-12;epsilon_si = 11.68;   % S.I. units
epsilon_sio2 = 3.9;
nm = 1e-9;                                      % nano meter -> m
Eg = 1.1;                                       % eV
k = 1.38e-23;                                   % S.I. units
q = 1.6e-19;                                    % S.I. units
T = 300;                                        % K

% p-type
doping_na = 1e16;           %(cm^-3)
doping_na = doping_na*1e6;  %(m^-3)
n_i = 1.5e10;               %(cm^-3)
n_i = n_i*1e6;              %(m^-3)

np0 = n_i^2/doping_na;
pp0 = doping_na;
epsilon = epsilon0*epsilon_si;
epsilon_ox=epsilon_sio2*epsilon0;

% Assumed Surface Potential
Vs = 0.4;

%% Empirical Solution
W_empirical = sqrt(2*Vs*epsilon/(q*doping_na));

%% Numerical Solution

% Defining the region and parameters inside it
x = 0:0.1:1000;                     % 1um
x = x'*nm;                          % in m
theta = (x(2) - x(1));
N_A = doping_na*ones(size(x));      % p-type region

V = zeros(size(x));V(1) = Vs;

Error = 10; % Arbitrary High value
i = 0;
while Error > 10*eps
    i = i+1;
    d2V_by_dx2=(V(3:end) - 2*V(2:end-1) + V(1:end-2))/theta^2;
    rho = q*( - N_A(2:end-1) + ( -np0*exp(q*V(2:end-1)/(k*T)) + pp0*exp(-q*V(2:end-1)/(k*T))));
    R = d2V_by_dx2+rho/epsilon;

    Mj = 2/theta^2 + (q/epsilon)*((q/(k*T))*np0*exp(q*V(2:end-1)/(k*T)) + (q/(k*T))*pp0*exp(-q*V(2:end-1)/(k*T)));

    m = size(x,1);
    CM=sparse(1:m-2,1:m-2,Mj,m-2,m-2)...
        +sparse(1:m-2-1,2:m-2,(-1/theta^2)*ones(m-2-1,1),m-2,m-2)+...
        sparse(2:m-2,1:m-2-1,(-1/theta^2)*ones(m-2-1,1),m-2,m-2); 

    DV = CM\R;
    V(2:end-1)=V(2:end-1)+DV;
    Error=norm(DV,2)/sqrt(m);
end

%% Electron and Hole Densities
n_x = np0*exp(q*V/(k*T))*1e-6;   % (cm^-3)
h_x = pp0*exp(-q*V/(k*T))*1e-6;   % (cm^-3)

figure;
semilogy(x,n_x,'g','LineWidth',1.5);
hold on;
semilogy(x,h_x,'r','LineWidth',1.5);
xlabel('x (\mum) \rightarrow');
ylabel('ln(density(cm^-^3)) \rightarrow');
title('Electron & Hole Densities');
xlim([x(1),x(end)]);
plot([0.2e-6 0.2e-6],[n_x(end) h_x(end)],'--');
plot([0.22e-6 0.22e-6],[n_x(end) h_x(end)],'--');
plot([0.24e-6 0.24e-6],[n_x(end) h_x(end)],'--');
plot([0.26e-6 0.26e-6],[n_x(end) h_x(end)],'--');
plot([0.28e-6 0.28e-6],[n_x(end) h_x(end)],'--');
plot([0.3e-6 0.3e-6],[n_x(end) h_x(end)],'--');
text(0.3e-6,n_x(1),'\leftarrow x = 0.3\mum');
legend('n(x)','h(x)');
hold off;

%% Charge Density
% rho = 1e-6*q*( - N_A(1:end) + ( -np0*exp(q*V(1:end)/(k*T)) + pp0*exp(-q*V(1:end)/(k*T))));
% figure;
% plot(x,rho,'LineWidth',1.5);
% hold on;
% xlabel('x (\mum) \rightarrow');
% ylabel('\rho (C-cm^-^3) \rightarrow');
% title('Net Charge Density');
% xlim([x(1),x(end)]);
% hold off;