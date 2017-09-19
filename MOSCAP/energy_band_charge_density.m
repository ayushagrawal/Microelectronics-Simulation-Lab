close all;
%% NUMERICALLY SOLVING POISSON FOR N-MOSCAP

% Defining constants
epsilon0 = 8.854 * 10^-12;epsilon_si = 11.68;   % S.I. units
nm = 1e-9;                                      % nano meter -> m
Eg = 1.1;                                       % eV
k = 1.38e-23;                                   % S.I. units
q = 1.6e-19;                                    % S.I. units
T = 300;                                        % K

% p-type
doping_na = 1e18;           %(cm^-3)
doping_na = doping_na*1e6;  %(m^-3)
n_i = 1.5e10;               %(cm^-3)
n_i = n_i*1e6;              %(m^-3)


np0 = n_i^2/doping_na;
pp0 = doping_na;
epsilon = epsilon0*epsilon_si;


% Defining the region and parameters inside it
x = 0:0.1:100;                       % 1um
x = x'*nm;                          % in m
theta = (x(2) - x(1));
N_A = doping_na*ones(size(x));      % p-type region

% Assumed Surface Potential
Vs = -0.15;

V = zeros(size(x));V(1) = Vs;

figure;
plot(V);
hold on;
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
    plot(V);
end
title('Numerical Convergence');
xlabel('x (m) \rightarrow');
ylabel('V (volt) \rightarrow');

%% Plotting the final Potential Profile

% figure;
% plot(x,V,'LineWidth',1.5);
% hold on;
% xlabel('x (m) \rightarrow');
% ylabel('V (volt) \rightarrow');
% title('Potential Profile of the device');
% xlim([x(1),x(end)]);
% ylim([-abs(Vs),abs(Vs)]);
% line([0 0],[-abs(Vs) abs(Vs)],'Color','red','LineStyle','--');
% text(-0.5e-7,0,'n-type','HorizontalAlignment','center');
% text(0.5e-7,0,'p-type','HorizontalAlignment','center');

%% Plotting the Energy Band Diagram
Ev = -V-ones(size(V))*Eg;
x=x*1e6;

figure;
plot(x,-V,'LineWidth',2);
hold on;
plot(x,Ev,'LineWidth',2);
xlabel('x (\mum) \rightarrow');
ylabel('Energy (eV) \rightarrow');
title('Energy Band Diagram');
xlim([x(1),x(end)]);
hold off;

%% Electron and Hole Densities
% n_x = np0*exp(q*V/(k*T))*1e-6;   % (cm^-3)
% h_x = pp0*exp(-q*V/(k*T))*1e-6;   % (cm^-3)
% 
% figure;
% semilogy(x,n_x,'g','LineWidth',1.5);
% hold on;
% semilogy(x,h_x,'r','LineWidth',1.5);
% xlabel('x (\mum) \rightarrow');
% ylabel('ln(density(cm^-^3)) \rightarrow');
% title('Electron & Hole Densities');
% xlim([x(1),x(end)]);
% legend('n(x)','h(x)');
% hold off;

%% Charge Density
rho = q*( - N_A(1:end) + ( -np0*exp(q*V(1:end)/(k*T)) + pp0*exp(-q*V(1:end)/(k*T))));
figure;
plot(x,rho*1e-6,'g','LineWidth',1.5);
hold on;
xlabel('x (\mum) \rightarrow');
ylabel('\rho (C-cm^-^3) \rightarrow');
title('Net Charge Density');
xlim([x(1),x(end)]);
hold off;