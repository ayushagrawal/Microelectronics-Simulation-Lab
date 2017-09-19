close all;clear;
%% NUMERICALLY SOLVING POISSON FOR N-MOSCAP

% Defining constants
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
trap_den = [0,1e11,1e12]';  %(cm^-2 eV^-1)
trap_den = trap_den*1e4;    %(m^-2 eV^-1)

np0 = n_i^2/doping_na;
pp0 = doping_na;
epsilon = epsilon0*epsilon_si;
epsilon_ox=epsilon_sio2*epsilon0;

% Defining the region and parameters inside it
x = 0:.1:1000;                        % (nm)
x = x'*nm;                          % (m)
theta = (x(2) - x(1));
N_A = doping_na*ones(size(x));      % p-type region

% Constructing test and result vectors
Vs = [-0.2:0.001:0.9]';
length = size(Vs,1);
V_applied = zeros(length,1);
charge = zeros(length,1);
capacitor = zeros(length,1);

figure;
hold on;

for num=1:size(trap_den,1)
    for i=1:length

        V = zeros(size(x));
        V(1) = Vs(i);   % Defining the Boundary Conditions

        Error = 10;     % Arbitrary High value
        while Error > 10*eps
            %i = i+1;
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
            %plot(V);
        end

        % Calculating Actual Trap Density (Surface Trap Density)
        Ev = -Vs(i)*q;
        Ef = ((Eg*q)/2) - (k*T*log(pp0/n_i));
        step = (Ef-Ev)/1000;
        Et = [Ev:step:Ef]';
        gD = 2;
        F_SD = 1./(1+gD*exp((Et-Ef)/(k*T)));
        net_trap = -trap_den(num)*sum(F_SD,1)*(step/q);        % (m^-2)
        trap_sigma = net_trap*q;                               % (Cm^-2)

        % This portion is modified due to presence of interface charges
        E_si = -(V(2)-V(1))/theta;
        E_sio2 = ((E_si*epsilon)-trap_sigma)/epsilon_ox;
        thickness = 10*nm;
        V_applied(i) = V(1) + E_sio2*thickness;
        charge(i) = -E_sio2*epsilon_ox;
    end
    charge=charge*1e-4;
    for i=1:length-1
        capacitor(i) = (charge(i+1)-charge(i))/(V_applied(i+1) - V_applied(i));
    end
    capacitor(end) = capacitor(end-1);
    
    % Plotting
    plot(V_applied,-capacitor,'LineWidth',1.5);
end
%% Plotting Charge Density
% figure;
% semilogy(Vs,abs(charge),'LineWidth',1.5);
% title('Charge Density Vs Surface Potential');
% xlabel('V_s \rightarrow');
% ylabel('\rho (cm^{-3}) \rightarrow');
% 
% figure;
% semilogy(V_applied,abs(charge),'LineWidth',1.5);
% title('Charge Density Vs Applied Potential');
% xlabel('V_{applied} \rightarrow');
% ylabel('\rho (cm^{-3}) \rightarrow');
% xlim([-2 3]);

%% Plotting the LFCV Curve
title('LFCV Curve');
xlabel('V_{applied} (V) \rightarrow');
ylabel('C (F/cm^2)\rightarrow');
%xlim([-2 3]);
%ylim([0 1.2]);
legend('D_t = 0 cm^{-2}eV^{-1}','D_t = 10^{11} cm^{-2}eV^{-1}','D_t = 10^{12} cm^{-2}eV^{-1}');