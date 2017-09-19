close all;clear;
%% NUMERICALLY SOLVING POISSON FOR N-MOSCAP

% Defining constants
epsilon0 = 8.854 * 10^-12;epsilon_si = 11.68;   % S.I. units
epsilon_sio2 = 3.9;
nm = 1e-9;                                      % nano meter -> m
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

% Defining the region and parameters inside it
x = 0:.1:1000;                        % (nm)
x = x'*nm;                          % (m)
theta = (x(2) - x(1));
N_A = doping_na*ones(size(x));      % p-type region

% Constructing test and result vectors
Vs = [-0.3:0.001:1]';
length = size(Vs,1);

V_applied1 = zeros(length,1);
V_applied2 = zeros(length,1);
V_applied3 = zeros(length,1);
V_applied4 = zeros(length,1);
V_applied5 = zeros(length,1);

charge = zeros(length,1);

capacitor1 = zeros(length,1);
capacitor2 = zeros(length,1);
capacitor3 = zeros(length,1);
capacitor4 = zeros(length,1);
capacitor5 = zeros(length,1);

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
    
    %rho = q*( - N_A(1:end) + ( -np0*exp(q*V(1:end)/(k*T)) + pp0*exp(-q*V(1:end)/(k*T))));
    E_si = -(V(2)-V(1))/theta;
    E_sio2 = (E_si*epsilon)/epsilon_ox;
    
    % ############################################
    thickness1 = 10*nm;                         %#
    V_applied1(i) = V(1) + E_sio2*thickness1;   %#
    thickness2 = 20*nm;                         %#
    V_applied2(i) = V(1) + E_sio2*thickness2;   %#
    thickness3 = 40*nm;                         %#
    V_applied3(i) = V(1) + E_sio2*thickness3;   %#
    thickness4 = 60*nm;                         %#
    V_applied4(i) = V(1) + E_sio2*thickness4;   %#
    thickness5 = 100*nm;                        %#
    V_applied5(i) = V(1) + E_sio2*thickness5;   %#
    %#############################################
    
    charge(i) = -E_sio2*epsilon_ox;
end
charge(301)=0;
charge=charge*1e-4;

for i=1:length-1
    capacitor1(i) = (charge(i+1)-charge(i))/(V_applied1(i+1) - V_applied1(i));
end
capacitor1(end) = capacitor1(end-1);

for i=1:length-1
    capacitor2(i) = (charge(i+1)-charge(i))/(V_applied2(i+1) - V_applied2(i));
end
capacitor2(end) = capacitor2(end-1);

for i=1:length-1
    capacitor3(i) = (charge(i+1)-charge(i))/(V_applied3(i+1) - V_applied3(i));
end
capacitor3(end) = capacitor3(end-1);

for i=1:length-1
    capacitor4(i) = (charge(i+1)-charge(i))/(V_applied4(i+1) - V_applied4(i));
end
capacitor4(end) = capacitor4(end-1);

for i=1:length-1
    capacitor5(i) = (charge(i+1)-charge(i))/(V_applied5(i+1) - V_applied5(i));
end
capacitor5(end) = capacitor5(end-1);

%% Plotting the LFCV Curve
figure;
hold on;
plot(V_applied1,capacitor1/capacitor1(1),'LineWidth',1.5);
plot(V_applied2,capacitor2/capacitor2(1),'LineWidth',1.5);
plot(V_applied3,capacitor3/capacitor3(1),'LineWidth',1.5);
plot(V_applied4,capacitor4/capacitor4(1),'LineWidth',1.5);
plot(V_applied5,capacitor5/capacitor5(1),'LineWidth',1.5);
title('LFCV Curve');
xlabel('V_{applied} \rightarrow');
ylabel('C/C_i \rightarrow');
xlim([-2 3]);
ylim([0 1.2]);
legend('10 nm','20 nm','40 nm', '60 nm', '100 nm');
hold off;