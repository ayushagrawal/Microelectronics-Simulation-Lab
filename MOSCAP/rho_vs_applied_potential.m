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
V_applied = zeros(length,1);
charge = zeros(length,1);
capacitor = zeros(length,1);

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
    thickness = 10*nm;
    V_applied(i) = V(1) + E_sio2*thickness;
    charge(i) = -E_sio2*epsilon_ox;
end
charge(301)=0;
charge=charge*1e-4;
for i=1:length-1
    capacitor(i) = (charge(i+1)-charge(i))/(V_applied(i+1) - V_applied(i));
end
capacitor(end) = capacitor(end-1);

%% Plotting Charge Density
figure;
semilogy(Vs,abs(charge),'LineWidth',1.5);
line([0 0], [1e-10 1e-5],'Color','red','LineWidth',1);
phi_bp = (k*T/q)*log(doping_na/n_i);
[temp, index] = min(abs(Vs - 2*phi_bp));
text(-0.2,1e-9,'Accumulation','HorizontalAlignment','center');
text(Vs(index)/2, 1e-9, 'Depletion','HorizontalAlignment','center');
text(0.9, 1e-9, 'Inversion ','HorizontalAlignment','center');
line([Vs(index) Vs(index)], [1e-10 1e-5],'Color','red','LineWidth',1);
title('Charge Density Vs Surface Potential');
xlabel('V_s \rightarrow');
ylabel('\rho (C-cm^{-2}) \rightarrow');

figure;
semilogy(V_applied,abs(charge),'LineWidth',1.5);
line([0 0], [1e-10 1e-6],'Color','red','LineWidth',1);
phi_bp = (k*T/q)*log(doping_na/n_i);
[temp, index] = min(abs(Vs - 2*phi_bp));
text(-1,1e-8,'Accumulation','HorizontalAlignment','center');
text(V_applied(index)/2, 1e-8, 'Depletion','HorizontalAlignment','center');
text(1.5 + V_applied(index)/2, 1e-8, 'Inversion ','HorizontalAlignment','center');
line([V_applied(index) V_applied(index)], [1e-10 1e-6],'Color','red','LineWidth',1);
title('Charge Density Vs Applied Potential');
xlabel('V_{applied} \rightarrow');
ylabel('\rho (C-cm^{-2}) \rightarrow');
xlim([-2 3]);

%% Plotting the LFCV Curve
figure;
plot(V_applied,capacitor/capacitor(1),'LineWidth',1.5);
line([0 0], [0 1.2],'Color','red','LineWidth',1);
phi_bp = (k*T/q)*log(doping_na/n_i);
[temp, index] = min(abs(Vs - 2*phi_bp));
text(-1,0.6,'Accumulation','HorizontalAlignment','center');
text(V_applied(index)/2, 0.6, 'Depletion','HorizontalAlignment','center');
text(1.5 + V_applied(index)/2, 0.6, 'Inversion ','HorizontalAlignment','center');
line([V_applied(index) V_applied(index)], [0 1.2],'Color','red','LineWidth',1.2);
title('LFCV Curve');
xlabel('V_{applied} \rightarrow');
ylabel('C/C_i \rightarrow');
xlim([-2 3]);
ylim([0 1.2]);