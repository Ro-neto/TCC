function [t,y] = concentracoes

%% Input
tspan = 500; % Time window in s
n_H = 0;
n_D = 1e20;
n_T = 1e20;
n_He3 = 0;
n_B11 = 0;
T = 3;

%% Calculate time evolution
disp(' --- Calculating plasma burning evolution ---')
opts    = odeset('RelTol',1e-2,'AbsTol',1e-2);
[t,y]   = ode45(@nuclear,[0 tspan],[0; n_H; n_D; n_T; n_He3; 0; n_B11; T],opts);
n_H     = y(:,2);
n_D     = y(:,3);
n_T     = y(:,4);
n_He3   = y(:,5);
n_alpha = y(:,6);
n_B11   = y(:,7);
T       = y(:,8);

%% Calculating Q factor
e    = 1.602176634e-19;
s_DT = zeros(length(T),1);
s_DD = s_DT; s_DHe3 = s_DT; s_HB = s_DT; s_TT = s_DT; s_THe3 = s_DT; s_He3He3 = s_DT;
for ii = 1:length(T)
    out = get_T(T(ii));
    s_DT(ii)     = out.DT;
    s_DD(ii)     = out.DD;  
    s_DHe3(ii)   = out.DHe3;  
    s_HB(ii)     = out.HB; 
    s_TT(ii)     = out.TT;  
    s_THe3(ii)   = out.THe3;  
    s_He3He3(ii) = out.He3He3;
end
%Se   = (3.5e3*n_D.*n_T.*s_DT + 4.85e3*n_D.^2.*s_DD + 18.3e3*n_D.*n_He3.*s_DHe3 + 41/100*14.32e3*n_T.*n_He3.*s_THe3)*e*1e3;
Se   = 3.5e3*n_D.*n_T.*s_DT*e*1e3;
Saux = 3*0.65e6*(tanh((t - 30)/10) + 1)/2;
Q = Se./Saux;

%% Plotting
figure(1)
clf
hs(1) = subplot(3,1,1);
hold on
plot(t,n_H/1e20,'color',[0 0.5 0],'linewidth',3)
plot(t,n_D/1e20,'r','linewidth',3)
plot(t,n_T/1e20,'b','linewidth',3)
plot(t,n_He3/1e20,'m','linewidth',3)
plot(t,n_alpha/1e20,'k','linewidth',3)
plot(t,n_B11/1e20,'y','linewidth',3)
hold off
title('Evolução temporal das concentração')
ylabel('n ( 1\times10^{20} m^{-3} )')
legend('H','D','T','^3He','\alpha','^{11}B')
hs(2) = subplot(3,1,2);
plot(t,T,'r','linewidth',3)
title('Evolução temporal da temperatura do plasma')
ylabel('T ( keV )')
hs(3) = subplot(3,1,3);
plot(t,Se/1e6,'b','linewidth',3)
hold on
plot(t,Saux/1e6,'r','linewidth',3)
plot(t,Q,'k','linewidth',3)
hold off
title('Evolução do Ganho de Fusão Q')
legend('S_e ( MW/m^3 )','S_{aux} ( MW/m^3 )','Q')
xlabel('Tempo ( s )')
linkaxes(hs,'x')
drawnow
disp('Done')

end
