function dydt = nuclear(t,y)

%% Characteristic loss time
tau_n = 1;
tau_p = 1e10;
tau_E = 1;
t0 = 40;
dt = 40;
dnhdt_aux = 0;
dnddt_aux = 1e17;
dntdt_aux = 1e17;
dnhe3dt_aux = 0;
dnbdt_aux = 0;
Saux = 3*0.65e6;

%% Getting particle densities
n_n     = y(1);
n_H     = y(2);
n_D     = y(3);
n_T     = y(4);
n_He3   = y(5);
n_alpha = y(6);
n_B     = y(7);
T       = y(8);

%% Getting cross sections at temperature T
out      = get_T(T);
s_DT     = out.DT;
s_DD     = out.DD;  
s_DHe3   = out.DHe3;  
s_HB     = out.HB; 
s_TT     = out.TT;  
s_THe3   = out.THe3;  
s_He3He3 = out.He3He3;

%% Reaction rates
dnndt     = n_D*n_T*s_DT + n_D^2*s_DD/2 + 59*n_T*n_He3*s_THe3/100 + 2*n_T^2*s_TT - n_n/tau_n;
dnhdt     = n_D^2*s_DD/2 + n_D*n_He3*s_DHe3 + 2*n_He3^2*s_He3He3 + 59*n_T*n_He3*s_THe3/100 - n_H*n_B*s_HB - n_H/tau_p + dnhdt_aux*(tanh((t - t0)/dt) + 1)/2;
dnddt     = 41*n_T*n_He3*s_THe3/100 - n_D*n_T*s_DT - 2*n_D^2*s_DD - n_D*n_He3*s_DHe3 - n_D/tau_p + dnddt_aux*(tanh((t - t0)/dt) + 1)/2;
dntdt     = n_D^2*s_DD/2 - n_D*n_T*s_DT - n_T*n_He3*s_THe3 - 2*n_T^2*s_TT - n_T/tau_p + dntdt_aux*(tanh((t - t0)/dt) + 1)/2;
dnhe3dt   = n_D^2*s_DD/2 - n_D*n_He3*s_DHe3 - 2*n_He3^2*s_He3He3 - n_T*n_He3*s_THe3 - n_He3/tau_p + dnhe3dt_aux*(tanh((t - t0)/dt) + 1)/2;
dnalphadt = n_D*n_T*s_DT + n_D*n_He3*s_DHe3 + n_He3^2*s_He3He3 + n_T*n_He3*s_THe3 + n_T^2*s_TT + 3*n_H*n_B*s_HB - n_alpha/tau_p;
dnbdt     = - n_H*n_B*s_HB - n_B/tau_p + dnbdt_aux*(tanh((t - t0)/dt) + 1)/2;
dnidt     = dnhdt + dnddt + dntdt + dnhe3dt + dnalphadt + dnbdt;
dnedt     = dnhdt + dnddt + dntdt + 2*dnhe3dt + 2*dnalphadt + 5*dnbdt;
dn0dt     = dnidt + dnedt;

%% Particle densities
ni = n_H + n_D + n_T + n_He3 + n_alpha + n_B;
ne = n_H + n_D + n_T + 2*n_He3 + 2*n_alpha + 5*n_B;
n0 = ni + ne;

%% Energy equation
e    = 1.602176634e-19;
Se   = (3.5e3*n_D*n_T*s_DT + 4.85e3*n_D^2*s_DD + 18.3e3*n_D*n_He3*s_DHe3 + 41/100*14.32e3*n_T*n_He3*s_THe3)*e*1e3;
dTdt = 2*(Se + Saux*(tanh((t - 30)/10) + 1)/2)/n0/3/e/1e3 - T*dn0dt/n0 - 2*T/tau_E/3;

%% Output
dydt = [dnndt;
        dnhdt;
        dnddt;
        dntdt;
        dnhe3dt;
        dnalphadt;
        dnbdt;
        dTdt];

end