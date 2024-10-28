function out = get_T(T)
% T must be given in energy units of keV

%% Defining physical constants
e    = 1.602176634e-19;
u    = 1.66053906660e-27;
mH   = 1.007276466620409*u;
mD   = 2.014*u;
mT   = 3.01604928*u;
mHe3 = 3.016029*u;
mHe4 = 4.002602*u;
mB11 = 11.009306*u;

%% D + T -> He4 + n reaction
s = load('D_T_-_a_n.txt');
E = s(:,1)*1e3; % Energy in keV
s = s(:,2)*1e-28; % Cross section in m^2
mr = mD*mT/(mD+mT); % reduced mass
out.DT = sqrt(1e3*e)*4/sqrt(2*pi*mr)/T^1.5*trapz(E,s.*E.*exp(-E/T));
%out.DT = sqrt(1e3*e)*sqrt(8/pi)*(mr/T)^1.5/mD^2*trapz(E,s.*E.*exp(-mr*E/T/mD));

%% D + D -> T + p reaction
s = load('D_D_-_T_p.txt');
E = s(:,1)*1e3; % Energy in keV
s = s(:,2)*1e-28; % Cross section in m^2
mr = mD*mD/(mD+mD); % reduced mass
out.DD_p = sqrt(1e3*e)*4/sqrt(2*pi*mr)/T^1.5*trapz(E,s.*E.*exp(-E/T));

%% D + D -> He-3 + n reaction
s = load('D_D_-_3He_n.txt');
E = s(:,1)*1e3; % Energy in keV
s = s(:,2)*1e-28; % Cross section in m^2
mr = mD*mD/(mD+mD); % reduced mass
out.DD_n = sqrt(1e3*e)*4/sqrt(2*pi*mr)/T^1.5*trapz(E,s.*E.*exp(-E/T));
out.DD = out.DD_n + out.DD_p;

%% D + He-3 -> He-4 + p reaction
s = load('D_3He_-_4He_p-endf.txt');
E = s(:,1)*1e3; % Energy in keV
s = s(:,2)*1e-28; % Cross section in m^2
mr = mD*mHe3/(mD+mHe3); % reduced mass
out.DHe3 = sqrt(1e3*e)*4/sqrt(2*pi*mr)/T^1.5*trapz(E,s.*E.*exp(-E/T));

%% p + B-11 -> 3 alpha reaction
s = load('p_11B_-_3a.txt');
E = s(:,1)*1e3; % Energy in keV
s = s(:,2)*1e-28; % Cross section in m^2
mr = mH*mB11/(mH+mB11); % reduced mass
out.HB = sqrt(1e3*e)*4/sqrt(2*pi*mr)/T^1.5*trapz(E,s.*E.*exp(-E/T));

%% T + T -> He-4 + 2n reaction
s = load('T_T_-_4He_n_n.txt');
E = s(:,1)*1e3; % Energy in keV
s = s(:,2)*1e-28; % Cross section in m^2
mr = mT*mT/(mT+mT); % reduced mass
out.TT = sqrt(1e3*e)*4/sqrt(2*pi*mr)/T^1.5*trapz(E,s.*E.*exp(-E/T));

%% T + He3 -> n + p + He-4 reaction
s = load('T_3He_-_n_p_4He.txt');
E = s(:,1)*1e3; % Energy in keV
s = s(:,2)*1e-28; % Cross section in m^2
mr = mT*mHe3/(mT+mHe3); % reduced mass
out.THe3a = sqrt(1e3*e)*4/sqrt(2*pi*mr)/T^1.5*trapz(E,s.*E.*exp(-E/T));

%% T + He3 -> D + He-4 reaction
s = load('T_3He_-_D_4He.txt');
E = s(:,1)*1e3; % Energy in keV
s = s(:,2)*1e-28; % Cross section in m^2
mr = mD*mHe4/(mD+mHe4); % reduced mass
out.THe3b = sqrt(1e3*e)*4/sqrt(2*pi*mr)/T^1.5*trapz(E,s.*E.*exp(-E/T));
out.THe3 = (59*out.THe3a + 41*out.THe3b)/100;

%% He-3 + He3 -> 2p + He-4 reaction
s = load('3He_3He_-_p_p_4He.txt');
E = s(:,1)*1e3; % Energy in keV
s = s(:,2)*1e-28; % Cross section in m^2
mr = mHe3*mHe3/(mHe3+mHe3); % reduced mass
out.He3He3 = sqrt(1e3*e)*4/sqrt(2*pi*mr)/T^1.5*trapz(E,s.*E.*exp(-E/T));

end