%% Normalization
mu0 = 4*pi*1e-7;
mi = 1.672621637e-27;
me = 9.10938215e-31;
e = 1.6021766208e-19;
% kB = 1.380649e-23;

B0 = Bc;
R0 = Rc;
P0 = B0^2/(2*mu0);
VA0 = B0/sqrt(mu0*IonMass*mi*ni_c);
v0 = VA0;
L0 = R0;
T0 = B0^2/(2*mu0)*L0^3;

n0 = L0^-3;
ni_c_norm = ni_c/n0;
m0 = mi*ni_c/n0;
mi_norm = mi/m0;
q0 = L0^2*sqrt(mi*ni_c/mu0);
e_norm = e/q0;
1/(e_norm*ni_c_norm)/2;

eta0 = mu0*L0*v0;
L0_R0 = L0/R0;
t0 = L0/v0;
gcon0 = L0^-2;
gcov0 = L0^2;
J0 = L0^3;
va_2_0 = B0^-2;
Jp_B_0 = 1/(mu0*L0);
Bcon0 = B0/L0;
bkcon0 = L0^-2;

delta0 = L0^-1;
curr0 = B0*L0;

phi0 = B0*v0*L0;
w0 = phi0/v0^2/L0^2;

Omega_i = e*B0/mi;
rw = VA0/(2*L0*Omega_i);

% dissipate0 = VA0*L0; % careful
dissipate0 = L0; % careful

psi_pn_max = psi_pn(end);

Pc = more_pres(1);
beta_c = Pc(1)/P0;
Tc = Pc/(ni_c*e); %eV

CHI0 = psi_p(1,1);
CHI1 = psi_p(end,1);
chimax = psi_pm;

va0 = VA0;

IonBeta = 0;
AlphaBeta = 0;
BeamBeta = 0;

if IonType~=3
    if IonType==1
        IonBeta = polyval(fliplr(nipoly), rho(1,1))*polyval(fliplr(Tipoly), rho(1,1))*1.0e19*1000*1.6021766208e-19;
    elseif IonType==2
        IonBeta = polyval(fliplr(Pipoly), rho(1,1));
    end
end

if AlphaType~=3
    if AlphaType==1
        AlphaBeta = polyval(fliplr(napoly), rho(1,1))*polyval(fliplr(Tapoly), rho(1,1))*1.0e19*1000*1.6021766208e-19;
    elseif AlphaType==2
        AlphaBeta = polyval(fliplr(Papoly), rho(1,1));
    end
end

if BeamType~=3
    if BeamType==1
        BeamBeta = polyval(fliplr(nbpoly), rho(1,1))*polyval(fliplr(Tbpoly), rho(1,1))*1.0e19*1000*1.6021766208e-19;
    elseif BeamType==2
        BeamBeta = polyval(fliplr(Pbpoly), rho(1,1));
    end
end

IonBeta = IonBeta/P0;
AlphaBeta = AlphaBeta/P0;
BeamBeta = BeamBeta/P0;

OutputFileName = strcat(outputPath,'Normalization.mat');
save(OutputFileName,'B0','L0','VA0','CHI0','CHI1','IonBeta','AlphaBeta','BeamBeta');

