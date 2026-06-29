%% Normalization

mu0 = 4*pi*1e-7;
mi = 1.672621637e-27;
me = 9.10938215e-31;
e = 1.6021766208e-19;

B0 = B00;
L0 = L00;
VA0 = B0/sqrt(mu0*IonMass*mi*ni_c);
P0 = B0^2/(2*mu0);
T0 = B0^2/(2*mu0)*L0^3;
N0 = L0^-3;
gcon0 = L0^-2;
gcov0 = L0^2;
J0 = L0^3;
Jp_B_0 = 1/(mu0*L0);
Bcon0 = B0/L0;

RHO0 = rho(1,1);
RHO1 = rho(end,1);
PSITMAX = psit_max;

gridNx = N_rho;
gridNy = N_theta;
NFP = 1;

IonBeta = 0;
AlphaBeta = 0;
BeamBeta = 0;

if IonType~=3
    if IonType==1
        IonBeta = ni(1,1)*Ti(1,1)*1e19*1000*1.6021766208e-19/P0;
    elseif IonType==2
        IonBeta = Pi(1,1)/P0;
    end
end

if AlphaType~=3
    if AlphaType==1
        AlphaBeta = na(1,1)*Ta(1,1)*1e19*1000*1.6021766208e-19/P0;
    elseif AlphaType==2
        AlphaBeta = Pa(1,1)/P0;
    end
end

if BeamType~=3
    if BeamType==1
        BeamBeta = nb(1,1)*Tb(1,1)*1e19*1000*1.6021766208e-19/P0;
    elseif BeamType==2
        BeamBeta = Pb(1,1)/P0;
    end
end

OutputFileName = strcat(outputPath,'normalization2D.mat');
save(OutputFileName,'gridNx','gridNy','NFP','B0','L0','VA0','RHO0','RHO1','PSITMAX','IonBeta','AlphaBeta','BeamBeta');

