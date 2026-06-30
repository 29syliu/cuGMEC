%%
%先得到所有平衡量在PEST坐标下的导数，再转换为Shifted Metric坐标下的导数，并进行归一化

transfer2xy_1d = @(f, df_drho, f0) {f/f0, df_drho.*Drho/f0};
transfer2xy_3d = @(f, df_drho, df_dtheta, df_dphi, f0) {f/f0, df_drho.*Drho/f0, (df_dtheta+df_dphi.*q)/f0, df_dphi/f0};
transfer2xy_3d2 = @(f, df_drho, df_dtheta, df_dphi, df_drho2, df_drhodtheta, df_drhodphi, df_dtheta2, df_dthetadphi, df_dphi2, f0) ...
    {f/f0, ...
    df_drho.*Drho/f0, (df_dtheta+df_dphi.*q)/f0, df_dphi/f0, ...
    df_drho2.*Drho.^2/f0, (df_drhodtheta.*Drho+df_dphi.*dq_drho.*Drho+df_drhodphi.*Drho.*q)/f0, df_drhodphi.*Drho/f0, ...
    (df_dtheta2+2.*df_dthetadphi.*q+df_dphi2.*q.^2)/f0, (df_dthetadphi+df_dphi2.*q)/f0, df_dphi2/f0};

output_ = @(f) reshape(permute(f,[3,2,1]),[],1);

metric_seq_ = @(g) {g{1,1}, g{1,2}, g{2,2}, g{1,3}, g{2,3}, g{3,3}};
output_metric_ = @(g,dg_dx,dg_dy,dg_dz) reshape(mbind_col(metric_seq_,g,dg_dx,dg_dy,dg_dz),[],1);

zero3d = zeros(size(rho));

%% 1d

%electron

ne = 1.0e19*ne;
dne_dr = 1.0e19*dne_dr;
ne_out = transfer2xy_1d(ne,dne_dr,N0);

Te = Te*1000*1.6021766208e-19;
dTe_dr = dTe_dr*1000*1.6021766208e-19;
Te_out = transfer2xy_1d(Te,dTe_dr,T0);

Pe = ne.*Te;
dPe_dr = dne_dr.*Te+dTe_dr.*ne;
Pe_out = transfer2xy_1d(Pe,dPe_dr,P0);



%Ion

if IonType~=3

    ni = 1.0e19*ni;
    dni_dr = 1.0e19*dni_dr;
    ni_out = transfer2xy_1d(ni,dni_dr,N0);

    if IonType==1

        Ti_out = transfer2xy_1d(Ti,dTi_dr,1);
  
        Pi = ni.*Ti*1000*1.6021766208e-19;
        dPi_dr = (dni_dr.*Ti+dTi_dr.*ni)*1000*1.6021766208e-19;
        Pi_out = transfer2xy_1d(Pi,dPi_dr,P0);

    elseif IonType==2

        Ti_out = transfer2xy_1d(Ti,dTi_dr,VA0);
        
        Pi_out = transfer2xy_1d(Pi,dPi_dr,P0);

    end

else

    ni_out = transfer2xy_1d(zero3d,zero3d,N0);
    Ti_out = transfer2xy_1d(zero3d,zero3d,1);
    Pi_out = transfer2xy_1d(zero3d,zero3d,P0);

end



%Alpha

if AlphaType~=3

    na = 1.0e19*na;
    dna_dr = 1.0e19*dna_dr;
    na_out = transfer2xy_1d(na,dna_dr,N0);

    if AlphaType==1

        Ta_out = transfer2xy_1d(Ta,dTa_dr,1);

        Pa = na.*Ta*1000*1.6021766208e-19;
        dPa_dr = (dna_dr.*Ta+dTa_dr.*na)*1000*1.6021766208e-19;
        Pa_out = transfer2xy_1d(Pa,dPa_dr,P0);

    elseif AlphaType==2

        Ta_out = transfer2xy_1d(Ta,dTa_dr,VA0);

        Pa_out = transfer2xy_1d(Pa,dPa_dr,P0);

    end

else

    na_out = transfer2xy_1d(zero3d,zero3d,N0);
    Ta_out = transfer2xy_1d(zero3d,zero3d,1);
    Pa_out = transfer2xy_1d(zero3d,zero3d,P0);

end



%Beam

if BeamType~=3

    nb = 1.0e19*nb;
    dnb_dr = 1.0e19*dnb_dr;
    nb_out = transfer2xy_1d(nb,dnb_dr,N0);

    if BeamType==1

        Tb_out = transfer2xy_1d(Tb,dTb_dr,1);

        Pb = nb.*Tb*1000*1.6021766208e-19;
        dPb_dr = (dnb_dr.*Tb+dTb_dr.*nb)*1000*1.6021766208e-19;
        Pb_out = transfer2xy_1d(Pb,dPb_dr,P0);

    elseif BeamType==2

        Tb_out = transfer2xy_1d(Tb,dTb_dr,VA0);

        Pb_out = transfer2xy_1d(Pb,dPb_dr,P0);

    end

else

    nb_out = transfer2xy_1d(zero3d,zero3d,N0);
    Tb_out = transfer2xy_1d(zero3d,zero3d,1);
    Pb_out = transfer2xy_1d(zero3d,zero3d,P0);

end



% output

ne_out = mbind_cell_col(output_,ne_out)';
Te_out = mbind_cell_col(output_,Te_out)';
Pe_out = mbind_cell_col(output_,Pe_out)';

ni_out = mbind_cell_col(output_,ni_out)';
Ti_out = mbind_cell_col(output_,Ti_out)';
Pi_out = mbind_cell_col(output_,Pi_out)';

na_out = mbind_cell_col(output_,na_out)';
Ta_out = mbind_cell_col(output_,Ta_out)';

nb_out = mbind_cell_col(output_,nb_out)';
Tb_out = mbind_cell_col(output_,Tb_out)';

%% Jacobian

Jxyz = JPEST*Drho;
dJxyz_drho = dJPEST_drho*Drho;
dJxyz_dtheta = dJPEST_dtheta*Drho;
dJxyz_dphi = dJPEST_dphi*Drho;

Jxyz_out = transfer2xy_3d(Jxyz, dJxyz_drho, dJxyz_dtheta, dJxyz_dphi, J0);
Jxyz_out = mbind_cell_col(output_,Jxyz_out)';

%% Bny

Bny = 2*psit_max*Drho*rho./q./Jxyz;
dBny_drho = 2*psit_max*Drho*(1./q./Jxyz-rho.*dq_drho./Jxyz./q.^2-rho.*dJxyz_drho./q./Jxyz.^2);
dBny_dtheta = 2*psit_max*Drho*(-rho.*dJxyz_dtheta./q./Jxyz.^2);
dBny_dphi = 2*psit_max*Drho*(-rho.*dJxyz_dphi./q./Jxyz.^2);

Bny_out = transfer2xy_3d(Bny, dBny_drho, dBny_dtheta, dBny_dphi, B0/L0);
Bny_out = mbind_cell_col(output_,Bny_out)';

%% Metric

gconSFT_out = mbind_cell_col(output_,output_metric_(gconSFT,dgconSFT_dx,dgconSFT_dy,dgconSFT_dz))'./gcon0;
gcovSFT_out = mbind_cell_col(output_,output_metric_(gcovSFT,dgcovSFT_dx,dgcovSFT_dy,dgcovSFT_dz))'./gcov0;

gconSFA_out = mbind_cell_col(output_,output_metric_(gconSFA,dgconSFA_dx,dgconSFA_dy,dgconSFA_dz))'./gcon0;
gcovSFA_out = mbind_cell_col(output_,output_metric_(gcovSFA,dgcovSFA_dx,dgcovSFA_dy,dgcovSFA_dz))'./gcov0;

%% FLR 

dPi_drho = dPi_dr;
dni_drho = dni_dr;

Rho = sqrt(IonMass*mi.*Pi./ni)/e./B;
dRho_drho = sqrt(IonMass*mi)/e*(0.5*(Pi./ni).^(-0.5).*(dPi_drho./ni-dni_drho.*Pi./ni.^2).*B-dB_drho.*sqrt(Pi./ni))./B.^2;     
dRho_dtheta = -sqrt(IonMass*mi.*Pi./ni)/e./B.^2.*dB_dtheta;
dRho_dphi = -sqrt(IonMass*mi.*Pi./ni)/e./B.^2.*dB_dphi;

Rho_out = transfer2xy_3d(Rho, dRho_drho, dRho_dtheta, dRho_dphi, L0);
Rho_out = mbind_cell_col(output_,Rho_out)';

%% Va

Va = B./sqrt(mu0*IonMass*mi*ni);
dVa_drho = dB_drho./sqrt(mu0*IonMass*mi*ni) - 0.5*B.*(mu0*IonMass*mi*ni).^(-1.5)*mu0*IonMass*mi.*dni_drho;
dVa_dtheta = dB_dtheta./sqrt(mu0*IonMass*mi*ni);
dVa_dphi = dB_dphi./sqrt(mu0*IonMass*mi*ni);

Va_out = transfer2xy_3d(Va, dVa_drho, dVa_dtheta, dVa_dphi, VA0);
Va_out = mbind_cell_col(output_,Va_out)';

%% RZ

transfer2xy_RZ = @(f1, f2, f0) {f1/f0, f2/f0};

RZ_out = transfer2xy_RZ(R, Z, L0);
RZ_out = mbind_cell_col(output_,RZ_out)';

%% B

Bs_out = transfer2xy_3d2(B, dB_drho, dB_dtheta, dB_dphi, dB_drho2, dB_drhodtheta, dB_drhodphi, dB_dtheta2, dB_dthetadphi, dB_dphi2, B0);
Bs_out = mbind_cell_col(output_,Bs_out)';

%% JpB

jp_B_out = transfer2xy_3d(jp_B, djp_B_drho, djp_B_dtheta, djp_B_dphi, Jp_B_0);
jp_B_out = mbind_cell_col(output_,jp_B_out)';

%% q

q_out = transfer2xy_1d(q, dq_drho, 1);
q_out = mbind_cell_col(output_,q_out)';

%% psip

% dpsip_drho = dpsip_dr;

% psip_out = transfer2xy_1d(psip, dpsip_drho, B0*L0*L0);
% psip_out = mbind_cell_col(output_,psip_out)';

psip_out = q_out;

%% Output

qtheta = output_(qtheta)';

output1d_data = [q_out,psip_out,ni_out,Ti_out,Pi_out,ne_out,Te_out,Pe_out,na_out,Ta_out,nb_out,Tb_out];
output2d_data = [gconSFT_out,gcovSFT_out,Jxyz_out,Bny_out,jp_B_out,Rho_out,Va_out,RZ_out,Bs_out];

Output_data = [qtheta,gconSFA_out,gcovSFA_out,output1d_data,output2d_data];

fip = fopen([outputPath, equilibriumName],'wb');
fwrite(fip,Output_data,'double');
fclose(fip);

%% Phase-space mapping variables

MP = mi;
QE = e;
B = B / B0;
J = Jxyz / J0;
psip = q / (B0*L0^2);
SFAcovyz = gcovSFA{2,3} / L0^2;
SFAcovzz = gcovSFA{3,3} / L0^2;
SFAcovyz = SFAcovyz(:,1+ghost:size(B,2)+ghost,:);
SFAcovzz = SFAcovzz(:,1+ghost:size(B,2)+ghost,:);
theta = theta_pest;

normalizationFile = fullfile(outputPath, 'normalization3D.mat');
save(normalizationFile, ...
    'MP', 'QE', 'B', 'J', 'psip', 'SFAcovyz', 'SFAcovzz', ...
    'q', 'R', 'Z', 'rho', 'theta', 'phi', 'NFP', 'fullTorus', '-append');

%% function

function out = mbind_col(f,varargin)
    % f:: x -> [a]
    % mbind_col:: [varargin] >>= f
    out = f(varargin{1});
    for i=2:length(varargin)
        out = [out;f(varargin{i})];
    end
end

function out = mbind_cell_col(f,C)
    % f:: x -> [a]
    % mbind_col:: [varargin] >>= f
    out = f(C{1});
    for i=2:length(C)
        out = [out;f(C{i})];
    end
end
