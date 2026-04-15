
transfer2xy_1d = @(f, df_dpsi, f0) {f/f0, Dpsi.*df_dpsi/f0};
transfer2xy_3d = @(f, df_dpsi, df_dth, f0) {f/f0, df_dpsi.*Dpsi/f0, df_dth/f0};
transfer2xy_3d2 = @(f, df_dpsi, df_dth, df_dpsi2, df_dpsidth, df_dth2, f0) ...
    {f/f0, ...
    df_dpsi.*Dpsi/f0, df_dth/f0, ...
    df_dpsi2.*Dpsi.^2/f0, Dpsi.*df_dpsidth/f0, df_dth2/f0};

output_ = @(f) reshape(f',[],1);

metric_seq_ = @(g) {g{1,1}, g{1,2}, g{2,2}, g{1,3}, g{2,3}, g{3,3}};
output_metric_ = @(g,dg_dx,dg_dy) reshape(mbind_col(metric_seq_,g,dg_dx,dg_dy),[],1);


%% 1d

%electron

ne = zeros(N_rho,N_theta);
dne_dr = zeros(N_rho,N_theta);
dne_dpsi = zeros(N_rho,N_theta);
for i = 1:N_rho
    temp_rho = rho(i,1);
    for j = 1:1:length(nepoly)
        ne(i,:) = ne(i,:) + nepoly(j)*temp_rho^(j-1);
    end
    for j = 2:1:length(nepoly)
        dne_dr(i,:) = dne_dr(i,:) + (j-1)*nepoly(j)*temp_rho^(j-2);
    end
end
dne_dpsi = dne_dr.*dpsi_p_dr.^(-1);
ne = 1.0e19*ne;
dne_dpsi = 1.0e19*dne_dpsi;
ne_out = transfer2xy_1d(ne,dne_dpsi,n0);

Te = zeros(N_rho,N_theta);
dTe_dr = zeros(N_rho,N_theta);
dTe_dpsi = zeros(N_rho,N_theta);
for i=1:N_rho
    temp_rho = rho(i,1);
    for j = 1:1:length(Tepoly)
        Te(i,:) = Te(i,:) + Tepoly(j)*temp_rho^(j-1);
    end
    for j = 2:1:length(Tepoly)
        dTe_dr(i,:) = dTe_dr(i,:) + (j-1)*Tepoly(j)*temp_rho^(j-2);
    end
end
dTe_dpsi = dTe_dr.*dpsi_p_dr.^(-1);

Pe = ne.*Te*1000*1.6021766208e-19;
dPe_dpsi = (dne_dpsi.*Te+dTe_dpsi.*ne)*1000*1.6021766208e-19;
Pe_out = transfer2xy_1d(Pe,dPe_dpsi,P0);

Te = Te*1000*1.6021766208e-19;
dTe_dpsi = dTe_dpsi*1000*1.6021766208e-19;
Te_out = transfer2xy_1d(Te,dTe_dpsi,T0);



%Ion

ni = zeros(N_rho,N_theta);
dni_dr = zeros(N_rho,N_theta);
dni_dpsi = zeros(N_rho,N_theta);
for i = 1:N_rho
    temp_rho = rho(i,1);
    for j = 1:1:length(nipoly)
        ni(i,:) = ni(i,:) + nipoly(j)*temp_rho^(j-1);
    end
    for j = 2:1:length(nipoly)
        dni_dr(i,:) = dni_dr(i,:) + (j-1)*nipoly(j)*temp_rho^(j-2);
    end
end
dni_dpsi = dni_dr.*dpsi_p_dr.^(-1);
ni = 1.0e19*ni;
dni_dpsi = 1.0e19*dni_dpsi;
ni_out = transfer2xy_1d(ni,dni_dpsi,n0);

Ti = zeros(N_rho,N_theta);
dTi_dr = zeros(N_rho,N_theta);
dTi_dpsi = zeros(N_rho,N_theta);

Pi = zeros(N_rho,N_theta);
dPi_dr = zeros(N_rho,N_theta);
dPi_dpsi = zeros(N_rho,N_theta);

if IonType==1

    for i=1:N_rho
        temp_rho = rho(i,1);
        for j = 1:1:length(Tipoly)
            Ti(i,:) = Ti(i,:) + Tipoly(j)*temp_rho^(j-1);
        end
        for j = 2:1:length(Tipoly)
            dTi_dr(i,:) = dTi_dr(i,:) + (j-1)*Tipoly(j)*temp_rho^(j-2);
        end
    end
    dTi_dpsi = dTi_dr.*dpsi_p_dr.^(-1);
    Ti_out = transfer2xy_1d(Ti,dTi_dpsi,1);

    Pi = ni.*Ti*1000*1.6021766208e-19;
    dPi_dpsi = (dni_dpsi.*Ti+dTi_dpsi.*ni)*1000*1.6021766208e-19;
    Pi_out = transfer2xy_1d(Pi,dPi_dpsi,P0);

elseif IonType==2

    for i=1:N_rho
        temp_rho = rho(i,1);
        for j = 1:1:length(Pipoly)
            Pi(i,:) = Pi(i,:) + Pipoly(j)*temp_rho^(j-1);
        end
        for j = 2:1:length(Pipoly)
            dPi_dr(i,:) = dPi_dr(i,:) + (j-1)*Pipoly(j)*temp_rho^(j-2);
        end
    end
    dPi_dpsi = dPi_dr.*dpsi_p_dr.^(-1);
    Pi_out = transfer2xy_1d(Pi,dPi_dpsi,P0);

    Ti_out = transfer2xy_1d(Ti,dTi_dpsi,1);

end



%Alpha

na = zeros(N_rho,N_theta);
dna_dr = zeros(N_rho,N_theta);
dna_dpsi = zeros(N_rho,N_theta);

Ta = zeros(N_rho,N_theta);
dTa_dr = zeros(N_rho,N_theta);
dTa_dpsi = zeros(N_rho,N_theta);

Pa = zeros(N_rho,N_theta);
dPa_dr = zeros(N_rho,N_theta);
dPa_dpsi = zeros(N_rho,N_theta);

if AlphaType~=3

    na = zeros(N_rho,N_theta);
    dna_dr = zeros(N_rho,N_theta);
    dna_dpsi = zeros(N_rho,N_theta);
    for i = 1:N_rho
        temp_rho = rho(i,1);
        for j = 1:1:length(napoly)
            na(i,:) = na(i,:) + napoly(j)*temp_rho^(j-1);
        end
        for j = 2:1:length(napoly)
            dna_dr(i,:) = dna_dr(i,:) + (j-1)*napoly(j)*temp_rho^(j-2);
        end
    end
    dna_dpsi = dna_dr.*dpsi_p_dr.^(-1);
    na = 1.0e19*na;
    dna_dpsi = 1.0e19*dna_dpsi;
    na_out = transfer2xy_1d(na,dna_dpsi,n0);

    if AlphaType==1

        Ta = zeros(N_rho,N_theta);
        dTa_dr = zeros(N_rho,N_theta);
        dTa_dpsi = zeros(N_rho,N_theta);
        for i=1:N_rho
            temp_rho = rho(i,1);
            for j = 1:1:length(Tapoly)
                Ta(i,:) = Ta(i,:) + Tapoly(j)*temp_rho^(j-1);
            end
            for j = 2:1:length(Tapoly)
                dTa_dr(i,:) = dTa_dr(i,:) + (j-1)*Tapoly(j)*temp_rho^(j-2);
            end
        end
        dTa_dpsi = dTa_dr.*dpsi_p_dr.^(-1);
        Ta_out = transfer2xy_1d(Ta,dTa_dpsi,1);

        Pa = na.*Ta*1000*1.6021766208e-19;
        dPa_dpsi = (dna_dpsi.*Ta+dTa_dpsi.*na)*1000*1.6021766208e-19;
        Pa_out = transfer2xy_1d(Pa,dPa_dpsi,P0);

    elseif AlphaType==2

        for i=1:N_rho
            temp_rho = rho(i,1);
            for j = 1:1:length(Papoly)
                Pa(i,:) = Pa(i,:) + Papoly(j)*temp_rho^(j-1);
            end
            for j = 2:1:length(Papoly)
                dPa_dr(i,:) = dPa_dr(i,:) + (j-1)*Papoly(j)*temp_rho^(j-2);
            end
        end
        dPa_dpsi = dPa_dr.*dpsi_p_dr.^(-1);
        Pa_out = transfer2xy_1d(Pa,dPa_dpsi,P0);

        Ta_out = transfer2xy_1d(Ta,dTa_dpsi,1);

    end

else

    na_out = transfer2xy_1d(na,dna_dpsi,n0);
    Ta_out = transfer2xy_1d(Ta,dTa_dpsi,1);
    Pa_out = transfer2xy_1d(Pa,dPa_dpsi,P0);

end


%Beam

nb = zeros(N_rho,N_theta);
dnb_dr = zeros(N_rho,N_theta);
dnb_dpsi = zeros(N_rho,N_theta);

Tb = zeros(N_rho,N_theta);
dTb_dr = zeros(N_rho,N_theta);
dTb_dpsi = zeros(N_rho,N_theta);

Pb = zeros(N_rho,N_theta);
dPb_dr = zeros(N_rho,N_theta);
dPb_dpsi = zeros(N_rho,N_theta);

if BeamType~=3

    nb = zeros(N_rho,N_theta);
    dnb_dr = zeros(N_rho,N_theta);
    dnb_dpsi = zeros(N_rho,N_theta);
    for i = 1:N_rho
        temp_rho = rho(i,1);
        for j = 1:1:length(nbpoly)
            nb(i,:) = nb(i,:) + nbpoly(j)*temp_rho^(j-1);
        end
        for j = 2:1:length(nbpoly)
            dnb_dr(i,:) = dnb_dr(i,:) + (j-1)*nbpoly(j)*temp_rho^(j-2);
        end
    end
    dnb_dpsi = dnb_dr.*dpsi_p_dr.^(-1);
    nb = 1.0e19*nb;
    dnb_dpsi = 1.0e19*dnb_dpsi;
    nb_out = transfer2xy_1d(nb,dnb_dpsi,n0);

    if BeamType==1

        Tb = zeros(N_rho,N_theta);
        dTb_dr = zeros(N_rho,N_theta);
        dTb_dpsi = zeros(N_rho,N_theta);
        for i=1:N_rho
            temp_rho = rho(i,1);
            for j = 1:1:length(Tbpoly)
                Tb(i,:) = Tb(i,:) + Tbpoly(j)*temp_rho^(j-1);
            end
            for j = 2:1:length(Tbpoly)
                dTb_dr(i,:) = dTb_dr(i,:) + (j-1)*Tbpoly(j)*temp_rho^(j-2);
            end
        end
        dTb_dpsi = dTb_dr.*dpsi_p_dr.^(-1);
        Tb_out = transfer2xy_1d(Tb,dTb_dpsi,1);

        Pb = nb.*Tb*1000*1.6021766208e-19;
        dPb_dpsi = (dnb_dpsi.*Tb+dTb_dpsi.*nb)*1000*1.6021766208e-19;
        Pb_out = transfer2xy_1d(Pb,dPb_dpsi,P0);

    elseif BeamType==2
        
        for i=1:N_rho
            temp_rho = rho(i,1);
            for j = 1:1:length(Pbpoly)
                Pb(i,:) = Pb(i,:) + Pbpoly(j)*temp_rho^(j-1);
            end
            for j = 2:1:length(Pbpoly)
                dPb_dr(i,:) = dPb_dr(i,:) + (j-1)*Pbpoly(j)*temp_rho^(j-2);
            end
        end
        dPb_dpsi = dPb_dr.*dpsi_p_dr.^(-1);
        Pb_out = transfer2xy_1d(Pb,dPb_dpsi,P0);

        Tb_out = transfer2xy_1d(Tb,dTb_dpsi,1);

    end

else

    nb_out = transfer2xy_1d(nb,dnb_dpsi,n0);
    Tb_out = transfer2xy_1d(Tb,dTb_dpsi,1);
    Pb_out = transfer2xy_1d(Pb,dPb_dpsi,P0);

end



% output

ni_out = mbind_cell_col(output_,ni_out)';
Ti_out = mbind_cell_col(output_,Ti_out)';
Pi_out = mbind_cell_col(output_,Pi_out)';

ne_out = mbind_cell_col(output_,ne_out)';
Te_out = mbind_cell_col(output_,Te_out)';
Pe_out = mbind_cell_col(output_,Pe_out)';

na_out = mbind_cell_col(output_,na_out)';
Ta_out = mbind_cell_col(output_,Ta_out)';

nb_out = mbind_cell_col(output_,nb_out)';
Tb_out = mbind_cell_col(output_,Tb_out)';

%% Output Jacobi and Bny

Jxyz = JB*Dpsi;
dJxyz_dpsi = dJB_dpsi*Dpsi;
dJxyz_dth = dJB_dth*Dpsi;
dJxyz_dxi = dJB_dxi*Dpsi;
Jxyz_out = transfer2xy_3d(Jxyz, dJxyz_dpsi, dJxyz_dth, J0);
Jxyz_out = mbind_cell_col(output_,Jxyz_out)';

%%

Bny = Dpsi./Jxyz;
dBny_dpsi = -Dpsi./Jxyz.^2.*dJxyz_dpsi;
dBny_dth = -Dpsi./Jxyz.^2.*dJxyz_dth;
dBny_dxi = -Dpsi./Jxyz.^2.*dJxyz_dxi;
Bny_out = transfer2xy_3d(Bny, dBny_dpsi, dBny_dth, Bcon0);
Bny_out = mbind_cell_col(output_,Bny_out)';

%% Output Shift Metric、Align Metric

gcon_s_out = mbind_cell_col(output_,output_metric_(gcon_s,dgcon_s_dx,dgcon_s_dy))'./gcon0;
gcov_s_out = mbind_cell_col(output_,output_metric_(gcov_s,dgcov_s_dx,dgcov_s_dy))'./gcov0;

SFAcon_s_out = mbind_cell_col(output_,output_metric_(SFAcon_s,dSFAcon_s_dx,dSFAcon_s_dy))'./gcon0;
SFAcov_s_out = mbind_cell_col(output_,output_metric_(SFAcov_s,dSFAcov_s_dx,dSFAcov_s_dy))'./gcov0;

%%

mu0 = 4*pi*1e-7;
mi = 1.672621637e-27;
e = 1.6021766208e-19;

Rho = sqrt(IonMass*mi.*Pi./ni)/e./B;

dRho_dpsi = sqrt(IonMass*mi)/e*(0.5*(Pi./ni).^(-0.5).*(dPi_dpsi./ni-dni_dpsi.*Pi./ni.^2).*B-dB_dpsi.*sqrt(Pi./ni))./B.^2;     

dRho_dth = -sqrt(IonMass*mi.*Pi./ni)/e./B.^2.*dB_dth;


Rho_out = transfer2xy_3d(Rho, dRho_dpsi, dRho_dth, L0);
Rho_out = mbind_cell_col(output_,Rho_out)';

%%

Va = B./sqrt(mu0*IonMass*mi*ni);

dVa_dpsi = dB_dpsi./sqrt(mu0*IonMass*mi*ni) - 0.5*B.*(mu0*IonMass*mi*ni).^(-1.5)*mu0*IonMass*mi.*dni_dpsi;

dVa_dth = dB_dth./sqrt(mu0*IonMass*mi*ni);


Va_out = transfer2xy_3d(Va, dVa_dpsi, dVa_dth, va0);
Va_out = mbind_cell_col(output_,Va_out)';

%% Output RZ

transfer2xy_RZ = @(f1, f2, f0) {f1/f0, f2/f0};

RZ_out = transfer2xy_RZ(R, Z, L0);
RZ_out = mbind_cell_col(output_,RZ_out)';

%% Output B

Bs_out = transfer2xy_3d2(B, dB_dpsi, dB_dth, dB_dpsi2, dB_dpsidth, dB_dth2, B0);
Bs_out = mbind_cell_col(output_,Bs_out)';

%% Output Jp_B

jp_B_out = transfer2xy_3d(jp_B, djp_B_dpsi, djp_B_dth, Jp_B_0);
jp_B_out = mbind_cell_col(output_,jp_B_out)';

%% Output PIC

q_out = output_(q)';

i_curr_out = transfer2xy_1d(i_curr,di_curr_dpsi,curr0);
i_curr_out = mbind_cell_col(output_,i_curr_out)';

g_curr_out = transfer2xy_1d(g_curr,dg_curr_dpsi,curr0);
g_curr_out = mbind_cell_col(output_,g_curr_out)';

K = -(i_curr.*gcon{1,2} + g_curr.*gcon{1,3})./gcon{1,1}./q;
dK_dth = (-(i_curr.*dgcon_dth{1,2} + g_curr.*dgcon_dth{1,3}).*gcon{1,1} + dgcon_dth{1,1}.*(i_curr.*gcon{1,2} + g_curr.*gcon{1,3}))./q./gcon{1,1}.^2;
transferK = @(f, df_dth, f0) {f/f0, df_dth/f0};
K_out = transferK(K,dK_dth,delta0);
K_out = mbind_cell_col(output_,K_out)';

PIC_data = [q_out,i_curr_out,g_curr_out,K_out];


%% Output

qtheta = output_(qtheta)';

output1d_data = [ni_out,Ti_out,Pi_out,ne_out,Te_out,Pe_out,na_out,Ta_out,nb_out,Tb_out];

output3d_data = [gcon_s_out,gcov_s_out,Jxyz_out,Bny_out,jp_B_out,Rho_out,Va_out,RZ_out,Bs_out];

length(output3d_data)/N_rho/N_theta

Output_data = [qtheta,SFAcon_s_out,SFAcov_s_out,PIC_data,output1d_data,output3d_data];

%% Out Equilibrium

if MHDstaggered == 0
    fip = fopen([outputPath, equilibriumName],'wb');
    fwrite(fip,Output_data,'double');
    fclose(fip);
end