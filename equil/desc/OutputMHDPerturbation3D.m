%% Disturb

if leftN==rightN

    leftN = leftN*tube;

    disturb_N = leftN;
    disturb_index = radialIndex;

    disturb_q = q(disturb_index,1,1);
    disturb_M = round(disturb_N*round(disturb_q*10)/10)
    disturb_x = (rho(:,1) - rho(1,1))/Drho;
    disturb_y = (0.5:N_theta-0.5)/N_theta*2*pi-pi;
    disturb_z = ((0.5:N_phi-0.5)/N_phi*2*pi-pi)/tube;

    [Z3D,X3D,Y3D] = ndgrid(disturb_z,disturb_x,disturb_y);

    sigma = width;
    center_x = disturb_x(disturb_index);
    Gaussianx_ = @(x)exp(-(x-center_x).^2/(2*sigma^2));

    dw0 = amplitude;
    dw = dw0*Gaussianx_(X3D).*cos(disturb_M*Y3D - disturb_N*Z3D);
    % dw = dw + dw0*Gaussianx_(X3D).*cos(disturb_M*Y3D - disturb_N*Z3D);

    outputd_ = @(f) reshape(f,[],1);
    Output_disturb = outputd_(dw);

elseif leftN~=rightN

    leftN = leftN*tube;
    rightN = rightN*tube;

    disturb_N = rightN;
    disturb_index = radialIndex;

    disturb_q = q(disturb_index,1,1);
    disturb_M = round(disturb_N*round(disturb_q*10)/10)
    disturb_x = (rho(:,1) - rho(1,1))/Drho;
    disturb_y = (0.5:N_theta-0.5)/N_theta*2*pi-pi;
    disturb_z = ((0.5:N_phi-0.5)/N_phi*2*pi-pi)/tube;

    [Z3D,X3D,Y3D] = ndgrid(disturb_z,disturb_x,disturb_y);

    sigma = width;
    center_x = disturb_x(disturb_index);
    Gaussianx_ = @(x)exp(-(x-center_x).^2/(2*sigma^2));

    dw0 = amplitude;
    % dw = dw0*Gaussianx_(X3D).*cos(disturb_M*Y3D - disturb_N*Z3D);
    dw = dw + dw0*Gaussianx_(X3D).*cos(disturb_M*Y3D - disturb_N*Z3D);

    for disturb_n = leftN:tube:(rightN-tube)

        disturb_m = round(disturb_n*round(disturb_q*10)/10)
        disturb_y = (0.5:N_theta-0.5)/N_theta*2*pi-pi+2*pi/disturb_M*disturb_m;
        disturb_z = ((0.5:N_phi-0.5)/N_phi*2*pi-pi)/tube+2*pi/disturb_N*disturb_n;

        [Z3D,X3D,Y3D] = ndgrid(disturb_z,disturb_x,disturb_y);
        dw = dw + dw0*Gaussianx_(X3D).*cos(disturb_m*Y3D - disturb_n*Z3D);

    end

    outputd_ = @(f) reshape(f,[],1);
    Output_disturb = outputd_(dw);

end

fip = fopen([outputPath, perturbationName],'wb');
fwrite(fip,Output_disturb,'double');
fclose(fip);

figure;
mesh(squeeze(dw(N_phi/2,:,:)));
title('perturbation \zeta=0');

figure;
mesh(squeeze(dw(:,:,N_theta/2)));
title('perturbation \theta=0');

figure;
mesh(squeeze(dw(:,disturb_index,:)));
title('perturbation x=rho(radialIndex)');

