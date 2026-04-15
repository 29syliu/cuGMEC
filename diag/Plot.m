%%

%{

进行任何操作前，请仔细阅读这一段注释。
进行任何操作前，请仔细阅读这一段注释。
进行任何操作前，请仔细阅读这一段注释。

请按照以下步骤进行操作:

1. 将输出的result文件夹复制进plot.m所在的文件夹。复制文件夹，不是复制文件夹里的文件。
2. 生成cuGMEC时用的MATLAB脚本指定了当时的inputPath和outputPath，在这两个路径下找到
   collocated.mat，plot.mat，Normalization.mat，并复制进plot.m所在的文件夹。
3. 代码块(1)设置参数并运行。
4. 代码块(2)中一键替换double或者double，取决于cuGMEC模拟时的参数，然后运行。
5. 运行代码块(3)画能量增长率，(4)画频率，(5)画2D扰动，(6)画FFT。按照每一步提示设置参数。   
%}

%% (1)

addpath('../lib/BSI')

%设置参数：nx,ny,nz,tube,leftN,rightN,dt。

nx=256;
ny=32;
nz=96;
tube=6;
leftN=1;
rightN=6;
dt=0.02;

ghost=2;
numN=rightN-leftN+1;
leftX=0;
rightX=nx-1;
numX=rightX-leftX+1;

load("standard2D.mat")
load("plot2D.mat")
load("normalization2D.mat")

%% (2)

Var_all = zeros(nz,nx,10*ny);

Filename = ['Phi.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
Phi = reshape(Array,nz,nx,[]);
Var_all(:,:,0*ny+1:1*ny)=Phi(:,:,1:ny);

Filename = ['A.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
A = reshape(Array,nz,nx,[]);
Var_all(:,:,1*ny+1:2*ny)=A(:,:,1:ny);

Filename = ['dNe.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dNe = reshape(Array,nz,nx,[]);
Var_all(:,:,2*ny+1:3*ny)=dNe(:,:,1:ny);

Filename = ['dTe.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dTe = reshape(Array,nz,nx,[]);
Var_all(:,:,3*ny+1:4*ny)=dTe(:,:,1:ny);

Filename = ['dPe.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dPe = reshape(Array,nz,nx,[]);
Var_all(:,:,4*ny+1:5*ny)=dPe(:,:,1:ny);

Filename = ['dJpB.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dJpB = reshape(Array,nz,nx,[]);
Var_all(:,:,5*ny+1:6*ny)=dJpB(:,:,1:ny);

Filename = ['w.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
w = reshape(Array,nz,nx,[]);
Var_all(:,:,6*ny+1:7*ny)=w(:,:,1:ny);

Filename = ['dPi.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dPi = reshape(Array,nz,nx,[]);
Var_all(:,:,7*ny+1:8*ny)=dPi(:,:,1:ny);

Filename = ['dPa.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');      
dPa = reshape(Array,nz,nx,[]);
Var_all(:,:,8*ny+1:9*ny)=dPa(:,:,1:ny);            

Filename = ['dPb.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dPb = reshape(Array,nz,nx,[]);
Var_all(:,:,9*ny+1:10*ny)=dPb(:,:,1:ny);

Filename = ['frequency.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
Frequency = reshape(Array,numX,[]);

Filename = ['IonDensity.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
IonDensity = reshape(Array,numX,[]);

Filename = ['AlphaDensity.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
AlphaDensity = reshape(Array,numX,[]);

Filename = ['BeamDensity.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
BeamDensity = reshape(Array,numX,[]);

Filename = ['AlphaDiffusivity.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
AlphaDiffusivity = reshape(Array,numX,[]);

Filename = ['BeamDiffusivity.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
BeamDiffusivity = reshape(Array,numX,[]);

Filename = ['IonDiffusivity.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
IonDiffusivity = reshape(Array,numX,[]);

%% (3)
va0 = VA0;

Filename = ['diagnose.bin'];
fip = fopen(['./',Filename],'rb');
[Energy,~]=fread(fip,inf,'double');

% 画能量随时间变化
figure;
plot(log(Energy));


% 根据这张图手动选择计算线性增长率的起点和终点
grow0 = 45000;
grow1 = 50000;

growth = log(Energy(grow0:grow1))-log(Energy(grow0-1:grow1-1));


% 输出增长率
mean(growth)/2*(1/dt)
% mean(growth)/2*(1/dt)

figure;
plot(growth/2*(1/dt))



%% (4)

Filename = ['frequency.bin'];
fip = fopen(['./',Filename],'rb');
[frequency,~]=fread(fip,inf,'double');
frequency = reshape(frequency,numX,[]);

% 选择在哪个径向坐标诊断频率
diagX = 275;
frequency = squeeze(frequency(diagX,:));

% 画这一点电势随时间变化
figure;
plot(log(abs(frequency)),'b');
title(diagX);

% 根据这张图手动选择计算频率的起点和终点，以及图上相应的周期数
fre0 = 3650;
fre1 = 5988;
period = 4;

% 输出频率
1/((fre1-fre0)/period*2*dt*L0/VA0)

%%

Filename = ['Epara.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
Epara = reshape(Array,numX,[]);
Epara = Epara';

Filename = ['EparaES.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
EparaES = reshape(Array,numX,[]);
EparaES = EparaES';

diagT = 5000;

figure;
plot(rho(:,1),Epara(diagT,:),'r','LineWidth',2);
hold on;
plot(rho(:,1),EparaES(diagT,:),'b','LineWidth',2);
plot(rho(:,1),EparaES(diagT,:)-Epara(diagT,:),'g','LineWidth',2);
% D:\ITER\ITER\0.1-0.9(512 64)\SingleN\N=22 dt0.02 8000\1
%%

figure;
plot(squeeze(abs(Phi(1,:,32))));
hold on;
plot(squeeze(abs(A(1,:,32))));

%% (5)

load("standard2D.mat")
load("plot2D.mat")
load("Normalization2D.mat")

zori = -pi/tube;
zrange = 2*pi/tube;

z_grid=((0.5:nz-0.5)/nz*2*pi-pi)/tube;
y_grid=(0.5:ny-0.5)/ny*2*pi-pi;
rho_grid = rho(:,1); 

[Nx_plot, Ny_plot] = size(qplot);

thgrid=(0.5:Ny_plot-0.5)/Ny_plot*2*pi-pi;
thplot = zeros(Nx_plot,Ny_plot);
for i=1:Nx_plot
    thplot(i,:)=thgrid;
end

qtheta = q.*theta_pest;
% rhoplot = repmat(rhoplot,[Ny_plot,1])';
% qplot = repmat(qplot,[Ny_plot,1])';
qthetaplot = qplot.*thplot;


% 画2D扰动图，ival=1:1就是画Phi，ival=1:2就是画Phi和A，以此类推。

name = ["Phi","A","dNe","dTe","dPe","dJpB","w","dPi","dPa","dPb"];

figure; 

for ival=1:1

    A_shift = Var_all(:,:,ny*(ival-1)+1:ny*ival);
    A_shift = A_shift/max(max(max(A_shift)));

    A_nonshift = zeros(nz,nx,ny);
    for ix = 1:nx
        for iy = 1:ny
            z_grid_new = z_grid+qtheta(ix,iy);
            A_nonshift(:,ix,iy) = bspline(uint64(4),true,{[z_grid(1),z_grid(1)+2*pi/tube]},A_shift(:,ix,iy),z_grid_new',uint64(0));
        end
    end


    tagz = 0.4-qthetaplot;
    coor_3d = [tagz(:),rhoplot(:),thplot(:)];
    derivative = uint64([0,0,0])';
    is_periodic = [true,false,false]';
    range_zxy = {[z_grid(1),z_grid(1)+2*pi/tube];[rho_grid(1),rho_grid(end)];[y_grid(1),y_grid(1)+2*pi]};
    Result = bspline(uint64(4),is_periodic,range_zxy,A_nonshift,coor_3d,derivative);
    Result = reshape(Result,Nx_plot,[]);

    Vars_interp{ival} = Result;

    surf(Rplot,Zplot,Result);

    shading interp;

    a = min(min(Result));
    b = max(max(Result));
    % 
    % x = [a,2/3*a,1/3*a,0,1/3*b,2/3*b,b];
    % colors = [0,1,1;0.3,0.3,0.8;0,0,1;1,1,1;1,0.2,0.2;1,0.3,0.6;1,1,0];
    % color_positons = [1,2,3,4,5,6,7];
    % xq = linspace(a,b,256);
    % colormap_custom = interp1(x,colors(color_positons,:),xq,'linear');
    % colormap(colormap_custom);

    % x = [a,1/2*a,0,1/2*b,b];
    % colors = [0.3,0.3,0.8;0,0,1;1,1,1;1,0.2,0.2;1,0.3,0.6];
    % color_positons = [1,2,3,4,5];
    % xq = linspace(a,b,256);
    % colormap_custom = interp1(x,colors(color_positons,:),xq,'linear');
    % colormap(colormap_custom);

    x = [a,1/2*a,0,1/2*b,b];
    colors = [0,1,1;0,0,1;1,1,1;1,0.2,0.2;1,1,0];
    color_positons = [1,2,3,4,5];
    xq = linspace(a,b,256);
    colormap_custom = interp1(x,colors(color_positons,:),xq,'linear');
    colormap(colormap_custom);

    % x = [a,0,b];
    % colors = [0,0,1;1,1,1;1,0,0];
    % color_positons = [1,2,3];
    % xq = linspace(a,b,256);
    % colormap_custom = interp1(x,colors(color_positons,:),xq,'linear');
    % colormap(colormap_custom);

    view(0,90);
    axis equal;
    xlabel('R');
    ylabel('Z');
    title(name(ival));
    grid off;
    hold on;

    % plot(Rplot(1,:),Zplot(1,:),'b','LineWidth',1.5)
    % plot(Rplot(end,:),Zplot(end,:),'b','LineWidth',1.5)

    plot3(Rplot(1,:),Zplot(1,:),ones(Nx_plot,Ny_plot)*10,'b','LineWidth',1.5)
    plot3(Rplot(end,:),Zplot(end,:),ones(Nx_plot,Ny_plot)*10,'b','LineWidth',1.5)
    plot3(Rplot(933,:),Zplot(933,:),ones(Nx_plot,Ny_plot)*10,'b','LineWidth',1.5)
    % plot3(Rplot(314,:),Zplot(314,:),ones(Nx_plot,Ny_plot)*10,'k-','LineWidth',1.5)
    % plot3(Rplot(253,:),Zplot(253,:),temp,'k-.','LineWidth',1.5)
    % plot3(Rplot(514,:),Zplot(514,:),temp,'k:','LineWidth',1.5)


    xlim([min(min(Rplot))-0.5,max(max(Rplot))+0.5]);
    ylim([min(min(Zplot))-0.5,max(max(Zplot))+0.5]);

    axisSize = 12;
    labelSize = 14;
    titleSize = 14;
    legendSize = 10;
    barSize = 12;
    lineWidth = 1.5;
    markerSize = 12;
    
    set(gca,'FontName','Times New Roman','FontSize',axisSize,'LineWidth',lineWidth);
    hold on;

    clim([a b]);
    set(gca,'FontName','Times New Roman','FontSize',axisSize,'LineWidth',1.2);
    set(gca,'xgrid','on','ygrid','on','GridLineWidth',1.0,'GridAlpha',0.3, ... 
    'MinorGridAlpha',0.3,'GridColor','k','MinorGridColor','k','GridLineStyle','-','MinorGridLineStyle','-');

    xlabel('$R/\mathrm{m}$','interpreter','latex','FontSize',labelSize);
    ylabel('$Z/\mathrm{m}$','interpreter','latex','FontSize',labelSize);
    box on;

end

%% (6)

load("collocated.mat")
load("plot.mat")

issave = true;
frontsize = 15;
legendsize = 11;
plotLineWidth = 1;
axisLineWidth = 1;
Linewidth = 1;

rhoplot = rhoplot';
is_rho = true;
qplot = qplot';
n = tube;

%%

ival = 1;
Result2 = Vars_interp{ival};

Plot_threshold = 1e-8;
minFFT = 0;

xrang = [rho(1,1) rho(N_rho,1)];

fftstartnum = 1;
fftendnum = 20;
fftbegin = 1;

clear h h2 lgd lgd2

% Result_fft = zeros(Nx_plot,Nth_plot);
% 
% for ipsi = 1:Nx_plot
%     temp = fft(Result2(ipsi,:));
%     Result_fft(ipsi,:) = temp;
% end

Result_fft = fft(Result2,[],2);
% Result_fft-Result_fft2

if is_rho
    plot_x = rhoplot;
    x_name = '$\sqrt{s}$';
else
    plot_x = rhoplot;
    x_name = '$\psi_p$';
end

R = abs(Result_fft(fftbegin:end,fftbegin:end));% max_fft = max(max(R(:,fftstartnum:fftendnum)));
R_max = max(R);
R_max_max = max(R_max);
[~,plot_index] = find(R_max>R_max_max*Plot_threshold);
plot_index=plot_index(plot_index<fftendnum+1);
plot_index=plot_index(plot_index>fftstartnum-1);
index_mid = plot_index(round(length(plot_index)/2)+1)
N_left = index_mid - plot_index(1);

% line
dfftnum = floor((plot_index(end)-plot_index(1))/20) + 1
fft_m = plot_index(1):plot_index(end);
ratio_q = fft_m/n;
ratio_x = interp1(qplot,plot_x,ratio_q,'spline');

p0 = figure;
set(p0,'Position',[1000,500,800,600]);
hold on

% axes('position',[0.72,0.15,0.25,0.77])
% axes('position',[0.366296296296296,0.082+0.565-0.05,0.25,0.36])
% hold on

set(gca,'LooseInset',get(gca,'TightInset'))
title('FFT','interpreter','none','FontSize',1);
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1.5);
set(gca,'xgrid','on','ygrid','on','GridLineWidth',1.0,'GridAlpha',0.3,'MinorGridAlpha',0.3,'GridColor','k','MinorGridColor','k','GridLineStyle','-','MinorGridLineStyle','-');
maxR = max(max(R(:,plot_index(1)+1:plot_index(end)+1)))

index = 1;
for i = plot_index(1):plot_index(end)
    % if mod(i,2) == 0
    %     h(index)=plot(plot_x(fftbegin:end),R(:,i+1)./maxR,'-','linewidth',Linewidth,'Color','#FF0000');  
    % else
    %     h(index)=plot(plot_x(fftbegin:end),R(:,i+1)./maxR,'-','linewidth',Linewidth,'Color','#000000'); 
    % end
    
    if(max(R(:,i+1)./maxR)>minFFT)
        h(index)=plot(plot_x(fftbegin:end),R(:,i+1)./maxR,'-','linewidth',2); 
        lgd(index)={['m=',num2str(i)]};
        index = index + 1;
    end
end
hlgd=legend(h,lgd,'location','northwest');% legend('boxoff');
set(hlgd, 'Fontsize', legendsize);
xlim(xrang)
ylim([0,1])
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1.5);
set(gca, 'LineWidth',axisLineWidth)
xlabel([x_name],'Interpreter','latex','fontsize',20);
ylabel('$\delta\phi$','Interpreter','latex','fontsize',20);
box on;

% for i=1:N_left-1
%     h=plot([ratio_x(i),ratio_x(i)],[0,1],'--');
%     set(h,'handlevisibility','off');
%     if mod(i,dfftnum)==0
%         text(ratio_x(i),1.05*1,num2str(fft_m(i)),HorizontalAlignment="center")
%     end
% end

% axesNew = axes('position',get(gca,'position'),'visible','off','Xlim',get(gca,'Xlim'),'Ylim',get(gca,'Ylim'));
% xlim(xrang)
% hold on
% index = 1;
% for i=index_mid:plot_index(end)
%     % if mod(i,2) == 0
%     %     h2(index)=plot(plot_x(fftbegin:end),R(:,i+1)./maxR,'-','linewidth',Linewidth,'Color','#FF0000');  
%     % else
%     %     h2(index)=plot(plot_x(fftbegin:end),R(:,i+1)./maxR,'-','linewidth',Linewidth,'Color','#000000'); 
%     % end
%     h2(index)=plot(plot_x(fftbegin:end),R(:,i+1)./maxR,'-','linewidth',Linewidth); 
%     lgd2(index)={['m=',num2str(i)]};
%     index = index + 1;
% end
% 
% 
% 
% % xlabel('$\psi$','Interpreter','latex','fontsize',frontsize);
% % % title(['\omega_{AE}=',num2str(omega_huge(ind_plot))]);
% % ylabel('|U_m|','fontsize',20);
% % grid on;a
% hlgd2=legend(h2,lgd2,'location','northeast');% legend('boxoff');
% 
% % Line
% for i=1:length(ratio_x)
%     h=plot([ratio_x(i),ratio_x(i)],[0,1],'--');
%     set(h,'handlevisibility','off');
%     if mod(i,1)==0
%         text(ratio_x(i),1.05*1,num2str(fft_m(i)),HorizontalAlignment="center")
%     end
% end
% set(hlgd2, 'Fontsize', legendsize);


% if issave
%   %  save_normal(p0,[OutputName,'_PhiFFT'],400)
% %     plot_latex(p0,'GMEC_phi_fft',400)
%     saveas(p0,'GMEC_phi_fft.fig');
%     saveas(p0,'GMEC_phi_fft.png');
% end

xlim([0,1])
xticks(0:0.2:1)
yticks(0:0.25:1)

%%



axisSize = 12;
labelSize = 14;
titleSize = 14;
legendSize = 10;
barSize = 12;
lineWidth = 1.5;
markerSize = 12;


Filename = ['amplitude.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
Amplitude = reshape(Array,numN,numX,[]);
% Amplitude = A500;

% Amplitude = A40000;
% Amplitude = AT;

% Amplitude = A4;
% 
% red: #FF0000
% green: #00800
% pink: #ff13a6
% blue: #0000FF
% blue: #02a8a8
% blue: #00BFBF

% colors = [
%     '#FF0000';  % 亮红色
%     '#0000FF';  % 深蓝色
%     '#00FF00';  % 亮绿色
%     '#FFA500';  % 橙色
%     '#800080';  % 紫色
%     '#00FFFF';  % 青色
%     '#FFFF00';  % 黄色
%     '#FFC0CB';  % 粉色
%     '#008000'   % 深绿色
% ];

% myHex = ['#FF0000';'#008000';'#ff13a6';'#0000FF';'#02a8a8';];
% myColor = zeros(size(myHex,1),3);
% for i = 1:size(myHex,1)
%     myColor(i,:) = sscanf(myHex(i,2:end),'%2x%2x%2x',[1 3])/255;
% end


figure;
% figure; set(gcf,'position',[350 50 400 260]);
% subplot('position',[0.15 0.19 300/400 300*0.5/0.8/260])
% hold on; 
% set(gca,'ColorOrder',myColor,'NextPlot','replacechildren');
hold on; 
set(gca,'FontName','Times New Roman','FontSize',axisSize,'LineWidth',1);
set(gca,'xgrid','on','ygrid','on','GridLineWidth',1.0,'GridAlpha',0.3, ... 
    'MinorGridAlpha',0.3,'GridColor','k','MinorGridColor','k','GridLineStyle','-','MinorGridLineStyle','-');

diagX = 80
AmplitudeResult = [];
for diagN = 1:1:rightN-leftN+1
    AmplitudeResult = [AmplitudeResult;log10(B0*L0*VA0/(35*1000)*squeeze(Amplitude(diagN,diagX,:)))'];
end

 % log10(B0*L0*VA0/(35*1000)*squeeze(Amplitude(diagN,diagX,:)))'
 % log10(B0*L0*VA0/(29.5*1000)*squeeze(Amplitude(diagN,diagX,:)))'
 % log(squeeze(Amplitude(diagN,diagX,:)))'

t = 1:size(AmplitudeResult,2)-1;
t=t;
t = t*dt;
t = t*L0/VA0*1000;
% t = t+15000*dt*L0/VA0*1000;

% [0,2,3,6,29,30]+1

for n = 1:1:rightN-leftN+1
    plot(t(1:10:end),AmplitudeResult(n,1:10:end-1),'LineWidth',1.5,'DisplayName',['n=',num2str((leftN+n-1)*tube)]);
end

h1 = legend('NumColumns',2);
h1.ItemTokenSize = [10,5];
box on;
xlim([0,t(end)])
% xlim([0,0.3])
% ylim([-10,-1]) 

[maxValue, linearIndex] = max(Amplitude(:));
[NN,XX,TT] = ind2sub(size(Amplitude),linearIndex)

title('$e\delta\phi/T_e$','interpreter','latex','FontSize',titleSize);
xlabel('$t/\mathrm{ms}$','interpreter','latex','FontSize',labelSize);
% 
ax = gca;
% ax.YTick = [-6,-5,-4,-3,-2,-1,0];
% ax.YTickLabel = {'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'};

grow0 = 3000;
grow1 = 4000;
% % 
growth = AmplitudeResult(1,grow0:grow1)-AmplitudeResult(1,grow0-1:grow1-1);
% % 
% 
% mean(growth)*(1/dt)

% 
% mean(growth)*(1/dt) 
mean(growth)*(1/dt)*VA0/L0
mean(growth)*(1/dt)
% mean(growth)*(1/dt)*VA0/L0*R835/cs
% % % 
figure;
plot((growth)*(1/dt)*VA0/L0);


%%

figure;plot(squeeze(A_nonshift(1,175:180,:))')
