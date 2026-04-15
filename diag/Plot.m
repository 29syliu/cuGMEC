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
4. 代码块(2)中一键替换double或者float，取决于cuGMEC模拟时的参数，然后运行。
5. 运行代码块(3)画能量增长率，(4)画频率，(5)画2D扰动，(6)画FFT。按照每一步提示设置参数。   
%}

%% (1)

addpath('../lib/BSI')

%设置参数：nx,ny,nz,tube,leftN,rightN,dt。

nx=512;
ny=64;
nz=288;
tube=2;
leftN=0;
rightN=18;
dt=0.02;

ghost=2;
numN=rightN-leftN+1;
leftX=0;
rightX=nx-1;
numX=rightX-leftX+1;

load("collocated.mat")
load("plot.mat")
load("Normalization.mat")

%% (2)

Var_all = zeros(nz,nx,10*ny);

Filename = ['/result/Phi.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
Phi = reshape(Array,nz,nx,[]);
Var_all(:,:,0*ny+1:1*ny)=Phi(:,:,1:ny);

Filename = ['/result/A.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
A = reshape(Array,nz,nx,[]);
Var_all(:,:,1*ny+1:2*ny)=A(:,:,1:ny);

Filename = ['/result/dNe.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dNe = reshape(Array,nz,nx,[]);
Var_all(:,:,2*ny+1:3*ny)=dNe(:,:,1:ny);

Filename = ['/result/dTe.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dTe = reshape(Array,nz,nx,[]);
Var_all(:,:,3*ny+1:4*ny)=dTe(:,:,1:ny);

Filename = ['/result/dPe.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dPe = reshape(Array,nz,nx,[]);
Var_all(:,:,4*ny+1:5*ny)=dPe(:,:,1:ny);

Filename = ['/result/dJpB.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dJpB = reshape(Array,nz,nx,[]);
Var_all(:,:,5*ny+1:6*ny)=dJpB(:,:,1:ny);

Filename = ['/result/w.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
w = reshape(Array,nz,nx,[]);
Var_all(:,:,6*ny+1:7*ny)=w(:,:,1:ny);

Filename = ['/result/dPi.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dPi = reshape(Array,nz,nx,[]);
Var_all(:,:,7*ny+1:8*ny)=dPi(:,:,1:ny);

Filename = ['/result/dPa.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');      
dPa = reshape(Array,nz,nx,[]);
Var_all(:,:,8*ny+1:9*ny)=dPa(:,:,1:ny);            

Filename = ['/result/dPb.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
dPb = reshape(Array,nz,nx,[]);
Var_all(:,:,9*ny+1:10*ny)=dPb(:,:,1:ny);

Filename = ['/result/energy.bin'];
fip = fopen(['./',Filename],'rb');
[Energy,~]=fread(fip,inf,'double');

Filename = ['/result/frequency.bin'];
fip = fopen(['./',Filename],'rb');
[Array,~]=fread(fip,inf,'double');
Frequency = reshape(Array,numX,[]);


%% (3)
 
% 画能量随时间变化
figure;
plot(log(Energy));


% 根据这张图手动选择计算线性增长率的起点和终点
grow0 = 3000;
grow1 = 4000;

growth = log(Energy(grow0:grow1))-log(Energy(grow0-1:grow1-1));


% 输出增长率
mean(growth)/2*(1/dt)*va0/L0
% mean(growth)/2*(1/dt)



%% (4)

% 选择在哪个径向坐标诊断频率
diagX = 115;
frequency = squeeze(Frequency(diagX,:));

% 画这一点电势随时间变化
figure;
plot(log(abs(frequency)),'b');
title(diagX);

% 根据这张图手动选择计算频率的起点和终点，以及图上相应的周期数
fre0 = 200;
fre1 = 100;
period = 2;

% 输出频率
1/((fre1-fre0)/period*2*dt*L0/va0)*2*pi



%% (5)

z_grid=((0.5:nz-0.5)/nz*2*pi-pi)/tube;
y_grid=(0.5:ny-0.5)/ny*2*pi-pi;
psigrid = psi_p(:,1); 

Nx_plot = 1024;
Ny_plot = 1024;

thgrid=(0.5:Ny_plot-0.5)/Ny_plot*2*pi-pi;
thplot = zeros(Nx_plot,Ny_plot);
for i=1:Nx_plot
    thplot(i,:)=thgrid;
end

qtheta = q.*theta_b;
psiplot = repmat(psiplot,[Ny_plot,1])';
qplot = repmat(qplot,[Ny_plot,1])';
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

    tagz = 1-qthetaplot;
    coor_3d = [tagz(:),psiplot(:),thplot(:)];
    derivative = uint64([0,0,0])';
    is_periodic = [true,false,false]';
    range_zxy = {[z_grid(1),z_grid(1)+2*pi/tube];[psigrid(1),psigrid(end)];[y_grid(1),y_grid(1)+2*pi]};
    Result = bspline(uint64(4),is_periodic,range_zxy,A_nonshift,coor_3d,derivative);
    Result = reshape(Result,Nx_plot,[]);
    Vars_interp{ival} = Result;

    surf(Rplot,Zplot,Result);

    shading interp;

    a = min(min(Result))
    b = max(max(Result))

    x = [a,2/3*a,1/3*a,0,1/3*b,2/3*b,b];
    colors = [0,1,1;0.3,0.3,0.8;0,0,1;1,1,1;1,0.2,0.2;1,0.3,0.6;1,1,0];
    color_positons = [1,2,3,4,5,6,7];
    xq = linspace(a,b,256);
    colormap_custom = interp1(x,colors(color_positons,:),xq,'linear');
    colormap(colormap_custom);

    view(0,90);
    axis equal;
    xlabel('R');
    ylabel('Z');
    title(name(ival));
    grid off;
    hold on;

    plot(Rplot(1,:),Zplot(1,:),'b','LineWidth',1.5)
    plot(Rplot(end,:),Zplot(end,:),'b','LineWidth',1.5)


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

