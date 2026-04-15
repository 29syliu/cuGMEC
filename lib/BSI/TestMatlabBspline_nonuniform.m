clear all
%%
mex COMPFLAGS='$COMPFLAGS /std:c++20' ...
    '-ID:\Backup\Document\C++\Bspline zq\BSplineInterpolation-dev\src\include\BSplineInterpolation' ...
   Test_nonuniform.cpp

%% uniform
N_point = 100;
order = uint64(4);
isperiodic = false;
range = [0,1];
derivative = uint64(0);
x_index = (0:N_point)'/N_point;
range_cell = {range};
f_val = x_index.^2;
coor1d = (0:50)'/50;
f_ref = coor1d.^2;

%  result1d = Test_nonuniform(order,isperiodic,range_cell,f_val(1:end),coor1d,derivative);
result1d = bspline(order,isperiodic,range_cell,f_val(1:end),coor1d,derivative);

figure
hold on
plot(x_index,f_val)
plot(coor1d,result1d)

figure
plot(coor1d,result1d-f_ref)


%% non-uniform
% f1d_ = @(x) x.^2;
% df1d_ = @(x) 2*x;

f1d_ = @(x) exp(-x.^2/0.3^2);
df1d_ = @(x) (-2*x/0.3^2).*exp(-x.^2/0.3^2);

order = uint64(4);
isperiodic = false;
derivative = uint64(0); derivative1 = uint64(1);

x_index = sort(rand(100,1),1);
range_cell = {x_index};
f_val = f1d_(x_index);

coor1d = (0:50)'/50;
f_ref = f1d_(coor1d);
df_ref = df1d_(coor1d);

% Test_nonuniform()

%  result1d = Test_nonuniform(order,isperiodic,range_cell,f_val(1:end),coor1d,derivative);
result1d = bspline(order,isperiodic,range_cell,f_val(1:end),coor1d,derivative);
result1d_dx = bspline(order,isperiodic,range_cell,f_val(1:end),coor1d,derivative1);
% result1d = bspline(order,isperiodic,range_cell,f_val(1:end));

figure
subplot 121
hold on
plot(x_index,f_val)
plot(coor1d,result1d)
plot(coor1d,f_ref)
subplot 122
plot(coor1d,result1d-f_ref)

figure
subplot 121
hold on
plot(coor1d,result1d_dx)
plot(coor1d,df_ref)
subplot 122
plot(coor1d,result1d_dx-df_ref)



%%
mex COMPFLAGS='$COMPFLAGS /std:c++20' ...
    '-IC:\Users\ALFVEN\Desktop\BSI\BSplineInterpolation' ...
   bspline.cpp

%% 2D non-uniform + periodic

isperiodic = [false, true]';
f_ = @(x,y) x.^2.* sin(y);
df_dx_ = @(x,y) 2*x.* sin(y);
nx_grid = 96; ny_grid = 48;
nx_interp = 64; ny_interp = 128;
x_grid = sort(rand(100,1),1)*2 - 1;
% x_grid = (0:nx_grid)'/nx_grid*2 - 1;
y_grid = (0.5:ny_grid-0.5)'/ny_grid*2*pi-pi;
x_interp = (0:nx_interp)'/nx_interp*2 - 1;
y_interp = (0.5:ny_interp-0.5)/ny_interp*2*pi-pi;

[XX_grid,YY_grid]=ndgrid(x_grid,y_grid);
[XX_interp,YY_interp]=ndgrid(x_interp,y_interp);
coor_interp = [XX_interp(:),YY_interp(:)];

order = uint64(4);
range_2d = {x_grid; [y_grid(1),y_grid(1)+2*pi]};
% range_2d = [-1,1;y_grid(1),y_grid(1)+2*pi];
FF = f_(XX_grid,YY_grid);
FF_dx = df_dx_(XX_grid,YY_grid);
derivative = uint64([0,0])';
derivative_dx= uint64([1,0])';
FF_interp_exact = f_(XX_interp,YY_interp);
FF_interp_dx_exact = df_dx_(XX_interp,YY_interp);

FF_interp = reshape(bspline(order,isperiodic,range_2d,FF,coor_interp,derivative),nx_interp+1,[]);
FF_interp_dx = reshape(bspline(order,isperiodic,range_2d,FF,coor_interp,derivative_dx),nx_interp+1,[]);

figure
subplot 131
mesh(XX_grid,YY_grid,FF)
subplot 132
mesh(XX_interp,YY_interp,FF_interp)
subplot 133
mesh(XX_interp,YY_interp,FF_interp -FF_interp_exact)

figure
subplot 131
mesh(XX_grid,YY_grid,FF_dx)
subplot 132
mesh(XX_interp,YY_interp,FF_interp_dx)
subplot 133
mesh(XX_interp,YY_interp,FF_interp_dx -FF_interp_dx_exact)

%% Test

[RRho, Tth] = ndgrid(rho,theta);

nr = length(rho);

figure
mesh(RRho,Tth,nu)

order = uint64(5);
isperiodic = [false, true]';
derivative= uint64([0,0])';
derivative_dx= uint64([0,1])';
range_2d = {[rho(1),rho(end)]; [theta(1),theta(1)+2*pi]};

coor_interp = [RRho(:),Tth(:)];

nu_interp = reshape(bspline(order,isperiodic,range_2d,nu,coor_interp,derivative),nr,[]);
nu_interp_dt = reshape(bspline(order,isperiodic,range_2d,nu,coor_interp,derivative_dx),nr,[]);


figure
subplot 131
mesh(RRho,Tth,nu_interp_dt)
subplot 132
mesh(RRho,Tth,nu_t)
subplot 133
mesh(RRho,Tth,(nu_t-nu_interp_dt)/mean(mean(abs(nu_t))))


