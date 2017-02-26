%%
P = Polyhedron('lb',[-1 -1],'ub',[1 1]);
%P = Polyhedron('lb',-1,'ub',1);
xmax = 3;
xmin = -3;
ymax = 3;
ymin = -3;
dx = 0.05;
dy = 0.05;

dim = size(P.A,2);

grid_x = xmin:dx:xmax;
grid_y = ymin:dy:ymax;
tic;

% load data if it exits
if(exist('2d_data.mat')>0)
'loading data'
%	load('2d_data.mat');
end


%% get signed dist
if(exist('dist_array_xy') &&  sum(size(dist_array_xy) == size(meshgrid(grid_x,grid_y))) == dim)
    
'already have data'
else
keyboard;
    'computing dist'
dist_array_xy = zeros(numel(grid_x),numel(grid_y));
[msg_x,msg_y] = meshgrid(grid_x,grid_y);
if(exist('Staliro')>0)
   [msg_x,msg_y] = meshgrid(grid_x,grid_y);
   a = arrayfun(@(x,y) SignedDist([x;y],P.A,P.b), msg_x(:), msg_y(:));
   dist_array_xy_vec = reshape(a,numel(grid_x),numel(grid_y)); 
   dist_array_xy = dist_array_xy_vec;  
end
if(dim==2)
    figure;
   mesh(grid_x,grid_y,dist_array_xy_vec); 
end

end
%% get coefficients
'getting coeffs'
E_dash = double(getcondvects(dim)); %set of vertices of [0,1]^dim
E = E_dash(2:end,:); % set of non-zero vertices in E_dash
k_min = -10;
k_max = +10;
j_min = 0;
j_max = 3;

E = E_dash; %test
if( (size(C_ejk,1) == size(E,1)) && (size(C_ejk,2) == numel(j_min:j_max)) && ...
	(size(C_ejk,3)==numel(k_min:k_max)^2)) %if C_ejk already exists
'coeffs loaded'
else
'computing coeffs'
C_ejk = getCoefficientsVector(grid_x,dist_array_xy,dx,j_min,j_max,k_min,k_max,E,0);
end
%% test on grid
'visualizing approximation on a smaller range'
grid_x_half = xmin/2:dx:xmax/2;
grid_y_half = ymin/2:dy:ymax/2;

[msg_x,msg_y] = meshgrid(grid_x_half,grid_y_half);

%computed distance as well
if(exist('Staliro')>0)
    a = arrayfun(@(x,y) SignedDist([x;y],P.A,P.b), msg_x(:), msg_y(:));
   dist_array_xy_restrict = reshape(a,numel(grid_x_half),numel(grid_y_half));
    
end


if(1)
   [msg_x,msg_y] = meshgrid(grid_x_half,grid_y_half);
   a = arrayfun(@(x,y) getWavApprox_vector([x;y],C_ejk,k_min,k_max,j_min,j_max,E), msg_x(:), msg_y(:));
   fhat_array = reshape(a,numel(grid_x_half),numel(grid_y_half)); 
   
end
if(dim==2)
    figure;
   mesh(grid_x_half,grid_y_half,fhat_array); 
   figure;
   mesh(grid_x_half,grid_y_half,dist_array_xy_restrict-fhat_array);
end

%% save
wavparams.C_ejk = C_ejk;
wavparams.E = E;
wavparams.k_min = k_min;
wavparams.k_max = k_max;
wavparams.j_min = j_min;
wavparams.j_max = j_max;
save('2d_data.mat','C_ejk','j_min','j_max','k_min','k_max','xmin','xmax',...
    'dist_array_xy','dist_array_xy_restrict','E','dx','wavparams','P','dim',...
    'grid_x','fhat_array','dim','grid_y','grid_x_half','grid_y_half');
