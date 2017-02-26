function wavparams = WavSignedDistVector(P,xmin,xmax,dx,viz,optional)

if(nargin==4)
    viz = 0;
    optional.filename = [];
    optional.savefile = 'newdata.mat';
elseif(nargin==5)
    optional.filename = [];
    optional.savefile = 'newdata.mat';
end
%%
ymax = xmax;
ymin = xmin;
dy = dx;

dim = size(P.A,2);

grid_x = xmin:dx:xmax;
grid_y = ymin:dy:ymax;
tic;

% load data if it exits
if(exist(optional.filename)>0)
    'loading data'
    load(optional.filename);
end

%% get signed dist
if(exist('dist_array_xy') &&  sum(size(dist_array_xy) == size(meshgrid(grid_x,grid_y))) == dim)
    
    'already have data'
else
    %keyboard;
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
j_max = 2;

if(exist('C_00k','var') && exist('D_ejk','var')) %if some params loaded
    
    if( (size(C_00k,3) == numel(k_min:k_max)^2) && (size(D_ejk,1) == size(wavparams.E_dash,1)) ...
            && (size(D_ejk,2) == numel(j_min:j_max)) && (size(D_ejk,3) == numel(k_min:k_max)^2))
        'params loaded'
        
    end
    
else
    'computing coeffs'
    [C_00k,D_ejk] = alt_getCoefficientsVector_genable(grid_x,dist_array_xy,dx,j_min,j_max,k_min,k_max,E_dash,0);
end
%% test on grid
if(viz)
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
        a = arrayfun(@(x,y) alt_getWavApprox_vector([x;y],C_00k,D_ejk,k_min,k_max,j_min,j_max,E_dash), msg_x(:), msg_y(:));
        fhat_array = reshape(a,numel(grid_x_half),numel(grid_y_half));
        
    end
    
end
%%
if(dim==2 && viz)
    figure;
    mesh(grid_x_half,grid_y_half,fhat_array);
    figure;
    mesh(grid_x_half,grid_y_half,dist_array_xy_restrict-fhat_array);
end

%% save
wavparams.D_ejk = D_ejk;
wavparams.C_00k = C_00k;
wavparams.E = E;
wavparams.E_dash = E_dash;
wavparams.k_min = k_min;
wavparams.k_max = k_max;
wavparams.j_min = j_min;
wavparams.j_max = j_max;


save(optional.savefile,'D_ejk','C_00k','j_min','j_max','k_min','k_max','xmin','xmax',...
    'dist_array_xy','dist_array_xy_restrict','E','dx','wavparams','P','dim',...
    'grid_x','fhat_array','dim','grid_y','grid_x_half','grid_y_half');
