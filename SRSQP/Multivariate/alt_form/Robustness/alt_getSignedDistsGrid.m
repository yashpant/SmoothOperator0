function dist_array_xy =  alt_getSignedDistsGrid(grid_x,P)
    dim = size(P.A,2);
    M1 = zeros(numel(grid_x)^dim,dim);
    parfor i = 1:(numel(grid_x))^dim
     M1(i,:) = permn2(grid_x,dim,i); 
    end
    %M1 = permn(grid_x,dim);	
	%C = num2cell(M1,2);
	%size(C)
	%keyboard;	
    a_cell = zeros((numel(grid_x))^dim,1);
    'prelim work done, now getting dist'
    parfor i=1:(numel(grid_x))^dim
	% i/((numel(grid_x))^dim)
       a_cell(i) = SignedDist(M1(i,:)',P.A,P.b);
    end
    
    if(dim~=1)
	dist_array_xy_vec = reshape(a_cell, numel(grid_x)*ones(1,dim));
	dist_array_xy = dist_array_xy_vec;		
    else
       dist_array_xy = a_cell; 
    end
