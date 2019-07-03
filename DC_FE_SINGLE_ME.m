clear;clc;
tic;

load DC_64;   % Coefficient Matrix of 63 freedom electrodes.

%% Region Optimization
%%% file needs absolute address %%%
file = '\HC001_AAL.nii';   % atlas
V = spm_vol( file );
atlas = spm_read_vols( V );
clear file;
%%% file needs absolute address %%%
file(1,:) = '\HC001_brain_seg_1.nii';  % gray matter
file(2,:) = '\HC001_brain_seg_2.nii';  % white matter
V1 = spm_vol( file );
data1 = spm_read_vols( V1 );
Vm = sum( data1(:) );                  % brain = gray + white

coor = round( A{1}(:,[2,1,3])' ) + 1;  % A; Coeffient Matrix
                                       % including x,y,z,and value of E
IDX = sub2ind( V.dim , coor(1,:) , coor(2,:) , coor(3,:) );
data = atlas( IDX );                   % coor mapping into atlas

IDX = find( data == 1 );               % 4 R SFG; 1 L PreCG; 52 R OG

original = V.mat \ [0;0;0;1];          % calculate the center of brain
original = original(1:3)';

IDX_ex = 1 : size( A{1}(:,1:3) , 1 );  
IDX_ex(IDX) = [];                      % Index of corr outside selected ROI

Vm_in = length( find( atlas == 1 ) ) * 1e-9 / length( IDX );  % 1e-9, mm3->m3
Vm_out = ( Vm - Vm_in ) * 1e-9 / length( IDX_ex );

% normailized radial direction vector of each coor in A
radial_norm = ( A{1}(:,1:3) - repmat( original , size( A{1} , 1 ) , 1 ) ) ./ repmat( sqrt( sum( ( A{1}(:,1:3) - repmat( original , size( A{1} , 1 ) , 1 ) ) .^ 2 , 2 ) ) , 1 , 3 );

% Calculation of W
W = zeros( 1 , length( A ) );
for i = 1 : length( IDX )
    
    J = zeros( 3 , length( A ) );
    for j = 1 : length( A )
        J(:,j) = A{j}(IDX(i),4:6)';
    end
    W = W + radial_norm(IDX(i),:) * J;
    
end
save WS_PreCG_ME W;

% Calculation of Q
Q = zeros( length( A ) , length( A ) );
for i = 1 : length( IDX_ex )
    
    J = zeros( 3 , length( A ) );
    for j = 1 : length( A )
        J(:,j) = A{j}(IDX_ex(i),4:6)';
    end
    Q = Q + J' * J * Vm_out;
    
end
save QS_PreCG Q;

% 100 times inter-points methods
c = zeros( length( A ) , 100 ); % c current
fval = zeros( 100 , 1 );        % fval maximum intensity
flag = fval;                    % flag: 2 corresponds to successful solving the optimization problem
pow = flag;                     % pow: power outside selected ROI
options = optimset( 'Tolcon' , 1e-8 , 'Tolfun' , 1e-10 , 'Display' , 'off' , 'MaxIter' , inf , 'MaxFunEval' , inf );  
parfor i = 1 : 100 % parallel
    [c(:,i),fval(i),flag(i)] = fmincon( @fun_ME_S , randn( length( A ) , 1 ) , [] , [] ,[] , [] , -2 * ones( length( A ) , 1 ) , 2 * ones( length( A ) , 1 ) , @con_ME_S , options );
    pow(i) = c(:,i)'*Q*c(:,i);
end

% genetic algorithm
% options = optimset( 'Tolcon' , 1e-8 , 'Tolfun' , 1e-10 );  
% [c,fval,flag] = ga( @fun3 , length( A ) , [] , [] ,[] , [] , -2 * ones( 1 , length( A )  ) , 2 * ones( 1, length( A ) ) , 'fun22' , options );
% pow = c*Q*c';

IDXX = find( flag == 2 );
ID = find( fval == min( fval(IDXX) ) );
cc = c(:,IDXX);

%% E norm
c = c(:,ID);
if c( abs(c) == max( abs(c) ) ) < 0
    c = -c;
end
fval = fval(ID);
pow = pow(ID);

save c_PreCG_ME c fval pow;

toc;