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

original = V.mat \ [0;0;0;1];          % calculate the center of brain
original = original(1:3)';                    
 
IDX = find( data == 1 | data == 4 | data == 52 ); % 4 R SFG; 1 L PreCG; 52 R OG
IDX1{1} = find( data == 1 );
IDX1{2} = find( data == 4 );
IDX1{3} = find( data == 52 );

IDX_ex = 1 : size( A{1}(:,1:3) , 1 );  
IDX_ex(IDX) = [];                      % Index of coor outside selected ROI 

Vm_in = length( find( atlas == 1 | atlas == 4 | atlas == 52 ) ) * 1e-9 / length( IDX );  % 1e-9, mm3->m3
Vm_out = ( Vm - Vm_in ) * 1e-9 / length( IDX_ex );

% normailized radial direction vector of each coor in A
radial_norm = ( A{1}(:,1:3) - repmat( original , size( A{1} , 1 ) , 1 ) ) ./ repmat( sqrt( sum( ( A{1}(:,1:3) - repmat( original , size( A{1} , 1 ) , 1 ) ) .^ 2 , 2 ) ) , 1 , 3 );

% Using Eq.22 calculating beta
load c_PreCG;  
fv(1) = fval;
load c_SFG;
fv(2) = fval;
load c_OG;
fv(3) = fval;

beta(1) = fv(1) / sum( fv ); 
beta(2) = fv(2) / sum( fv ); 
beta(3) = fv(3) / sum( fv ); 

% Calculating W for each ROI
W1 = zeros( 1 , length( A ) );
for i = 1 : length( IDX1{1} )
    
    J = zeros( 3 , length( A ) );
    for j = 1 : length( A )
        J(:,j) = A{j}(IDX1{1}(i),4:6)';
    end
    W1 = W1 + radial_norm(IDX1{1}(i),:) * J;
    
end

W2 = zeros( 1 , length( A ) );
for i = 1 : length( IDX1{2} )
    
    J = zeros( 3 , length( A ) );
    for j = 1 : length( A )
        J(:,j) = A{j}(IDX1{2}(i),4:6)';
    end
    W2 = W2 + radial_norm(IDX1{2}(i),:) * J;
    
end

W3 = zeros( 1 , length( A ) );
for i = 1 : length( IDX1{3} )
    
    J = zeros( 3 , length( A ) );
    for j = 1 : length( A )
        J(:,j) = A{j}(IDX1{3}(i),4:6)';
    end
    W3 = W3 + radial_norm(IDX1{3}(i),:) * J;
    
end

save WME_1 W1 W2 W3;

% S = [0.3,0.6,0.3]
W1 = W1 / beta(1) * 0.3;
W2 = W2 / beta(2) * 0.6;
W3 = W3 / beta(3) * 0.3;
save WME_M W1 W2 W3;

% Calculation of Q
Q = zeros( length( A ) , length( A ) );
for i = 1 : length( IDX_ex )
    
    J = zeros( 3 , length( A ) );
    for j = 1 : length( A )
        J(:,j) = A{j}(IDX_ex(i),4:6)';
    end
    Q = Q + J' * J * Vm_out;
    
end
save QS_ME_M Q;

% 100 times inter-points methods
c = zeros( length( A ) , 100 ); % c current
fval = zeros( 100 , 1 );        % fval maximum intensity
flag = fval;                    % flag: 2 corresponds to successful solving the optimization problem
pow = flag;                     % pow: power outside selected ROI
options = optimset( 'Tolcon' , 1e-9 , 'Tolfun' , 1e-10 , 'Display' , 'off' , 'MaxIter' , 10000 , 'MaxFunEval' , inf );  
parfor i = 1 : 100
    [c(:,i),fval(i),flag(i)] = fmincon( @fun_ME_M , randn( length( A ) , 1 ) , [] , [] ,[] , [] , -2 * ones( length( A ) , 1 ) , 2 * ones( length( A ) , 1 ) , @con_ME_M , options );
    pow(i) = c(:,i)'*Q*c(:,i);
end

IDXX = find( flag == 2 );
ID = find( fval == min( fval(IDXX) ) );
cc = c(:,IDXX);

%% current norm
c = c(:,ID);
if c( abs(c) == max( abs(c) ) ) < 0
    c = -c;
end
fval = fval(ID);
pow = pow(ID);

save c_ME_M_W13 cc c fval pow;

toc;