function [ c , ceq ] = con_MLS_M( x )

    c(1) = norm( [ x ; -sum(x) ] , 1 ) - 4;
    
    load('QS_MLS_M.mat');
    c(2) = x' * Q * x - 1e-6; 
    
    ceq = [];

end

