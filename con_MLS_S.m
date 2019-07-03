function [ c , ceq ] = con_MLS_S( x )

    load('QS_PreCG.mat');
    c(1) = norm( [ x ; -sum(x) ] , 1 ) - 4;
    c(2) = x' * Q * x - 1e-6;
    ceq = [];

end

