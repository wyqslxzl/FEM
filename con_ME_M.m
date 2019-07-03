function [ c , ceq ] = con_ME_M( x )

    c(1) = norm( [ x ; -sum(x) ] , 1 ) - 4;
    
    load( 'QS_ME_M.mat' );
    c(2) = x' * Q * x - 1e-6;    
    
    load( 'WME_M.mat' );
    c(3) = std( abs( [W1*x,W2*x,W3*x] ) ) - 0.1;
    
    ceq = [];

end

