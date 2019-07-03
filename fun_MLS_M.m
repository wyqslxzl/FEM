function f = fun_MLS_M( x )
 
    load('WMLS_M.mat');
    f = x' * WW * x - 2 * YW * x + YY;

end
