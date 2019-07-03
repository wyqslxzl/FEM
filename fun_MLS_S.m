function f = fun_MLS_S( x )
 
    load('WS_PreCG_MLS.mat');
    f = x' * WW * x - 2 * YW * x + YY;

end
