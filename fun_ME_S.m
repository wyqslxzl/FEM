function f = fun_S( x )
 
    load('WS_PreCG_ME.mat');
    f = -abs( W * x );

end

