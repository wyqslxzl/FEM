function f = fun_ME_M( x )
 
    load('WME_M.mat');
    f = -abs( (W1+W2+W3) * x );

end

