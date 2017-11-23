% Funzione algoritmo di bisezione

function [xvect,it] = bisez(a,b,toll,fun)

xvect = [];
it = 0;

nmax = ceil(log2((b-a)/toll)-1);
err = toll+1;
while (err>toll || it<nmax)
    
    it = it+1;
    xk = (b+a)/2;
    
    if (fun(xk) == 0)
        zero = x;        
        return
        
    elseif (fun(a)*fun(xk)<0)
        a=a;
        b=xk;
        
    else
        a=xk;
        b=b;
    end
    
    err = abs(fun(xk));
    xvect = [xvect;xk];
end

    
    

