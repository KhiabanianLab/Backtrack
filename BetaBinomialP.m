function logP = BetaBinomialP (m, Dm, n, Dn)

minX = 50;

logP = log ((Dn+1)/(Dm+Dn+1));

if (m>=minX && m<(Dm-minX)) 
    logP = logP + ((Dm+0.5)*log(Dm) - (Dm-m+0.5)*log(Dm-m) - (0.5*log(2*pi()) + (m+0.5)*log(m)));
elseif (m>0 && m<minX)
    logP = logP + ((Dm+0.5)*log(Dm) - (Dm-m+0.5)*log(Dm-m) - log(factorial(m)) - m);
elseif (m~=0 && m~=Dm)
    logP = logP + ((Dm+0.5)*log(Dm) - (m+0.5)*log(m) - log(factorial(Dm-m)) - (Dm-m));
end

if (n>=minX && n<(Dn-minX))
    logP = logP + ((Dn+0.5)*log(Dn) - (Dn-n+0.5)*log(Dn-n) - (0.5*log(2*pi()) + (n+0.5)*log(n)));
elseif (n>0 && n<minX)
    logP = logP + ((Dn+0.5)*log(Dn) - (Dn-n+0.5)*log(Dn-n) - log(factorial(n)) - n);
elseif (n~=0 && n~=Dn)
    logP = logP + ((Dn+0.5)*log(Dn) - (n+0.5)*log(n) - log(factorial(Dn-n)) - (Dn-n));
end

t=m+n;
Dt=Dm+Dn;

if (t>=minX && t<(Dt-minX))
    logP = logP - ((Dt+0.5)*log(Dt) - (Dt-t+0.5)*log(Dt-t) - (0.5*log(2*pi()) + (t+0.5)*log(t)));
elseif (t>0 && t<minX)
    logP = logP - ((Dt+0.5)*log(Dt) - (Dt-t+0.5)*log(Dt-t) - log(factorial(t)) - t);
elseif (t~=0 && t~=Dt)
    logP = logP - ((Dt+0.5)*log(Dt) - (t+0.5)*log(t) - log(factorial(Dt-t)) - (Dt-t));
end
