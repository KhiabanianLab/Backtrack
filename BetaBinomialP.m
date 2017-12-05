function logP = BetaBinomialP (m, Dm, n, Dn)

logP = log ((Dn+1)/(Dm+Dn+1));

if (m>0 && m<Dm) 
    logP = logP + ((Dm+0.5)*log(Dm) - (Dm-m+0.5)*log(Dm-m) - (0.5*log(2*pi()) + (m+0.5)*log(m)));
end

if (n>0 && n<Dn)
    logP = logP + ((Dn+0.5)*log(Dn) - (Dn-n+0.5)*log(Dn-n) - (0.5*log(2*pi()) + (n+0.5)*log(n)));    
end

logP = logP - ((Dm+Dn+0.5)*log(Dm+Dn) - (Dm+Dn-m-n+0.5)*log(Dm+Dn-m-n) - (0.5*log(2*pi()) + (m+n+0.5)*log(m+n)));
