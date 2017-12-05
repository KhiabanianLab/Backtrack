function Backtrack_BetaBinomialP (file_in, file_out)

%f ile_in: tab delimited columns of variant allele depth and total depth across Nd samples

p_cut = 1e-4;
max_freq = 0.20;

Data = dlmread (file_in);

Nd = size(Data, 1);
p_out = zeros(1, Nd);
cutoffP = zeros(1, Nd);

for k=1:Nd
    m = Data(k,1);
    Dm = Data(k,2);
    n = sum(Data(1:Nd~=k & (Data(:,1)./Data(:,2)<max_freq)', 1));
    Dn = sum(Data(1:Nd~=k & (Data(:,1)./Data(:,2)<max_freq)', 2));
    
    x = 0:Dm; 
    p = zeros(1,length(x));
    for i=1:length(x)
        p(i) = BetaBinomialP(x(i), Dm, n, Dn);
    end

    x = 0:min(max(m+1, 2000), Dm);
    cumsumP = zeros(1,length(x));
    for i=1:length(x)
        cumsumP(i) = sum(exp(p(x>x(i))));
    end

    p_out(k) = min(sum(exp(p(x>m))), 1);
    
    if (Dn == 0) 
       cutoffP(k) = binoinv(1-p_cut, Dm, m/Dm);      
    else
       cutoffP(k) = min(x(cumsumP<(p_cut)));
    end  
end

p_fdr = mafdr (p_out, 'BHFDR', 'true');

fid = fopen (file_out, 'wt');
for k=1:Nd
    fprintf (fid, '%e\t%e\t%i\n', p_out(k), p_fdr(k), cutoffP(k));
end
fclose(fid);

