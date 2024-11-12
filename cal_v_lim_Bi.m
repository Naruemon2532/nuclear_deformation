function velocity_lim_Bi = cal_v_lim_Bi(N,M,f,kon,koff,p,q,delta)

b=q/p;
a= factorial(N)*factorial(M);
kons = kon/M;

for c=1:N+1
    
    pp0(c) = (a/(factorial(c-1)*factorial(N-(c-1))*factorial(M-(c-1))))*(kons/koff)^(c-1);
    
end


P0= 1/sum(pp0);


for s=1:N+1
    Prob_Bi_lim(s)=(1/factorial(s-1))*(kons/koff)^(s-1)*(factorial(N)/factorial(N-(s-1)))*...
        (factorial(M)/factorial(M-(s-1)))*P0;
end

for i =1:N
    Pn(i)=Prob_Bi_lim(i+1);
    top(i) = p*(1-(exp(f))*(b^i))*(1-b);
    bottom(i) = (exp(f*delta))*(1-b)+ (exp(f)*(b-b^i));
    sub_v_lim_bi(i) = (Pn(i)/(1-P0))*(top(i)/bottom(i));
end

velocity_lim_Bi  = sum(sub_v_lim_bi);
end