function y=I1(a,b,Li,lT,sigmaT2);
y=0.5*sqrt(pi)*sigmaT2*lT* exp(-b.*b/(lT*lT)) .* (erf((Li+a)/lT)-erf(a/lT));