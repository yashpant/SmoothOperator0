function psi = RickerWavelet(x,sigma)

if(nargin==1)
   sigma = .25; 
end

psi = (2/(sqrt(3*sigma)*(pi^(1/4))))*(1-((x^2)/(sigma^2)))* ...
    exp((-x^2)/(2*sigma^2));