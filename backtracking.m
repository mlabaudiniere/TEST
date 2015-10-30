function [alpha] = backtracking(A,grad,P,mu)
%backtracking algorithm

rho=.9;      %backtracking coefficient
c1=.1;       %backtracking coefficient
alpha=1;   %initial alpha

x=A(1);
y=A(2);

for i=1:25
    
    fplus=(1-(x+alpha*P(1)))^2+100*((y+alpha*P(2))-((x+alpha*P(1))^2))^2+mu/2*(1-(x+alpha*P(1))^2-(y+alpha*P(2))^2)^2;
    f=(1-x)^2+100*(y-x^2)^2+mu/2*(1-x^2-y^2)^2;
    del=c1*alpha*grad'*P;
    
    if fplus <= f+del
        break
    else
        alpha=rho*alpha;
    end
end

