function [A,Convergence] = QuadPenalty(A,Convergence,tol,maxk)
%Constrained Optimization 

mu=1;

for i=1:maxk

x=A(1,i);       %x
y=A(2,i);       %y

grad(1,i)=-400*x*(y-x^2)+2*x-2-mu*2*x*(1-x^2-y^2);     %compute Fx
grad(2,i)=200*(y-x^2)-mu*2*y*(1-x^2-y^2);              %compute Fy

Q=[-400*y+1200*x^2+2+6*mu*x^2-2*mu+2*mu*y^2 -400*x+4*x*y*mu; -400*x+4*x*y*mu 200+6*mu*y^2-2*mu+2*mu*x^2];  %Hessian

Convergence(1,i)=i;
Convergence(2,i)=sqrt(grad(1,i)^2+grad(2,i)^2); %Residual

P(:,i)=-inv(Q)*grad(:,i); %Newton Search Direction
 
[alpha]=backtracking(A(:,i),grad(:,i),P(:,i),mu); %compute alpha with backtracking

A(:,i+1)=A(:,i)+alpha*P(:,i); %compute next points


 if Convergence(2,i)<=tol
     
     if mu>5000
        break
     else
         mu=mu*1.5 %if converged, increase mu unless large
         A(:,i+1)
     end
 end


end


end

