function [A,Convergence] = quasi(A,Convergence,tol,maxk)
%Quasi-Newton Method

x=A(1,1);       %x
y=A(2,1);       %y    

grad(1,1)=-400*x*(y-x^2)+2*x-2;     %compute fx
grad(2,1)=200*(y-x^2);              %compute fy

Q=[-400*y+1200*x^2+2 -400*x; -400*x 200];  %Hessian


for i=1:maxk    

Convergence(1,i)=i;
Convergence(2,i)=sqrt(grad(1,i)^2+grad(2,i)^2);

 if Convergence(2,i)<=tol
        break
 end
    
P(:,i)=-Q*grad(:,i);


[alpha]=backtracking(A(:,i),grad(:,i),P(:,i));

A(:,i+1)=A(:,i)+alpha*P(:,i);

dA=A(:,i+1)-A(:,i);

x=A(1,i+1);       %x
y=A(2,i+1);       %y      
    
grad(1,i+1)=-400*x*(y-x^2)+2*x-2;     %compute fx
grad(2,i+1)=200*(y-x^2);              %compute fy

dgrad=grad(:,i+1)-grad(:,i);

Q=Q+(1+(dgrad'*Q*dgrad)/(dA'*dgrad))*(dA*dA')/(dA'*dgrad)-((Q*dgrad*dA')+(Q*dgrad*dA')')/(dA'*dgrad);

end

