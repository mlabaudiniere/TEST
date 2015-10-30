function [A,Convergence] = newton(A,Convergence,tol,maxk)
%Newton method

for i=1:maxk

x=A(1,i);       %x
y=A(2,i);       %y

grad(1,i)=-400*x*(y-x^2)+2*x-2;     %compute Fx
grad(2,i)=200*(y-x^2);              %compute Fy

Q=[-400*y+1200*x^2+2 -400*x; -400*x 200];  %Hessian

Convergence(1,i)=i;
Convergence(2,i)=sqrt(grad(1,i)^2+grad(2,i)^2);

 if Convergence(2,i)<=tol
        break
 end

P(:,i)=-inv(Q)*grad(:,i);
 
[alpha]=backtracking(A(:,i),grad(:,i),P(:,i));

A(:,i+1)=A(:,i)+alpha*P(:,i);

end
end

