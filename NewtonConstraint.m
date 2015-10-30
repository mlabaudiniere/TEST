function [A,lambda,Convergence] = NewtonConstraint(A,lambda,Convergence,tol,maxk)
%Constrained Newton method

for i=1:maxk

x=A(1,i);       %x
y=A(2,i);       %y

c=1-x^2-y^2;        %calculate Constraint

a(1,1)=-2*x;        %Compute GradC
a(1,2)=-2*y;

grad(1,i)=-400*x*(y-x^2)+2*x-2;     %compute Fx
grad(2,i)=200*(y-x^2);              %compute Fy

dxL=grad(:,i)-a'*lambda(:,i);       %Compute GradF-A*Lambda

B=[-dxL;-c];                        %Ax=B (RHS)

dldx2=[-400*y+1200*x^2+2+2*lambda(:,i), -400*x; -400*x, 200+2*lambda(:,i)];     %Compute Lxx

Fprime=[dldx2, -a';a, 0];       %A in Ax=B



Convergence(1,i)=i;         %Compute Residual
Convergence(2,i)=norm(dxL); 


%Check Convergence
if Convergence(2,i)<=tol
        break
 end

P(:,i)=Fprime\B;        %Compute Search Direction
Px(:,i)=[P(1,i);P(2,i)];    %Split in Search direction for x and lambda
Pl(:,i)=[P(3,i)];

[alpha]=backtracking(A(:,i),[grad(:,i);c],P(:,i),lambda);       %Calculate alpha with backtracking function
A(:,i+1)=A(:,i)+alpha*Px(:,i);      %Next Point
lambda(:,i+1)=lambda(:,i)+alpha*Pl(:,i);    %Next Lambda

end
end

