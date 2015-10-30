function [A,Convergence] = conjugate(A,Convergence,tol,maxk)
%Conjugate Gradient Approach

x=A(1,1);       %Initial x
y=A(2,1);       %Initial y
    
grad(1,1)=-400*x*(y-x^2)+2*x-2;     %compute fx
grad(2,1)=200*(y-x^2);              %compute fy
  
%Q=[-400*y+1200*x^2+2 -400*x; -400*x 200];  %Hessian
   
P(1,1)=-grad(1,1);          %initial search direction
P(2,1)=-grad(2,1);          %initial search direction

sigma(1)=0;

for i=1:maxk
    
   % Q=[-400*y+1200*x^2+2 -400*x; -400*x 200];  %Hessian
    
    [alpha]=backtracking(A(:,i),grad(:,i),P(:,i));
    
    A(:,i+1)=A(:,i)+alpha*P(:,i);
    
    x=A(1,i+1);
    y=A(2,i+1);
    
    grad(1,i+1)=-400*x*(y-x^2)+2*x-2;
    grad(2,i+1)=200*(y-x^2); 
    
    Convergence(1,i)=i;
    Convergence(2,i)=sqrt(grad(1,i+1)^2+grad(2,i+1)^2);
    
    if Convergence(2,i)<=tol
        break
    end
        
    sigma(i+1)=transpose(grad(:,i+1))*(grad(:,i+1)-grad(:,i))/(transpose(P(:,i))*(grad(:,i+1)-grad(:,i)));  %Heines-Steiffel
    P(:,i+1)=-grad(:,i+1)+sigma(i+1)*P(:,i);
        
    end

end



