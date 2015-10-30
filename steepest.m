function [A,Convergence] = steepest(alpha,A,Convergence,tol,maxk)


for i=1:maxk
    
    x=A(1,i);       %what is x?
    y=A(2,i);       %what is y?
    
    gradx=-400*x*(y-x^2)+2*x-2;     %compute fx
    grady=200*(y-x^2);              %compute fy
    
    A(1,i+1)=A(1,i)-alpha*gradx;    %compute next value
    A(2,i+1)=A(2,i)-alpha*grady;    %compute next value
    
    Convergence(1,i)=i;
    Convergence(2,i)=sqrt(gradx^2+grady^2);      %compute magnitude of gradient
    
    if Convergence(2,i)<=tol   %close enough to 0?
        break
    end
    
    

end
end

