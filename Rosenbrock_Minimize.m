clear
clc

%%Minimize Rosenbrock Function

%% f(x,y)= (1-x)^2+100(y-x^2)^2

%using steepest descent
%x(k+1)=x(k)+alpha(k)*P(k)

%P(k)=-grad(F(k))

%% initialize variables

tol=1e-10; %define tolerance
maxk=100000; %define max iterations
constraints=1; %define number of constraints

A=zeros(2,1); %intialize array of results
Convergence=zeros(2,1); %initialize convergence array

A(1,1)=1.5;   %initialize first x
A(2,1)=1.5;   %initialize first y
B=A;
C=A;
D=A;

L=100*ones(constraints,1);


for scheme=1:4
%     if scheme==1
%         alpha=.001;
%         [A,ConvergenceA]=steepest(alpha,A,Convergence,tol,maxk);
%     end
    
%     if scheme==2
%         [B,ConvergenceB]=conjugate(B,Convergence,tol,maxk);
%     end
    
%     if scheme==3
%         [C,L,ConvergenceC]=NewtonConstraint(C,L,Convergence,tol,maxk);
%     end
    
     if scheme==4
         [D,ConvergenceD]=QuadPenalty(D,Convergence,tol,maxk);
     end
end

figure;
%ConvergenceA(1,:),ConvergenceA(2,:)
semilogy(ConvergenceD(1,:),ConvergenceD(2,:));
legend('Quadratic Penalty: Newton Method','location','southwest')
xlabel('Iteration')
ylabel('Residual')
title('Convergence of Residual for constrained Rosenbrock Function')



[xplot,yplot]=meshgrid(-2:.2:2,-2:.2:2);
Z=(1-xplot).^2+100.*(yplot-xplot.^2).^2;

% figure;
% contour(xplot,yplot,Z,40);
% axis([-1 1 -1 1]);
% hold on;
% plot(A(1,:),A(2,:),...
%     'MarkerEdgeColor','r',...
%     'Marker','.',...
%     'LineStyle','-',...
%     'Color','r')
% legend('Rosenbrock Contour','Iteration Path','location','southwest')
% xlabel('x')
% ylabel('y')
% title('Steepest Descent Iteration for Rosenbrock Function')

% 
% figure;
% contour(xplot,yplot,Z,40);
% axis([-1 1 -1 1]);
% hold on;
% plot(D(1,:),D(2,:),...
%     'MarkerEdgeColor','r',...
%     'Marker','.',...
%     'LineStyle','-',...
%     'Color','r')
% legend('Rosenbrock Contour','Iteration Path','location','southwest')
% xlabel('x')
% ylabel('y')
% title('Conjugate Gradient Method Iteration for Rosenbrock Function')

% figure;
% contour(xplot,yplot,Z,80);
% axis([-2 2 -2 2]);
% hold on;
% plot(C(1,:),C(2,:),...
%     'MarkerEdgeColor','r',...
%     'Marker','.',...
%     'LineStyle','-',...
%     'Color','r')
% legend('Rosenbrock Contour','Iteration Path','location','southwest')
% xlabel('x')
% ylabel('y')
% title('Newton Method Iteration for Constrained Rosenbrock Function')

 figure;
 contour(xplot,yplot,Z,40);
 axis([-2 2 -2 2]);
 hold on;
 plot(D(1,:),D(2,:),...
     'MarkerEdgeColor','r',...
     'Marker','.',...
     'LineStyle','-',...
     'Color','r')
 legend('Rosenbrock Contour','Iteration Path','location','southwest')
 xlabel('x')
 ylabel('y')
 title('Quadtratic Penalty Newton Method Iteration for Rosenbrock Function')
