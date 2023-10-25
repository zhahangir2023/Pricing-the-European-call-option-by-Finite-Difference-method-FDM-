


function [Umat,S_axis,dtau] = Vanilla_call_option_example_FDM_for_thesis(Nx,Nt)

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PRE-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this example, we consider the European Vanilla call option
%  V_t = \sigma^2/2* V_xx + (r-\sigma^2/2)*V_x -r*V
% I.C.: V(x,0) = max(exp(x)*S_int -K,0)
% B.C.: V(0,t) = 0; V(S_max,t) = exp(x)*S_int -K;


%%%%%%% Data for the problem



T = 1;



sigma = 0.2;
K_tilda = 100;
kappa = 1;
r = 0.05;
Sint = 100;
% Spatial mesh information
% Asset S-coordinate
% xmin = 0;
% xmax = 1;


% Number and size of temporal discretization
%
% Time integration parameter
% t-coordinate
% Nt = 12800;
t0 = 0;
t1 = T;
tau0 = T-t1;
tau1 = T-t0;
dtau = (tau1-tau0)/Nt;
theta = 1/2;



% Nx = 1024;
%
% Construct relevant matrices
% Dxx matrix : Mass matrix
K = Matrix_Dxx(Nx);
%

% Dx matrix :
N = Matrix_Dx(Nx);
%

% D matrix : Stiffness matrix
M = Matrix_D(Nx);

%%%%%%%%%%%%%%%%%%%%
Xmin = -2;
Xmax = 2;

dx = (Xmax - Xmin)/Nx;
x = (Xmin:dx:Xmax)';

Uinit = max(exp(x)*Sint - K_tilda , 0 ) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Umat = zeros(Nx+1,Nt+1);


% Initial time ste

% Boundary contributions
S_axis = exp(x)*Sint;

uFixed = [0 exp(x(end,1))*Sint - K_tilda]';

Umat(1,1) = uFixed(1,1);
Umat(2:end-1,1) = Uinit(2:end-1,1);
Umat(end,1) = uFixed(2,1);

% plot(S_axis,Umat(:,1))

% - for the CB PDE: use fully implicit scheme

% Block matrices for the time dependent problem
A11  = M - theta*dtau*(sigma^2/2/dx^2)*K - theta*dtau*( r -sigma^2/2)/(2*dx)*N + theta*dtau*r*M;
At11 = M + (1-theta)*dtau*(sigma^2/2/dx^2)*K + (1-theta)*dtau*(r-sigma^2/2)/(2*dx)*N - (1-theta)*dtau*r*M;



% Initial solution in the interior
U = Uinit(2:(end-1),1);
% Initial boundary condition
U0 = Umat(1,1);
Un = Umat(end,1);
% Boundary condition vectors
% - for the CB PDE:
fbcU0 = zeros(Nx-1,1);
fbcU0(1,1)   = (1/2*(sigma^2/dx^2) - (r - sigma^2/2)/(2*dx))*U0;
fbcU0(end,1) = (1/2*(sigma^2/dx^2) + (r - sigma^2/2)/(2*dx))*Un;

%
% Set up some vectors and matrices
% Boundary condition vectors at the new time level
fbcU1 = zeros(Nx-1,1);

for i = 1:Nt
    %
    tau = tau0 + i*dtau;
    %
    % (i) Crank-Nicolson (CN) Ste



    % Boundary conditions:
    fbcU1(1,1)   = (1/2*(sigma^2/dx^2) - (r - sigma^2/2)/(2*dx))*uFixed(1,1);
    fbcU1(end,1) = (1/2*(sigma^2/dx^2) + (r - sigma^2/2)/(2*dx))*uFixed(2,1);


    % Solve the CB PDE for U:
    fU = At11*U + theta*dtau*fbcU1 + (1-theta)*dtau*fbcU0;

    Unew = A11\(fU);


    U = Unew;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Switching roles of U0,Un,Bo,Bn,U
    % Update boundary condition vectors at the old time level m using m+1
    Umat(1,i+1) = uFixed(1,1);
    Umat(2:end-1,i+1) = Unew;
    Umat(end,i+1) = uFixed(2,1);
    fbcU0(1,1)   = (1/2*(sigma^2/dx^2) - (r - sigma^2/2)/(2*dx))*Umat(1,i+1);
    fbcU0(end,1) = (1/2*(sigma^2/dx^2) + (r - sigma^2/2)/(2*dx))*Umat(end,i+1);
    %




end
return


disp([num2str(toc),'  SOLVING THE SYSTEM'])
% 
% [xx, yy] = meshgrid(T:-dtau:0,S_axis);
% size(Umat)
% size(xx)
% size(yy)
% surf(xx,yy,Umat)
% [Vexact, S_exact,Delta_exact,Gamma_exact] = Vanilla_Option_Call_Exact_BS_inS(sigma,dtau);
% 
% figure
% 
% plot(S_axis, Umat(:,end),'bs',S_axis,Umat(:,1),'g--',S_exact,Vexact,'r');
% axis([50 300 -20 200])
% 
% legend('FDM','Pay-off','Linear B-S Exact solution')
% title('Vanilla European Call option price by IGA(NURBS) FEM vs Exact solution')
% open the file with write permission
disp([num2str(toc),'  SOLVING THE SYSTEM'])
S_position = Nx/2 + 1;
Umat(S_position,end)
% Exact solution  = 10.4506
%
% Results by FDM without Rannacher but with Crank Nicolson
%   10.441556227122483
%   10.448330078981066
%   10.450020406588441
%   10.450442793751337
%   10.450548378394123
%   10.450574773725130
%   10.450581372478545
%   10.450583022195358
%   10.450583434227791 corr. to Nx = 51 200
%   10.450583520828317 corr. to Nx = 102 400










%             %%%%%%%%%%%%%%%%%%% GREEKS FDM - Starts %%%%%%%%%%%%%%%%%%%
% 
%         Delta_FDM = zeros(Nx-1,Nt+1);
%         Gamma_FDM = zeros(Nx-1,Nt+1);
%         Theta_FDM = zeros(Nx-1,Nt+1);
%         % Theta = zeros(nE-1,Nt-1);
% 
%         for i = 2:Nx
%             Delta_FDM(i-1,:) = 1/(S_axis(i,1)).*((Umat(i+1,:) - Umat(i-1,:))/(2*dx));
%         end
%         % for i = 1:nE
%         %     Delta(i-,:) = ((U(i+1,:) - U(i,:))/(dx))*(1/Sint*exp(-x(i,1)));
%         % end
% %         [xdelta,ydelta] = meshgrid(S_axis(2:end-1,1),0:dtau:T);
% %         figure;
% %         mesh(xdelta,ydelta,Delta_FDM')
% % %         axis([0 400 0 5 -0.5 1.5])
% %         xlabel('Stock price')
% %         ylabel('Time')
% %         zlabel('Delta')
% %         % axis([0 250 0 7 80 250])
% %         % legend('P1 at $t = 0$')
% %         title('Delta by FDM')
%         for i = 2:Nx
%             Gamma_FDM(i-1,:) =  1/(S_axis(i,1))^2.*((Umat(i+1,:) - 2*Umat(i,:) + Umat(i-1,:))/(dx^2) -  (Umat(i+1,:) - Umat(i-1,:))/(2*dx) ) ;
%         end
% %         hold on;
% %         figure;
% %         mesh(xdelta,ydelta,Gamma_FDM')
% % %         axis([0 300 0 5 -0.01 0.45])
% %         xlabel('Stock price')
% %         ylabel('Time')
% %         zlabel('Gamma ')
% %         %  axis([0 300 0 5])
% %         % legend('P1 at $t = 0$')
% %         title('Gamma by FDM')
% %%%%% GREEKS comparison
% plot(S_exact(1,2:end-1), Delta_exact(2:end-1,1),'r',S_axis(2:end-1,:),Delta_FDM(:,end),'b--')   
%   axis([20 200  0 1])
%  plot(S_exact, Gamma_exact,S_axis(2:end-1,:),Gamma_FDM(:,end),'*')   
% %      axis([20 200  0 0.025])

