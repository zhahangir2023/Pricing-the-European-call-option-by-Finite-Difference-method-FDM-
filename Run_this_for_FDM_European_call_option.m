

clear all
close all
format long



refinement = 10;
JJ = size(refinement);
Nt = [2^refinement];
N_nodes = 2.^refinement;
for i = 1
[U_FDM,S_FDM,dtau] = Vanilla_call_option_example_FDM_for_thesis(N_nodes(i,1),Nt(i,1));




end

[Vexact, S_exact, Delta_exact, Gamma_exact] = Vanilla_Option_Call_Exact_BS_inS(0.2,0.01);

figure;

plot(S_FDM,U_FDM(:,end),'o-',S_exact,Vexact,'b-',S_FDM,U_FDM(:,1),'g--')
xlabel('Stock')
ylabel('Option price')
title('European call option')

legend('FDM','Exact','Pay-off')

axis([40 140 0 60])
% 
[x_FDM,y_FDM] = meshgrid(S_FDM,1:-dtau:0);
figure;
surf(x_FDM,y_FDM,U_FDM')
axis([40 200 0 1 -0.2 120])
xlabel('Stock')
ylabel('Time')
zlabel('Option')
title('European call option by FDM')
