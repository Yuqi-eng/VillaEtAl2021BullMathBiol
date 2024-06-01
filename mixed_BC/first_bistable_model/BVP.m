clc
close all

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')
%% Numerical set up
par.K = 100;                    % Number of spatial grid cells
par.L = 1;                      % Domain length
x = linspace(0,par.L,par.K);    % Discretise spatial domain

%% Create initial guess
solinit = bvpinit(x, @guess);

%% Solve
options = bvpset(Stats="on");
sol = bvp4c(@bvpfcn, @bcfcn, solinit, options);

%%% reshape output
n = sol.y(1,:);
p = sol.y(3,:);
u = sol.y(5,:);
sig = sol.y(6,:);
% plot(sol.x,[n;p;u]);
figure('Units','normalized','Position',[0 0 0.5 0.45])
subplot(1,4,1)
plot(sol.x,n)
title('$n(t,x)$')
axis square
subplot(1,4,2)
plot(sol.x,p)
title('$\rho(t,x)$')
axis square
subplot(1,4,3)
plot(sol.x,u)
title('$u(t,x)$')
axis square
subplot(1,4,4)
plot(sol.x,sig)
title('$\sigma(t,x)$')
axis square
pic_name = ['BVP_res' '.png'];
saveas(gcf,pic_name);

%% Main function implementing the model
function f = bvpfcn(x,y)
    %%% Parameter values
    a = 4;
    D2 = 0.1;
    d2 = 1;
    E1 = 1.4e8;
    tau1 = 0.4;
    
    %%% hill function traction force term
    N0 = 1.25;
    k2 = 4;
    fn = (y(1)^k2)/(N0^k2 + y(1)^k2);

    sig0 = 0.3;
    k1 = 4;
    fsig = (y(6)^k1)/(sig0^k1 + y(6)^k1);

    % y(1) = n, y(2) = n_x, y(3) = rho, y(4) = rho_x, y(5) = u, y(6) =
    % sigma
    f = [y(2)
        y(1)-1-a*fsig
        y(4)
        (1/D2)*(d2*y(3)-y(1))
        (1/E1)*(y(6)-tau1*y(3)*fn)
        0];
end

%% Helper functions

function res = bcfcn(ya,yb) % boundary conditions
res = [yb(1)-1
       ya(2)
       yb(3)-1
       ya(4)
       ya(5)
       yb(5)];
end

function g = guess(x)
g = [3
    0
    3
    0
    0
    2];
end

% function g = guess(x)
% g = [2
%      0
%      2
%      0
%      0
%      5];
% end