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
% plot(sol.x,[n;p;u]);
figure('Units','normalized','Position',[0 0 0.5 0.45])
subplot(1,3,1)
plot(sol.x,n)
title('$n(t,x)$')
axis square
subplot(1,3,2)
plot(sol.x,p)
title('$\rho(t,x)$')
axis square
subplot(1,3,3)
plot(sol.x,u)
title('$u(t,x)$')
axis square
pic_name = ['BVP_res' '.png'];
saveas(gcf,pic_name);

%% Main function implementing the model
function f = bvpfcn(x,y)
    %%% prameter definitions
    E = 1;
    s = 1;
    D = 0.01;
    Dp = 1e-3;
    dn = 0.5;
    an = 0.5;
    m = 0.1;
    dp = 0.1;
    tau = 0.5;
    a1 = 0.5;
    
    %%% hill function traction force term
    n0 = 1.45;
    k2 = 5;
    h1 = (y(1)^k2)/(n0^k2 + y(1)^k2);

    sig0 = 0.18;
    k1 = 5;
    fsig = (y(6)^k1)/(sig0^k1 + y(6)^k1);

    % y(1) = n, y(2) = n_x, y(3) = rho, y(4) = rho_x, y(5) = u, y(6) =
    % sigma
    f = [y(2)
        (1/D)*(dn*y(1)-an-a1*fsig)
        y(4)
        (1/Dp)*(dp*y(3)-m*y(1))
        (1/E)*(y(6)-tau*y(3)*h1)
        s*y(3)*y(5)];
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
g = [1
    0
    0
    0
    0
    0
    ];
end

% function g = guess(x)
% g = [1.5
%      0
%      1.5
%      0
%      0
%      0.3];
% end