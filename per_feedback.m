clc;
% clear;
close all;
%% Initialization
xf = 1.5;
x0 = 1;
tf = 25;
t0 = 0;


%% iLQR parameters
N = 10000;
dt = (tf - t0)./N;
alpha = 0.5;


Q = 2;
R = 2;
x_nom = zeros(1,N+1);
x_nom(:,1) = x0;
x_new = zeros(1,N+1);
x_new(:,1) = x0;
u_nom = randn(1,N);
dx = zeros(1,N+1);
du = zeros(1,N);
S = zeros(1,N+1);
S(:,N+1) = 0;
P = zeros(1,N+1);
P(:,N+1) = 0;
V = zeros(1,N+1);
V(:,N+1) = 1;
U = zeros(1,N+1);
U(:,N+1) = 1;
w = zeros(1,N+1);
w(:,N+1) = 0;
K = zeros(1,N);
K_u = zeros(1,N);
k = zeros(1,N);

eps = 0;



criteria = true;
iter = 0;
max_iter = 10000;

while iter < max_iter

    iter = iter + 1;

    for i = 1:1:N
        x_nom(:,i+1) = x_nom(:,i) + dt*(-x_nom(:,i)^3 + u_nom(:,i));
    end

    d_psi = xf - x_nom(:,N+1);
    % d_psi = dx(:,N+1);
    for i = N:-1:1
        A = (1 + dt*(-3*(x_nom(:,i)^2)));
        B = dt;
        Quu = 2 + B'*S(:,i+1)*B;
        if i == N
            w(:,i+1) = -inv(B)*2*u_nom(:,N);
        end

        S(:,i) = A'*S(:,i+1)*A + Q - A'*S(:,i+1)*B*inv(Quu)*B'*S(:,i+1)*A;
        V(:,i) = A'*V(:,i+1) - A'*S(:,i+1)*B*inv(Quu)*B'*V(:,i+1);
        w(:,i) = Q*x_nom(:,i) + A'*w(:,i+1) - A'*S(:,i+1)*B*...
            inv(Quu)*(B*w(:,i+1) + R*u_nom(:,i));
        U(:,i) = U(:,i+1)*(A - B*inv(Quu)*B'*S(:,i+1)*A);
        P(:,i) = P(:,i+1) + U(:,i+1)*B*(-inv(Quu)*B'*V(:,i+1));
        K(:,i) = inv(Quu)*B'*V(:,i+1)*inv(P(:,i))*U(:,i) - ...
            inv(Quu)*B'*S(:,i+1)*A;
        K_u(:,i) = -inv(Quu)*B'*V(:,i+1)*inv(P(:,i));
        k(:,i)   = -inv(Quu)*2*u_nom(:,i) -inv(Quu)*B'*w(:,i+1);
        % K_u(:,i) = dt;
        
    end

    for i = 1:1:N
        du(:,i) = K(:,i)*dx(:,i) + K_u(:,i)*d_psi + alpha*k(:,i);
        u_nom(:,i) = u_nom(:,i) + du(:,i);
        x_new(:,i+1) = x_new(:,i) + dt*(-x_new(:,i)^3 + u_nom(:,i));
        dx(:,i+1) = x_new(:,i+1) - x_nom(:,i+1);
    end
    
    % x_nom = x_new;
    
    if abs(d_psi)<1e-10
        break;
    end
    terminal_error(iter) = d_psi;
    state_error(iter) = norm(dx,inf);
    % iter
end

% u_nom = [u_nom(:,1), u_nom];
iter
% J = (sum(u_nom.^2) + sum(x_nom.^2))*dt;
lambda = zeros(1,N+1);
mu = -inv(P(:,1))*(U(:,1)*dx(:,1) - d_psi);
for i = N+1:-1:1
    lambda(i) = S(:,i)*dx(:,i) + V(:,i)*mu + w(:,i);
    lambda(i) = lambda(i)*dt;
end
%% Plots

% figure(1)
% plot(0:dt:tf-dt,K)
% xlabel('Time')
% ylabel('State feedback')
% grid on
% hold on
% 
% figure(2)
% plot(0:dt:tf-dt,K_u)
% xlabel('Time')
% ylabel('Term. offset feedback')
% grid on
% hold on
% 
figure(3)
plot(0:dt:tf,x_nom)
xlabel('Time')
ylabel('State(x)')
grid on
hold on

figure(4)
plot(0:dt:tf,lambda)
xlabel('Time')
ylabel('\lambda')
grid on
% hold on
% plot(0:dt:tf,lambda)
% plot(0:dt:tf,w.*dt)
%%
% sig = 0.6;
% for k = 1:1:100
%     x = zeros(1,N+1);
%     err = normrnd(0,sig,[1,1]);
%     x(:,1) = x0 + err;
%     noise = zeros(size(x));
%     noise(:,1) = err;
%     for i = 1:1:N
%         % err = randn(1,1);
%         x_b(:,i) = x(:,i);
%         for j = i:1:N
%             x_b(:,j+1) = x_b(:,j) + dt*(-x_b(:,j)^3 + u_nom(:,j));
%         end
%         d_psi_n = xf - x_b(:,N+1);
%         u_new(:,i) = u_nom(:,i) + 1*K(:,i)*(x(:,i) - x_nom(:,i)) + ...
%             1*K_u(:,i)*d_psi_n;
%         x(:,i+1) = x(:,i) + dt*(-x(:,i)^3 + u_new(:,i));
%         if i<(N-10)
%                 err = normrnd(0,sig,[1,1]);
%                 x(:,i+1) = x(:,i+1) + 1*err;
%                 noise(:,i+1) = err;
%         end
%     end
% 
% 
%     error = x-x_nom;
%     final_error(k) = error(end);
%     figure(1)
%     plot(0:dt:tf,noise)
%     hold on
% end
% 
% 
% figure(1)
% plot(0:dt:tf,3*sig*ones(1,N+1),'--r')
% plot(0:dt:tf,-3*sig*ones(1,N+1),'--r')
% 
% figure(2)
% plot(1:1:100,(final_error./xf)*100,'o')
% ylabel('Final Error as % of final state')
% xlabel('Iterations')
% 
% x10 = x;
% u10 = u_new;
% figure
% plot(0:dt:tf,x_nom)
% hold on
% plot(0:dt:tf,x)
% legend('Opt. Sol.','Feedback sol.')
% 
% figure
% plot(0:dt:tf-dt,u_nom)
% hold on
% plot(0:dt:tf-dt,u_new)
% legend('Opt. Sol.','Feedback sol.')
% 
% x_opt = x_nom;
% u_opt = u_nom;

    
% figure
% plot(0:dt:tf,P)
% xlabel('Time')
% ylabel('P')
% grid on
% 
% figure
% plot(0:dt:tf,V)
% xlabel('Time')
% ylabel('V')
% grid on


%% Continuous-Time results


% Forward marching

% x_init = [x0;lambda(1)];
% traj = x_init;
% sim_diff_fac = 100;
% N1 = N*sim_diff_fac;
% dt1 = dt/sim_diff_fac;
% 
% for i = 1:1:N1
%     [t,x] = ode113(@(t,x) odefun(t,x), [(i-1)*dt1, i*dt1], x_init);
%     x_init = x(end,:)';
%     traj = [traj, x_init];
% end
% 
% frac = 1;
% t_plot = tf*frac;
% n = N*frac+ 1;
% n1 = N1*frac + 1;
% 
% figure(1)
% plot(0:dt1:t_plot,traj(1,1:n1));
% hold on
% plot(0:dt:t_plot,x_nom(:,1:n));
% 
% figure(2)
% plot(0:dt1:t_plot,traj(2,1:n1));
% hold on
% plot(0:dt:t_plot,lambda(:,1:n));

% Backward Marching
% x_init = [x_nom(1,end);lambda(1,end)];
% traj = x_init;
% sim_diff_fac = 100;
% N1 = N*sim_diff_fac;
% dt1 = dt/sim_diff_fac;
% 
% for i = N1:-1:1
%     [t,x] = ode113(@(t,x) odefun(t,x), [(i)*dt1, (i-1)*dt1], x_init);
%     x_init = x(end,:)';
%     traj = [x_init, traj];
% end
% 
% 
% frac = 1/1;
% t_plot = tf - tf*frac;
% n = N - N*frac + 1;
% n1 = N1 - N1*frac + 1;
% 
% figure(1)
% plot(t_plot:dt1:tf,traj(1,n1:end));
% hold on
% plot(t_plot:dt:tf,x_nom(:,n:end));
% 
% figure(2)
% plot(t_plot:dt1:tf,traj(2,n1:end));
% hold on
% plot(t_plot:dt:tf,lambda(:,n:end));
% 
% function xout = odefun(t,y)
% 
%  x = y(1);
%  lambda = y(2);
% 
%  xout = zeros(2,1);
%  xout(1) = -x.^3 - lambda./2;
%  xout(2) = -2*x + 3*lambda.*(x.^2);
% end
