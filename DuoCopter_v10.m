clear all
close all
clc


%% CONFIGURATION

% Geometry
config.m = 1; % mass
config.g = 9.8; % gravitational constant
config.L = 1; % length
config.I = 1/12*config.m*config.L^2; % moment of intertia

% Controler
config.Ncontrollers = 50; % # of duocopters
config.lqr_use_nonlinear = false; % toggle to re-calculate LQR K at different thetas (better performance)
config.lqr_theta_vec = linspace(-pi,pi, 31);  % number of LQR interpolation points non nonlinear LQR

% Simulation
config.dt = .01; % simulation time step
config.Nsteps = 5000; % number of steps
config.t = 0:config.dt:config.dt*(config.Nsteps-1); % time vector (do not change)
config.xLims = [-10 10]; % x boundaries
config.yLims = [-10 10]; % y boundaries
config.targetall = make_target(config.t); % [Nsteps x 6] target vector repretsenting desired {x, y, theta, vx, vy, w} for every step

% Plotting and Animation
config.play_animation = true; % toggle to play animation
config.frametime = .01; % frametime of animation (does not affect simulation)
config.play_every_n_frames = 2; % play 1 of every n frames (to speed up animation)
config.mycolors = lines(config.Ncontrollers); % [Ncontrollers x 3] rgb for desired colors for every duocopter 




%% SIMULATION

[config, XU] = set_up_controllers(config);
[XU, config] = run_simulation(config, XU);
make_static_plot(XU, config);
plot_duocopter(XU, config)













%% FUNCTIONS





function [XU, config] = run_simulation(config, XU)

for istep = 2:config.Nsteps
    config.target = config.targetall(istep,:);
    
    for ic = 1:config.Ncontrollers
        XU(istep-1, 7:8, ic) = controller(squeeze(XU(istep-1,:,ic)), config, ic);
    end
    
    XU(istep,:,:) = duocopter_dynamics(squeeze(XU(istep-1,:,:)), config);
end

disp('done running simulation');

end 



function [config, XU] = set_up_controllers(config)
rng(0);
X0 = zeros(8,config.Ncontrollers);
XU = zeros(config.Nsteps,8, config.Ncontrollers);
XU(1, :,:) = X0;

for ic = 1:config.Ncontrollers
    Q = eye(6).*([1 1 1 10 10 100] + 1000*rand(1,6));
    R = eye(2).*rand().*[100 100];
    
    config.K{ic} = calc_k_LQR([0 0 0 0 0 0], config, Q, R);
    
    for itheta = length(config.lqr_theta_vec):-1:1
        KNL(itheta,:, :) = calc_k_LQR([0 0 config.lqr_theta_vec(itheta) 0 0 0], config, Q, R);     
    end
    config.KNL{ic} = KNL;
    
end
end 


function T = make_target(t)


y = .5*[10 10 0 0 0 0
    -10 10 0 0 0 0
    -10 -10 0 0 0 0
    10 -10 0 0 0 0
    10 10 0 0 0 0];

x = linspace(0,t(end),size(y,1));
T = interp1(x', y, mod(2*t, t(end)), 'previous');






end



function hfig = make_static_plot(XU, config)

mycolors = config.mycolors;
colormap(mycolors)

hfig = figure();

hax = subplot(3,2,2);
hold on
h1 = plot(config.t, squeeze(XU(:,1,:)));
hax.ColorOrderIndex = 1;
h2 = plot(config.t, squeeze(XU(:,2,:)), '--');
xlabel('t'); ylabel('position');
legend([h1(1) h2(1)], {'x', 'y'})
hold off


hax = subplot(3,2,4);
hold on
h1 = plot(config.t,squeeze(XU(:,7,:))/(config.m*config.g));
hax.ColorOrderIndex = 1;
h2 = plot(config.t,squeeze(XU(:,8,:))/(config.m*config.g), '--');
xlabel('t');
ylabel('force/mg')
legend([h1(1) h2(1)],{'F_1', 'F_2'})
hold off

subplot(3,2,6);
hold on
plot(config.t,squeeze(XU(:,3,:))*180/pi)
ylabel('\theta (deg)')
xlabel('t');
hold off

subplot(3,2,[1 3 5]);
hold on
plot(squeeze(XU(:,1,:)), squeeze(XU(:,2,:)))
plot(squeeze(XU(end,1,:)), squeeze(XU(end,2,:)), 'kx')
title('y vs. x')
hold off
axis equal

drawnow

end



function K = calc_k_LQR(XU, config, Q, R)


%XU(1,:) = x
%XU(2,:) = y
%XU(3,:) = theta
%XU(4,:) = vx
%XU(5,:) = vy
%XU(6,:) = w
%XU(7,:) = F1
%XU(8,:) = F2

F1 = config.m*config.g/2/cos(XU(3));
F2 = config.m*config.g/2/cos(XU(3));
% F1 = config.m*config.g/2;
% F2 = config.m*config.g/2;

A = [0 0 0 1 0 0
    0 0 0 0 1 0
    0 0 0 0 0 1
    0 0 -cos(XU(3))*(F1+F2)/config.m 0 0 0
    0 0 -sin(XU(3))*(F1+F2)/config.m 0 0 0
    0 0 0 0 0 0];

B = [0 0
    0 0
    0 0
    -sin(XU(3))/config.m -sin(XU(3))/config.m
    cos(XU(3))/config.m cos(XU(3))/config.m
    -config.L/(2*config.I) config.L/(2*config.I)];

%Q = eye(6).*[1 1 1 1 1 10];
%R = eye(2).*[.1 .1];

[K,~,~] = lqr(A, B, Q, R);
end


%%

function FF = controller(XU, config, ic)


if config.lqr_use_nonlinear
    
    
    I = find(config.lqr_theta_vec>= XU(3), 1);
    
    if isempty(I)
        K = squeeze(interp1(config.lqr_theta_vec, config.KNL{ic}, XU(3), 'nearest', 'extrap'));
    else
        tmp = config.KNL{ic}(I, :, :);
        K(:, :) = tmp(1, :, :);
        %K = squeeze([config.KNL{ic}(I, :, :)]);
    end
    
else
    K = config.K{ic};
end

X = XU(1:6);

FF = -K*(X-config.target)';

%FF = FF + config.m*config.g/2/cos(XU(3));
FF = FF + config.m*config.g/2;

FF(FF>config.m*config.g) = config.m*config.g;
FF(FF<-config.m*config.g) = -config.m*config.g;

end



function X2 = duocopter_dynamics(X, config)

%X(1,:) = x
%X(2,:) = y
%X(3,:) = theta
%X(4,:) = vx
%X(5,:) = vy
%X(6,:) = w
%X(7,:) = F1
%X(8,:) = F2

Fx = -sin(X(3,:)).*(X(7,:)+ X(8,:));
Fy = cos(X(3,:)).*(X(7,:)+ X(8,:)) - config.m*config.g;
M = config.L/2*(X(8,:)-X(7,:));

ddx = Fx./config.m;
ddy = Fy./config.m;
ddtheta = M./config.I;

X2 = X;
X2(1,:)  = X(1,:) + X(4,:)*config.dt + .5*ddx*config.dt^2;
X2(2,:)  = X(2,:) + X(5,:)*config.dt + .5*ddy*config.dt^2;
X2(3,:)  = X(3,:) + X(6,:)*config.dt + .5*ddtheta*config.dt^2;
X2(4,:)  = X(4,:) + ddx*config.dt;
X2(5,:)  = X(5,:) + ddy*config.dt;
X2(6,:)  = X(6,:) + ddtheta*config.dt;
X2(7,:)  = X(7,:);
X2(8,:)  = X(8,:);

X2(1, X2(1,:)>config.xLims(2)) = config.xLims(2);
X2(1, X2(1,:)<config.xLims(1)) = config.xLims(1);
X2(2, X2(2,:)>config.yLims(2)) = config.yLims(2);
X2(2, X2(2,:)<config.yLims(1)) = config.yLims(1);
% X2(3, X2(3,:)>pi())  = X2(3, X2(3,:)>pi()) - 2*pi();
% X2(3, X2(3,:)<-pi())  = X2(3, X2(3,:)<-pi()) + 2*pi();
end




function plot_duocopter(XU, config)

if config.play_animation
    
    Ncontrollers = config.Ncontrollers;
    handles.hfig = figure();
    handles.hax = axes();
    axis equal
    handles.hax.XLim=config.xLims;
    handles.hax.YLim=config.yLims;
    handles.mycolors = config.mycolors;
    
    hold on
    
    geom.pX = [-1 -1;    1   1;    -1 1;    -1.4 -.6;    .6 1.4;    0 0];
    geom.pY = [0 -1;    0 -1;    0 0;    .2 .2;    .2 .2;    0 0];
    
    handles.ship = cell(1,Ncontrollers);
    for ic = 1:Ncontrollers
        hl(1) = plot([nan nan], [nan nan], 'Color', [.7 .7 .7], 'Parent', handles.hax, 'LineWidth', 4);
        hl(2) = plot([nan nan], [nan nan], 'Color', [.7 .7 .7], 'Parent', handles.hax, 'LineWidth', 4);
        hl(3) = plot([nan nan], [nan nan], 'Color', handles.mycolors(ic,:), 'Parent', handles.hax, 'LineWidth', 2);
        hl(4) = plot([nan nan], [nan nan], 'Color',handles.mycolors(ic,:), 'Parent', handles.hax, 'LineWidth', 2);
        hl(5) = plot([nan nan], [nan nan], 'Color',handles.mycolors(ic,:), 'Parent', handles.hax, 'LineWidth', 2);
        hl(6) = plot([nan nan], [nan nan], 'Color',handles.mycolors(ic,:), 'Marker', 'd', 'Parent', handles.hax, 'LineWidth', 2, 'MarkerFaceColor', handles.mycolors(ic,:));
        handles.ship{ic} = hl;
    end
    
    hold off
    
    for it = 1:config.play_every_n_frames:size(XU, 1)
        tic
        handles.hfig.Visible = 0;
        for ic = 1:Ncontrollers
            x = XU(it,1, ic);
            y = XU(it,2, ic);
            theta =XU(it,3, ic);
            F1 = XU(it,7, ic);
            F2 = XU(it,8, ic);
            
            pY_ad = geom.pY;
            pY_ad(1,2) = -F1/(config.m*config.g);
            pY_ad(2,2) = -F2/(config.m*config.g);
            
            xr = geom.pX*cos(theta) - pY_ad*sin(theta) + x;
            yr = geom.pX*sin(theta) + pY_ad*cos(theta) + y;
            xvalues = num2cell(xr, 2);
            yvalues = num2cell(yr, 2);
            [handles.ship{ic}.XData] =  xvalues{:};
            [handles.ship{ic}.YData] = yvalues{:};
            
        end
        
        dt_frame = toc;
        if dt_frame<config.frametime
            pause(config.frametime-dt_frame);
        end
        handles.hfig.Visible = 1;
        drawnow
        
    end
end

end



