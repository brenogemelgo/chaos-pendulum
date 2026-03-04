clc; clearvars; close all

%-------------------------------------------------------------------------%

m1 = 1;
m2 = 1;

l1 = 3;
l2 = 2;

th1 = pi;
th2 = pi;

w1 = 0;
w2 = 0;

while true
    vt = str2double(input('Enter the simulation time (in seconds): ', 's'));

    if isfinite(vt) && (vt > 0)
        break
    end

    disp('Error: Input value must be a finite positive number!')
    disp(' ')
end

g = 9.81;
m = [m1 m2];
l = [l1 l2];
y0 = [th1 th2 w1 w2];

fprintf('Simulation Time: %.2f seconds \n', vt)

%-------------------------------------------------------------------------%

sim_fps = 400;
render_fps = 60;
speed = 1.0; % 1.0 real-time, 0.5 half-speed, 2.0 double-speed

nSim = max(4000, round(vt * sim_fps) + 1);
tSim = linspace(0, vt, nSim);

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

disp('Simulating...')
[tSol, ySol] = ode45(@(t, y) pendulum(t, y, m, l, g), tSim, y0, opts);
disp('Rendering...')

nR = max(2, round(vt * render_fps) + 1);
tR = linspace(0, vt, nR);

yR = interp1(tSol, ySol, tR, 'linear');

th1v = yR(:, 1);
th2v = yR(:, 2);

s1 = sin(th1v); c1 = cos(th1v);
s2 = sin(th2v); c2 = cos(th2v);

x1 = l(1) * s1;
y1p = -l(1) * c1;

x2 = x1 + l(2) * s2;
y2p = y1p - l(2) * c2;

%-------------------------------------------------------------------------%

screen = get(0, 'ScreenSize');
figure_width = 1200;
figure_height = 600;
figure_left = (screen(3) - figure_width) / 2;
figure_bottom = (screen(4) - figure_height) / 2;

L = sum(l);
axlim = 1.2 * L;

h = figure('Color', 'k', 'Name', 'Pendulums Simulation', 'NumberTitle', 'off');
set(h, 'Position', [figure_left figure_bottom figure_width figure_height]);

ax = axes('Parent', h);
axis(ax, 'equal');
axis(ax, [-axlim axlim -axlim axlim]);
set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
hold(ax, 'on');

rod = plot(ax, [0 x1(1) x2(1)], [0 y1p(1) y2p(1)], 'o-', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'w');

traj1 = plot(ax, nan, nan, 'm', 'LineWidth', 1);
traj2 = plot(ax, nan, nan, 'r', 'LineWidth', 1);

title(ax, sprintf('Time: %.2f s', tR(1)), 'Color', 'w');

traj_x = nan(2, nR);
traj_y = nan(2, nR);

dt = tR(2) - tR(1);

tic

for i = 1:nR

    if ~ishandle(h)
        break
    end

    traj_x(:, i) = [x1(i); x2(i)];
    traj_y(:, i) = [y1p(i); y2p(i)];

    set(rod, 'XData', [0 x1(i) x2(i)], 'YData', [0 y1p(i) y2p(i)]);
    set(traj1, 'XData', traj_x(1, 1:i), 'YData', traj_y(1, 1:i));
    set(traj2, 'XData', traj_x(2, 1:i), 'YData', traj_y(2, 1:i));

    title(ax, sprintf('Time: %.2f s', tR(i)), 'Color', 'w');
    drawnow limitrate nocallbacks

    target = (i - 1) * (dt / speed);

    while toc < target
        pause(min(0.002, target - toc));
    end

end

%-------------------------------------------------------------------------%

disp('Plotting full trajectory...')
figure('Color', 'k', 'Name', 'Full Trajectory', 'NumberTitle', 'off');
plot(x1, y1p, 'm', 'LineWidth', 1); hold on
plot(x2, y2p, 'r', 'LineWidth', 1); hold off
axis equal
axis([-axlim axlim -axlim axlim])
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
title({sprintf('Time: %.2f s', vt), 'Full Trajectory'}, 'Color', 'w')
legend({'Pendulum 1 Trajectory', 'Pendulum 2 Trajectory'}, 'TextColor', 'white', 'Location', 'southoutside')

%-------------------------------------------------------------------------%

function dydt = pendulum(~, y, m, l, g)
    th1 = y(1);
    th2 = y(2);
    w1 = y(3);
    w2 = y(4);

    delta = th2 - th1;
    s = sin(delta);
    c = cos(delta);

    den1 = (m(1) + m(2)) * l(1) - m(2) * l(1) * c * c;
    den2 = (l(2) / l(1)) * den1;

    dth1dt = w1;
    dth2dt = w2;

    dw1dt = (m(2) * l(1) * w1 * w1 * s * c ...
        + m(2) * g * sin(th2) * c ...
        + m(2) * l(2) * w2 * w2 * s ...
        - (m(1) + m(2)) * g * sin(th1)) / den1;

    dw2dt = (-m(2) * l(2) * w2 * w2 * s * c ...
        + (m(1) + m(2)) * (g * sin(th1) * c - l(1) * w1 * w1 * s - g * sin(th2))) / den2;

    dydt = [dth1dt; dth2dt; dw1dt; dw2dt];
end
