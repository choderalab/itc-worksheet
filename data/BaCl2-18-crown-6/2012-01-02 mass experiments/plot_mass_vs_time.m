% Plot 18-crown-6 mass vs. time.

% Read time/mass data.
[time, mass] = textread('mass-vs-time.txt', '%f %f', 'headerlines', 4);

mass_error = 0.01 * ones(size(mass));

% Plot.
clf;
%plot(time, mass, '.', 'markersize', 10);
errorbar(time, mass, mass_error, '.', 'markersize', 10);
xlabel('time / min');
ylabel('mass / mg');
title('18-crown-6 mass at 70.1F (21.2C) on 2 Jan 2012 (Mettler-Toledo XA105 DualRange)');
axis([min(time) max(time) min(mass) max(mass)]);

% Fit exponential.
expfun = @(P, time) P(1) + P(2)*exp(-P(3)*time);
objfun = @(P) sum((expfun(P, time) - mass).^2);

P0 = [mass(end) (mass(end) - mass(1)) 1.0/time(end)];
P = fminsearch(objfun, P0);
P
hold on;
timevec = linspace(min(time), max(time), 1000);
plot(timevec, expfun(P, timevec), 'r-'); 


