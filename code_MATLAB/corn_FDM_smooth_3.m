 clear; clf;

Nt = 100; % # of time steps
T = 1/365; % Near the maturity
dt = T/Nt; % For around the maturity
E = 100; % strike price
L = 300; % sufficiently large value(computational domain)
sig = 0.5; % volatility
r = 0.03; % riskless interest rate
dx = 0.1; % spatial step

% grid construct - nonuniform
x=[0 0.5*dx:dx:L-0.5*dx];
Nx = length(x);
h = diff(x);
h = [h, h(end)];

% digital option(cach-or-nothing)
cash = 100;
u = zeros(Nx, Nt+1);

%%%%%%% version 2 %%%%%%%
equd = 0.5;
rev = (x-E)./(equd);

y3 = zeros(1, Nx);
w = 0.7; % weight
for i = 1:Nx
    if  rev(i) >= -2 && rev(i) <= -1
        y3(i) = (1 - w) * ( 0.25*(2 + rev(i))^4 - (1/3)*(2 + rev(i))^3 );
    elseif rev(i) >= -1 && rev(i) <= 0
        y3(i) = (1 - w) * ( 0.25*(1+rev(i))^4 + (1/3)*(1+rev(i))^3 + 0.5*(1+rev(i))^2 - (1/12) )...
            + w*(0.5*(1+rev(i))^2) ;
    elseif rev(i) >= 0 && rev(i) <= 1
        y3(i) = 1 - ((1 - w) * ( 0.25*(1-rev(i))^4 + (1/3)*(1-rev(i))^3 + 0.5*(1-rev(i))^2 - (1/12) )...
            + w*(0.5*(1-rev(i))^2));        
    elseif rev(i) >= 1 && rev(i) <= 2
        y3(i) = 1 - (1 - w) * ( 0.25*(2-rev(i))^4 - (1/3)*(2-rev(i))^3);
    elseif rev(i) < -2
        y3(i) = 0;
    else
        y3(i) = 1;
    end
end

u(:, 1) = cash*y3; % payoff
%%%%%%%%%%%%%%%%%%%%

% tridiagonal matrix
d = zeros(1,Nx-1); c = zeros(1,Nx-1); a = zeros(1,Nx-1);
for i = 1:Nx-1
    d(i) = 1 + dt*( ((sig*x(i+1))^2 - r*x(i+1)*(h(i+1)-h(i))) / (h(i)*h(i+1)) + r );
    c(i) = dt*( -(sig*x(i+1))^2 - r*x(i+1)*h(i)) / (h(i+1)*(h(i)+h(i+1)) );
    a(i) = dt*( -(sig*x(i+1))^2 + r*x(i+1)*h(i+1)) / (h(i)*(h(i)+h(i+1)) );
end
% linear boundary condition
d(Nx-1) = d(Nx-1) + 2*c(Nx-1);
a(Nx-1) = a(Nx-1) - c(Nx-1);

% time loop
for n = 1:Nt
    b = u(2:Nx,n);
    u(2:Nx,n+1) = thomas(a, d, c, b);
end


% exact option price
d1 = (log(x/E) + (r + 0.5*sig^2)*T) / (sig*sqrt(T));
d2 = d1 -(sig*sqrt(T));
exc = cash * exp(-r * T)*normcdf(d2);

% error
bc = 0.8; ec = 1.2;
bidx = find(x > bc*E, 1 );
eidx = find(x < ec*E, 1, 'last' );

RMSE = sqrt(mean( (u(bidx:eidx,end)-exc(bidx:eidx)').^2 ) );
maxerr = max(abs(u(:,end) - exc'));
fprintf('RMSE : %.8f\n', RMSE);
fprintf('maxerr : %.8f\n\n', maxerr);

% plot
FigHandle = figure(1);
set(FigHandle, 'Position', [100, 100, 750, 750]);
set(gca,'fontsize',20);
xlabel('x'); ylabel('u(x,t)');
axis([90 110 -5 105])
hold on;
grid on;
plot(x,u(:,1)','k-')
plot(x, exc,'r*-')
plot(x,u(:,end)','ko-')

yyaxis right
plot(x, abs(u(:,end) - exc'), 'b--' )

title('corn\_FDM\_smooth\_3.m')
legend('Payoff', 'Exact', 'FDM solution', 'error', 'location','northwest')