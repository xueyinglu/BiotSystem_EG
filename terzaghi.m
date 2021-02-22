% Terzaghi analytical solution
k = 1e-3TT;
E = 1e5;
nu = 0.2;
alpha = 1;
inv_M = 0.1;
lambda = E*nu./(1+nu)./(1-2*nu);
mu = E./2./(1+nu);
K_bulk = lambda + 2./3*mu;
K_u = K_bulk + alpha*alpha*inv_M;
c_f = inv_M * k *(K_bulk +4./3*mu)./ (K_u + 4./3 *mu);

F =1000;

t_interval = [0.1:0.1:1];
[y,t] = meshgrid(0:0.01/256:0.01,0.1:0.1:1);

pressure = zeros(size(y));
s = zeros(size(t_interval));
n_terms = 100000;
alpha * F * inv_M ./(K_u + 4./3 *mu)

for i=0:n_terms
    M = pi *(2*i+1)/2;
    pressure = pressure+ alpha * F * inv_M /(K_u + 4./3 *mu) * 2/M .* sin (M .*y) .* exp(-M*M*c_f.*t);
    %pressure = pressure+ 2/M .* sin (M .*y) .* exp(-M*M*c_f.*t);
    s = s+ 2/M/M .* exp(-M*M*c_f.*t_interval);
end
figure
surf(y,t,pressure)
xlabel("y")
ylabel("t")
zlabel("p")
figure
plot(t_interval,s)