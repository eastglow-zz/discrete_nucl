
length_scale = 1e-9;
eVpJ = 6.24150934e+18;

gamma = 0.406*eVpJ*length_scale^2;
alpha = 2.9444;

delta = 0.04;

kappa = 3*gamma*delta/alpha;
w = 3*alpha*gamma/delta;


disp(kappa);
disp(w);