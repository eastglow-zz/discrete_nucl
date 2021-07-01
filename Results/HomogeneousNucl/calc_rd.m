
% computes diffusion radius
rN = 2; %nm

CN_eq = 1;
CB_eq = 0;
Ci = 0.001;


rd = rN * sqrt(1+(CN_eq-CB_eq)/(Ci-CB_eq));

disp(rd);
