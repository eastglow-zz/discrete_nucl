

A_b = 270.5;
B_b = 32;
eps = 0.01;

A_v = 1000;
ceq_v = 1;

cp = 0.0001;

fc_b_cp = A_b * cp /sqrt(cp^2 + eps^2) + B_b*cp;

cp_v = fc_b_cp/A_v + ceq_v;

disp(cp_v);
