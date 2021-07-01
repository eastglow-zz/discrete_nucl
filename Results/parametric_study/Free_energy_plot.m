clear;
close all;


Va = 0.02024; %atoms per nm^3
kb = 8.6173303e-5; 
T = 300; %K
length_scale = 1e-9;
time_scale = 1e-6;

cv = 0.0:0.0001:1.0;

Ef = 0.52; %eV
fbulk = Ef/Va*cv + kb*T/Va*(cv.*log(cv) + (1-cv).*log(1-cv));


fbulk_200 = 0.200/Va*cv + kb*T/Va*(cv.*log(cv) + (1-cv).*log(1-cv));
fbulk_300 = 0.300/Va*cv + kb*T/Va*(cv.*log(cv) + (1-cv).*log(1-cv));
fbulk_400 = 0.400/Va*cv + kb*T/Va*(cv.*log(cv) + (1-cv).*log(1-cv));
fbulk_500 = 0.500/Va*cv + kb*T/Va*(cv.*log(cv) + (1-cv).*log(1-cv));
fbulk_600 = 0.600/Va*cv + kb*T/Va*(cv.*log(cv) + (1-cv).*log(1-cv));
fbulk_700 = 0.700/Va*cv + kb*T/Va*(cv.*log(cv) + (1-cv).*log(1-cv));
fbulk_800 = 0.800/Va*cv + kb*T/Va*(cv.*log(cv) + (1-cv).*log(1-cv));



% tol = 1e-3;
% plog = log(tol) + (cv - tol)/tol - ((cv - tol).^2)/(2*tol*tol) + ((cv-tol).^3)/(3*tol*tol*tol);
% fbulk_plog = Ef/Va*cv + kb*T/Va*(cv.*plog + (1-cv).*plog);

plot(cv,fbulk_200,'linewidth',2)
hold on;
plot(cv,fbulk_300, 'linewidth',2);
plot(cv,fbulk_400, 'linewidth',2);
plot(cv,fbulk_500, 'linewidth',2);
plot(cv,fbulk_600, 'linewidth',2);
plot(cv,fbulk_700, 'linewidth',2);
plot(cv,fbulk_800, 'linewidth',2);
plot(cv,fbulk, '-k','linewidth',3);
ylim([-1,60]);
% xlim([0,1]);
set(gca,'fontsize',18)
xlabel('c_v (molar fraction)')
ylabel('Free energy (eV/nm^3)')
% legend('f_B','f_{d}','f_{V}','location','north')
legend('E_f = 0.2 eV','E_f = 0.3 eV','E_f = 0.4 eV','E_f = 0.5 eV','E_f = 0.6 eV','E_f = 0.7 eV','E_f = 0.8 eV','E_f = 0.52 eV (Original)','location','northwest')
legend boxoff



