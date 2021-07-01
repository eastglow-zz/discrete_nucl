clear;
close all

fname = './Em_vary_out_samples_0001.csv';
data_Em_sample = table2array(readtable(fname));
fname = './Em_vary_out_results_0001.csv';
data_Em_result = table2array(readtable(fname));
fname = './Ef_vary_out_samples_0001.csv';
data_Ef_sample = table2array(readtable(fname));
fname = './Ef_vary_out_results_0001.csv';
data_Ef_result = table2array(readtable(fname));
fname = './flux_vary_out_samples_0001.csv';
data_flux_sample = table2array(readtable(fname));
fname = './flux_vary_out_results_0001.csv';
data_flux_result = table2array(readtable(fname));
fname = './sig_vary_out_samples_0001.csv';
data_sig_sample = table2array(readtable(fname));
fname = './sig_vary_out_results_0001.csv';
data_sig_result = table2array(readtable(fname));

Em = data_Em_sample(:,1);
Ef = data_Ef_sample(:,2);
log_flux = data_flux_sample(:,3);
sig = data_sig_sample(:,4);

time_Em = data_Em_result(:,2);
maxcv_Em = data_Em_result(:,1);
time_Ef = data_Ef_result(:,2);
time_flux = data_flux_result(:,2);
time_sig = data_sig_result(:,2);


figure(1)
plot(Em, time_Em, 'k*-',  'LineWidth', 2);
hold on;
xlim([0.055 0.205]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',20);
xlabel('E_m (eV)', 'FontSize', 24);
ylabel('Time (microseconds)', 'FontSize', 24);
hold off;

figure(2)
plot(Em, maxcv_Em, 'k*-',  'LineWidth', 2);
hold on;
xlim([0.055 0.205]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',20);
xlabel('E_m (eV)', 'FontSize', 24);
ylabel('c_{v,max} (mole fraction)', 'FontSize', 24);
hold off;

figure(3)
plot(Ef, time_Ef, 'k*-',  'LineWidth', 2);
hold on;
xlim([0.1 0.4]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',20);
xlabel('E_f (eV)', 'FontSize', 24);
ylabel('Time (microseconds)', 'FontSize', 24);
hold off;

figure(4)
plot(log_flux, time_flux, 'k*-',  'LineWidth', 2);
hold on;
xlim([-5.0 1]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',20);
set(gca,'YScale', 'log');
xlabel('log_{10}(Flux) (eV)', 'FontSize', 24);
ylabel('Time (microseconds)', 'FontSize', 24);
hold off;

figure(5)
plot(log_flux, time_flux, 'k*-',  'LineWidth', 2);
hold on;
xlim([-5.0 1]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',20);
% set(gca,'YScale', 'log');
xlabel('log_{10}(Flux) (eV)', 'FontSize', 24);
ylabel('Time (microseconds)', 'FontSize', 24);
hold off;

figure(6)
plot(sig, time_sig, 'k*-',  'LineWidth', 2);
hold on;
xlim([0.1 0.4]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',20);
xlabel('Sigma (J/m^2)', 'FontSize', 24);
ylabel('Time (microseconds)', 'FontSize', 24);
hold off;

figure(7)
plot(sig, time_sig, 'k*-',  'LineWidth', 2);
hold on;
xlim([0.1 0.4]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',20);
set(gca,'YScale', 'log');
xlabel('Sigma (J/m^2)', 'FontSize', 24);
ylabel('Time (microseconds)', 'FontSize', 24);
hold off;

