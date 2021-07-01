

r = [0.02                 0.1              0.15              0.175             0.2              2];
Fi = [0.0093153422029699  0.23288355504762 0.52398799885721  0.71320588733332  0.93153422019047 9.315342e+01];
Ff = [0.042761595117328   0.38677393742063 0.82616794246561  1.0038255406563   1.2177427500954  115.17074773059];
vol =[1^2                 5^2              7.5^2             8.75^2            10^2             100^2];

df = (Ff - Fi)./vol;

disp(df);

scatter(r, df, 'k*', 'MarkerSize', 5);

% xlim([0 2]);


a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',20);
xlabel('Radius (nm)', 'FontSize', 24);
ylabel('\Deltag', 'FontSize', 24);