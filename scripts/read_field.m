i_theta = 8;
i_chi = 1;

if i_theta <= size(Field,2)/2
    disp('WARNING: mu < 0! Setting i_theta to the first emerging direction:')    

    i_theta = size(Field,2)/2 + 1;

    disp(i_theta)

end

% load_nu_grid;
nu_grid = nu_grid(:,1);

mu  = mu_grid(i_theta);
chi = chi_grid(i_chi);

xmin = 4225;
xmax = 4228;

subplot(2,2,1)
set_plot_defaults
plot(nu_grid,Field{1,i_theta,i_chi})
xlim([xmin, xmax])
hold on
title(['$(\mu,\chi) = ($', num2str(mu,3), ', ', num2str(chi,3),')'])
ylabel('$I$')
set_plot_defaults

subplot(2,2,2)
set_plot_defaults
plot(nu_grid,Field{2,i_theta,i_chi})
xlim([xmin, xmax])
hold on
ylabel('$Q/I[\%]$')
set_plot_defaults

subplot(2,2,3)
set_plot_defaults
plot(nu_grid,Field{3,i_theta,i_chi})
xlim([xmin, xmax])
hold on
ylabel('$U/I[\%]$')
set_plot_defaults

subplot(2,2,4)
set_plot_defaults
plot(nu_grid,Field{4,i_theta,i_chi})
xlim([xmin, xmax])
hold on
ylabel('$V/I[\%]$')
set_plot_defaults


