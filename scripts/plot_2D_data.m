% input path where parse.sh is 
disp('Set correct filename in parse.sh and LOS')
% disp('Reading input files...')
% 
% % Start the Bash script
% system('./parse.sh'); 
% disp('Read completed!');

% load data 
I = load("I.txt");
Q = load("Q.txt");
U = load("U.txt");
V = load("V.txt");

N = sqrt(length(I));

I = reshape(I,N,N);
Q = reshape(Q,N,N);
U = reshape(U,N,N);
V = reshape(V,N,N);

colormap(hot)

subplot(2,2,1)
set_plot_defaults
imagesc(I)
colorbar
title('$I$')

subplot(2,2,2)
set_plot_defaults
imagesc(Q)
colorbar
title('$Q/I[\%]$')

subplot(2,2,3)
set_plot_defaults
imagesc(U)
colorbar
title('$U/I[\%]$')

subplot(2,2,4)
set_plot_defaults
imagesc(V)
colorbar
title('$V/I[\%]$')

% title
load('nu_grid_96.mat')
nu = nu_grid(48);
sgtitle(['$\mu=,\chi=,\nu = $ ', num2str(nu)]);

