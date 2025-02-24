close all

q1 = readmatrix("AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_V0_KQ_quad1_cont.txt");
q2 = readmatrix("AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_V0_KQ_quad2_cont.txt");

cm2 = readmatrix("cont_map_2048_8.csv");
cm1 = readmatrix("cont_map_2048_5.csv");

cm2 = cm2';
cm1 = cm1';

q2 = q2(1:64,:);
q1 = q1(1:64,:);

figure

subplot(1,2,1)
plot_rel_error(cm1, q1, "Err QUAD 1");

subplot(1,2,2)
plot_rel_error(cm2, q2, "Err QUAD 2");

figure

subplot(2,2,1)
plot_map(q1, "Hanle HQ Q1")

subplot(2,2,2)
plot_map(q2, "Hanle HQ Q2")

subplot(2,2,3)
plot_map(cm1, "TRIP Q1")

subplot(2,2,4)
plot_map(cm2, "TRIP Q2")

function ee = plot_rel_error(cm, q, tit)
    e = log10(100 * abs(cm - q) ./ abs(q));
    h2 = pcolor(e);
    set(h2, 'EdgeColor', 'none');
    colormap(gca, 'turbo');
    colorbar;
    title(tit)
end

function ee = plot_map(m, tit)
    h2 = pcolor(m);
    set(h2, 'EdgeColor', 'none');
    colormap(gca, 'turbo');
    colorbar;
    title(tit)
end


