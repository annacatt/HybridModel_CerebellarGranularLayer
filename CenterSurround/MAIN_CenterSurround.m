% Please cite **** TO BE ADDED ****

clear all
close all

grid_res = 65; % Grid resolution of GrC population, namely number of nodes per edge

disp(['Using GrC grid step of ' mat2str(grid_res)])
Hybrid_CenterSurround(grid_res,1) % Control
Hybrid_CenterSurround(grid_res,0.03) % Inhibition blocked
[E,I]=MakeFigure_EI_Fig4_func(grid_res);
eval(['save mat_EI_',num2str(grid_res),'.mat E I'])



