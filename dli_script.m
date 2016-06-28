% This script runs various DLI usages

%% run and save - 1D chain, default values
dli = DLI_1D_CHAIN();
dli = dli.run_and_plot(); %#ok<NASGU>
save(['Results',filesep,'1D_chain_default'],'dli')

%% run and save - 1D chain, default values
dli = DLI_1D_CHAIN('BC','sharp','n_cells',2);
dli = dli.run_and_plot(); %#ok<NASGU>
save(['Results',filesep,'1D_chain_2cells'],'dli')

%% load and plot - 1D chain, default values
ws = load(['Results',filesep,'1D_chain_default'],'dli');
dli = ws.dli;
clear ws
dli.plot_all()

%% run without plotting, printing message or saving - 1D chain, default values
dli = DLI_1D_CHAIN();
dli = dli.run_and_plot('skip_plot','skip_print'); %#ok<NASGU>

%% run and save - 2D hexagonal lattice, 12 cells, l_max=6
clc
dli = DLI_2D_HEX('solver',@ode45,'v_cells',3,'h_cells',4,'l_max',6,'t_max',2);
dli = dli.run_and_plot(); %#ok<NASGU>
save(['Results',filesep,'2D_hex_cells12_lmax6'],'dli')

%% run and save - 2D hexagonal lattice, 12 cells, l_max=18, no SS
close all,clc
dli = DLI_2D_HEX('solver',@ode45,'v_cells',3,'h_cells',4,'l_max',18,'t_max',100);
dli = dli.run_and_plot('skip_ss','skip_plot'); %#ok<NASGU>
save(['Results',filesep,'2D_hex_cells12_lmax6'],'dli')

%% run without save - 2D hexagonal lattice, 12 cells, l_max=42, no SS
% close all,clc
dli = DLI_2D_HEX('solver',@ode45,'v_cells',3,'h_cells',4,'l_max',6*7,'t_max',150,'IC','one_higher');
dli = dli.run_and_plot();%'skip_ss','skip_plot'); %#ok<NASGU>
% dli.plot_distributions();
save(['Results',filesep,'2D_hex_cells12_lmax42'],'dli')