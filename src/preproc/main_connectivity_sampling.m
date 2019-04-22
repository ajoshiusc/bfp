
clc;clear all;close all;

load /deneb_disk/SCD_BOLDdata/atlas.mat
full_graph=atlas_data';
sample_rate=.1;
para_sigma=1;
nSurf=32492;
v_conn_ext_list='conn';
tic
full_graph_left=full_graph(1:nSurf,:);
[ sample_graph_left, sample_index_left ] = connectivity_sampling(full_graph_left, sample_rate, para_sigma );
toc

full_graph_right=full_graph(1+nSurf:2*nSurf,:);
[ sample_graph_right, sample_index_right ] = connectivity_sampling(full_graph_right, sample_rate, para_sigma );
toc
full_graph_subc=full_graph(2*nSurf+1:end,:);
[ sample_graph_subc, sample_index_subc ] = connectivity_sampling( full_graph_subc, sample_rate, para_sigma );
toc


