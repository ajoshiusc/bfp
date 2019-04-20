function [ sample_graph, sample_index ] = IPCS( full_graph, sample_rate, para_sigma, v_conn_ext_list )
%TEST_MAX_SAMPLE Summary of this function goes here
%   Detailed explanation goes here
%%
[size_1, size_2] = size(full_graph) ;
if size_1 ~= size_2
    T = size_2;
    full_graph = (full_graph * full_graph') ./ T ;
end
rho_threshold = sqrt(1 / (size_2 - 1)) ;
[V_full, ~] = size(full_graph) ;
if sample_rate <= 1
    sample_rate = sample_rate * V_full ;
end
V_d = round(sample_rate) ;

mask_graph = ones(V_full) - eye(V_full);
if exist('v_conn_ext_list', 'var')
    for current_node = 1 : V_full
        current_list = v_conn_ext_list{current_node} ;
        mask_graph(current_node, current_list) = 0 ;
        mask_graph(current_list, current_node) = 0 ;
    end
end

%%
sample_index = [] ;
tmp_graph = abs(full_graph) ;
for current_sample = 1 : V_d
    mask_graph(abs(tmp_graph) < rho_threshold) = 0 ;
    tmp_graph = exp((tmp_graph - 1) ./ para_sigma) ;
    tmp_graph = tmp_graph .* mask_graph ;
    [~, max_node] = max(sum(tmp_graph, 1)) ;
    if ismember(max_node, sample_index)
        break
    end
    sample_index(current_sample) = max_node ;
    if rcond(full_graph(sample_index, sample_index)) > 1e-9
        tmp_graph = abs(full_graph - full_graph(:, sample_index) * (full_graph(sample_index, sample_index) \ full_graph(sample_index, :))) ;
    else
        break
    end
end
if length(sample_index) < V_d
    current_sample = length(sample_index) ;
    tmp_index = (1 : V_full) ;
    tmp_index(sample_index) = [] ;
    rand_index = randperm(length(tmp_index), V_d - current_sample) ;
    rand_sample = tmp_index(rand_index) ;
    sample_index(current_sample + 1 : V_d) = rand_sample;
end
sample_graph = full_graph(:, sample_index) ;
end

