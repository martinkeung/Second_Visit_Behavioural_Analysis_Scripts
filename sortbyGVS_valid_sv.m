% Sort behaviour results in the GVS stim order of 1–9
%
% INPUT
%   - input_ds - structural variable containing all the experiment information
%                (results are in a pseudo-random order)
%
% @author: Martin Keung
% @last update: 2024.06.17

function res_sorted = sortbyGVS_valid_sv(res)

pink_stim_ind = find(res.gvs_order == 1); % find indices of pink stim 
sham_stim_ind = find(res.gvs_order == 2); % find indices of sham
best_stim_ind = find(res.gvs_order == 3); % find indices of best stim

sorted_gvs_ordering = [pink_stim_ind, sham_stim_ind, best_stim_ind];

% [~, sorted_gvs_ordering] = sort(res.gvs_order);

res_sorted.response_threshold = res.response_threshold(sorted_gvs_ordering);
res_sorted.peak_time = res.peak_time(sorted_gvs_ordering);
res_sorted.peak_loc = res.peak_loc(sorted_gvs_ordering);
res_sorted.peak_time_base = res.peak_time_base;
res_sorted.peak_pressure = res.peak_pressure(sorted_gvs_ordering);
res_sorted.peak_pressure_base = res.peak_pressure_base;
res_sorted.sqrwd_success = res.sqrwd_success(sorted_gvs_ordering);
res_sorted.sqrwd_prem = res.sqrwd_prem(sorted_gvs_ordering);
res_sorted.sqrwd_vlate = res.sqrwd_vlate(sorted_gvs_ordering);
res_sorted.sqrwd_delay = res.sqrwd_delay(sorted_gvs_ordering);
res_sorted.sqrwd_nosq = res.sqrwd_nosq(sorted_gvs_ordering);
res_sorted.catchtrials = res.catchtrials(sorted_gvs_ordering);
res_sorted.others = res.others(sorted_gvs_ordering);
res_sorted.goodtrials = res.goodtrials(sorted_gvs_ordering);
res_sorted.sqtype = res.sqtype(sorted_gvs_ordering);
res_sorted.reward = res.reward(sorted_gvs_ordering);
res_sorted.cue_times = res.cue_times(sorted_gvs_ordering);
res_sorted.sqdata_all = res.sqdata_all(sorted_gvs_ordering);
res_sorted.markers_all = res.markers_all(sorted_gvs_ordering);
res_sorted.gvs_order = 1:3;


