function ds = create_ds_valid_v1_sv(input_ds,sqtper)
% This code should be used for participants PSPD007 and earlier.
% Create a structural variable containing threshold, response time,
% success, and bad trial information based on
% 1) squeeze timeseries (combined_data_matrix)
% 2) flags (flags_matrix)
%
% INPUT
%   - input_ds - structural variable containing all the experiment information
%                (i.e., input_ds = load('experiment_file.mat'))
%
% @author: Varsha Sreenivasan
% @last update: 2024.03.19

warning('off', 'signal:findpeaks:largeMinPeakHeight')
% Find out how many blocks were completed
total_block = size(input_ds.go_times_matrix, 1);
total_trials = size(input_ds.go_times_matrix, 2);

% Note for matrix dimension
% Row = trials
% Column = GVS conditions
fs = 1000;
dsize = [total_block, 30];
base_squeeze = input_ds.baseline_squeeze_during_go_x_matrix;
sqtype = zeros(dsize);

base_squeeze_array = zeros(1, 20); % Hard code 20 trials because it's always going to be 20 trials

if size(base_squeeze, 2) < 20
    base_squeeze_array = 0.5*ones(1,20); % If there's no baseline squeeze data, assign it to be all ones
else
    for i = 1:20
        base_squeeze_array(i) = max(base_squeeze(:, i));
    end
end

mpeak = nanmean(base_squeeze_array); % Computes the mean of the base squeeze array

% Preallocate res variables
res.bad_trials = zeros(total_block, total_trials);
res.response_threshold = nan(total_block, total_trials);
res.cue_times = nan(total_block, total_trials, 3);

% Initialize cell arrays to hold each trial's data
trials = cell(total_block, total_trials);
trials_markers = cell(total_block, total_trials);


for i = 1:total_block % Run through each block
    sqdata = input_ds.combined_data_matrix{i}(:, 1); % first channel is squeeze data
    markers = input_ds.combined_data_matrix{i}(:, 5); % fifth channel is marker data

    % Precompute trial start and end indices
    trialstart_inds = find(markers == 11 | markers == 12 | markers == 13);
    trialend_inds = find(markers == 51);
    markers_positive = find(markers > 0);


    for j = 1:total_trials
        start_idx = trialstart_inds(j); % finds the index of the cue flag
        end_idx = trialend_inds(j); % finds the index of the iti flag
        trial_data = sqdata(start_idx:end_idx); % takes the data between the cue and iti flag
        trials{i, j} = trial_data; % put the data into the cell corresponding to it's trial and block number
        bs = nanmean(trials{i, j}); % gets the mean of the squeeze data but without the NaNs
        start_marker_idx = find(markers_positive == start_idx); % get index of cue flag in trial
        end_marker_idx = find(markers_positive == end_idx); % get index of iti flag in trial
        if sum(markers_positive(start_marker_idx:end_marker_idx)) > start_idx % this will make sure we don't get a 0 value for markers_positive
            trials_markers{i, j} = markers_positive(start_marker_idx:end_marker_idx) - start_idx; % get time positions of markers in each trial
        else
            trials_markers{i, j} = 99999; % this will let us know if theres' an issue
        end
        markerinds = trials_markers{i, j} / fs; % Converting into seconds

        % duration for 1) reward presentation, 2) anticipation, and 3) go
        markerdiff = diff(markerinds);

        % Find bad trials
        if length(trials_markers{i, j}) ~= 5
            % number of markers are not proper
            res.bad_trials(i, j) = 1;  % Indexed by i and j
            trials_markers{i, j} = nan(5, 1);
        end

        % Find threshold used
        res.response_threshold(i, j) = input_ds.go_times_matrix(j); % feedback - go flag

        % Save cue duration
        res.cue_times(i, j, 1) = markerdiff(1); % 1) reward presentation
        res.cue_times(i, j, 2) = markerdiff(2); % 2) anticipation
        res.cue_times(i, j, 3) = markerdiff(3); % 3) go

        % Find peak and compute response time
        if sum(isnan(markerinds)) == 0
            % check for squeeze
            [pks, locs] = findpeaks(trials{i, j}, 'MinPeakHeight', max(trials{i, j}) * 0.8, 'MinPeakDistance', 5000);

            if sum(pks - bs > sqtper*mpeak)>0 % If a full squeeze in the Go period

                % There may be trials where there was a peak before Go and
                % one was initiated in Go in which case we will have a flag 41 but findpeaks may find
                % the first peak instead. We need to account for that case
                % (see under second elseif)
                if length(locs) == 1 && locs > markerinds(3)*fs && (locs/fs - markerinds(3)) <= 1.2*res.response_threshold(i,j) % peak in the Go + 0.2*Go time period
                    res.peak_time(i,j) = locs/fs - markerinds(3);
                    res.peak_loc(i,j) = locs;
                    res.peak_pressure(i,j) = pks - bs;
                    sqtype(i,j) = 1; % successful squeeze
                    txt = 'Successful squeeze';
                elseif length(locs) == 1 && locs > markerinds(3)*fs && (locs/fs - markerinds(3)) > 1.2*res.response_threshold(i,j) % peak outside Go + 0.2*Go time period
                    res.peak_time(i,j) = locs/fs - markerinds(3);
                    res.peak_loc(i,j) = locs;
                    res.peak_pressure(i,j) = pks - bs;
                    sqtype(i,j) = 2; % valid late squeeze
                    txt = 'Valid late squeeze';
                elseif length(locs) == 1 && locs <= markerinds(3)*fs
                    % If there is a peak before Go and one was initiated in Go and findpeaks found the first (premature peak), we findpeaks again
                    [pks_pre,locs_pre] = findpeaks(x,'MinPeakHeight',max(x)*0.8,...
                        'MinPeakDistance',1000); % find peaks again, this time a little more liniently
                    if length(locs_pre) == 2 && locs_pre(1) <= markerinds(3)*fs && locs_pre(2) > markerinds(3)*fs % If we get two peaks and first peak is before Go but the second after Go

                        % check if second peak is a full squeeze
                        if (pks_pre(2) - bs > sqtper*mpeak)
                            % Take the second squeeze
                            if (locs_pre(2)/fs - markerinds(3)) <= 1.2*res.response_threshold(i,j) % peak in the Go + 0.2*Go time period
                                res.peak_time(i,j) = locs_pre(2)/fs - markerinds(3);
                                res.peak_loc(i,j) = locs_pre(2);
                                res.peak_pressure(i,j) = pks_pre(2) - bs;
                                sqtype(i,j) = 1; % successful squeeze
                                txt = 'Successful squeeze';
                            elseif (locs_pre(2)/fs - markerinds(3)) > 1.2*res.response_threshold(i,j) % peak outside Go + 0.2*Go time period
                                res.peak_time(i,j) = locs_pre(2)/fs - markerinds(3);
                                res.peak_loc(i,j) = locs_pre(2);
                                res.peak_pressure(i,j) = pks_pre(2) - bs;
                                sqtype(i,j) = 2; % valid late squeeze
                                txt = 'Valid late squeeze';
                            end
                        else
                            res.peak_time(i,j) = nan;
                            res.peak_loc(i,j) = nan;
                            res.peak_pressure(i,j) = nan;
                            sqtype(i,j) = 4; % no squeeze
                            txt = 'No squeeze';
                        end

                    elseif length(locs_pre) == 1 && sum(pks_pre(1) > sqtper*mpeak) > 0% If only one squeeze before Go and full squeeze
                        % Reject as premature squeeze
                        res.peak_time(i,j) = nan;
                        res.peak_loc(i,j) = nan;
                        res.peak_pressure(i,j) = nan;
                        sqtype(i,j) = 3; % premature squeeze
                        txt = 'Premature squeeze';
                    elseif length(locs_pre) > 2
                        txt = 'other'; % sqtype will be 0 in this case, which will be regarded as other subsequently
                    else
                    end

                elseif length(locs) == 2 && locs(1) <= markerinds(3)*fs && locs(2) > markerinds(3)*fs

                    % Take the second squeeze
                    if (locs(2)/fs - markerinds(3)) <= 1.2*res.response_threshold(i,j) % peak in the Go + 0.2*Go time period
                        res.peak_time(i,j) = locs(2)/fs - markerinds(3);
                        res.peak_loc(i,j) = locs(2);
                        res.peak_pressure(i,j) = pks(2) - bs;
                        sqtype(i,j) = 1; % successful squeeze
                        txt = 'Successful squeeze';
                    else % peak outside Go + 0.2*Go time period
                        res.peak_time(i,j) = locs(2)/fs - markerinds(3);
                        res.peak_loc(i,j) = locs(2);
                        res.peak_pressure(i,j) = pks(2) - bs;
                        sqtype(i,j) = 2; % valid late squeeze
                        txt = 'Valid late squeeze';
                    end
                elseif length(locs) == 2 && locs > markerinds(3)*fs
                    % Take the first squeeze
                    if (locs(1)/fs - markerinds(3)) <= 1.2*res.response_threshold(i,j) % peak in the Go + 0.2*Go time period
                        res.peak_time(i,j) = locs(1)/fs - markerinds(3);
                        res.peak_loc(i,j) = locs(1);
                        res.peak_pressure(i,j) = pks(1) - bs;
                        sqtype(i,j) = 1; % successful squeeze
                        txt = 'Successful squeeze';
                    else % peak outside Go + 0.2*Go time period
                        res.peak_time(i,j) = locs(1)/fs - markerinds(3);
                        res.peak_loc(i,j) = locs(1);
                        res.peak_pressure(i,j) = pks(1) - bs;
                        sqtype(i,j) = 2; % valid late squeeze
                        txt = 'Valid late squeeze';
                    end

                elseif length(locs) > 2
                    % Reject as premature squeeze
                    res.peak_time(i,j) = nan;
                    res.peak_loc(i,j) = nan;
                    res.peak_pressure(i,j) = nan;
                    sqtype(i,j) = 3; % premature squeeze
                    txt = 'Premature squeeze';
                end
            elseif sum(pks - bs > sqtper*mpeak)==0 % not full squeeze then this is a no squeeze
                res.peak_time(i,j) = nan;
                res.peak_loc(i,j) = nan;
                res.peak_pressure(i,j) = nan;
                sqtype(i,j) = 4; % no squeeze
                txt = 'No squeeze';
            end

            % mark successful squeeze trials
            if input_ds.reward_matrix(i, j) >= 0 && sqtype(i, j) == 1
                res.sqrwd_success(i, j) = 1; % successful squeeze
            else
                res.sqrwd_success(i, j) = 0; % other
            end

            % mark valid late squeeze trials
            if input_ds.reward_matrix(i, j) >= 0 && sqtype(i, j) == 2
                res.sqrwd_vlate(i, j) = 1; % late squeeze
                res.sqrwd_delay(i, j) = res.peak_time(i, j) - res.response_threshold(i, j); % delay (how much later)
            else
                res.sqrwd_vlate(i, j) = 0; % other
                res.sqrwd_delay(i, j) = 0; % no delay
            end

            % mark premature squeeze trials
            if input_ds.reward_matrix(i, j) >= 0 && sqtype(i, j) == 3
                res.sqrwd_prem(i, j) = 1; % premature squeeze
            else
                res.sqrwd_prem(i, j) = 0; % other
            end

            % mark no squeeze trials
            if input_ds.reward_matrix(i, j) >= 0 && sqtype(i, j) == 4
                res.sqrwd_nosq(i, j) = 1; % no squeeze
            else
                res.sqrwd_nosq(i, j) = 0; % other
            end

            % mark catch trials
            if input_ds.reward_matrix(i, j) < 0
                res.catchtrials(i, j) = 1; % catch trial
            else
                res.catchtrials(i, j) = 0; % other
            end

            % mark other weird invalid trials
            if sqtype(i, j) == 0
                res.others(i, j) = 1;
                res.peak_time(i, j) = nan;
                res.peak_pressure(i, j) = nan;
            else
                res.others(i, j) = 0;
            end

            % mark overall bad trials
            if res.sqrwd_prem(i, j) == 1 || res.sqrwd_nosq(i, j) == 1 || res.catchtrials(i, j) == 1 || res.others(i, j) == 1
                res.goodtrials(i, j) = false;
            else
                res.goodtrials(i, j) = true;
            end
        end
    end
end


% Extra step to repmove any trials that are weird
fpp = isnan(res.peak_pressure);
fpt = isnan(res.peak_loc);

if ~isequal(fpp,fpt)
    fn = find(fpp & ~fpt);
    [tr,bl] = ind2sub(fn,size(fpp));
    for ii = 1:length(tr)
        res.peak_loc(tr(ii),bl(ii)) = nan;
        res.goodtrials(tr(ii),bl(ii)) = 0;
        res.others(tr(ii),bl(ii)) = 1;
    end
end


ds.response_threshold = res.response_threshold';
ds.peak_time = res.peak_time';
ds.peak_loc = res.peak_loc';
ds.peak_time_base = input_ds.starting_response_time_array;
ds.peak_pressure_base = base_squeeze_array;
ds.peak_pressure = res.peak_pressure';
ds.sqrwd_success = res.sqrwd_success';
ds.sqrwd_prem = res.sqrwd_prem';
ds.sqrwd_vlate = res.sqrwd_vlate';
ds.sqrwd_delay = res.sqrwd_delay';
ds.sqrwd_nosq = res.sqrwd_nosq';
ds.catchtrials = res.catchtrials';
ds.others = res.others';
ds.goodtrials = res.goodtrials';
ds.sqtype = sqtype';
ds.reward = input_ds.reward_matrix';
ds.cue_times = res.cue_times;
ds.sqdata_all = trials';
ds.markers_all = trials_markers';
ds.gvs_order = input_ds.true_stim_order';

