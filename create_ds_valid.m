% This code should be used for participants aftyer PSPD007
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


function ds = create_ds_valid(input_ds,sqtper)

warning('off', 'signal:findpeaks:largeMinPeakHeight')

% Find out how many blocks were completed
total_block = size(input_ds.go_times_matrix,1);
total_trials = size(input_ds.go_times_matrix,2);

% Note for matrix dimension
% Row = trials
% Column = GVS conditions
fs = 1000;
dsize = [total_block,total_trials];
sqdata_all = squeeze(input_ds.combined_data_matrix_processed(:,1,:,:));
markers_all = squeeze(input_ds.combined_data_matrix_processed(:,5,:,:));
base_squeeze = input_ds.baseline_squeeze_during_go_x_matrix;
sqtype = zeros(dsize);

% get the baseline squeeze data from each of the 20 trials into an array to calculate the mean afterwards
base_squeeze_array = zeros(1, 20);
for i = 1:2*total_trials
    base_squeeze_array(i) = max(base_squeeze(:, i)) - base_squeeze(1, i);
end

mpeak = nanmean(base_squeeze_array);

for i = 1:dsize(1)
    sqdata = squeeze(sqdata_all(:,i,:));
    markers = squeeze(markers_all(:,i,:));
    x = []; m = [];
    for j = 1:dsize(2) % run through each trial
        x = sqdata(:,j);
        bs = nanmean(x);

        m = markers(:,j);
        markerinds = find(m > 0 & ~isnan(m))/fs; % we want non-zero elements + NaNs, seems like processed code already divides by the sampling rate so removed fs as well
        % duration for 1) reward presentation, 2) anticipation, and 3) go
        markerdiff = diff(markerinds);

        if length(markerinds) ~= 5 && ~((i == 3 && j == 10) || (i == 6 && j == 10) || (i == 9 && j == 10))
            % any trial that is not a trial before the REST should have 5 markers
            res.bad_trials(i,j) = 1;
            markerinds = nan(5,1);
            markerdiff = diff(markerinds);
        end

        % Find threshold used
        res.response_threshold(i,j) = input_ds.go_times_matrix(i,j); %)markerinds(4) - markerinds(3);

        % Save cue duration
        res.cue_times(j,i,1) = markerdiff(1); % 1) reward presentation
        res.cue_times(j,i,2) = markerdiff(2); % 2) anticipation
        res.cue_times(j,i,3) = markerdiff(3); % 3) go

        xp = (1:length(x))'/fs;
        slope = gradient(x)./gradient(xp);

        % Find peak and compute response time
        if sum(isnan(markerinds)) == 0

            % check for squeeze
            [pks,locs] = findpeaks(x,'MinPeakHeight',max(x)*0.8,...
                'MinPeakDistance',5000); % x = sqdata

            % if sum(ismember(m,41))&& sum(pks - bs > sqtper*mpeak)>0 % If a full squeeze in the Go period
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
                % elseif sum(ismember(m,42)) && sum(pks - bs > sqtper*mpeak)>0 && locs > markerinds(3)*fs% If full squeeze after the Go period
                %     res.peak_time(i,j) = locs(1)/fs - markerinds(3);
                %     res.peak_loc(i,j) = locs(1);
                %     res.peak_pressure(i,j) = pks(1) - bs;
                %     sqtype(i,j) = 2.5; % invalid late squeeze
                %     txt = 'Invalid late squeeze';
            elseif sum(pks - bs > sqtper*mpeak)==0 % not full squeeze then this is a no squeeze
                res.peak_time(i,j) = nan;
                res.peak_loc(i,j) = nan;
                res.peak_pressure(i,j) = nan;
                sqtype(i,j) = 4; % no squeeze
                txt = 'No squeeze';
                % f = find(x > sqtval); % in general, find where squeeze was above threshold. This could be full squeeze or just a blip, which would be considered as a no squeeze
                % if ~isempty(f) && f(1) < markerinds(3)*fs && % atleast one point (which means the threshold was crossed) and before the Go signal and full squeeze
                %     res.peak_time(i,j) = nan;
                %     res.peak_loc(i,j) = nan;
                %     res.peak_pressure(i,j) = nan;
                %     sqtype(i,j) = 3; % premature squeeze
                %     txt = 'Premature squeeze';
                % else % no squeeze
                %     res.peak_time(i,j) = nan;
                %     res.peak_loc(i,j) = nan;
                %     res.peak_pressure(i,j) = nan;
                %     sqtype(i,j) = 4; % no squeeze
                %     txt = 'No squeeze';
                % end
            end
        end

        % mark successful squeeze trials - when participant squeezed on reward >=0
        % trials and squeezed within the allotted time
        if input_ds.reward_matrix(i,j) >= 0 && sqtype(i,j) == 1
            res.sqrwd_success(i,j) = 1; % successful squeeze
        else
            res.sqrwd_success(i,j) = 0; % other
        end

        % mark valid late squeeze trials - when participant initiated squeeze within the Go period on reward >=0
        % trials
        if input_ds.reward_matrix(i,j) >= 0 && sqtype(i,j) == 2
            res.sqrwd_vlate(i,j) = 1; % late squeeze
            res.sqrwd_delay(i,j) = res.peak_time(i,j) - res.response_threshold(i,j); % delay (how much later)
        else
            res.sqrwd_vlate(i,j) = 0; % other
            res.sqrwd_delay(i,j) = 0; % no delay
        end

        % % mark invalid late squeeze trials - when participant initiated
        % % squeeze beyond the Go period on reward >=0 trials
        % if input_ds.reward_matrix(i,j) >= 0 && sqtype(i,j) == 2.5
        %     res.sqrwd_ivlate(i,j) = 1; % late squeeze
        % else
        %     res.sqrwd_ivlate(i,j) = 0; % other
        % end

        % mark premature squeeze trials - when participant squeezed on reward >=0
        % trials but before the go period
        if input_ds.reward_matrix(i,j) >= 0 && sqtype(i,j) == 3
            res.sqrwd_prem(i,j) = 1; % premature squeeze
        else
            res.sqrwd_prem(i,j) = 0; % other
        end


        % mark no squeeze trials - when participant didn't squeeze on reward >=0
        % trials
        if input_ds.reward_matrix(i,j) >= 0 && sqtype(i,j) == 4
            res.sqrwd_nosq(i,j) = 1; % no squeeze
        else
            res.sqrwd_nosq(i,j) = 0; % other
        end

        % mark catch trials - reward <=0 trials, squeeze or not
        if input_ds.reward_matrix(i,j) < 0
            res.catchtrials(i,j) = 1; % catch trial
        else
            res.catchtrials(i,j) = 0; % other
        end

        % mark other weird invalid trials (several peaks, peaks but no threshold crossed, etc)
        if sqtype(i,j) == 0
            res.others(i,j) = 1; %
            res.peak_time(i,j) = nan;
            res.peak_loc(i,j) = nan;
            res.peak_pressure(i,j) = nan;
        else
            res.others(i,j) = 0;
        end

        % mark overall bad trials - bad trials are those where the participant
        % squeezed prematurely (before the Go) or didn't squeeze on >=0 reward trials or had invalid late squeezes. We also use this to exclude
        % the $-5 trials (whether the participant squeezed or not),
        % although catch trial are not technically "bad" trials.
        % Good trials are successful squeezes, valid late squeezes,
        if res.sqrwd_prem(i,j) == 1 || res.sqrwd_nosq(i,j) == 1 || res.catchtrials(i,j) == 1 || res.others(i,j) == 1
            res.goodtrials(i,j) = false;
        else
            res.goodtrials(i,j) = true;
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
% ds.sqrwd_ivlate = res.sqrwd_ivlate';
ds.sqrwd_delay = res.sqrwd_delay';
ds.sqrwd_nosq = res.sqrwd_nosq';
ds.catchtrials = res.catchtrials';
ds.others = res.others';
ds.goodtrials = res.goodtrials';
ds.sqtype = sqtype';
ds.reward = input_ds.reward_matrix';
ds.cue_times = res.cue_times;
ds.sqdata_all = permute(sqdata_all, [1,3,2]);
ds.markers_all = permute(markers_all, [1,3,2]);
ds.gvs_order = input_ds.true_stim_order;
