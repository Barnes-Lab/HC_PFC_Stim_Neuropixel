function all_SP_filtered = load_and_filter_spikes(baseDir, rats, ages)
%Load SP_good.mat per rat and filter each unit to
% good spikes only (via IX_of_good_spikes). Attaches rat ID, age, and good-
% spike count metadata to each unit.
%

% OUTPUT
%   all_SP_filtered 
% SS 2026

all_SP_filtered = [];

for r = 1:length(rats)
    rat = rats{r};
    age = ages{r};
    fprintf('Loading rat %s (age: %s)...\n', rat, age);

    spFile = fullfile(baseDir, rat, 'Denoised', 'SP_good.mat');
    if ~exist(spFile, 'file')
        warning('SP file not found for rat %s', rat); continue;
    end
    SP = load(spFile, 'SP_good').SP_good;

    for u = 1:length(SP)
        if isfield(SP(u), 'IX_of_good_spikes')
            good_idx = logical(SP(u).IX_of_good_spikes);
        else
            warning('No IX_of_good_spikes for rat %s unit %d, using all spikes', rat, u);
            good_idx = true(size(SP(u).t_uS));
        end

        SP_filt                  = SP(u);
        SP_filt.t_uS             = SP(u).t_uS(good_idx);
        SP_filt.rat              = rat;
        SP_filt.age              = age;
        SP_filt.n_good_spikes    = sum(good_idx);
        SP_filt.n_total_spikes   = length(good_idx);
        SP_filt.pct_good_spikes  = 100 * sum(good_idx) / length(good_idx);
        SP_filt.good_spike_idx   = good_idx;

        all_SP_filtered = [all_SP_filtered, SP_filt];
    end
end

n_young = sum(strcmp({all_SP_filtered.age}, 'Y'));
n_old   = sum(strcmp({all_SP_filtered.age}, 'O'));
fprintf('Loaded %d neurons total (Young: %d, Old: %d)\n', ...
    length(all_SP_filtered), n_young, n_old);
end
