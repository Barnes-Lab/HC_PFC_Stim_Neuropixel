function export_I80_table_clean(I80_DATA, output_path, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exports a flat I80 table from I80_DATA.mat, 
% SS 2026

%% PARAMETERS TO CHANGE
PRM_overwrite = false;
Extract_varargin;

if exist(output_path, 'file') && ~PRM_overwrite
    fprintf('I80 table CSV already exists, skipping export: %s\n', output_path);
    return;
end

rats     = fieldnames(I80_DATA);
rats     = rats(startsWith(rats, 'Rat'));
sess_nms = {'iHC', 'vHC'};

n_max    = numel(rats) * numel(sess_nms);
Rat       = strings(n_max, 1);
Stim_Site = strings(n_max, 1);
Layer     = strings(n_max, 1);
I80       = nan(n_max, 1);
n         = 0;

for ir = 1:numel(rats)
    ratFld = rats{ir};
    rat_id = erase(ratFld, 'Rat');
    for is = 1:numel(sess_nms)
        sn = sess_nms{is};
        if ~isfield(I80_DATA.(ratFld), sn); continue; end
        d = I80_DATA.(ratFld).(sn);
        if ~isfield(d, 'I80'); continue; end
        n            = n + 1;
        Rat(n)       = rat_id;
        Stim_Site(n) = sn;
        Layer(n)     = "L5";
        I80(n)       = d.I80;
    end
end

T = table(Rat(1:n), Stim_Site(1:n), Layer(1:n), I80(1:n), ...
          'VariableNames', {'Rat', 'Stim_Site', 'Layer', 'I80'});

[out_dir, ~, ~] = fileparts(output_path);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

writetable(T, output_path);
fprintf('Exported I80 table (%d rows): %s\n', n, output_path);

end
