function plot_waveform_overlay(individual_wfs, mean_wf, color, neuron_type)
n_samples = size(individual_wfs, 2);
if strcmp(neuron_type, 'INT'); individual_linewidth = 1;
else;                          individual_linewidth = 0.75;
end
for i = 1:size(individual_wfs, 1)
    plot(1:n_samples, individual_wfs(i, :), 'Color', [color, 0.1], 'LineWidth', individual_linewidth);
end
plot(1:n_samples, mean_wf, 'Color', color, 'LineWidth', 4);
xlim([1, n_samples]); ylim([-1.1, 0.6]);
end
