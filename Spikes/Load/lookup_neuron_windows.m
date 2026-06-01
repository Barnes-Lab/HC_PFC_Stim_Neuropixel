function [mono, poly] = lookup_neuron_windows(W, rat, sess, region)
% LOOKUP_NEURON_WINDOWS  Mono/poly windows for one (rat, sess, region).
% Returns [] if no row found.


m = strcmp(W.Rat, rat) & strcmp(W.Sess, sess) & strcmp(W.Region, region);
if ~any(m); mono = []; poly = []; return; end
mono = [W.Mono_Start(m), W.Mono_End(m)];
poly = [W.Poly_Start(m), W.Poly_End(m)];
end
