function cm = cmap_white_to_red(n, gamma)
% N-step white-to-red colormap, gamma controls midtone
% taper 

if nargin < 1; n     = 256; end
if nargin < 2; gamma = 1;   end

t  = linspace(0, 1, n)';
cm = [ones(n, 1), (1 - t).^gamma, (1 - t).^gamma];
end
