function plot_stim(onset, duration, color)
% add patches on plot to indicate stimuli
% INPUT
% onset[vector], index of stim onset
% duration[scalar/vector], stim duration
% color, color of patch

if nargin < 3; color = [0.6,0.6,0.6]; end

onset = onset(:);
duration = duration(:);
nstim = size(onset,1);

if nstim > 1 && size(duration,1) == 1
    duration = repmat(duration,nstim,1);
end

if nstim > 1 && size(color,1) == 1
    color = repmat(color,nstim,1);
end

for i = 1:nstim
    hold on
    a = get(gca,'ylim');
    x = [onset(i) onset(i)+duration(i) onset(i)+duration(i) onset(i)];
    y = [a(1) a(1) a(2) a(2)];
    patch(x,y,color(i,:),'facealpha',0.3,'edgecolor','none')
    hold off

end