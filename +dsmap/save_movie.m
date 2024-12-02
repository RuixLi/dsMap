function save_movie(DSresults,ops)
% save the movie of average responses
% DSresults: the results of DS map
% ops: the parameters of DS map


vidObj = VideoWriter(fullfile(ops.save_dir,'00_averaged_response'),'MPEG-4');
set(vidObj,'FrameRate',30);
set(vidObj, 'Quality',100);
open(vidObj);
avg = reshape(DSresults.trialAvgData,ops.d1,ops.d2,[]);
T = size(avg,3);
clim = [quantile(avg(1:1e5),0.01),quantile(avg(1:1e5),0.99)];
cm = gray(256);
h = figure('Position',[600,300,380,400],'Color','k','MenuBar','none','ToolBar','none','Visible','off');

trialLength = sum(ops.ON_OFF_frames);
dirList = repmat(ops.dir_list,trialLength,1);
dirList = dirList(:);
for i = 1:T
    % try figure(h); catch; break; end
    imagesc(avg(:,:,i),clim);
    axis image off
    set(gca,'Position',[0,0,1,0.95])
    colormap(cm)
    title(sprintf('%d deg : %03d/%d', dirList(i),mod(i,trialLength)+1,trialLength),'color',[1,1,1])
    writeVideo(vidObj,getframe(h));
end
close(h)
close(vidObj);
fprintf('movie saved to %s\n',ops.save_dir)