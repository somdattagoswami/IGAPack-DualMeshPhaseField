
locationFolder = pwd;
filePattern = fullfile(locationFolder,'*.PNG');% Get a directory listing
imageFiles = dir(filePattern);

writerObj = VideoWriter('shearMovie.avi');% Open the video writer object.
writerObj.FrameRate = 55;
open(writerObj);

for iframe = 1 : length(imageFiles)
    filename = sprintf('Loadstep%d.png',iframe);
    frame = imread(filename);
    writeVideo(writerObj,frame);% Write this frame out to the AVI file
end

close(writerObj);% Close down the video writer object to finish the file
