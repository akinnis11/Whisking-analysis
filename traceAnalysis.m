function traceAnalysis_c(filename,plotFig) %Use filename + .mat extension

close all

%Loading all the necessary files
vidfile = sprintf('%s',filename(1:end-4),'.mp4');
disp('Loading video file...')
vidobj = VideoReader(vidfile); %Load the video file

measurefile = sprintf('%s',filename(1:end-4),'.measurements');
disp('Loading measurements file...')
measurements = LoadMeasurements(measurefile);

if exist(filename,'file') ~= 2 %Create a mat file if doesn't already exist and run the 'setup' analysis
    analyzeWhiskersSetup(filename,vidobj);
end
disp('Loading mat file...')
load(filename)

frame = zeros(1,size(measurements,1));
for i = 1:size(measurements,1)
    frame(i) = measurements(i).fid;
end

minLength = 80;
%whiskerFrames(1) = struct('cdata',zeros(vidHeight,vidWidth,'uint8'),...
%    'colormap',[]);

whiskMov = struct('cdata',zeros(vidHeight,vidWidth,'uint8'),...
    'colormap',[]);

IRledSignal = zeros(3,nFrames);

z = 1;
%plotFig = 'y';

if isequal(plotFig,'y')
    h = figure(1);
    set(0,'CurrentFigure',h)
end

percDone = 0;
disp(sprintf('%s',num2str(percDone),'% done...'))

%%%%%%%%% Analyze whisker tracing information by frame %%%%%%%
for i = 1:nFrames
    
    percDone = floor(100*(i/nFrames));
    percDoneLast = floor(100*((i-1)/nFrames));
    if isequal(percDone,percDoneLast) == 0
        disp(sprintf('%s',num2str(percDone),'% done...'))
    end
    whiskMov.cdata = read(vidobj,i);
    
    temp = whiskMov.cdata(IRledLocation(2),IRledLocation(1),:);
    IRledSignal(:,i) = reshape(temp,3,1);
    
    if isequal(plotFig,'y')
        %Plot the current frame
        figure(1)
        image(whiskMov.cdata)
        hold on
        title(sprintf('%s','Frame ',num2str(i-1)))
        line([xThresh1 xThresh2],[yThresh1 yThresh1])
        line([xThresh1 xThresh1],[yThresh1 yThresh2])
        line([xThresh2 xThresh2],[yThresh1 yThresh2])
        line([xThresh1 xThresh2],[yThresh2 yThresh2])
        plot(faceEdgeX,faceEdgeY,'-y')
        plot(IRledLocation(1),IRledLocation(2),'.y','MarkerSize',20)
    end
    indList = find(frame == (i-1)); %Find indices in 'whiskers' that correspond to current frame
    
    whiskAngle = [];
    
    %Going through each object in this frame for analysis
    for j = 1:length(indList)
        t = indList(j);
        
        %Find the 'follicle' point, which is the closest point on the
        %whisker to the face edge
        follicleX = measurements(t).follicle_x;
        follicleY = measurements(t).follicle_y;
        whiskerTipX = measurements(t).tip_x;
        whiskerTipY = measurements(t).tip_y;
        minFollicleDistance = findFollicle_b(follicleX,follicleY,faceEdgeX,faceEdgeY);
        whiskerCurve = measurements(t).curvature;
        whiskAngle =  faceAngle + abs(measurements(t).angle) + 90;
        
        %Check each traced object and determine if potential whisker or
        %not; returns '0' if not whisker and '1' if it is
        isWhisker(j) = checkTrace_b(follicleX,follicleY,whiskerTipX,whiskerTipY,xThresh1,yThresh1,xThresh2,yThresh2,...
            minLength,faceEdgeX,faceEdgeY,minFollicleDistance,faceAngle,whiskAngle);
        
        if isequal(isWhisker(j),0)
            if isequal(plotFig,'y')
                plot([follicleX whiskerTipX],[follicleY whiskerTipY],'-r','MarkerSize',20) %Rejected whisker objects
            end
            measurements(t).fid = NaN;
            measurements(t).follicle_x = [];
            measurements(t).follicle_y = [];
            measurements(t).tip_x = [];
            measurements(t).tip_y = [];
            
            measurements(t).angle = NaN;
            whiskerAngles(j) = NaN;
            noWhiskerInd(z) = t; %Saving a list of the objects that are not whiskers, in order to delete later
            z = z+1;
            
        else
            if isequal(plotFig,'y')
                plot([follicleX whiskerTipX],[follicleY whiskerTipY],'-g','MarkerSize',20)
                plot([follicleX whiskerTipX],[follicleY whiskerTipY],'.b','MarkerSize',20)
            end
            
            whiskerAngles(j) = whiskAngle;
        end
        
    end
    
    %Save info for all objects in current frame
    temp = find(isnan(whiskerAngles));
    whiskerAngles(temp) = [];
    whiskerPosition(i) = mean(whiskerAngles);
    whiskerPosition_median(i) = median(whiskerAngles);
    whiskerCurvature(i) = mean(whiskerCurve);
    
    if isequal(plotFig,'y')
        figure(1)
        title(sprintf('%s','Frame ',num2str(i-1),': Median whisker angle = ',num2str(whiskerPosition_median(i)),' degrees'))
        hold off
    end
end
whiskersAll = measurements(1:t);
whiskersAll(noWhiskerInd) = [];

disp('Saving .mat file...')
save(filename,'whiskersAll','frame','noWhiskerInd','whiskerPosition','whiskerCurvature','whiskerPosition_median','IRledSignal','-append','-v7.3')

end


