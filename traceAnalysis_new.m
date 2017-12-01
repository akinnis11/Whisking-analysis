function traceAnalysis_new(filename,parameters) %Use filename

% December 1, 2017 AK: This is the most recent version of program to
% analyze the output from 'whiski' to calculate whisking from mouse whisker video. Before running
% this, run the analyzeWhiskersSetup program to obtain the appropriate
% parameters.

%Input arguments:
%Filename - enter the filename (with or without the extension) for the whisker
%video. All variables will get saved in a .mat file with the same name

%parameters: should have already run the "analyzeWhiskersSetup" which will
%generate a .mat file containing a structure with necessary parameters. Load
%this file and pass the parameters structure as an input. 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
global nFrames vidobj faceEdgeX faceEdgeY faceAngle measurements xThresh1 yThresh1 xThresh2 yThresh2

%Initialize variables
ind=1;
frame = 0;
whiskerPosition_median = zeros(1,nFrames);
IRledSignal = zeros(1,nFrames);

%Remove any filename extension if exists
if isequal(filename(end-3),'.')
    filename = filename(1:end-4);
end

%Create a mat file if doesn't already exist and run the 'setup' analysis
if ~exist([filename '.mat'],'file')
    analyzeWhiskersSetup(filename);
end

%Load necessary variables
loadFiles()
loadParams(parameters)

%Analyze whisker tracking information by frame
while hasFrame(vidobj)
    
    frame = frame+1; percCounter(frame,1,nFrames) %keeps track of progress
    
    %Get the current frame
    im = readFrame(vidobj);
    
    %Get all traced objects in the current frame
    [indList,ind] = getInds(frame-1,ind);
    
    %Goes through each object in this frame for analysis; returns median
    %angle per frame & coordinates of objects for plotting
    [whiskerPosition_median(frame),isWhisker,fX,fY,tX,tY] = checkWhiskerObjs(indList);
    
    %Plot frame output (at least for first few hundred frames)
    plotOutput(isWhisker,fX,fY,tX,tY)
    
    %Save the mean value of the frame to track IR led flashes
    IRledSignal(frame) = readIRLedpixel(im);
    
end

disp('Saving .mat file...')
save(filename,'IRledSignal','whiskerPosition_median','-append')

plotFinalOutput(whiskerPosition_median)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function loadParams(p)
        xThresh1 = p.xThresh1;
        xThresh2 = p.xThresh2;
        yThresh1 = p.yThresh1;
        yThresh2 = p.yThresh2;
        faceEdgeX = p.faceEdgeX;
        faceEdgeY = p.faceEdgeY;
        faceAngle = p.faceAngle;
        
        nFrames = vidobj.Duration*vidobj.FrameRate;
    end

    function loadFiles()
        %Loading video file
        disp('Loading video file...')
        vidobj = VideoReader([filename '.mp4']); %Load the video file
        %Load measurements file
        measurements = LoadMeasurements([filename '.measurements']);

    end

    function [indList,b] = getInds(f,b)
        first = b;
        while measurements(b).fid == f && (b < length(measurements))
            b = b+1;
        end
        last = b-1;
        indList = first:last;
    end

    function [whiskerPosition_median,isWhisker,fX,fY,tX,tY] = checkWhiskerObjs(iList)
        
        whiskerAngles = zeros(1,length(iList));
        isWhisker = ones(1,length(iList));
        fX = zeros(1,length(iList));
        fY = zeros(1,length(iList));
        tX = zeros(1,length(iList));
        tY = zeros(1,length(iList));
        
        for j = 1:length(iList)
            t = iList(j);
            
            %Find the 'follicle' point, which is the closest point on the
            %whisker to the face edge
            fX(j) = measurements(t).follicle_x;
            fY(j) = measurements(t).follicle_y;
            tX(j) = measurements(t).tip_x;
            tY(j) = measurements(t).tip_y;
            minFollicleDistance = findFollicle_b(fX(j),fY(j),faceEdgeX,faceEdgeY);
            whiskerAngles(j) =  faceAngle + abs(measurements(t).angle) + 90;
            
            %Check each traced object and determine if potential whisker or
            %not; returns '0' if not whisker and '1' if it is
            isWhisker(j) = checkTrace_c(fX(j),fY(j),tX(j),tY(j),xThresh1,yThresh1,xThresh2,yThresh2,...
                faceEdgeX,faceEdgeY,minFollicleDistance,faceAngle,whiskerAngles(j));
        end
        %Just return the median angle per frame to keep
        whiskerPosition_median = median(whiskerAngles(isWhisker == 1));
    end

    function plotOutput(isWhisker,fx,fy,tx,ty)
      
        if frame < 201
            h = figure(1);
            set(0,'CurrentFigure',h)
            image(im)
            hold on
            title(sprintf('%s','Frame ',num2str(frame-1)))
            w = xThresh1 - xThresh2;
            h = yThresh2 - yThresh1;
            rectangle('Position',[xThresh2 yThresh1 w h],'EdgeColor','b')
            plot(faceEdgeX,faceEdgeY,'-y')
            
            line([fx(isWhisker == 1)' tx(isWhisker == 1)']',[fy(isWhisker == 1)' ty(isWhisker == 1)']','Color','g') %Kept whisker objects
            plot([fx(isWhisker == 1)' tx(isWhisker == 1)']',[fy(isWhisker == 1)' ty(isWhisker == 1)']','.b','MarkerSize',20)
            line([fx(isWhisker == 0)' tx(isWhisker == 0)']',[fy(isWhisker == 0)' ty(isWhisker == 0)']','Color','r') %Rejected whisker objects
            
            drawnow
            hold off
        end
    end

    function signal = readIRLedpixel(im2)
        signal = mean(im2(:));
    end

    function plotFinalOutput(whiskerPosition_median)
        figure(1)
        ax(1) = subplot(2,1,1);
        plot(whiskerPosition_median,'-b')
        axis tight
        ylabel('Median whisker angle per frame (deg)')
        
        subplot(2,1,2)
        hist(whiskerPosition_median,50)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','b','EdgeColor','w')
        axis tight
        title('histogram of median whisker angle')
    end
end


