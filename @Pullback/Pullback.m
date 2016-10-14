classdef Pullback < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = private)
        pbfn;
        fileInfo;
        numOfFrames;
        %data;
        frames;
        centerline;
       
        attributes;
        
    end
    properties % Public access
       
    end
    
    methods
        function obj = Pullback(fn)    % Constructor
            
            obj.centerline = [];
            obj.pbfn = fn;
            obj.fileInfo = imfinfo(fn);
            obj.numOfFrames = size(obj.fileInfo, 1); % get the number of frames in the file  
            
            %firstFrame = imread(fn, 1);  % read the first frame              
            %obj.data = uint8(zeros([size(firstFrame) 1 obj.numOfFrames])); % allocate a 4D matrix to hold the movie            
            %obj.data(:, :, :, 1) = firstFrame;            
            
            data(1:obj.numOfFrames) = Frame(); % initialize array of classes for the r-theta pullback
            convert2XYasYouLoad = true;
           
            % now, read the rest of the frames
            h = waitbar(0,'Loading pullback...');
            for i=1:obj.numOfFrames
                %obj.data(:, :, :, i) = imread(fn, i);
                data(i) = Frame(i, obj.fileInfo(1).Width, obj.fileInfo(1).Height, fn);   
                
                waitbar(i/obj.numOfFrames,h,'Loading pullback...');
            end
            obj.frames = data;
            close(h);
        end %constructor
        % -----------------------------------------------------------------
        function compareCartesianImages(obj, pbfn, frameList)
            for i=frameList                
                obj.frames(i).compareCartesianImages(pbfn);
            end
        end
        % -----------------------------------------------------------------
        function savePullbackAsImages(obj, dir, frameList)
            
            for i=frameList                
                obj.frames(i).saveAsTiff(dir);
            end
        end
        % -----------------------------------------------------------------
        % Perform pre-processing on the data
        % (cath correction->removeSpeckle->backBorder ID->subtract baseline->lumenDetect->pixelShift-> & clip to zero->ln(I+1))
        function preProcess(obj, frameList, pbfn, ui)          
            
            ppFlags = ui.getPreProcessingFlags();
            gFlags = ui.getGlobalFlags();
           
             %% Check if any global operation is required prior to processing
            if bitget(gFlags, Globals.REMOVE_GUIDE_WIRE)
                obj.removeGuideWire(frameList);
            end                
            
            for frameNum = frameList
                if bitget(ppFlags, Globals.REMOVE_SPECKLE_BIT)
                    obj.frames(frameNum).removeSpeckle('Enhanced Lee', pbfn, 0);
                    %obj.removeSpeckle('Enhanced Lee', frameNum, pbfn, 0); % 'Frost'
                end
                if bitget(ppFlags, Globals.REMOVE_BASELINE_BIT)
                    obj.frames(frameNum).removeBaseline(0);
                    %obj.removeBaseline(frameList, false); % THIS CAUSES THE IMAGE
                    % TO DARKEN.  SOLVE IT !! (the minimum number in the image
                    % before baseline removal is zero.
                end
                if bitget(ppFlags, Globals.CATH_CORRECT_BIT);
                    obj.frames(frameNum).correctForCatheter(0);
                    %obj.correctForCatheter(frameList, false); %This causes for stripes to appear.  Need to check why !!!!!!!!!!!!
                end
                if bitget(ppFlags, Globals.EXTRACT_LUMEN_BIT)                    
                    obj.frames(frameNum).computeLumenBorder(0);
                    %obj.computeLumenBorder(frameList, 0);
                end
                if bitget(ppFlags, Globals.PIX_SHIFT_BIT)
                    I = obj.frames(frameNum).pixShift(0);
                    %obj.pixShift(frameList, false);
                end
                if bitget(ppFlags, Globals.EXTRACT_BACKB_BIT)
                    obj.frames(frameNum).computeBackBorder(0);
                    %obj.computeBackBorder(frameList, false);
                end                
            end
            
           
        end        
        % -----------------------------------------------------------------
        % This function removes the guidewire images.  I replaced the 
        % original image data with the file so
        % that any follow up operation on the image will not have the
        % guidewire and thus will affect the results.  I modeified the
        % function guidewire_longi() from Hong's code
        function [gw_pos, border, Img_longi] = removeGuideWire(obj, frameList)
            
            %[dir, name, ext] = fileparts(obj.pbfn);
            %[gw_pos, border, Img_longi] = guidewire_longi([dir '\'], [name ext], frameList(1), frameList(end));
            
            M= obj.fileInfo(1).Height;
            start_frame = frameList(1);
            end_frame = frameList(end);
            Img_longi=zeros(M,end_frame-start_frame+1);
            
            for frame_num=start_frame:end_frame
                longi = obj.frames(frame_num).scaleImageData();                
                Img_longi(:,frame_num) = longi;
            end
            
            Img_accumulate = Img_longi';
            add_lines = size(Img_longi,1);% changed from 300 to size(Img_longi,1), Apr.17 2014, Hong Lu
            Img_accumulate = [Img_accumulate(:,end-add_lines+1:end) Img_accumulate Img_accumulate(:,1:add_lines)];
            inten_dif1 = zeros(size(Img_accumulate));
            for j=1:size(inten_dif1,2)
                inten_dif1(:,j) = mean(Img_accumulate(:,max(1,j-10):j),2)-mean(Img_accumulate(:,j:min(j+10,size(inten_dif1,2))),2);
            end
            inten_dif2 = -inten_dif1;
            connectivity = 11; % changed from 21 to 31, Apr.17 2014, Hong Lu
            [border1,countour_index1] = DP_fast(inten_dif1', connectivity); 
            countour_index1 = countour_index1+int32(countour_index1<1);
            border1 = border1';
            im_local = inten_dif2(:,min(countour_index1):min(max(countour_index1)+200,size(inten_dif1,2)));
            [border2_local] = DP_fast(im_local', connectivity);% changed from 21 to 31, Apr.17 2014, Hong Lu
            border2_local = border2_local';
            border2 = zeros(size(Img_accumulate));
            border2(:,min(countour_index1):min(max(countour_index1)+200,size(inten_dif1,2)))=border2_local;
            border12 = border1|border2;
 
            % modified to handles cases with more than 360 degree change in guide wire, Apr. 17 2014, Hong Lu
            border = border12(:,1:add_lines)| border12(:,add_lines+1:2*add_lines)|border12(:,2*add_lines+1:end);
            
            % commented out on Apr. 17 2014, Hong Lu
            % border12(:,1:2*add_lines)=temp; border12(:,end-2*add_lines+1:end)=temp;
            % border=border12(:,add_lines+1:end-add_lines);
            if sum(border(:))<1
                kao=zeros(size(Img_longi'));
                kao(:,1:add_lines)=border12(:,end-add_lines+1:end);
                border=kao;
            end
            
            border=border';
            gw_pos=zeros(size(border,2),2);
            for j=1:size(border,2)
                ind=find(border(:,j));
                gw_pos(j,1)=min(ind);
                gw_pos(j,2)=max(ind);
            end
            
            % set the mask of the guidewire in case I want to remove it            
            for i=start_frame:end_frame
                obj.frames(i).setGuidewireMask(gw_pos(i, :));
            end
            
            dispResults = 0;
            if dispResults
                figure, imshow(Img_longi,[]);
            end
            
        end
        % -----------------------------------------------------------------
        % Perform post-processing on the data
        
        function postProcess(obj, frameList, pbfn)
            obj.edgeDetect(frameList, false); % "true" to show the results (works only for a single frame - for development purposes)
            
            
        end
        % -------------------------------------------------------
        function show3D(obj)
            %process a bunch of TIFF files into a volume
            %the original images were had-colored by region to one of 7 colors
            %so the scalar value on the output runs from from 0 to 7
            %file setup
            %filebase = 'C:\Users\rys.ADS\Documents\Projects\TissueCharacterization\OCT_Data\ValidationData\Vessel 52\3D Data\Cartesian\';
            %filebase = 'C:\Users\rys.ADS\Documents\Projects\TissueCharacterization\OCT_Data\ValidationData\Vessel 52\OCT Labels\Vessel52_Segment1_labels.00';
            filebase = 'C:\Users\Erez\Documents\Projects\TissueCharacterization\OCT_Data\ValidationData\Prabu\Vessel 52\OCT Labels\Vessel52_Segment1_labels.00';
            
            
            startFrame = 1;
            endFrame = startFrame + 64;
            %read frames, reduce size, show frames, and build volume
            for i=startFrame:endFrame
                %filename=[filebase, num2str(i) '.tiff'];
                filename=[filebase, num2str(i-1,'%02d') '.tif'];
                %filename=[filebase, num2str(i) 'Calcium' '.tiff'];
                temp=double(imresize(imread(filename), 0.5));
                stack(:,:,i) = temp; %temp~=0; % (temp(:,:,1)==255) + 2*(temp(:,:,2)==255) + 4*(temp(:,:,3)==255);
                imagesc(stack(:,:,i));
                colormap('gray')
                drawnow
            end
            method = 1;
            if method ==1
%                figure
                load mri; % just so I can use their colormap
                %stack = squeeze(D);
%                 colormap('prism')
%                 image_num = 8;
%                 image(stack(:,:,image_num))
%                 axis image;
%                 
%                 x = xlim;
%                 y = ylim;
%                 
%                 cm = brighten(jet(length(map)),-.5);
%                 figure
%                 colormap(cm)
%                 contourslice(stack,[],[],image_num)
%                 axis ij
%                 xlim(x)
%                 ylim(y)
%                 daspect([1,1,1])
%                 
%                 figure
%                 colormap(cm)
%                 contourslice(stack,[],[],[1,12,19,27],8);
%                 view(3);
%                 axis tight                
                
                figure
                %colormap('gray')
                colormap(Globals.cmMatrix);
                Ds = smooth3(stack);
                hiso = patch(isosurface(Ds,1),... % (Ds,5)
                    'FaceColor',[1,.75,.65],...
                    'EdgeColor','none');
                isonormals(Ds,hiso)
                
                hcap = patch(isocaps(stack,1),...
                    'FaceColor','interp',...
                    'EdgeColor','none');
                
                view(35,30)
                axis tight
                daspect([1,1,.4]);
                
                lightangle(45,30);
                lighting gouraud
                hcap.AmbientStrength = 0.6;
                hiso.SpecularColorReflectance = 0;
                hiso.SpecularExponent = 50;

            else
                %patch smoothing factor
                rfactor = 0.125;
                %isosurface size adjustment
                level = .8;
                %useful string constants
                c2 = 'facecolor';
                c1 = 'edgecolor';
                
                %build one isosurface for each of 7 different levels
                %The "slice" matrix takes on one of 7 integer values,
                %so each of the following converts the data to a binary
                %volume variable, then computes the isosurface between the
                %1 and 0 regions
                
                p=patch(isosurface(smooth3(stack==1),level)); % 1 is the lumen border
                reducepatch(p,rfactor)
                set(p,c2,[1,0,0],c1,'none', 'FaceAlpha', 0.25);
                material shiny;
                
                p=patch(isosurface(smooth3(stack==2),level)); % 2 is fiber
                reducepatch(p,rfactor)
                set(p,c2,[0,1,0],c1,'none', 'FaceAlpha', 0.5, 'SpecularColorReflectance', 1.0);
                material shiny;
                
                p=patch(isosurface(smooth3(stack==3),level)); % 3 is lipid
                reducepatch(p,rfactor)
                set(p,c2,[0.427,0.43,0.168],c1,'none', 'FaceAlpha', 0.3);
                material shiny;
                
                p=patch(isosurface(smooth3(stack==5),level)); % 4 is clacium
                reducepatch(p,rfactor)
                set(p,c2,[1,1,1],c1,'none', 'FaceAlpha', 1.0, 'SpecularColorReflectance', 0.7);
                material shiny;
                
                %lighting/image control
                set(gca,'projection','perspective')
                box on
                lighting phong
                light('position',[1,1,1])
                light('position',[-1,-1,-1])
                light('position',[-1, 1,-1])
                light('position',[ 1,-1,-1])
                
                %set relative proportions
                daspect([1,1,1])
                axis on
                set(gca,'color',[1,1,1]*.6)
                
                view(-30,30)
                
                rotate3d on
            end
        end
        % -----------------------------------------------------------------
        function adventitiaEdgeDetect(obj, frameList, dispResultFlag)
             for i = frameList
                obj.frames(i).adventitiaEdgeDetect(dispResultFlag);
            end 
        end
        % -----------------------------------------------------------------
        function computeCentroidLine(obj, frameList, dispResultFlag)
            
            saveToFile_flag = 0; % maybe add this to the gui options
            root = 'C:\Users\rys.ADS\Documents\Projects\TissueCharacterization\OCT_Data\ValidationData\Vessel 52\3D Data\Centroid\';
            
            centroid = [];
            for i = frameList
                c = obj.frames(i).computeCentroidLine(dispResultFlag);                
            end
        end
        % -----------------------------------------------------------------
        function [segmentsAttributes, segmentsImgMask] = isolateEdges(obj, frameList, dispResultFlag)
            segmentsAttributes = {};
            segmentsImgMask = {};
            counter = 1;
            for i = frameList
                [segmentsAttributes{counter}, segmentsImgMask{counter}] = obj.frames(i).isolateEdges(dispResultFlag);
                %figure, imshow(segmentsImgMask{counter},[]);
                %colormap('hot');
            end
        end
        
        % -----------------------------------------------------------------
        function segmentNormal(obj, frameList, dispResultFlag)
            for i = frameList
                obj.frames(i).segmentNormal(dispResultFlag);
            end
        end
        % -----------------------------------------------------------------
        function segmentTCFA(obj, frameList, dispResultFlag)
            for i = frameList
                obj.frames(i).segmentTCFA(dispResultFlag);
            end
        end
        % -----------------------------------------------------------------
        function dispLmode(obj) % source: C:\Users\Erez\Documents\Projects\Old OCTivat Matlab\Source Code\plaque_detection.m start row 
            %longViewFile=[FileName(1:end-4) ' segmentation results\longiView.mat'];
            xyImg.width = 512;
            M = xyImg.width;
            num_rows = obj.fileInfo(1).Height;
            N = obj.fileInfo(1).Height;
            angle_index = 359;
            
            
            fileHandle.totalNumFrames = size(obj.fileInfo, 1);
            fileHandle.num_rows = obj.fileInfo(1).Height;
            fileHandle.num_cols = obj.fileInfo(1).Width;
            fileHandle.start_frame = 1;
            fileHandle.ending_frame = size(obj.fileInfo, 1);
            fileHandle.frame_index=1;
            [fileHandle.PathName, fileHandle.FileName, ext] = fileparts(obj.pbfn);
            fileHandle.FileName = ['\' fileHandle.FileName, ext];
            %fileHandle.PathName = pathToFile;
            %fileHandle.FileName = fileName;
            
            row_longi_first=ceil(M*angle_index/360);
            row_longi_first = min(row_longi_first, N);
            row_longi_first = max(1,row_longi_first);
            [im_longi_first,slider_im_first,~]=extract_img_longitudinal(row_longi_first,fileHandle);
            
            slider_im_first = slider_im_first(:,:,1);
            im_longi_first_display = imresize(slider_im_first,[size(slider_im_first,1),30*size(slider_im_first,2)]);
            figure('Position', [100, 100, 540, 10]);
            
            imshow(im_longi_first_display);
            axis 'off';
            colormap(gray);

%             if ~exist(longViewFile,'file')
%                 row_longi=ceil(M*angle_index/360);
%                 row_longi=min(row_longi,num_rows);row_longi=max(1,row_longi);
%                 [im_longi,slider_im,slider_im_color]=extract_img_longitudinal(row_longi,handles);
%                 longiView=cat(3,im_longi,slider_im,slider_im_color);
%                 save(longViewFile,'longiView');
%             else
%                 longiView=load(longViewFile);
%                 longiView=longiView.longiView;
%                 im_longi=longiView(:,:,1);
%                 slider_im=longiView(:,:,2:4);
%                 slider_im_color=longiView(:,:,5:end);
%             end
%             im_longi_color=goldenmap(im_longi);
%             if mydata.ColorMap==1
%                 imagesc(slider_im_color,'Parent',handles.axes2);set(handles.axes2,'YTick',[]);
%             else
%                 imagesc(slider_im,'Parent',handles.axes2);colormap(gray);set(handles.axes2,'YTick',[]);
%             end
%             set(handles.axes2,'XTick',mydata.longi_axes);drawnow expose
 
        end
        % -----------------------------------------------------------------
        % inputs: frameNum: I use it for debugging and development: if
        % frameNum is a number, then it corrects only that frame number.  if
        % frameNum='all" I correct the whole pullback
        function correctForCatheter(obj, frameList, dispResultFlag)
            
            for i = frameList        
                obj.frames(i).correctForCatheter(dispResultFlag);
            end
        end
        % -----------------------------------------------------------------
        % inputs: 
        % filterType: either Enhanced Lee of the one published by Schmitt
        % frameNum: I use it for debugging and development: if
        % frameNum is a number, then it correct only that frame number.  if
        % frameNum='all" I correct the whole pullback
        function filteredPullback = removeSpeckle(obj, filterType, frameList, pbfn, dispResultFlag)
            
            filteredPullback = [];
            
            for i = frameList
                    obj.frames(i).removeSpeckle(filterType, pbfn, dispResultFlag);
            end
        end
        % -----------------------------------------------------------------
        % inputs: 
        % filterType: either Enhanced Lee of the one published by Schmitt
        % frameNum: I use it for debugging and development: if
        % frameNum is a number, then it correct only that frame number.  if
        % frameNum='all" I correct the whole pullback
        function computeLumenBorder(obj, frameList, dispResultFlag)           
       
            for i = frameList
               obj.frames(i).computeLumenBorder(dispResultFlag);
               
            end              
            
        end
        % -----------------------------------------------------------------
        function improvedComputeLumenBorder(obj, frameList, dispResultFlag)           
            for i = frameList
                %obj.frames(i).improvedComputeLumenBorder(dispResultFlag);
                 obj.frames(i).borderDetectDP(dispResultFlag);
            end       
        end
        % ----------------------------------------------------------------
        function loadExpertMarking(obj, frameList, dispResultFlag)           
            % Now, import the marking te expert did on this pullback for
            % analysis later on
            [path, name, ext] = fileparts(obj.pbfn);
            marking_fn = '.\Data\XOL-036-00M-LAD-PRE ANALYZED.xlsx'; % marking_fn = [path name '.xlsx']; % '.\Data\XOL-036-00M-LAD-PRE ANALYZED.xlsx';
            %excelData = xlsread(marking_fn);
            excelData = readtable(marking_fn);
            excelData = excelData(5:end, :); % skeep the header
            labelColumn = excelData.(5); % extract the 5th column which is where the tissue type is specified
            if strcmp(labelColumn{1}, 'LABEL')
                for i=2:size(labelColumn,1)
                    disp(labelColumn{i});
                end                
            end
        end            
        % -----------------------------------------------------------------
        function removeBaseline(obj, frameList, dispResultFlag)            
            for i = frameList
                obj.frames(i).removeBaseline(dispResultFlag);
            end           
        end
        % -----------------------------------------------------------------
        % inputs:         
        % frameNum: I use it for debugging and development: if
        % frameNum is a number, then it correct only that frame number.  if
        % frameNum=[min frame: max frame] I correct the frame range
        function computeBackBorder(obj, frameList, dispResultFlag)        
            for i = frameList
                obj.frames(i).computeBackBorder(dispResultFlag);                
            end                    
        end
        % ----------------------------------------------------------------------
        function verifyCoords(obj, coords, frameList)
            for i = frameList
                obj.frames(i).verifyCoords(coords);
            end 
        end
        % -----------------------------------------------------------------
        function pixShift(obj, frameList, dispResultFlag)            
            for i = frameList
                I = obj.frames(i).pixShift(dispResultFlag);
            end           
        end        
       
        % -----------------------------------------------------------------
        function unknownDataset = getAttributes(obj, guiOptions, unknownOptions, frameList, dispFlag)
              
            box = guiOptions.getBoxSize();
            [thetaStride, rStride] = guiOptions.getStride();
            
            imageWidth = obj.fileInfo(1).Width;
            imageHeight = obj.fileInfo(1).Height;
            numOfFrames = obj.numOfFrames;            
            
            tic;
            for i = frameList
                [unknownDataset, ~, ~, ~] = obj.computeAttributes(i, box, thetaStride, rStride, dispFlag);
            end
            fprintf('%s %3.3f %s\n', 'Unknown feature calculations took', toc, ' seconds');
            saveToFile_flag = true; % maybe add this to the gui options
            if saveToFile_flag
                csvwrite(unknownOptions.dataset_fn, unknownDataset);
            end           
        end
        % -----------------------------------------------------------------
        % I assume here that the attributes file exists and that the first
        % 3 columns are the coordinates
       
        function [labelImg, Ytest] = readValidationLabels(obj, ui, frameList, fns, false)
                   
             for i=1:numel(frameList)             
                att_fn = fns.attributes{i};
                imglabel_fn = fns.trueLabels{i};
                validationDataset = load(att_fn); % read the coordinates from the attribute file
                coords = validationDataset(:, [1 2 3]); % (r, theta, frameNumber)
                [labelImg, Ytest] = obj.frames(frameList(i)).readValidationLabels(imglabel_fn, coords);                
                
                csvwrite(fns.labels_fn{i}, Ytest);
                fprintf('Labels dadaset written to : %s\n', fns.labels_fn{i});
             end            
        end
        % -----------------------------------------------------------------
        function [Xnew] = computeAttributes(obj, ui, frameList, fns, dispFlag)        
           
            Xnew = [];
          
            box = ui.getBoxSize();
            [thetaStride, rStride] = ui.getStride();            
            frameIndexInBuffer = 1;
            for fNumber = frameList                  
                frame = obj.frames(fNumber);
                [imageW, imageH] = frame.getImgSize();
                %frame.writeFrameToDisk('MLframe97.tif');
                lumenVector = frame.getLumenBorder(); % get the location of the border (arranged in a vector, each element represnts the r value and the entry line number is the theta index)
                backBorderVector = frame.getBackBorder();
               
                % to use the lumen and backborder in VC, I store them is a
                % text file and then read them from the VC implementation
                               
                %csvwrite([getenv('OCTDATA_DIR') 'ValidationData\Prabu\Vessel 64\Borders\' 'lumen_' num2str(fNumber) '.txt'], lumenVector-1); 
                %csvwrite([getenv('OCTDATA_DIR') 'ValidationData\Prabu\Vessel 64\Borders\' 'backB_' num2str(fNumber) '.txt'], backBorderVector-1);fprintf('wrote %d\n', fNumber);
                frameImg = frame.getData();
                if dispFlag
                    
                    frame.displayImage(frameImg, 'rt', 512);hold on;
                    title('Original r-theta view');
                    y_lumen = 1:numel(lumenVector); % Globals.NUM_LINES_PER_FRAME;
                    plot(lumenVector, y_lumen, 'r-'); hold on;
                    %frame.computeBackBorder(false);
                    plot(backBorderVector, y_lumen, 'y-'); hold on;
                    % There is a rectangle() command I should use:rec = rectangle('Position', [box.origin(1) box.origin(2) box.width box.height], 'EdgeColor', 'r')
                    %                 rec = plot([box.origin(2) box.width box.width box.origin(2) box.origin(2)],...
                    %                     [box.height box.height box.origin(1) box.origin(1) box.height], 'r-');
                end
                
                h = waitbar(0,'Computing validation attributes...');
                
                %best_mutMap= zeros(size(lumenVector,1)-box.height, obj.fileInfo(1).Width-box.width);
                %avg_mutMap= zeros(size(lumenVector,1)-box.height, obj.fileInfo(1).Width-box.width);
                med_mutMap= zeros(size(lumenVector,1)-box.height, obj.fileInfo(1).Width-box.width);
                
                %% Now, compute the attributes
                %for theta=1+floor(box.height/2):thetaStride:size(lumenVector,1)-floor(box.height/2) % imageHeight
                
                pixelCounter = 0;
                for theta=1:thetaStride:size(lumenVector,1)
                    waitbar(theta/size(lumenVector,1), h, 'Computing unknown data...');                    
                    for r=lumenVector(theta): rStride:backBorderVector(theta)                   
                        pixelCounter = pixelCounter + 1;
                        box.origin = [theta; r; fNumber]; % = [theta index, r index, frame number]
                        minTheta = max(theta-floor(box.height/2),1);
                        maxTheta = min(theta+floor(box.height/2), imageH);
                        minR = max(r-floor(box.width/2),1);
                        maxR = min(r+floor(box.width/2), imageW);
                        if dispFlag
                            %mBox.paintMbox(pbfn);drawnow;
                            %delete(rec); rec = plot([r r+box.width r+box.width r r], [theta+box.height theta+box.height theta theta theta+box.height], 'r-');
                            h = patch([minR maxR maxR minR], [minTheta minTheta maxTheta maxTheta], 'r');
                            set(h,'FaceAlpha',0.25);               
                        end                                          
                        boxImg = frameImg(minTheta:maxTheta, minR:maxR);
                        mBox = obj.createMBox(box, size(lumenVector,1));
                                              
                        [mut, logInot, Ihat] = mBox.computeOpticalAttributes('LSQM');                                              
                        %best_mutMap(theta, r) = mut.best;                       
                        %avg_mutMap(theta, r) = mut.avg;
                        med_mutMap(theta, r) = mut.med;                                
                        [dLumen, dCath, dDepth] = mBox.computeGeometricAttributes(r, lumenVector(theta), backBorderVector(theta));    
                        if bitget(ui.getPreProcessingFlags(), Globals.FLATTEN_BIT)
                            mBox.flattenData();
                        end
                        texture = mBox.computeTextureAttributes('mine', nan);
                        glcm = mBox.computeGLCMattributes(boxImg); %Gray Level Co-occurrence Matrix
%                         pullbackFileName = [pullbackFileName; obj.pbfn];
%                         frameNumber = [frameNumber; fNumber];
%                         rt_coordinates = [rt_coordinates;  [theta r]];                        
                       
                        Xnew(end+1, :) = [r, theta, fNumber,...
                                          mut.avg, mut.std, mut.med, mut.se, mut.iqr,...
                                          logInot.avg, logInot.std, logInot.med, logInot.se,...
                                          nan, Ihat,...
                                          dLumen, dCath,dDepth...
                                          texture.smoothness, texture.uniformity,...
                                          texture.homogeneity, texture.contrast, texture.entropy,...
                                          glcm.Contrast, glcm.Correlation, glcm.Energy, glcm.Homogeneity,...
                                          glcm.maxProbability, glcm.entropy];

                    end
                end               
                fprintf('%d pixels in frame %d\n', pixelCounter, fNumber);
                close(h)                
                
%                 saveToFile_flag = true; % maybe add this to the gui options
%                 
%                 if saveToFile_flag                  
                    csvwrite(fns.attributes{frameIndexInBuffer}, Xnew);                   
                    fprintf('Validation dadaset written to : %s\n', fns.attributes{frameIndexInBuffer});
%                end
                frameIndexInBuffer = frameIndexInBuffer + 1;
                Xnew = [];
                
                
                %imwrite(best_mutMap, 'bestmut.tif');
                %imwrite(avg_mutMap, [num2str(ui.getCaseNum()), '_avgmut.tif']);
                %imwrite(med_mutMap, [num2str(ui.getCaseNum()), '_medmut_F' num2str(fNumber) '.tif']);
                %
                %frame.display(best_mutMap, 'both');hold on; title('rmse-based best mut image');
                %frame.display(avg_mutMap, 'both');hold on;title('Average mut image');
                %frame.display(med_mutMap, 'both');hold on;title('Median mut image');
            end          
        end
            function [Xnew] = computeAttributes_old(obj, ui, frameList, fns, dispFlag)%(obj, fNumber, box, thetaStride, rStride, dispFlag)        
           
            Xnew = [];
          
            box = ui.getBoxSize();
            [thetaStride, rStride] = ui.getStride();            
            frameIndexInBuffer = 1;
            for fNumber = frameList                  
                frame = obj.frames(fNumber);
                [imageW, imageH] = frame.getImgSize();
                %frame.writeFrameToDisk('MLframe97.tif');
                lumenVector = frame.getLumenBorder(); % get the location of the border (arranged in a vector, each element represnts the r value and the entry line number is the theta index)
                backBorderVector = frame.getBackBorder();
               
                % to use the lumen and backborder in VC, I store them is a
                % text file and then read them from the VC implementation
                               
                %csvwrite([getenv('OCTDATA_DIR') 'ValidationData\Prabu\Vessel 64\Borders\' 'lumen_' num2str(fNumber) '.txt'], lumenVector-1); 
                %csvwrite([getenv('OCTDATA_DIR') 'ValidationData\Prabu\Vessel 64\Borders\' 'backB_' num2str(fNumber) '.txt'], backBorderVector-1);fprintf('wrote %d\n', fNumber);
                frameImg = frame.getData();
                if dispFlag
                    
                    frame.displayImage(frameImg, 'rt', 512);hold on;
                    title('Original r-theta view');
                    y_lumen = 1:numel(lumenVector); % Globals.NUM_LINES_PER_FRAME;
                    plot(lumenVector, y_lumen, 'r-'); hold on;
                    %frame.computeBackBorder(false);
                    plot(backBorderVector, y_lumen, 'y-'); hold on;
                    % There is a rectangle() command I should use:rec = rectangle('Position', [box.origin(1) box.origin(2) box.width box.height], 'EdgeColor', 'r')
                    %                 rec = plot([box.origin(2) box.width box.width box.origin(2) box.origin(2)],...
                    %                     [box.height box.height box.origin(1) box.origin(1) box.height], 'r-');
                end
                
                h = waitbar(0,'Computing validation attributes...');
                
                %best_mutMap= zeros(size(lumenVector,1)-box.height, obj.fileInfo(1).Width-box.width);
                %avg_mutMap= zeros(size(lumenVector,1)-box.height, obj.fileInfo(1).Width-box.width);
                med_mutMap= zeros(size(lumenVector,1)-box.height, obj.fileInfo(1).Width-box.width);
                
                %% Now, compute the attributes
                %for theta=1+floor(box.height/2):thetaStride:size(lumenVector,1)-floor(box.height/2) % imageHeight
                
                pixelCounter = 0;
                for theta=1:thetaStride:size(lumenVector,1)
                    waitbar(theta/size(lumenVector,1), h, 'Computing unknown data...');                    
                    for r=lumenVector(theta): rStride:backBorderVector(theta)                   
                        pixelCounter = pixelCounter + 1;
                        box.origin = [theta; r; fNumber]; % = [theta index, r index, frame number]
                        minTheta = max(theta-floor(box.height/2),1);
                        maxTheta = min(theta+floor(box.height/2), imageH);
                        minR = max(r-floor(box.width/2),1);
                        maxR = min(r+floor(box.width/2), imageW);
                        if dispFlag
                            %mBox.paintMbox(pbfn);drawnow;
                            %delete(rec); rec = plot([r r+box.width r+box.width r r], [theta+box.height theta+box.height theta theta theta+box.height], 'r-');
                            h = patch([minR maxR maxR minR], [minTheta minTheta maxTheta maxTheta], 'r');
                            set(h,'FaceAlpha',0.25);               
                        end                                          
                        boxImg = frameImg(minTheta:maxTheta, minR:maxR);
                        mBox = obj.createMBox(box, size(lumenVector,1));
                                              
                        [mut, logInot, Ihat] = mBox.computeOpticalAttributes('LSQM');                                              
                        %best_mutMap(theta, r) = mut.best;                       
                        %avg_mutMap(theta, r) = mut.avg;
                        med_mutMap(theta, r) = mut.med;                                
                        [dLumen, dCath, dDepth] = mBox.computeGeometricAttributes(r, lumenVector(theta), backBorderVector(theta));    
                        if bitget(ui.getPreProcessingFlags(), Globals.FLATTEN_BIT)
                            mBox.flattenData();
                        end
                        texture = mBox.computeTextureAttributes('mine', nan);
                        
%                         pullbackFileName = [pullbackFileName; obj.pbfn];
%                         frameNumber = [frameNumber; fNumber];
%                         rt_coordinates = [rt_coordinates;  [theta r]];                        
                       
                        Xnew(end+1, :) = [r, theta, fNumber,...
                                          mut.avg, mut.std, mut.med, mut.se, mut.iqr,...
                                          logInot.avg, logInot.std, logInot.med, logInot.se,...
                                          nan, Ihat,...
                                          dLumen, dCath,dDepth...
                                          texture.smoothness, texture.uniformity, texture.homogeneity, texture.contrast, texture.entropy];
                    end
                end               
                fprintf('%d pixels in frame %d\n', pixelCounter, fNumber);
                close(h)                
                
%                 saveToFile_flag = true; % maybe add this to the gui options
%                 
%                 if saveToFile_flag                  
                    csvwrite(fns.attributes{frameIndexInBuffer}, Xnew);                   
                    fprintf('Validation dadaset written to : %s\n', fns.attributes{frameIndexInBuffer});
%                end
                frameIndexInBuffer = frameIndexInBuffer + 1;
                Xnew = [];
                
                
                %imwrite(best_mutMap, 'bestmut.tif');
                %imwrite(avg_mutMap, [num2str(ui.getCaseNum()), '_avgmut.tif']);
                %imwrite(med_mutMap, [num2str(ui.getCaseNum()), '_medmut_F' num2str(fNumber) '.tif']);
                %
                %frame.display(best_mutMap, 'both');hold on; title('rmse-based best mut image');
                %frame.display(avg_mutMap, 'both');hold on;title('Average mut image');
                %frame.display(med_mutMap, 'both');hold on;title('Median mut image');
            end          
        end
        % -----------------------------------------------------------------
                % Parameters I needed to assume to implement Ughi's approach to
        % automated plaque classification.  In addition, I specify here
        % things he did differently than what I did.
        % 1. mut estimates - he uses something that sounds like RANASC
        % while varying layers, which yields also the number of layers the
        % a-line goes through.  Does not look very reliable since an a-line
        % is very noisy, yert he does not specifiy that he does nay
        % flterring to mut's other than ignoring negative values.
        % 2. what is the Kernel, K ?
        % 2.5 He mentions in sec 3.5 that w_min = 255 um.  Is that the
        % window of the filter? Also tt=25.5um . what is tt?
        % 3. He pre-process the image with low-pass Gaussian (I think he
        % introduced a bias on the mut estimates).  He mentiones the
        % filter's parameters as h_g=7 and sd_g=2.5
        % 4. He assumes d_max to be 1.5mm
        % 5. He does not mention if the ROI is rectangular or free-hand.  I
        % assumed "rectangular"
        %
        % It seems like Ughi computes the slop for each line separately
        % without any regard to a region.  He simply does it within the
        % ROI.  so, to mimick his ROI, I run a window.  His average
        % intensity was computed for a kernel (I assume, the filter's
        % kernel), so his kernel is like an ROI enclosing the pixel of interest 
        function [Xnew] = computeUghiAttributes(obj, ui, frameList, fns, dispFlag)%(obj, fNumber, box, thetaStride, rStride, dispFlag)
           
            Xnew = [];
          
            box = ui.getBoxSize();
            [thetaStride, rStride] = ui.getStride();            
            frameIndexInBuffer = 1;
            for fNumber = frameList                  
                frame = obj.frames(fNumber);
                I = frame.getData();
                lumenVector = frame.getLumenBorder(); % get the location of the border (arranged in a vector, each element represnts the r value and the entry line number is the theta index)
                backBorderVector = lumenVector + 300;               
                
                if dispFlag
                    frame.display(frame.getData(), 'rt');hold on;
                    title('Original r-theta view');
                    y_lumen = 1:numel(lumenVector); % Globals.NUM_LINES_PER_FRAME;
                    plot(lumenVector, y_lumen, '.r'); hold on;
                    %frame.computeBackBorder(false);
                    plot(backBorderVector, y_lumen, '.y'); hold on;
                    % There is a rectangle() command I should use:rec = rectangle('Position', [box.origin(1) box.origin(2) box.width box.height], 'EdgeColor', 'r')
                    %                 rec = plot([box.origin(2) box.width box.width box.origin(2) box.origin(2)],...
                    %                     [box.height box.height box.origin(1) box.origin(1) box.height], 'r-');
                end
                
                h = waitbar(0,'Computing validation - Ughi Method - attributes...');
                
                roiSize = 7; % that i sthe ROI
                kernel = 7;
                sig = 2.5;
                
                pixelCounter = 0;
                for theta=1:round(roiSize/2):size(lumenVector,1)-roiSize % window origin is the upper left hand corner
                    waitbar(theta/size(lumenVector,1), h, 'Computing Ughi''s attributes...');
                    for r=lumenVector(theta): round(roiSize/2):backBorderVector(theta)-roiSize
                        pixelCounter = pixelCounter + 1;
                        box.origin = [theta; r; fNumber]; % = [theta index, r index, frame number]
                        mBox = obj.createMBox(box, size(lumenVector,1));
                        
                        mask = fspecial ('gaussian', kernel, sig);
                        I_filtered = imfilter(I, mask, 'replicate'); % Ughi does not specify, so I use 'replicate' which is the most common
                        winROI = I_filtered(theta:theta+roiSize-1, r:r+roiSize-1);
                        Iavg = mean(winROI(:));
                        [mut, logInot, Ihat] = mBox.computeOpticalAttributes('LSQM');
                        
                        
                        [dLumen, dCath, dDepth] = mBox.computeGeometricAttributes(r, lumenVector(theta), backBorderVector(theta));
                        
                        texture = mBox.computeTextureAttributes('ughi', winROI);
                         
                        texture.diffMoments = nan;
                        texture.shade = nan;
                        
                        Xnew(end+1, :) = [r, theta, fNumber,...
                                          mut.med, Ihat, double(dCath),...
                                          texture.maxProbability, texture.contrast, texture.shade, texture.energy,...
                                          texture.entropy, texture.homogeneity,texture.correlation, texture.diffMoments];
                    end
                end               
                close(h)                                
                saveToFile_flag = true; % maybe add this to the gui options                
                if saveToFile_flag                 
                    csvwrite(fn, Xnew);                   
                    fprintf('Validation dadaset written to : %s\n', fn);
                end
                frameIndexInBuffer = frameIndexInBuffer + 1;
                Xnew = [];
            end          
        end   
        % -----------------------------------------------------------------
        function [X, y] = computeSectorAttributes(obj, ui, frameList, fns)
            
            y = [];
            X = [];
            
            sec = ui.getSectorMeasures();
            tic;
            for fNumber = frameList(1):frameList(2) % frameList 
                frame = obj.frames(fNumber);
                %frame.writeFrameToDisk('MLframe97.tif');
                lumenVector = frame.getLumenBorder(); % get the location of the border (arranged in a vector, each element represnts the r value and the entry line number is the theta index)
                backBorderVector = frame.getBackBorder();   
                %%%%annotatedFrame = frame.getAnnotatedFrame(ui);
              
                h = waitbar(0,'Computing regional attributes...');                
                % Compute frame attributes (note that these are only
                % CURRENT frame attributes, not all frames which the 3D
                % sector spans !!!!!!!)
                lumenShape = frame.computeShapeAttributes();
                [wallLayers, intimaLayerProperties] = frame.computeLayerness();    
                
                for theta=ceil(sec.height/2):sec.stride:obj.fileInfo(1).Height
                    waitbar(theta/size(lumenVector,1), h, 'Computing regional attributes...');
                      
                    % get the annotated label for the sector
                    sector = obj.createSector(fNumber, theta, ui, lumenVector, backBorderVector);
% % % %                     if ~isempty(annotatedFrame)
% % % %                         [y_i, y_i_sequence] = sector.computeAnnotation(theta, obj.fileInfo(1).Height, lumenVector, backBorderVector, annotatedFrame.rt);
% % % %                     end
                    % compute the sectors' attributes
                    texture = sector.computeTextureAttributes();
                    geom = sector.computeGeometricAttributes(theta, obj.fileInfo(1).Height, lumenVector, wallLayers, backBorderVector);
                    
                    % David, this function calculates the edges, their
                    % direction, leading/trailing etc.  It is all of what
                    % you need. You just need to interpret it correctly
                    [edgeSegmentsAttributes, edgeSegmentsImg] = obj.isolateEdges(fNumber, 0);
                    
                    X(end+1, :) = [double(fNumber), double(theta),...
                        lumenShape.eccentricity, lumenShape.area, lumenShape.circularity,...
                        texture.avg, texture.smoothness, texture.uniformity, texture.homogeneity, texture.contrast, texture.entropy,...
                        geom.d_depth, geom.intimalArea, geom.intimalDepth, geom.mediaArea, geom.mediaDepth, geom.adventitiaArea, geom.adventitiaDepth];
                    
                    %%%%y = [y; y_i]; % y is an array of classes
                    
                end
                %%%%uniqueLabels = unique(y);
                close(h);
                fprintf('%s %3.3f %s\n', 'One frame took ', toc, ' seconds');
            end
            %%%%obj.saveFile(X, y, ui);
            %attributeNames = ['FrameNum', 'Theta' Globals.REGIONALATTRIBUTENAMES];
            attributeNames = Globals.REGIONALATTRIBUTENAMES;
            T_x = array2table(X, 'VariableNames', attributeNames);            
            writetable(T_x, 'X.csv');   
            
        end
        % ------------------------------------------------------
        function saveFile(obj, X, y, ui)
            
            fn = obj.createAttributeFilenames(ui);
            
            %attributeNames = {'FrameNum' 'Theta' 'AvgI' 'smoothness' 'uniformity' 'homogeneity' 'contrast' 'entropy' 'BeamDepth' 'initimalArea' 'MediaDepth' 'mediaArea' 'mediaDepth' 'AdArea' 'AdDepth'};
            attributeNames = ['FrameNum', 'Theta' Globals.REGIONALATTRIBUTENAMES];
            T_x = array2table(X, 'VariableNames', attributeNames);            
            writetable(T_x, fn.X);            
           
            T_y = table(y, 'VariableNames', {'label'});
            writetable(T_y, fn.y);
        end
        % -----------------------------------------------------------------
        function fn = createAttributeFilenames(obj, ui)
            tr = ui.getTrainingOptions();
            [~, name, ~] = fileparts(obj.pbfn);
            
            fn.X = [tr.folder name '_X' tr.ext];
            fn.y =  [tr.folder name '_Y' tr.ext];             
        end           
        % -----------------------------------------------------------------
        function viewExpertAnnotations(obj, frameList, fns)            
            trueLabels_fns = fns.trueLabels;
            offset = fns.offset;
            
            for i = 1: numel(frameList)
                frameNum = frameList(i);                
                frame = obj.frames(frameNum);
                trueLables_fn = trueLabels_fns{i}; % trueLabels_fns{frameNum-offset+1}; % case 47 is special and I use trueLabels_fns{i};
                [trueLabelsImg, ~] = frame.convertLabelImage(trueLables_fn, 'rt');
                saveCartesianImage = 0; % for use with Amira to generate 3D of the vessel
%                 if saveCartesianImage
%                     dir = ui.get3Ddir(fns.pb);
%                     frame.saveAsTiff(dir);
%                 end
            end
        end
         % ------------------------------------------------------
        % "loc" is the theta value of the center line of the sector
        function sector = createSector(obj, fNumber, theta, ui, lumenVector, backBorderVector)
            
            secMeas = ui.getSectorMeasures();
            secMeas.origin =  [theta; 1; fNumber]; % [theta, r, frame]
            sector = Sector(secMeas);
            lineCollection = [];
            
            for fn = fNumber-floor(secMeas.depth/2):fNumber+floor(secMeas.depth/2)
                f = obj.frames(fn);
                lowTheta = theta-floor(secMeas.height/2);
                hiTheta = theta+floor(secMeas.height/2);
                
                imgWidth = obj.fileInfo(1).Width;
                for t=max(1, lowTheta):min(hiTheta, obj.fileInfo(1).Height)
                    lineDepth = backBorderVector(t)-lumenVector(t);                    
                    oneLine = f.getSingleLineIntensities(t, lumenVector(t), lineDepth, imgWidth);
                    lineCollection = [lineCollection oneLine];
                end                
            end            
            sector.setLineAggregation(lineCollection);
        end
     end
end
