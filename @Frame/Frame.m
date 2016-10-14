classdef Frame < handle
    
    properties  (GetAccess=private, SetAccess=private)
        index;  % index of the frame in the pullback;
        ImgWidth;
        ImgHeight;
        fn;     % pullback file name this frame belongs to
        data;   % actual frame values - changes as I process the frame's data
        originalData; %I keep a copy of the original, unprocessed data in case I need to display it
        gwMask; % mask of the guidewire
        guidewire; % the guidewire is lines 1-guidewire(1) and guidewire(2)-imageHeight
        lumenBorder;
        iemBorder; %Intima-Media border
        eemBorder; % Media-Adventitia border
        adventitiaBorder;
        backBorder;
        aLineLabel; % a vector: each row represents the theta and the value is the label of the complete A-line
        aLineLabelProb; % a vector: each row represents the theta and the value is the probability that the a-line is of label in prediyedLabel
        textonResultImgMask; % a mask image the same size as "data" - result of classification using textons
        
        % Variables for the gui display
        fig;
        displayAxis;
        movingWindowSpecPanel;
        windowWidth;
        wText;
        aspectRatio;
        arText;        
        strideValue;
        sText;
        updatePB;
        movingWindowOptions;
       
        recHdl; % handle for the window

        
    end
    % =====================================================================
    methods
        function obj = Frame(idx,w, h, fn) % Constructor
            if nargin > 0            
               
                obj.index = idx;
                obj.fn = fn;
                obj.guidewire = [];
                obj.ImgWidth = w;
                obj.ImgHeight = h;
                obj.data = imread(fn, idx);   % 16 bit data

                obj.originalData = obj.data;                
                obj.backBorder = nan(h,1);
                obj.eemBorder = zeros(h,1);
                % set initial value for the moving window
                obj.movingWindowOptions.width = 10;
                obj.movingWindowOptions.aspectRatio = 1;
                obj.movingWindowOptions.stride = 1;
                obj.movingWindowOptions.origin = [0; 0];
                
            end
        end % Constructor
        % -----------------------------------------------------------------
        function [w,h] = getImgSize(obj)
            w = obj.ImgWidth;
            h = obj.ImgHeight;
        end
        % -----------------------------------------------------------------
        function setClassificationImage(obj, classificationImg)
            obj.textonResultImgMask = classificationImg;
        end       
        % -----------------------------------------------------------------
        function windowLevel(obj, I, brightness, contrast)
            % do automatic window level (assume maximum contrast of 100) on
            % the image I, which is currently displayed
            %contrast = 80;
            %brightness = 90;
            maxImage = max(max(I));
            minImage = min(min(I));
            level = (1-brightness/100)*(maxImage - minImage) + minImage;
            width = (1-contrast/100)*(maxImage - minImage);
            
            caxis([level-(width/2) level+(width/2)])
 
        end
        % ------------------------------------------------------        
        function longi = scaleImageData(obj)
            
            Img=log10(double(obj.data+1));
            Img=(Img-min(Img(:)))./(max(Img(:))-min(Img(:)));
            longi = sum(Img,2);
        end        
        % ------------------------------------------------------
        % David please note: The circularity is defined as 
        % "compactness" (see Sonka, 3rd edition, p. 382 and p. 356)
        % = (border length^2/area) where
        % most compact region is a circle.  NOTE THAT THIS IS AN
        % ATTRIBUTE OF THE COMPLETE FRAME. NOT A SECTOR
        function shape = computeShapeAttributes(obj)
            % First, create a binary image of lumen and everything else
            bw = ones(size((obj.data)));
            
            for i=1:numel(obj.lumenBorder)
                bw(i, obj.lumenBorder(i):end) = 0;
            end
            % Convert the resulting image to Cartesian coordinates
            I_xy = rectopolar_fast(im2double(bw)', 512);
            
            % Now, I can compute the statistics I need from it
            cc = bwconncomp(logical(I_xy), 8);
            measurements  = regionprops(cc,'centroid', 'Eccentricity', 'Perimeter', 'Area');
            largestIdx = 1;
            if numel(measurements)>1
                largestIdx = 0;
                a = 0;
                for i=1:numel(measurements)
                    if measurements(i).Area > a
                        largestIdx = i;
                        a = measurements(i).Area;
                    end
                end
            end
            shape.centroid = measurements(largestIdx).Centroid;
            shape.eccentricity = measurements(largestIdx).Eccentricity;
            shape.area = measurements(largestIdx).Area; % * Globals.RTHETA_PIXEL_SIZE^2; % results in microns^2
            
            % Calculate the "circularity" (AKA "compactness")
            p = measurements(largestIdx).Perimeter;
            area = measurements(largestIdx).Area;
            if area > 0
                shape.circularity = p^2/area;
            else
                fprintf('Problem: zero lumen area in frame number %d', obj.index)
            end
            
            dispFlag = 0;
            if dispFlag
                imshow(I_xy)
                hold on
                plot(shape.centroid(:,1), shape.centroid(:,2), 'b*')
                hold off
            end
        end
        % ------------------------------------------------------
         function [layers, intimaLayerProperties] = computeLayerness(obj)
           
             adventitia = obj.adventitiaEdgeDetect(0);
             [iem, intimaLayerProperties, intimaMask] = obj.segmentTCFA(0);
             eem = obj.eemEdgeDetect(0);
             
             layers.iem = iem; % I find the iem in the TCFA segmentation
             layers.eem = eem;
             layers.adventitia = adventitia;             
             
             dispResultFlag = 0;
             if dispResultFlag
                thickness = 3;
                obj.displayImage(obj.originalData, 'rt'); hold on;                
                %plot(obj.lumenBorder, [1:numel(obj.lumenBorder)], 'r-', 'LineWidth', thickness);
                plot(layers.iem+10, [1:numel(layers.iem)], 'g-', 'LineWidth', thickness);
                plot(layers.eem, [1:numel(layers.eem)], 'b-', 'LineWidth', thickness);
                plot(layers.adventitia, [1:size(layers.adventitia,1)], 'm-', 'LineWidth', thickness);
                plot(obj.backBorder, [1:size(obj.backBorder,1)], 'y-', 'LineWidth', thickness);
             end             
         end
        % -----------------------------------------------------------------        
        % This is a copy of getIndividualImages()'s first part
        %
        function binaryImages = separateClassificationMasks(obj, labelImg, coords, predicted, showResults)
           
            % Separate the predicted labels into separate binary images
            
            [m, n] = size(labelImg.rt);
            binaryImages = zeros(m, n, 4);
            for j = 1:size(coords,1)
                r = coords(j, 1);
                theta = coords(j, 2);
                switch predicted.labels(j)
                    case Globals.CALCIUM
                        binaryImages(theta, r, Globals.CALCIUM) = 1;
                    case Globals.LIPID
                        binaryImages(theta, r, Globals.LIPID) = 1;
                    case Globals.FIBER
                        binaryImages(theta, r, Globals.FIBER) = 1;
                    case Globals.BKGD
                        binaryImages(theta, r, Globals.BKGD) = 1;                            
                    otherwise
                        binaryImages(theta, r, Globals.OTHER) = 1;                    
                end                
            end
            
           
            % Note that inspite of the fact that my labelSet has only 3
            % labels, the end results has 5 labels, so I need to show them
            % all
            if strcmp(showResults, 'xy')
                %% need tgo do it
            elseif strcmp(showResults, 'rt')
                actualLabelSet = unique(predicted.labels);
                % First, display the expert-annotated label image
                figure;
                img_handle = imshow(labelImg.rt,[]);
                colormap(Globals.cmMatrix);
                set(img_handle,'CDataMapping','direct');
                % Now, show the separated binary images
                for i=actualLabelSet' % 1:numel(labelSet)
                    %plaqueNum = actualLabelSet(i);
                    figure;
                    fig_handle = imshow(binaryImages(:,:,i),[]);
                    binaryColorMap = [0 0 0; Globals.RGBCOLORS(i,:)];
                    colormap(binaryColorMap);
                    title(Globals.plaqueLabels{i});
                    drawnow;
                end
            end
        end
        % -----------------------------------------------------------------
        function rv = fun(obj, x)
            rv = mode(x(:));
        end
        % -----------------------------------------------------------------
        
        % -----------------------------------------------------------------
        % I should remember (prety much what I do below using the flag
        % combinedBinaryImageProcessing): best results are by following 
        % the following: first, seperate to binary images (each 
        % is the results of one label)
        % second, remove all small classifications (smaller than
        % user-defined size) who do not have neighboring voxels (i.e. run a
        % majority filter and remove all those sreas whose neighberhood is
        % background
        % third, combine the binary images to a single classification image
        % (i.e it will be 4 colors)
        % run a closing operation on the combined image to remove holes and fjords
        %
        function [y_predicted_revised, binaryImages] = cleanClassificationRestuls(obj, coords, binaryImages, predicted)
            % Apply a morphological operations on the fiber binary image
            % "openning" is erosion followed by dialation - removes peninsulas and islands,
            % "closing" is dilation followed by erosion -  removes holes and fjords)
            %
            % TO DO: NEED TO TRY MEDIAN FILTER ON THE COMBINED IMAGES (I.E.
            % NOT JUST THE BINARY IMAGE OF EACH PLAQUE TYPE. THIS WAY
            % SINGLE PIXELS OF ONE TYPE IN A NEIGHBERHOOD OF ANOTER TYPE 
            % WOULD BE SOMEHOW SMOOTHED.  IT IS EASIER TO EXPLAIN IN THE 
            % TEXT 
             
            
            
            % First, make sure everything behind the guidewire, inside the
            % lumen or behind the back-border is eliminated. 
            [~, mask] = obj.extractBloodVesselMask(); % since I called, removeGuideWire(), the I here is already the image w/o guidewire and without data beyond borders before
            %figure, imshow(log(double(I+1)),[])
            for k=1:size(binaryImages,3)
                binaryImages(:,:,k) = uint16(binaryImages(:,:,k)) .* mask;
            end

                
            combinedBinaryImageProcessing = 1;
            if combinedBinaryImageProcessing
                a = rand(5);
                nhood = logical((a+a')/2>.5);
                se = strel('arbitrary', nhood);                
                %fun = @(x) mode(x(:));
                predictedImag_final = zeros(size(binaryImages,1), size(binaryImages,2));
                predictedImag_final = double(predictedImag_final);
                for k=unique(predicted.labels)'
                    % Clean the binary image
                    bw = imerode(binaryImages(:,:,k), se);
                    %bw = imdilate(bw, se);
                    %bw =  imdilate(bw, se);
                    %bw = imerode(bw, se);
                    binaryImages(:,:,k) = bw;
                    % combine the cleaned binary image to the combined image
                    g = binaryImages(:,:,k);
                    gMask = logical(g);
                    predictedImag_final(gMask) = k;
                end
                %ni = nlfilter(predictedImag_final, [9, 5], fun); 
                ni = nlfilter(predictedImag_final, [9, 5], @obj.fun); 
                ni(ni==0) = Globals.OTHER;
                % update the classification results
                y_predicted_revised.labels = zeros(size(predicted.labels,1),1);
                
                for j = 1:size(coords,1)
                    r = coords(j, 1);
                    theta = coords(j, 2);
                    y_predicted_revised.labels(j) = ni(theta, r);
                end
%                 figure;
%                 img_handle = imshow(floor(ni),[]);
%                 colormap(Globals.cmMatrix);
%                 set(img_handle,'CDataMapping','direct');
            else % else work on the individual images
                cleanMethod = 'morphology'; % 'morphology' 'byArea'
                
                if strcmp(cleanMethod, 'byArea')
                    subMethod = 'byLargestArea'; % 'byNumOfPixels' 'byLargestArea'
                    for k=1:size(binaryImages,3)
                        g = logical(binaryImages(:,:,k));
                        if strcmp(subMethod, 'byNumOfPixels')
                            conn = 8;
                            CC = bwconncomp(g,conn); % after Connected Components here I can do area analysis etc...
                            S = regionprops(CC,'Area');
                            areaInPixels = cellfun(@numel,CC.PixelIdxList);
                            toRemove = find(areaInPixels < 40);
                            for i = toRemove
                                g(CC.PixelIdxList{i})=0;
                            end
                            binaryImages(:,:,k) = +g; % the "+" converts the logical to double
                            %figure, imshow(g,[]);
                        else
                            % retain 5 largest components
                            if unique(g)>0
                                g2 = bwareafilt(g, 4, 'largest');
                                binaryImages(:,:,k) = +g2; % the "+" converts the logical to double
                            end
                            %imshowpair(g,g2,'montage')
                        end
                        
                    end
                else
                    %se = strel('disk', 3, 0);
                    actualLabelSet = unique(predicted.labels);
                    for i=actualLabelSet'
                        binaryColorMap = [0 0 0; Globals.RGBCOLORS(i,:)];
                        operation = 'openThenClose'; %'open' 'close' 'openThenClose'
                        %
                        % +++ checkout: http://www.peterkovesi.com/matlabfns/
                        % he seem to have better morphological cleaning !!!
                        %[seg, Am, mask] = tcUtils.mcleanupregions(binaryImages(:,:,i), 0);
                        %
                        % bw2 = bwmorph(binaryImages(:,:,i), operation) ;
                        % figure, imshow(bw2);
                        % colormap(binaryColorMap);
                        %                 operation = 'close';
                        %                 bw = bwmorph(binaryImages(:,:,i), operation) ;
                        %                 figure, imshow(bw);
                        %                 colormap(binaryColorMap);
                        a = rand(5);
                        nhood = logical((a+a')/2>.5);
                        se = strel('arbitrary', nhood);
                        if strcmp(operation, 'open')
                            bw = imerode(binaryImages(:,:,i), se);
                            bw2 = imdilate(bw, se);
                        elseif  strcmp(operation, 'close')
                            bw =  imdilate(binaryImages(:,:,i), se);
                            bw2 = imerode(bw, se);
                        elseif  strcmp(operation, 'openThenClose')
                            bw = imerode(binaryImages(:,:,i), se);
                            bw2 = imdilate(bw, se);
                            bw3 =  imdilate(bw2, se);
                            bw4 = imerode(bw3, se);
                            bw2 = bw4;
                        end
                        %figure, imshow(bw2);
                        %colormap(binaryColorMap);
                        %drawnow;
                        binaryImages(:,:,i) = bw2;
                    end
                end
                %% now update the acuracy numbers to reflect the improved images
                y_predicted_revised.labels = zeros(size(predicted.labels,1),1);
                
                for j = 1:size(coords,1)
                    r = coords(j, 1);
                    theta = coords(j, 2);
                    if binaryImages(theta, r, Globals.CALCIUM) == 1
                        y_predicted_revised.labels(j) = Globals.CALCIUM;
                    elseif binaryImages(theta, r, Globals.LIPID) == 1
                        y_predicted_revised.labels(j) = Globals.LIPID;
                    elseif binaryImages(theta, r, Globals.FIBER) == 1
                        y_predicted_revised.labels(j) = Globals.FIBER;
                    elseif binaryImages(theta, r, Globals.OTHER) == 1
                        y_predicted_revised.labels(j) = Globals.OTHER;
                    else
                        y_predicted_revised.labels(j) = Globals.OTHER;
                    end
                end
                
            end
            y_predicted_revised.probs = predicted.probs;
         end
        % -----------------------------------------------------------------
        function animateTextonScanning(obj, trueLabelsImg)
            
            animateView = 'rtheta'; % 'xy' 'rtheta'
            xyImageSize = 512;
              
            pause on;
            aniFN = 'C:\Users\rys.ADS\Documents\PhD\Dissertation\Presentation\algoAnimationNoPatch.gif'; % animated gif
            createAnimatedGif = 1;
            siSize.w = 50; % size in r-theta view in pixels
            siSize.h = 50;
            lumenVector = obj.getLumenBorder(); % I cal the function since in there I check if the lumen exist
            backBorderVector = obj.getBackBorder();
                       
            %% extract the image data (r-theta) which is only between the borders and
            % 0 elsewhere
            [Imasked, mask] = obj.extractBloodVesselMask();  
            if strcmp(animateView, 'rtheta')
                %obj.displayImage(trueLablesImg.rt, 'rt');hold on;
                obj.displayImage(Imasked, 'rt', xyImageSize);hold on;                
                y_lumen = 1:numel(lumenVector); % Globals.NUM_LINES_PER_FRAME;
                plot(lumenVector, y_lumen, 'r-'); hold on;
                %frame.computeBackBorder(false);
                plot(backBorderVector, y_lumen, 'y-'); hold on;
             elseif strcmp(animateView, 'xy')              
                obj.displayImage(Imasked, 'xy', xyImageSize);hold on;
                %I_xy = rectopolar_fast(im2double(Imasked)',xyImageSize);                       
                % Plot the back border
                x=zeros(size(backBorderVector,1), 1);
                y=zeros(size(backBorderVector,1), 1);
                for i=1:size(backBorderVector,1)
                    xy_pt = obj.convertRTpointToCartesian([obj.backBorder(i); i], 1, xyImageSize) ;        
                    x(i) = xy_pt(1);
                    y(i) = xy_pt(2);
                end                                
                plot(x, y, 'y-', 'MarkerSize', 4);
                
                % plot the lumen
                x=zeros(size(lumenVector,1), 1);
                y=zeros(size(lumenVector,1), 1);
                for i=1:size(lumenVector,1)
                    xy_pt = obj.convertRTpointToCartesian([obj.lumenBorder(i); i], 1, xyImageSize) ;        
                    x(i) = xy_pt(1);
                    y(i) = xy_pt(2);
                end                
                plot(x, y, 'r-', 'LineWidth', 4); hold on;
            end
            options.fillColor = 'y';
            options.drawSectorLines = true;
            options.sectorLineColor = 'r';
            shapeHandle = obj.plotSI([lumenVector(1) 1], siSize, animateView, xyImageSize, options); % 
            
            pixelCounter = 0;
            for theta=1:round(siSize.h/2):size(lumenVector,1)-siSize.w % window origin is the upper left hand corner
                for r=lumenVector(theta): round(siSize.w/2):backBorderVector(theta)-siSize.w
                    if createAnimatedGif
                        f = getframe(1);
                        im = frame2im(f);
                        [imind,cm] = rgb2ind(im,256);
                        if pixelCounter==0
                            imwrite(imind,cm,aniFN,'gif', 'Loopcount',inf);
                        else
                            imwrite(imind,cm,aniFN,'gif','WriteMode','append');
                        end
                    end
                    % STILL NEED TO IGNOR THE GUIDEWIRE - i HAVE NOT
                    % COMPLETED IT YET !
                    pixelCounter = pixelCounter + 1;
                    delete(shapeHandle);                    
                    shapeHandle = obj.plotSI([r theta], siSize, animateView, xyImageSize, options); % 
                    % ensure I am not behind the guidewire
                    %winThetaRange = theta:theta+siSize-1;
                    %winRrange = r:r+siSize-1;
                    %winThetaRange = obj.resizeForGuidewire(winThetaRange);
                    drawnow;
                    pause(0.02);
                end
            end           
        end

        % ------------------------------------------------------
        % Used by texton-based calcium classification
        % assumption: binaryMaskImg is the r-theta classification result (1
        % is calcium, 0 is non-calcium)
        function dispBinaryClassificationImage(obj, Imasked, binaryMaskImg, trueLablesImg, displayMode)            
            binaryColorMap = [0 0 0; Globals.RGBCOLORS(Globals.CALCIUM,:)];
            if strcmp(displayMode, 'xy')
                % image
                xyImageSize = size(trueLablesImg.xy,1);
                figure;
                img_axH = subplot(221);
                I_xy = rectopolar_fast(double(Imasked)', xyImageSize);                    
                imshow(log(I_xy+1),[]);
                hold on; plot(floor(xyImageSize/2), floor(xyImageSize/2), '+y');
                % true labels
                img_axH = subplot(223);                
                img_handle = imshow(trueLablesImg.xy,[]);
                colormap(img_axH, Globals.cmMatrix); 
                set(img_handle,'CDataMapping','direct');
                %hold on; plot(floor(xyImageSize/2), floor(xyImageSize/2), '+y');
                % Calcium segmentation image
                img_axH = subplot(224);
                I_calcium_xy = round(rectopolar_fast(double(binaryMaskImg)', xyImageSize));  
                img_handle = imshow(I_calcium_xy,[]);
                colormap(img_axH, binaryColorMap);                
                %hold on; plot(floor(xyImageSize/2), floor(xyImageSize/2), '+y');                
            elseif strcmp(displayMode, 'rt')
                figure;
                img_axH = subplot(221);
                imshow(log(double(Imasked)+1), []); 
                I_trueLabelsIdx = trueLablesImg.rt==Globals.CALCIUM;                
                img_axH = subplot(223);
                img_handle = imshow(trueLablesImg.rt,[]);
                colormap(img_axH, Globals.cmMatrix); 
                set(img_handle,'CDataMapping','direct');               
                img_axH = subplot(224);
                imshow(binaryMaskImg,[]);               
                colormap(img_axH, binaryColorMap);              
                drawnow;
            end
            
        end
        % ------------------------------------------------------
        % input:
        % ulhc: origin of window in rt view as [r; theta] point
        % siSize: size of subimage side
        % xyImageSize: size of Cartesian image to display
        function h = plotSI(obj, ulhc, siSize, mode, xyImageSize, options)            
            
            if strcmp(mode, 'rtheta')
                x = ulhc(1);
                t = ulhc(2);
                h = plot([x x+siSize.w x+siSize.w x x], [t+siSize.h t+siSize.h t t t+siSize.h], 'y-', 'LineWidth', 3);
            else
                nIncrements = 10;
                delta = siSize.h/nIncrements;
                farArc = zeros(nIncrements,2);
                nearArc = zeros(nIncrements,2);
                
                farArc(1,1) = ulhc(1)+siSize.w; farArc(1,2) = ulhc(2);
                nearArc(1,:) = ulhc;
                
                for k=2:nIncrements
                    farArc(k,:) = [farArc(1,1), farArc(k-1,2)+delta];
                    nearArc(k,:) = [nearArc(1,1), nearArc(k-1,2)+delta];
                end
                
                for k=1:nIncrements
                    pt = obj.convertRTpointToCartesian(farArc(k,:)', 1, xyImageSize);
                    farArc(k,:) = pt(1:2)';
                    pt = obj.convertRTpointToCartesian(nearArc(k,:)', 1, xyImageSize);
                    nearArc(k,:) = pt(1:2)';
                end                
                
                farArc = flipud(farArc);
                X = [nearArc(:,1); farArc(:,1)];
                Y = [nearArc(:,2); farArc(:,2)];
                
                h = patch(X, Y, options.fillColor);
                %set(h,'FaceAlpha',0.25);
             
%                 if options.drawSectorLines
%                     hold on;
%                     plot(options.sectorLineColor);
%                 end
            end
        end
        % ------------------------------------------------------
        function newThetaRange = resizeForGuidewire(obj, winThetaRange)
            newThetaRange = winThetaRange;
            
            if winThetaRange(1) < obj.guidewire(1) && winThetaRange(end) < obj.guidewire(2) && winThetaRange(end) > obj.guidewire(1)
                newThetaRange(winThetaRange>obj.guidewire(1)) = [];
            end
            disp('');
        end
        % -----------------------------------------------------------------
        % zero out all pixels beyond the borders
        function [I, mask] = extractBloodVesselMask(obj)
            
            fun = @(A, B) ~bitxor(A, B);         
            
            lumenVector = obj.getLumenBorder(); % I call the function directly since in there I check if the lumen exist
            backBorderVector = obj.getBackBorder();
            %% extract the image data which is only between the borders and
            % 0 elsewhere
            mask_left = bsxfun(@ge, 1:size(obj.data, 2), lumenVector); % compare row vector 1:size(obj.data, 2) to each element of lumenVector and return 1 when it is greater or equal (ge)
            %mask = uint16(mask);
            %I = obj.data .* mask;
            mask_right = bsxfun(@le, 1:size(obj.data, 2), backBorderVector); % compare row vector 1:size(I, 2) to each element of backborderVector and return 1 when it is less then or equal (le)
            mask = bsxfun(fun, mask_left, mask_right);
            mask = uint16(mask);
            
            % combine with the guidewire mask (must call removeGuideWire()
            % first
            mask(obj.gwMask==0) = 0; % combine with the guidewire mask
            
            I = obj.data .* mask;
        end
        % -----------------------------------------------------------------
        function [assignedLabel, prob] = assignEdgeLabel(obj, edgeMask)
            edgesPresentInWindow = unique(edgeMask);
            
            if numel(edgesPresentInWindow) > 3
                assignedLabel = 1;
                prob = 0.9;
            else
                assignedLabel = -1;
                prob = 0.9;
            end
        end
        
        % -----------------------------------------------------------------
        function [X, coords, y] = loadRegionalValidationDatasets(obj, ui)
           
                % I did it for case 6 just temporarily to make sure everything
                % works
                tr = ui.getTrainingOptions(); % I use the training data til Prabhu gives me more data
                [~, name, ~] = fileparts(obj.fn);
                
                fn.X = [tr.folder name '_X' tr.ext];
                fn.y =  [tr.folder name '_Y' tr.ext];
                
                Xtable = readtable(fn.X);
                X = table2array(Xtable);
                coords = X(:, 1:2); % [frameNum, theta] in the regional dataset
                X = X(:, 3:end);
                
                yTable = readtable(fn.y);
                y = yTable.label;
         
        end
        % ------------------------------------------------------
        % Here I assume the dataset exists.  If not, it is generated at the
        % pullback level.
        function [Xtest, coords, Ytest] = loadValidationDataset(obj, ui, attributes_fn, labels_fn) 
            
           validationDataOptions = ui.getValidationOptions();    
           
           if strcmp(validationDataOptions.source, 'existing')
               if strcmp(ui.getAnalysisType(), 'regional')
                %                    [Xtest, coords, Ytest] = obj.loadRegionalValidationDatasets(ui);
                %                    [Xtest, muTest, sigmaTest] = classifier.featureNormalize(Xtest); % DO NOT FORGET TO NORMALIZE AS DONE IN THE TRAINING AND MODEL CREATION
                %                    labelImg.rt = zeros(obj.ImgHeight, obj.ImgWidth);
                %                    labelImg.xy = zeros(512, 512);
               else % voxel based
                   vdc = ui.getValidationDataLocations();
                   validationDataset = load(attributes_fn);
                   coords = validationDataset(:, [1 2 3]); % (r, theta, frameNumber)
                   Xtest = validationDataset(:, vdc); % [4 5 6];
                   % Now, load the labels
                   if ~isempty(labels_fn)
                       %trueLabelImage_fn = fns.trueLabels{fileNameIdx};
                       %[labelImg, Ytest] = obj.readValidationLabels(ui, trueLabelImage_fn, coords); % the validation values created by the expert
                       Ytest = load(labels_fn);
                   else
                       Ytest = zeros(size(Xtest,1), 1);
                   end
               end
           elseif strcmp(validationDataOptions.source, 'create')
               % Easier to put the computeAttribute() at the Pullback class
               % level (as oppose to the Frame() class level since it
               % involves adjacent frames
               [labelImg, Ytest] = pb.readValidationLabels(ui, frameList, fns, false, '');               
           end            
        end
        % ------------------------------------------------------
        function [convertedLabel_img, labels] = readValidationLabels(obj, fn, coords,  dispResultsFlag)
            
            [convertedLabel_img, originalMembers] = obj.convertLabelImage(fn, dispResultsFlag);            
            
            notMembers = setdiff(unique(convertedLabel_img.rt), originalMembers);            
                        
            labels = zeros(size(coords,1), 1);
            for j = 1:size(coords,1)
                r = coords(j, 1);
                theta = coords(j, 2);                
                labels(j) = convertedLabel_img.rt(theta, r);
            end
            % ==== Verify I did not currupt the data =======
            verify = 0;
            if verify
                [m, n] = size(convertedLabel_img.rt);
                binaryImage = ones(m, n) * Globals.BKGD;
                for j = 1:size(coords,1)
                    r = coords(j, 1);
                    theta = coords(j, 2);
                    switch labels(j)
                        case Globals.CALCIUM
                            binaryImage(theta, r) = Globals.CALCIUM;
                        case Globals.LIPID
                            binaryImage(theta, r) = Globals.LIPID;
                        case Globals.FIBER
                            binaryImage(theta, r) = Globals.FIBER;
                        case Globals.NORMAL
                            binaryImage(theta, r) = Globals.NORMAL;
                        case Globals.OTHER
                            binaryImage(theta, r) = Globals.OTHER;
                        case Globals.BKGD
                            binaryImage(theta, r) = Globals.BKGD;
                    end
                end                
                figure;
                img_handle = imshow(binaryImage,[]);
                colormap(Globals.cmMatrix); set(img_handle,'CDataMapping','direct');
                colorbar('Ticks',[Globals.CALCIUM:Globals.BKGD]+0.5,'YTickLabel', Globals.plaqueLabels);
                title('Converted labels image : r-theta');
                figure;
                img_handle = imshow(convertedLabel_img.rt,[]);
                colormap(Globals.cmMatrix); set(img_handle,'CDataMapping','direct');
                title('Final Label Image (r-\theta) view');
            end            
            % ============= end of verification ============
            
            for k=1:numel(notMembers)
                loc = (labels==notMembers(k));
                labels(loc) = Globals.OTHER; % for now, leave arbitrary assignment
            end
            % unique(labels)
            r = coords(:,1);
            theta = coords(:, 2);
        end
        % -----------------------------------------------------------------
        % Read the label image David Prabhu marked.  The marking is:
        % For frame 267 in pullback "prox LAD.oct":
        % 0 - background
        % 1 - lumen (i.e. blood/flush (contrast) agent)
        % 2 - fibrous
        % 3 - lipid 
        % 4 - Calcium
        % 5 - Normal
        % 6 - Other
        % I do not have the Calcium in this image
        
        function [convertedLabel_img, originalMembers] = convertLabelImage(obj, fn, dispResultsFlag)
                        % label_img is the original XY label image the expert gave me (Dave
            % marked on the XY (of size 1024x1024 but I prefer to work on a 512x512 image) image not the r-t)
            label_img = imread(fn); % Remember, I saved the label images one per frame
            mem = unique(label_img);
            originalMembers = zeros(size(mem));
            originalMembers(mem==0) = Globals.BKGD; % Globals.OTHER; % he marked background as 0
            originalMembers(mem==1) = Globals.BKGD; % Globals.OTHER; % he marked lumen as 1
            originalMembers(mem==2) = Globals.FIBER;
            originalMembers(mem==3) = Globals.LIPID;
            originalMembers(mem==4) = Globals.CALCIUM;
            originalMembers(mem==5) = Globals.BKGD; % Globals.OTHER;
            originalMembers(mem==6) = Globals.OTHER;
            
            
     
            scale = 512/size(label_img, 1);
            % img_xy is the label image after I resized it to 512x512    
            img_xy = imresize(label_img, scale);
            
            % Since the expert gave me labels with numbers which are different
            % from what I use, I created convertedLabel_img to be the lable
            % image after conversion to my label values.
            convertedLabel_img.xy = img_xy; %ones(size(img_xy)) * Globals.BKGD;
            
            convertedLabel_img.xy(img_xy==0) = 0; % Globals.OTHER; % he marked background as 0
            convertedLabel_img.xy(img_xy==1) = 0; % Globals.OTHER; % he marked lumen as 1
            convertedLabel_img.xy(img_xy==2) = Globals.FIBER;
            convertedLabel_img.xy(img_xy==3) = Globals.LIPID;
            convertedLabel_img.xy(img_xy==4) = Globals.CALCIUM;
            convertedLabel_img.xy(img_xy==5) = 0; % Globals.OTHER; % he marked normal as 5 
            convertedLabel_img.xy(img_xy==6) = Globals.OTHER;
            
            % Note, only if I rotate the XY image, I get the r-theta image
            % in the right orientation.  Since I do calculations on the
            % r-theta, I will rotate the XY, but need to keep in mind that
            % for diplay purposes, I should not rotate it since it will
            % apprear 90 degrees rotated relative to the image provided by
            % you (the application of Amr gives the same orientation as
            % yours, so I think that the conversion program we use does
            % some rotation)
            convertedLabel_img.xy = imrotate(convertedLabel_img.xy, -90);
            
            cleanImage = 0;
            if cleanImage
                g = zeros(size(convertedLabel_img.xy));
                g2 = zeros(size(convertedLabel_img.xy));
                for k=1:1 % clean the expert annotation to remove conversion residue (cab repeat several time for better appearance)
                    g(convertedLabel_img.xy==Globals.CALCIUM) = 1;
                    g = medfilt2(g);
                    g2(g==1)= Globals.CALCIUM;
                    g = zeros(size(convertedLabel_img.xy));
                    
                    g(convertedLabel_img.xy==Globals.LIPID) = 1;
                    g = medfilt2(g);
                    g2(g==1)= Globals.LIPID;
                    g = zeros(size(convertedLabel_img.xy));
                    
                    g(convertedLabel_img.xy==Globals.FIBER) = 1;
                    g = medfilt2(g);
                    g2(g==1)= Globals.FIBER;
                    g = zeros(size(convertedLabel_img.xy));
                end
                convertedLabel_img.xy = g2;
            end
            convertedLabel_img.rt = polartorect_fast(double(convertedLabel_img.xy), size(obj.data,2), size(obj.data,1));            
            convertedLabel_img.rt = uint8(round(convertedLabel_img.rt)');
            
            if strcmp(dispResultsFlag, 'xy')
                %% First show the original image
                figure;
                I_xy = log(double(obj.data+1)); % if I first convert and then take the log(Ixy) since it is too dark otherwise)
                I_xy = rectopolar_fast(im2double(I_xy)',512);                
                imshow(I_xy, []);
                title('Original Image (x,y) view');                
                %% Now, show the annotated image
                figure;
                %I1 = imrotate(convertedLabel_img.xy, 90);
                %img_handle = imshow(I1,[]);
                % +++ need to try and work with ind2rgb(X,map) function +++
                img_handle = imshow(convertedLabel_img.xy, []);
                colormap(Globals.cmMatrix);
                set(img_handle,'CDataMapping','direct');
                %colorbar('Ticks',[Globals.CALCIUM:Globals.BKGD],'YTickLabel', Globals.plaqueLabels);
                title('Expert marking labels image : XY');
                drawnow;
                %lcolorbar(Globals.plaqueLabels,'fontweight','bold');
            elseif strcmp(dispResultsFlag, 'rt')
                figure;
                colormap('gray');
                imshow(log(double(obj.data+1)), []);
                title('Original Image (r,t) view');                
                %% Now, show expwert annotation
                figure;
                convertedLabel_img.rt = medfilt2(convertedLabel_img.rt);
                img_handle = imshow(convertedLabel_img.rt,[]);
                colormap(Globals.cmMatrix);
                set(img_handle,'CDataMapping','direct');
                title('Final Label Image (r-\theta) view');
                %[cmin, cmax] = caxis;
                %colorbar('Ticks',[cmin:cmax],'YTickLabel', Globals.plaqueLabels(c1));
            end
        end
        % ------------------------------------------------------
        % in case the guidewire is split at the end of the image,
        % gw(1) stores the first line of the guidewire is 
        % gw(2) is the second line(s.t lines1-gw(1) is the 
        % guidewire and so is and gw(2)-imageHeight  
        function setGuidewireMask(obj, gw)
            
            im1 = im2double(obj.data);       
            obj.guidewire = gw; % this way I will have access to the guidewire location
            mask = ones(size(im1));
           
            if gw(2)-gw(1)<(size(im1,1)/2)
                mask(gw(1):gw(2),:)=0;
            else                
                mask(1:gw(1),:)=0;
                mask(gw(2):obj.ImgHeight,:)=0;
            end
            
            obj.gwMask = uint16(mask);
            
            % update the image data
            obj.data = obj.data .* obj.gwMask;
%             offset = 1;
%             for lineNum=1:obj.guidewire(1)
%                 obj.data(lineNum, :) = obj.data(obj.guidewire(1)+offset, :);
%             end
%             for lineNum=obj.guidewire(2):obj.ImgHeight
%                 obj.data(lineNum, :) = obj.data(obj.guidewire(2)-offset, :);
%             end
            
            
            dispResults = 0;
            if dispResults
                maskedImage = obj.data .* obj.gwMask;
                figure, imshow(obj.gwMask,[]);
                figure, imshow(log10(double(obj.data+1)), []); title('Before');  
                figure, imshow(log10(double(maskedImage+1)), []); title('After');       
            end
        end
        % ------------------------------------------------------
        function updateFrameData(obj, data_in)
            obj.data = data_in;
            obj.ImgHeight = size(data_in,1);
        end
        % -------------------------------------------------------------
        function compareCartesianImages(obj, pbfn)
            if isempty(obj.lumenBorder)
                [~] = obj.computeLumenBorder(false);
            end
            if isnan(obj.backBorder)
                obj.computeBackBorder(false);
            end
            
            [dir, name, ext] = fileparts(pbfn);   
           
            
            size = 512;
            
            sjImg_pbfn = [dir '\' name '_XY_' num2str(size) 'x' num2str(size) '.tif'];  % location of St Jude image
            
            sjImg = imread(sjImg_pbfn, obj.index);
            
            temp = log(double(obj.data+1));
               I_xy = rectopolar_fast(im2double(temp)', size);
               I_xy = imrotate(I_xy,+90);
               grayImg = uint8(255*mat2gray(I_xy));
               
               figure, imshow(grayImg,[]);
               figure, imshow(sjImg,[]);
               
               sjImg(sjImg==255) = 1;
               sjImg(sjImg==0) = 1;
               %% Note that St Jude put 255 values all arounf the perimeter of the image
               quotientImg = grayImg./sjImg;
               figure, imshow(quotientImg,[]);
               
               figure, imhist(quotientImg)
               
        end
        % -------------------------------------------------------------

        function saveAsTiff(obj, root3d_xy)
%             fn_xy = ['Cartesian\' num2str(obj.index) '_704' '.tiff'];
%             fn_rt = ['Polar\' num2str(obj.index) '.tiff'];
%             if ~exist('Cartesian', 'dir')  
%                 mkdir('Cartesian');
%             end
%             if ~exist('Polar', 'dir')  
%                 mkdir('Polar');
%             end
            
           size = 512;
           
           saveTextonResultImage = 0;
           if saveTextonResultImage               
               I_xy = rectopolar_fast(double(obj.textonResultImgMask)', size);
               fn_xy = [root3d_xy 'TexonResult\' num2str(obj.index) '.tiff'];
               imwrite(I_xy, fn_xy); % note that I saving the masked image
               fprintf('Wrote texton result mask [%d]-> %s\n', obj.index, fn_xy);
           end
           
           saveDataFlag = 1;
           if saveDataFlag
               [I, ~] = extractBloodVesselMask(obj);               
               temp = log(double(I+1));
               I_xy = rectopolar_fast(im2double(temp)', size);
               fn_xy = [root3d_xy 'Vessel\' num2str(obj.index) '.tiff'];
               grayImg = uint8(255*mat2gray(I_xy));
               imwrite(grayImg, fn_xy); % note that I saving the masked image
               fprintf('Wrote frame %d -> %s\n', obj.index, fn_xy);
           %imwrite(uint8(uint8(255*mat2gray(obj.data))), fn_rt);
           end
           
           saveLumenFlag = 0;
           if saveLumenFlag
               lumenBinaryImage = obj.computeLumenBorder(false);
               lumenBinaryImage = rectopolar_fast(im2double(lumenBinaryImage)', size);
               lumenBinaryImage = logical(lumenBinaryImage);
               lumen_fn_xy = [root3d_xy '\Lumen\' num2str(obj.index) 'lumen' '.tiff'];
               imwrite(lumenBinaryImage, lumen_fn_xy);
           end
           
           saveBBflag = 0;
           if saveBBflag
               bbBinaryImage = obj.computeBackBorder(false);
               bbBinaryImage = rectopolar_fast(im2double(bbBinaryImage)', size);
               bbBinaryImage = logical(bbBinaryImage);
               bb_fn_xy = [root3d_xy '\Backborder\' num2str(obj.index) 'bb' '.tiff'];
               imwrite(bbBinaryImage, bb_fn_xy);
           end
           
           saveAdventitiaFlag = 0;
           if saveAdventitiaFlag
               obj.adventitiaEdgeDetect(0);
               adBinaryImage = zeros(504,968);
               for i=1:numel(obj.eemBorder)
                   adBinaryImage(i, obj.eemBorder(i)) = 1;
               end
               adBinaryImage = rectopolar_fast(im2double(adBinaryImage)', size);
               adBinaryImage = logical(adBinaryImage);
               ad_fn_xy = [root3d_xy '\Adventitia\' num2str(obj.index) 'ad' '.tiff'];
               imwrite(adBinaryImage, ad_fn_xy);
           end
        end
        % ***********************Post-processing *************************
        % Here I try to find the adventitia layer.  I use the Watershed
        % Transform as a basis because (as mentioned on p. 536):
        % The watershed approach is particularly attractive because:
        % 1. It produces closed, well-defined region boundaries, 
        % 2. It behaves in a global manner, and 
        % 3. It provides a framework in which a priori knowledge can 
        %     be utilized to improve segmentation results
        function border = adventitiaEdgeDetect(obj, dispResultFlag)
            
            if isempty(obj.lumenBorder)
                obj.computeLumenBorder(false);                
            end            
            if isnan(obj.backBorder)
                obj.computeBackBorder(false);                
            end
            
            % Here I use marker-controlled watershed segmentation: see 
            % "Digital Image Processing Using Matlab" book page 593 (watershed example)
            % I consider as "internal" the area to the right of the lumen (i.e. what's
            % behind the media layer)
            
            h = fspecial('sobel');
            % Compute the gradient image 
            g = sqrt(imfilter(double(obj.data), h, 'replicate') .^ 2 + imfilter(double(obj.data), h', 'replicate') .^ 2);
            %rm = imregionalmin(g);
            %rm = imregionalmin(obj.data); % compute the location of all regional minima in an image
            heightThresh = 10;  % set a threshold
            internalMarkers = imextendedmin(obj.data, heightThresh); % im is an image whose foreground pixels mark the locations of the deep regional minima
            % superimpose the data into the original image
            fim = obj.data;
            fim(internalMarkers) = 175;
            Lim = watershed(bwdist(internalMarkers));
            em = (Lim == 0); % em is a logical matrix with zeros and 1's
            %             g2 = imimposemin(g, im | em);  % modify the gradient image using a "minima imposition" (see Soille, 2003)
            %             %g2(g2==-inf) = 0;
            %             L2 = watershed(g2);
            %             f2 = obj.data;
            %             f2(L2 == 0) = 255;
            
            % Now, to extract the back border I scan the rows from the
            % right and mark as border the first slope greater than 1
            lineCounter = 1;
            border = zeros(size(em,1),1);
            for row=1:size(em,1)
                slopes = diff(em(row, :));
                largestSlopeLocation = (slopes == 1);
                ind = find(largestSlopeLocation, 1, 'last'); % by 'last' I in fact, take the last location along the vector with slope==1
                if ~isempty(ind)                    
                    border(lineCounter) = ind;                    
                end
                lineCounter = lineCounter + 1;
            end
            
            obj.adventitiaBorder = obj.correctForGuidewire(border, 10);
            
            
            if dispResultFlag
                % First, display the frame as a surface plot
                x = 1:size(obj.data, 1);
                y = 1:size(obj.data, 2);
                [X,Y] = meshgrid(x,y);
                surf(X,Y,double(obj.data)');
                camlight left;
                lighting phong;
                view([73 58]);
                shading interp;
                %% surface plot end
                obj.displayImage(obj.originalData, 'rt', nan); hold on;
                y_lumen = 1:size(obj.adventitiaBorder,1);
                plot(obj.adventitiaBorder, [1:size(obj.adventitiaBorder,1)], 'm-', 'LineWidth', 3);
                plot(obj.backBorder, [1:size(obj.backBorder,1)], 'y-', 'LineWidth', 3);
                plot(obj.lumenBorder, [1:size(obj.lumenBorder,1)], 'r-', 'LineWidth', 3);
                
                title('(r,\theta) Image with back border overlay');                
                
                obj.displayImage(obj.originalData, 'xy', 704); hold on;
                %[x,y,cp] = obj.convertRTpointToCartesian(obj.lumenBorder, 1, 704)                   
                
                %% Display the XY images
                % Convert lumen
                x=zeros(size(obj.lumenBorder,1), 1);
                y=zeros(size(obj.lumenBorder,1), 1);
                for i=1:size(obj.lumenBorder,1)
                    xy_pt = obj.convertRTpointToCartesian([obj.lumenBorder(i); i], 1, 704) ;        
                    x(i) = xy_pt(1);
                    y(i) = xy_pt(2);
                end
                  centroid.x = sum(x(:))/size(x,1);
                centroid.y = sum(y(:))/size(y,1);           
                plot(x, y, 'r-', 'LineWidth', 3); hold on;
                %plot(centroid.x, centroid.y, '+r', 'MarkerSize', 5);
                
                % Convert the back border
                x=zeros(size(obj.backBorder,1), 1);
                y=zeros(size(obj.backBorder,1), 1);
                for i=1:size(obj.backBorder,1)
                    xy_pt = obj.convertRTpointToCartesian([obj.backBorder(i); i], 1, 704);      
                    x(i) = xy_pt(1);
                    y(i) = xy_pt(2);
                end
                            
                plot(x, y, 'y-', 'LineWidth', 3); hold on;
                
                % Convert Adventitia border
                x=zeros(size(obj.adventitiaBorder,1), 1);
                y=zeros(size(obj.adventitiaBorder,1), 1);
                for i=1:size(obj.adventitiaBorder,1)
                    xy_pt = obj.convertRTpointToCartesian([obj.eemBorder(i); i], 1, 704) ;         
                    x(i) = xy_pt(1);
                    y(i) = xy_pt(2);
                end
                             
                plot(x, y, 'm-', 'LineWidth', 3); hold on;
                
                
            end
        end
        % ------------------------------------------------------
        % I will scan the image between the lumen and the adventitia
        % and use dynamic programing to try and detect a maximum line that
        % will be the eem (external elastic membrane) - the border between
        % the media and adventitia
        function border = eemEdgeDetect(obj, dispResultFlag)
            
            % create a mask of the region between the lumen and the
            % adventitia            
            regionMask  = uint16(zeros(size(obj.data)));
            
%             for t=1:size(obj.data,1)
%                 for r=obj.lumenBorder(t):obj.adventitiaBorder(t)
%                     regionMask(t, r) = 1;
%                 end
%             end

            startColBeyondIEM = 10;
            bandwidth = 5;
            [M, N] = size(obj.data);
            slopeImg = zeros(M, N);
            border = zeros(M,1);
            halfW = 5;
            
            for t = 1:M
                %for r =obj.lumenBorder(t):obj.adventitiaBorder(t)
                for r =obj.iemBorder(t):obj.adventitiaBorder(t)
                    minTheta = max(1, t-halfW);
                    maxTheta = min(t+halfW, M);
                    %slope = mean(mean(obj.data(minTheta:maxTheta, r:min(r+bandwidth,N)),2) - mean(obj.data(minTheta:maxTheta, r-bandwidth:r),2));                  
                   slope = mean(obj.data(t, r:min(r+bandwidth,N)),2) - mean(obj.data(t, max(r-bandwidth,1):r),2);                   
                    slopeImg(t, r) = max(slope, 0);                    
                end
                %startCol = obj.lumenBorder(t)+startColBeyondLumen;
                startCol = obj.iemBorder(t)+startColBeyondIEM;
                [~, border(t)] = max(slopeImg(t, startCol:end), [], 2);
                border(t) = border(t) + startCol;
            end
            
            % create a binary image of the result so I can connect the
            % border
            tempImg = zeros(M,N);
            for j=1:numel(border)
                tempImg(j,border(j))=1;
            end
            
           [tempImg, connectedBorder] = connect_border(tempImg);
           
           obj.eemBorder = connectedBorder;
            
%             catheterPos = 70;
%                 lumen_connectivity_rect = 15;
%                 bandwidth = 60; 
%             borderFrame = lumen_segmentation_DP(f,catheterPos,lumen_connectivity_rect,bandwidth); 
            
            
            if dispResultFlag
                figure, imshow(log(double(slopeImg+1)), []);hold on;
                plot(obj.lumenBorder, [1:M], 'r-', 'LineWidth', 1);
                plot(obj.iemBorder, [1:M], 'g-', 'LineWidth', 1);
                plot(obj.eemBorder, [1:numel(obj.eemBorder)], 'b-', 'LineWidth', 2);
                
            end
            
        end
        % ------------------------------------------------------
        function centroid = computeCentroidLine(obj, dispResultFlag)
            if isempty(obj.lumenBorder)
                obj.computeLumenBorder(false);                
            end   
            
            % Convert lumen
            x=zeros(size(obj.lumenBorder,1), 1);
            y=zeros(size(obj.lumenBorder,1), 1);
            for i=1:size(obj.lumenBorder,1)
                xy_pt = obj.convertRTpointToCartesian([obj.lumenBorder(i); i], 1, 512) ;
                x(i) = xy_pt(1);
                y(i) = xy_pt(2);
            end
            centroid.x = sum(x(:))/size(x,1);
            centroid.y = sum(y(:))/size(y,1);
            centroid.z = obj.index;
                           
            root3d_xy = ['C:\Users\rys.ADS\Documents\Projects\TissueCharacterization\OCT_Data\ValidationData\Prabu\Vessel 52\3D Data\'];
            
            saveCentroidFlag = 1;
            if saveCentroidFlag
                centroidBinaryImage = zeros(512);
                centroidBinaryImage(round(centroid.x), round(centroid.y)) = 1;
                se = strel('disk',10, 0);
                centroidBinaryImage = imdilate(centroidBinaryImage, se);
                centroidBinaryImage = logical(centroidBinaryImage);
                centroid_filename = [root3d_xy '\Centroid\' num2str(obj.index) 'centroid' '.tiff'];
                imwrite(centroidBinaryImage, centroid_filename);
            end
            
     
           
             
            
            if dispResultFlag
                obj.displayImage(obj.data, 'xy', 704); hold on;
                plot(x, y, 'r-', 'LineWidth', 3); hold on;
                plot(centroid.x, centroid.y, '+y', 'MarkerSize', 5);
            end
            
        end
        % ------------------------------------------------------------------
        % Return value: CC
        % CC is initialized using connected components and then I add more
        % fields using regionprops().  The regions of the connected
        % compnents are, in fact, the edges found.
        % Then I add more features such as:
        % CC.Orientation : See regionprops() 
        % 1) CC.Extrema : See regionprops()  
        % 2) CC.PixelValues : See regionprops() 
        % 3) CC.PixelListCartesian : See regionprops() 
        % 4) CC.gradientMagnitude: gradient at the edge pixel - it's an image
        % 5) CC.gradientDir gradient dir at the edge pixel - it's an image
        % 6) Curvature is d(angle)/d(arc) = rate of angle change
        % w.r.t. the arc it spans.  I calculate the curvature
        % at the end points of each edge (i.e. for each edge, I
        % have 2 curvature values
        % 7) Continuity: This will go along with the values of curvature of
        % the edge's end points (see my presentation:
        % C:\Users\rys.ADS\Documents\Projects\Lab Presentations\3D Data Extraction.pptx
        % I am looking to have C0 & C1 curvatues 
        function [CC, edgeMaskImg, g] = isolateEdges(obj, dispResultFlag)
            
            CC = [];
            
            method = 'automated'; % 'automated'; % 'interactive';
            
            if strcmp(method, 'automated')
                
                %% 1) Try Sobel
%                 cathLoc = 60;
%                 T = 0.000883; %0.000883                
%                 [g, t] = edge(obj.data(:,cathLoc:end), 'sobel', T, 'both');
%                 %fprintf('Threshold used = %f\n', t);
%                 figure, imshow(g,[]); hold on;
                %% 2) try entropy of gradients
%                 dx = 10;
%                 dy = 10;
%                 % Looking at the whole image is too dense, so subsample
%                 f = obj.data(1:dy:end,1:dx:end);
%                 
%                 [Gmag,Gdir] = imgradient(f);
%                 x = 1:size(f,2); % 1:obj.ImgWidth;
%                 y = 1:size(f,1); % 1:obj.ImgHeight;
%                 z = f; % obj.data;
%                 figure, surf(x,y,z,Gdir);
%                 colorbar;
%                 figure, surf(x,y,z,Gmag);
%                 colorbar;
%                 e = entropy(f);
%                 figure, imshow(e,[]);
                %% 3) Try Canny Edge detector
                
                [I, bloodVesselMask] = obj.extractBloodVesselMask();
               
                if dispResultFlag
                    figure, imshow(log(double(obj.data)+1),[]);
                    figure, imshow(log(double(I)+1),[]);
                end
                T = []; %[170, 400]; %[] will let the algo compute the T values. T = [T_low, T_hi] where T_low < T_hi < 1.0
                sigma = 10;
                [g, t] = edge(I, 'canny', T, sigma);
                %[gGradient, tGradient] = edge(I, 'sobel', T, 'both'); % sobel returns a binary image w/ the edges at those points where the gradient of I is maximum
                [Gmag,Gdir] = imgradient(I, 'sobel');
                % Since g is a binary image I can use connected components 
                conn = 8;
                CC = bwconncomp(g,conn);
                numEdges = CC.NumObjects;
                % Create a mask with the edge numbers for later use
                edgeMaskImg = zeros(size(g));
                for edgeNum=1:numEdges
                  edgeMaskImg(CC.PixelIdxList{edgeNum}) = edgeNum;  
                end
                if dispResultFlag
                    %figure, imshow(g,[]);                    
                    % keep this code to see how to access the data efficiently                    
                    %% +++(keep to know how to access the data)
                    %numPixels = cellfun(@numel,CC.PixelIdxList);
                    %[biggest,idx] = max(numPixels);
                    %rgbImg = zeros(size(g));
                    %rgbImg(CC.PixelIdxList{idx}) = 255; %  make the largest line white
                    %%rtrgb = ind2rgb(rgbImg, Globals.cmMatrix);
                    %%figure , imshow(rtrgb,[]);
                    %%+++
                    redChannel = zeros(size(g,1), size(g,2), 'uint8');
                    greenChannel = zeros(size(g,1), size(g,2), 'uint8');
                    blueChannel = zeros(size(g,1), size(g,2), 'uint8');
                    imshow(log(double(I)+1),[]); hold on;
                    %% Run through the lines and create an rgb image of the edges
                    %visboundaries(g,'Color','r');
                    for edgeNum=1:numEdges
                        rng; % set a seed to that results are repeatable
                        colorValues = randi(255,3,1);
                        redChannel(CC.PixelIdxList{edgeNum}) = colorValues(1);
                        greenChannel(CC.PixelIdxList{edgeNum}) = colorValues(2);
                        blueChannel(CC.PixelIdxList{edgeNum}) = colorValues(3);
                        coloredImage = cat(3, redChannel, greenChannel, blueChannel);
                        imshow(coloredImage); hold on;
                        %hold on, plot(CC.PixelIdxList{nol}, '.r'); hold on;
                        drawnow;
                    end
                    % Show images side by side
                    figure, imshowpair(I,coloredImage,'montage');
                    title('Images side by side using imshowpair()');
                    % Fuse the images together
                    C = imfuse(I,coloredImage);
                    figure, imshow(C,[])
                    title('Fused images'); % using imfuse()');
                end
                
            else % interactive methods
                edgeMaskImg = [];
                if isempty(obj.lumenBorder)
                    obj.computeLumenBorder(false);
                end
                %% Select the ROI
                imshow(obj.data,[]);
                h = imrect;
                pos = getPosition(h);
                setColor(h,'r');
                woi.origin = [pos(1);pos(2)];
                woi.width = pos(3);
                woi.height = pos(4);
                woi.data = obj.data(woi.origin(2):woi.origin(2)+woi.height, woi.origin(1):woi.origin(1)+woi.width);
                woi.lumen = obj.lumenBorder(woi.origin(2):woi.origin(2)+woi.height);
                woi.lumen = woi.lumen - woi.origin(1);
                
                figure, imshow(woi.data,[]);
                hold on;
                plot(woi.lumen, [1:numel(woi.lumen)], 'r-','LineWidth', 4);
                
                %% Now, analysis is done only on the Window of Interest
                
                %% Try to see the gradient image
                [Gmag,Gdir] = imgradient(woi.data);
                
                dx = 1;
                dy = 1;
                %% Lookimg at the whole image is too dense, so subsample using dx & dy
                f = woi.data(1:dy:end,1:dx:end);
                
                [Gx,Gy] = imgradientxy(f);
                row = 1:size(f, 2);
                X = repmat(row, size(f,1),1 );
                col = [1:size(f, 1)]';
                Y  = repmat(col, 1, size(f,2) );
                
                
                if dispResultFlag
                    figure;
                    imshow(f,[]); hold on;
                    contour(X,Y,f, 'w-');hold on;
                    quiver(X, Y, Gx, Gy, 'y-');
                    
                end
            end
        end
        
        % -----------------------------------------------------------------
        function segmentNormal(obj, dispResultFlag)
            
            w = imread('XOL-002-00M-LAD-PRE_normal_template.tif'); % load the template
            
            g = abs(normxcorr2(w, obj.data));
            figure, imshow(g, []);
            gT = (g== max(g(:))); % gT is a logical array.
            % Find out how many peaks there are.
            idx = find(gT == 1); % We use idx again later.
            numel(idx);
            gT = imdilate(gT, ones(7)); % since one point is hard to see
            figure, imshow(gT, []);
        end
        % *********************** Pre-processing *****************************
        function filteredFrame = removeSpeckle(obj, filterType, pbfn, dispResultFlag)
            
            switch filterType 
                case 'Enhanced Lee'
                    filteredFrame = obj.enhancedLee(pbfn);
                case 'RGMAP' % Refined Gamma Maximum-A-Posteriori (RGMAP) filter
                case 'Frost'
                    %I = double(imread('..\TCFA-S-1.oct', 57));                    
                    filteredFrame = fcnFrostFilter(obj.data, getnhood(strel('disk',3,0)));
                otherwise
                    disp('No Filterring done\n');
            end
            
            if dispResultFlag
                obj.displayImage(obj.data, 'both', 704);
                title('Image - before speckle removal');                
                obj.displayImage(filteredFrame, 'both', 704);
                title('Denoised image - after speckle removal');  
            end
            
            obj.data = filteredFrame; % I am not keeping the original data.  Will have to remove this and return the denoised image to the calling fcn
            
        end       
        % ----------------------------------------------------------------
        % Just extracted from Zaho's code for now.  I need to improve
        % results and incorporate into my code (naming, hard coded values
        % etc............).  After excusion of this function the variable 
        % obj.iemBorder will include the border around the whole intima layer
            % (both sides the lumen and the far side from the lumen).  But
            % the function's return value will be the same format as I do
            % with the lumen and the other borders: each entry of a column
            % border will include the r value of the border of the iem
        function [border, cap, intimaMask] = segmentTCFA(obj, dispResultFlag)
            
            band_width=round(0.5/(4.8/obj.ImgWidth));
            band_width = 20;
            inten_width=35;
            %cordX = [157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219];
            cordX = 1:obj.ImgHeight;
            if isempty(obj.lumenBorder)
                    obj.computeLumenBorder(false);
            end
            lumen_index =   obj.lumenBorder;
            Img_gw_remove = obj.data;
            
            fibrous_cap=TCFA_boundary_extra_signelRegion(Img_gw_remove,cordX,lumen_index,band_width,inten_width);
            % Total area in r-theta view is just the number of pixels 
            % whose value is 1 in the variable fibrous_cap times 
            % the area of a single pixel
            intimaMask = fibrous_cap;
            cap.area = sum(sum(fibrous_cap)); % * Globals.RTHETA_PIXEL_SIZE^2; % result in micros^2
            cap.intimaPerimeter = bwperim(fibrous_cap,8); % assign output value
            obj.iemBorder = obj.extractiem(cap.intimaPerimeter); % update class variable
            obj.iemBorder(obj.iemBorder==0) = 1; % because Iuse the values as index and I cannot have zero index (in eemEdgeDetect())
            border = obj.iemBorder; % assign function output value
             
            if dispResultFlag
                %% show the r-theta results
                %imshow(obj.originalData,[]);
                obj.displayImage(obj.originalData, 'rt', nan); hold on;
                hold on;
                [r1,c1] = find(cap.intimaPerimeter);
                plot(c1,r1,'.y', 'MarkerSize', 0.25);
                plot(obj.iemBorder, [1:numel(obj.iemBorder)], 'w-', 'MarkerSize', 2);
                % Plot the other borders
                 plot(obj.adventitiaBorder, [1:size(obj.adventitiaBorder,1)], 'm-', 'LineWidth', 1);
                 if isempty(obj.lumenBorder)
                     obj.computeLumenBorder(false);
                 end
                 if isnan(obj.backBorder)
                     obj.computeBackBorder(false);
                 end
                 plot(obj.backBorder, [1:size(obj.backBorder,1)], 'y-', 'LineWidth', 1);
                 plot(obj.lumenBorder, [1:size(obj.lumenBorder,1)], 'r-', 'LineWidth', 2);
                 
                 %% Show, XY view
                 [~,imLU_r,imLU_t]=rectopolar_fast(im2double(obj.data)', 512);
                 
                 fibrous_cap_polar=rectopolar_fast_LU(double(fibrous_cap)', imLU_r, imLU_t);
                 fibrous_cap_polar=logical(fibrous_cap_polar);
                 cap.xyPerimeter = bwperim(fibrous_cap_polar,8);
                 boundaryIndex=find(cap.xyPerimeter);
                 [r,c] = find(cap.xyPerimeter);                 
                 
                 im_xy = rectopolar_fast(im2double(obj.data)',512);
                 figure, imshow(log(double(im_xy+1)),[])
                 hold on;
                 plot(c,r,'.y', 'MarkerSize', 0.25);
                 
            end
        end
        % ----------------------------------------------------------------
        % Just extracted from Zaho's code for now.  I need to improve
        % results and incorporate into my code (naming, hard coded values
        % etc............)
        function segmentTCFA_ORIGINAL(obj, dispResultFlag)
            
            band_width=round(0.5/(4.8/obj.ImgWidth)); %=968
            inten_width=35;
            cordX = [157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219];
            
            if isempty(obj.lumenBorder)
                    obj.computeLumenBorder(false);
            end
            lumen_index =   obj.lumenBorder;
            Img_gw_remove = obj.data;
            
            fibrous_cap=TCFA_boundary_extra_signelRegion(Img_gw_remove,cordX,lumen_index,band_width,inten_width);
            % Toral area in r-theta view is just the number of pixels in
            % the variable fibrous_cap times ghe area of a single pixel
            cap.rtBoundary = bwperim(fibrous_cap,8);
            [~,imLU_r,imLU_t]=rectopolar_fast(im2double(obj.data)', 512);           
            
            fibrous_cap_polar=rectopolar_fast_LU(double(fibrous_cap)', imLU_r, imLU_t);
            fibrous_cap_polar=logical(fibrous_cap_polar);
            cap.xyBoundary = bwperim(fibrous_cap_polar,8);
            boundaryIndex=find(cap.xyBoundary);
            [r,c] = find(cap.xyBoundary);
           
            im_xy = rectopolar_fast(im2double(obj.data)',512);
            
            % show the r-theta results
            imshow(obj.data,[]);
            hold on;
            [r1,c1] = find(cap.rtBoundary);
            plot(c1,r1,'.y', 'MarkerSize', 0.25);
            % Show, XY view
            
            figure, imshow(im_xy,[])
            hold on;            
            plot(c,r,'.y', 'MarkerSize', 0.25);
        end
        % -------------------------------------------------------
        function iem = extractiem(obj, perimeter)
            
            % Find the last occurance of 1 in each row
            %[row, col] = find(perimeter);
            % with accumarray take the maximum column index for every row
            %iem = accumarray(row,col,[],@max);
            iem = zeros(obj.ImgHeight,1);
            for row=1:obj.ImgHeight
                oneLine = perimeter(row,:);
                loc = find(oneLine, 1, 'last');
                if loc
                    iem(row) = loc;
                end
            end
            
            % Since "perimeter" does not contain values in all image rows,
            % I will simply assign the border value of the closest non-zero
            % location
            idx = ~(iem==0); % => there will be a 0 where there is no value in the iem vector
            prev = 1;
            while prev<=numel(iem)
                next = find(idx(prev:end), 1, 'first') + prev - 1;
                if next
                    if prev<numel(iem)
                        for i=prev:next-1 % assign the lumen border value if the value of iem(next) is too far (I could do iem(prev:next-1) = iem(next);  but it would make it a constant
                            iem(i) = min(obj.lumenBorder(i), iem(next));
                        end
                        prev = next;
                    elseif prev==numel(iem)
                        iem(prev) = iem(prev-1);
                        break;
                    end
                    while idx(prev) && prev<numel(iem)% progress with the prev pointer up to the next 0
                        prev = prev + 1;
                    end
                else
                    break;
                end
            end
            
            
            dispResultFlag=0;
            if dispResultFlag
                % show the r-theta results
                imshow(log(double(obj.data+1)),[]);
                hold on;
                plot(iem, [1:numel(iem)], 'w-', 'MarkerSize', 0.25);
            end
        end
        % -----------------------------------------------------------------
        function filteredFrame = enhancedLee(obj, pbfn)
            
            %window size for (may want to add this as a GUI option           
            ws.w = 7;%window width (in r direction) for the sigma filter
            ws.h = 11;%window height for the sigma filter
            w_half.w=floor(ws.w/2);
            w_half.h=floor(ws.h/2);            
            
            fileInfo = imfinfo(pbfn);
            %disp(fileInfo(1));
            if numel(fileInfo) > 1 % check that I have more than a single image (I use single images when I analyze averaged stationary images)
                firstFrameNum = fileInfo(1).PageNumber(1) + 1;
                lastFrameNum = fileInfo(1).PageNumber(2);
                
                % MAKE SURE I AM NOT OUT OF BOUNDS WITHIN THE PULLBACK
                if (obj.index < lastFrameNum) && (obj.index > firstFrameNum)
                    Iprev = double(imread(pbfn, obj.index-1));
                    I_original = obj.data; % double(imread(pbfn, obj.index));
                    Ipost = double(imread(pbfn, obj.index+1));
                elseif obj.index == lastFrameNum
                    Iprev = double(imread(pbfn, obj.index-1));
                    I_original = obj.data; % double(imread(pbfn, obj.index));
                    Ipost = zeros(fileInfo(1).Height, fileInfo(1).Width);
                elseif obj.index == firstFrameNum
                    Iprev = zeros(fileInfo(1).Height, fileInfo(1).Width);
                    I_original = obj.data; % double(imread(pbfn, obj.index));
                    Ipost = double(imread(pbfn, obj.index+1));
                end
                % Concatenate rows from the images before and after the image of interest
                I = [Iprev(504-w_half.h+1:504,:); I_original; Ipost(1:w_half.h,:)];
            elseif numel(fileInfo) == 1
                I = double(imread(pbfn));
            end
             
            filteredFrame = enhancedLeeFilterC(double(I), ws.w);
            %filteredFrame = uint16(filteredFrame);
            %filteredFrame = obj.enhancedLeeFilterM(ws);  
            filteredFrame = filteredFrame(w_half.h+1:end-w_half.h, 1:end);
            
            dispFlag = false;
            if dispFlag
                figure, imshow(log(double(obj.data)+1),[]); title('Noisy image')
                figure, imshow(log(double(filteredFrame)+1),[]);title('Denoised image');
                %figure, imshow(log(double(filteredFrame2)+1),[]);title('MATLAB: Denoised image');
            end
            
        end
        % -----------------------------------------------------------------
        function filteredFrame = enhancedLeeFilterM(obj, ws)
            
            pbfn = obj.fn;
            fileInfo = imfinfo(pbfn);
            %disp(fileInfo(1));
            firstFrameNum = fileInfo(1).PageNumber(1) + 1;
            lastFrameNum = fileInfo(1).PageNumber(2);
            
            % MAKE SURE I AM NOT OUT OF BOUNDS WITHIN THE PULLBACK
            frameNum = obj.index;
            if (frameNum < lastFrameNum) && (frameNum > firstFrameNum)
                Iprev = double(imread(pbfn, frameNum-1));
                I_original = double(imread(pbfn, frameNum));
                Ipost = double(imread(pbfn, frameNum+1)); 
            elseif frameNum == lastFrameNum
                Iprev = double(imread(pbfn, frameNum-1));
                I_original = double(imread(pbfn, frameNum));
                Ipost = zeros(fileInfo(1).Width, fileInfo(1).Height);
            elseif frameNum == firstFrameNum
                Iprev = zeros(fileInfo(1).Height, fileInfo(1).Width);
                I_original = double(imread(pbfn, frameNum));
                Ipost = double(imread(pbfn, frameNum+1)); 
            end
            
            
            w_half.w=floor(ws.w/2);
            w_half.h=floor(ws.h/2);
            
            I = [Iprev(504-w_half.h+1:504,:); I_original; Ipost(1:w_half.h,:)];
            
            fr = zeros(size(I));
            L = 3; % Number of looks (input parameter)
           
            
            D = 1; % Damping factor (input parameter) - Specifies the damping factor to define the extent of smoothing (make 1.0 as default)
            Cu = 1/sqrt(L); % Noise variation coefficient
            Cmax = sqrt (1.0 + 2.0 / L); % Maximum noise variation coefficient)
            
            for iter=1:1  
                for row=1+w_half.h:(size(I,1)-w_half.h)
                    for col=1+w_half.w:(size(I,2)-w_half.w)
                        localWin = I(row-w_half.h:row+w_half.h,col-w_half.w:col+w_half.w);
                        Lm = mean(localWin(:));% Lm = Local mean of filter window
                        SD = std(double(localWin(:)));
                        Ci = SD/Lm; % Image variation coeficient      
                        P_c = I(row, col); % Center  pixel value of the window
%                         if ( (row==w_half_s+1) && (col==w_half_s+1) )
%                             fprintf('MATLAB: [Lm, SD, Ci, P_c] = [%f, %f, %f, %f]\n', Lm, SD, Ci, P_c);
%                             disp(localWin);
%                         end
                        K = exp(-D*(Ci-Cu)/(Cmax-Ci));
                        
                        if Ci <= Cu
                            fr(row, col) = Lm;
                        elseif (Ci > Cu) && (Ci < Cmax)
                            fr(row, col) = Lm*K + P_c*(1-K);
                        elseif (Ci >= Cmax)
                            fr(row, col) = P_c;
                        end
                     end
                end
            end
            filteredFrame = fr(w_half.h+1:end-w_half.h, 1:end);
           
            
%              figure;
%              subplot(1,2,1), imshow(log(double(I)+1),[]); title('Noisy image')
%              subplot(1,2,2), imshow(log(double(fr)+1),[]);title('Denoised image');
        end
        % -----------------------------------------------------------------
        function filteredFrame = enhancedLeeFilterM_original(obj, ws)
            
            
            Iprev = double(imread(obj.fn, obj.index-1));
            I = obj.data;
            Ipost = double(imread(obj.fn, obj.index+1));
            w_half_s = fix(ws/2);
            
            I = [Iprev(504-w_half_s:504,:);I;Ipost(1:w_half_s,:)];
            
            fr = zeros(size(I));
            L = 3; % Number of looks (input parameter)
           
            
            D = 1; % Damping factor (input parameter) - Specifies the damping factor to define the extent of smoothing (make 1.0 as default)
            Cu = 1/sqrt(L); % Noise variation coefficient
            Cmax = sqrt (1.0 + 2.0 / L); % Maximum noise variation coefficient)
            
            for iter=1:1  %preforming 3 iterations
                for row=1+w_half_s:(size(I,1)-w_half_s)
                    for col=1+w_half_s:(size(I,2)-w_half_s)
                        localWin = I(row-w_half_s:row+w_half_s,col-w_half_s:col+w_half_s);
                        Lm = mean(localWin(:));% Lm = Local mean of filter window
                        SD = std(double(localWin(:)));
                        Ci = SD/Lm; % Image variation coeficient      
                        P_c = I(row, col); % Center  pixel value of the window
                        
                        K = exp(-D*(Ci-Cu)/(Cmax-Ci));
                        
                        if Ci <= Cu
                            fr(row, col) = Lm;
                        elseif (Ci > Cu) && (Ci < Cmax)
                            fr(row, col) = Lm*K + P_c*(1-K);
                        elseif (Ci >= Cmax)
                            fr(row, col) = P_c;
                        end
                     end
                end
            end
            filteredFrame = fr;
            
%              figure;
%              subplot(1,2,1), imshow(log(double(I)+1),[]); title('Noisy image')
%              subplot(1,2,2), imshow(log(double(fr)+1),[]);title('Denoised image');
        end
        % --------------------------------------------------------------------
        % Correct for catheter for the whole image in one shot (as oppose to a single a-line)
        %
        % Inputs intensities=the a-line intensities vector to be corrected
        %           rVec = vector of the absolute indices of the a-line
        function Icorrected = correctForCatheter(obj, dispResultFlag)
            epsilon = 1e-8; % to avoid divide-by-zero            
           
            Zc = Globals.Zc;
            Zw = Globals.Zw;
            Z0 = Globals.Z0;
            Zr = Globals.Zr;
            
            T = @(r) 1./(epsilon + sqrt(((r-Z0)./Zr).^2+1)); % Confocal function (Van Soest paper)
            S = @(r) exp(-1.0*(((r-Zc)./Zw).^2));               % Roll-off function (Van Soest paper)
            
            % get a 504x968 matrix compsed of r repeated for each row
            rSingleLine = 1:size(obj.data, 2);
            r = rSingleLine' * diag(eye(size(obj.data,1)))'; 
            r = r' .* Globals.RTHETA_PIXEL_SIZE;  % Convert the r vector to micro-meters
            
            % Now, correct the intensities
            Icorrected = double(obj.data)./T(r);
            Icorrected = Icorrected./S(r);
            
            if dispResultFlag
                obj.displayImage(obj.data, 'rt', nan);
                title('Image before catheter correction');
                obj.displayImage(Icorrected, 'rt', nan);
                title('Image after catheter correction');            
            end
            
            obj.data = uint16(Icorrected);
             
        end
        % ----------------------------------------------------------------     
        % I added the return value of the frame just so I can use it to
        % save the lumen image for 3D vieweing.  I did not use the retun
        % value anywhere but in saveAsTiff() function and I did not check
        % if adding this value created a problem !!!
        function lumenBorderFrame = computeLumenBorder(obj, dispResultFlag)
            % Extract the lumen border in the frame and store it as a
                % 504x1 matrix (i.e index is theta number and the value is the r value
                catheterPos = 80;
                lumen_connectivity_rect = 25;
                bandwidth = 60;                    
                
                lumenBorderFrame = lumen_segmentation_DP(obj.data,catheterPos,lumen_connectivity_rect,bandwidth);                
                
                lb = [];
                %lb2=[];
                for row=1:size(lumenBorderFrame, 1)
                     rInd = find(lumenBorderFrame(row,:), 1, 'first');
                     %rInd2 = find(lumenBorderFrame2(row,:), 1, 'first');
                     lb(end+1) = rInd;
                     %lb2(end+1) = rInd2;
                end
                
                obj.lumenBorder = lb';
                obj.lumenBorder = obj.correctForGuidewire(obj.lumenBorder, 2);
                
%                 if ~isempty(obj.guidewire) % if the obj.guidewre is not empty in means that I wanted to remove the guidewire
%                     offset = 1;                    
%                     if obj.guidewire(2)-obj.guidewire(1)<(size(obj.data,1)/2)
%                         obj.lumenBorder(obj.guidewire(1):obj.guidewire(2)) = obj.lumenBorder(obj.guidewire(1)-offset);
%                     else
%                         obj.lumenBorder(1:obj.guidewire(1)+offset) = obj.lumenBorder(obj.guidewire(1)+offset);
%                         obj.lumenBorder(obj.guidewire(2)-offset:end) = obj.lumenBorder(obj.guidewire(2)-offset);
%                     end                    
%                 end
%                 
                if dispResultFlag
                    obj.displayImage(lumenBorderFrame, 'rt', nan);
                    % imtool(lumenBorderFrame); % just so I can see precise locations for analysis
                    title('Binary lumen border');
                    obj.displayImage(obj.data, 'rt', nan); hold on;                 
                    y_lumen = 1:numel(obj.lumenBorder); % Globals.NUM_LINES_PER_FRAME;
                    plot(obj.lumenBorder, y_lumen, 'r-', 'LineWidth', 4);
                    %plot(lb2, y_lumen, '.y', 'MarkerSize', 4);
                    title('Image on which lumen detection was implemented with lumen border overlay');
                    hold off;
                end
        end         
        % ----------------------------------------------------------------
        function correctedBorder = correctForGuidewire(obj, border, offset)
            
            correctedBorder = border;
            
            if ~isempty(obj.guidewire) % if the obj.guidewre is not empty in means that I wanted to remove the guidewire
              
                if obj.guidewire(2)-obj.guidewire(1)<(size(obj.data,1)/2)
                    correctedBorder(obj.guidewire(1)-offset:obj.guidewire(2)+offset) = border(obj.guidewire(1)-offset);
                else
                    correctedBorder(1:obj.guidewire(1)+offset) = border(obj.guidewire(1)+offset);
                    correctedBorder(obj.guidewire(2)-offset:end) = border(obj.guidewire(2)-offset);
                end
            end
            
            dispResultFlag = 0;
            if dispResultFlag
                obj.displayImage(obj.data, 'rt', 704); hold on;
                y_lumen = 1:numel(border);
                plot(border, y_lumen, 'r-', 'LineWidth', 4);
                plot(correctedBorder, y_lumen, 'y-', 'LineWidth', 4);
            end
            
        end
        % ----------------------------------------------------------------        
        function improvedComputeLumenBorder(obj, dispResultFlag)
            % Extract the lumen border in the frame and store it as a
                % 504x1 matrix (i.e index is theta number and the value is the r value
                catheterPos = 70;
                lumen_connectivity_rect = 30;
                bandwidth = 60;                    
                
                lumenBorderFrame1 = lumen_segmentation_DP(obj.data,catheterPos,lumen_connectivity_rect,bandwidth);
                
                 
                % take the first rIndex in each row (some rows have more
                % than one r index)
                lb = [];
                
                for row=1:size(lumenBorderFrame1, 1)
                     rInd = find(lumenBorderFrame1(row,:), 1, 'first');
                     %rInd2 = find(lumenBorderFrame2(row,:), 1, 'first');
                     lb(end+1) = rInd;
                     %lb2(end+1) = rInd2;
                end
                obj.lumenBorder = lb';
%                 lb_temp = bwboundaries(lumenBorderFrame);
%                 lb = lb_temp{1,1}(:,2);
%                 obj.lumenBorder = lb';
%                 y_lumen = lb_temp{1,1}(:,1);
                
                if dispResultFlag
                    obj.displayImage(lumenBorderFrame1, 'rt', 704);
                    % imtool(lumenBorderFrame); % just so I can see precise locations for analysis
                    title('Binary lumen border');
                    obj.displayImage(obj.data, 'rt', 704); hold on;                 
                    y_lumen = 1:numel(obj.lumenBorder); % Globals.NUM_LINES_PER_FRAME;
                    plot(obj.lumenBorder, y_lumen, 'r-', 'LineWidth', 3);
                    %plot(lb2, y_lumen, '.y', 'MarkerSize', 4);
                    title('Image on which lumen detection was implemented with lumen border overlay');
                    hold off;
                end
        end         
        % -----------------------------------------------------------------
        function removeBaseline(obj, dispResultFlag)
            
            baseline = 6.7; % (true value we found was 6.7, need to put this in the Globals
            obj.data = obj.data - baseline;
            obj.data(obj.data<0.0) = 0;     %THIS CAUSES THE IMAGE TO DARKEN.  SOLVE IT !!      
            
            if dispResultFlag
                obj.displayImage(obj.data, 'rt', 512); hold on;                           
            end
        end
        % ----------------------------------------------------------------        
        function backBorderFrame = computeBackBorder(obj, dispResultFlag)
            % Extract the back border in the frame and store it as a
            % 504x1 matrix (i.e each row is the theta index and its value is the r index)
            
            %method = 'simple addition';                  
            %method = 'image processing';
            method = 'dynamic programming';
            
            switch method
                case 'simple addition'                    
                    % Just add 400 pixels (~2mm=2000 microns)
                    obj.backBorder = obj.lumenBorder + 400;
                case 'dynamic programming'          
                    backBorderFrame = obj.dpa(); % Dynamic Programming Algorithm
                    % Now, translate the binary frame to a vector
                    border = [];                    
                    for row=1:size(backBorderFrame, 1)
                        rInd = find(backBorderFrame(row,:), 1, 'first');
                        border(end+1) = rInd;
                    end
                    %figure, imshow(backBorderFrame, []);
                    obj.backBorder = border';
                case 'image processing'
                    % Here I use marker-controlled watershed segmentation: see EBME book page 593 (watershed example)
                    % I consider as "internal" the area to the right of the lumen (i.e. what's
                    % behind the media layer)                    
                    h = fspecial('sobel');
                    g = sqrt(imfilter(double(obj.data), h, 'replicate') .^ 2 + imfilter(double(obj.data), h', 'replicate') .^ 2);
                    %rm = imregionalmin(g);
                    %rm = imregionalmin(obj.data); % compute the location of all regional minima in an image
                    heightThresh = 10;  % set a threshold
                    im = imextendedmin(obj.data, heightThresh); % im is an image whose foreground pixels mark the locations of the deep regional minima
                    % superimpose the data into the original image
                    fim = obj.data;
                    fim(im) = 175;
                    Lim = watershed(bwdist(im));
                    em = (Lim == 0); % em is a logical matrix with zeros and 1's
                    %             g2 = imimposemin(g, im | em);  % modify the gradient image using a "minima imposition" (see Soille, 2003)
                    %             %g2(g2==-inf) = 0;
                    %             L2 = watershed(g2);
                    %             f2 = obj.data;
                    %             f2(L2 == 0) = 255;
                    
                    % Now, to extract the back border I scan the rows from the
                    % right and mark as border the first slope greater than 1
                    lineCounter = 1;
                    obj.backBorder = zeros(size(em,1),1);
                    for row=1:size(em,1)
                        slopes = diff(em(row, :));
                        largestSlopeLocation = (slopes == 1);
                        ind = find(largestSlopeLocation, 1, 'last');
                        if isempty(ind)
                            lineCounter = lineCounter + 1;
                        else
                            obj.backBorder(lineCounter) = ind;
                            lineCounter = lineCounter + 1;
                        end
                    end
                    
                    obj.backBorder = obj.backBorder + 50; % arbitrarily shift everything to the right by 50 pixels
            end
           
            obj.backBorder(isnan(obj.backBorder)) = obj.ImgHeight; % obj.backBorder is initialized to nan.  Need to make sure non is left
             
            if dispResultFlag
%                 figure, imshow(im, []); hold on; title('Extended Minima Transform - set of " low spots" in the image that are deeper (by a certain height threshold) than their immediate surroundings');
%                 figure, imshow(fim, []); hold on; title('Extended minima locations superimposed on the original image');
%                 figure, imshow(Lim, []); hold on; title('Watershed transform of the distance transform of the Extended minima image');
%                 figure, imshow(em, []); hold on; title('Truncated-to-zero all watershed areas to mark the foreground');
                
                
                obj.displayImage(obj.originalData, 'rt', 704); hold on;
                %%%%%%%%%%%%obj.displayImage(obj.data, 'rt', 704); hold on;
                y_lumen = 1:size(obj.data,1);
                
                %plot(obj.backBorder, y_lumen, 'y-', 'LineWidth', 2);
                plot(obj.backBorder, y_lumen, 'y-', 'LineWidth', 2);
                plot(obj.lumenBorder, y_lumen, 'r-', 'LineWidth', 2);
                %title('(r,\theta) Image with back border overlay');                
                
               obj.displayImage(obj.data, 'xy', 704); hold on;
               % gXY = imread('C:\Users\rys.ADS\Documents\Projects\TissueCharacterization\OCT_Data\ValidationData\Prabu\Vessel 64\OCT Data\prox LAD_XY_704x704.tif', 267);
                %gXY = imrotate(gXY,-90);
                %gXY = imresize(gXY, 0.98);
                %figure, imshow(gXY,[]); hold on;
                % Use Sepia colormap (I use a fcn from the file exchange
%      AdvancedColormap('bone'); 
%     AdvancedColormap('invert'); 
%     AdvancedColormap('reverse');
%     AdvancedColormap('kryw');
                
                
                %[x,y,cp] = obj.convertRTpointToCartesian(obj.lumenBorder, 1, 704)                   
                
                x=zeros(size(obj.lumenBorder,1), 1);
                y=zeros(size(obj.lumenBorder,1), 1);
                for i=1:size(obj.lumenBorder,1)
                    xy_pt = obj.convertRTpointToCartesian([obj.lumenBorder(i); i], 1, 704) ;        
                    x(i) = xy_pt(1);
                    y(i) = xy_pt(2);
                end
                centroid.x = sum(x(:))/size(x,1);
                centroid.y = sum(y(:))/size(y,1);               
               
                plot(x, y, 'r-', 'LineWidth', 4); hold on;
                %plot(centroid.x, centroid.y, '+r', 'MarkerSize', 5);
                
                x=zeros(size(obj.backBorder,1), 1);
                y=zeros(size(obj.backBorder,1), 1);
                for i=1:size(obj.backBorder,1)
                    xy_pt = obj.convertRTpointToCartesian([obj.backBorder(i); i], 1, 704) ;        
                    x(i) = xy_pt(1);
                    y(i) = xy_pt(2);
                end
                                
                plot(x, y, 'y-', 'MarkerSize', 4);
                title('(x,y) Image with back border overlay');
            end
            
        end
        % -----------------------------------------------------------------
        function verifyCoords(obj, coords)
            
            mode = 'rt'; % 'xy'
            %% display the original image
            obj.displayImage(obj.data, mode, 704); hold on;   
            %% Display the lumen border
            x=zeros(size(obj.lumenBorder,1), 1);
            y=zeros(size(obj.lumenBorder,1), 1);
            if isempty(obj.lumenBorder)
                obj.computeLumenBorder(false);                
            end
            for i=1:size(obj.lumenBorder,1)
                if strcmp(mode, 'xy')
                    xy_pt = obj.convertRTpointToCartesian([obj.lumenBorder(i); i], 1, 704) ;
                    x(i) = xy_pt(1);
                    y(i) = xy_pt(2);
                else
                    x(i) = obj.lumenBorder(i);
                    y(i) = i;
                end
            end
            plot(x, y, '.r', 'MarkerSize', 2); hold on;
            %% Display the centroid location
            if strcmp(mode, 'xy')
                centroid.x = sum(x(:))/size(x,1);
                centroid.y = sum(y(:))/size(y,1);
                plot(centroid.x, centroid.y, '+r', 'MarkerSize', 5);
            end
            %% Display the back-border
             if isnan(obj.backBorder)
                obj.computeBackBorder(false);                
            end
            x=zeros(size(obj.backBorder,1), 1);
            y=zeros(size(obj.backBorder,1), 1);
            for i=1:size(obj.backBorder,1)
                if strcmp(mode, 'xy')
                    xy_pt = obj.convertRTpointToCartesian([obj.backBorder(i); i], 1, 704) ;
                    x(i) = xy_pt(1);
                    y(i) = xy_pt(2);
                else
                    x(i) = obj.backBorder(i);
                    y(i) = i;                    
                end
            end            
            plot(x, y, '.y', 'MarkerSize', 2);
            
            %% Now, display the coordinates of the classified image
            x=zeros(size(coords,1), 1);
            y=zeros(size(coords,1), 1);
            for i=1:size(coords,1)
                rt_x = coords(i,1);
                rt_y = coords(i,2);
                if strcmp(mode, 'xy')
                    xy_pt = obj.convertRTpointToCartesian([rt_x; rt_y] , 1, 704) ;
                    x(i) = xy_pt(1);
                    y(i) = xy_pt(2);
                else
                     x(i) = rt_x;
                     y(i) = rt_y;
                 end
                plot(x(i), y(i), '.g', 'MarkerSize', 1);
            end
            
            plot(x, y, '.g', 'MarkerSize', 1);
        end
        % -----------------------------------------------------------------
        function I = pixShift(obj, dispResultFlag)
            
            if isempty(obj.lumenBorder)
                obj.computeLumenBorder(false);                
            end
            
            leftAnchor = min(obj.lumenBorder); % THAT WILL BE THE COLUN TO WHICH ALL a-LINES WILL BE SHIFTED TO
            
            if isnan(obj.backBorder)
                obj.computeBackBorder(false);                
            end
            
            % shift to the left and pad right side with zeros
            I = double(zeros(numel(obj.lumenBorder), size(obj.data,2)-leftAnchor));            
            for i=1:numel(obj.lumenBorder)
                aLineI = obj.data(i, obj.lumenBorder(i):obj.backBorder(i));
                I(i, 1:numel(aLineI)) = aLineI;                
            end            
            
            % show only the part between the borders
            I2 = double(zeros(size(obj.data)));
            %I2 = obj.data(:, obj.lumenBorder:obj.backBorder);
            for i=1:numel(obj.lumenBorder)
                low = obj.lumenBorder(i);
                hi = obj.backBorder(i);
                aLineI = obj.data(i,low:hi);
                I2(i, low:hi) = aLineI;                
            end 
            
            %imwrite(I, 'temp.tif');
            
            if dispResultFlag
                fig;
                imshow(log(double(obj.data)+1),[]); hold on;title('Original');
                y_lumen = 1:numel(obj.lumenBorder); % Globals.NUM_LINES_PER_FRAME;
                plot(obj.lumenBorder, y_lumen, '.r', 'MarkerSize', 4);
                plot(obj.backBorder, y_lumen, '.y', 'MarkerSize', 4);
                fig;
                imshow(log(double(I)+1), []); hold on;title('Shifted');
                fig;
                imshow(log(double(I2)+1), []); hold on;title('only ROI');
            end            
        end
        % -----------------------------------------------------------------
        % The Sonka DP implementation is based on the Sonka book, 3rd edition (p. 210)
        % it takes into account only the 8-neighborhood..  Zaho's DP
        % implementation is more flexible and lets me select neighber=hood
        % size.  Therefore, the entgropy partial cost works better with
        % Zaho's DP implementation of DP.  But, Zaho's cost function is
        % also better and more stable than  my entropy partial cost
        % approach.
        function border = borderDetectDP(obj, dispResultFlag)
            %     lumenBorderFrame = lumen_segmentation_DP(obj.data,catheterPos,lumen_connectivity_rect,bandwidth);
            
            % First, assign the cost (entropy of local neighborhood defined
            % by nhood) for each pixel
            
            catheterPos = 70;
            connectivity = 11;            
            add_lines = 20;
            
            method = 'entropy';
            %method = 'mean'; % Zaho's method
            %method = 'gradient';
            %method = 'invertedPixelGrayValues';
            
            im = obj.data;
            im = [im(end-add_lines+1:end,:);im;im(1:add_lines,:)];
            [M, N] = size(im);
            costImg = zeros(M,N);
            
            switch method
                case 'entropy' % Note that Sonka's DP minimizes cost, yet Zaho maximizes cost
                    ww = 11;  % window width (r direction)
                    wh = 5; % window height
                    %nhood = true(ww, wh);
                    se = strel('rectangle', [wh, ww]);
                    se = translate(se,[0 -floor(ww/2)]); % shift it to the left so that I use only the entropy to the left of the point
                    nhood = getnhood(se);
                    costImg = entropyfilt(im, nhood);
                    %costImg_complement = imcomplement(costImg);
                    costImg(1:obj.ImgHeight, 1:catheterPos) = -max(max(costImg)); % to make sure that the catheter and guide wire are out of consideration
                    %border = sonka.dpboundary( -double(costImg) ); % Sonka implementatiopn
                    [borderBinaryImg, index] = DP_fast(costImg',connectivity); % Notice that DP_fast requires a rotated image
                    borderBinaryImg = borderBinaryImg';
                    borderBinaryImg = connect_border(borderBinaryImg);
                    border = index;
                case 'gradient'                    
                    [Gmag,Gdir] = imgradient(im); % get the gradient and gradient direction images
                    costImg = Gmag; 
                    %costImg = Gdir;
                    costImg(1:obj.ImgHeight, 1:catheterPos) = 0; % to make sure that the catheter and guide wire are out of consideration
                    
                    %border = sonka.dpboundary( -double(costImg) );
                    [borderBinaryImg, index] = DP_fast(costImg',connectivity); % Notice that DP_fast requires a rotated image
                    borderBinaryImg = borderBinaryImg';
                    borderBinaryImg = connect_border(borderBinaryImg);
                    border = index;
                    
                case 'mean' %Zaho's method
                    band_width = min(40,catheterPos-1);                    
                    for j=catheterPos:N-band_width
                        costImg(:,j)=mean(im(:,j:min(j+band_width,N)),2)-mean(im(:,j-band_width:j),2);
                    end      
                    costImg(1:obj.ImgHeight, 1:catheterPos) = 0; % to make sure that the catheter and guide wire are out of consideration
                    
                    %% Now, implement the dynamic programming
                    [borderBinaryImg, index] = DP_fast(costImg',connectivity); % Notice that DP_fast requires a rotated image
                    borderBinaryImg = borderBinaryImg';
                    
                    %            borderBinaryImg = borderBinaryImg(add_lines+1:end-add_lines,:,:);
                    %            index = index(add_lines+1:end-add_lines);
                    borderBinaryImg = connect_border(borderBinaryImg);
                    border = index;
                case 'invertedPixelGrayValues'
                    costImg = im;
                    costImg(1:obj.ImgHeight, 1:catheterPos) = -65536;
                    %border = sonka.dpboundary( -double(costImg) );
                    border = sonka.dpboundary( -double(costImg) );
                    
            end
            
            
            
            %% Plot stuff
            if dispResultFlag
                figure, imshow(log(double(im)),[]); %original Image
                hold on, plot(border, [1:numel(border)], 'r-');
                %figure, imshow(log(double(costImg)),[]); % partial cost image
                %figure, imshow(log(double(borderBinaryImg)),[]); % binary image of the border
            end
        end
         % -----------------------------------------------------------------
        % Computes the border 
        function border = dpa(obj)
            
            neighberhood = 9; 
            add_lines = 20;% To guarantee the contour is closed, there may be better ways to do it ( I borrowed Zaho's method for now, but need to improve it!!)
            
            f = obj.computePartialCostMatrix(add_lines);    % Compute the individual pixel cost values matrix            
            [border,idx]=dpaC(f', neighberhood);                % Compute the Cost function using dynamic programming
            border=border';
            
            border = border(add_lines+1:end-add_lines,:,:);
            idx = idx(add_lines+1:end-add_lines);
            
            border = connect_border(border);                    % Connect the border points (i.e. back track the lowest cost path)
            
            %figure, obj.displayImage(border, 'rt', 704);            
        end
        % -----------------------------------------------------------------
%         function backBorder = dpa_mine(obj)
%             
%             %f = @(theta,r,ws) mean(obj.data(theta, r-floor(ws/2):r+floor(ws/2)));
%             wh = 10; % window height (theta)
%             ww = 20; % window width (r direction)
%             
%             if isempty(obj.lumenBorder)
%                 obj.computeLumenBorder(false);                
%             end
%             
%             
%            
%             % Pre allocated 2D matrix of classes (remember: the constractor
%             % needs to have the option of zero nargin)
%             grid = DPCell.empty(0, size(obj.data, 2));
%             
% %             h = waitbar(0,'Loading pullback...');
% %             %creates  GRID of cells with 'DPCell' class object
% %             for i=1:size(obj.data,1)
% %                 for j=1:size(obj.data,2)
% %                     gridCell = DPCell(i, j);
% %                     grid(i,j) = gridCell;
% %                     clear gridCell;
% %                 end
% %                 waitbar(i/size(obj.data,1),h,'Creating empty grid...');
% %             end
% %             close(h);
%             
%             h = waitbar(0,'Creating Cost Grid...');
%             methods = {'mean', 'absolute mean', 'second diff', 'entropy'};
%             costMap = nan(size(obj.data));
%             firstRow = ceil(wh/2);
%             firstColumn = min(obj.lumenBorder); % ceil(ww/2); %min(obj.lumenBorder);
%             lastRow = obj.ImgHeight-ceil(wh/2); %  for debugging just scan 4th height : floor(obj.ImgHeight/4); %
%             
%            
%             for theta=firstRow:lastRow;
%                 costHist = [];
%                 for r=firstColumn:obj.ImgWidth-ceil(ww/2) % for r=obj.lumenBorder(theta):obj.ImgWidth-ceil(ww/2)
%                     newCell = DPCell(theta, r);
%                     grid(theta,r) = newCell;
%                     cost = obj.dpCost(theta,r,wh,ww, methods{2});
%                     %cost = cost  * obj.dpCost(theta,r,wh,ww, methods{4});
%                     if theta==firstRow          
%                         grid(theta,r).set(0, [0,0]);                        
%                         costMap(theta, r) = 0;
% %                         grid(theta,r).set(cost, [0,0]);                        
% %                         costMap(theta, r) = cost;
%                     elseif r==firstColumn                      
%                         c1 = costMap(theta-1, r);
%                         c2 = costMap(theta-1, r+1);                       
%                         
%                         [minValue, idx] = min([c1, c2]);                        
%                         if idx==1
%                             predecessor = [theta-1, r];
%                         elseif idx==2
%                             predecessor = [theta-1, r+1];
%                         end
%                         grid(theta,r).set(minValue + cost, predecessor);
%                         costMap(theta, r) = minValue + cost;
%                     else
%                         c1 = costMap(theta-1, r-1);
%                         c2 = costMap(theta-1, r);
%                         c3 = costMap(theta-1, r+1);                       
%                         
%                         [minValue, idx] = min([c1, c2, c3]);
%                         
%                         if idx==1
%                             predecessor = [theta-1, r-1];
%                         elseif idx==2
%                             predecessor = [theta-1, r];
%                         elseif idx==3
%                             predecessor = [theta-1, r+1];
%                         end
%                         grid(theta,r).set(minValue + cost, predecessor);                        
%                         costMap(theta, r) = minValue + cost;
%                     end                    
%                
%                     costHist(end+1) = cost;
%                 end
%                 waitbar(theta/lastRow,h,'Creating Cost grid...');
%                 %figure,  plot(costHist, '-r'); title('Cost Function - Neighberhood Entropy');
%               
%                 %costMap(theta,:) = costHist;
%             end
%             close(h);
%             
%             % Now, back track and get the best path
%             lowestLoc = [lastRow, firstColumn];
%             lowestCost = grid(lowestLoc(1), lowestLoc(2)).optimalCost;
%             costRow = [];
%            for i=1:size(grid,2)
%                costRow = [costRow grid(lastRow, i).optimalCost];
%            end
%            lowestY = lastRow;
%            [lowestValue, lowestX] = min(costRow);
% 
%             
%             backBorder = zeros(lastRow,1);
%             prevX = lowestX;
%             for theta=lowestY:-1:firstRow;
%                try
%                 backBorder(theta) = grid(theta, prevX).predecessor(2);             
%                catch
%                    disp('');
%                end
%                 prevX = grid(theta, prevX).predecessor(2);
%             end
%             %obj.displayImage(costMap, 'rt');
%             figure, obj.displayImage(obj.data, 'rt'); hold on; 
%             plot (backBorder, [1:size(backBorder,1)], 'r-');
%             title('Backborder');            
%         end
        % -----------------------------------------------------------------
        % This function fills up a matrix whose enetries are the cost
        % values for each pixel.  IT'S THE ORIGINAL ONE WHICH WOKS WELL.  
        % THE NEXT ONE (FOLLOWING FCN) IS ONE i USE TO PLAY WITH VARIOUS
        % COST FUNCTIONS
%         function f = computePartialCostMatrixOld(obj, add_lines)            
%             
%             methods = {'mean', 'absolute mean', 'second diff', 'entropy'};
%             wh = 10; % window height (theta direction)
%             ww = 20; % window width (r direction)
%             
%             if isempty(obj.lumenBorder)
%                 obj.computeLumenBorder(false);
%             end
%             firstColumn = min(obj.lumenBorder) + 40; % ceil(ww/2); %min(obj.lumenBorder);
%             firstRow = ceil(wh/2);            
%             lastRow = obj.ImgHeight-ceil(wh/2);   
%            
%             % pad the image with "add_lines" on top and bottom            
%             im = [obj.data(end-add_lines+1:end,:); obj.data; obj.data(1:add_lines,:)];
%             
%             M = size(im,1);  
%             N = size(im,2);
%             
%             f = zeros(M,N);
%             
%             % This is Zaho's cost function for lumen detection
% %             bandwidth = 40;
% %             for j=firstColumn:N-bandwidth
% %                 f(:,j)=mean(im(:,j:min(j+bandwidth,N)),2)-mean(im(:,j-bandwidth:j),2);
% %             end
%             
%             bandwidth = 120;
%             %for theta=firstRow:lastRow;
%                 for r=firstColumn:obj.ImgWidth-ceil(ww/2)                  
%                     %cost = obj.dpCost(theta,r,wh,ww, methods{2});
%                     %f(theta,r) = cost * obj.dpCost(theta,r,wh,ww, methods{4});
%                     %f(:,r) = abs(mean(im(:,r:min(r+bandwidth,N)),2) - mean(im(:,r-bandwidth:r),2));
%                     %f(:,r) = abs(mean(im(:,min(r+bandwidth, firstColumn):r),2) - mean(im(:,r:max(r-bandwidth,N)),2));
%                     
%                    
%                     leftSide = im(:,max(r-bandwidth, firstColumn):r);
%                     rightSide = im(:, r:min(r+bandwidth,N));
%                     f(:,r) = abs(mean(leftSide, 2) - mean(rightSide, 2));
%                     
% %                     leftSide = im(:,max(r-bandwidth, firstColumn):r);
% %                     rightSide = im(:, r:min(r+bandwidth,N));
% %                     pLeft = imhist(leftSide);
% %                     pRight = imhist(rightSide);
% %                     
% %                     eLeft = sum(pLeft .* log2(pLeft+1));
% %                     el = eLeft * ones(size(f,1),1);
% %                     entropyColRight = sum(pRight .* log2(pRight+1));
% %                     entropyColLeft = eRight * ones(size(f,1),1);                 
%                     
%                     if isa(f, 'uint16') || isa(f, 'double')
%                         f(:,r) = (2^16-1)*ones(size(f,1),1) - f(:,r);
%                     elseif isa(obj.data, 'int8')
%                        f(:,r) = (2^8-1)*ones(size(f,1),1) - f(:,r);
%                     else
%                         disp(['image data type is' class(obj.data)]);
%                     end                    
%                 end
%             %end
%         end
        % -----------------------------------------------------------------
        % This function fills up a matrix whose enetries are the cost
        % values for each pixel.  I use this function TO PLAY WITH VARIOUS
        % COST FUNCTIONS
        function f = computePartialCostMatrix(obj, add_lines)            
            
            methods = {'mean', 'absolute mean', 'second diff', 'entropy'};
            wh = 5; % window height (theta direction)
            ww = 50; % window width (r direction)
            
            if isempty(obj.lumenBorder)
                obj.computeLumenBorder(false);
            end
            firstColumn = min(obj.lumenBorder) + 40; % ceil(ww/2); %min(obj.lumenBorder);
            firstRow = ceil(wh/2);            
            lastRow = obj.ImgHeight-ceil(wh/2);   
           
            % pad the image with "add_lines" on top and bottom (BETTER TO PADD WITH PREV AND POST IMAGE DATA)        
            im = [obj.data(end-add_lines+1:end,:); obj.data; obj.data(1:add_lines,:)];
            
            M = size(im,1);  
            N = size(im,2);
            
            f = zeros(M,N);           
           
            bandwidth = 80;
            lineEntropyVec = [];
            %fig;
            %for theta=firstRow:lastRow;
                for r=firstColumn:obj.ImgWidth-ceil(ww/2)                  
                    %cost = obj.dpCost(theta,r,wh,ww, methods{2});
                    %f(theta,r) = cost * obj.dpCost(theta,r,wh,ww, methods{4});
                    %f(:,r) = abs(mean(im(:,r:min(r+bandwidth,N)),2) - mean(im(:,r-bandwidth:r),2));
                    %f(:,r) = abs(mean(im(:,min(r+bandwidth, firstColumn):r),2) - mean(im(:,r:max(r-bandwidth,N)),2));
                    
                   % method 1: mean of window to left - mean window to
                   % right (window width of 20 works nice)
                    leftSide = im(:,max(r-bandwidth, firstColumn):r);
                    rightSide = im(:, r:min(r+bandwidth,N));
                    f(:,r) = abs(mean(leftSide, 2) - mean(rightSide, 2)); % computed for a complete column
                     
                    % method 2: compute entropy to the right and use it as
                    % the partial cost (need to uncomment the theta loop
                    % and comment the below since it is implemented in the dpCost fcn)
                    %[f(theta, r), originalPartialCost] = obj.dpCost(theta, r, wh, ww, N, methods{4});
                    %lineEntropyVec(end+1) = f(theta,r); % originalPartialCost;
                    if isa(f, 'uint16') || isa(f, 'double')
                        f(:,r) = (2^16-1)*ones(size(f,1),1) - f(:,r);
                    elseif isa(obj.data, 'int8')
                       f(:,r) = (2^8-1)*ones(size(f,1),1) - f(:,r);
                    else
                        disp(['image data type is' class(obj.data)]);
                    end                    
                end % of r
                %plot(lineEntropyVec);
            %end % of theta loop
            %lineEntropyVec = [];
        end

        % ------------------------------------------------------------------
        function [partialCost, originalPartialCost] = dpCost(obj, theta, r, wh, ww, N, method)
            
            switch method
                case 'mean'
                    thetaRange = [theta-floor(wh/2)+1:theta+floor(wh/2)];
                    rRange = [r-floor(ww/2)+1:r+floor(ww/2)];
                    
                    window = obj.data(thetaRange, rRange);
                    partialCost = mean(window(:));
                case 'absolute mean'
                    thetaRange = [theta-floor(wh/2)+1:theta+floor(wh/2)];
                    rRange = [r-floor(ww/2)+1:r+floor(ww/2)];
                    
                    window = obj.data(thetaRange, rRange);
                    partialCost = abs(mean(window(:)));                    
                case 'second diff'
                    h = 0.001;
                    slopes = diff(obj.data(theta, :))/h; %first derivative                       
                    partialCost = min(diff(slopes)/h); % second derivavtive
                case 'entropy'
                    % if I want the currect location to be at the center of
                    % the ROI                    
                    %thetaRange = [theta-floor(wh/2)+1:theta+floor(wh/2)];
                    %rRange = [r-floor(ww/2)+1:r+floor(ww/2)];
                    
                    % if I want the ROI to be to the right of my current
                    % location
                   thetaRange = [theta-floor(wh/2)+1:theta+floor(wh/2)];
                   % thetaRange = theta;
                    rRange = [r:min(r+ww, N)];
                    
                    % definition of entropy: -sum [over all possible values
                    % of x] of (p(x==i) .* log2(p(x==i)) where p is the
                    % probability that x equals i (i.e. 0<p<1)
                   window = obj.data(thetaRange, rRange);                   
                   [n, ~] = hist(window);
                   p = n/sum(n); % scale the distributions properly to transform into probabilites.
                   e = -sum(p .* log2(p+0.001));                    
                   %partialCost = entropy(window);
                   partialCost = e;
                   originalPartialCost = e;
            end
            
            % I subtract the value from 255 because my dynamic programming maximizes (not minimizes)
            if isa(obj.data, 'uint16') || isa(obj.data, 'double')
                partialCost = (2^16-1)-abs(mean(window(:)));
            elseif isa(obj.data, 'int8')
                partialCost = (2^8-1)-abs(mean(window(:)));
            else
                disp(['image data type is' class(obj.data)]);
            end
                    
        end
        %------------------------------------------------------------------
        function [muts, Inots, mse] = preComputeAttributes(obj, showResultsFlag)
            
            slidingWinSize = 7;
            for theta=1:obj.ImgHeight
                [muts, Inots, mse] = obj.leastSquares(theta, slidingWinSize, showResultsFlag);
            end
            
        end
        %-----------------------------------------------------------------
        % return values:
        %        gof : goodness of fit measure (for now I use mse)
        % 
        function [mutVector, InotVector, mseVector] = leastSquares(obj, thetaIdx, ws, plotResultsFlag)

            % Init output values
            mutVector = nan(1, obj.backBorder(thetaIdx)-obj.lumenBorder(thetaIdx));
            InotVector = nan(size(mutVector));
            mseVector = nan(size(mutVector));
            
            %method = 'Least Squares';
            method = 'Linear Regression';
            switch method
                case 'Linear Regression' % this is the least square method Wilson sent me a short document about (not sure why the name)
               %tic;
                     for r=obj.lumenBorder(thetaIdx):obj.backBorder(thetaIdx)
                        r_idx = r:r+ws;                        
                        r_mm = r_idx * Globals.RTHETA_PIXEL_SIZE;
                        windowIntensities = log(double(obj.data(thetaIdx, r_idx)) + 1);
                        
                        n = size(r_mm,2);
                        r_mm_mean = mean(r_mm);
                        windowIntensities_mean = mean(windowIntensities);
                        SS_xx = (r_mm.^2) * ones(n, 1) - n*r_mm_mean^2;
                        SS_xy =  (r_mm - r_mm_mean) * (windowIntensities-windowIntensities_mean)';
                        
                        mutVector(r) = (-1.0) * (SS_xy/SS_xx); % multiply by -1 since the algo computes the slope and mu_t = -slope
                        InotVector(r) = exp(windowIntensities_mean - mutVector(r)*r_mm_mean); % check the version in Aline.m, I think I need to change sign
                        
                     end
                     %fprintf('LSQ: took %3.3f seconds\n', toc);
                case 'Least Squares'
                %tic;
                    for r=obj.lumenBorder(thetaIdx):obj.backBorder(thetaIdx)
                        r_idx = r:r+ws;
                        r_mm = r_idx * Globals.RTHETA_PIXEL_SIZE;
                        windowIntensities = log(double(obj.data(thetaIdx, r_idx)) + 1);
                        
                        [p, s] = polyfit(r_mm, windowIntensities, 1);
                        mutVector(r) = (-1.0) * p(1); % multiply by -1 since the fcn computes the slope and mu_t = -slope
                        logInot = p(2);
                        InotVector(r) = exp(logInot);
                        
                        % Compute the mse
                        estimated = polyval(p, r_mm);
                        observed = windowIntensities;
                        mseVector(r) = mean(((estimated - observed).^2)); % remember: rms computes the square of the differences from the mean, yet residual is the difference from the fit
                        if plotResultsFlag
                            plot(r_mm, observed, 'o', r_mm, estimated,'-');
                        end
                    end
                    %fprintf('LSQ: took %3.3f seconds\n', toc);
            end
        end
        % -----------------------------------------------------------------
        function rv = getSingleLineIntensities(obj, theta, rStart, depth, imgWidth)
            rv = obj.data(theta, [rStart: min((rStart + depth-1), imgWidth)]); 
        end
        % -----------------------------------------------------------------
        function rv = getLineIntensities(obj, varargin)
            
            if nargin==3 %if I call it with a rectangular box shape
                thetaRange = varargin{1};
                rRange = varargin{2};              
                rv = obj.data(thetaRange, rRange);               
            elseif nargin==4 %if I call it with different aline lengths I return a matrix where I pad the missing values with zeros
                thetaRange = varargin{1};
                actualBorder = varargin{2};
                box =  varargin{3};
                
                rightSide = max( max(actualBorder) + box.width, box.origin(2) + box.width); % NOT GOOD, THIS WILL MODIFY THE BOX SIZE
                leftSide = min( min(actualBorder), box.origin(2));
                maxRangeVec = [leftSide:rightSide];
                rv = ones(size(thetaRange,2), size(maxRangeVec, 2)) * (-1); % MAKE IT ZEROS NOT -1
                for i=1:size(thetaRange,2)
                    v = obj.data(thetaRange(i), [actualBorder(i):max(maxRangeVec)]);
                    rv(i, size(maxRangeVec,2)-length(v)+1:end) = v;
                end
            else
                fprintf('function is called with wrong number of parameters');
                disp('');
            end
        end        
        % ************* Accessors *********************************
        function [labels, probabilities] = getAlinePredictions(obj)
            labels = obj.aLineLabel;
            probabilities = obj.aLineLableProb;
        end
        % ----------------------------------------------------------------
        function rv = getData(obj)
            rv = obj.data;
        end
        % -----------------------------------------------------------------
        function [h, w] = getImageDim(obj)
            
            [h, w] = size(obj.data);
        end
        % -----------------------------------------------------------------        
        % Output: A vector containing the column number (in pixels)
        % of the first nonzero border pixel of the lumen 
        % in each row.  
        function column = getLumenBorder(obj)
            if isempty(obj.lumenBorder)
                obj.computeLumenBorder(false);                
            end
            column = obj.lumenBorder;
        end
        % -----------------------------------------------------------------        
        function column = getBackBorder(obj)            
            if isnan(obj.backBorder)
                obj.computeBackBorder(false);
            end
            column = obj.backBorder;
        end
        % -----------------------------------------------------------------       
        function writeFrameToDisk(obj, fn)  
            imwrite(obj.data, fn);
            fprintf('Wrote %s to disk\n', fn);
        end
        % -----------------------------------------------------------------       
        % Note: the colors displayed by the St-Jude OCT machine is called
        % Sepia.  I read somewhere that sepia tone colorization of grayscale
        % photographs can be achived using Matlab's built-in pink, i.e.  colormap(pink)
        function [hdl_rt, hdl_xy] = displayImage(obj, I, what, s) 
            figure;
            switch what
                case 'rt'                    
                    hdl_rt = imshow(log(double(I)+1),[]);
                    hdl_xy = nan;
                case 'xy'
                    I_xy = rectopolar_fast(double(I)',s);                    
                    hdl_xy = imshow(log(I_xy+1),[]);
                    hdl_rt = nan;
                case 'both'                    
                    I_xy = rectopolar_fast(im2double(I)',s);
                    subplot(121), hdl_rt = imshow(log(double(I)+1),[]);
                    subplot(122), hdl_xy = imshow(log(double(I_xy)+1),[]);
            end
            datacursormode on;
        end                
        % ----------------------------------------------------------------
        function plotWindow(obj, origin)
            hold on;
            w = obj.movingWindowOptions.width;
            h = w * obj.movingWindowOptions.aspectRatio;
            
            x = [origin(1) origin(1)+w origin(1)+w origin(1) origin(1)];
            y = [origin(2) origin(2) origin(2)+h origin(2)+h origin(2)];
            obj.recHdl = plot(x,y,'-','Color','r','LineWidth',1);
         end   
        
        % ----------------------------------------------------------------
        function doLayout(obj)
            %Redo the layout of the figure.
            %There is a 150 pixel wide area on the left for gui components. The
            %rest of the figure is used for the image
            pos = get(obj.fig, 'Position');
            fig_width = pos(3) - 210;
            fig_height = pos(4) - 60;
            
            %set(obj.displayAxis, 'Position', [180 (pos(4) - fig_size - 30) fig_size fig_size]);
            set(obj.displayAxis, 'Position', [300 30 976*0.8 504*0.75]);
            top = pos(4) - 200;
            
            set(obj.movingWindowSpecPanel,'Position', [5/pos(3) 0.3 200/pos(3) 150/pos(4)]); 
           
            top = 100; 
            set(obj.wText, 'Position',[5 top 100 20]); 
            set(obj.windowWidth,'Position',[125 top 70 20]);top = top - 30;
            set(obj.arText, 'Position',[5 top 100 20]); 
            set(obj.aspectRatio,'Position',[125 top 70 20]);top = top - 30;
            set(obj.sText, 'Position',[5 top 100 20]); 
            set(obj.strideValue,'Position',[125 top 70 20]);top = top - 60;
            
            set(obj.updatePB, 'Position',[50 30 100 30]);
           
        end
        % -----------------------------------------------------------------       
        function xy_pt = convertRTpointToCartesian(obj, rt_pt, magnification, xy_size)  
            % Inputs:
           % rt_pt(1) = x index value in (r, theta) view
           % rt_pt(2) = y index value in (r, theta) view (i.e theta index)
           % In this example I assume the XY size is 704 (and is square) - I may want to use a variable
           % magnification: if no magnification required just use 1.0
           % outputs: the x,y pixel coordinates in XY (2 seperate columns)
                 
            %xy_size = 512; % 704 is the XY view image size (assumed square)
            
            m = magnification;
            dx = xy_size/2;
            dy = xy_size/2;
            dtheta = (2*pi)/obj.ImgHeight;
            
            conversionRatio = 0.5 * (xy_size/obj.ImgWidth);  
            
            
            rt_x = rt_pt(1);
            rt_y = rt_pt(2);
            if size(rt_pt,1)==3
                rt_z = rt_pt(3);
            end
           
            xy_radius = rt_x * conversionRatio;
            
            rt_theta = (rt_y-1) * dtheta;         % theta in (r,theta) view
            xy_theta = rt_theta - (pi/2.0); % theta in XY view           
             
            % Projection 
            Mprojection = [cos(xy_theta) 0 0; 0 sin(xy_theta) 0; 0 0 1];
            % Scaling
            Mscaling = [m 0 0; 0 m 0; 0 0 1];
            % Reflection of Y axis about the X axis
            Mreflection = [1 0 0; 0 -1 0; 0 0 1];
            % Translation to center of image
            Mtranslation = [1 0 dx; 0 1 dy; 0 0 1];
            
            xy_pt = Mtranslation * Mreflection * Mscaling * Mprojection * [xy_radius; xy_radius; 1];
            
             if size(rt_pt,1)==3
                 xy_pt(3) = obj.index * rt_pt(2)/obj.ImgHeight;
                 
             end
           
            
%             out = [];
%             for i=1:obj.ImgHeight
%                 theta = i * (2*pi)/obj.ImgHeight;
%                 
%                 pt = M*[r(i);r(i);1];   
%                 out(:, end+1) = pt;
%             end           
%             x = out(1,:);
%             y = out(2,:);
%             centerPoint(1) = sum(x(:))/size(x, 2);
%             centerPoint(2) = sum(y(:))/size(y, 2);
            
        end
        %------------------------------------------------------------------
        function [centerX_polar,centerY_polar,rho,theta_polar] = rect_to_polar_coordinate(obj, x,y,x_to_rho_factor,dimY)
            
            % convert coordinate in (r,theta) view image to coordinate in (x,y) view
            % image
            % rect means (r,theta) view image, polar means(x,y) view image
            
            % input: [x,y] is the pixel coordiate in (r,theta) view
            % output: [centerX_polar,centerY_polar] is the pixel coordinate in (x,y)
            % view
            % [rho,theta_polar] is the (rho,theta) cooridate in (x,y) view.
            % x_to_rho_factor = dimP/2/dimX, dimP is the dimension of (x,y) view image,
            % e.g. 1024, dimX is the number of columns of the (r,theta) view image,
            % e.g. 970
            % dimY is the number of rows of the (r,theta) image
            
            % Example
            % dimP = size(Img_polar,1);
            % [dimY, dimX] = size(Img_rect);
            % x_to_rho_factor = dimP/2/dimX;
            % [centerX_polar,centerY_polar,rho,theta_polar] = rect_to_polar_coordinate(x,y,x_to_rho_factor,dimY)
            
            rho = x*x_to_rho_factor; % convert x coordinate in (r,theta) view to radius (in number of pixels) in (x,y) view
            theta_rect = 2*pi/(dimY-1)*(y-1); % calculate theta value in (r,theta)view
            
            theta_polar = theta_rect - pi/2; % 0 in (r,theta) view corresponds -0.5pi in (x,y) view
            [centerX_polar,centerY_polar] = pol2cart(theta_polar,rho);
            
            %% origin is at (512,512) in the image
            centerX_polar = centerX_polar + 512;
            centerY_polar = 512 - centerY_polar;
        end
        % -----------------------------------------------------------
        % predicts the label of an a-line: it finds the label of the a-line
        % and the probability of that prediction.  The  
        % A-line's probability is selected to be the median of the
        % respective label's probability.  
        % Inputs: 
        % fCoords: m x 2 pairs of (r, theta) coordinates coming from the
        % classifier's output
        % y - label vector of the pixels
        function rv = aLinePredict(obj, fCoords, y, yProb, numOfLables,  typeOfLabeling)            
            
            
            %preferedMethod = 'label order';
            preferedMethod = 'heighest probability';
            
            
            switch preferedMethod
                case 'heighest probability' % in case I want a  label based on the heighest probability label in the line
                    rv = obj.probVote(fCoords, y, yProb, numOfLables);
                case 'label order' % if I want to set line a-label based on the order of labels along the a-line
                    %rv = struct('startPixels',{},'endPixels',{},'labelOrder',{});
                    rv = {};
                    for imgRow=1:obj.ImgHeight
                        rowIdx = (fCoords(:,2)==imgRow);
                        if sum(rowIdx) > 0
                            yOfRow = y(rowIdx)';
                            rv{imgRow} = obj.pickWidestLabelsAlongLine(yOfRow);
                        end
                    end
            end
            obj.aLineLabel = rv.label;
            obj.aLineLabelProb = rv.prob;
        end
        % -----------------------------------------------------------
        function rv = probVote(obj, fCoords, y, yProb, numOfLables)
            rv.label = zeros(obj.ImgHeight, 1);
            rv.prob = zeros(obj.ImgHeight, 1);
            
            for imgRow=1:obj.ImgHeight
                rowIdx = (fCoords(:,2)==imgRow);
                if sum(rowIdx) > 0
                    yOfRow = y(rowIdx)'; % make it a row just to emphasize that each element is a class of a pixel in the row
                    probOfRow = yProb(rowIdx, :)'; %row 1 is prob of label 1 etc
                    [rv.label(imgRow), rv.prob(imgRow)] = obj.determineAlineLabel(yOfRow, probOfRow, numOfLables, 'prediction');
                end
            end            
        end
        % -----------------------------------------------------------
        % given a row of labels (corresponding to a-line pixels), this
        % function resturns a strcture containg:
        % a vector specifying the order in which the labels appear along
        % the-line. and the start and end points of the corresponding
        % segments
        
        % This function finds the labes alongthe a-line and determines the
        % start and end pixel of each segment. 
        % (note that the function pickWidestLabelsAlongLine() does the same
        % operation only returns the widest one.
        
        function aLineLabels = pickWidestLabelsAlongLine(obj, labelsOfAline)
           
            labels = unique(labelsOfAline);
            labelOrder.s = zeros(1, length(labels));
            labelOrder.e = zeros(1, length(labels));
            
            for j=1:length(labels)
                iVals = find(labelsOfAline==labels(j));
                sVals = iVals([true diff(iVals)~=1]);
                eVals = iVals([diff(iVals)~=1 true]);
                dVals = eVals-sVals;
                [dValsSorted, dIdx] = sort(dVals);  % sort to find which segment is the largest within the plaque type
                                                                    % find only the areas were the plaque is the widest since I assume all of the
                                                                    % rest are noise or too short to be accounted for
                labelOrder.s(j) = sVals(dIdx(end));
                labelOrder.e(j) = eVals(dIdx(end));
                
            end
            [~, orderIdx] = sort(labelOrder.s); % sort to find the order of the labels along the line
            aLineLabels.labelOrder = labels(orderIdx);
            aLineLabels.startPixels = labelOrder.s(orderIdx);
            aLineLabels.endPixels = labelOrder.e(orderIdx);
            
        end
        % -----------------------------------------------------------
        % this function determines the label of an a-line based 
        % on the heighest median probablitlity of all lables in the line
        % It does it for two types of labels: 
        %   1. for the predicted values of the classifier based 
        % on the predicted labels and label probabilities.
        %   2. for the label image given by the expert:
         
        function [predictionLabel, predictionProb] = determineAlineLabel(obj, yOfRow, probOfRow, numOfLables, source)
            
            labels = unique(yOfRow);  % unique() also returns the index of the first occurance only
            switch source
                case 'prediction'
                    perLabelProb = zeros(length(labels),1);
                    for k=1:length(labels)
                        lableLoc = (yOfRow==labels(k));
                        %[prediction.label(imgRow), numOfAccurances] = mode(yOfRow(lableLoc)); % plurality vode
                        perLabelProb(k) = median(probOfRow(labels(k), lableLoc)); % consider taking the mean
                    end
                    [predictionProb, predictionLabel]  = max(perLabelProb);
                case 'expert'
                    nLabels = [labels(labels==Globals.CALCIUM) labels(labels==Globals.LIPID) labels(labels==Globals.FIBER)];
                    % method 1: take the majority vote out of the three : calcium,
                    % lipid, fibrous (not a good option)
                    % count the number of accurances of each label
                    
                    count = histc(yOfRow, nLabels); % get the count of labels
                    [~, maxVidx] = max(count);
                    predictionLabel = labels(maxVidx);
                    % and the probability of that specific plaque for all three : calcium,
                    % lipid, fibrous (not a good option)
                    predictionProbOfAll = count./sum(count); % compute the fraction (NOT PROBABILITY !!!) of the 3 main plaques
                    predictionProb = predictionProbOfAll(maxVidx);
             end
        end
        
        % -----------------------------------------------------------------
        % Here I separate the labels such that each label has its own
        % logical overlay images.  All in r-theta view
        function [binaryImages, y_predicted_revised] = getIndividualImages(obj, labelImg, coords, y_predicted, labelSet, applyMorphologicalFilteringFlag)
                    
%             y_predicted_revised = y_predicted; % [];
%             binaryImages = [];
%             return;
 
            actualLabelSet = unique(y_predicted);

            y_predicted_revised = [];
            % Separate the predicted labels into separate binalry
            % images
            
            [m, n] = size(labelImg.rt);
            binaryImages = zeros(m, n, Globals.NUMOFVOXELLABELS);
            for j = 1:size(coords,1)
                r = coords(j, 1);
                theta = coords(j, 2);
                switch y_predicted(j)
                    case Globals.CALCIUM
                        binaryImages(theta, r, Globals.CALCIUM) = 1;
                    case Globals.LIPID
                        binaryImages(theta, r, Globals.LIPID) = 1;
                    case Globals.FIBER
                        binaryImages(theta, r, Globals.FIBER) = 1;
                    case Globals.OTHER
                        binaryImages(theta, r, Globals.OTHER) = 1;
                    case Globals.NORMAL
                        binaryImages(theta, r, Globals.NORMAL) = 1;                    
                    case Globals.BKGD
                        binaryImages(theta, r, Globals.BKGD) = 1;
                end                
            end
             % Note that inspite of the fact that my labelSet has only 3
            % labels, the end results has 5 labels, so I need to show them
            % all
            %actualLabelSet = unique(y_predicted);
            showResult = 0;
            
            if showResult
                % First, display the expert-annotated label image
                figure;
                img_handle = imshow(labelImg.rt,[]);
                colormap(Globals.cmMatrix); set(img_handle,'CDataMapping','direct');
                % Now, show the separated binary images
                for i=actualLabelSet' % 1:numel(labelSet) 
                    %plaqueNum = actualLabelSet(i);
                    figure;
                    fig_handle = imshow(binaryImages(:,:,i),[]);
                    binaryColorMap = [0 0 0; Globals.RGBCOLORS(i,:)];
                    colormap(binaryColorMap);
                    title(Globals.plaqueLabels{i});
                    drawnow;
                end
            end
            
            % Apply a morphological operations on the fiber binary image
             % "openning" is erosion followed by dialation - removes peninsulas and islands, 
            % "closing" is dilation followed by erosion -  removes holes and fjords)
            %applyMorphologicalFilteringFlag = 0;
            % TO DO: NEED TO TRY MEDIAN FILTER ON THE COMBINED IMAGES (I.E.
            % NOT JUST THE BINARY IMAGE OF EACH PLAQUE TYPE. THIS WAY
            % SINGLE PIXELS OF ONE TYPE IN A NEIGHBERHOOD OF ANOTER TYPE 
            % WOULD BE SOMEHOW SMOOTHED.  IT IS EASIER TO EXPLAIN IN THE 
            % TEXT 
            if ~applyMorphologicalFilteringFlag
                y_predicted_revised = y_predicted;
            else 
                %se = strel('disk', 3, 0);                 
                for i=actualLabelSet'
                    binaryColorMap = [0 0 0; Globals.RGBCOLORS(i,:)];
                    operation = 'openThenClose'; %'open' 'close' 'openThenClose'
                    %
                    % +++ checkout: http://www.peterkovesi.com/matlabfns/
                    % he seem to have better morphological cleaning !!!
                    %[seg, Am, mask] = tcUtils.mcleanupregions(binaryImages(:,:,i), 0);
                    %
                    % bw2 = bwmorph(binaryImages(:,:,i), operation) ;
                    % figure, imshow(bw2);
                    % colormap(binaryColorMap);
                    %                 operation = 'close';
                    %                 bw = bwmorph(binaryImages(:,:,i), operation) ;
                    %                 figure, imshow(bw);
                    %                 colormap(binaryColorMap);
                    a = rand(5);
                    nhood = logical((a+a')/2>.5);
                    se = strel('arbitrary', nhood);
                    if strcmp(operation, 'open')
                        bw = imerode(binaryImages(:,:,i), se);
                        bw2 = imdilate(bw, se);
                    elseif  strcmp(operation, 'close')
                        bw =  imdilate(binaryImages(:,:,i), se);
                        bw2 = imerode(bw, se);
                    elseif  strcmp(operation, 'openThenClose')
                        bw = imerode(binaryImages(:,:,i), se);
                        bw2 = imdilate(bw, se);                        
                        bw3 =  imdilate(bw2, se);
                        bw4 = imerode(bw3, se);
                        bw2 = bw4;
                    end
                    figure, imshow(bw2);
                    colormap(binaryColorMap);
                    drawnow;
                    binaryImages(:,:,i) = bw2;
                    
                end
                
                %% now update the acuracy numbers to reflect the improved images
                y_predicted_revised = zeros(size(y_predicted,1),1);
                for j = 1:size(coords,1)
                    r = coords(j, 1);
                    theta = coords(j, 2);
                    if binaryImages(theta, r, Globals.CALCIUM) == 1
                        y_predicted_revised(j) = Globals.CALCIUM;
                    elseif binaryImages(theta, r, Globals.LIPID) == 1
                        y_predicted_revised(j) = Globals.LIPID;
                    elseif binaryImages(theta, r, Globals.FIBER) == 1
                        y_predicted_revised(j) = Globals.FIBER;
                    elseif binaryImages(theta, r, Globals.NORMAL) == 1
                        y_predicted_revised(j) = Globals.NORMAL;
                    elseif binaryImages(theta, r, Globals.OTHER) == 1
                        y_predicted_revised(j) = Globals.OTHER;
                    elseif binaryImages(theta, r, Globals.BKGD) == 1
                        y_predicted_revised(j) = Globals.BKGD;
                    end
                end
                
                %             bw = imerode(binaryImages(:,:,Globals.FIBER), se);
                %             bw2 = imdilate(bw, se);
                %             figure, imshow(bw2);
                %             binaryColorMap = [0 0 0; Globals.RGBCOLORS(Globals.FIBER,:)];
                %             colormap(binaryColorMap);
            end
            
                       
            obj.showSectorPluralityVote(actualLabelSet, binaryImages);
            
        end
        % -----------------------------------------------------------------
        function computeSectorStat(obj, binaryImages, trueLablesImg)
            
            numOfSectors = 4;
           
            sectorThetaDelta = floor(obj.ImgHeight/numOfSectors);
            xyImageSize = 512;
            
            %labelSet = [1 2 3]; % I deliberatly do not count the "other" since it's always the majority
            
            estimatedNumOfVoxels = zeros(numOfSectors, size(binaryImages,3));
            trueNumOfVoxels = zeros(numOfSectors, size(binaryImages,3));
            pctOverlap = zeros(numOfSectors, size(binaryImages,3));
            dice = zeros(numOfSectors, size(binaryImages,3));
            
            if obj.index==240 % when I used vessel 21, I wanted only one frame's result (255)
                
                for i=1:size(binaryImages,3) % numel(labelSet)
                    I = binaryImages(:,:,i);
                    trueLabelsI = trueLablesImg.rt==i;
                    figure;
                    imshow(I,[]);
                    binaryColorMap = [0 0 0; Globals.RGBCOLORS(i,:)];
                    colormap(binaryColorMap);
                    title(Globals.plaqueLabels{i});
                    figure;
                    for k=1:5
                        trueLabelsI= medfilt2(trueLabelsI);
                    end
                    imshow(trueLabelsI(:,:,i),[]);
                    colormap(binaryColorMap);
                    drawnow;
                    for s=1:numOfSectors
                        sectorMask = bsxfun(@eq, I((s-1)*sectorThetaDelta+1:s*sectorThetaDelta, :), 1);
                        estimatedNumOfVoxels(s, i) = sum(sum(sectorMask));
                        trueLabelsMask = trueLabelsI((s-1)*sectorThetaDelta+1:s*sectorThetaDelta, :);
                        trueNumOfVoxels(s, i) = sum(sum(trueLabelsMask));                        
                        dice(s) = 2*nnz(sectorMask&trueLabelsMask)/(nnz(sectorMask) + nnz(trueLabelsMask));
                        if estimatedNumOfVoxels(s,i)
                            pctOverlap(s, i) =  trueNumOfVoxels(s,i)/estimatedNumOfVoxels(s,i);
                        else
                            pctOverlap(s, i) = 0;
                        end
                        hold on; plot([1 984], [s*sectorThetaDelta s*sectorThetaDelta], '--y');                        
                    end                    
                end
                fprintf('\t Calcium \t\t Lipid \t\t Fiber \t\t Other\n');
                disp(estimatedNumOfVoxels);
                disp(pctOverlap);
                disp(dice);
                
                %% Now, show results within a range of r (imitate: if clinician want to see all calcium thicker than some value
                showResultInRange = 0;
                if showResultInRange
                    dt = 1.9; % desired thickness in mm in XY view
                    % convert to r-theta
                    dr = floor((dt/Globals.XY_PIXEL_SIZE)/2); % the distance in r-theta view in number of pixels (roughly half of XY view)
                    [Imasked, mask] = obj.extractBloodVesselMask();
                    % modify the mask to show only the thicker part of the
                    % calcium
                    lumenVector = obj.getLumenBorder(); % I call the function since in there I check if the lumen exist
                    calciumMask = binaryImages(:,:,1);
                    rightMask = bsxfun(@ge, 1:size(calciumMask, 2), min((lumenVector+dr), obj.ImgWidth)); % compare row vector 1:size(mask, 2) to each element of lumenVector and return 1 when it is greater or equal (ge)
                    revisedMask = calciumMask .* uint16(rightMask);
                    obj.displayImage(Imasked, 'xy', xyImageSize);hold on;
                    plot(floor(xyImageSize/2), floor(xyImageSize/2), '+r');
                    obj.displayImage(revisedMask, 'xy', xyImageSize);
                    binaryColorMap = [0 0 0; Globals.RGBCOLORS(1,:)];
                    colormap(binaryColorMap);
                end
                
                %% Now, show the plurality vote using colored sector arcs
                [Imasked, mask] = obj.extractBloodVesselMask();
                %Imasked(Imasked==0) = 65535;
                obj.displayImage(Imasked, 'xy', xyImageSize);hold on;
                %obj.displayImage(obj.data, 'xy', xyImageSize); hold on;
                hold on;
                plot(floor(xyImageSize/2), floor(xyImageSize/2), '+y');
                
                siSize.w = 70;
                siSize.h = sectorThetaDelta;
                radius = obj.ImgWidth - siSize.w - 1;
                ulhc = [radius 1]; % [r,theta] in r-theta view
                animateView = 'xy';
                tempC = 'bbbgggbbbw';
                for s=1:numOfSectors
                    [~, majorityPlaque] = max(estimatedNumOfVoxels(s,:));
                    options.fillColor = Globals.COLORS(majorityPlaque);
                    %options.fillColor = tempC(s);
                    shapeHandle = obj.plotSI(ulhc, siSize, animateView, xyImageSize, options); %
                    hold on;
                    ulhc = [radius s*sectorThetaDelta];
                    hold on;
                    %plot([256 352], [256 161], '--r')
                end
            end
        end
        % -----------------------------------------------------------------
        function showSectorPluralityVote(obj, binaryImages)
            
            numOfSectors = 8;
           
            sectorThetaDelta = floor(obj.ImgHeight/numOfSectors);
            xyImageSize = 512;
            
            %labelSet = [1 2 3]; % I deliberatly do not count the "other" since it's always the majority
            
            numOfVoxels = zeros(numOfSectors, size(binaryImages,3));
            
            for i=1:size(binaryImages,3) % numel(labelSet)
                I = binaryImages(:,:,i);
                for s=1:numOfSectors                    
                    sectorMask = bsxfun(@eq, I((s-1)*sectorThetaDelta+1:s*sectorThetaDelta, :), 1);
                    numOfVoxels(s, i) = sum(sum(sectorMask));
                    
%                     figure;
%                     fig_handle = imshow(binaryImages(:,:,i),[]);
%                     binaryColorMap = [0 0 0; Globals.RGBCOLORS(i,:)];
%                     colormap(binaryColorMap);
%                     title(Globals.plaqueLabels{i});
%                     drawnow;
                end                
            end
            fprintf('\t Calcium \t\t Lipid \t\t Fiber \t\t Other\n');
            disp(numOfVoxels);
            
            [Imasked, mask] = obj.extractBloodVesselMask();
            %Imasked(Imasked==0) = 65535;
            obj.displayImage(Imasked, 'xy', xyImageSize);hold on;            
            %obj.displayImage(obj.data, 'xy', xyImageSize); hold on;
            hold on;
            plot(floor(xyImageSize/2), floor(xyImageSize/2), '+y');
            
            siSize.w = 70;
            siSize.h = sectorThetaDelta;
            radius = obj.ImgWidth - siSize.w - 1;
            ulhc = [radius 1]; % [r,theta] in r-theta view
            animateView = 'xy'; 
            tempC = 'bbbgggbbbw';
            for s=1:numOfSectors    
                [~, majorityPlaque] = max(numOfVoxels(s,:));
                options.fillColor = Globals.COLORS(majorityPlaque); 
                %options.fillColor = tempC(s);
                shapeHandle = obj.plotSI(ulhc, siSize, animateView, xyImageSize, options); % 
                hold on;
                ulhc = [radius s*sectorThetaDelta]; 
                hold on; 
                %plot([256 352], [256 161], '--r')
            end
        end        
        % --------------------------------------------------------
        function outImg = fillROIs(obj, xyImg, xy, bwMask, numOfROIs)
            
            figure;           
            img_handle1 = imshow(xyImg, []); hold on;
            colormap(Globals.cmMatrix); set(img_handle1,'CDataMapping','direct');
            %title('Original cleaned image');
%             subplot(2, 2, 2);
%             img_handle2 = imshow(bwMask{1}); hold on;
%             title('Mask');
%             subplot(2, 2, 3);
%             img_handle3 = imshow(bwMask{2}); hold on;
%             title('Mask');
            border = {};
            for i=1:numOfROIs
                structBoundaries = bwboundaries(bwMask{i});
                border{i}=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
                hold on; % Don't blow away the image.
                plot(border{i}(:, 2), border{i}(:, 1), 'w-', 'LineWidth', 1);
            end
            
            
            % Burn the marked regions into the image
            outImg = xyImg;
            outImg(bwMask{1}) = Globals.CALCIUM;
            outImg(bwMask{2}) = Globals.CALCIUM;
            outImg(bwMask{3}) = Globals.FIBER;
            figure;
            img_handle2 = imshow(outImg, []); hold on;
            colormap(Globals.cmMatrix); set(img_handle2,'CDataMapping','direct');
            %title('Image with burned VOIs');
            
            % Now, I need to add noise to the roi's selected by the
            % freehand
 
            for n=1:numOfROIs
                for iter=1:10
                    for i = 1:size(border{n},1)
                        x = border{n}(i,2) + randi([-10,10],[1,1]); % generate a vector of random integers between -10 and 10.
                        y = border{n}(i,1) + randi([-15,20],[1,1]); % generate a vector of random integers between -15 and 20.
                        if n==3
                            outImg(y, x) = Globals.FIBER;
                        else
                            outImg(y, x) = Globals.CALCIUM;
                        end
                    end
                end
            end
            figure;
            img_handle3 = imshow(outImg, []); hold on;
            colormap(Globals.cmMatrix); set(img_handle3,'CDataMapping','direct');
            
            for i=1:numOfROIs
                hold on; % Don't blow away the image.
                plot(border{i}(:, 2), border{i}(:, 1), 'w-', 'LineWidth', 1);
            end
        end
        %--------------------------------------------------------
        function finalLabel = determineSingleLabel(obj, aLineLabelVector, singleLine)
            % method 1: take the majority vote out of the three : calcium,
            % lipid, fibrous (not a good option)
                       
            % count the number of accurances of each label
            labels = unique(singleLine);  % gives the index of the first occurance only
            nLabels = [labels(labels==Globals.CALCIUM) labels(labels==Globals.LIPID) labels(labels==Globals.FIBER)];
            count = histc(singleLine, nLabels); % get the count of labels
            [~, maxVidx] = max(count);
            finalLabel = labels(maxVidx);
            
            % method 2: return the probability for all three : calcium,
            % lipid, fibrous (not a good option)
            finalLabelProb = count./sum(count); % compute the fraction of the 3 main plaques
        end
        % -----------------------------------------------------------------
        % given a row of labels (corresponding to a-line pixels), this
        % function returns a strcture containg:
        % a vector specifying the order in which the labels appear along
        % the-line. and the start and end points of the corresponding
        % segments        
        % (note that the function pickWidestLabelsAlongLine() does the same
        % operation only returns the widest one.
        function aLineLabels = evaluateLabelsLengths(obj, aLine)
            
            aLineLabels = struct('startPix',{},'endPix',{},'label',{});
            
            % count the number of accurances of each label
            labels = unique(aLine);  % gives the index of the first occurance only
            count=histc(aLine,labels); % get the count of elements
            
            counter = 1;
            for i=1:length(labels)
                iVals = find(aLine==labels(i));               
                sVals = iVals([true diff(iVals)~=1]);
                eVals = iVals([diff(iVals)~=1 true]);
                labelVals = ones(1, size(sVals,2)) * labels(i);
                
                for j=1:length(sVals)
                    aLineLabels(counter).startPix = sVals(j);
                    aLineLabels(counter).endPix = eVals(j);
                    aLineLabels(counter).label = labelVals(j);  
                    counter = counter + 1;
                end                
            end
        end
        % --------------------------------------------------------
        % inputs:
        % predictedLabels : a cell array of the a-line predicition (i.e
        % a-line-based prediction) 
        % trueLabels: a cell array of the a-line expert labeling I used to
        % compare my prediction)
        % the two vectors have different structure depending on the value of
        % typeOfLabeling
        %
        % ASSUMPTION: THE TWO VECTORS HAVE THE SAME LENGTH WHICH IS THE
        % HEIGHT OF THE FRAME !!
        %
        % Notice that the predictedLabels will contain only the labels that
        % appear in the training data, yet the trueLabels will contain all
        % labels marked by the expert.  In the case of typeOfLabeling='single
        % value' then I do it  in determineAlineLabel()           
           
        function matchingLines = analyzeAlinePrediction(obj, predictedLabels, trueLabels, typeOfLabeling)            
           
            matchingLines = zeros(size(trueLabels,2), 1);
            labelsOfInterest = [Globals.CALCIUM Globals.LIPID Globals.FIBER];
            
            id = false(1, length(labelsOfInterest));            
               
            switch typeOfLabeling
                case 'label order'
                    useMethod = 1;
                    if useMethod == 1
                        % one approach: look what labels appear on the
                        % trueLabels and see if I have the same labels in the
                        % same order in the predicted.
                        for i=1:size(trueLabels,2)
                            % I have no idea how the following 7 lines work
                            A = trueLabels{i}.labelOrder;
                            B = predictedLabels{i}.labelOrder;
                            n = length(A);
                            m = length(B);
                            K = (1:n-m+1)' * ones(1,m) + ones(n-m+1,1) * (0:m-1);
                            K(:) = A(K);
                            J = find(all(K'==(B(:)*ones(1,n-m+1))));
                            if ~isempty(J)
                                matchingLines(i) = true;
                            end                            
                        end                    
                        % Second approach
                    elseif useMethod==2
                        TP = 0;
                        for i=1:size(trueLabels,2)
                            id = [];
                            for j=1:length(labelsOfInterest)
                                id(j, :) = trueLabels{i}.labelOrder == labelsOfInterest(j);
                            end
                            idx = false(1, size(id,2));
                            for k=1:size(id,1)
                                idx(1, :) = idx(1,:) | id(k,:);
                            end
                            if sum(idx)>1
                                disp('');
                            end
                            tl = trueLabels{i}.labelOrder(idx);
                            if strfind(tl, predictedLabels{i}.labelOrder)
                                TP = TP + 1;
                            end                            
                        end
                    end
                case 'single value'
                    
                    
            end
        end
        % ---------------------------------------------------------
    end % methods
    
end

