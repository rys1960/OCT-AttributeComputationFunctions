classdef Sector < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here    
     
    properties %(GetAccess=private, SetAccess=private)  
   
        origin;     % [theta, r, frame]
        height;
        depth;
        lineCollection;
    end
    methods
        function obj = Sector(secMeas)    % Constructor           
            obj.origin = secMeas.origin;
            obj.height = secMeas.height;
            obj.depth = secMeas.depth;
            
        end %constructor        
         % ---------------------------------------------------------------  
         function [y, ySequence] = computeAnnotation(obj, theta, h, lumenVector, backBorderVector, labelImg)
             % Define the sector's min and max theta in the frame
             lowTheta = max(1, theta-floor(obj.height/2));
             hiTheta = min(h, theta+floor(obj.height/2));
             % acc will hold the number of occurances of each label
             acc = zeros(numel(Globals.plaqueLabels), 1);
             
             numRows = hiTheta-lowTheta+1;
             rowCounter = 1;
             
             %Ys = zeros(numRows, 1); 
             Ys = []; 
             YsSequenceVector = {}; 
             for t=lowTheta: hiTheta                 
                 oneRow = labelImg(t, lumenVector(t):backBorderVector(t));
                 % fix transformation errors: make sure the first layer is fibrous
                 oneRow(1:20) = Globals.FIBER;                 
                             
                 or = obj.cleanRow(oneRow); % For row label, I remove labels that are too short
                 
                 % get the order of label occurances to generate the line's
                 % class (based on
                 % http://blogs.mathworks.com/loren/2009/11/26/unique-values-without-rearrangement/)     
                 [Xs, ix] = sort(or(:));
                 uv(ix) = ([1; diff(Xs)] ~= 0);
                 order = or(uv);
               
                 [Ys(end+1), YsSequenceVector{end+1}] = obj.findSingleRowSequence(order);                
                 
                 % Accumulate the number of pixels of each plaque
                 % in case I need to use it later on in the annotation
                 %vv = tabulate(oneRow);
                 values = unique(oneRow);
                 [instances, binIdx] = histc(oneRow(:),values);
                 for j=1:numel(values)
                     acc(values(j)) = acc(values(j)) + instances(j);
                 end                 
                 rowCounter = rowCounter + 1;
             end
             
             y = mode(Ys); 
             % Find which letter sequence this label is
             idx = find(Ys==y);
             ySequence = YsSequenceVector(idx(1));
             dispResultsFlag = 0;
             if dispResultsFlag
                 hold on;
                 imshow(labelImg,[]);
                 hold on;
                 plot(backBorderVector, [1:size(backBorderVector,1)], 'y-', 'LineWidth', 3);
                 plot(lumenVector, [1:size(lumenVector,1)], 'r-', 'LineWidth', 3);
                 plot(lumenVector(lowTheta):backBorderVector(lowTheta), lowTheta, 'g-')
                 plot(lumenVector(hiTheta):backBorderVector(hiTheta), hiTheta, 'r-')
             end
         end
         % -----------------------------------------------------         
         function [y, ySequence] = findSingleRowSequence(obj, order)
             
             y = 0;
             ySequence = char.empty();
             allClasses = [1:Globals.NUMOFVOXELLABELS];
             p = perms(allClasses);
             ind = p(:,1)~=3;
             p(ind,:) = [];             
             
             for i=1:size(p,1)
                 for j=1:size(p,2)
                     if isequal(order, p(i,1:j))                        
                         %y = num2str(order);
                         %y = order;
                         ySequence = Globals.LETTERS(order);                        
                         break;
                     end
                 end                 
             end             
             
             lastOne = numel(Globals.REGIONCLASSES);
             if strcmp(ySequence, Globals.REGIONCLASSES{lastOne})
                 y = lastOne;
             elseif strcmp(ySequence, Globals.REGIONCLASSES{lastOne-1})
                 y = lastOne - 1;
             else                 
                 for j=1:lastOne-2
                     if strfind(ySequence, Globals.REGIONCLASSES{j})
                         y = j;
                         break; % just to save time
                     end
                 end
             end
             if ~y
                 disp('Label not defined for sequence'); disp(order);
             end
         end
         % -----------------------------------------------------
         function or = cleanRow(obj, oneRow)
             
             
             or = oneRow;
             originalMembers = unique(oneRow);
             notMembers = setdiff([1:Globals.NUMOFVOXELLABELS], originalMembers);
             %subplot(131);plot(or); ylim([0 6]);
             
             win = 30;
             temp=zeros(1, floor(win/2)+ numel(or));
             temp(1:floor(win/2)) = or(1);
             temp(floor(win/2)+1:end) = or;
             temp = medfilt1(temp, 30);
             or = temp(floor(win/2)+1:end);
             or = round(or);
             %subplot(132); plot(or); ylim([0 6]);
             % Now, remove those values that were not in the vector
             % originaly
             for k=1:numel(notMembers)
                 loc = find(or==notMembers(k));
                 or(loc) = or(min(loc+floor(win/2), end)); % need to assign the closest one, but for now leave arbitrary assignment
             end
             %subplot(133); plot(or); ylim([0 6]);
         end
         % -----------------------------------------------------
         function setLineAggregation(obj, in_lineCollection)
             obj.lineCollection = in_lineCollection;
         end
         % --------------------------------------------------------
        function statTexture = computeTextureAttributes(obj) 
           textureMeasures = Texture(obj.lineCollection, 'mine');
           statTexture = textureMeasures.getMeasures('mine');            
         end
        % ------------------------------------------------------
        % See Sonka 3rd edition in champter 8 for region based features'
        % definition
        function geom = computeGeometricAttributes(obj, theta, h, lumenVector, wallLayers, backBorderVector)
            % First define the sector's min and max theta in the specific
            % frame for which I have the lumenShape, wallLayers and intimaLayerProperties
            lowTheta = max(1, theta-floor(obj.height/2));
            hiTheta = min(h, theta+floor(obj.height/2));
            
            %% Now, calculate all of the features
            % Calculate the average beam penetration depth
            % over the lines of the sectors of all layers.  The area is
            % calculated as the sum of all the pixels between the
            % respective borders (each pixel is 1)
            
            e = ones(hiTheta-lowTheta+1,1);
            numElements = numel(e);
            u = lumenVector(lowTheta:hiTheta);
            v = backBorderVector(lowTheta:hiTheta);
            w1 = wallLayers.iem(lowTheta:hiTheta);
            w2 = wallLayers.eem(lowTheta:hiTheta);
            w3 = wallLayers.adventitia(lowTheta:hiTheta);
            
            % Area are calculated in units of pixels
            geom.d_depth = ( (v-u)'*e )/numElements;
            geom.intimalArea = (w1-u)'*e;
            geom.intimalDepth = geom.intimalArea/numElements;
            geom.mediaArea = (w2-w1)'*e;
            geom.mediaDepth = geom.mediaArea/numElements;
            geom.adventitiaArea = (w3-w2)'*e;
            geom.adventitiaDepth = geom.adventitiaArea/numElements;
            
           
         end        
        %------------ Accessors ----------------------------------------
        function rv = getNumberOfAlines(obj)
            rv = obj.numberOf_ALines;           
        end
        %-------------------------------------------------------------------
        function rv = getSliceNumbers(obj)
            rv = obj.sliceNumbers;           
        end
        %-------------------------------------------------------------------
        function rv = getNumOfSubvolumes(obj)
            rv = size(obj.sliceNumbers, 3);           
        end
        % ----------------------------------------------------------------
        function rv = getDlumen(obj)
                rv = obj.Dlumen;                
        end
        % ----------------------------------------------------------------
        function rv = getDcath(obj)
                rv = obj.Dcath;                
        end        
        %-------------------- Modifiers -------------------------------
        function setSliceNumbers(obj, vec)
            obj.sliceNumbers = vec;
        end
        % ---------------------------------------------------------------
    end
end
