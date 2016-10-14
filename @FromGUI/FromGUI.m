classdef FromGUI < handle    
       
    properties (Constant)
        C1 = 9; % just for testing
        C2 = 10;
    end
    % ---------------------------------------------------------------
    properties (GetAccess=private)        % (GetAccess=private, SetAccess=private)  
        analysisType;
        sector; % in r-t view
        
    end
    % =====================================
    methods
        function obj = FromGUI()    % Constructor         
            
        end % Constructor        
        % -----------------------------------------------------------------
        function setGUIvalues(obj) % set all values at once
            obj.analysisType = 'sector'; % 'voxel' 'sector' 'HOG' 'Texton'
            % sector sized for regional-based analysis
            obj.sector.origin = [1;1;1]; % [theta, r, frame] of ulhc of the region in r-theta view
            obj.sector.height = 63; % number of lines (in r-theta view) the sector covers
            obj.sector.depth = 1; % number of frames along Z.  Can be even or odd
            obj.sector.stride = obj.sector.height+1;
        end
        
        % -----------------------------------------------------------------
        function rv = getSectorMeasures(obj)
            rv = obj.sector;
        end
        % ------------------------------------------------------
        function rv = getAnalysisType(obj)
            rv = obj.analysisType;
        end
        % -----------------------------------------------------------------  
        % fns : structure containing the file names
        % frameList : a vector of teh frames of interest (good when
        % debugging so I do not process all frames every run.
        function [fns, frameList] = getFileNames(obj)
            
            datadir = getenv('OCTDATA_DIR'); % need to set environment variable indicating where data is located
            
            fns.pb = [datadir 'ValidationData\Prabu\Vessel 52\OCT Raw\distalLAD.oct'];
            frameList = [346:436+46]; %: [346:346+64];  % the first frame Prabhu gave me is 346
            fns.trueLabels = cell(numel(frameList), 1); % the file of the true label image
            fns.attributes = cell(numel(frameList), 1);
            
            trueLabelsRoot = [datadir 'ValidationData\Prabu\Vessel 52\OCT Labels\Vessel52_Segment1_labels.'];
            
            fns.offset = 346;        % the offset between the first frame prabhu annotated and the file numbering
            
        end
    end
end
