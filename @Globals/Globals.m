classdef Globals < handle
    
    
%      enumeration
%          CC, LL, FIB, NOR
%      end
    properties (Constant)
        
        CALCIUM = 1; LIPID = 2;  FIBER = 3;
        OTHER = 4; NORMAL = 5;   BKGD = 6; 
        
        RGBCOLORS = [[1.0, 0.0, 0.0]; [0.0, 1.0, 0.0]; [0.0, 0.0, 1.0]; [0.5, 0.5, 0.5]; [0.0, 0.0, 0.0]; [0.0, 0.0, 0.0]];
        cmMatrix = [[0 0 0];...
                    Globals.RGBCOLORS(1,:);... % red % REMEMBER THAT 0 MAPS TO THE FIRST ENTRY, SO i INITIALIZE AN IMAGE WITH BKGD
                    Globals.RGBCOLORS(2,:);... % green
                    Globals.RGBCOLORS(3,:);... % blue
                    Globals.RGBCOLORS(4,:);... % Other
                    Globals.RGBCOLORS(5,:);... % magenta
                    Globals.RGBCOLORS(6,:)];% background                                                
                          
        plaqueLabels = {'Calcium','Lipid','Fibrous','Other', 'Normal', 'Background'};
        
        % Catheter parameters (as derived by Madhu-names are per Van Soest paper)
        Zc = 0.0;            % A machine-specific value taken from the manual
        Zw = 12.0;         % A machine-specific value taken from the manual
        Z0 = 1.0565;      % Computed using non-lin fit using the water pullback (see mainFindCathParams.m)
        Zr = 0.59023;     % Computed using non-lin fit using the water pullback

        %BASELINE = [6.0+3.6092, 6.0+3.6092, 6.0+3.6092]; % analyzed in a seperate utility mainPreProcessingAnalysis.m
        BASELINE = 6.0;
        BASELINE_STD = 3.6092;
             
        SCAN_WIDTH = 300; % 400 INDEX VALUES WHICH IS 2mm beyond the lumen
        FRAME_THICKNESS = 54/271; % ~200 microns
        COLORS = 'rgbckmy';
        SHAPES = '+ov^sd*';        
        
        % Some machine constants
        RTHETA_PIXEL_SIZE = 4.9587 * 1e-3; % 5 microns converted to mm
        XY_PIXEL_SIZE = 4.8/512; % in mm
            
        NUM_LINES_PER_FRAME = 504;
        NUM_COL_PER_FRAME = 976;        
                
        REGIONALATTRIBUTENAMES = {'FrameNum',...
            'Theta',...
            'lumenEccentricity',...
            'lumenArea',...
            'lumenCircularity',...
            'AvgI',...
            'smoothness',...
            'uniformity',...
            'homogeneity',...
            'contrast',... 
            'entropy',...
            'BeamDepth',...
            'initimalArea',...
            'MediaDepth',...
            'mediaArea',...
            'mediaDepth',...
            'AdArea',...
            'AdDepth'};
    end   
end
