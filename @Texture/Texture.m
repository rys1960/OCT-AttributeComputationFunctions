classdef Texture < handle 
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties  %(GetAccess=private, SetAccess=private)  
        avg;
        stdDev;
        smoothness;
        uniformity;
        homogeneity;
        contrast;
        entropy;
        
        % additional variables when using co-occurance matrix
        maxProbability;
        correlation;
        energy;
        
    end
    
    methods
        function obj = Texture(v, which, winROI) % Constructor
            
            if strcmp(which, 'mine')
                v = v(:); % vectorize the data
                obj.avg = mean(v);
                obj.stdDev = std(double(v));
                
                [nn, xx] = hist(double(v), 5);
                p = nn/sum(nn);
                
                obj.smoothness = 1 - 1/(1+obj.stdDev^2); % Measure of smoothness (zero for constant intensity and 1 for regions for large excursion)
                obj.uniformity = sum(p .* p);
                obj.homogeneity = obj.stdDev/obj.avg;
                obj.contrast = obj.stdDev;    % std as a measure of contrast
                obj.entropy = -sum(p .* log2(p+0.001)); % entropy(v);
            
            elseif strcmp(which, 'ughi')
                %obj.correctForCatheter();
                %bb = obj.getBoundingBox('pixels'); % vector [min_theta, max_theta, min_r, max_r]
                subI = uint16(winROI);
                g = graycomatrix(subI, 'NumLevels', 64); % Ughi says he selected 64
                g_n = g/sum(g(:)); % Normalized matrix.
                statTexture = graycoprops(g, 'all'); % Descriptors.
                obj.contrast = statTexture.Contrast;
                obj.energy = statTexture.Energy;
                obj.correlation = statTexture.Correlation;
                obj.homogeneity = statTexture.Homogeneity;
                statTexture.maxProbability = max(g_n(:));
                obj.maxProbability = statTexture.maxProbability;
                for i = 1:size(g_n,1)
                    sumcols(i) = sum(-g_n(i,1:end).*log2(g_n(i,1:end) + eps));                    
                end
                statTexture.entropy = sum(sumcols);
                obj.entropy = statTexture.entropy;
                
            end
        end
        % -------------------------------------------------------
        function rv = getMeasures(obj, which)
            if strcmp(which, 'mine')
                rv.smoothness = obj.smoothness;
                rv.uniformity = obj.uniformity;
                rv.homogeneity = obj.homogeneity;
                rv.contrast = obj.contrast;    
                rv.entropy = obj.entropy;
                rv.avg = obj.avg;
            elseif strcmp(which, 'ughi')
                rv.maxProbability = obj.maxProbability;
                rv.entropy = obj.entropy;
                rv.contrast = obj.contrast;
                rv.energy = obj.energy;
                rv.correlation = obj.correlation;
                rv.homogeneity = obj.homogeneity;
            end
        end
    end
    
end

