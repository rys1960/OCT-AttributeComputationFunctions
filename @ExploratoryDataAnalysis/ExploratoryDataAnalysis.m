classdef ExploratoryDataAnalysis <handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    % ---------------------------------------------------------------------
    methods
        function obj = ExploratoryDataAnalysis()
        end % constructor
        %------------------------------------------------------------------
        function checkFeatureDistribution(obj, X, guiOptions)
            % See whether the features are normally distribution
           % (If the features aren't normally distributed, I shouldn't use
           % discriminant analysis)
           figure;
           if strcmp(guiOptions.getAnalysisType(), 'voxel')
               for i = 1:size(X, 2)
                   subplot(3,4,i)
                   %The purpose of a normal probability plot is to graphically assess whether the data in X could come from a
                   % normal distribution. If the data are normal the plot will be linear. Other distribution types will introduce curvature in the plot.
                   normplot(double(X(:,i))) %
                   title(guiOptions.getAttName(i));
                   drawnow;
               end
           elseif strcmp(guiOptions.getAnalysisType(), 'sector')
               for i = 1:size(X, 2)
                   subplot(4,5,i)
                   %The purpose of a normal probability plot is to graphically assess whether the data in X could come from a
                   % normal distribution. If the data are normal the plot will be linear. Other distribution types will introduce curvature in the plot.
                   normplot(double(X(:,i))) %
                   title(Globals.REGIONALATTRIBUTENAMES{i});
               end
           end
           
        end
        % -----------------------------------------------------------------
        function checkFeatureCorrelation(obj, X, guiOptions)
            % See whether the features are correlated with one another.
            % (If the features are highly correlated we shouldn't use
            % Naive Bayes Classifier) nor should I use the higly correlated
            % features since they contribute nothing
            covmat = corrcoef(double(X));
            
            figure
            x = size(X, 2);
            imagesc(covmat);
            set(gca,'XTick',1:x);
            set(gca,'YTick',1:x);
            if strcmp(guiOptions.getAnalysisType(), 'voxel')
                set(gca,'XTickLabel',guiOptions.getAttName('all'));
                set(gca,'YTickLabel',guiOptions.getAttName('all'));
            elseif strcmp(guiOptions.getAnalysisType(), 'regional')
                set(gca,'XTickLabel',Globals.REGIONALATTRIBUTENAMES);
                set(gca,'YTickLabel',Globals.REGIONALATTRIBUTENAMES);
            end
            title('Correlation among features');
            axis([0 x+1 0 x+1]);
            grid;
            colorbar;
            
        end
        % -----------------------------------------------------------------
        function showLabelImage(obj, y)
            figure;
            imagesc(y); 
            title('Class labels');            
        end
        % -----------------------------------------------------------------
        function showFeatureImage(obj, X)
            figure;            
            imagesc(X); 
            title('Attributes');
        end
        % *****************************************************************
    end
end

