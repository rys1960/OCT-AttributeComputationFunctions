function mainComputeFeatures()
   close all; close all hidden;   
   
   ui = FromGUI(); % for now, it is all hardcoded, later need to build a GUI
   ui.setGUIvalues();
   %    numOfCPUs = getenv('NUMBER_OF_PROCESSORS');
   %    if matlabpool('size') == 0
   %        matlabpool open
   %    end   
  
   
   globals = Globals;   % Initilaize the Globals
   
   [fns, frameList] = ui.getFileNames(); % get all of the file names to work with including the frames of interest
   pb = Pullback(fns.pb);
   X = pb.computeSectorAttributes(ui, frameList, fns);
   
   %% --------- Perform Exploratory Data Analysis to check assumptions ----------
   doExploratoryAnalysis = 1;
   if doExploratoryAnalysis
       exploratoryAnalysis = ExploratoryDataAnalysis();
       exploratoryAnalysis.checkFeatureDistribution(X, ui);
       exploratoryAnalysis.checkFeatureCorrelation(X, ui);      
   end

   disp('....Done!');
end