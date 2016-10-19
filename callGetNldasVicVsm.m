% This script will call the getNldasVsm function
% P. Shellito
% 7/29/16

clear all
close all

% -------------------------------------------------------------------------
% In this example, there are sites with latitude and longitude in the
% following input file:
inFile = './inFile.txt';
% Date range requested
qStart = [2015,3,31];
% qStart = 'apnd';
qEnd = [2015,4,1];

% -------------------------------------------------------------------------
% Read the input data from the text file
% Open the input file
fid = fopen(inFile);
% Read the data in the input file
data = textscan(fid,'%s\t%f\t%f', 'headerlines', 1);
% Close the input file
fclose(fid);

% -------------------------------------------------------------------------
% Organize input data and create the needed vectors to pass into the function
% A cell array of strings
qNames = data{1,1};
% Latitude of the sites in qNames
qLat = data{1,2};
% Longitude of the sites in qNames
qLon = data{1,3};
% Directory to hold the output text files
outDir = './nldasVsm';

% -------------------------------------------------------------------------
% Record what time is is before the function is called
disp('Starting the script at')
disp(datetime)
startTime = datetime;

% -------------------------------------------------------------------------
% Call the function
outDirectory = getNldasVsm(qNames, qLat, qLon, qStart, qEnd, outDir);

% -------------------------------------------------------------------------
% Report where the data are held and how long the script took to run
disp('Finished! Site data can be found here:')
disp(outDirectory)
disp('Start and finish times were:')
disp(startTime)
disp(datetime)