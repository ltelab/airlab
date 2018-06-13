clearvars;
tic;
% CP : I stopped the processing at 1315
% imageFile='D:\gpm\2DSimage\2DS.20151201_1.V.cdf';
imageFile='/ltedata/MASC/OAP/OAP_flight_data/20151112_WF/HVPS/HVPS.20151112.V.cdf';
hour=ncread(imageFile,'hour');
minute=ncread(imageFile,'minute');
second=ncread(imageFile,'second');
% saveDir='D:\gpm\HabitClassification\case\20151202_0039-0054\Alexis classification scheme\2DS\inputDir\';
saveDir='/ltedata/MASC/OAP/OAP_flight_data/20151112_WF/HVPS/19h19-19h31/';
%saveDir='/data/mcfarq/a/saisai/GPM/HabitClassification/case/20151201_230300-230359/2DS/inputDir/';
%%
startTime='191900';%hh:mm:ss
endTime='193100';
process_all = false;

if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

start_hh=str2double(startTime(1:2));
start_mm=str2double(startTime(3:4));
start_ss=str2double(startTime(5:6));
end_hh=str2double(endTime(1:2));
end_mm=str2double(endTime(3:4));
end_ss=str2double(endTime(5:6));

startFrame=min(find(hour==start_hh & minute==start_mm & second==start_ss));
endFrame=max(find(hour==end_hh & minute==end_mm & second==end_ss));

if process_all
    startFrame = 1;
    endFrame = numel(hour);
end

parfor i=startFrame:endFrame % <--- I stopped at 39'500 last time... should be enough :)
    image_buffer(imageFile,saveDir,i);
    fprintf('%u / %u \n',i,endFrame-startFrame+1)
end

finished = toc;
display(finished);

