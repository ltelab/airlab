tic;
% imageFile='D:\gpm\2DSimage\2DS.20151201_1.V.cdf';
imageFile='D:\gpm\HVPSimage\HVPS.20151201_1.V.cdf';
hour=ncread(imageFile,'hour');
minute=ncread(imageFile,'minute');
second=ncread(imageFile,'second');
% saveDir='D:\gpm\HabitClassification\case\20151202_0039-0054\Alexis classification scheme\2DS\inputDir\';
saveDir='D:\gpm\HabitClassification\case\20151202_0039-0054\Alexis classification scheme\HVPS\inputDir\';
%saveDir='/data/mcfarq/a/saisai/GPM/HabitClassification/case/20151201_230300-230359/2DS/inputDir/';
%%
startTime='004600';%hh:mm:ss
endTime='004659';
start_hh=str2double(startTime(1:2));
start_mm=str2double(startTime(3:4));
start_ss=str2double(startTime(5:6));
end_hh=str2double(endTime(1:2));
end_mm=str2double(endTime(3:4));
end_ss=str2double(endTime(5:6));

startFrame=min(find(hour==start_hh & minute==start_mm & second==start_ss));
endFrame=max(find(hour==end_hh & minute==end_mm & second==end_ss));
for i=startFrame:endFrame
    image_buffer(imageFile,saveDir,i);
    i
end

toc;