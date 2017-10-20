function image_buffer(imageFile,saveDir,FrameNo)% generate matrix of 1's and 0's from buffer
%     imageFile='D:\gpm\2DSimage\2DS.20151201_1.V.cdf';
%     FrameNo=1000;
    buf=ncread(imageFile,'data',[1 1 FrameNo],[8 1700 1])';
    hh=double(ncread(imageFile,'hour',FrameNo,1));
    mm=double(ncread(imageFile,'minute',FrameNo,1));
    ss=double(ncread(imageFile,'second',FrameNo,1));
    hhmmss=hh*10000+mm*100+ss;
    %%
    boundaryInd = []; % index of the particle boundary (series of 8 consecutive '43690' values)
    boundary = [43690 43690 43690 43690 43690 43690 43690 43690];  
    boundaryTime = 0;
    j = 1;
    i=1;
    while buf(j,1) ~= -1 && j+1 <= 1700
        if isequal(buf(j,:),boundary) && isequal(buf(j+1,1),boundaryTime) % particle boundary
            boundaryInd(i)=j;
            i=i+1;
        end
        j = j + 1;
    end
%%
    boundaryData = repmat([2 2 1 1],1, 32); % alternate 1's and 2's for boundary slice (white & cyan pixels)
    buf(buf==-1) = 0; % change invalid values to 0 (unshadowed segment)
    buf = 65535 - buf; % 0: shadowed; 1: unshadowed

    % convert decimal to binary (8 image words for each slice)
    imageData = ones(1700,128); % set up image buffer (1's mean unshadowed pixels)

    for x=1:1700
        tempBuf = [dec2bin(buf(x,1),16); dec2bin(buf(x,2),16); dec2bin(buf(x,3),16);
                            dec2bin(buf(x,4),16); dec2bin(buf(x,5),16); dec2bin(buf(x,6),16);
                            dec2bin(buf(x,7),16); dec2bin(buf(x,8),16)];
        sliceBuf = [];
        for y =0:127
            temp=str2num(tempBuf(floor(y/16)+1,mod(y,16)+1));
            sliceBuf=[sliceBuf;temp];
%             append(sliceBuf,temp);
        end
%         sliceBuf = np.asarray(sliceBuf, dtype='int');
        imageData(x,:) = sliceBuf;
    end
    for i=1:length(boundaryInd)
        imageData(boundaryInd(i),:) = boundaryData;
    end
    if boundaryInd(end)+1 <= 128
        imageData(boundaryInd+1,:) = 1.;
    else
        imageData(boundaryInd(1:end-1)+1,:) = 1.;
    end
    
    particle=cell(length(boundaryInd),1);
    particle{1}=imageData(1:boundaryInd(1)-1,:);
    for i=1:length(boundaryInd)-1
        particle{i+1}=imageData(boundaryInd(i)+2:boundaryInd(i+1)-1,:);
    end
    title=cell(length(boundaryInd),1);
%     saveDir='D:\gpm\HabitClassification\case\20151202_0039-0054\Alexis classification scheme\2DS\inputDir\';
    for i=1:length(boundaryInd)
        title{i}=sprintf('%s%d_%d_%d.png',saveDir,hhmmss,FrameNo,i);
        particle_mask=particle{i}';
        if ~isempty(particle_mask)
            % save(title{i},'particle_mask'); <- to save MatFile
            imwrite(particle_mask,title{i},'png','BitDepth',8);
        end
    end
%     return imageData;
end
    