%% small script to identify 2DS/HVPS/CPI images and isolate a selected time interval dataset for validation purposes

%CPI
datadir = '/ltedata/MASC/OAP/OAP_flight_data/20151112_WF/CPI/all_proc';
targetdir = '/ltedata/MASC/OAP/OAP_flight_data/20151112_WF/CPI/19h20-19h30_proc';

if ~exist(targetdir,'dir')
    mkdir(targetdir);
end


t_str_date = '20151112';
t_str_date2 = '20151112';
t_str_start = '20151112192000';
t_str_stop  = '20151112193000';
t_num_start = datenum(t_str_start,'yyyymmddHHMMSS');
t_num_stop = datenum(t_str_stop,'yyyymmddHHMMSS');


all_img = dir(fullfile(datadir,'**','part*__*'));
all_folder = {all_img.folder}';
all_img = {all_img.name}';

k=0;
ntot = numel(all_img);
for i=1:ntot
    
    fprintf('%u/%u\n',i,ntot); 
    current_img = all_img{i};
    s = strfind(current_img,t_str_date);
    if ~isempty(s)
        cat_string = current_img(s:s+14);
    else
        s = strfind(current_img,t_str_date2);
        if ~isempty(s)
            cat_string = current_img(s:s+14);
        else
            fprintf('WARNING : no date found in image %s \n!',current_img);
        end
    end
    
    tmp_datenum = datenum(cat_string,'yyyymmdd_HHMMSS');
    if tmp_datenum > t_num_start && tmp_datenum <= t_num_stop
        k=k+1;
        fprintf('%u image(s) found yet! \n',k);
        copyfile(fullfile(all_folder{i},all_img{i}),fullfile(targetdir,all_img{i}));
    end
        
    
end
