% small script to check sub-classification of AG to AG-BR and AG-O for the time period 00h45-50 for F2 (needed for the paper)
clearvars; close all;

%% copy AGG only to a different folder
if 0

    datadir = '/ltedata/MASC/OAP/OAP_flight_data/20151201_WF/CPI/00h45-00h50_proc';
    targetdir = '/ltedata/MASC/OAP/OAP_flight_data/20151201_WF/CPI/00h45-00h50_only_AGG';

    [X, Xlab, Xname, Xt, Xfullprob] = load_processed_2DS_labels(datadir,[],[],true);
    idx = find(X(:,1) == 1); % select AGG

    for i=1:numel(idx)

        source = fullfile(datadir,Xname{idx(i)});
        target = fullfile(targetdir,Xname{idx(i)});
        copyfile(source,target);

    end

end
    
% labels = X(:,1);
% Dmax = X(:,3);
% area = X(:,2);
% MatP = Xfullprob;
        
%% compute ratio of AG-BR and AG-O for the event
datadir = '/ltedata/MASC/OAP/OAP_flight_data/20151201_WF/CPI/00h45-00h50_only_AGG';
[X, Xlab, Xname, Xt, Xfullprob] = load_processed_2DS_labels(datadir,[],[],true);