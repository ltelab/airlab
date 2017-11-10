% some folders :
%'/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/processed_05-Feb-2016/DATA/GOOD';
%'/media/praz/Masc-Data/APRES3_2015/PROCESSED_20160119/DATA/GOOD';


function [X,Xlab,Xname,Xt] = load_processed_2DS_data(dirname,t_str_start,t_str_stop,feat_vec)

    fprintf('Load processed data...');

    if nargin > 1
    
        tmin = datenum(t_str_start,'yyyymmddHHMMSS');
        tmax = datenum(t_str_stop,'yyyymmddHHMMSS');
        
    else
        
        tmin = datenum('20000101000000','yyyymmddHHMMSS');
        tmax = datenum('20300101000000','yyyymmddHHMMSS');
        
    end

    file_list = dir(fullfile(dirname,'*.mat'));
    file_only_list = {file_list.name};
    file_list = fullfile(dirname,file_only_list);

    % loop over the snowflakes and save each feature in a vector of the main
    % structure data
    n_flakes = 0;
    X = [];
    Xlab = {};
    Xname = {};
    Xt = [];
    
    i = 1;
    for k=1:length(file_list)

        load(file_list{k});
        % load only pictures within the time interval and large enough to process
        if (~isfield(roi,'tnum') || (roi.tnum >= tmin && roi.tnum <= tmax)) %&& ~roi.is2small

            j = 1;
            n_flakes = n_flakes + 1;
            
            Xname{i} = file_only_list{i};
            %disp(Xname{i});
            if isfield(roi,'tnum')
                Xt(i) = roi.tnum;
            else 
                Xt = NaN;
            end

            % dimension related features 
            X(i,j) = get_struct_field(roi,'area'); Xlab{j} = 'area'; j=j+1; %1
            X(i,j) = get_struct_field(roi,'Dmax'); Xlab{j} = 'Dmax'; j=j+1; %2
            X(i,j) = get_struct_field(roi,'width'); Xlab{j} = 'width'; j=j+1; %3
            X(i,j) = get_struct_field(roi,'height'); Xlab{j} = 'height'; j=j+1; %4
            X(i,j) = get_struct_field(roi,'perim'); Xlab{j} = 'perim'; j=j+1; %5
            X(i,j) = get_struct_field(roi,'eq_radius'); Xlab{j} = 'eq radius'; j=j+1; %6
            X(i,j) = get_struct_field(roi,'area_porous'); Xlab{j} = 'porous area'; j=j+1; %7
            X(i,j) = (get_struct_field(roi,'area') - get_struct_field(roi,'area_porous'))/get_struct_field(roi,'area'); Xlab{j} = 'porous ratio'; j=j+1; %8

            % ellipse related features
            X(i,j) = get_struct_field(roi,'E','a'); Xlab{j} = 'ellipse fit A dim'; j=j+1; %9
            X(i,j) = get_struct_field(roi,'E','b'); Xlab{j} = 'ellipse fit B dim'; j=j+1; %10
            X(i,j) = get_struct_field(roi,'E','theta'); Xlab{j} = 'orientation'; j=j+1; %11
            X(i,j) = get_struct_field(roi,'E_in','a'); Xlab{j} = 'ellipse in A dim'; j=j+1; %12
            X(i,j) = get_struct_field(roi,'E_in','b'); Xlab{j} = 'ellipse in B dim'; j=j+1; %13
            X(i,j) = get_struct_field(roi,'E_out','a'); Xlab{j} = 'ellipse out A dim'; j=j+1; %14
            X(i,j) = get_struct_field(roi,'E_out','b'); Xlab{j} = 'ellipse out B dim'; j=j+1; %15
            X(i,j) = pi*get_struct_field(roi,'E_in','a')*get_struct_field(roi,'E_in','b'); Xlab{j} = 'ellipse in area'; j=j+1; %16
            X(i,j) = pi*get_struct_field(roi,'E_out','a')*get_struct_field(roi,'E_out','b'); Xlab{j} = 'ellipse out area'; j=j+1; %17
            X(i,j) = pi*get_struct_field(roi,'E','a')*get_struct_field(roi,'E','b'); Xlab{j} = 'ellipse fit area'; j=j+1; %18
            X(i,j) = get_struct_field(roi,'E','a')/get_struct_field(roi,'E_out','a'); Xlab{j} = 'ratio ellipse A fit/out'; j=j+1; %19
            X(i,j) = get_struct_field(roi,'E','b')/get_struct_field(roi,'E_out','b'); Xlab{j} = 'ratio ellipse B fit/out'; j=j+1; %20 
            X(i,j) = get_struct_field(roi,'E_in','a')/get_struct_field(roi,'E_out','a'); Xlab{j} = 'ratio ellipse A in/out'; j=j+1; %21
            X(i,j) = get_struct_field(roi,'E_in','b')/get_struct_field(roi,'E_out','b'); Xlab{j} = 'ratio ellipse B in/out'; j=j+1; %22
            X(i,j) = get_struct_field(roi,'E_in','a')/get_struct_field(roi,'E','a'); Xlab{j} = 'ratio ellipse A in/fit'; j=j+1; %23
            X(i,j) = get_struct_field(roi,'E_in','b')/get_struct_field(roi,'E','b'); Xlab{j} = 'ratio ellipse B in/fit'; j=j+1; %24
            X(i,j) = (get_struct_field(roi,'E_in','a')*get_struct_field(roi,'E_in','b')) / (get_struct_field(roi,'E_out','a')*get_struct_field(roi,'E_out','b')); Xlab{j} = 'ratio ellipse area in/out'; j=j+1; %25
            X(i,j) = (get_struct_field(roi,'E_in','a')*get_struct_field(roi,'E_in','b')) / (get_struct_field(roi,'E','a')*get_struct_field(roi,'E','b')) ; Xlab{j} = 'ratio ellipse area in/fit'; j=j+1; %26
            X(i,j) = (get_struct_field(roi,'E','a')*get_struct_field(roi,'E','b')) / (get_struct_field(roi,'E_out','a')*get_struct_field(roi,'E_out','b')); Xlab{j} = 'ratio ellipse area fit/out'; j=j+1; %27

            % Skeleton
            X(i,j) = get_struct_field(roi,'skel','p_ratio'); Xlab{j} = 'skel/perim ratio'; j=j+1; %28
            X(i,j) = get_struct_field(roi,'skel','A_ratio'); Xlab{j} = 'skel/area ratio'; j=j+1; %29
            if ~isnan(get_struct_field(roi,'skel','N_ends'))
                X(i,j) = roi.skel.N_ends; Xlab{j} = 'skel N ends'; j=j+1; %30  
            else
                X(i,j) = 0; Xlab{j} = 'skel N ends'; j=j+1; %30  
            end
            X(i,j) = get_struct_field(roi,'skel','N_junctions'); Xlab{j} = 'skel N junctions'; j=j+1; %31

            % Fractal
            X(i,j) = get_struct_field(roi,'F'); Xlab{j} = 'boxcounting fractal dim'; j=j+1; %32
            X(i,j) = get_struct_field(roi,'F_jac'); Xlab{j} = 'theoretical fractal dim'; j=j+1; %33

            % Hull
            X(i,j) = get_struct_field(roi,'hull','solidity'); Xlab{j} = 'solidity'; j=j+1; %34
            X(i,j) = get_struct_field(roi,'hull','convexity'); Xlab{j} = 'convexity'; j=j+1; %35
            X(i,j) = length(get_struct_field(roi,'hull','xh')); Xlab{j} = 'Hull nb angles'; j=j+1; %36
               
            % Others
            X(i,j) = get_struct_field(roi,'roundness'); Xlab{j} = 'roundness'; j=j+1; %37
            X(i,j) = get_struct_field(roi,'compactness'); Xlab{j} = 'compactness'; j=j+1; %38
            if get_struct_field(roi,'E','a') >= get_struct_field(roi,'E','b')
                X(i,j) = roi.E.b/roi.E.a; Xlab{j} = 'aspect ratio'; j=j+1; %39
                X(i,j) = sqrt(1-roi.E.b/roi.E.a); Xlab{j} = 'eccentricity'; j=j+1; %40
            else
                X(i,j) = 1; Xlab{j} = 'aspect ratio'; j=j+1; %39
                X(i,j) = 0; Xlab{j} = 'eccentricity'; j=j+1; %40
            end
            
            if get_struct_field(roi,'nb_holes') > 0
               X(i,j) = 1; Xlab{j} = 'has holes'; j=j+1; %41
            else
               X(i,j) = 0; Xlab{j} = 'has holes'; j=j+1; %41
            end
            X(i,j) = get_struct_field(roi,'nb_holes'); Xlab{j} = 'N holes'; j=j+1; %42
            
            % New descriptors based on min rectangle box
            X(i,j) = get_struct_field(roi,'Rect','A_ratio'); Xlab{j} = 'rectangularity'; j=j+1; %43
            X(i,j) = get_struct_field(roi,'Rect','p_ratio'); Xlab{j} = 'rect perim ratio'; j=j+1; %44
            X(i,j) = get_struct_field(roi,'Rect','aspect_ratio'); Xlab{j} = 'rect aspect ratio'; j=j+1; %45
            X(i,j) = get_struct_field(roi,'Rect','eccentricity'); Xlab{j} = 'rect eccentricity'; j=j+1; %46
            
            % New descriptor based on circumscribed circle
            X(i,j) = (2*pi*get_struct_field(roi,'C_out','r'))/get_struct_field(roi,'perim'); Xlab{j} = 'circle out perim ratio'; j=j+1; %47  
            X(i,j) = get_struct_field(roi,'complex'); Xlab{j} = 'complex'; j=j+1; %48

            % New descriptors which are probably useless
            tmp_vec = [];
            X(i,j) = get_struct_field(roi,'Sym','P0'); Xlab{j} = 'fft P0'; j=j+1; %49
            X(i,j) = get_struct_field(roi,'Sym','P1'); tmp_vec(end+1) = X(i,j); Xlab{j} = 'fft P1'; j=j+1; %50
            X(i,j) = get_struct_field(roi,'Sym','P2'); tmp_vec(end+1) = X(i,j); Xlab{j} = 'fft P2'; j=j+1; %51
            X(i,j) = get_struct_field(roi,'Sym','P3'); tmp_vec(end+1) = X(i,j); Xlab{j} = 'fft P3'; j=j+1; %52
            X(i,j) = get_struct_field(roi,'Sym','P4'); tmp_vec(end+1) = X(i,j); Xlab{j} = 'fft P4'; j=j+1; %53
            X(i,j) = get_struct_field(roi,'Sym','P5'); tmp_vec(end+1) = X(i,j); Xlab{j} = 'fft P5'; j=j+1; %54
            X(i,j) = get_struct_field(roi,'Sym','P6'); tmp_vec(end+1) = X(i,j); Xlab{j} = 'fft P6'; j=j+1; %55
            X(i,j) = get_struct_field(roi,'Sym','P7'); tmp_vec(end+1) = X(i,j); Xlab{j} = 'fft P7'; j=j+1; %56
            X(i,j) = get_struct_field(roi,'Sym','P8'); tmp_vec(end+1) = X(i,j); Xlab{j} = 'fft P8'; j=j+1; %57
            X(i,j) = get_struct_field(roi,'Sym','P9'); tmp_vec(end+1) = X(i,j); Xlab{j} = 'fft P9'; j=j+1; %58
            X(i,j) = get_struct_field(roi,'Sym','P10'); tmp_vec(end+1) = X(i,j); Xlab{j} = 'fft P10'; j=j+1; %59
            X(i,j) = get_struct_field(roi,'Sym','mean'); Xlab{j} = 'Sym mean'; j=j+1; %60
            X(i,j) = get_struct_field(roi,'Sym','std'); Xlab{j} = 'Sym std'; j=j+1; %61
            X(i,j) = get_struct_field(roi,'Sym','std')/get_struct_field(roi,'Sym','mean'); Xlab{j} = 'Sym std/mean'; j=j+1; %962
            X(i,j) = get_struct_field(roi,'Sym','P6')/nanmax(tmp_vec); Xlab{j} = 'fft P6/Pmax'; j=j+1; %63
            [~,idx_max] = nanmax(tmp_vec); 
            X(i,j) = idx_max-1; Xlab{j} = 'fft P# max'; j=j+1; %64 <--- check why is there idx_max-1 here, why -1 ?
            
            % new descriptor based on truncated ellipse fitted (useful to detect truncated Spheres)
            X(i,j) = get_struct_field(roi,'E2','fit_OA'); Xlab{j} = 'truncated ellipse fitted OA'; j=j+1; %65
            X(i,j) = get_struct_field(roi,'touch_edge'); Xlab{j} = 'is touching image edge'; j=j+1; %66
            X(i,j) = get_struct_field(roi,'touch_edge_ratio'); Xlab{j} = 'ratio of perimeter touching the edge'; j=j+1; %67
            X(i,j) = get_struct_field(roi,'n_true'); Xlab{j} = 'total number of 1 pixels'; j=j+1; %68
            X(i,j) = get_struct_field(roi,'ntrue_ratio'); Xlab{j} = 'ratio 1 pixels within mask / total 1 pixels'; j=j+1; % 69
            
            % new descriptors based on IC-PCA
            X(i,j) = get_struct_field(roi,'icpca','b','a_skaalattu_poly'); Xlab{j} = 'a_skaalattu_poly'; j=j+1; %70
            X(i,j) = get_struct_field(roi,'icpca','b','roundness2'); Xlab{j} = 'roundness2'; j=j+1; %71
            X(i,j) = get_struct_field(roi,'icpca','b','roundness2max'); Xlab{j} = 'roundness2max'; j=j+1; %72
            X(i,j) = get_struct_field(roi,'icpca','b','rectangularity'); Xlab{j} = 'icpca-rectangularity'; j=j+1; %73
            X(i,j) = get_struct_field(roi,'icpca','b','asprat2'); Xlab{j} = 'asprat2'; j=j+1; %74
            X(i,j) = get_struct_field(roi,'icpca','b','cornerarea'); Xlab{j} = 'cornerarea'; j=j+1; %75
            X(i,j) = get_struct_field(roi,'icpca','b','holearea'); Xlab{j} = 'holearea'; j=j+1; %76
            X(i,j) = get_struct_field(roi,'icpca','b','number_of_corners'); Xlab{j} = 'num of corners'; j=j+1; %77
            X(i,j) = get_struct_field(roi,'icpca','b','mean_lengths'); Xlab{j} = 'mean lengths'; j=j+1; %78
            X(i,j) = get_struct_field(roi,'icpca','b','std_lengths'); Xlab{j} = 'std lengths'; j=j+1; %79
            jana_ka = get_struct_field(roi,'icpca','b','jana_ka');
            if ~isnan(jana_ka)
                X(i,j) = jana_ka(1); Xlab{j} = 'jana_ka1'; j=j+1; %80
                X(i,j) = jana_ka(2); Xlab{j} = 'jana_ka2'; j=j+1; %81
                X(i,j) = jana_ka(3); Xlab{j} = 'jana_ka3'; j=j+1; %82
            else
                X(i,j) = NaN; Xlab{j} = 'jana_ka1'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_ka2'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_ka3'; j=j+1;
            end
            jana_kov = get_struct_field(roi,'icpca','b','jana_kov');
            if ~isnan(jana_kov)
                X(i,j) = jana_kov(1); Xlab{j} = 'jana_kov1'; j=j+1; %83
                X(i,j) = jana_kov(2); Xlab{j} = 'jana_kov2'; j=j+1; %84
                X(i,j) = jana_kov(3); Xlab{j} = 'jana_kov3'; j=j+1; %85
                X(i,j) = jana_kov(4); Xlab{j} = 'jana_kov4'; j=j+1; %86
                X(i,j) = jana_kov(5); Xlab{j} = 'jana_kov5'; j=j+1; %87
                X(i,j) = jana_kov(6); Xlab{j} = 'jana_kov6'; j=j+1; %88
                X(i,j) = jana_kov(7); Xlab{j} = 'jana_kov7'; j=j+1; %89
                X(i,j) = jana_kov(8); Xlab{j} = 'jana_kov8'; j=j+1; %90
                X(i,j) = jana_kov(9); Xlab{j} = 'jana_kov9'; j=j+1; %91
                X(i,j) = jana_kov(10); Xlab{j} = 'jana_kov10'; j=j+1; %92
                X(i,j) = jana_kov(11); Xlab{j} = 'jana_kov11'; j=j+1; %93
                X(i,j) = jana_kov(12); Xlab{j} = 'jana_kov12'; j=j+1; %94
            else
                X(i,j) = NaN; Xlab{j} = 'jana_kov1'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_kov2'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_kov3'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_kov4'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_kov5'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_kov6'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_kov7'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_kov8'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_kov9'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_kov10'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_kov11'; j=j+1;
                X(i,j) = NaN; Xlab{j} = 'jana_kov12'; j=j+1;
            end
            
            X(i,j) = get_struct_field(roi,'D90','Dmax_90'); Xlab{j} = 'Dmax 90'; j=j+1; % 95
            X(i,j) = get_struct_field(roi,'D90','AR'); Xlab{j} = 'D90 AR'; j=j+1; % 96
            X(i,j) = get_struct_field(roi,'D90','angle'); Xlab{j} = 'D90 angle'; j=j+1; % 97
            X(i,j) = get_struct_field(roi,'frame_fraction'); Xlab{j} = 'perim touching frame/frame length'; j=j+1; %98
            
            
            % textural descriptors
            if strcmp(roi.probe,'CPI')
                X(i,j) = get_struct_field(roi,'mean_intens'); Xlab{j} = 'mean intens'; j=j+1; %99
                X(i,j) = get_struct_field(roi,'max_intens'); Xlab{j} = 'max intens'; j=j+1; %100 
                X(i,j) = get_struct_field(roi,'range_intens'); Xlab{j} = 'range intens'; j=j+1; %101
                X(i,j) = get_struct_field(roi,'focus'); Xlab{j} = 'focus'; j=j+1; %102
                X(i,j) = get_struct_field(roi,'std'); Xlab{j} = 'global std'; j=j+1; %103
                X(i,j) = get_struct_field(roi,'local_std'); Xlab{j} = 'avg local std 3x3'; j=j+1; %104
                X(i,j) = get_struct_field(roi,'lap'); Xlab{j} = 'lap energy'; j=j+1; %105
                X(i,j) = get_struct_field(roi,'wavs'); Xlab{j} = 'wavelet desc'; j=j+1; % 106
                X(i,j) = get_struct_field(roi,'hist_entropy'); Xlab{j} = 'hist_entropy'; j=j+1; %107
                X(i,j) = get_struct_field(roi,'H','Contrast'); Xlab{j} = 'Haralick contrast'; j=j+1; %108
                X(i,j) = get_struct_field(roi,'H','Correlation'); Xlab{j} = 'Haralick correlation'; j=j+1; %109
                X(i,j) = get_struct_field(roi,'H','Energy'); Xlab{j} = 'Haralick energy'; j=j+1; %110
                X(i,j) = get_struct_field(roi,'H','Homogeneity'); Xlab{j} = 'Haralick homogeneity'; j=j+1; %111
            end
    
            i = i+1;

        end

    end
    
    % 4th argument is list of features wanted
    if nargin==4 && ~isempty(X)

        X = X(:,feat_vec);
        Xlab = Xlab(feat_vec);

    end
    
    
    fprintf('   Done! %u images found and processed. %u descriptors computed. \n',n_flakes,size(X,2));
    Xlab = Xlab';
    Xname = Xname';
    Xt = Xt';
    
end