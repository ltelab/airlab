function glint_idx = compute_glint_idx(data_bs,data_bin,dist_thresh,illustration,illu_path)

    %filename = 'part134__CPI__20151112_192003_033_.mat';
    %load(fullfile(datadir,filename));
    if ~exist('dist_thresh','var')
        dist_thresh = 3; % pixels
    end
    
    if ~exist('illustration','var')
        illustration = false;
    end
    
    M = data_bs;
    M = (M(:,:,1) + M(:,:,2) + M(:,:,3))./3; %M = 1-M;
    M = double(M);
    Mbin = imfill(data_bin,'holes');
    M(Mbin == 0) = NaN;
    
    % calculation of D, a matrix of distances from the particle centroid
    [h,w] = size(M);
    c = [h/2, w/2];
    
    % compute the "intensity weighted" centroid in the vincinity of the geometrical centroid
    if h>=30 && w>=30
        vinc = 2;
    elseif h>=25 && w>=25
        vinc = 1;
    else 
        vinc = 0;
    end
    
    if vinc > 0
        c = ceil(c);
        subM = M(c(1)-vinc:c(1)+vinc, c(2)-vinc:c(2)+vinc);
        [idx1, idx2] = find(subM == max(subM(:)));
        if ~isempty(idx1) && ~isempty(idx2)
            c(1) = c(1) + idx1(1) - (vinc+1);
            c(2) = c(2) + idx2(1) - (vinc+1);
        end
    end
    y = 1:h;
    x = 1:w;
    D = sqrt((y.' - c(1)) .^ 2 + (x - c(2)) .^ 2); 
    D(Mbin == 0) = NaN;
    
    % translate into vectors
    dist = D(:);
    intens = M(:);
    %dist(intens == 0 ) = NaN;
    %intens(intens == 0) = NaN;

    if 0.25 * max(dist) > dist_thresh
    
        % compute intensity at the center (where the bright spot is supposed to be)
        intens_in = quantile(intens(dist < dist_thresh),0.5);

        % compute average intensity within the "black" part of the particle
        D_min = 0.25 * max(dist);
        D_max = 0.75 * max(dist);
        intens_out = quantile(intens(dist >= D_min & dist <= D_max),0.5);
        
    else
        
        D_min = dist_thresh;
        D_max = 0.75 * max(dist);
        intens_in = quantile(intens(dist < dist_thresh),0.5);
        intens_out = quantile(intens(dist >= dist_thresh & dist <= D_max),0.5);
      
    end
    
    intens_out_std = nanstd(intens(dist >= D_min & dist <= D_max));
    

    % compute the glint index in [-1;1]
    glint_idx = (intens_in - intens_out)/min(intens_in,intens_out);
    
    if glint_idx > 0 && intens_out_std > 0.04
        glint_idx = 0;
    end
    
    if illustration
        
        c1 = [65,171,93]./255;
        c2 = [227,74,51]./255;
        c3 = [115,115,115]./255;

        fig1 = figure('Outerposition',[100 1000 1000 400]); hold on; box on;
        subplot(121); hold on; title(sprintf('%s=%2.2f','\gamma_{BR}',glint_idx));% std=%2.2f',glint_idx,intens_out_std));
        hold on;
        imshow(data_bs);
        %set(gca,'ydir','normal');
        
        angle = linspace(0,2*pi,100);
        x1 = c(2) + dist_thresh .* sin(angle);
        y1 = c(1) + dist_thresh .* cos(angle);
        x2 = c(2) + D_min .* sin(angle);
        y2 = c(1) + D_min .* cos(angle);
        x3 = c(2) + D_max .* sin(angle);
        y3 = c(1) + D_max .* cos(angle);
        p1 = plot(c(2),c(1),'rx','linewidth',2); l1 = 'C^*';
        plot(x1,y1,'k-','color',c1,'linewidth',2);
        plot(x2,y2,'k-','color',c2,'linewidth',2);
        plot(x3,y3,'k-','color',c2,'linewidth',2);
        legend(p1,l1); set(gca,'Fontsize',12);
        
        
        subplot(122); hold on; grid on; box on;
        plot(dist,intens,'ko','markerfacecolor',c3);
        plot([dist_thresh dist_thresh],[0 1],'r--','color',c1,'linewidth',1.5);
        plot([D_min D_min],[0 1],'b--','color',c2,'linewidth',1.5);
        plot([D_max D_max],[0 1],'b--','color',c2,'linewidth',1.5);
        plot([0 dist_thresh],[intens_in intens_in],'r-','color',c1,'linewidth',2.5);
        plot([D_min D_max],[intens_out intens_out],'b-','color',c2,'linewidth',2.5);
        xlabel('Distance to centroid C^* [pixels]');
        ylabel('Intensity [0-1]');
        set(gca,'Fontsize',12);
        disp(illu_path);
        export_fig(illu_path,'-r400');
        
        
    end
    
end