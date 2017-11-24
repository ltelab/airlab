function [x,y,amax,area_no_holes,name,cout,centroid,adiff,imsize,edgepx,xtremepx] = get_perimeter(filePath, outputDir, plots, debug, enableFilters, BW_mask_in)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% (c) Hannakaisa Lindqvist, Hanne Hakkarainen & Jussi Tiira
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Tätä aliohjelmaa kutsuu classify.m ohjelma. Aliohjelma irrottaa reunakäyrät (x,y) kuvista.
%

%% Config
threshcoef = 1; % jääkiteille sopiva on 1, omenoille 0.8 (riippuu kuvasta; jos on tikku -> 1)
threshbracket = 0.20; % coef+-bracket, 0 for no threshold bracketing
kaanteinen = true; % false (tausta tummempi kuin siluetti: omenoille) tai true (tausta vaaleampi kuin siluetti: jääkiteille)

plot_overlay = true; % Overlay perimeter curve on output image.
crop_frame = true; % Causes a blinking popup window, but does the job.
%% Reunakäyrän irrottaminen kuvasta.

% Tiedostonimen eri osat.
[~,name,ext] = fileparts(filePath);

% disp([sprintf('\nComputing parameters for ') name '...'])

newFileName = [name '_out' ext];
fullFilePath_output = fullfile(outputDir,newFileName);

if ~exist('BW_mask_in','var')
 
    % Ladataan kuva, muutetaan 0-255 kokonaislukudata 0-1 liukuluvuksi
    % (double) ja varmistutaan, että kuva on harmaasävykuva. 
    info = imfinfo(filePath);
    switch info.ColorType
        case 'grayscale'
            im = imread(filePath);
        case 'truecolor'
            im = rgb2gray(im2double(imread(filePath)));
        case 'indexed'
            [X, map] = imread(filePath);
            im = im2double(ind2gray(X, map));
        otherwise
            error('invalid image type')
    end

    % Pienennetään kuva, jotta se on kevyempi käsitellä (lähinnä omenoille). Jos
    % kuitenkin kuva valmiiksi pienempi niin ei muuteta kokoa.
    if info.Width >= info.Height && info.Width > 1024
        im=imresize(im,[NaN 1024]);
    elseif info.Width < info.Height && info.Height > 1024
        im=imresize(im,[1024 NaN]);
    end
    
    % Jos kuvassa objekti tummempi, käännetään kuvan värit siten, että objekti on valkoinen.
    if kaanteinen
        im = imcomplement(im);
    end
    
else
    
    im = BW_mask_in;
    
end

% filtteri reunakäyrän epätasaisuuksien poistamiseksi
im = conv2(double(im),double(ones(5)/9),'full'); % by default : 'valid'

% Luodaan maski käyttäen thresholdia ja sopivaa kerrointa (inputparametri).
mask=im2bw(im,graythresh(im)*threshcoef);

% Test threshold sensitivity
adiff = 0;
if threshbracket ~= 0 && enableFilters
    avgbrightness=mean(mean(im));
    for i=[0 -1]
        props=regionprops(logical(im2bw(im,graythresh(im)*(threshcoef+i*threshbracket*(1-avgbrightness)))),'Area');
        if adiff<0
            k=-1;
        else
            k=1;
        end
        adiff = adiff - k*max([props.Area]);
    end
end

% Määritetään suurin objekti ja käsitellään vain sitä 
% (näin päästään eroon mahdollisista ylimääräisistä, pienemmistä kappaleista).
label = bwlabel(mask);
props=regionprops(label,'FilledArea');
area = [props.FilledArea];

props2=regionprops(label,'Area');
area_noholes = [props2.Area];

% Suurimman objektin pinta-ala.
amax = max(area);
area_no_holes = max(max(area_noholes));
idx = find(area == max(area));
BW = ismember(label, idx);
label2 = logical(BW);

% Määritetään reunat.
B = bwboundaries(BW,'noholes');

%cout=corner(im,B,1.0,165,3,0.35,0,1,1);
cout = corner(im,B,1.0,165,3,0.35,0,1,1);

if (isempty(cout) == 1)
  cout = [0,0;0,0];
elseif (size(cout) == [1,2])
  cout = [0,0;0,0];
elseif (size(cout) == [2,1])
  cout = [0,0;0,0];
end

% Objektin painopiste.
stats=regionprops(label2,'Centroid');
centroid = cat(1, stats.Centroid);

% Poimitaan reunan x- ja y- koordinaatit.
boundaries = B{1,1};
N = length(boundaries);
I = 1:N;
y = boundaries(I,1);
x = boundaries(I,2);

% Image size and boundaries touching image edges
imsize = size(BW);
edgepx = 0; % pixel count on image edges
xtremepx = 0;
for i=[1 2]
    edgepx = edgepx + sum(boundaries(:,i)==0) + sum(boundaries(:,i)==imsize(i));
    xtremepx = xtremepx + sum(boundaries(:,i)==min(boundaries(:,i))) + sum(boundaries(:,i)==max(boundaries(:,i)));
end

% Reunakäyräkuvien piirtäminen yksittäisistä kiteistä output-kansioon.
if (plots || debug)

  if debug
      visible = 'on';
  else
      visible = 'off';
  end
  f = figure('Visible',visible);
  set(0,'defaultaxesfontsize',22,'defaultaxeslinewidth',2,...
        'defaultlinelinewidth',2,'defaultpatchlinewidth',2,...
        'defaulttextfontsize',22);
  if plot_overlay
    imshow(im, 'Border', 'tight');
  end
  hold on
  plot(x,y,'r-')
  plot(cout(:,2),cout(:,1),'gs','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',4)
  plot(centroid(1),centroid(2),'mh','MarkerEdgeColor','k','MarkerFaceColor','m','MarkerSize',8)
  axis equal
  grid off
  if debug||crop_frame
      axesH = gca;
      set(axesH,'units','pixels') % set the axes units to pixels
      set(f,'units','pixels') % set the figure units to pixels
      axpos = get(axesH, 'position');
      figpos = get(f, 'position');
      set(f, 'position', [figpos(1) figpos(2) axpos(3) axpos(4)])
      set(axesH,'units','normalized','position',[0 0 1 1])
      F = getframe(f);
      imwrite(F.cdata,fullFilePath_output)
  else
      saveas(gcf,fullFilePath_output)
  end
  if debug
      pause on
      disp 'DEBUG MODE: Paused. Press any key to continue.'
      pause
  end
  close(gcf)
end % plottaukset
end

