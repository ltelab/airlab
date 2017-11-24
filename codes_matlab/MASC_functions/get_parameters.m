function out_params = get_parameters(filename, inputDir, paramfile, namesfile, plots, aste, x, y, amax, area_no_holes, name, cout, centroid,adiff,imsize,edgepx,xtremepx,enableFilters,masterPath)
% [maxdist_orig,area_ratio,area_ratio_filled,asprat2,number_of_corners,th,flags]<---initialoutput
% GET_PARAMETERS
% A subroutine for calculating particle parameters
%
% (c) Hannakaisa Lindqvist, Hanne Hakkarainen & Jussi Tiira
% This is a subroutine called by classify.m
% The subroutine computes needed parameters from x,y perimeter curve.

%% CONFIG

% Choose the parameters to be calculated.
% K�ytt�j�n t�ytyy valita JOKO objektin_kaanto TAI pystyasento.
% My�s JOKO maxdist_skaala TAI convex_skaala.
objektin_kaanto = true; % Rotate object so that the maximum dimension is along x-axis
janat = true;           % Kahden pisteen v�lisen janan pituuden keskiarvo ja autokovarianssi. jana_ka, jana_kov
normaalit = false; 		% Ulkonormaalien v�lisen kulman muutoksen keskiarvo ja autokovarianssi. normaali_ka, normaali_kov
tahkopisteet = false;	% Reunak�yr�n approksimoiminen tahkoilla. Sopii kulmikkaille kappaleille.
plottaukset = false; 	% Yksitt�isten reunak�yr�kuvaajien piirt�minen output-kansioon.

% Misc
px_per_unit = 43.0; % pixels per unit distance (or angle). For CPI images 43px = 1 micron.
reunat_tiedostoon = false;	% Reunakäyrien kirjoitus tiedostoon 'reunapisteet.out'.
verbose = false;

% Filter config
skip_filtered = false; % if set to false, parameters will be calculated also for the filtered particles
skip_small = false;
filter_dirs = true; % move/copy filtered images to their own subfolders
    move = false; % Move images instead of copying. See option filter_dirs
small_dirs = true; % move/copy small crystals
    move_small = false;
small_limit_um = 100; % Size limit for small crystals in human readable units (um for ice crystals).

% Filters
filter_small = false;
filter_cropped = true;
filter_blurry = true;
filter_rect = true;

%% autoconfig
if strcmp(masterPath,pwd)
    masterOut = false;
else
    masterOut = true;
end

%% Reunakäyrän pituuden normittaminen 2*pi:hin ja kuutiosplini-interpoloinnit.

[~,kharea1] = convhull(x,y);
area2 = polyarea(x,y);
a_skaalattu_poly = area2./kharea1;

number_of_corners = size(cout,1);
area_inside_corners = polyarea(cout(:,1),cout(:,2));
cornerarea = area_inside_corners ./ area2;
holearea = 1- (area_no_holes./amax);
holearea = holearea(1);

N = length(x);
reunan_pituus = 1:N;
reunan_pituus(1) = 0.0; 
for i = 2:N
  reunan_pituus(i) = reunan_pituus(i-1) + sqrt((x(i)-x(i-1)).^2 + (y(i)-y(i-1)).^2);
end

% Parametrissä alkup_reuna on tallessa reunakäyrän alkuperäinen pituus.
alkup_reuna = reunan_pituus(N) + sqrt((x(N)-x(1)).^2 + (y(N)-y(1)).^2);

% Normitus.
const = (2*pi)/alkup_reuna;  
norm_pituudet = reunan_pituus.*const; 
amax_scaled = amax.*(const.^2);
area2_scaled = area2.*(const.^2);
area_no_holes_scaled = area_no_holes.*(const.^2);
small_limit = small_limit_um*0.01*px_per_unit*const; % Analyze only large crystals (D_max > small_limit)

% Kuutiosplini-interpoloinnit, joiden tarkoituksena on saada reunakäyrältä tasavälein 
% pisteet x ja y (eli yhtä pitkin välimatkoin). Erikseen x:lle ja y:lle.
% Input-parametri aste ilmoittaa resoluution, eli millä väleillä halutaan pisteitä 
% käyrältä (1.0 = 1 asteen välein = yhteensä 360 pistettä).
T2 = (0:aste*(pi./180.0):(2*pi)-(aste*(pi./180.0)))';
y2 = spline(norm_pituudet,y,T2);
x2 = spline(norm_pituudet,x,T2);
x2 = x2.*const;
y2 = y2.*const;
centroid = centroid.*const;

% Normitus vielä uudestaan.
x2length=length(x2);
reunan_pituus2 = 1:x2length;
reunan_pituus2(1) = 0.0; 
for i = 2:x2length
  reunan_pituus(i) = reunan_pituus(i-1) + sqrt((x2(i)-x2(i-1)).^2 + (y2(i)-y2(i-1)).^2);
end
reunan_pituus(x2length+1) = reunan_pituus(x2length) + sqrt((x2(x2length)-x2(1)).^2 + (y2(length(y2))-y2(1)).^2);
const2 = (2*pi)/reunan_pituus(x2length+1); 
x2 = x2.*const2;
y2 = y2.*const2;
centroid = centroid.*const2;
small_limit = small_limit*const2;
amax_scaled = amax_scaled.*(const2.^2);
area2_scaled = area2_scaled.*(const2.^2);
area_no_holes_scaled = area_no_holes_scaled.*(const2.^2);

%% Reunakäyrän kahden pisteen välisen suurimman mahdollisen etäisyyden etsiminen.

[maxdist,rivi,sarake] = diameter(x2,y2);

maxdist_orig = (maxdist./(const2.*const)).*(100.0./px_per_unit);
max_diameter = maxdist_orig;
area_ratio = area_no_holes_scaled./(pi.*(0.5.*maxdist).^2);

% Pisteet, joilla maksimiarvo saadaan.
xp = [x2(rivi), x2(sarake)]';
yp = [y2(rivi), y2(sarake)]';

%% Filters
flags=''; % initialize
if enableFilters
    flags=''; % initialize
    filter_code = '';
    isSmall = maxdist < small_limit;
    
    if (filter_small && isSmall)
        filter_code = ['size_below_' num2str(small_limit_um) '_microns'];
        if small_dirs
            filter2dirs(filename, inputDir, filter_code, false, move_small);
        end
    end
    if (filter_cropped && outOfBounds(imsize,edgepx))
        filter_code = 'bad_crop';
        flags = filter_actions(filter_dirs, filename, inputDir, filter_code, plots, move, flags);
    end
    if (filter_blurry && blurryBoundaries(adiff,amax))
        filter_code = 'blurry';
        flags = filter_actions(filter_dirs, filename, inputDir, filter_code, plots, move, flags);
    end
    if (filter_rect && rectBounds(imsize,xtremepx))
        filter_code = 'rectangular';
        flags = filter_actions(filter_dirs, filename, inputDir, filter_code, plots, move, flags);
    end
    
    if ~isempty(filter_code)
        disp(['Particle filtered: ' filter_code])
        if (skip_filtered||(skip_small&&isSmall))
            area_ratio_filled=NaN;
            asprat2 = NaN;
            return
        end
    end
end

%% Object rotation
% siten, että kahden pisteen välinen pisin mahdollinen etäisyys 
% laitetaan uudeksi vaaka-akseliksi (sopii esim. jääkiteille).

if objektin_kaanto

  apu = 0;
  xx=zeros(1,x2length+1);
  yy=xx;
% Järjestetään x,y parit siten, että pisteet alkavat tulevalta x-akselilta ja kiertävät pos.suuntaan.
  for j = rivi:x2length 
    apu = apu + 1; 
    xx(apu) = x2(j);
    yy(apu) = y2(j);
  end
  for j = 1:rivi
    apu = apu + 1;
    xx(apu) = x2(j);
    yy(apu) = y2(j);
  end

% Move and rotate the coordinate system. Check where the new x axis is
% pointing.
  if (yp(2)-yp(1) >= 0) % 1st and 2nd quarter
    th = kulma([(xp(2)-xp(1)), (yp(2)-yp(1))]',[1, 0]');  
  elseif (yp(2)-yp(1) < 0) % 3rd and 4th quarter
    th = 2*pi-(kulma([(xp(2)-xp(1)), (yp(2)-yp(1))]',[1, 0]')); 
  end
  xk=zeros(1,x2length); % Initialize xk and yk.
  yk=xk;
  for j = 1:x2length
    xx(j) = xx(j) - (xp(1)+xp(2))./2; % move
    yy(j) = yy(j) - (yp(1)+yp(2))./2;
    xy = [xx(j),yy(j)]*[cos(th), -sin(th); sin(th), cos(th)]; % rotate
    xk(j) = xy(1); 
    yk(j) = xy(2);
  end
  centroid(1) = centroid(1) - (xp(1)+xp(2))./2; % move
  centroid(2) = centroid(2) - (yp(1)+yp(2))./2;
  centroidxy = [centroid(1),centroid(2)]*[cos(th), -sin(th); sin(th), cos(th)]; % rotate
  centroid(1) = centroidxy(1); 
  centroid(2) = centroidxy(2);
  for j = 1:x2length
    xk(j) = xk(j) - centroid(1); % siirretään
    yk(j) = yk(j) - centroid(2);
  end

end % objektin_kaanto

roundprec = 50; % rounding distance. 50 seems to work, might not be optimal.
width = 0;
for k = unique(round(xk*roundprec))
    ind = round(xk*roundprec)==k;
    len = max(yk(ind))-min(yk(ind));
    if len>width, width = len; end;
end

roundness2 = area2_scaled./(pi.*((width./2).^2));
roundness2max = area2_scaled./(pi.*((maxdist./2).^2));
area_ratio_filled = amax_scaled./(pi.*((maxdist./2).^2));
asprat2 = maxdist./width;
rectangularity = amax_scaled./(width.*maxdist);

if verbose
    disp(['Max diameter: ' num2str(max_diameter)])
    disp(['Area ratio: ' num2str(area_ratio)])
    disp(['Area ratio with holes filled: ' num2str(area_ratio_filled)])
end

%% sivutahkojen väliset piste- ja ristitulot 
% (tai oikeastaan cos ja sin, koska nämä on jaettu vektorien pituuksien tulolla):

[dots,crosses,lengths] = vec_ops(cout);
crosses = crosses(:,3);
kosinit = zeros(1,5);
sinit = zeros(1,2);

lengths = lengths./alkup_reuna;

std_lengths = std(lengths);
mean_lengths = mean(lengths);

for j = 1:length(dots)
  if (dots(j) < -0.77), kosinit(1) = kosinit(1) + 1; end
  if (dots(j) > -0.77 && dots(j) < -0.26), kosinit(2) = kosinit(2) + 1; end
  if (dots(j) > -0.26 && dots(j) < 0.26), kosinit(3) = kosinit(3) + 1; end
  if (dots(j) > 0.26 && dots(j) < 0.77), kosinit(4) = kosinit(4) + 1; end
  if (dots(j) > 0.77), kosinit(5) = kosinit(5) + 1; end
  if (crosses(j) < 0.0), sinit(1) = sinit(1) + 1; end
  if (crosses(j) > 0.0), sinit(2) = sinit(2) + 1; end
end


%% Kahden pisteen välisen janan pituuden autokovarianssi ja keskiarvo.
% Tässä oletetaan, että reunakäyrän kulmaresoluutio on 1 aste! 
% Parametreista:
% dgamma = se invariantin kulman arvo, jolta väliltä mitataan janan pituus
% deltaphi = kahden eri janan aloituspisteiden välinen kulmaero 
% dstep = askeleen suuruus, eli kuinka suurin välein kierretään koko 360 asteen kierros

if janat

  dstep = 3; 
  apu = 0;
  apu2 = 0;

% Koordinaatit kaksi täyttä kierrosta, jotta ympäri meneminen on helpompaa.  
  xk2 = [xk(1:(length(xk))) xk(1:(length(xk)))];  
  yk2 = [yk(1:(length(yk))) yk(1:(length(yk)))];
  jana_ka=zeros(1,3);
  d1=zeros(1,length(xk));
  jana_kov=zeros(1,12);
  for dgamma = [30 90 180];
    apu2 = apu2 + 1;
% d-vektorin elementeissä on tieto siitä, kuinka pitkä etäisyys on deltagamman päässä
% olevaan käyrän pisteeseen (tämä siis tulee samalla resoluutiolla kuin xk).
% Huom. jana_kov(1) = varianssi!
    for ii = 1:length(xk)
      d1(ii) = sqrt((xk2(ii)-xk2(ii+dgamma)).^2 + (yk2(ii)-yk2(ii+dgamma)).^2);
    end
    for dphi = [0 30 90 180]; % orig: 0,5,10,30,90,180
      apu = apu + 1;
      [dcov1,dmean1] = autokovis(d1,dphi,dstep);
      jana_kov(apu) = dcov1;
    end
    jana_ka(apu2) = dmean1;
  end

end % janat

%% Pinnan ulkonormaalin suunnan muutoksien keskiarvo ja autokovarianssi.

if normaalit

  for j3 = 1:length(xk)-1
    nvec1 = [(-yk(j3+1)+yk(j3)),(xk(j3+1)-xk(j3)),0];
    nvec(j3,:) = nvec1./sqrt((-yk(j3+1)+yk(j3)).^2+(xk(j3+1)-xk(j3)).^2);
  end

% Matriisissa 'nvec' on ulkonormaalivektorit riveittäin. 
% Kovarianssianalyysin parametrien selitys: ks. "janat". 

  dstep = 3; 
  apu = 0;
  apu2 = 0;
  nvec = [nvec;nvec];
  for dgamma = [5 10 30 90 180];
    apu2 = apu2 + 1;

% Vektorin 'angle' elementeissä on tieto siitä, kuinka suuri kulman muutos on 
% deltagamman päässä olevaan käyrän ulkonormaaliin. 
    for ii = 1:length(xk)
      angle1(ii) = kulma(nvec(ii,:),nvec(ii+dgamma,:));
    end
    for dphi = [0 10 30 90 180];
      apu = apu + 1;
      [angcov1,angmean1] = autokovis(angle1,dphi,dstep);
      normaali_kov(apu) = angcov1;
    end
    normaali_ka(apu2) = angmean1;
  end

end % normaalit

%% Tunnuslukujen kirjoitus tiedostoon.

% the parameters
savetus = [a_skaalattu_poly,roundness2,roundness2max,rectangularity,asprat2,cornerarea,holearea,number_of_corners,a_skaalattu_poly,asprat2,mean_lengths,std_lengths,kosinit,sinit,jana_ka,jana_kov];
out_params.b.a_skaalattu_poly = a_skaalattu_poly;
out_params.b.roundness2 = roundness2;
out_params.b.roundness2max = roundness2max;
out_params.b.rectangularity = rectangularity;
out_params.b.asprat2 = asprat2;
out_params.b.cornerarea = cornerarea;
out_params.b.holearea = holearea;
out_params.b.number_of_corners = number_of_corners;
out_params.b.mean_lengths = mean_lengths;
out_params.b.std_lengths = std_lengths;
out_params.b.kosinit = kosinit;
out_params.b.sinit = sinit;
out_params.b.jana_ka = jana_ka;
out_params.b.jana_kov = jana_kov;

out_params.a.maxdist_orig = maxdist_orig;
out_params.a.area_ratio = area_ratio;
out_params.a.area_ratio_filled = area_ratio_filled;
out_params.a.asprat2 = asprat2;
out_params.a.number_of_corners = number_of_corners;
out_params.a.th = th;

%%%%%% I don't want to save the text file, hence commented part %%%%%
% save (paramfile, 'savetus', '-ASCII','-append')
% 
% fid = fopen(namesfile, 'a');
% fprintf(fid, '%s\n', name);
% fclose(fid);
% 
% if masterOut
%     masterfile = fullfile(masterPath,'class_parameters.out');
%     save (masterfile, 'savetus', '-ASCII','-append');
% 
%     masterfile2 = fullfile(masterPath, 'filenames.out');
%     mid = fopen(masterfile2, 'a+');
%     fprintf(mid, '%s\n', name);
%     fclose(mid);
% end %masterOut

%% Reunakäyrän koordinaattien kirjoitus tiedostoon. Ensin x-koordinaatit, sitten y.  

if reunat_tiedostoon
  reunakoordinaatit = [xk yk];
  filenimi = fullfile(outputDir,'reunapisteet.out');
  save (filenimi, 'reunakoordinaatit', '-ASCII','-append')
end % reunat_tiedostoon

%% Reunakäyräkuvien piirtäminen yksittäisistä kiteistä output-kansioon.

if plottaukset
  figure('Visible','off');
  if tahkopisteet
    plot(xk,yk,'k--',kulmapisteet(:,1),kulmapisteet(:,2),'bo')
  else
    plot(xk,yk,'k--')
  end
  axis equal
  grid on
  saveas(gcf,fullFilePath_output)
  close(gcf)
  clear all
end % plottaukset
end

%% Subfunctions

function [autokovarianssi,keskiarvo] = autokovis(d,deltaphi,askel)
%AUTOKOVIS [autokovarianssi,keskiarvo] = autokovis(d,deltaphi,askel)
%   Detailed explanation goes here
keskiarvo = real(mean(d)); % < d > eli keskiarvo
d11 = [d d]; % helpompi mennä "ympäri"

num = 0;

dd = zeros(1,length(d)./askel); % preallocate
for j = 1:askel:length(d);
  num = num + 1;
  dd(num) = d11(j).*d11(j+deltaphi);
end

autokovarianssi = real(((sum(dd))./num)-keskiarvo.^2);
end


function ang = kulma(v1,v2)
%KULMA function ang = kulma(v1,v2)
%   laskee vektorien v1 ja v2 välisen kulman (rad.)
ang = acos(dot(v1,v2)./((norm(v1)).*(norm(v2))));
end


function [dots,crosses,sivujen_pituudet] = vec_ops(coords)

sivujen_pituudet = sqrt((coords(1:(end-1),1)-coords(2:end,1)).^2+(coords(1:(end-1),2)-coords(2:end,2)).^2);
sivujen_pituudet(length(sivujen_pituudet)+1) = sqrt((coords(end,1)-coords(1,1)).^2+(coords(end,2)-coords(1,2)).^2);

vec1 = [(coords(1,1)-coords(end,1)),(coords(1,2)-coords(end,2)),0];
vec2 = [(coords(2,1)-coords(1,1)),(coords(2,2)-coords(1,2)),0];
dots(1) = dot(vec1,vec2)./(norm(vec1).*norm(vec2));
crosses(1,:) = cross(vec1,vec2)./(norm(vec1).*norm(vec2));

for j = 2:length(coords)-1
    vec1 = [(coords(j,1)-coords(j-1,1)),(coords(j,2)-coords(j-1,2)),0];
    vec2 = [(coords(j+1,1)-coords(j,1)),(coords(j+1,2)-coords(j,2)),0];
    dots(j) = dot(vec1,vec2)./(norm(vec1).*norm(vec2));
    crosses(j,:) = cross(vec1,vec2)./(norm(vec1).*norm(vec2));
end

vec1 = [(coords(end,1)-coords(end-1,1)),(coords(end,2)-coords(end-1,2)),0];
vec2 = [(coords(1,1)-coords(end,1)),(coords(1,2)-coords(end,2)),0];
dots(length(coords)+1) = dot(vec1,vec2)./(norm(vec1).*norm(vec2));
crosses(length(coords)+1,:) = cross(vec1,vec2)./(norm(vec1).*norm(vec2));
end


function [dmax,ip1,ip2] = diameter(x,y)
dmax=0;
if numel(x)>4 
    points = convhull(x,y)';
else
    points = 1:numel(x);
end
j = 1;
for i = points
    j = j + 1;
    x0 = x(i);
    y0 = y(i);
    for ii = points(j:end-1)
        d = sqrt((x0-x(ii)).^2+(y0-y(ii)).^2);
        if d>dmax
            dmax=d;
            ip1=i;
            ip2=ii;
        end;
    end
end
end


function filter2dirs(filename, inputDir, filter_code, plots, move)
if nargin<4
    plots = false;
end
[~,name,ext]=fileparts(filename);
filterDir = fullfile(inputDir,filter_code);
if(~exist(filterDir,'dir'))
    mkdir(filterDir);
end;
sourcepath = fullfile(inputDir,filename);
if move
    movefile(sourcepath,filterDir)
else
    copyfile(sourcepath,filterDir)
end
if plots
    plotpath = fullfile(inputDir,'classes',[name '_out' ext]);
    if move
        movefile(plotpath,filterDir)
    else
        copyfile(plotpath,filterDir)
    end
end
end


function oob = outOfBounds(imsize,edgepx)
edgepx_relat = edgepx/(2*sum(imsize));
oob = (edgepx > 55 && edgepx_relat > 0.10) || edgepx > 75 || edgepx_relat > 0.15;
end


function bb = blurryBoundaries(adiff,amax)
%disp(['adiff: ' num2str(adiff) ', adiff/amax: ' num2str(adiff/amax)])
bb = adiff/amax > 0.2 && adiff > 100;
end


function rb = rectBounds(imsize,xtremepx)
xtremepx_relat = xtremepx/(2*sum(imsize));
%rb = (xtremepx > 80 && xtremepx_relat > 0.15) || (xtremepx > 120) || (xtremepx_relat > 0.35);
rb = (xtremepx > 75 && xtremepx_relat > 0.12) || (xtremepx > 110) || (xtremepx_relat > 0.30);
end


function flags = filter_actions(filter_dirs, filename, inputDir, filter_code, plots, move, flags)
if filter_dirs && ~(~isempty(flags) && move)
    filter2dirs(filename, inputDir, filter_code, plots, move);
end
flags = flags_append(flags,filter_code);
end


function flags = flags_append(flags,newflag)
delimiter=';';
if isempty(flags)
    flags = newflag;
else
    flags = [flags delimiter newflag];
end
end

