% script to read .olympex flight data which contain information about
% aircraft speed, trajectory and environmental conditions
%
% some interesting entries:
% 1  = time in seconds from midnight, day aircraft flight started
% 2  = air temperature
% 6  = pressure altitude
% 8  = static pressure
% 10 = RH
% 33 = LWC based in King Probe
% 34 = Total Water Content based on Nevzorov
% 35 = LWC based on Nevzorov
% 36 = IWC based on Nevzorov
% 39 = Number concentration of droplets based on cloud droplet probe #/cc
% 40 = LWC based on cloud droplet probe
% 48 = 2D-S Horizontal Total Normalize Particle Concentration > 105 microns (11 pixels) #/m3
% 49 = 2D-S Horizontal Total Normalize Particle Concentration all bins #/m3
% 50 = same as 48 for Vertical channel
% 51 = same as 49 for Vertical channel
% 52 = HVPS3 horizontal total normalized particle concentration for all bins #/m3
% 53 = HVPS3 vertical total normalized particle concentration for all bins #/m3
%
% ... 


clearvars; close all;

filename = '../../OAP_flight_data/20151201_WF/Env/15_12_01_22_29_17.olympex';
date_ref = [2015 11 12 00 00 00];

fid = fopen(filename);
data = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',70,'Delimiter','\t');  

dt =  data{1};
dt_duration = seconds(dt);

flight.time = datetime(date_ref) + dt_duration;
flight.alt = data{6};

env.T = data{2};
env.p = data{8};
env.RH = data{10};
env.LWC_king = data{33};
env.LWC_nev = data{35};
env.LWC_clp = data{40};

% replacing N/A values by NaNs
env.LWC_king(env.LWC_king > 1e5) = NaN;
env.LWC_nev(env.LWC_nev > 1e5) = NaN;
env.LWC_clp(env.LWC_clp > 1e5) = NaN;

OAP.DSH_filt = data{48};
OAP.DSH_all = data{49};
OAP.DSV_filt = data{50};
OAP.DSV_all = data{51};
OAP.HVPSH = data{52};
OAP.HVPSV = data{53};

OAP.Dmean_cprobe = data{40};
OAP.Dmean_2DC = data{45};

% replacing N/A values by NaNs
OAP.Dmean_2DC(OAP.Dmean_2DC >= 1e3) = NaN;

% illustration
figure; subplot(311); hold on; grid on; box on;
title('Flight #2 : 2015.12.01');
yyaxis left;
plot(flight.time,flight.alt,'-','linewidth',2);
ylabel('Altitude [m agl]');

yyaxis right;
plot(flight.time,env.T,'-','linewidth',2);
ylabel('Temp. [°C]');


subplot(312); hold on; grid on; box on;
yyaxis left;
plot(flight.time,env.RH,'-');
ylabel('Rel. Hum. [%]');

yyaxis right;
plot(flight.time,env.LWC_king,'-');
ylabel('LWC [g/m^3]');

subplot(313); hold on; grid on; box on;

yyaxis right;
%plot(flight.time,OAP.Dmean_cprobe,'-','linewidth',1.5);
plot(flight.time,OAP.Dmean_2DC,'-','linewidth',0.5);
ylabel('Dmean [um]');
legend('Dmean');

yyaxis left;
plot(flight.time,OAP.DSV_filt,'r-','linewidth',1.5);
plot(flight.time,OAP.HVPSV,'b-','linewidth',1.5);
ylabel('part. concentration [#/m3]');
set(gca,'yscale','log','ylim',[1e0 1e7]);
legend('2D-S','HVPS');





% figure;
% subplot(121);
% scatter(OAP.HVPSV,OAP.DSH_all);
% subplot(122);
% scatter(OAP.HVPSV,OAP.DSH_filt);
% 
% 
% figure; hold on; grid on; box on;
% plot(flight.time,env.LWC_king,'k-');
% plot(flight.time,env.LWC_nev,'b-');
% plot(flight.time,env.LWC_clp,'r-');
% xlabel('Time UTC');
% ylabel('LWC [g/m^3]');













