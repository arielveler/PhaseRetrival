clear; clc; close all;
% Example - fused silica:
material = 'Fused silica';

%lambda0 = 0.21e-6:0.01e-6:3.71e-6;
lambda = 0.5e-6:0.001e-6:3.71e-6;
%Sellmeier's coefficients - fused silica:
B1 = 0.6961663;
B2 = 0.4079426;
B3 = 0.8974794;
C1 = 0.0684043e-6^2; %[um^2]
C2 = 0.1162414e-6^2; %[um^2]
C3 = 9.896161e-6^2; %[um^2]

syms lam 
c = 2.99792458e8; %[m/s]

%Sellmeier's formula:
%n=sqrt(1+0.6961663./(1-(0.0684043./x).^2)+0.4079426./(1-(0.1162414./x).^2)+0.8974794./(1-(9.896161./x).^2))
n(lam) = sqrt(B1*lam.^2./(lam.^2-C1)+B2*lam.^2./(lam.^2-C2)+B3.*lam.^2./(lam.^2-C3)+1);
%double(n(850e-9))
vg(lam) = c./(n-diff(n,lam,1)*lam); % pulse average speed in medium
% disspersion by GVD:
GVD(lam) = (lam^3/(2*pi*c^2))*diff(n, lam, 2); %GVD coefficient
%double(GVD(850e-9)) %check - 32.228 [fs2/mm] for fused silica at 850nm (32.228e-27)
D(lam) = -(lam/c)*diff(n, lam, 2); %Dissperssion coefficient D 
%double(D(850e-9)) %check - -84.024 ps/(nm km) for fused silica at 850nm (-84.024e-6)
lambda_zd = lambda(abs(double(D(lambda)))==min(abs(double(D(lambda))))); %zero-disspersion wavelength

figure(1);
subplot(1,3,1); plot(lambda,double(n(lambda)),'linewidth',2);
title('n(\lambda_0)'); xlabel('Wavelenght'); ylabel('n');
subplot(1,3,2); plot(lambda,double(GVD(lambda)),'linewidth',2);
title('GVD(\lambda_0)'); xlabel('Wavelenght'); ylabel('GVD [fs^2/mm]');
subplot(1,3,3); plot(lambda,double(D(lambda)),'linewidth',2);
title('D_{Material}(\lambda_0)'); xlabel('Wavelenght'); ylabel('D [ps/(nm*km)]');
sgtitle(['n, GVD, D as function of \lambda_0 from 3 term sellmeir eqation for ',material]);

%% calculation for specific fibers:
% fibers course broadening at 850[nm]:
% parameters:
NA = 0.11; % of F-CPL-M18850
core_diam = 9e-6; % of F-CPL-M18850
lam_0 = 1550e-9; % [m] central wavelength of laser NPL52C
%d_lam = 1.2e-9; % Spectral bandwidth of NPL52C
d_lam=0.2e-9;
% lam_0 = 1030e-9; % [m] central wavelength of Pharos
% d_lam = 5e-9; % Spectral bandwidth of Pharos fot 300[fs] pulse
%d_lam = linspace(lam_0-delta_lambda/2,lam_0+delta_lambda/2,10); % [m] FWHM of pulse for lam_0=500[nm]
%d_lam = linspace(0.028e-9,0.030e-9,10); % [m] FWHM of pulse for lam_0=850[nm]
%d_lam = linspace(0.035e-9,0.05e-9,10); % [m] FWHM of pulse for
%lam_0=1030[nm]

fibers_num = 6;
delay = 2; %how many "quiet" pulses there will be between 2 pulses
L = 1; % [L] initial fiber

%check if enougth modes:
V = (2*pi/lam_0)*core_diam*NA; % V number
M = 4*V^2/pi^2; % approximation of no. of modes

t_pulse_FWHM = zeros(numel(d_lam),1);

for j=1:numel(d_lam)
    t_pulse_FWHM(j) = dlam2dt(d_lam(j),lam_0,'gaussian');
    
    % Simulation of the ideal fibers number and length according to dispersion:
    L_fiber = zeros(fibers_num,1);
    t_pulse_tot = zeros(fibers_num,1);

    L_fiber(1) = L; %length of the first fiber
    for i=1:fibers_num
       %[dist_max,dt_NA] = calc_NA_dispersion(NA,L_fiber(i),n,lam_0,vg); 
       t_pulse_1 = calc_GVD_dispersion(double(GVD(lam_0)),L_fiber(i),d_lam(j),lam_0);
       %t_pulse_2 = calc_GVD_dispersion(double(GVD(lam_0)),dist_max,d_lam(j),lam_0);
       %t_pulse_tot(i) = dt_NA+t_pulse_1/2+t_pulse_2/2;
       t_pulse_tot(i) = t_pulse_1;
       if i==fibers_num 
           break;
       end
       %L_fiber(i+1)=L_fiber(i)+t_pulse_tot(i)*double(vg(lam_0))*delay;
       L_fiber(i+1)=ceil(L_fiber(i)+t_pulse_tot(i)*double(vg(lam_0))*delay);
    end
    L_max(j) = L_fiber(end); %longest fiber length
    t_pulse_max(j) = t_pulse_tot(end);
    factor(j) = t_pulse_max(j)/t_pulse_FWHM(j);
    figure(1)
    plot([t_pulse_FWHM(j);t_pulse_tot],[0;L_fiber],'-o','DisplayName',num2str(d_lam(j)))
    hold on
end
ylabel('fibers length [m]');
xlabel('pulse times [s]');

% giving a score to that laser:
% figure(2)
% plot(d_lam.*1e9,factor,'-o',t_pulse_FWHM*1e12,factor,'-o','DisplayName',['\lambda_0=',num2str(lam_0)]);
% xlabel('${\Delta}{\lambda}/{\Delta}{t} - \textrm{FWHM [nm]/[ps]}$','Interpreter','latex','FontSize',16);
% ylabel('$\textrm{Factor} - L_{max}{\cdot}{\frac{{\Delta}t-Pulse_{max}}{{\Delta}t-Pulse_{initial}}}$','Interpreter','latex','FontSize',16);
% legend({'\delta\lambda [pm]','\deltat [ps]'},'FontSize',12)
% title(['\lambda_0 - ',num2str(lam_0*1e9), '[nm]'],'FontSize',16);

%% attempts
lam_0=520e-9;
d_lam_FWHM=0.00006e-9;
L=1;
%gvd = double(GVD(lam_0));
%t_pulse=calc_GVD_dispersion(gvd,L,d_lam_FWHM,lam_0)
t_pulse_FWHM = dlam2dt(d_lam_FWHM,lam_0,'gaussian')
%% Modal disperssion for splitting transvers modes:
clear; clc; close all;
material_core = 'SILICA';
material_clad = 'SILICACLAD';
%L_fiber = 100; % [m]
lam0 = 1030e-9; %central wavwlength
factor = 1e6; %factor for refractive index function
k0 = 2*pi/lam0;
c = 2.99792458e8; %[m/s]

to_plot = 4; % 1 - plot b(v). 2 - plot dvb/dv(v). 0 - do not plot

syms lambda

N1 = refrIndex(material_core);
N2 = refrIndex(material_clad)/0.9994;
n1 = double(N1(lam0));
n2 = double(N2(lam0));

NA = sqrt(n1^2-n2^2); %NA of fiber

d = 20e-6; % core diameter 
a = d/2; % core radius

V = (2*pi/lam0)*a*NA; % V number
V_cutoff = @(m,l)m*pi+(l-3/2)*(pi/2); 
M = V^2/4; % approximation of no. of modes

l_num = 9;
m_range = 12;
del = (n1-n2)/n2;
if to_plot
    figure('units','normalized','outerposition',[0 0 1 1]);
%     y_lim = 1*to_plot;
%     ylim([0 y_lim]);
    hold on
end
UC = load('UC.mat'); % loading zeros of bessele(l-1,v)
UC = UC.UC;

for l=0:l_num
    syms v
    uc = UC{l+1};
    for m=1:numel(uc)
        if uc(m) >= m_range
            continue
        end
        w(v) = sqrt(v^2-uc(m)^2);
        k(v) = 1-(w^2+l^2+1)^(-1/2);
        s = sqrt(uc(m)^2-l^2-1);
        if l==0 && m==1
            u(v) = (1+sqrt(2))*v/(1+(4+v^4)^0.25);
        else
            u(v) = uc(m)*exp((asec(v/s)-asec(uc(m)/s))/s);
        end
        b(v) = 1-(u(v)/v)^2;
        if to_plot == 1
            fplot(b,[0 m_range]);
            str = ['LP_ ',[num2str(l) ',' num2str(m)]];
            dim = [0.1333+((m_range-(m_range-uc(m))/2)/m_range)*0.75,0.125+(double(b(m_range-(m_range-uc(m))/2)))*0.8125, 0 , 0];
            annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','white');
            xlabel('Fiber parameter - v'); ylabel('Normalized propagation parameter b(v)');
        end
        
        dvb_dv(v) = 1-(u(v)/v)^2*(1-2*k(v));
        if to_plot == 2
            fplot(dvb_dv,[uc(m) m_range]);
            str = ['LP_ ',[num2str(l) ',' num2str(m)]];
            dim = [0.1333+((m_range-(m_range-uc(m))/2)/m_range)*0.75,0.125+(double(dvb_dv(m_range-(m_range-uc(m))/2)))*0.8125/y_lim, 0 , 0];
            annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','white');
            xlabel('Fiber parameter - v'); ylabel('Normalized group delay d(vb)/dv (v)');
        end
        
        Neff =  N2.*(double(b(V)).*((N1-N2)./N2)+1);
        gvd_tmp = (lambda^3/(2*pi*c^2))*diff(Neff, lambda, 2);

        if isreal(double(w(V)))
            tau_modal(l+1,m) = n2*del*(double(dvb_dv(V)));
            gvd(l+1,m) = double(gvd_tmp(lam0));
            neff(l+1,m) = n2*(double(b(V))*del+1);
        end
    end
    if l == 0
        tau_modal = [tau_modal; NaN(l_num-1,numel(tau_modal))];
        gvd = [gvd ; NaN(l_num-1,numel(gvd))];
        neff = [neff ; NaN(l_num-1,numel(neff))];
    end 
end
tau_material = N2-lambda*diff(N2, lambda, 1); %d_(nk)/d_k = n(lam)-lam*d_n/d_lam

vg_total = c./(tau_modal+double(tau_material(lam0)));

% calculations:
syms L_fiber
aa = vg_total(:); 
aa(isnan(aa)) = 0;
[aa,indx] = sort(aa,'descend');
bb = aa(aa>0);
bb = 1e12./bb - 1e12./max(bb); %delay in [ns/km] from the first mode

tau_min(L_fiber) = L_fiber.*1e-3.*min(bb(2:end)-bb(1:end-1)); % minimum time separating the closest modes 

% disspersion by GVD:
GVD(lambda) = (lambda^3/(2*pi*c^2))*diff(N1, lambda, 2); %GVD coefficient
D(lambda) = -(lambda/c)*diff(N1, lambda, 2); %Dissperssion coefficient D 

d_lam=0.01e-9;
t_pulse_FWHM_before = dlam2dt(d_lam,lam0,'gaussian');
t_pulse_FWHM_after(L_fiber) = calc_GVD_dispersion(double(GVD(lam0)),L_fiber,d_lam,lam0).*1e9; % pulse width in [ns] after fiber
L_range = 1000;


if to_plot == 3
    fplot(tau_min,[0 L_range]);
    fplot(t_pulse_FWHM_after,[0 L_range]);
    xlabel('Fiber Length - L'); ylabel('Time - t [ns]');
    legend('time difference between closest modes','pulse length at end of fiber');
end
if to_plot == 4
    plot(1:numel(bb),bb);
    xlabel('number of mode'); ylabel('Delay [ns/km] from 1st mode');
    title('Delay in [ns/km] from the first mode');
    %legend('time difference between closest modes','pulse length at end of fiber');
end
%%
% function b_lm = find_LP_modes(l, n1, n2, lam0, a)
% 
% k0 = 2*pi/lam0;
% syms n_eff
% beta(n_eff) = n_eff*k0;
% h = sqrt(n1^2*k0^2-beta^2);
% q = sqrt(beta^2-n2^2*k0^2);
% eqn = h*besselj(l+1,h*a)/besselj(l,h*a) == q*besselk(l+1,q*a)/besselk(l,q*a);
% 
% b_lm = findzeros(eqn,[n2 n1]);
% b_lm = (b_lm-n2)/(n1-n2)
% end
%     
%% calculation for lens, BS etc. of system:
clear;clc; close all;
% parameters to control:
pix = [2048 2048]; % in acA2440-35um - Basler ace there are 2448x2048
pixsize = [3.45e-6 3.45e-6]; % for acA2440-35um - Basler ace
margins = 0.1;
res = 135e-6 %resolution  on object
spots = [5 5];
overlap = 0.7;

lam0 = 1030e-9;

%rmax = pix.*pixsize.*(1-margins)./spots/(2*(1+margins));
rmax = pix.*pixsize.*(1-margins)^2./spots/2;
f = min(rmax)*res/lam0

dkmin = pixsize./lam0./f;

angle = 2*atand((spots-1).*rmax./(1-margins)./f) % total angle of the BS component
% FOV1 is collimated ray diameter on the object
FOV1 = 1./dkmin./2;   % divide by 2 - because of intensity, so a square on the camera may see only half of the spacial frequencies
% FOVtot is total area covered by the rays on the object:
Fovtot = FOV1.*(1+(spots-1).*(1-overlap));

%% Code to save besselj zeros for determined l,m 
clear; clc; close all;
l_num = 10; % set to determind the number of besselj function
m_range = 20; % set to determind the maximum number for 
syms v
for l=0:l_num  
    disp(['number ',num2str(l),'/',num2str(l_num)]);
    eqn1 = besselj(l-1,v)==0;
    if l==0
        uc = findzeros(besselj(l-1,v),[0 m_range]);
    else
        uc = findzeros(besselj(l-1,v),[eps*1e4 m_range]);
    end
    UC{l+1} = uc;
end
save('uc.mat','UC')
%%
% GVD dispersion calculation
function t_pulse_FWHM = dlam2dt(d_lam_FWHM,lam_0,method)
% Help function for converting a pulse from spectral width (in [m], FWHM)
% into temporal width (in [s], FWHM).
% input:
% d_lam - spectral bandwidth [m]
% lam_0 - mean wavelength [m]
% method - type of pulse ('gaussian','sech',lorenzian')
% output:
% t_pulse_FWHM - temporal pulse width [s]
    c = 2.99792458e8; %[m/s]
    d_omega = d_lam_FWHM*2*pi*c/lam_0^2; % FWHM bandwidth
    d_nu = d_omega/(2*pi);
    switch method
        case 'gaussian'
            time_bandwidth_const = 2*log(2)/pi; % for a gaussian
        case 'sech'
            time_bandwidth_const = 4*log(sqrt(2)+1)^2/pi^2;  %for a sech^2
        case 'lorentzian'
            time_bandwidth_const = log(2)*sqrt(sqrt(2)-1)/pi; %for a lorentzian
        otherwise
            time_bandwidth_const = 2*log(2)/pi; % gaussian as default
    end
    t_pulse_FWHM = time_bandwidth_const/d_nu; % [s] temporal FWHM of pulse

end


function t_FWHM_new=calc_GVD_dispersion(GVD,L,d_lam_FWHM,lam_0)
% function for calculation of the GVD dispersion
%input:
% D - dispersion coefficient [ps/nm*km]
% L - fiber length [m]
% d_lam_FWHM - spectral bandwidth FWHM [m]
% t_init_FWHM - pulse width FWHM [s]
% output:
% dt - addition time for the pulse due to material dispersion (GVD) 
    %t_pulse = t_init_FWHM + abs(D)*L*d_lam;
    
    %FWHM broadening:
%     c = 2.99792458e8; %[m/s]
%     d_omega_FWHM = d_lam_FWHM*2*pi*c/lam_0^2; % FWHM bandwidth
    fac = 2*sqrt(log(2)); % FWHM = fac*sqrt(2)*std - for a gaussian (sqrt(2)*std=1/e amplitude)
    %std broadening:
    std_pulse = dlam2dt(d_lam_FWHM,lam_0,'gaussian')/fac;
    std_pulse_new = sqrt(std_pulse^2+((GVD*L)/std_pulse)^2); %see https://www.ece.rutgers.edu/~orfanidi/ewa/ch03.pdf
    t_FWHM_new = std_pulse_new*fac;
end

% NA dispersion calculation ("Modal dispersion")
function [L_max,dt_NA] = calc_NA_dispersion(NA,L,n,lam_0,vg)
    factor_dist = cos(asin(NA./double(n(lam_0))));
    L_max = L./factor_dist; % maximum fiber length according to wavelenth
    
    %option a - without difference in ligth speed
    t_max = L_max./double(vg(lam_0)); %time of fligth maximum according to wavelenth
    t_min = L./double(vg(lam_0)); %
    dt_NA = t_max-t_min;
    
    % %option b - with difference in ligth speed
    % lam_min = lam_0-d_lam/2; lam_max = lam_0+d_lam/2; % min and max wavelength.
    % c_min = c./double(n(lam_min)); c_max = c./double(n(lam_max));  
    % t_max = L_max./c_min; %time of fligth maximum according to wavelenth
    % t_min = L./c_max; %
    % dt_NA = t_max-t_min
end
