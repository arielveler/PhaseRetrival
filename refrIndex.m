function n = refrIndex(mat, lambda)
% Returns refractive index of some known materials at a given wavelength
% Temperature deviation from the "standard" temperature (Attention! this 
% may be different for different materials) can be specified together with
% the material name, e.g. 'silica&+5.4K'.
%
% The material can be specified 
%   - by case-insensitive name, like 'silica' or 'air'
%   - or by the value of refractive index, in which case dispersion is
%     ignored and temperature deviation is not allowed
%
% Lambda is in nanometers (>50) or microns (<50)
%
% Currently known materials: SILICA, SILICACLAD, AIR, WATER, SM800CORE,
% ACETONE, ETHANOL, METHANOL, CARGILLE(133, 134, 135, 136, 137, 1.38, 1.39)
%
% Temperature coefficients supported for: SILICA, ACETONE, WATER, ETHANOL, 
% METHANOL, CARGILLE.
%
% EXAMPLE
%   RefrIndex('silica', 900)
%   RefrIndex('acetone&-0.5K')
%
% See also contents

% Copyright: 
% (cc-by) Fibres team, Meschede Research Group, Uni Bonn, 2010-2011
% kotya.karapetyan@gmail.com

if nargin == 1
    returnFunction = true; % wavelength is not specified, so return result as a function
    syms lambda
elseif nargin == 2
    returnFunction = false;
    if sum(lambda > 50) == length(lambda) % lambda is >50 so probably in nm
        lambda = lambda / 1000; % convert to microns
    else % wavelength is apparently in um, all lambdas should be <= 50 
        assert(sum(lambda > 50) == 0, 'Wavelength should be either in nm or um');
    end
else
    error('%s: Ivalid number of input parapeters\n', upper(mfilename));
end

% Parameters can be passed to refrIndex together with the material name.
% Currently, the only supported parameter is temperature deviation from the
% standard value. 
dt = 0;
if ~isa(mat, 'double')
    matWithParam = regexp(mat, '&', 'split');
    if numel(matWithParam) > 1 % parameters given
        param2string = matWithParam{2};
        assert(strcmpi(param2string(end), 'K'), ...
            'Invalid material specification');
        param2string(param2string == ',') = '.';
        try
            dt = eval(param2string(1:end-1));
        catch ME
            if strcmpi(ME.message, 'Error: Unexpected MATLAB expression.')
                error('Failed to obtain temperature difference from the material name\n');
            else
                rethrow(ME)
            end
        end
    end
    mat = matWithParam{1};
end

% For the case if the index n has to be different from 'silica',
% 'silicaclad' or 'air'.
if isa(mat, 'double') % in case the refractive index is specified by a double value
    if mat >= 1
        n = mat;
    else
        error('%s: Invalid refractive index specified\n', upper(mfilename));
    end
else % material or refr. index is specified by a string
    switch upper(mat)
        case 'SILICA'
            if returnFunction
                n(lambda) =  silica(lambda, dt);
            else
                n = silica(lambda, dt);
            end
        case 'SILICACLAD'
            if returnFunction
                n(lambda) =  silicaclad(lambda, dt);
            else
                n = silicaclad(lambda, dt);
            end
        case 'SM800CORE'
            if returnFunction
                n(lambda) =  sm800core(lambda, dt);
            else
                n = sm800core(lambda, dt);
            end
        case 'AIR'
            if returnFunction
                n(lambda) =  1;
            else
                n = 1 * ones(size(lambda));
            end
        case 'WATER'
            if returnFunction
                n(lambda) =  water(lambda, dt);
            else
                n = water(lambda, dt);
            end
        case 'ACETONE'
            if returnFunction
                n(lambda) =  acetone(lambda, dt);
            else
                n = acetone(lambda, dt);
            end
        case 'ETHANOL'
            if returnFunction
                n(lambda) =  ethanol(lambda, dt);
            else
                n = ethanol(lambda, dt);
            end
        case 'METHANOL'
            if returnFunction
                n(lambda) =  methanol(lambda, dt);
            else
                n = methanol(lambda, dt);
            end
        case 'CARGILLE130'
            if returnFunction
                n(lambda) =  cargille130(lambda, dt);
            else
                n = cargille130(lambda, dt);
            end
        case 'CARGILLE131'
            if returnFunction
                n(lambda) =  cargille131(lambda, dt);
            else
                n = cargille131(lambda, dt);
            end
        case 'CARGILLE132'
            if returnFunction
                n(lambda) =  cargille132(lambda, dt);
            else
                n = cargille132(lambda, dt);
            end
        case 'CARGILLE133'
            if returnFunction
                n(lambda) =  cargille133(lambda, dt);
            else
                n = cargille133(lambda, dt);
            end
        case 'CARGILLE134'
            if returnFunction
                n(lambda) =  cargille134(lambda, dt);
            else
                n = cargille134(lambda, dt);
            end
        case 'CARGILLE135'
            if returnFunction
                n(lambda) =  cargille135(lambda, dt);
            else
                n = cargille135(lambda, dt);
            end
        case 'CARGILLE136'
            if returnFunction
                n(lambda) =  cargille136(lambda, dt);
            else
                n = cargille136(lambda, dt);
            end
        case 'CARGILLE137'
            if returnFunction
                n(lambda) =  cargille137(lambda, dt);
            else
                n = cargille137(lambda, dt);
            end
        case 'CARGILLE138'
            if returnFunction
                n(lambda) =  cargille138(lambda, dt);
            else
                n = cargille138(lambda, dt);
            end
        case 'CARGILLE139'
            if returnFunction
                n(lambda) =  cargille139(lambda, dt);
            else
                n = cargille139(lambda, dt);
            end
        case 'HEPTANE'
            if returnFunction
                n(lambda) =  heptane(lambda, dt);
            else
                n = heptane(lambda, dt);
            end
        case 'OSCILLATOR'
            if returnFunction
                n(lambda) =  oscillator(lambda, dt);
            else
                n = oscillator(lambda, dt);
            end
        otherwise % expect refractive index as a string
            if isnan(str2double(mat)) % check if parseable
                error('%s: Cannot understand "%s" as material or refractive index', upper(mfilename), mat);
            end
            if returnFunction 
                n(lambda) =  str2double(mat);
            else 
                n = str2double(mat) * ones(size(lambda));
            end
    end
end

% assert(sum(n > 0) == numel(n(~isnan(n))), 'Refractive index < 0');

end % main function

%% silica
function n = silica(lambda, dt)
% LAMBDA is in meters. Calculates refractive index for fused silica.
%
% Coefficients taken from: Tong, Lou, Mazur. Single-mode guiding properties
% of subwavelength-diameter silica and silicon wire waveguides. Opt.
% Express 12 1025 (2004). http://dx.doi.org/10.1364/OPEX.12.001025
%
% Error not more than 0.00001 from 300 nm and beyond, compared to Melles
% Griot Synthetic Fused Silica Materials properties, and not more than
% 0.00003 from 400 nm compared to CVI Index of refraction Technical notes. 
%
% (KK 2009-03-07) Coefficients checked with Handbook of optics, p.
% II-33.69, table 23. The wavelength range given in HoO is 0.21...3.71 um.

if isa(class(lambda),'double') && lambda < 0.21e-6
    error('%s: The result cannot be guaranteed for such short wavelength (%g)\n', upper(mfilename), lambda);
end

x = lambda;

n = sqrt(1 + (0.6961663 * x.^2) ./ (x.^2 - 0.0684043e-6^2) +...
    (0.4079426 .* x.^2) ./ (x.^2 - 0.1162414e-6^2) +...
    (0.8974794 * x.^2) ./ (x.^2 - 9.896161e-6^2));

n = n + dt * 12E-6; % http://www.sciner.com/Opticsland/FS.htm, http://en.wikipedia.org/wiki/Fused_quartz
end % silica function

%% silica clad (F-doped silica)
function n = silicaclad(lambda, dt)
% LAMBDA is in the same units as in SILICA
% The base dispersion equation is for fused silica, and the
% SILICA function is used. The cladding is known to have lower refractive
% index by 0.36% [Corning, for telecom fibres:
% http://www.corning.com/WorkArea/showcontent.aspx?id=15535]

x = lambda;

n = silica(x, dt) / 1.0036;
end % silicaclad function

%% sm800core
function n = sm800core(lambda, dt)
% Refractive index of the core material of Fibercore SM800 fibre. This
% fibre has a nominal NA of 0.14 and, according to John Wooler (see email 
% to Uli Wiedemann from 2009-08-21), is within 0.138 to 0.145 for the range
% from ~370 to 1750 nm (0.138 to 0.140 for the whole NIR range), see the
% first picture in the email. 

x = lambda;
NA = 0.14;

n = sqrt(silica(x, dt) .^2 + NA^2) ;
end % sm800core function

%% water
function n = water(lambda, dt)
% LAMBDA is in meters. Calculates refractive index for water at 21.5 ן¿½C.
%
% Coefficients taken from Daimon 2007, Applied Optics Vol. 46, No. 18, p. 3814.

if isa(class(lambda),'double') && lambda < 0.18
    error('%s: The result cannot be guaranteed for such short wavelength (%g)\n', upper(mfilename), lambda);
end

x = lambda;

n = sqrt(1 + (0.5689093832 * x.^2) ./ (x.^2 - 0.005110301794e-6^2) +...
    (0.1719708856 .* x.^2) ./ (x.^2 - 0.018251801554e-6^2) +...
    (0.02062501582 .* x.^2) ./ (x.^2 - 0.02624158904e-6^2) +...    
    (0.1123965424 * x.^2) ./ (x.^2 - 10.67505178e-6^2));

n = n + dt * (-0.8E-4); % 10.1364/AO.12.001584
end % water function

%% acetone
function n = acetone(lambda, dt)
% LAMBDA is in microns. 
%
% Coefficients taken from http://refractiveindex.info and there from
% doi:10.1088/0957-0233/8/6/003 for 476.5...830 nm at 20 degrees

x = lambda;
C1 = 1.34979; C2 = 0.00306; C3 = 0.00006;
n = C1 + C2./x.^2 + C3./x.^4;

n = n + dt * (-5E-4); % 10.1364/AO.12.001584
end

%% ethanol
function n = ethanol(lambda, dt)
% LAMBDA is in microns. 
%
% Coefficients taken from http://refractiveindex.info and there from
% doi:10.1088/0957-0233/8/6/003 for 476.5...830 nm

x = lambda;
C1 = 1.35265; C2 = 0.00306; C3 = 0.00002;
n = C1 + C2./x.^2 + C3./x.^4;

n = n + dt * (-4E-4); % Condensed Matter Physics 2006, Vol. 9, No 2, p. 385
end

%% methanol
function n = methanol(lambda, dt)
% LAMBDA is in microns. 
%
% Coefficients taken from http://refractiveindex.info and there from
% doi:10.1016/S0921-4526(99)00856-X for 400...800 nm

x = lambda;
C1 = 1.294611; C2 = 12706.403e-6;
n = C1 + C2./x.^2;

n = n + dt * (-4E-4); % 10.1364/AO.12.001584
end

%% cargille 1.35
function n = cargille135(lambda, dt)
% Refractive index for refractive index oil Cargille Labs AAA 1.35.
% LAMBDA is in microns. 
%
% Data from SeriesAAA_13500.pdf sent to Kotya by cargillelabs@aol.com on
% 2011-07-07.
%
% Temp. Coefficient: -0.000339 dnD/dt ( 15 - 35 °C )

x = lambda * 1E4; % equation is written for angstroms
A = 1.3432154;
B = 237036;
C = -4.943692E+10;
n = A + B ./ x.^2 + C ./ x.^4;  % at 25 degrees

n = n + dt * (-3.39E-4); % ( 15 - 35 °C )
end

%% cargille 1.36
function n = cargille136(lambda, dt)
% Refractive index for refractive index oil Cargille Labs AAA 1.36.
% LAMBDA is in microns. 
%
% Data from SeriesAAA_13600.pdf sent to Kotya by cargillelabs@aol.com on
% 2011-07-07.
%
% Temp. Coefficient: -0.000341 dnD/dt ( 15 - 35 °C )

x = lambda * 1E4; % equation is written for angstroms
A = 1.3527514;
B = 254675;
C = -1.024360E+11;
n = A + B ./ x.^2 + C ./ x.^4; % at 25 degrees

n = n + dt * (-3.41E-4); 
end

%% cargille 1.37
function n = cargille137(lambda, dt)
% Refractive index for refractive index oil Cargille Labs AAA 1.36.
% LAMBDA is in microns. 
%
% Data from SeriesAAA_13700.pdf sent to Kotya by cargillelabs@aol.com on
% 2011-07-07.
%
% Temp. Coefficient: -0.000342 dnD/dt ( 15 - 35 °C )

x = lambda * 1E4; % equation is written for angstroms
A = 1.3622874;
B = 272314;
C = -1.554351E+11;
n = A + B ./ x.^2 + C ./ x.^4; % at 25 degrees

n = n + dt * (-3.42E-4); 
end

%% heptane (not implemented)
function n = heptane(lambda, dt) %#ok<STOUT,INUSD>
% LAMBDA is in microns. 
%
% Data taken from doi:10.1007/978-3-540-75291-2_253, for 20 degrees

error('REFRINDEX: Not yet implemented for HEPTANE');

x = lambda; %#ok<UNRCH>

lambda1 = 589.3; 
n1 = mean([1.3873 1.3872 1.38791 1.38776 1.38776 1.3876]);
lambda2 = 632.8; 
n2 = mean([1.3865 1.3865]);

% must still fit C1 and C2 
n = C1 + C2./x.^2;
end

%% cargille 1.33
function n = cargille133(lambda, dt)
% Refractive index for refractive index oil Cargille Labs AAA 1.36.
% LAMBDA is in microns. 
%
% Data from SeriesAAA_13300.pdf sent to Kotya by cargillelabs@aol.com on
% 2011-07-07.

x = lambda * 1E4; % equation is written for angstroms
A = 1.3241434;
B = 201757;
C = 5.656121E+10;
n = A + B ./ x.^2 + C ./ x.^4; % at 25 degrees

n = n + dt * (-3.37E-4); % within 15-35 degrees
end

%% cargille 1.34
function n = cargille134(lambda, dt)
% Refractive index for refractive index oil Cargille Labs AAA 1.36.
% LAMBDA is in microns. 
%
% Data from SeriesAAA_13400.pdf sent to Kotya by cargillelabs@aol.com on
% 2011-07-07.


x = lambda * 1E4; % equation is written for angstroms
A = 1.3336794;
B = 219396;
C = 3.562146E+09;
n = A + B ./ x.^2 + C ./ x.^4; % at 25 degrees

n = n + dt * (-3.38E-4); % ( 15 - 35 degrees C )
end

%% cargille 1.39
function n = cargille139(lambda, dt)
% Refractive index for refractive index oil Cargille Labs AAA 1.36.
% LAMBDA is in microns. 
%
% Data from SeriesAAA_13900.pdf sent to Kotya by cargillelabs@aol.com on
% 2011-07-07.

x = lambda * 1E4; % equation is written for angstroms
A = 1.3813595;
B = 307592;
C = -2.614332E+11;
n = A + B ./ x.^2 + C ./ x.^4; % at 25 degrees

n = n + dt * (-3.44E-4); % ( 15 - 35 degrees C )
end

%% cargille 1.38
function n = cargille138(lambda, dt)
% Refractive index for refractive index oil Cargille Labs AAA 1.36.
% LAMBDA is in microns. 
%
% Data from SeriesAAA_13800.pdf sent to Kotya by cargillelabs@aol.com on
% 2011-07-07.

x = lambda * 1E4; % equation is written for angstroms
A = 1.3718235;
B = 289953;
C = -2.084341E+11;
n = A + B ./ x.^2 + C ./ x.^4; % at 25 degrees

n = n + dt * (-3.43E-4); % ( 15 - 35 degrees C )
end

%% cargille 1.30
function n = cargille130(lambda, dt)
% Refractive index for refractive index oil Cargille Labs AAA 1.36.
% LAMBDA is in microns. 
%
% Data from SeriesAAA_13000.pdf sent to Kotya by cargillelabs@aol.com on
% 2011-07-07.

x = lambda * 1E4; % equation is written for angstroms
A = 1.2955353;
B = 148840;
C = 2.155584E+11;
n = A + B ./ x.^2 + C ./ x.^4; % at 25 degrees

n = n + dt * (-3.33E-4); % ( 15 - 35 degrees C )
end

%% cargille 1.31
function n = cargille131(lambda, dt)
% Refractive index for refractive index oil Cargille Labs AAA 1.36.
% LAMBDA is in microns. 
%
% Data from SeriesAAA_13100.pdf sent to Kotya by cargillelabs@aol.com on
% 2011-07-07.

x = lambda * 1E4; % equation is written for angstroms
A = 1.3050713;
B = 166479;
C = 1.625594E+11;
n = A + B ./ x.^2 + C ./ x.^4; % at 25 degrees

n = n + dt * (-3.34E-4); % ( 15 - 35 degrees C )
end

%% cargille 1.32
function n = cargille132(lambda, dt)
% Refractive index for refractive index oil Cargille Labs AAA 1.36.
% LAMBDA is in microns. 
%
% Data from SeriesAAA_13200.pdf sent to Kotya by cargillelabs@aol.com on
% 2011-07-07.

x = lambda * 1E4; % equation is written for angstroms
A = 1.3146073;
B = 184118;
C = 1.095603E+11;
n = A + B ./ x.^2 + C ./ x.^4; % at 25 degrees

n = n + dt * (-3.36E-4); % ( 15 - 35 degrees C )
end
