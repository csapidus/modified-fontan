%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original File by M. De Luca & A. Veneziani 2010                         %
% Modified by Mahdi Al-Husseini, 11 October 2017 [1.4],                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ribint,Ribext, Rpvext, pvspan, RlwextL, RlwextS, lvspan, ...
    offset, extensionspan, hmesh]= BypassPreproc(Ribint, Ribext, Rpvext, ...
    pvspan, RlwextL, RlwextS, lvspan, offset, extensionspan, hmesh)

%bypassPreproc prepares the geometrical data for bypass to generate
%parametric geometries of an idealized coronary bypass
% Some check on the consistency of the data is performed and bypass is run
% upon request of the user.
%
%% Input Parameters
% Miscellaneous **********************************************************
% hmesh = max size of the dsired mesh
% offset = lateral off-set between SVC and IVC
% Superior Vena Cava (SVC) **************************************************************
% Rib = radius of curvature of SVC
% Ribext = external radius of SVC
% Ribint = internal radius of SVC
% extensionspan = span of SVC extension
% Pulmonary Artery (PA) ******************************************************* 
% Rpvext = external radius of PA
% pvspan = span of PA
% Inferior Vena Cava (IVC) ********************************************************
% RlvextL = larger radius of the IVC
% RlvextS = smaller radius of the IVC
% lvspan = span of the IVC

% BypassPreproc(Ribint, Ribext, Rpvext, ...
%    pvspan, RlwextL, RlwextS, lvspan, offset, extensionspan, hmesh)
% As Requested by Dr.Corno: BypassPreproc(4.9, 15, 5, 70, 10, 5, 40, 10, 10, 0.01)

%% Input Validation
counter = 0;
if (offset - RlwextS > (pvspan/2) +1), ...
    disp('ERROR: offset is beyond bounds of pulmonary artery');
    counter = counter +1;
end

if (Ribint > Ribext), ...
    disp('ERROR: internal torus radius cannot be larger than external torus radius');
    counter = counter +1;
end

if (RlwextS > Rpvext), ...
    disp('ERROR: smaller radius of inferior vena cava may not have a larger radius than pulmonary artery');
    counter = counter +1;
end

if (Ribext > pvspan/2), ...
    disp('ERROR: pulmonary artery span too short for connecting geometry');
    counter = counter +1;
end

if (counter > 0)
    return;
end

%% Bypass.m Initiation
R=input('Do you want to run bypass? (y/n)','s');

if (R=='y')
    condition = true;
    while(condition == true)
        condition = false;
        F = input('Would you like a fontan model or the alternative? (f/a)','s');
        if (F == 'f')
            Bypass('bypass-fontan',Ribint, Ribext, Rpvext, pvspan, RlwextL, RlwextS, lvspan, ...
                offset, extensionspan, hmesh, 'fontan');
        elseif (F == 'a')
            Bypass('bypass-alternative',Ribint, Ribext, Rpvext, pvspan, RlwextL, RlwextS, lvspan, ...
                offset, extensionspan, hmesh, 'alternative');
        else
            O = input('That is not an acceptable input. Would you like to try again? (y/n)', 's');
            if (O == 'y')
                condition = true;
            end
        end
    end
end

end