###########################################################
# Idealized Coronary Bypass Generated by Bypass.m
# Original File By M.R. De Luca and A. Veneziani - 2010
# Modified By Mahdi Al-Husseini for Coronary Bypass - 2017
###########################################################

algebraic3d

#################### Pulmonary Artery ####################
solid pa = cylinder (-35.000000,0,0;35.000000,0,0;5.000000)
and plane (-34.000000, 0, 0; -1.000000, 0, 0) 
and plane (34.000000, 0, 0; 1.000000, 0, 0); 

#################### Superior Vena Cava (SVC) ###################
solid svc = torus (0,15.090000,0;0,0,1;15.000000;4.900000)
and plane (0, 0, 0; -1, 0, 0) 
and plane (0, 15.090000, 0; 0, 1, 0); 

#################### SVC Extension ####################
solid svcE = cylinder (15.000000,15.090000,0;15.000000,25.090000,0;4.900000)
and plane (0, 25.090000, 0; 0, 1, 0) 
and plane (0, 15.090000, 0; 0, -1, 0); 
##################### Inferior Vena Cava (IVC) #####################
solid ivc = cone (-10.000000,0,0;5.000000;-10.000000,-40.000000,0;10.000000)
and plane (-10.000000, -10.000000, 0; 0, -1, 0) 
and plane (-10.000000, 0.1, 0; 0, 1, 0); 

solid lumen = svcE or pa or svc or ivc;
tlo lumen; -maxh=0.010000
