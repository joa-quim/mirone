function Angles = dec2deg(vector_decimal_angles, opt)
%Converts decimal angles to angles in degrees
%
%SYNOPSIS:  Angle_degrees = dec2deg(Angle_decimal,opt)
%      where Angle_degrees is the output Matrix in degrees, minutes
%      and seconds and Angle_decimal is the input vector of decimal angles.
%      OPT is a flag (its value is irrelevant) to signal that output  
%      is an array with 3 columns [deg minutes seconds]. Otherwise, the
%      output is an array of strings with the form deg:min:sec

if (nargin == 1)
    opt = [];
    Angles = zeros(length(vector_decimal_angles),3);
end
for i = 1:length(vector_decimal_angles)  
   decimal_angle = vector_decimal_angles(i) ;
   degrees = fix(decimal_angle) ;
   minutes = fix(60*(decimal_angle - degrees)) ;
   seconds = 3600*(decimal_angle - degrees) - 60*minutes ;
   if seconds > 60 , seconds = seconds - 60 ;
            minutes = minutes + 1;
   end  
   diff = 60 - seconds ;
   if diff < 1e-10
       seconds = 0;
        minutes = minutes +  1;
   end
   if (~isempty(opt))
       Angles(i,:) = [degrees abs(minutes) abs(seconds)];
   else
       Angles(i,:) = [num2str(abs(degrees)) ':' num2str(abs(minutes)) ':' num2str(abs(seconds),'%.3f')];
   end
end
