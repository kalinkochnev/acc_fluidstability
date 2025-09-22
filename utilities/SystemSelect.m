% Indicate which system matrix to use to solve the system for.
% RADKO_SHEAR_INSTABILITY - 4x4 matrix
% RADKO_GRL - 5x5 matrix
classdef SystemSelect
   enumeration
       RADKO_SHEAR_INSTABILITY, RADKO_GRL, UNSELECTED
   end
end
