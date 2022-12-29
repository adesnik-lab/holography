function [ Sutter ] = function_Sutter_Start( Setup )
%FUNCTION_SUTTER_START Summary of this function goes here
%   Detailed explanation goes here

Sutter.obj = sutterMP285(Setup.Sutterport);
[stepMult, currentVelocity, vScaleFactor] = getStatus(Sutter.obj);
setVelocity(Sutter.obj, 50, vScaleFactor); 
Sutter.obj.verbose=0;
end

