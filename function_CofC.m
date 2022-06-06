function [ Get ] = function_CofC( refAsk,refGet,Ask )
% Ref ask Ref get, size N by 2 or N by 3

[LN LP] = size(refAsk);


if LP == 1
modelterms = [0 ; 1 ; 2 ];     %XY spatial calibration model for C_Of_C between SLM and true space
FitX =  polyfitn(refAsk,refGet(:,1),modelterms);


GetX = polyvaln(FitX,Ask);
Get = [GetX];

elseif LP == 2
modelterms = [0 0 ; 1 0 ; 0 1 ; 1 1];     %XY spatial calibration model for C_Of_C between SLM and true space
FitX =  polyfitn(refAsk,refGet(:,1),modelterms);
FitY =  polyfitn(refAsk,refGet(:,2),modelterms);

GetX = polyvaln(FitX,Ask);
GetY = polyvaln(FitY,Ask);
Get = [GetX GetY];

elseif LP == 3;

modelterms = [0 0 0 ; 1 0 0 ; 0 1 0; 0 0 1;  1 1 0 ; 1 0 2];     %XYZ spatial calibration model for C_Of_C between SLM and true space
FitX =  polyfitn(refAsk,refGet(:,1),modelterms);
FitY =  polyfitn(refAsk,refGet(:,2),modelterms);
FitZ =  polyfitn(refAsk,refGet(:,3),modelterms);
GetX = polyvaln(FitX,Ask);
GetY = polyvaln(FitY,Ask);
GetZ = polyvaln(FitZ,Ask);

Get = [GetX GetY GetZ];    
   
else
    disp('you idfiot')
end
end