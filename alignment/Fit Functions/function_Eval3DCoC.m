function Get = function_Eval3DCoC(CoC,Ask)

FitX = CoC.FitX;
FitY = CoC.FitY;
FitZ = CoC.FitZ;


GetX = polyvaln(FitX,Ask);
GetY = polyvaln(FitY,Ask);
GetZ = polyvaln(FitZ,Ask);

Get = [GetX GetY GetZ]; 