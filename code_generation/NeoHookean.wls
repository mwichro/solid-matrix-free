#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*\[FilledSquare] 3D Element (compressible neo-hookean)*)


(* ::Input:: *)
(*<<AceGen`;*)


dim = 2;
ScalarProduct2[U_, V_]=Tr[U . Transpose[V]];


energy:=
(
{Ee,\[Nu]}\[DoubleRightTee]{2500,0.35};
\[Lambda]\[DoubleRightTee](Ee*\[Nu])/((1+\[Nu])(1-2\[Nu]));
\[Mu]\[DoubleRightTee]Ee/(2(1+\[Nu]));
\[DoubleStruckCapitalH]=gradU;
Ii\[DoubleRightTee]IdentityMatrix[dim];
\[DoubleStruckCapitalF]e\[DoubleRightTee]SMSFreeze[Ii+\[DoubleStruckCapitalH]];
\[DoubleStruckCapitalC]e\[DoubleRightTee]Transpose[\[DoubleStruckCapitalF]e] . \[DoubleStruckCapitalF]e;
Je\[DoubleRightTee]SMSSqrt[Det[\[DoubleStruckCapitalC]e]];
bx=0;
by= 1000;
bz= 0;
Force= If[dim==2, {bx,by}, {bx,by,bz}];
\[Psi]el\[DoubleRightTee]\[Mu]/2 (Tr[\[DoubleStruckCapitalC]e]-dim-2Log[Je])+\[Lambda] (Log[Je])^2 + u . Force; 
)


residual :=
  (
    Ru \[DoubleRightTee] SMSD[\[Psi]el, u];
    Rgradu \[DoubleRightTee] SMSD[\[Psi]el, gradU];
  )


tangent :=
  (
    residual;
    Ku \[DoubleRightTee] SMSD[Ru, u];
    Kgradu \[DoubleRightTee] SMSD[Ru, gradU];
    Hu \[DoubleRightTee] SMSD[Rgradu, u];
    Hgradu \[DoubleRightTee] SMSD[Rgradu, gradU];
  )


SMSInitialize["Residual", "Language" -> "C++", "VectorLength" -> 100,
   "Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["Residual", Real[uIn$$[dim], graduIn$$[dim,dim], valueOut$$[dim], gradientOut$$[
  dim,dim]]];
u \[RightTee] Table[SMSReal[uIn$$[i]], {i, 1, dim}];

gradU \[RightTee] Table[SMSReal[graduIn$$[i,j]], {i, 1, dim}, {j, 1, dim}];

energy;

Ru \[DoubleRightTee] SMSD[\[Psi]el, u];
Rgradu \[DoubleRightTee] SMSD[\[Psi]el, gradU];
SMSExport[Ru, valueOut$$];
SMSExport[Rgradu, gradientOut$$];
SMSWrite[];



SMSInitialize["Tangent", "Language" -> "C++", "VectorLength" -> 6000, 
  "Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["Tangent", Real[uIn$$[dim], graduIn$$[dim, dim], duIn$$[dim
  ], gradduIn$$[dim, dim], valueOut$$[dim], gradientOut$$[dim, dim]]];

u \[RightTee] Table[SMSReal[uIn$$[i]], {i, 1, dim}];

gradU \[RightTee] Table[SMSReal[graduIn$$[i, j]], {i, 1, dim}, {j, 1, dim}];

\[Delta]u \[RightTee] Table[SMSReal[duIn$$[i]], {i, 1, dim}];

grad\[Delta]U \[RightTee] Table[SMSReal[gradduIn$$[i, j]], {i, 1, dim}, {j, 1, dim}];

energy;

tangent;
SMSExport[Ku . \[Delta]u + 
 Table[ Sum[Kgradu[[i,k,l]] grad\[Delta]U[[k,l]], {k,1, dim}, {l,1,dim }], {i,1,dim}], valueOut$$];
(* TODO: clean the mess below....*) 
SMSExport[Hu . \[Delta]u + Table[ Sum[Hgradu[[i,j,k,l]] grad\[Delta]U[[k,l]], {k,1, dim}, {l,1,dim }], {i,1,dim}, {j, 1, dim}], gradientOut$$];
SMSWrite[];
