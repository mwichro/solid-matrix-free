#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Solid Mechanics code generation*)


(* ::Input::Initialization:: *)
<<AceGen`;


(* ::Input::Initialization:: *)
dim = 3;
ScalarProduct2[U_, V_]=Tr[U . Transpose[V]];
cachedim=If[dim==3,21+6,6+3];
vectordim=If[dim==3,6,3];


(* ::Input::Initialization:: *)
initialization:=
(
\[DoubleStruckCapitalH]=gradU;
Ii\[DoubleRightTee]IdentityMatrix[dim];
\[DoubleStruckCapitalF]\[DoubleRightTee]SMSFreeze[Ii+\[DoubleStruckCapitalH]];

sq2=Sqrt[2];
If[dim==3
,
Tensor2Vector={#[[1,1]],#[[2,2]],#[[3,3]],sq2 #[[2,3]],sq2 #[[1,3]],sq2 #[[1,2]]}&;
Vector2Tensor={{#[[1]],#[[6]]/sq2,#[[5]]/sq2},{#[[6]]/sq2,#[[2]],#[[4]]/sq2},{#[[5]]/sq2,#[[4]]/sq2,#[[3]]}}&;
,
Tensor2Vector={#[[1,1]],#[[2,2]],sq2 #[[1,2]]}&;
Vector2Tensor={{#[[1]],#[[3]]/sq2},{#[[3]]/sq2,#[[2]]}}&;
];
)


(* ::Subsection:: *)
(*Energies*)


(* ::Subsubsection:: *)
(*Davydov*)


(* ::Input:: *)
(*energy:=( *)
(**)
(*\[DoubleStruckCapitalC]temp\[DoubleRightTee]Transpose[\[DoubleStruckCapitalF]] . \[DoubleStruckCapitalF];*)
(*SMSFreeze[\[DoubleStruckCapitalC]V,\[DoubleStruckCapitalC]temp//Tensor2Vector];*)
(*\[DoubleStruckCapitalC]\[DoubleRightTee](\[DoubleStruckCapitalC]V//Vector2Tensor);*)
(**)
(*J\[DoubleRightTee]SMSSqrt[Det[\[DoubleStruckCapitalC]]];*)
(*\[Psi]el\[DoubleRightTee]\[Mu]/2 (Tr[\[DoubleStruckCapitalC]]-dim-2Log[J])+\[Lambda] (Log[J])^2 ; *)
(*)*)


(* ::Subsubsection:: *)
(*neo-Hookean with volumetric-isochoric split*)


(* ::Input:: *)
(*energy:=( *)
(*\[Kappa]=\[Lambda]+2\[Mu]/3;*)
(**)
(**)
(*\[DoubleStruckCapitalC]temp\[DoubleRightTee]Transpose[\[DoubleStruckCapitalF]] . \[DoubleStruckCapitalF];*)
(*SMSFreeze[\[DoubleStruckCapitalC]V,\[DoubleStruckCapitalC]temp//Tensor2Vector];*)
(**)
(*\[DoubleStruckCapitalC]\[DoubleRightTee](\[DoubleStruckCapitalC]V//Vector2Tensor);*)
(*I3 \[RightTee]Det[\[DoubleStruckCapitalC]];*)
(*(*SMSFreeze[I3,Det[\[DoubleStruckCapitalC]]];*)*)
(*If [dim==3, rtI3\[RightTee]SMSCbrt[I3];*)
(*SMSDefineDerivative[rtI3,I3,1/3  rtI3/I3];,*)
(*rtI3\[RightTee]SMSSqrtF[I3];*)
(*SMSDefineDerivative[rtI3,I3,1/2  rtI3/I3];*)
(*];*)
(*\[DoubleStruckCapitalC]bar\[RightTee]\[DoubleStruckCapitalC]/rtI3;*)
(*IC1\[DoubleRightTee]Tr[\[DoubleStruckCapitalC]bar];*)
(**)
(*\[Psi]el\[DoubleRightTee]1/2 \[Mu](IC1-dim)+1/2 \[Kappa](1/2(I3-1)-1/2Log[I3]);*)
(*)*)


(* ::Subsection:: *)
(*Code Generation*)


(* ::Input::Initialization:: *)
residual :=  (
Ru \[DoubleRightTee] SMSD[\[Psi]el, u];
(* Rgradu \[DoubleRightTee] SMSD[\[Psi]el, Flatten[gradU]]; *)
\[DoubleStruckCapitalS]V\[DoubleRightTee]2SMSD[\[Psi]el, \[DoubleStruckCapitalC]V];
\[DoubleStruckCapitalS]\[DoubleRightTee](\[DoubleStruckCapitalS]V//Vector2Tensor);
\[DoubleStruckCapitalP]\[DoubleRightTee]\[DoubleStruckCapitalF] . \[DoubleStruckCapitalS];
Rgradu\[DoubleRightTee]Flatten[\[DoubleStruckCapitalP]];
)


(* ::Input::Initialization:: *)
tangent :=
  (
    residual;

    Ku \[DoubleRightTee] SMSD[Ru, u];
    Kgradu \[DoubleRightTee] SMSD[Ru, gradU];
    Hu \[DoubleRightTee] SMSD[Partition[Rgradu,dim], u];
    Hgradu \[DoubleRightTee] SMSD[Rgradu, Flatten[gradU]];
  )


(* ::Subsubsection:: *)
(*Residual*)


(* ::Input::Initialization:: *)
SMSInitialize["Residual", "Language" -> "C++", "VectorLength" -> 100,"Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["residual", Real[uIn$$[dim], graduIn$$[dim,dim],valueOut$$[dim], gradientOut$$[dim,dim], mu$$ , lambda$$]];

u \[RightTee] Table[SMSReal[uIn$$[i]], {i, 1, dim}];
gradU \[RightTee] Table[SMSReal[graduIn$$[i,j]], {i, 1, dim}, {j, 1, dim}];

\[Mu]\[DoubleRightTee]SMSReal[mu$$];
\[Lambda]\[DoubleRightTee]SMSReal[lambda$$];

initialization;

energy;

residual;

SMSExport[Ru, valueOut$$];
SMSExport[Partition[Rgradu,dim], gradientOut$$];
SMSWrite[];


(* ::Subsubsection:: *)
(*Tangent*)


(* ::Input::Initialization:: *)
SMSInitialize["Tangent", "Language" -> "C++", "VectorLength" -> 6000, "Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["tangent", Real[uIn$$[dim], graduIn$$[dim, dim], duIn$$[dim], gradduIn$$[dim, dim], valueOut$$[dim], gradientOut$$[dim, dim], mu$$ , lambda$$]];

u \[RightTee] Table[SMSReal[uIn$$[i]], {i, 1, dim}];

gradU \[RightTee] Table[SMSReal[graduIn$$[i, j]], {i, 1, dim}, {j, 1, dim}];

\[Delta]u \[RightTee] Table[SMSReal[duIn$$[i]], {i, 1, dim}];

grad\[Delta]U \[RightTee] Table[SMSReal[gradduIn$$[i, j]], {i, 1, dim}, {j, 1, dim}];

\[Mu]\[DoubleRightTee]SMSReal[mu$$];
\[Lambda]\[DoubleRightTee]SMSReal[lambda$$];

initialization;

energy;

tangent;

(* Hgg\[DoubleRightTee]Hgradu.Flatten[grad\[Delta]U]; *)
x\[DoubleRightTee]SMSFictive[];
Hgg\[DoubleRightTee]SMSD[Rgradu,x,"Dependency"->{Flatten[gradU],x,Flatten[grad\[Delta]U]}];

SMSExport[Ku . \[Delta]u + Table[ Sum[Kgradu[[i,k,l]] grad\[Delta]U[[k,l]], {k,1, dim}, {l,1,dim }], {i,1,dim}], valueOut$$];
SMSExport[Hu . \[Delta]u + Partition[Hgg,dim],gradientOut$$];
SMSWrite["OptimizingLoops"-> 5];


(* ::Subsubsection:: *)
(*Reference configuration Cache*)


(* ::Input::Initialization:: *)
SMSInitialize["CacheReference", "Language" -> "C++", "VectorLength" -> 5100,"Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["cache", Real[uIn$$[dim], graduIn$$[dim,dim],valueOut$$[dim], gradientOut$$[dim,dim],cacheOut$$[cachedim], mu$$ , lambda$$]];

u \[RightTee] Table[SMSReal[uIn$$[i]], {i, 1, dim}];
gradU \[RightTee] Table[SMSReal[graduIn$$[i,j]], {i, 1, dim}, {j, 1, dim}];
\[Mu]\[DoubleRightTee]SMSReal[mu$$];
\[Lambda]\[DoubleRightTee]SMSReal[lambda$$];
initialization;

energy;

residual;

\[DoubleStruckCapitalE]\[DoubleRightTee]2SMSD[\[DoubleStruckCapitalS]V,\[DoubleStruckCapitalC]V];
cache=Join[
Table[\[DoubleStruckCapitalE][[i,j]],{j,1,vectordim},{i,j,vectordim}]//Flatten,
Table[\[DoubleStruckCapitalS][[i,j]],{j,1,dim},{i,j,dim}]//Flatten
];

SMSExport[Ru, valueOut$$];
SMSExport[Partition[Rgradu,dim], gradientOut$$];
SMSExport[cache, cacheOut$$];
SMSWrite[];


(* ::Subsubsection:: *)
(*Deformed configuration Cache*)


(* ::Input::Initialization:: *)
SMSInitialize["CacheDeformed", "Language" -> "C++", "VectorLength" -> 5100,"Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["cache", Real[uIn$$[dim], graduIn$$[dim,dim],valueOut$$[dim], gradientOut$$[dim,dim],cacheOut$$[cachedim], mu$$ , lambda$$]];

u \[RightTee] Table[SMSReal[uIn$$[i]], {i, 1, dim}];
gradU \[RightTee] Table[SMSReal[graduIn$$[i,j]], {i, 1, dim}, {j, 1, dim}];
\[Mu]\[DoubleRightTee]SMSReal[mu$$];
\[Lambda]\[DoubleRightTee]SMSReal[lambda$$];
initialization;

energy;

residual;

\[DoubleStruckCapitalE]\[DoubleRightTee]2SMSD[\[DoubleStruckCapitalS]V,\[DoubleStruckCapitalC]V];
cache=Join[
Table[\[DoubleStruckCapitalE][[i,j]],{j,1,vectordim},{i,j,vectordim}]//Flatten,
Table[\[DoubleStruckCapitalS][[i,j]],{j,1,dim},{i,j,dim}]//Flatten
];

SMSExport[Ru, valueOut$$];
SMSExport[Partition[Rgradu,dim], gradientOut$$];
SMSExport[cache, cacheOut$$];
SMSWrite[];


(* ::Subsubsection:: *)
(*CachedReferenceTangent*)


(* ::Input::Initialization:: *)
SMSInitialize["CachedReferenceTangent", "Language" -> "C++", "VectorLength" -> 6000, "Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];
(* I guess we can remove valueOut$$
Should I remove  uIn$$ and  graduIn$$ *)
SMSModule["Tangent", Real[uIn$$[dim], graduIn$$[dim, dim], duIn$$[dim], gradduIn$$[dim, dim],cacheIn$$[cachedim], valueOut$$[dim], gradientOut$$[dim, dim]]];

u \[RightTee] Table[SMSReal[uIn$$[i]], {i, 1, dim}];

gradU \[RightTee] Table[SMSReal[graduIn$$[i, j]], {i, 1, dim}, {j, 1, dim}];

\[Delta]u \[RightTee] Table[SMSReal[duIn$$[i]], {i, 1, dim}];

grad\[Delta]U \[RightTee] Table[SMSReal[gradduIn$$[i, j]], {i, 1, dim}, {j, 1, dim}];

cache \[RightTee] Table[SMSReal[cacheIn$$[i]], {i, 1, cachedim}];

initialization;

Ru \[DoubleRightTee] Table[0,{dim}];
Rgradu \[DoubleRightTee] Table[0,{dim*dim}];
tangent;

If[dim==3
,
\[DoubleStruckCapitalE]\[DoubleRightTee](cache[[1;;21]]//{
{#[[1]],#[[2]],#[[3]],#[[4]],#[[5]],#[[6]]},
{#[[2]],#[[7]],#[[8]],#[[9]],#[[10]],#[[11]]},
{#[[3]],#[[8]],#[[12]],#[[13]],#[[14]],#[[15]]},
{#[[4]],#[[9]],#[[13]],#[[16]],#[[17]],#[[18]]},
{#[[5]],#[[10]],#[[14]],#[[17]],#[[19]],#[[20]]},
{#[[6]],#[[11]],#[[15]],#[[18]],#[[20]],#[[21]]}
}&);
\[DoubleStruckCapitalS]\[DoubleRightTee](cache[[22;;27]]//{{#[[1]],#[[2]],#[[3]]},{#[[2]],#[[4]],#[[5]]},{#[[3]],#[[5]],#[[6]]}}&);
,
\[DoubleStruckCapitalE]\[DoubleRightTee](cache[[1;;6]]//{{#[[1]],#[[2]],#[[3]]},{#[[2]],#[[4]],#[[5]]},{#[[3]],#[[5]],#[[6]]}}&);
\[DoubleStruckCapitalS]\[DoubleRightTee](cache[[7;;9]]//{{#[[1]],#[[2]]},{#[[2]],#[[3]]}}&);
];

\[DoubleStruckCapitalG]X2\[DoubleRightTee]grad\[Delta]U . \[DoubleStruckCapitalS];
FX\[DoubleRightTee]Transpose[\[DoubleStruckCapitalF]] . grad\[Delta]U;
FXsymV\[DoubleRightTee]((FX+Transpose[FX])/2//Tensor2Vector);
\[DoubleStruckCapitalA]V\[DoubleRightTee]\[DoubleStruckCapitalE] . FXsymV;
\[DoubleStruckCapitalA]\[DoubleRightTee](\[DoubleStruckCapitalA]V//Vector2Tensor);
\[DoubleStruckCapitalB]X\[DoubleRightTee]\[DoubleStruckCapitalF] . \[DoubleStruckCapitalA];
\[DoubleStruckCapitalG]X\[DoubleRightTee]\[DoubleStruckCapitalB]X+\[DoubleStruckCapitalG]X2;

SMSExport[Ku . \[Delta]u + Table[ Sum[Kgradu[[i,k,l]] grad\[Delta]U[[k,l]], {k,1, dim}, {l,1,dim }], {i,1,dim}], valueOut$$];
SMSExport[Hu . \[Delta]u + \[DoubleStruckCapitalG]X,gradientOut$$];
SMSWrite[];


(* ::Subsubsection:: *)
(*Cached Deformed Tangent*)


(* ::Input::Initialization:: *)
SMSInitialize["CachedDeformedTangent", "Language" -> "C++", "VectorLength" -> 6000, "Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];
(* I guess we can remove valueOut$$ *)
SMSModule["Tangent", Real[ duIn$$[dim], gradduIn$$[dim, dim],cacheIn$$[cachedim], valueOut$$[dim], gradientOut$$[dim, dim]]];

u \[RightTee] Table[SMSReal[uIn$$[i]], {i, 1, dim}];

gradU \[RightTee] Table[SMSReal[graduIn$$[i, j]], {i, 1, dim}, {j, 1, dim}];

\[Delta]u \[RightTee] Table[SMSReal[duIn$$[i]], {i, 1, dim}];

grad\[Delta]U \[RightTee] Table[SMSReal[gradduIn$$[i, j]], {i, 1, dim}, {j, 1, dim}];

cache \[RightTee] Table[SMSReal[cacheIn$$[i]], {i, 1, cachedim}];

initialization;

Ru \[DoubleRightTee] Table[0,{dim}];
Rgradu \[DoubleRightTee] Table[0,{dim*dim}];
tangent;

If[dim==3
,
\[DoubleStruckCapitalE]\[DoubleRightTee](cache[[1;;21]]//{
{#[[1]],#[[2]],#[[3]],#[[4]],#[[5]],#[[6]]},
{#[[2]],#[[7]],#[[8]],#[[9]],#[[10]],#[[11]]},
{#[[3]],#[[8]],#[[12]],#[[13]],#[[14]],#[[15]]},
{#[[4]],#[[9]],#[[13]],#[[16]],#[[17]],#[[18]]},
{#[[5]],#[[10]],#[[14]],#[[17]],#[[19]],#[[20]]},
{#[[6]],#[[11]],#[[15]],#[[18]],#[[20]],#[[21]]}
}&);
\[DoubleStruckCapitalS]\[DoubleRightTee](cache[[22;;27]]//{{#[[1]],#[[2]],#[[3]]},{#[[2]],#[[4]],#[[5]]},{#[[3]],#[[5]],#[[6]]}}&);
,
\[DoubleStruckCapitalE]\[DoubleRightTee](cache[[1;;6]]//{{#[[1]],#[[2]],#[[3]]},{#[[2]],#[[4]],#[[5]]},{#[[3]],#[[5]],#[[6]]}}&);
\[DoubleStruckCapitalS]\[DoubleRightTee](cache[[7;;9]]//{{#[[1]],#[[2]]},{#[[2]],#[[3]]}}&);
];

\[DoubleStruckCapitalG]X2\[DoubleRightTee]grad\[Delta]U . \[DoubleStruckCapitalS];
FX\[DoubleRightTee]Transpose[\[DoubleStruckCapitalF]] . grad\[Delta]U;
FXsymV\[DoubleRightTee]((FX+Transpose[FX])/2//Tensor2Vector);
\[DoubleStruckCapitalA]V\[DoubleRightTee]\[DoubleStruckCapitalE] . FXsymV;
\[DoubleStruckCapitalA]\[DoubleRightTee](\[DoubleStruckCapitalA]V//Vector2Tensor);
\[DoubleStruckCapitalB]X\[DoubleRightTee]\[DoubleStruckCapitalF] . \[DoubleStruckCapitalA];
\[DoubleStruckCapitalG]X\[DoubleRightTee]\[DoubleStruckCapitalB]X+\[DoubleStruckCapitalG]X2;

SMSExport[Ku . \[Delta]u + Table[ Sum[Kgradu[[i,k,l]] grad\[Delta]U[[k,l]], {k,1, dim}, {l,1,dim }], {i,1,dim}], valueOut$$];
SMSExport[Hu . \[Delta]u + \[DoubleStruckCapitalG]X,gradientOut$$];
SMSWrite[];
