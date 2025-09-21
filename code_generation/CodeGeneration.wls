#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Hyperelastic material - code generation*)


(* ::Text:: *)
(*Before running define model and dim*)


(* ::Input::Initialization:: *)
model={"neoHooke","neoHookeSplit","Hencky1","Hencky2","Hencky3"}[[1]];
dim = {2,3}[[2]];


(* ::Subsection::Closed:: *)
(*Initialization*)


(* ::Input::Initialization:: *)
<<AceGen`;


(* ::Input::Initialization:: *)
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


(* ::Input::Initialization:: *)
residual :=  (
\[DoubleStruckCapitalS]V\[DoubleRightTee]2SMSD[\[CapitalPsi], \[DoubleStruckCapitalC]V];
\[DoubleStruckCapitalS]\[DoubleRightTee](\[DoubleStruckCapitalS]V//Vector2Tensor);
\[DoubleStruckCapitalP]\[DoubleRightTee]\[DoubleStruckCapitalF] . \[DoubleStruckCapitalS];
)


(* ::Input::Initialization:: *)
(* Kelvin notation *)
If[dim==3
,
i1={1,0,0};
i2={0,1,0};
i3={0,0,1};
v1=Outer[Times,i1,i1];
v2=Outer[Times,i2,i2];
v3=Outer[Times,i3,i3];
v4=1/Sqrt[2] (Outer[Times,i2,i3]+Outer[Times,i3,i2]);
v5=1/Sqrt[2] (Outer[Times,i1,i3]+Outer[Times,i3,i1]);
v6=1/Sqrt[2] (Outer[Times,i2,i1]+Outer[Times,i1,i2]);
V={v1,v2,v3,v4,v5,v6}
,
i1={1,0};
i2={0,1};
v1=Outer[Times,i1,i1];
v2=Outer[Times,i2,i2];
v3=1/Sqrt[2] (Outer[Times,i2,i1]+Outer[Times,i1,i2]);
V={v1,v2,v3}
];


(* ::Subsection::Closed:: *)
(*Energy*)


(* ::Input::Initialization:: *)
(* neoHookean as in Davydov - but with corrected \[Lambda]-term *)
If[ model=="neoHooke",
energy:=(
\[DoubleStruckCapitalC]temp\[DoubleRightTee]Transpose[\[DoubleStruckCapitalF]] . \[DoubleStruckCapitalF];
SMSFreeze[\[DoubleStruckCapitalC]V,\[DoubleStruckCapitalC]temp//Tensor2Vector];
\[DoubleStruckCapitalC]\[DoubleRightTee](\[DoubleStruckCapitalC]V//Vector2Tensor);
J2\[DoubleRightTee]Det[\[DoubleStruckCapitalC]];
\[CapitalPsi]\[DoubleRightTee]\[Mu]/2 (Tr[\[DoubleStruckCapitalC]]-dim-Log[J2])+\[Lambda]/2 (Log[J2]/2)^2 
);
];


(* ::Input::Initialization:: *)
(* neoHookean with volumetric-isochoric split *)
If[ model=="neoHookeSplit",
energy:=(
\[Kappa]=\[Lambda]+2\[Mu]/3;
\[DoubleStruckCapitalC]temp\[DoubleRightTee]Transpose[\[DoubleStruckCapitalF]] . \[DoubleStruckCapitalF];
SMSFreeze[\[DoubleStruckCapitalC]V,\[DoubleStruckCapitalC]temp//Tensor2Vector];
\[DoubleStruckCapitalC]\[DoubleRightTee](\[DoubleStruckCapitalC]V//Vector2Tensor);
J2 \[RightTee]Det[\[DoubleStruckCapitalC]];
I3 \[DoubleRightTee]J2;
If [dim==3, 
rtI3\[RightTee]SMSCbrt[I3];
SMSDefineDerivative[rtI3,I3,1/3  rtI3/I3]
,
rtI3\[RightTee]SMSSqrtF[I3];
SMSDefineDerivative[rtI3,I3,1/2  rtI3/I3]
];
\[DoubleStruckCapitalC]bar\[RightTee]\[DoubleStruckCapitalC]/rtI3;
IC1\[DoubleRightTee]Tr[\[DoubleStruckCapitalC]bar];
\[CapitalPsi]\[DoubleRightTee]1/2 \[Mu](IC1-dim)+1/2 \[Kappa](1/2(I3-1)-1/2Log[I3])
);
];


(* ::Input::Initialization:: *)
(* Hencky model with Pade approximant order 1 *)
If[ model=="Hencky1",
energy:=(
\[DoubleStruckCapitalC]temp\[DoubleRightTee]Transpose[\[DoubleStruckCapitalF]] . \[DoubleStruckCapitalF];
SMSFreeze[\[DoubleStruckCapitalC]V,\[DoubleStruckCapitalC]temp//Tensor2Vector];
\[DoubleStruckCapitalC]\[DoubleRightTee](\[DoubleStruckCapitalC]V//Vector2Tensor);
J2\[DoubleRightTee]Det[\[DoubleStruckCapitalC]];
(* (m,n)=(1,1) *)
\[DoubleStruckCapitalP]m\[DoubleRightTee]2(-Ii+\[DoubleStruckCapitalC]);
\[DoubleStruckCapitalQ]m\[DoubleRightTee]Ii+\[DoubleStruckCapitalC];
\[DoubleStruckCapitalQ]mi\[DoubleRightTee]SMSInverse[\[DoubleStruckCapitalQ]m];
\[DoubleStruckCapitalH]e\[DoubleRightTee]\[DoubleStruckCapitalP]m . \[DoubleStruckCapitalQ]mi/2;
\[CapitalPsi]\[DoubleRightTee]1/2\[Lambda] (Tr[\[DoubleStruckCapitalH]e])^2+\[Mu] Tr[\[DoubleStruckCapitalH]e . \[DoubleStruckCapitalH]e]
);
];


(* ::Input::Initialization:: *)
(* Hencky model with Pade approximant order 2 *)
If[ model=="Hencky2",
energy:=(
\[DoubleStruckCapitalC]temp\[DoubleRightTee]Transpose[\[DoubleStruckCapitalF]] . \[DoubleStruckCapitalF];
SMSFreeze[\[DoubleStruckCapitalC]V,\[DoubleStruckCapitalC]temp//Tensor2Vector];
\[DoubleStruckCapitalC]\[DoubleRightTee](\[DoubleStruckCapitalC]V//Vector2Tensor);
\[DoubleStruckCapitalC]2\[DoubleRightTee]\[DoubleStruckCapitalC] . \[DoubleStruckCapitalC];
J2\[DoubleRightTee]Det[\[DoubleStruckCapitalC]];
(* (m,n)=(1,1) *)
\[DoubleStruckCapitalP]m\[DoubleRightTee]3(-Ii+\[DoubleStruckCapitalC]2);
\[DoubleStruckCapitalQ]m\[DoubleRightTee]Ii+4\[DoubleStruckCapitalC]+\[DoubleStruckCapitalC]2;
\[DoubleStruckCapitalQ]mi\[DoubleRightTee]SMSInverse[\[DoubleStruckCapitalQ]m];
\[DoubleStruckCapitalH]e\[DoubleRightTee]\[DoubleStruckCapitalP]m . \[DoubleStruckCapitalQ]mi/2;
\[CapitalPsi]\[DoubleRightTee]1/2\[Lambda] (Tr[\[DoubleStruckCapitalH]e])^2+\[Mu] Tr[\[DoubleStruckCapitalH]e . \[DoubleStruckCapitalH]e]
);
];


(* ::Input::Initialization:: *)
(* Hencky model with Pade approximant order 3 *)
If[ model=="Hencky3",
energy:=(
\[DoubleStruckCapitalC]temp\[DoubleRightTee]Transpose[\[DoubleStruckCapitalF]] . \[DoubleStruckCapitalF];
SMSFreeze[\[DoubleStruckCapitalC]V,\[DoubleStruckCapitalC]temp//Tensor2Vector];
\[DoubleStruckCapitalC]\[DoubleRightTee](\[DoubleStruckCapitalC]V//Vector2Tensor);
\[DoubleStruckCapitalC]2\[DoubleRightTee]\[DoubleStruckCapitalC] . \[DoubleStruckCapitalC];
\[DoubleStruckCapitalC]3\[DoubleRightTee]\[DoubleStruckCapitalC]2 . \[DoubleStruckCapitalC];
J2\[DoubleRightTee]Det[\[DoubleStruckCapitalC]];
(* (m,n)=(1,1) *)
\[DoubleStruckCapitalP]m\[DoubleRightTee]-11Ii-27\[DoubleStruckCapitalC]+27\[DoubleStruckCapitalC]2+11\[DoubleStruckCapitalC]3;
\[DoubleStruckCapitalQ]m\[DoubleRightTee]3(Ii+9\[DoubleStruckCapitalC]+9\[DoubleStruckCapitalC]2+\[DoubleStruckCapitalC]3);
\[DoubleStruckCapitalQ]mi\[DoubleRightTee]SMSInverse[\[DoubleStruckCapitalQ]m];
\[DoubleStruckCapitalH]e\[DoubleRightTee]\[DoubleStruckCapitalP]m . \[DoubleStruckCapitalQ]mi/2;
\[CapitalPsi]\[DoubleRightTee]1/2\[Lambda] (Tr[\[DoubleStruckCapitalH]e])^2+\[Mu] Tr[\[DoubleStruckCapitalH]e . \[DoubleStruckCapitalH]e]
);
];


(* ::Subsection:: *)
(*Code Generation*)


(* ::Subsubsection::Closed:: *)
(*Residual*)


(* ::Input::Initialization:: *)
SMSInitialize[model<>"_"<>ToString[dim]<>"D_Residual", "Language" -> "C++", "VectorLength" -> 100,"Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["residual", Real[ gradU$$[dim,dim], out$$[dim,dim], mu$$, lambda$$]];

gradU \[RightTee] Table[SMSReal[gradU$$[i,j]], {i, 1, dim}, {j, 1, dim}];
\[Mu]\[DoubleRightTee]SMSReal[mu$$];
\[Lambda]\[DoubleRightTee]SMSReal[lambda$$];

initialization;
energy;
residual;

SMSExport[ \[DoubleStruckCapitalP], out$$];
SMSWrite[];


(* ::Subsubsection::Closed:: *)
(*Tangent*)


(* ::Input::Initialization:: *)
SMSInitialize[model<>"_"<>ToString[dim]<>"D_Tangent", "Language" -> "C++", "VectorLength" -> 6000, "Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["tangent", Real[gradU$$[dim, dim], gradDU$$[dim, dim], out$$[dim, dim], mu$$ , lambda$$]];

gradU \[RightTee] Table[SMSReal[gradU$$[i, j]], {i, 1, dim}, {j, 1, dim}];
grad\[CapitalDelta]U \[RightTee] Table[SMSReal[gradDU$$[i, j]], {i, 1, dim}, {j, 1, dim}];
\[Mu]\[DoubleRightTee]SMSReal[mu$$];
\[Lambda]\[DoubleRightTee]SMSReal[lambda$$];

initialization;
energy;
residual;

\[Xi]\[DoubleRightTee]SMSFictive[];
\[DoubleStruckCapitalG]\[DoubleRightTee]SMSD[ Flatten[\[DoubleStruckCapitalP]],\[Xi],"Dependency"->{Flatten[gradU],\[Xi],Flatten[grad\[CapitalDelta]U]}];

SMSExport[ Partition[\[DoubleStruckCapitalG],dim],out$$ ];
SMSWrite["OptimizingLoops"-> 5];


(* ::Subsubsection::Closed:: *)
(*Residual - cache - reference configuration*)


(* ::Input::Initialization:: *)
SMSInitialize[model<>"_"<>ToString[dim]<>"D_CacheReference_Residual", "Language" -> "C++", "VectorLength" -> 5100,"Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["cache", Real[gradU$$[dim,dim], out$$[dim,dim],cacheOut$$[cachedim], mu$$ , lambda$$]];

gradU \[RightTee] Table[SMSReal[gradU$$[i,j]], {i, 1, dim}, {j, 1, dim}];
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

SMSExport[ \[DoubleStruckCapitalP], out$$ ];
SMSExport[ cache, cacheOut$$ ];
SMSWrite[];


(* ::Subsubsection::Closed:: *)
(*Tangent - cache - reference configuration*)


(* ::Input::Initialization:: *)
SMSInitialize[model<>"_"<>ToString[dim]<>"D_CacheReference_Tangent", "Language" -> "C++", "VectorLength" -> 6000, "Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["Tangent", Real[gradU$$[dim, dim], gradDU$$[dim, dim],cacheIn$$[cachedim], out$$[dim, dim]]];

gradU \[RightTee] Table[SMSReal[gradU$$[i, j]], {i, 1, dim}, {j, 1, dim}];
grad\[CapitalDelta]U \[RightTee] Table[SMSReal[gradDU$$[i, j]], {i, 1, dim}, {j, 1, dim}];
cache \[RightTee] Table[SMSReal[cacheIn$$[i]], {i, 1, cachedim}];

initialization;

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

\[DoubleStruckCapitalG]2\[DoubleRightTee]grad\[CapitalDelta]U . \[DoubleStruckCapitalS];
FX\[DoubleRightTee]Transpose[\[DoubleStruckCapitalF]] . grad\[CapitalDelta]U;
FXsymV\[DoubleRightTee]((FX+Transpose[FX])/2//Tensor2Vector);
\[DoubleStruckCapitalA]V\[DoubleRightTee]\[DoubleStruckCapitalE] . FXsymV;
\[DoubleStruckCapitalA]\[DoubleRightTee](\[DoubleStruckCapitalA]V//Vector2Tensor);
\[DoubleStruckCapitalB]X\[DoubleRightTee]\[DoubleStruckCapitalF] . \[DoubleStruckCapitalA];
\[DoubleStruckCapitalG]\[DoubleRightTee]\[DoubleStruckCapitalB]X+\[DoubleStruckCapitalG]2;

SMSExport[ \[DoubleStruckCapitalG],out$$ ];
SMSWrite[];


(* ::Subsubsection::Closed:: *)
(*Residual - cache - deformed configuration*)


(* ::Text:: *)
(*Attention: here integration is to be performed over the deformed configuration*)


(* ::Input::Initialization:: *)
SMSInitialize[model<>"_"<>ToString[dim]<>"D_CacheDeformed_Residual", "Language" -> "C++", "VectorLength" -> 5100,"Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["cache", Real[gradU$$[dim,dim], out$$[dim,dim],cacheOut$$[cachedim], mu$$ , lambda$$]];

(* here, gradU is the gradient in the REFERENCE configuration *)
gradU \[RightTee] Table[SMSReal[gradU$$[i,j]], {i, 1, dim}, {j, 1, dim}];
\[Mu]\[DoubleRightTee]SMSReal[mu$$];
\[Lambda]\[DoubleRightTee]SMSReal[lambda$$];

initialization;
energy;
residual;

J\[DoubleRightTee]SMSSqrt[J2];

\[DoubleStruckCapitalE]kelv\[DoubleRightTee]2SMSD[\[DoubleStruckCapitalS]V,\[DoubleStruckCapitalC]V];
\[DoubleStruckCapitalE]tens=Table[Sum[V[[A,i,j]]*\[DoubleStruckCapitalE]kelv[[A,B]]*V[[B,k,l]],{A,1,vectordim},{B,1,vectordim}],{i,1,dim},{j,1,dim},{k,1,dim},{l,1,dim}];
J\[DoubleStruckCapitalC]tens=Table[Sum[\[DoubleStruckCapitalF][[i,I]]*\[DoubleStruckCapitalF][[j,J]]*\[DoubleStruckCapitalF][[k,K]]*\[DoubleStruckCapitalF][[l,L]]*\[DoubleStruckCapitalE]tens[[I,J,K,L]],{I,1,dim},{J,1,dim},{K,1,dim},{L,1,dim}]
,{i,1,dim},{j,1,dim},{k,1,dim},{l,1,dim}];
\[DoubleStruckCapitalC]tens=J\[DoubleStruckCapitalC]tens/J;
\[DoubleStruckCapitalC]kelv\[DoubleRightTee]Table[Sum[V[[A,i,j]]*\[DoubleStruckCapitalC]tens[[i,j,k,l]]*V[[B,k,l]],{i,1,dim},{j,1,dim},{k,1,dim},{l,1,dim}],{A,1,vectordim},{B,1,vectordim}];

\[Tau]\[DoubleRightTee]\[DoubleStruckCapitalP] . Transpose[\[DoubleStruckCapitalF]];
\[Sigma]\[DoubleRightTee]\[Tau]/J;

cache=Join[
Table[\[DoubleStruckCapitalC]kelv[[i,j]],{j,1,vectordim},{i,j,vectordim}]//Flatten,
Table[\[Sigma][[i,j]],{j,1,dim},{i,j,dim}]//Flatten
];

SMSExport[ \[Sigma], out$$ ];
SMSExport[ cache, cacheOut$$ ];
SMSWrite[];


(* ::Subsubsection::Closed:: *)
(*Tangent - cache - deformed configuration*)


(* ::Text:: *)
(*Attention: here integration necessarily is to be performed over the deformed configuration!*)


(* ::Input::Initialization:: *)
SMSInitialize[model<>"_"<>ToString[dim]<>"D_CacheDeformed_Tangent", "Language" -> "C++", "VectorLength" -> 6000, "Mode" -> "Optimal", "GlobalNames"-> {"acegen_scratch","i","b"}];

SMSModule["Tangent", Real[gradDU$$[dim, dim],cacheIn$$[cachedim], out$$[dim, dim]]];

(* here, grad\[CapitalDelta]U is the gradient in the CURRENT configuration *)
grad\[CapitalDelta]U \[RightTee] Table[SMSReal[gradDU$$[i, j]], {i, 1, dim}, {j, 1, dim}];
cache \[RightTee] Table[SMSReal[cacheIn$$[i]], {i, 1, cachedim}];

(* initiated for compatibility with other formulations - no effect *)
gradU=Table[0, {i, 1, dim}, {j, 1, dim}];

initialization;

If[dim==3
,
\[DoubleStruckCapitalC]kelv\[DoubleRightTee](cache[[1;;21]]//{
{#[[1]],#[[2]],#[[3]],#[[4]],#[[5]],#[[6]]},
{#[[2]],#[[7]],#[[8]],#[[9]],#[[10]],#[[11]]},
{#[[3]],#[[8]],#[[12]],#[[13]],#[[14]],#[[15]]},
{#[[4]],#[[9]],#[[13]],#[[16]],#[[17]],#[[18]]},
{#[[5]],#[[10]],#[[14]],#[[17]],#[[19]],#[[20]]},
{#[[6]],#[[11]],#[[15]],#[[18]],#[[20]],#[[21]]}
}&);
\[Sigma]\[DoubleRightTee](cache[[22;;27]]//{{#[[1]],#[[2]],#[[3]]},{#[[2]],#[[4]],#[[5]]},{#[[3]],#[[5]],#[[6]]}}&);
,
\[DoubleStruckCapitalC]kelv\[DoubleRightTee](cache[[1;;6]]//{{#[[1]],#[[2]],#[[3]]},{#[[2]],#[[4]],#[[5]]},{#[[3]],#[[5]],#[[6]]}}&);
\[Sigma]\[DoubleRightTee](cache[[7;;9]]//{{#[[1]],#[[2]]},{#[[2]],#[[3]]}}&);
];

Gg\[DoubleRightTee]grad\[CapitalDelta]U . \[Sigma];
gsV\[DoubleRightTee]((grad\[CapitalDelta]U+Transpose[grad\[CapitalDelta]U])/2//Tensor2Vector);
CgsV\[DoubleRightTee]\[DoubleStruckCapitalC]kelv . gsV;
Cgs\[DoubleRightTee](CgsV//Vector2Tensor);

SMSExport[ Gg+Cgs,out$$ ];
SMSWrite[];
