#!/ usr / bin / env wolframscript
(*::Package::*)

    (*Initilization and input *)
  << AceGen`;
dim = 3;

(**)
(* Using "acegen_scratch" allows ReplaceAceGenScratch.sh to find and replace all temporary variables *)
(**)

SMSInitialize["neoHooke_" < > ToString[dim]<> "D_Tangent",
                "Language"->"C++",
                "GlobalNames"->{"acegen_scratch", "i", "b"}];
SMSModule
  ["tangent",
   Real
     [gradU$$[dim, dim], gradDU$$[dim, dim], out$$[dim, dim], mu$$, lambda$$]];
gradU \[RightTee] Table[SMSReal[gradU$$[i, j]], {i, 1, dim}, {j, 1, dim}];
grad\[Delta] U \[RightTee] Table[SMSReal[gradDU$$[i, j]],
                                 {i, 1, dim},
                                 {j, 1, dim}];
\[Mu] \[RightTee] SMSReal[mu$$];
\[Lambda] \[RightTee] SMSReal[lambda$$];
\[DoubleStruckCapitalH] \[DoubleRightTee] ArrayPad[gradU, {0, 3 - dim}];


(*Strain energy function *)
\[DoubleStruckCapitalF] \[DoubleRightTee] IdentityMatrix[3] +
  [DoubleStruckCapitalH];
\[DoubleStruckCapitalC] \[DoubleRightTee] Transpose[\[DoubleStruckCapitalF]]
  . \[DoubleStruckCapitalF];
J2 \[DoubleRightTee] Det[\[DoubleStruckCapitalC]];
\[Psi] el \[DoubleRightTee] \[Mu] /
      2(Tr[\[DoubleStruckCapitalC]] - 3 - Log[J2]) +\[Lambda](Log[J2] / 2) ^
  2;


(*Tangent *)
\[DoubleStruckCapitalP] \[DoubleRightTee] SMSD
  [\[Psi] el, \[DoubleStruckCapitalH], "IgnoreNumbers"->True];
\[Xi] \[RightTee] SMSFictive[];
Kgrad\[Delta] U \[DoubleRightTee] SMSD[\[DoubleStruckCapitalP], \[Xi],
                                       "Dependency"->{Flatten[gradU],\[Xi],
                                                      Flatten[grad\[Delta] U]}];


(*Output *)SMSExport[Kgrad\[Delta] U[[1;; dim, 1;; dim]], out$$];
SMSWrite["OptimizingLoops"->2];
