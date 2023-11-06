#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
ClearAll["Global`*"]
K=4*^6 (*Carrying capacity*);

r=1.4 (*Growth Rate*);
  
\[Beta]1=0.7(*Cell death rate*);
\[Beta]2=0.9(*Cells infected with Mimivirus death rate*)(*Virophages\[LongDash]Known and Unknown Facts*);
\[Beta]3=0.8(*Cells infected with both Mimivirus and Virophage Death rate*)(*Virophages\[LongDash]Known and Unknown Facts*);

\[Gamma]=4.3*^-8(*Infection Rate*)(*partially derived*);
\[Delta]=2.2*^-8(*Composite formation rate*)(*partially derived*);

\[Phi]=3.2*^-2 (*Mimivirus Decay Rate*);
\[Lambda]=3.2*^-1 (*Virophage Decay Rate*);
\[Eta]=(\[Lambda]+\[Phi]) (*Composite Decay Rate*)(*\[Phi]+\[Lambda]*);
\[Psi]=0.001(*Mobilome Decay Rate*);

\[Sigma]=300 *\[Beta]2(*Mimivirus Outbreak*);
\[Omega]=1500*\[Beta]3(*Virophage Outbreak*);
\[Zeta]=30*\[Beta]3(*Composite outbreak*);


L=1+(A+Ai+Aiv+Ait+Aitv)/K(*Logistic Growth*);


(* ::Input::Initialization:: *)
ode[ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_,\[Epsilon]_,\[Kappa]_] := { (A    (r-(\[Beta]1  ( 1+(A  +Ai  +Aiv  +Ait  +Aitv  )/K))))-(\[Gamma] M    A  )-(\[Gamma]  Mpro    A  )-(\[Gamma]  Mt    A  )-(\[Gamma]  Mprot    A  )-(\[Gamma] A  MproV  )-(\[Gamma] A  MV  )-(\[Gamma] A  MVt  )-(\[Gamma] A  MtVt  )-(\[Gamma] A  MproVt  )-(\[Gamma] A  MprotV  )-(\[Gamma] A  MprotVt  ),

 (\[Gamma]  M    A  )-(\[Beta]2 ( 1+(A  +Ai  +Aiv  +Ait  +Aitv  )/K) Ai  ),(\[Eta] MV  )+(\[Eta] MVt  )+(\[Sigma] \[Epsilon] Ait   )+(\[Epsilon] \[Chi] \[Sigma]  Aitv  )+( \[Sigma]  Ai  )+( \[Chi] \[Sigma] Aiv  )+(\[Psi] Mpro  )+(\[Psi] Mt  )-(\[Gamma]  M    A  )-(\[Delta]   M    V  )-(\[Delta]   M  Vt  )-(\[Phi]  M  ),(\[Eta] MprotV  )+(\[Eta] MV  )+(\[Eta] MproV  )+(\[Eta] MtV  )+(\[Kappa] \[Omega]  Aitv  )+( \[Omega]  Aiv  )+(\[Psi] Vt  )-(\[Delta]   M    V  )-(\[Delta]  Mt    V  )-(\[Delta]  V    Mprot  )-(\[Delta]  V    Mpro  )-(\[Lambda]  V  ),
(\[Gamma]  Mpro    A  )-(\[Beta]3 ( 1+(A  +Ai  +Aiv  +Ait  +Aitv  )/K) Aiv  )+(\[Gamma] A  MproV  )+(\[Gamma] A  MV  ),

 (\[Eta] MproV )+(\[Eta] MproVt )+(\[Epsilon] \[Chi] \[Sigma]  Aitv )+(\[Chi] \[Omega]  Aiv )+(\[Psi] Mprot )-(\[Phi]  Mpro )-(\[Gamma]  Mpro   A )-(\[Psi] Mpro )-(\[Delta]    Vt   Mpro )-(\[Delta]  V   Mpro ),

 (\[Eta] MtVt )+(\[Eta] MtV )+( \[Sigma] \[Epsilon] Ait  )+( \[Epsilon] \[Chi] \[Sigma]  Aitv )+(\[Psi] Mprot )-(\[Gamma]  Mt   A )-(\[Delta]  Mt   Vt )-(\[Delta]    V   Mt )-(\[Phi]  Mt )-(\[Psi] Mt ),

 (\[Eta] MprotVt )+(\[Eta] MproVt )+(\[Eta] MtVt )+(\[Eta] MVt )+( \[Kappa] \[Omega]  Aitv )-(\[Delta]  M   Vt )-(\[Delta]  Mt   Vt )-(\[Delta]  Vt   Mpro )-(\[Delta]  Vt   Mprot )-(\[Lambda]  Vt )-(\[Psi] Vt ),

(\[Gamma]  Mt    A  )-(\[Beta]2 ( 1+(A  +Ai  +Aiv  +Ait  +Aitv  )/K) Ait  ),
 (\[Gamma]  Mprot    A  )-(\[Beta]3 ( 1+(A  +Ai  +Aiv  +Ait  +Aitv  )/K) Aitv  )+(\[Gamma] A  MVt  )+(\[Gamma] A  MtVt  )+(\[Gamma] A  MproVt  )+(\[Gamma] A  MprotV  )+(\[Gamma] A  MprotVt  )+(\[Gamma] A  MtV  ),

(\[Eta] MprotVt )+(\[Eta] MprotV )+( \[Epsilon] \[Sigma] \[Chi]  Aitv )-(\[Psi] Mprot )-(\[Psi] Mprot )-(\[Phi]  Mprot )-(\[Gamma]  Mprot   A )-(\[Delta]  V   Mprot )-(\[Delta]  Vt   Mprot ),

  (\[Delta]   V   Mt )-(\[Gamma] A MtV )-(\[Eta] MtV ),

  (\[Delta]   Vt   M )-(\[Gamma] A MVt )-(\[Eta] MVt ),
  (\[Delta]  Vt   Mt )-(\[Gamma] A MtVt )-(\[Eta] MtVt ),

(\[Delta]  Vt   Mpro )-(\[Gamma] A MproVt )-(\[Eta] MproVt ),

  (\[Delta]  V   Mpro )-(\[Gamma] A MproV )-(\[Eta] MproV ),
 (\[Delta]    V   M )-(\[Gamma] A MV )-(\[Eta] MV ),

 (\[Delta]  V   Mprot )-(\[Gamma] A MprotV )-(\[Eta] MprotV ),

(\[Delta]  Vt   Mprot )-(\[Gamma] A MprotVt )-(\[Eta] MprotVt )}



(* ::Text::Initialization:: *)
(*Bifurcation Analysis*)


(* ::Input::Initialization:: *)
findSteadyState[\[Chi]_,\[Epsilon]_,\[Kappa]_]:=NSolve[ode[ A,M,V,Ai,Aiv,Mpro,Mt,Vt,Ait,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,\[Chi],\[Epsilon],\[Kappa]] ==0 ,{A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt},Reals];


resultsTable=Table[solutions=findSteadyState[\[Chi]val,\[Epsilon]val,\[Kappa]val];Append[{\[Chi]val,\[Epsilon]val,\[Kappa]val},#]&/@solutions,{\[Chi]val,0,1,0.1},{\[Epsilon]val,0,1,0.1},{\[Kappa]val,0,1,0.1}]
flattenedTable=Flatten[resultsTable,1];
bi2=TableForm[flattenedTable,TableHeadings->None]
Export["bi2.csv",bi2]


(* ::Input::Initialization:: *)
 databi =Import["/home/kai-hua/OneDrive/Thesis/bi2.csv"]


(* ::Text::Initialization:: *)
(*Bifurcation Analysis Data processing*)


(* ::Input::Initialization:: *)
processedDatabi=ToExpression[databi];
processedDatabi2=Flatten[Flatten[processedDatabi,1],1];


(* ::Input::Initialization:: *)
filteredLists=Select[processedDatabi2,FreeQ[#[[4]],Rule[_,_?Negative]]&];


(* ::Input::Initialization:: *)
 fixedpointsbi=DeleteDuplicatesBy[filteredLists,ToString];


(* ::Input::Initialization:: *)
Export["fpbi.csv",fixedpointsbi]


(* ::Input::Initialization:: *)
fixedpointsbi=ToExpression[Import["/home/kai-hua/fpbi.csv"]]
