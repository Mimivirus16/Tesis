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


(* ::Text::Initialization:: *)
(*\[Chi] sensitivity analysis*)


(* ::Input::Initialization:: *)
eqA[ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_] :=  (A    (r-(\[Beta]1  L)))-(\[Gamma] M    A  )-(\[Gamma]  Mpro    A  )-(\[Gamma]  Mt    A  )-(\[Gamma]  Mprot    A  )-(\[Gamma] A  MproV  )-(\[Gamma] A  MV  )-(\[Gamma] A  MVt  )-(\[Gamma] A  MtVt  )-(\[Gamma] A  MproVt  )-(\[Gamma] A  MprotV  )-(\[Gamma] A  MprotVt  );

eqAi [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_] := (\[Gamma]  M    A  )-(\[Beta]2 L Ai  );

eqM [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_] := (\[Eta] MV  )+(\[Eta] MVt  )+(\[Sigma] \[Epsilon] Ait   )+(\[Epsilon] \[Chi] \[Sigma]  Aitv  )+( \[Sigma]  Ai  )+( \[Chi] \[Sigma] Aiv  )+(\[Psi] Mpro  )+(\[Psi] Mt  )-(\[Gamma]  M    A  )-(\[Delta]   M    V  )-(\[Delta]   M  Vt  )-(\[Phi]  M  );

eqV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_] := (\[Eta] MprotV  )+(\[Eta] MV  )+(\[Eta] MproV  )+(\[Eta] MtV  )+(\[Kappa] \[Omega]  Aitv  )+( \[Omega]  Aiv  )+(\[Psi] Vt  )-(\[Delta]   M    V  )-(\[Delta]  Mt    V  )-(\[Delta]  V    Mprot  )-(\[Delta]  V    Mpro  )-(\[Lambda]  V  );

eqAiv [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_] := (\[Gamma]  Mpro    A  )-(\[Beta]3 L Aiv  )+(\[Gamma] A  MproV  )+(\[Gamma] A  MV  );

eqMpro [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_] := (\[Eta] MproV )+(\[Eta] MproVt )+(\[Epsilon] \[Chi] \[Sigma]  Aitv )+(\[Chi] \[Omega]  Aiv )+(\[Psi] Mprot )-(\[Phi]  Mpro )-(\[Gamma]  Mpro   A )-(\[Psi] Mpro )-(\[Delta]    Vt   Mpro )-(\[Delta]  V   Mpro );

eqMt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_] := (\[Eta] MtVt )+(\[Eta] MtV )+( \[Sigma] \[Epsilon] Ait  )+( \[Epsilon] \[Chi] \[Sigma]  Aitv )+(\[Psi] Mprot )-(\[Gamma]  Mt   A )-(\[Delta]  Mt   Vt )-(\[Delta]    V   Mt )-(\[Phi]  Mt )-(\[Psi] Mt );

eqVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_] := (\[Eta] MprotVt )+(\[Eta] MproVt )+(\[Eta] MtVt )+(\[Eta] MVt )+( \[Kappa] \[Omega]  Aitv )-(\[Delta]  M   Vt )-(\[Delta]  Mt   Vt )-(\[Delta]  Vt   Mpro )-(\[Delta]  Vt   Mprot )-(\[Lambda]  Vt )-(\[Psi] Vt );

eqAit  [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_] :=(\[Gamma]  Mt    A  )-(\[Beta]2 L Ait  );

eqAitv [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_] := (\[Gamma]  Mprot    A  )-(\[Beta]3 L Aitv  )+(\[Gamma] A  MVt  )+(\[Gamma] A  MtVt  )+(\[Gamma] A  MproVt  )+(\[Gamma] A  MprotV  )+(\[Gamma] A  MprotVt  )+(\[Gamma] A  MtV  );

eqMprot [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_]:=(\[Eta] MprotVt )+(\[Eta] MprotV )+( \[Epsilon] \[Sigma] \[Chi]  Aitv )-(\[Psi] Mprot )-(\[Psi] Mprot )-(\[Phi]  Mprot )-(\[Gamma]  Mprot   A )-(\[Delta]  V   Mprot )-(\[Delta]  Vt   Mprot );

eqMtV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_]:=  (\[Delta]   V   Mt )-(\[Gamma] A MtV )-(\[Eta] MtV );

eqMVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_]:=  (\[Delta]   Vt   M )-(\[Gamma] A MVt )-(\[Eta] MVt );

eqMtVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_]:=  (\[Delta]  Vt   Mt )-(\[Gamma] A MtVt )-(\[Eta] MtVt );

eqMproVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_]:=  (\[Delta]  Vt   Mpro )-(\[Gamma] A MproVt )-(\[Eta] MproVt );

eqMproV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_]:=  (\[Delta]  V   Mpro )-(\[Gamma] A MproV )-(\[Eta] MproV );

eqMV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_]:= (\[Delta]    V   M )-(\[Gamma] A MV )-(\[Eta] MV );

eqMprotV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_]:=  (\[Delta]  V   Mprot )-(\[Gamma] A MprotV )-(\[Eta] MprotV );

eqMprotVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Chi]_]:=  (\[Delta]  Vt   Mprot )-(\[Gamma] A MprotVt )-(\[Eta] MprotVt );



(* ::Input::Initialization:: *)
sensitivity1={};
parametervalue={2.997,2.991,2.985,0.303,0.309,0.315};
For[\[Chi]=1,\[Chi]<=Length[parametervalue],\[Chi]++,(*Solve for the fixed point*)fp=NSolve[{eqA[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]] == 0, eqM[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0, eqV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0, eqAi[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0, eqAiv[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0, eqMpro[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0, eqMt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0, eqVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0, eqAit[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0, eqAitv[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0, eqMprot[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0,
eqMtV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0,
eqMVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0,
eqMtVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0,
eqMproVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0,eqMproV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0,
eqMV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0,eqMprotV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0,
eqMprotVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Chi]]]]==0},{A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt},Reals];AppendTo[sensitivity1,fp]];
Export["sensitivity1",sensitivity1]


(* ::Text::Initialization:: *)
(*\[Epsilon] sensitivity analysis*)


(* ::Input::Initialization:: *)
eqA[ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_] :=  (A    (r-(\[Beta]1  L)))-(\[Gamma] M    A  )-(\[Gamma]  Mpro    A  )-(\[Gamma]  Mt    A  )-(\[Gamma]  Mprot    A  )-(\[Gamma] A  MproV  )-(\[Gamma] A  MV  )-(\[Gamma] A  MVt  )-(\[Gamma] A  MtVt  )-(\[Gamma] A  MproVt  )-(\[Gamma] A  MprotV  )-(\[Gamma] A  MprotVt  );

eqAi [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_] := (\[Gamma]  M    A  )-(\[Beta]2 L Ai  );

eqM [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_] := (\[Eta] MV  )+(\[Eta] MVt  )+(\[Sigma] \[Epsilon] Ait   )+(\[Epsilon] \[Chi] \[Sigma]  Aitv  )+( \[Sigma]  Ai  )+( \[Chi] \[Sigma] Aiv  )+(\[Psi] Mpro  )+(\[Psi] Mt  )-(\[Gamma]  M    A  )-(\[Delta]   M    V  )-(\[Delta]   M  Vt  )-(\[Phi]  M  );

eqV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_] := (\[Eta] MprotV  )+(\[Eta] MV  )+(\[Eta] MproV  )+(\[Eta] MtV  )+(\[Kappa] \[Omega]  Aitv  )+( \[Omega]  Aiv  )+(\[Psi] Vt  )-(\[Delta]   M    V  )-(\[Delta]  Mt    V  )-(\[Delta]  V    Mprot  )-(\[Delta]  V    Mpro  )-(\[Lambda]  V  );

eqAiv [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_] := (\[Gamma]  Mpro    A  )-(\[Beta]3 L Aiv  )+(\[Gamma] A  MproV  )+(\[Gamma] A  MV  );

eqMpro [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_] := (\[Eta] MproV )+(\[Eta] MproVt )+(\[Epsilon] \[Chi] \[Sigma]  Aitv )+(\[Chi] \[Omega]  Aiv )+(\[Psi] Mprot )-(\[Phi]  Mpro )-(\[Gamma]  Mpro   A )-(\[Psi] Mpro )-(\[Delta]    Vt   Mpro )-(\[Delta]  V   Mpro );

eqMt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_] := (\[Eta] MtVt )+(\[Eta] MtV )+( \[Sigma] \[Epsilon] Ait  )+( \[Epsilon] \[Chi] \[Sigma]  Aitv )+(\[Psi] Mprot )-(\[Gamma]  Mt   A )-(\[Delta]  Mt   Vt )-(\[Delta]    V   Mt )-(\[Phi]  Mt )-(\[Psi] Mt );

eqVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_] := (\[Eta] MprotVt )+(\[Eta] MproVt )+(\[Eta] MtVt )+(\[Eta] MVt )+( \[Kappa] \[Omega]  Aitv )-(\[Delta]  M   Vt )-(\[Delta]  Mt   Vt )-(\[Delta]  Vt   Mpro )-(\[Delta]  Vt   Mprot )-(\[Lambda]  Vt )-(\[Psi] Vt );

eqAit  [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_] :=(\[Gamma]  Mt    A  )-(\[Beta]2 L Ait  );

eqAitv [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_] := (\[Gamma]  Mprot    A  )-(\[Beta]3 L Aitv  )+(\[Gamma] A  MVt  )+(\[Gamma] A  MtVt  )+(\[Gamma] A  MproVt  )+(\[Gamma] A  MprotV  )+(\[Gamma] A  MprotVt  )+(\[Gamma] A  MtV  );

eqMprot [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_]:=(\[Eta] MprotVt )+(\[Eta] MprotV )+( \[Epsilon] \[Sigma] \[Chi]  Aitv )-(\[Psi] Mprot )-(\[Psi] Mprot )-(\[Phi]  Mprot )-(\[Gamma]  Mprot   A )-(\[Delta]  V   Mprot )-(\[Delta]  Vt   Mprot );

eqMtV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_]:=  (\[Delta]   V   Mt )-(\[Gamma] A MtV )-(\[Eta] MtV );

eqMVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_]:=  (\[Delta]   Vt   M )-(\[Gamma] A MVt )-(\[Eta] MVt );

eqMtVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_]:=  (\[Delta]  Vt   Mt )-(\[Gamma] A MtVt )-(\[Eta] MtVt );

eqMproVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_]:=  (\[Delta]  Vt   Mpro )-(\[Gamma] A MproVt )-(\[Eta] MproVt );

eqMproV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_]:=  (\[Delta]  V   Mpro )-(\[Gamma] A MproV )-(\[Eta] MproV );

eqMV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_]:= (\[Delta]    V   M )-(\[Gamma] A MV )-(\[Eta] MV );

eqMprotV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_]:=  (\[Delta]  V   Mprot )-(\[Gamma] A MprotV )-(\[Eta] MprotV );

eqMprotVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Epsilon]_]:=  (\[Delta]  Vt   Mprot )-(\[Gamma] A MprotVt )-(\[Eta] MprotVt );



(* ::Input::Initialization:: *)
sensitivity2={};
parametervalue={0.99,0.97,0.95};
For[\[Epsilon]=1,\[Epsilon]<=Length[parametervalue],\[Epsilon]++,(*Solve for the fixed point*)fp=NSolve[{eqA[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]] == 0, eqM[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0, eqV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0, eqAi[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0, eqAiv[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0, eqMpro[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0, eqMt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0, eqVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0, eqAit[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0, eqAitv[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0, eqMprot[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0,
eqMtV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0,
eqMVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0,
eqMtVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0,
eqMproVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0,eqMproV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0,
eqMV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0,eqMprotV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0,
eqMprotVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Epsilon]]]]==0},{A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt},Reals];AppendTo[sensitivity2,fp]];
Export["sensitivity2",sensitivity2]


(* ::Text::Initialization:: *)
(*\[Kappa] sensitivity analysis*)


(* ::Input::Initialization:: *)
eqA[ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_] :=  (A    (r-(\[Beta]1  L)))-(\[Gamma] M    A  )-(\[Gamma]  Mpro    A  )-(\[Gamma]  Mt    A  )-(\[Gamma]  Mprot    A  )-(\[Gamma] A  MproV  )-(\[Gamma] A  MV  )-(\[Gamma] A  MVt  )-(\[Gamma] A  MtVt  )-(\[Gamma] A  MproVt  )-(\[Gamma] A  MprotV  )-(\[Gamma] A  MprotVt  );

eqAi [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_] := (\[Gamma]  M    A  )-(\[Beta]2 L Ai  );

eqM [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_] := (\[Eta] MV  )+(\[Eta] MVt  )+(\[Sigma] \[Epsilon] Ait   )+(\[Epsilon] \[Chi] \[Sigma]  Aitv  )+( \[Sigma]  Ai  )+( \[Chi] \[Sigma] Aiv  )+(\[Psi] Mpro  )+(\[Psi] Mt  )-(\[Gamma]  M    A  )-(\[Delta]   M    V  )-(\[Delta]   M  Vt  )-(\[Phi]  M  );

eqV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_] := (\[Eta] MprotV  )+(\[Eta] MV  )+(\[Eta] MproV  )+(\[Eta] MtV  )+(\[Kappa] \[Omega]  Aitv  )+( \[Omega]  Aiv  )+(\[Psi] Vt  )-(\[Delta]   M    V  )-(\[Delta]  Mt    V  )-(\[Delta]  V    Mprot  )-(\[Delta]  V    Mpro  )-(\[Lambda]  V  );

eqAiv [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_] := (\[Gamma]  Mpro    A  )-(\[Beta]3 L Aiv  )+(\[Gamma] A  MproV  )+(\[Gamma] A  MV  );

eqMpro [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_] := (\[Eta] MproV )+(\[Eta] MproVt )+(\[Epsilon] \[Chi] \[Sigma]  Aitv )+(\[Chi] \[Omega]  Aiv )+(\[Psi] Mprot )-(\[Phi]  Mpro )-(\[Gamma]  Mpro   A )-(\[Psi] Mpro )-(\[Delta]    Vt   Mpro )-(\[Delta]  V   Mpro );

eqMt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_] := (\[Eta] MtVt )+(\[Eta] MtV )+( \[Sigma] \[Epsilon] Ait  )+( \[Epsilon] \[Chi] \[Sigma]  Aitv )+(\[Psi] Mprot )-(\[Gamma]  Mt   A )-(\[Delta]  Mt   Vt )-(\[Delta]    V   Mt )-(\[Phi]  Mt )-(\[Psi] Mt );

eqVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_] := (\[Eta] MprotVt )+(\[Eta] MproVt )+(\[Eta] MtVt )+(\[Eta] MVt )+( \[Kappa] \[Omega]  Aitv )-(\[Delta]  M   Vt )-(\[Delta]  Mt   Vt )-(\[Delta]  Vt   Mpro )-(\[Delta]  Vt   Mprot )-(\[Lambda]  Vt )-(\[Psi] Vt );

eqAit  [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_] :=(\[Gamma]  Mt    A  )-(\[Beta]2 L Ait  );

eqAitv [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_] := (\[Gamma]  Mprot    A  )-(\[Beta]3 L Aitv  )+(\[Gamma] A  MVt  )+(\[Gamma] A  MtVt  )+(\[Gamma] A  MproVt  )+(\[Gamma] A  MprotV  )+(\[Gamma] A  MprotVt  )+(\[Gamma] A  MtV  );

eqMprot [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_]:=(\[Eta] MprotVt )+(\[Eta] MprotV )+( \[Epsilon] \[Sigma] \[Chi]  Aitv )-(\[Psi] Mprot )-(\[Psi] Mprot )-(\[Phi]  Mprot )-(\[Gamma]  Mprot   A )-(\[Delta]  V   Mprot )-(\[Delta]  Vt   Mprot );

eqMtV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_]:=  (\[Delta]   V   Mt )-(\[Gamma] A MtV )-(\[Eta] MtV );

eqMVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_]:=  (\[Delta]   Vt   M )-(\[Gamma] A MVt )-(\[Eta] MVt );

eqMtVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_]:=  (\[Delta]  Vt   Mt )-(\[Gamma] A MtVt )-(\[Eta] MtVt );

eqMproVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_]:=  (\[Delta]  Vt   Mpro )-(\[Gamma] A MproVt )-(\[Eta] MproVt );

eqMproV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_]:=  (\[Delta]  V   Mpro )-(\[Gamma] A MproV )-(\[Eta] MproV );

eqMV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_]:= (\[Delta]    V   M )-(\[Gamma] A MV )-(\[Eta] MV );

eqMprotV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_]:=  (\[Delta]  V   Mprot )-(\[Gamma] A MprotV )-(\[Eta] MprotV );

eqMprotVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_,\[Kappa]_]:=  (\[Delta]  Vt   Mprot )-(\[Gamma] A MprotVt )-(\[Eta] MprotVt );



(* ::Input::Initialization:: *)
sensitivity3={};
parametervalue={0.99,0.97,0.95};
For[\[Kappa]=1,\[Kappa]<=Length[parametervalue],\[Kappa]++,(*Solve for the fixed point*)fp=NSolve[{eqA[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]] == 0, eqM[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0, eqV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0, eqAi[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0, eqAiv[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0, eqMpro[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0, eqMt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0, eqVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0, eqAit[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0, eqAitv[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0, eqMprot[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0,
eqMtV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0,
eqMVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0,
eqMtVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0,
eqMproVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0,eqMproV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0,
eqMV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0,eqMprotV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0,
eqMprotVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt,parametervalue[[\[Kappa]]]]==0},{A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt},Reals];AppendTo[sensitivity3,fp]];
Export["sensitivity3",sensitivity3]


(* ::Text::Initialization:: *)
(*Sensitivity Analysis Data processing*)


(* ::Input::Initialization:: *)
s1 =Import["/home/kai-hua/OneDrive/Thesis/sensitivity1.csv"];
processedDatas1=ToExpression[s1];
selecteds11less=Select[processedDatas1[[1]],And@@Thread[Values[#]>=0]&];
selecteds13less=Select[processedDatas1[[2]],And@@Thread[Values[#]>=0]&];
selecteds15less=Select[processedDatas1[[3]],And@@Thread[Values[#]>=0]&];
selecteds11greater=Select[processedDatas1[[4]],And@@Thread[Values[#]>=0]&];
selecteds13greater=Select[processedDatas1[[5]],And@@Thread[Values[#]>=0]&];
selecteds15greater=Select[processedDatas1[[6]],And@@Thread[Values[#]>=0]&];
 fixedpointss11less=DeleteDuplicatesBy[selecteds11less,ToString];
 fixedpointss13less=DeleteDuplicatesBy[selecteds13less,ToString];
 fixedpointss15less=DeleteDuplicatesBy[selecteds15less,ToString];
 fixedpointss11greater=DeleteDuplicatesBy[selecteds11greater,ToString];
 fixedpointss13greater=DeleteDuplicatesBy[selecteds13greater,ToString];
 fixedpointss15greater=DeleteDuplicatesBy[selecteds15greater,ToString] ;


(* ::Input::Initialization:: *)
s2 =Import["/home/kai-hua/OneDrive/Thesis/sensitivity2.csv"];
processedDatas2=ToExpression[s2];
selecteds21=Select[processedDatas2[[1]],And@@Thread[Values[#]>=0]&];
selecteds23=Select[processedDatas2[[2]],And@@Thread[Values[#]>=0]&];
selecteds25=Select[processedDatas2[[3]],And@@Thread[Values[#]>=0]&];
 fixedpointss21=DeleteDuplicatesBy[selecteds21,ToString];
 fixedpointss23=DeleteDuplicatesBy[selecteds23,ToString];
 fixedpointss25=DeleteDuplicatesBy[selecteds25,ToString];


(* ::Input::Initialization:: *)
s3 =Import["/home/kai-hua/OneDrive/Thesis/sensitivity3.csv"]; (*Kappa*)
processedDatas3=ToExpression[s3];
selecteds31=Select[processedDatas3[[1]],And@@Thread[Values[#]>=0]&];(*1%*)
selecteds33=Select[processedDatas3[[2]],And@@Thread[Values[#]>=0]&];(*3%*)
selecteds35=Select[processedDatas3[[3]],And@@Thread[Values[#]>=0]&];(*5%*)
 fixedpointss31=DeleteDuplicatesBy[selecteds31,ToString];
 fixedpointss33=DeleteDuplicatesBy[selecteds33,ToString];
 fixedpointss35=DeleteDuplicatesBy[selecteds35,ToString];


(* ::Text::Initialization:: *)
(*Base Values Fixed Points*)


(* ::Input::Initialization:: *)
data =Import["/home/kai-hua/OneDrive/Thesis/steadystate.csv"];
processedData=ToExpression[data];
selected=Select[processedData,And@@Thread[Values[#]>=0]&];
 fixedpoints=DeleteDuplicatesBy[selected,ToString];


(* ::Text::Initialization:: *)
(*Sensitivity analysis  calculation*)


(* ::Input::Initialization:: *)
extractCoordinates[point_]:=Values[point]

euclideanDistance[point1_List,point2_List]:=EuclideanDistance[extractCoordinates[point1],extractCoordinates[point2]]
minDistance=Infinity;
minPair={};


closestPairss31={};
closestPairss33={};
closestPairss35={};
localsensitivityresults31={};
localsensitivityresults33={};
localsensitivityresults35={};


(* ::Text::Initialization:: *)
(*Euclidean distance approximation*)


(* ::Input::Initialization:: *)
For[i=1,i<=Length[fixedpointss31],i++,minDistance=Infinity;
minPair={};
For[j=1,j<=Length[fixedpoints],j++,distance=euclideanDistance[fixedpoints[[j]],fixedpointss31[[i]]];
If[distance<minDistance,minDistance=distance;
minPair={fixedpointss31[[i]],fixedpoints[[j]]}]];
 AppendTo[closestPairss31,minPair];]

(*Absolute local sensitivity coefficient  Numerical Approximation*)

For[i=1,i<=Length[closestPairss31],i++,approximation=(Values[closestPairss31[[i,1]]]-Values[closestPairss31[[i,2]]])/0.01;AppendTo[localsensitivityresults31,approximation]]




(* ::Input::Initialization:: *)
For[i=1,i<=Length[fixedpointss33],i++,minDistance=Infinity;
minPair={};
For[j=1,j<=Length[fixedpoints],j++,distance=euclideanDistance[fixedpoints[[j]],fixedpointss33[[i]]];
If[distance<minDistance,minDistance=distance;
minPair={fixedpointss33[[i]],fixedpoints[[j]]}]];
 AppendTo[closestPairss33,minPair];]

(*Absolute local sensitivity coefficient  Numerical Approximation*)

For[i=1,i<=Length[closestPairss33],i++,approximation=(Values[closestPairss33[[i,1]]]-Values[closestPairss33[[i,2]]])/0.03;AppendTo[localsensitivityresults33,approximation]]




(* ::Input::Initialization:: *)
For[i=1,i<=Length[fixedpointss35],i++,minDistance=Infinity;
minPair={};
For[j=1,j<=Length[fixedpoints],j++,distance=euclideanDistance[fixedpoints[[j]],fixedpointss35[[i]]];
If[distance<minDistance,minDistance=distance;
minPair={fixedpointss35[[i]],fixedpoints[[j]]}]];
 AppendTo[closestPairss35,minPair];]

(*Absolute local sensitivity coefficient  Numerical Approximation*)

For[i=1,i<=Length[closestPairss35],i++,approximation=(Values[closestPairss35[[i,1]]]-Values[closestPairss35[[i,2]]])/0.05;AppendTo[localsensitivityresults35,approximation]]


(* ::Input::Initialization:: *)
localsensitivityresults31
localsensitivityresults33
localsensitivityresults35
