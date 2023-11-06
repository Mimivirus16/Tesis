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

\[Chi]=0.3 (*Control Parameter (Virophage to Mimivirus*)(*Virophages\[LongDash]Known and Unknown Facts*);
\[Epsilon]=1(*Control Parameter (Transpoviron to Mimivirus*);
\[Kappa]=1(*Control Parameter (Transpoviron to Virophage*);


(* ::Input::Initialization:: *)
eqA[ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_] :=  (A    (r-(\[Beta]1  L)))-(\[Gamma] M    A  )-(\[Gamma]  Mpro    A  )-(\[Gamma]  Mt    A  )-(\[Gamma]  Mprot    A  )-(\[Gamma] A  MproV  )-(\[Gamma] A  MV  )-(\[Gamma] A  MVt  )-(\[Gamma] A  MtVt  )-(\[Gamma] A  MproVt  )-(\[Gamma] A  MprotV  )-(\[Gamma] A  MprotVt  );

eqAi [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_] := (\[Gamma]  M    A  )-(\[Beta]2 L Ai  );

eqM [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_] := (\[Eta] MV  )+(\[Eta] MVt  )+(\[Sigma] \[Epsilon] Ait   )+(\[Epsilon] \[Chi] \[Sigma]  Aitv  )+( \[Sigma]  Ai  )+( \[Chi] \[Sigma] Aiv  )+(\[Psi] Mpro  )+(\[Psi] Mt  )-(\[Gamma]  M    A  )-(\[Delta]   M    V  )-(\[Delta]   M  Vt  )-(\[Phi]  M  );

eqV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_] := (\[Eta] MprotV  )+(\[Eta] MV  )+(\[Eta] MproV  )+(\[Eta] MtV  )+(\[Kappa] \[Omega]  Aitv  )+( \[Omega]  Aiv  )+(\[Psi] Vt  )-(\[Delta]   M    V  )-(\[Delta]  Mt    V  )-(\[Delta]  V    Mprot  )-(\[Delta]  V    Mpro  )-(\[Lambda]  V  );

eqAiv [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_] := (\[Gamma]  Mpro    A  )-(\[Beta]3 L Aiv  )+(\[Gamma] A  MproV  )+(\[Gamma] A  MV  );

eqMpro [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_] := (\[Eta] MproV )+(\[Eta] MproVt )+(\[Epsilon] \[Chi] \[Sigma]  Aitv )+(\[Chi] \[Omega]  Aiv )+(\[Psi] Mprot )-(\[Phi]  Mpro )-(\[Gamma]  Mpro   A )-(\[Psi] Mpro )-(\[Delta]    Vt   Mpro )-(\[Delta]  V   Mpro );

eqMt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_] := (\[Eta] MtVt )+(\[Eta] MtV )+( \[Sigma] \[Epsilon] Ait  )+( \[Epsilon] \[Chi] \[Sigma]  Aitv )+(\[Psi] Mprot )-(\[Gamma]  Mt   A )-(\[Delta]  Mt   Vt )-(\[Delta]    V   Mt )-(\[Phi]  Mt )-(\[Psi] Mt );

eqVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_] := (\[Eta] MprotVt )+(\[Eta] MproVt )+(\[Eta] MtVt )+(\[Eta] MVt )+( \[Kappa] \[Omega]  Aitv )-(\[Delta]  M   Vt )-(\[Delta]  Mt   Vt )-(\[Delta]  Vt   Mpro )-(\[Delta]  Vt   Mprot )-(\[Lambda]  Vt )-(\[Psi] Vt );

eqAit  [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_] :=(\[Gamma]  Mt    A  )-(\[Beta]2 L Ait  );

eqAitv [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_] := (\[Gamma]  Mprot    A  )-(\[Beta]3 L Aitv  )+(\[Gamma] A  MVt  )+(\[Gamma] A  MtVt  )+(\[Gamma] A  MproVt  )+(\[Gamma] A  MprotV  )+(\[Gamma] A  MprotVt  )+(\[Gamma] A  MtV  );

eqMprot [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_]:=(\[Eta] MprotVt )+(\[Eta] MprotV )+( \[Epsilon] \[Sigma] \[Chi]  Aitv )-(\[Psi] Mprot )-(\[Psi] Mprot )-(\[Phi]  Mprot )-(\[Gamma]  Mprot   A )-(\[Delta]  V   Mprot )-(\[Delta]  Vt   Mprot );

eqMtV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_]:=  (\[Delta]   V   Mt )-(\[Gamma] A MtV )-(\[Eta] MtV );

eqMVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_]:=  (\[Delta]   Vt   M )-(\[Gamma] A MVt )-(\[Eta] MVt );

eqMtVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_]:=  (\[Delta]  Vt   Mt )-(\[Gamma] A MtVt )-(\[Eta] MtVt );

eqMproVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_]:=  (\[Delta]  Vt   Mpro )-(\[Gamma] A MproVt )-(\[Eta] MproVt );

eqMproV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_]:=  (\[Delta]  V   Mpro )-(\[Gamma] A MproV )-(\[Eta] MproV );

eqMV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_]:= (\[Delta]    V   M )-(\[Gamma] A MV )-(\[Eta] MV );

eqMprotV [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_]:=  (\[Delta]  V   Mprot )-(\[Gamma] A MprotV )-(\[Eta] MprotV );

eqMprotVt [ A_,M_,V_,Ai_,Aiv_,Mpro_,Mt_,Vt_,Ait_,Aitv_,Mprot_,MtV_,MVt_,MtVt_,MproVt_,MproV_,MV_,MprotV_,MprotVt_]:=  (\[Delta]  Vt   Mprot )-(\[Gamma] A MprotVt )-(\[Eta] MprotVt );



(* ::Input::Initialization:: *)
    steadystate=NSolve[{eqA[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt] == 0, eqM[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0, eqV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0, eqAi[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0, eqAiv[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0, eqMpro[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0, eqMt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0, eqVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0, eqAit[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0, eqAitv[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0, eqMprot[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0,
eqMtV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0,
eqMVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0,
eqMtVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0,
eqMproVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0,eqMproV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0,
eqMV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0,eqMprotV[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0,
eqMprotVt[   A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt]==0},{A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt},Reals]
Export["steadystate.csv",steadystate]
