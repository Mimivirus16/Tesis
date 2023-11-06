#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title::Initialization:: *)
(*In Silico Simulations*)


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


L=1+(A[t]+Ai[t]+Aiv[t]+Ait[t]+Aitv[t])/K(*Logistic Growth*);


(* ::Input::Initialization:: *)
\[Chi]=0.3 (*Control Parameter (Virophage to Mimivirus*)(*Virophages\[LongDash]Known and Unknown Facts*);
\[Epsilon]=1(*Control Parameter (Transpoviron to Mimivirus*);
\[Kappa]=1(*Control Parameter (Transpoviron to Virophage*);


(* ::Input::Initialization:: *)
zamilonsputnik1={A'[t]==(A[t]  (r-(\[Beta]1  L)))-(\[Gamma] M[t]  A[t])-(\[Gamma]  Mpro[t]  A[t])-(\[Gamma]  Mt[t]  A[t])-(\[Gamma]  Mprot[t]  A[t])-(\[Gamma] A[t]MproV[t])-(\[Gamma] A[t]MV[t])-(\[Gamma] A[t]MVt[t])-(\[Gamma] A[t]MtVt[t])-(\[Gamma] A[t]MproVt[t])-(\[Gamma] A[t]MprotV[t])-(\[Gamma] A[t]MprotVt[t]),

M'[t]==(\[Eta] MV[t])+(\[Eta] MVt[t])+( \[Sigma]  Ai[t])+( \[Chi] \[Sigma] Aiv[t])+(\[Sigma] \[Epsilon] Ait[t] )+(\[Epsilon] \[Chi] \[Sigma]  Aitv[t])+(\[Psi] Mpro[t])+(\[Psi] Mt[t])-(\[Gamma]  M[t]  A[t])-(\[Delta]   M[t]  V[t])-(\[Delta]   M[t]Vt[t])-(\[Phi]  M[t]),

V'[t]==(\[Eta] MprotV[t])+(\[Eta] MV[t])+(\[Eta] MproV[t])+(\[Eta] MtV[t])+(\[Kappa] \[Omega]  Aitv[t])+( \[Omega]  Aiv[t])+(\[Psi] Vt[t])-(\[Delta]   M[t]  V[t])-(\[Delta]  Mt[t]  V[t])-(\[Delta]  V[t]  Mprot[t])-(\[Delta]  V[t]  Mpro[t])-(\[Lambda]  V[t]),

Ai'[t]==(\[Gamma]  M[t]  A[t])-(\[Beta]2 L Ai[t]),

Aiv'[t]==(\[Gamma]  Mpro[t]  A[t])-(\[Beta]3 L Aiv[t])+(\[Gamma] A[t]MproV[t])+(\[Gamma] A[t]MV[t]),

Mpro'[t]==(\[Eta] MproV[t])+(\[Eta] MproVt[t])+(\[Epsilon] \[Chi] \[Sigma]  Aitv[t])+(\[Chi] \[Omega]  Aiv[t])+(\[Psi] Mprot[t])-(\[Phi]  Mpro[t])-(\[Gamma]  Mpro[t]  A[t])-(\[Psi] Mpro[t])-(\[Delta]    Vt[t]  Mpro[t])-(\[Delta]  V[t]  Mpro[t]),

Mt'[t]==(\[Eta] MtVt[t])+(\[Eta] MtV[t])+( \[Sigma] \[Epsilon] Ait[t] )+( \[Epsilon] \[Chi] \[Sigma]  Aitv[t])+(\[Psi] Mprot[t])-(\[Gamma]  Mt[t]  A[t])-(\[Delta]  Mt[t]  Vt[t])-(\[Delta]    V[t]  Mt[t])-(\[Phi]  Mt[t])-(\[Psi] Mt[t]),

Vt'[t]==(\[Eta] MprotVt[t])+(\[Eta] MproVt[t])+(\[Eta] MtVt[t])+(\[Eta] MVt[t])+( \[Kappa] \[Omega]  Aitv[t])-(\[Delta]  M[t]  Vt[t])-(\[Delta]  Mt[t]  Vt[t])-(\[Delta]  Vt[t]  Mpro[t])-(\[Delta]  Vt[t]  Mprot[t])-(\[Lambda]  Vt[t])-(\[Psi] Vt[t]),

Ait'[t]==(\[Gamma]  Mt[t]  A[t])-(\[Beta]2 L Ait[t]),

Aitv'[t]==(\[Gamma]  Mprot[t]  A[t])-(\[Beta]3 L Aitv[t])+(\[Gamma] A[t]MVt[t])+(\[Gamma] A[t]MtVt[t])+(\[Gamma] A[t]MproVt[t])+(\[Gamma] A[t]MprotV[t])+(\[Gamma] A[t]MprotVt[t])+(\[Gamma] A[t]MtV[t]),

Mprot'[t]==(\[Eta] MprotVt[t])+(\[Eta] MprotV[t])+( \[Epsilon] \[Sigma] \[Chi]  Aitv[t])-(\[Psi] Mprot[t])-(\[Psi] Mprot[t])-(\[Phi]  Mprot[t])-(\[Gamma]  Mprot[t]  A[t])-(\[Delta]  V[t]  Mprot[t])-(\[Delta]  Vt[t]  Mprot[t]),

MtV'[t]== (\[Delta]   V[t]  Mt[t])-(\[Gamma] A[t]MtV[t])-(\[Eta] MtV[t]),
MVt'[t]== (\[Delta]   Vt[t]  M[t])-(\[Gamma] A[t]MVt[t])-(\[Eta] MVt[t]),
MtVt'[t]== (\[Delta]  Vt[t]  Mt[t])-(\[Gamma] A[t]MtVt[t])-(\[Eta] MtVt[t]),
MproVt'[t]== (\[Delta]  Vt[t]  Mpro[t])-(\[Gamma] A[t]MproVt[t])-(\[Eta] MproVt[t]),
MproV'[t]== (\[Delta]  V[t]  Mpro[t])-(\[Gamma] A[t]MproV[t])-(\[Eta] MproV[t]),
MV'[t]== (\[Delta]    V[t]  M[t])-(\[Gamma] A[t]MV[t])-(\[Eta] MV[t]),
MprotV'[t]== (\[Delta]  V[t]  Mprot[t])-(\[Gamma] A[t]MprotV[t])-(\[Eta] MprotV[t]),
MprotVt'[t]== (\[Delta]  Vt[t]  Mprot[t])-(\[Gamma] A[t]MprotVt[t])-(\[Eta] MprotVt[t])};


(* ::Text::Initialization:: *)
(*Initial conditions*)


(* ::Input::Initialization:: *)
literature={A[0]==1*^6,M[0]== 10*^6,V[0]==0 ,Ai[0]== 0,Aiv[0]== 0,Mpro[0]== 0,Mt[0]== 0,Vt[0]== 10*^9,Ait[0]== 0,Aitv[0]==  0,Mprot[0]== 0,MtV[0]==0,MVt[0]==0,MtVt[0]==0,MproVt[0]==0,MproV[0]==0,MV[0]==0,MprotV[0]==0,MprotVt[0]==0};(*MOI=10 for Mimivirus,10^3 for 1 Mimivirus particle (Virophage with Transpoviron)*)
tmax=730;
tmin=0;

variables={A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt};


(* ::Subtitle::Initialization:: *)
(*Literature Simulation*)


(* ::Input::Initialization:: *)
simulation=NDSolve[{zamilonsputnik1,literature},{A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt},{t,tmax},Method->"StiffnessSwitching"];
Plot[Evaluate[{A[t],M[t],V[t],Mt[t],Vt[t]}/. simulation],{t,0,tmax},AxesLabel->Automatic,PlotStyle->{Red,Green,Blue,Purple,Yellow,Orange,Gray,Black,Cyan,Pink,Magenta},PlotLabels->{"A[t]","M[t]","V[t]","Mt[t]","Vt[t]"},PlotLegends->LineLegend[{Red,Green,Blue,Purple,Yellow,Orange,Gray,Black,Cyan,Pink,Magenta},{"A[t]","M[t]","V[t]","Mt[t]","Vt[t]"}]]

literaturesimulation=Table[Flatten[{t,Evaluate[{A[t],M[t],V[t],Ai[t],Aiv[t],Mpro[t],Mt[t],Vt[t],Ait[t],Aitv[t],Mprot[t],MtV[t],MVt[t],MtVt[t],MproVt[t],MproV[t],MV[t],MprotV[t],MprotVt[t]}/.simulation]}],{t,tmin,tmax}];
tableWithTitles=Insert[literaturesimulation,Prepend[Array[ToString[variables[[#]]]&,19],"Time"],1];

(*Export the table to a CSV file*)
Export["literaturesimulation.csv",tableWithTitles]


(* ::Subtitle::Initialization:: *)
(*Bifurcation Analysis Unstable oscillator Virophage Simulation*)


(* ::Input::Initialization:: *)
fixedpointsbi=ToExpression[Import["/home/kai-hua/fpbi.csv"]];
virophage=fixedpointsbi[[5841,4]];


(* ::Input::Initialization:: *)



(* ::Input::Initialization:: *)
virophage={A->465.478,M->3.710431561117332`*^6,V->1.181103015052778`*^6,Ai->82.67753334265693`,Aiv->314.97950513349065`,Mpro->1.1446329146778325`*^7,Mt->0,Vt->0,Ait->0,Aitv->0,Mprot->0,MtV->0,MVt->0,MtVt->0,MproVt->0,MproV->844907.720014968`,MV->273884.51182687125`,MprotV->0,MprotVt->0};
(*Change in -1 A*)
virophagefps=virophage/.A-> A[0] /.M->M[0]  /.V->V[0]  /.Ai-> Ai[0] /.Aiv->Aiv[0]  /.Mpro->Mpro[0]  /.Mt->Mt[0]  /.Vt->Vt[0]  /.Ait->Ait[0]  /.Aitv->Aitv[0]  /.Mprot->Mprot[0]  /.MtV->MtV[0]  /.MVt->MVt[0]  /.MtVt->MtVt[0]  /.MproVt-> MproVt[0] /.MproV->MproV[0]  /.MV->MV[0]  /.MprotV->MprotV[0] /.MprotVt->MprotVt[0] /.Rule->Equal;


(* ::Input::Initialization:: *)
\[Chi]=fixedpointsbi[[5841,1]](*Control Parameter (Virophage to Mimivirus*)(*Virophages\[LongDash]Known and Unknown Facts*)
\[Epsilon]=fixedpointsbi[[5841,2]](*Control Parameter (Transpoviron to Mimivirus*)
\[Kappa]=fixedpointsbi[[5841,3]](*Control Parameter (Transpoviron to Virophage*)



simulation2=NDSolve[{zamilonsputnik1,virophagefps},{A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt},{t,tmax},Method->"StiffnessSwitching"];
Plot[Evaluate[{A[t],M[t],V[t],Mt[t],Vt[t]}/. simulation2],{t,0,tmax},AxesLabel->Automatic,PlotStyle->{Red,Green,Blue,Purple,Yellow,Orange,Gray,Black,Cyan,Pink,Magenta},PlotLabels->{"A[t]","M[t]","V[t]","Mt[t]","Vt[t]"},PlotLegends->LineLegend[{Red,Green,Blue,Purple,Yellow,Orange,Gray,Black,Cyan,Pink,Magenta},{"A[t]","M[t]","V[t]","Mt[t]","Vt[t]"}]]

unstablevirophagesimulation=Table[Flatten[{t,Evaluate[{A[t],M[t],V[t],Ai[t],Aiv[t],Mpro[t],Mt[t],Vt[t],Ait[t],Aitv[t],Mprot[t],MtV[t],MVt[t],MtVt[t],MproVt[t],MproV[t],MV[t],MprotV[t],MprotVt[t]}/.simulation2]}],{t,tmin,tmax}];
tableWithTitles2=Insert[unstablevirophagesimulation,Prepend[Array[ToString[variables[[#]]]&,19],"Time"],1];

(*Export the table to a CSV file*)
Export["unstablevirophagesimulation.csv",tableWithTitles2]


(* ::Subtitle::Initialization:: *)
(*Sensitivity Analysis Transpoviron Simulation*)


(* ::Input::Initialization:: *)
\[Chi]=0.3;
\[Epsilon]=1;
\[Kappa]=0.99;
tmax=500;
s3 =Import["/home/kai-hua/OneDrive/Thesis/sensitivity3.csv"];
processedDatas3=ToExpression[s3];
selecteds31=Select[processedDatas3[[1]],And@@Thread[Values[#]>=0]&];
 fixedpointss31=DeleteDuplicatesBy[selecteds31,ToString];
transpoviron=fixedpointss31[[3]];

transpovironfps=transpoviron/.A-> A[0] /.M->M[0]  /.V->V[0]  /.Ai-> Ai[0] /.Aiv->Aiv[0]  /.Mpro->Mpro[0]  /.Mt->Mt[0]  /.Vt->Vt[0]  /.Ait->Ait[0]  /.Aitv->Aitv[0]  /.Mprot->Mprot[0]  /.MtV->MtV[0]  /.MVt->MVt[0]  /.MtVt->MtVt[0]  /.MproVt-> MproVt[0] /.MproV->MproV[0]  /.MV->MV[0]  /.MprotV->MprotV[0] /.MprotVt->MprotVt[0] /.Rule->Equal;


(* ::Input::Initialization:: *)
simulation3=NDSolve[{zamilonsputnik1,transpovironfps},{A ,M ,V ,Ai ,Aiv ,Mpro ,Mt ,Vt ,Ait ,Aitv,Mprot,MtV,MVt,MtVt,MproVt,MproV,MV,MprotV,MprotVt},{t,tmax},Method->"StiffnessSwitching"];
Plot[Evaluate[{A[t],M[t],V[t],Mt[t],Vt[t],Ai[t]}/. simulation3],{t,0,tmax},AxesLabel->Automatic,PlotStyle->{Red,Green,Blue,Purple,Yellow,Orange,Gray,Black,Cyan,Pink,Magenta},PlotLabels->{"A[t]","M[t]","V[t]","Mt[t]","Vt[t]","Ai[t]"},PlotLegends->LineLegend[{Red,Green,Blue,Purple,Yellow,Orange,Gray,Black,Cyan,Pink,Magenta},{"A[t]","M[t]","V[t]","Mt[t]","Vt[t]","Ai[t]"}]]

transpovironsimulation=Table[Flatten[{t,Evaluate[{A[t],M[t],V[t],Ai[t],Aiv[t],Mpro[t],Mt[t],Vt[t],Ait[t],Aitv[t],Mprot[t],MtV[t],MVt[t],MtVt[t],MproVt[t],MproV[t],MV[t],MprotV[t],MprotVt[t]}/.simulation3]}],{t,tmin,tmax}];
tableWithTitles3=Insert[transpovironsimulation,Prepend[Array[ToString[variables[[#]]]&,19],"Time"],1];

(*Export the table to a CSV file*)
Export["transpovironsimulation.csv",tableWithTitles3]
