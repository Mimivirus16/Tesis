#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
ClearAll["Global`*"]
K=4*^6 ;
r=1.4 ;
\[Beta]1=0.7;
\[Beta]2=0.9;
\[Beta]3=0.8;
\[Gamma]=4.3*^-8;
\[Delta]=2.2*^-8;
\[Phi]=3.2*^-2 ;
\[Lambda]=3.2*^-1 ;
\[Eta]=(\[Lambda]+\[Phi]) ;
\[Psi]=0.001 ;
\[Sigma]=300 *\[Beta]2;
\[Omega]=1500*\[Beta]3;
\[Zeta]=30*\[Beta]3;
L=1+(A[t]+Ai[t]+Aiv[t]+Ait[t]+Aitv[t])/K ;
fixedpointsbi=ToExpression[Import["/home/kai-hua/fpbi.csv"]];
vmatrix:={(A[t]  (r-(\[Beta]1  L)))-(\[Gamma] M[t]  A[t])-(\[Gamma]  Mpro[t]  A[t])-(\[Gamma]  Mt[t]  A[t])-(\[Gamma]  Mprot[t]  A[t])-(\[Gamma] A[t]MproV[t])-(\[Gamma] A[t]MV[t])-(\[Gamma] A[t]MVt[t])-(\[Gamma] A[t]MtVt[t])-(\[Gamma] A[t]MproVt[t])-(\[Gamma] A[t]MprotV[t])-(\[Gamma] A[t]MprotVt[t]),

(\[Eta] MV[t])+(\[Eta] MVt[t])+(\[Sigma] \[Epsilon] Ait[t] )+(\[Epsilon] \[Chi] \[Sigma]  Aitv[t])+( \[Sigma]  Ai[t])+( \[Chi] \[Sigma] Aiv[t])+(\[Psi] Mpro[t])+(\[Psi] Mt[t])-(\[Gamma]  M[t]  A[t])-(\[Delta]   M[t]  V[t])-(\[Delta]   M[t]Vt[t])-(\[Phi]  M[t]),

(\[Eta] MprotV[t])+(\[Eta] MV[t])+(\[Eta] MproV[t])+(\[Eta] MtV[t])+(\[Kappa] \[Omega]  Aitv[t])+( \[Omega]  Aiv[t])+(\[Psi] Vt[t])-(\[Delta]   M[t]  V[t])-(\[Delta]  Mt[t]  V[t])-(\[Delta]  V[t]  Mprot[t])-(\[Delta]  V[t]  Mpro[t])-(\[Lambda]  V[t]),

(\[Gamma]  M[t]  A[t])-(\[Beta]2 L Ai[t]),

(\[Gamma]  Mpro[t]  A[t])-(\[Beta]3 L Aiv[t])+(\[Gamma] A[t]MproV[t])+(\[Gamma] A[t]MV[t]),

(\[Eta] MproV[t])+(\[Eta] MproVt[t])+(\[Epsilon] \[Chi] \[Sigma]  Aitv[t])+(\[Chi] \[Omega]  Aiv[t])+(\[Psi] Mprot[t])-(\[Phi]  Mpro[t])-(\[Gamma]  Mpro[t]  A[t])-(\[Psi] Mpro[t])-(\[Delta]    Vt[t]  Mpro[t])-(\[Delta]  V[t]  Mpro[t]),

(\[Eta] MtVt[t])+(\[Eta] MtV[t])+( \[Sigma] \[Epsilon] Ait[t] )+( \[Epsilon] \[Chi] \[Sigma]  Aitv[t])+(\[Psi] Mprot[t])-(\[Gamma]  Mt[t]  A[t])-(\[Delta]  Mt[t]  Vt[t])-(\[Delta]    V[t]  Mt[t])-(\[Phi]  Mt[t])-(\[Psi] Mt[t]),

(\[Eta] MprotVt[t])+(\[Eta] MproVt[t])+(\[Eta] MtVt[t])+(\[Eta] MVt[t])+( \[Kappa] \[Omega]  Aitv[t])-(\[Delta]  M[t]  Vt[t])-(\[Delta]  Mt[t]  Vt[t])-(\[Delta]  Vt[t]  Mpro[t])-(\[Delta]  Vt[t]  Mprot[t])-(\[Lambda]  Vt[t])-(\[Psi] Vt[t]),

(\[Gamma]  Mt[t]  A[t])-(\[Beta]2 L Ait[t]),

(\[Gamma]  Mprot[t]  A[t])-(\[Beta]3 L Aitv[t])+(\[Gamma] A[t]MVt[t])+(\[Gamma] A[t]MtVt[t])+(\[Gamma] A[t]MproVt[t])+(\[Gamma] A[t]MprotV[t])+(\[Gamma] A[t]MprotVt[t])+(\[Gamma] A[t]MtV[t]),

(\[Eta] MprotVt[t])+(\[Eta] MprotV[t])+( \[Epsilon] \[Sigma] \[Chi]  Aitv[t])-(\[Psi] Mprot[t])-(\[Psi] Mprot[t])-(\[Phi]  Mprot[t])-(\[Gamma]  Mprot[t]  A[t])-(\[Delta]  V[t]  Mprot[t])-(\[Delta]  Vt[t]  Mprot[t]),

(\[Delta]   V[t]  Mt[t])-(\[Gamma] A[t]MtV[t])-(\[Eta] MtV[t]),
(\[Delta]   Vt[t]  M[t])-(\[Gamma] A[t]MVt[t])-(\[Eta] MVt[t]),
 (\[Delta]  Vt[t]  Mt[t])-(\[Gamma] A[t]MtVt[t])-(\[Eta] MtVt[t]),
 (\[Delta]  Vt[t]  Mpro[t])-(\[Gamma] A[t]MproVt[t])-(\[Eta] MproVt[t]),
(\[Delta]  V[t]  Mpro[t])-(\[Gamma] A[t]MproV[t])-(\[Eta] MproV[t]),
(\[Delta]    V[t]  M[t])-(\[Gamma] A[t]MV[t])-(\[Eta] MV[t]),
(\[Delta]  V[t]  Mprot[t])-(\[Gamma] A[t]MprotV[t])-(\[Eta] MprotV[t]),
(\[Delta]  Vt[t]  Mprot[t])-(\[Gamma] A[t]MprotVt[t])-(\[Eta] MprotVt[t])};
variables={A[t],M[t],V[t],Ai[t],Aiv[t],Mpro[t],Mt[t],Vt[t],Ait[t],Aitv[t],Mprot[t],MtV[t],MVt[t],MtVt[t],MproVt[t],MproV[t],MV[t],MprotV[t],MprotVt[t]};
jacobian={};
For[i=1,i<=Length[variables],i++,AppendTo[jacobian,D[vmatrix,variables[[i]]]]]
jacobian;
jacobianm=jacobian/. x_[t]:>x;
ji={};
eigen={};
ilength=Length[fixedpointsbi];
For[i=1,i<=ilength,i++,AppendTo[ji,jacobianm/.fixedpointsbi[[i,4]]]]


(* ::Input::Initialization:: *)
elength=Length[ji];


(* ::Text::Initialization:: *)
(*Bifurcation Analysis fixed point classification*)


(* ::Input::Initialization:: *)
eigen={};
For[i=1,i<=elength,i++,
\[Chi]=fixedpointsbi[[i,1]];\[Epsilon]=fixedpointsbi[[i,2]];\[Kappa]=fixedpointsbi[[i,3]];AppendTo[eigen,Eigenvalues[ji[[i]]]]]


(* ::Input::Initialization:: *)
dynamicalbehaviour={};
For[i=1,i<=5842,i++,If[AllTrue[eigen[[i]],Re[#]==0&&Im[#]==0&],systembehaviour="Neutral Point",
If[AllTrue[eigen[[i]],Re[#]==0&&Im[#]!=0&],systembehaviour="Neutral Center",If[AllTrue[eigen[[i]],Re[#]<0&&Im[#]==0&],systembehaviour="Stable Point",
If[AllTrue[eigen[[i]],Re[#]<0&&(Im[#]!=0||Im[#]==0)&],systembehaviour="Stable Oscillator",If[AllTrue[eigen[[i]],Re[#]>0&&Im[#]==0&],systembehaviour="Unstable Point",If[AllTrue[eigen[[i]],Re[#]>0&&Im[#]!=0&],systembehaviour="Unstable Oscillator",If[AllTrue[eigen[[i]],(Re[#]>0||Re[#]<0||Re[#]==0)&&Im[#]==0&],systembehaviour="Saddle Point",systembehaviour="Unstable Oscillator"
]]]]]]];AppendTo[dynamicalbehaviour,systembehaviour];]


(* ::Text::Initialization:: *)
(*Fixed point populations*)


(* ::Input::Initialization:: *)
dynamicalpopulation={};

For[i=1,i<=5842,i++,If[Positive[A/. fixedpointsbi[[i,4]]]&&Positive[M/.fixedpointsbi[[i,4]]]&&Positive[Ai/.fixedpointsbi[[i,4]]]&&Positive[V/.fixedpointsbi[[i,4]]]&&Positive[Aiv/.fixedpointsbi[[i,4]]]&&Positive[Aitv/.fixedpointsbi[[i,4]]]&&Positive[Ait/.fixedpointsbi[[i,4]]]&&Positive[Mt/.fixedpointsbi[[i,4]]]&&Positive[Vt/.fixedpointsbi[[i,4]]]&&Positive[Mprot/.fixedpointsbi[[i,4]]]&&Positive[Mpro/.fixedpointsbi[[i,4]]],systempopulation="AMPVT",If[Positive[A/. fixedpointsbi[[i,4]]]&&Positive[M/.fixedpointsbi[[i,4]]]&&Positive[Ai/.fixedpointsbi[[i,4]]]&&Positive[V/.fixedpointsbi[[i,4]]]&&Positive[Aiv/.fixedpointsbi[[i,4]]]&&Positive[Aitv/.fixedpointsbi[[i,4]]]&&Positive[Ait/.fixedpointsbi[[i,4]]]&&Positive[Mt/.fixedpointsbi[[i,4]]]&&Positive[Vt/.fixedpointsbi[[i,4]]],systempopulation="AMVT",If[Positive[A/. fixedpointsbi[[i,4]]]&&Positive[M/.fixedpointsbi[[i,4]]]&&Positive[Ai/.fixedpointsbi[[i,4]]]&&Positive[Ait/.fixedpointsbi[[i,4]]]&&Positive[Mt/.fixedpointsbi[[i,4]]],systempopulation="AMT",If[Positive[A/. fixedpointsbi[[i,4]]]&&Positive[M/.fixedpointsbi[[i,4]]]&&Positive[Ai/.fixedpointsbi[[i,4]]]&&Positive[V/.fixedpointsbi[[i,4]]]&&Positive[Aiv/.fixedpointsbi[[i,4]]]&&Positive[Mpro/.fixedpointsbi[[i,4]]],systempopulation="AMPV",If[Positive[A/. fixedpointsbi[[i,4]]]&&Positive[M/.fixedpointsbi[[i,4]]]&&Positive[Ai/.fixedpointsbi[[i,4]]]&&Positive[V/.fixedpointsbi[[i,4]]]&&Positive[Aiv/.fixedpointsbi[[i,4]]],systempopulation="AMV",If[Positive[A/. fixedpointsbi[[i,4]]]&&Positive[M/.fixedpointsbi[[i,4]]]&&Positive[Ai/.fixedpointsbi[[i,4]]],systempopulation="AM",If[Positive[A/. fixedpointsbi[[i,4]]],systempopulation="A",systempopulation="E"]]]]]]];AppendTo[dynamicalpopulation,systempopulation];]



(* ::Text::Initialization:: *)
(*R data export*)


(* ::Input::Initialization:: *)
Export["dynamicpopulation.csv",dynamicalpopulation];
Export["dynamicbehaviour.csv",dynamicalbehaviour];


(* ::Text::Initialization:: *)
(* R Data Import*)


(* ::Input::Initialization:: *)
bifurdata=Rest@Import["/home/kai-hua/Bifurcationdata.csv", HeaderLines -> 0]; (*Transformed CSV from R*)
bifurdata=MapAt[ToExpression,bifurdata,{All,4}]; (*Use ToExpression only for a list of Rules,Not Strings*)


(* ::Text::Initialization:: *)
(*\[Chi] Bifurcation DIagrams for base values of \[Epsilon] and \[Kappa]*)


(* ::Input::Initialization:: *)
chifixedpointsbi=Select[bifurdata,#[[2]]==1&&#[[3]]==1&];
Achibi=Values[chifixedpointsbi[[All,4,1]]];
Mchibi=Values[chifixedpointsbi[[All,4,2]]];
Vchibi=Values[chifixedpointsbi[[All,4,3]]];
chi=chifixedpointsbi[[All,1]];



(* ::Input::Initialization:: *)
bifurchangechiAMPV=Select[chifixedpointsbi,#[[6]]=="AMPV"&];
fitdataAchiAMPV=Values[bifurchangechiAMPV[[All, 4, 1]]];
fitdatachiAMPV=bifurchangechiAMPV[[All, 1]];
fitdataAchiAMPV=MapThread[List,{fitdatachiAMPV,fitdataAchiAMPV}];


bifurchangechiAM=Select[chifixedpointsbi,#[[6]]=="AM"&];
fitdataAchiAM=Values[bifurchangechiAM[[All, 4, 1]]];
fitdatachiAM=bifurchangechiAM[[All, 1]];
fitdataAchiAM=MapThread[List,{fitdatachiAM,fitdataAchiAM}];

bifurchangechiA=Select[chifixedpointsbi,#[[6]]=="A"&];
fitdataAchiA=Values[bifurchangechiA[[All, 4, 1]]];
fitdatachiA=bifurchangechiA[[All, 1]];
fitdataAchiA=MapThread[List,{fitdatachiA,fitdataAchiA}];

bifurchangechiE=Select[chifixedpointsbi,#[[6]]=="E"&];
fitdataAchiE=Values[bifurchangechiE[[All, 4, 1]]];
fitdatachiE=bifurchangechiE[[All, 1]];
fitdataAchiE=MapThread[List,{fitdatachiE,fitdataAchiE}];


(* ::Input::Initialization:: *)
formulaAchiAMPV=D[FindFormula[fitdataAchiAMPV,\[Chi]],\[Chi]];


(* ::Input::Initialization:: *)
Achiplot=Labeled[ListLinePlot[{fitdataAchiAMPV,fitdataAchiAM,fitdataAchiA,fitdataAchiE},PlotRange->All, ScalingFunctions->{None,"SignedLog"}, PlotLegends->{"AMPV","AM","A","E"},PlotStyle->{Directive[Blue,Dashed],Directive[Purple,Dashed],Directive[Orange,Thick],Directive[Red,Thick]},  LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]],{"\[Chi]","A Fixed Points"},{{Bottom,Center},Left},RotateLabel->True,LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]]



(* ::Input::Initialization:: *)
fitdataMchiAMPV=Values[bifurchangechiAMPV[[All, 4, 2]]];
fitdataMchiAMPV=MapThread[List,{fitdatachiAMPV,fitdataMchiAMPV}];
formulaMchiAMPV=D[FindFormula[fitdataMchiAMPV,\[Chi]],\[Chi]];


fitdataMchiAM=Values[bifurchangechiAM[[All, 4, 2]]];
fitdataMchiAM=MapThread[List,{fitdatachiAM,fitdataMchiAM}];


fitdataMchiA=Values[bifurchangechiA[[All, 4, 2]]];
fitdataMchiA=MapThread[List,{fitdatachiA,fitdataMchiA}];


fitdataMchiE=Values[bifurchangechiE[[All, 4, 2]]];
fitdataMchiE=MapThread[List,{fitdatachiE,fitdataMchiE}];


(* ::Input::Initialization:: *)
Mchiplot=Labeled[ListLinePlot[{fitdataMchiAMPV,fitdataMchiAM,fitdataMchiA,fitdataMchiE},PlotRange->All, ScalingFunctions->{None,"SignedLog"}, PlotLegends->{"AMPV","AM","A","E"},PlotStyle->{Directive[Blue,Dashed],Directive[Purple,Dashed],Directive[Orange,Thick],Directive[Red,Thick]},  LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]],{"\[Chi]","M Fixed Points"},{{Bottom,Center},Left},RotateLabel->True,LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]]


(* ::Input::Initialization:: *)
fitdataVchiAMPV=Values[bifurchangechiAMPV[[All, 4, 3]]];
fitdataVchiAMPV=MapThread[List,{fitdatachiAMPV,fitdataVchiAMPV}];
formulaVchiAMPV=D[FindFormula[fitdataVchiAMPV,\[Chi]],\[Chi]];


fitdataVchiAM=Values[bifurchangechiAM[[All, 4, 3]]];
fitdataVchiAM=MapThread[List,{fitdatachiAM,fitdataVchiAM}];


fitdataVchiA=Values[bifurchangechiA[[All, 4, 3]]];
fitdataVchiA=MapThread[List,{fitdatachiA,fitdataVchiA}];


fitdataVchiE=Values[bifurchangechiE[[All, 4, 3]]];
fitdataVchiE=MapThread[List,{fitdatachiE,fitdataVchiE}];





(* ::Input::Initialization:: *)
Vchiplot=Labeled[ListLinePlot[{fitdataVchiAMPV,fitdataVchiAM,fitdataVchiA,fitdataVchiE},PlotRange->All, ScalingFunctions->{None,"SignedLog"}, PlotLegends->{"AMPV","AM","A","E"},PlotStyle->{Directive[Blue,Dashed],Directive[Purple,Dashed],Directive[Orange,Thick],Directive[Red,Thick]},  LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]],{"\[Chi]","V Fixed Points"},{{Bottom,Center},Left},RotateLabel->True,LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]]


(* ::Text::Initialization:: *)
(*\[Epsilon] Bifurcation DIagrams for base values of \[Chi] and \[Kappa]*)


(* ::Input::Initialization:: *)
epsifixedpointsbi=Select[bifurdata,#[[1]]==0.3&&#[[3]]==1&];
Aepsibi=Values[epsifixedpointsbi[[All,4,1]]];
Mepsibi=Values[epsifixedpointsbi[[All,4,2]]];
Vepsibi=Values[epsifixedpointsbi[[All,4,3]]];
epsi=epsifixedpointsbi[[All,2]];

bifurchangeepsiAMPV=Select[epsifixedpointsbi,#[[6]]=="AMPV"&];
fitdataAepsiAMPV=Values[bifurchangeepsiAMPV[[All, 4, 1]]];
fitdataepsiAMPV=bifurchangeepsiAMPV[[All, 2]];
fitdataAepsiAMPV=MapThread[List,{fitdataepsiAMPV,fitdataAepsiAMPV}];


bifurchangeepsiAM=Select[epsifixedpointsbi,#[[6]]=="AM"&];
fitdataAepsiAM=Values[bifurchangeepsiAM[[All, 4, 1]]];
fitdataepsiAM=bifurchangeepsiAM[[All, 2]];
fitdataAepsiAM=MapThread[List,{fitdataepsiAM,fitdataAepsiAM}];

bifurchangeepsiA=Select[epsifixedpointsbi,#[[6]]=="A"&];
fitdataAepsiA=Values[bifurchangeepsiA[[All, 4, 1]]];
fitdataepsiA=bifurchangeepsiA[[All, 2]];
fitdataAepsiA=MapThread[List,{fitdataepsiA,fitdataAepsiA}];

bifurchangeepsiE=Select[epsifixedpointsbi,#[[6]]=="E"&];
fitdataAepsiE=Values[bifurchangeepsiE[[All, 4, 1]]];
fitdataepsiE=bifurchangeepsiE[[All, 2]];
fitdataAepsiE=MapThread[List,{fitdataepsiE,fitdataAepsiE}];


(* ::Input::Initialization:: *)
formulaAepsiAMPV=D[FindFormula[fitdataAepsiAMPV,\[Epsilon]],\[Epsilon]];


(* ::Input::Initialization:: *)
Aepsiplot=Labeled[ListLinePlot[{fitdataAepsiAMPV,fitdataAepsiAM,fitdataAepsiA,fitdataAepsiE},PlotRange->All, ScalingFunctions->{None,"SignedLog"}, PlotLegends->{"AMPV","AM","A","E"},PlotStyle->{Directive[Blue,Dashed],Directive[Purple,Dashed],Directive[Orange,Thick],Directive[Red,Thick]},  LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]],{"\[Epsilon]","A Fixed Points"},{{Bottom,Center},Left},RotateLabel->True,LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]]


(* ::Input::Initialization:: *)
fitdataMepsiAMPV=Values[bifurchangeepsiAMPV[[All, 4, 2]]];
fitdataMepsiAMPV=MapThread[List,{fitdataepsiAMPV,fitdataMepsiAMPV}];


fitdataMepsiAM=Values[bifurchangeepsiAM[[All, 4, 2]]];
fitdataMepsiAM=MapThread[List,{fitdataepsiAM,fitdataMepsiAM}];


fitdataMepsiA=Values[bifurchangeepsiA[[All, 4, 2]]];
fitdataMepsiA=MapThread[List,{fitdataepsiA,fitdataMepsiA}];


fitdataMepsiE=Values[bifurchangeepsiE[[All, 4, 2]]];
fitdataMepsiE=MapThread[List,{fitdataepsiE,fitdataMepsiE}];

formulaMepsiAMPV=D[FindFormula[fitdataMepsiAMPV,\[Epsilon]],\[Epsilon]];

Mepsiplot=Labeled[ListLinePlot[{fitdataMepsiAMPV,fitdataMepsiAM,fitdataMepsiA,fitdataMepsiE},PlotRange->All, ScalingFunctions->{None,"SignedLog"}, PlotLegends->{"AMPV","AM","A","E"},PlotStyle->{Directive[Blue,Dashed],Directive[Purple,Dashed],Directive[Orange,Thick],Directive[Red,Thick]},  LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]],{"\[Epsilon]","M Fixed Points"},{{Bottom,Center},Left},RotateLabel->True,LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]]




(* ::Input::Initialization:: *)
fitdataVepsiAMPV=Values[bifurchangeepsiAMPV[[All, 4, 3]]];
fitdataVepsiAMPV=MapThread[List,{fitdataepsiAMPV,fitdataVepsiAMPV}];


fitdataVepsiAM=Values[bifurchangeepsiAM[[All, 4, 3]]];
fitdataVepsiAM=MapThread[List,{fitdataepsiAM,fitdataVepsiAM}];


fitdataVepsiA=Values[bifurchangeepsiA[[All, 4, 3]]];
fitdataVepsiA=MapThread[List,{fitdataepsiA,fitdataVepsiA}];


fitdataVepsiE=Values[bifurchangeepsiE[[All, 4, 3]]];
fitdataVepsiE=MapThread[List,{fitdataepsiE,fitdataVepsiE}];

formulaVepsiAMPV=D[FindFormula[fitdataVepsiAMPV,\[Epsilon]],\[Epsilon]];


Vepsiplot=Labeled[ListLinePlot[{fitdataVepsiAMPV,fitdataVepsiAM,fitdataVepsiA,fitdataVepsiE},PlotRange->All, ScalingFunctions->{None,"SignedLog"}, PlotLegends->{"AMPV","AM","A","E"},Ticks->{Automatic,{#,NumberForm@#}&/@{1000,5000000.0}},PlotStyle->{Directive[Blue,Dashed],Directive[Purple,Dashed],Directive[Orange,Thick],Directive[Red,Thick]},  LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]],{"\[Epsilon]","V Fixed Points"},{{Bottom,Center},Left},RotateLabel->True,LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]]


(* ::Text::Initialization:: *)
(*\[Kappa] Bifurcation DIagrams for base values of \[Chi] and \[Epsilon]*)


(* ::Input::Initialization:: *)
kappafixedpointsbi=Select[bifurdata,#[[1]]==0.3&&#[[2]]==1&];
Akappa=Values[kappafixedpointsbi[[All,4,1]]];
Mkappa=Values[kappafixedpointsbi[[All,4,2]]];
Vkappa=Values[kappafixedpointsbi[[All,4,3]]];
kappa=kappafixedpointsbi[[All,3]];


(* ::Input::Initialization:: *)
bifurchangekappaAMPV=Select[kappafixedpointsbi,#[[6]]=="AMPV"&];
fitdataAkappaAMPV=Values[bifurchangekappaAMPV[[All, 4, 1]]];
fitdatakappaAMPV=bifurchangekappaAMPV[[All, 3]];
fitdataAkappaAMPV=MapThread[List,{fitdatakappaAMPV,fitdataAkappaAMPV}];


bifurchangekappaAM=Select[kappafixedpointsbi,#[[6]]=="AM"&];
fitdataAkappaAM=Values[bifurchangekappaAM[[All, 4, 1]]];
fitdatakappaAM=bifurchangekappaAM[[All, 3]];
fitdataAkappaAM=MapThread[List,{fitdatakappaAM,fitdataAkappaAM}];

bifurchangekappaA=Select[kappafixedpointsbi,#[[6]]=="A"&];
fitdataAkappaA=Values[bifurchangekappaA[[All, 4, 1]]];
fitdatakappaA=bifurchangekappaA[[All, 3]];
fitdataAkappaA=MapThread[List,{fitdatakappaA,fitdataAkappaA}];

bifurchangekappaE=Select[kappafixedpointsbi,#[[6]]=="E"&];
fitdataAkappaE=Values[bifurchangekappaE[[All, 4, 1]]];
fitdatakappaE=bifurchangekappaE[[All, 3]];
fitdataAkappaE=MapThread[List,{fitdatakappaE,fitdataAkappaE}];

formulaAkappaAMPV=D[FindFormula[fitdataAkappaAMPV,\[Kappa]],\[Kappa]];

Akappaplot=Labeled[ListLinePlot[{fitdataAkappaAMPV,fitdataAkappaAM,fitdataAkappaA,fitdataAkappaE},PlotRange->All, ScalingFunctions->{None,"SignedLog"}, PlotLegends->{"AMPV","AM","A","E"},PlotStyle->{Directive[Blue,Dashed],Directive[Purple,Dashed],Directive[Orange,Thick],Directive[Red,Thick]},  LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]],{"\[Kappa]","A Fixed Points"},{{Bottom,Center},Left},RotateLabel->True,LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]]


(* ::Input::Initialization:: *)
fitdataMkappaAMPV=Values[bifurchangekappaAMPV[[All, 4, 2]]];
fitdataMkappaAMPV=MapThread[List,{fitdatakappaAMPV,fitdataMkappaAMPV}];


fitdataMkappaAM=Values[bifurchangekappaAM[[All, 4, 2]]];
fitdataMkappaAM=MapThread[List,{fitdatakappaAM,fitdataMkappaAM}];


fitdataMkappaA=Values[bifurchangekappaA[[All, 4, 2]]];
fitdataMkappaA=MapThread[List,{fitdatakappaA,fitdataMkappaA}];


fitdataMkappaE=Values[bifurchangekappaE[[All, 4, 2]]];
fitdataMkappaE=MapThread[List,{fitdatakappaE,fitdataMkappaE}];

formulaMkappaAMPV=D[FindFormula[fitdataMkappaAMPV,\[Kappa]],\[Kappa]];

Mkappaplot=Labeled[ListLinePlot[{fitdataMkappaAMPV,fitdataMkappaAM,fitdataMkappaA,fitdataMkappaE},PlotRange->All, ScalingFunctions->{None,"SignedLog"}, PlotLegends->{"AMPV","AM","A","E"},PlotStyle->{Directive[Blue,Dashed],Directive[Purple,Dashed],Directive[Orange,Thick],Directive[Red,Thick]},  LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]],{"\[Kappa]","M Fixed Points"},{{Bottom,Center},Left},RotateLabel->True,LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]]


(* ::Input::Initialization:: *)
fitdataVkappaAMPV=Values[bifurchangekappaAMPV[[All, 4, 3]]];
fitdataVkappaAMPV=MapThread[List,{fitdatakappaAMPV,fitdataVkappaAMPV}];


fitdataVkappaAM=Values[bifurchangekappaAM[[All, 4, 3]]];
fitdataVkappaAM=MapThread[List,{fitdatakappaAM,fitdataVkappaAM}];


fitdataVkappaA=Values[bifurchangekappaA[[All, 4, 3]]];
fitdataVkappaA=MapThread[List,{fitdatakappaA,fitdataVkappaA}];


fitdataVkappaE=Values[bifurchangekappaE[[All, 4, 3]]];
fitdataVkappaE=MapThread[List,{fitdatakappaE,fitdataVkappaE}];

formulaVkappaAMPV=D[FindFormula[fitdataVkappaAMPV,\[Kappa]],\[Kappa]];


Vkappaplot=Labeled[ListLinePlot[{fitdataVkappaAMPV,fitdataVkappaAM,fitdataVkappaA,fitdataVkappaE},PlotRange->All, ScalingFunctions->{None,"SignedLog"}, PlotLegends->{"AMPV","AM","A","E"},Ticks->{Automatic,{#,NumberForm@#}&/@{1000,5000000.0}},PlotStyle->{Directive[Blue,Dashed],Directive[Purple,Dashed],Directive[Orange,Thick],Directive[Red,Thick]},  LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]],{"\[Kappa]","V Fixed Points"},{{Bottom,Center},Left},RotateLabel->True,LabelStyle->Directive[FontSize->14,FontFamily->"Helvetica",Black]]


(* ::Text::Initialization:: *)
(*3x3 Bifurcation Diagrams Plot*)


(* ::Input::Initialization:: *)
GraphicsGrid[{{ Achiplot, Aepsiplot, Akappaplot},{ Mchiplot, Mepsiplot, Mkappaplot},{ Vchiplot, Vepsiplot,Vkappaplot}},ImageSize->1000,AspectRatio->Full]



(* ::Text::Initialization:: *)
(*2x2 Bifurcation Diagrams Plot*)


(* ::Input::Initialization:: *)
GraphicsGrid[{{ Achiplot, Mchiplot},{ Vchiplot, Aepsiplot}},ImageSize->1000,AspectRatio->Full]


(* ::Input::Initialization:: *)
Export["/home/kai-hua/OneDrive/Article pictures/bifurcationdiagram.tiff",%744,"TIFF"]


(* ::Text::Initialization:: *)
(*Derivative Sensitivity Analysis*)


(* ::Input::Initialization:: *)
Grid[{{"AMPV","\[Chi] Derivate "},{"A",formulaAchiAMPV},{"M",formulaMchiAMPV},{"V",formulaVchiAMPV}},Frame->All,Alignment->{Left,Center}]
