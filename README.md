#  Testing...

<b>R code (rcode folder) and input files (data folder) for analyses of <i>Calanus hyperboreus</i> distribution and energy requirements.</b>

The <b>rcode</b> folder contains the following R scripts:<br>
<i>gams_chyp.r</i> : Fit GAMs of Chyp abundance data<br> 
<i>c3_seasonal.r</i> : Plot seasonal variation in Chyp C3 concentrations against growth season<br>
<i>energy_budget.r</i> : Estimate Chyp energy requirements to develop from N3-C3<br>
... plus additional files to plot polar-centered maps, draw color image scales and calculate leave-one-year-out cross validation<br>


The <b>data</b> folder contains the following R data files:<br>
<i>chyp_orig.rda</i>: Original Chyp abundance data in ind.m2 or ind.m3 (see <i>c3_seasonal.r</i>)<br>
<i>chyp_depthint.rda</i>: Standardized depth integrated Chyp abundance data in ind.m2 (see <i>gams_chyp.r</i>)<br>
<i>carb_weight_stages.rda</i>: Stage-specific Chyp body weights (see <i>energy_budget.r</i>)<br>
<i>filtrations.rda</i>: Stage-specific Chyp filtration rates (see <i>energy_budget.r</i>)<br>
<i>phyto_daymean.rda</i>: Averaged daily phytoplankton/phyto+microzooplankton concentrations (see <i>energy_budget.r</i>/<i>c3_seasonal.r</i>)<br>
<i>predframe.rda</i>: Dataframe for GAM predictions (see <i>gams_chyp.r</i>)<br> 
As well as <i>chyp_data_sources.xlsx</i>, a full list of data sources for Chyp abundance data.