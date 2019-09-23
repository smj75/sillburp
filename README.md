# sillburp
**Greenhouse gas emissions from igneous sills.**

Sillburp calculates the temperature and gas production history when an igneous sheet intrudes sedimentary rock.  The primary reference is, "Jones, S. M., Hoggett, M., Greene, S. E. & Dunkley Jones, T.,  Large Igneous Province thermogenic greenhouse gas flux could have initiated Paleocene-Eocene Thermal Maximum climate change, Nature Communications, in press Sept 2019," referred to below as Jea19; see the Methods section.

Make using "make".  Type "./sillburp" to see simple instructions.  Run GMT scripts in the order described to process model results.

## *sillburp.cc*, *sillburp.h*

Source code for sillburp.

## *burp.cc* , *burp.h*

Additional source code common with LIPburp code.

## *run_sillburp_jea19.sh*

Example script that runs *sillburp* to produce greenhouse gas emissions histories for sills of various thicknesses and depths.  These were the runs used in Jea19.  These calculations were done in separate directories for each depth of emplacement.

## *parameterize_flux_[D]km.gmt*

Scripts to fit tapered power law parameterisation to sillburp emissions results.  [D] is the emplacement depth.  These calculations were done in separate directories for each depth of emplacement.

## *collate_aureole.gmt*, *collate_decay_time.gmt*, *collate_p.gmt*

Scripts to collate, contour and plot the 3 power law parameters.  Also writes out the collated values in a format that can be used by *LIPburp/dim2mtp.awk*.

## *sill_labile_cracking_aureoles_direct.xyz*, *sill_refractory_aureoles_direct.xyz*, *dless_decay_time_oil.xyz*, *dless_decay_time_ref.xyz*, *p_oil.xyz*, *p_ref.xyz*

Collated results of dimensionless aureole thickness (Y), dimensionless decay time (tau_decay) and power law exponent (p) used to plot Jea19 Supplmentary Fig. 1.

