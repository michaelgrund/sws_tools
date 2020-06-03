% Generate different plots for synthetic splitting distributions for a
% dipping layer with different dip angles.

% see functions <<< split_model >>> and  <<< SWS_modelling_stereoplot_THEO_dipping >>>

%%

incangle=11

dips=0
split_model(incangle,dips)

SWS_modelling_stereoplot_THEO_dipping(dips)

dips=20
split_model(incangle,dips)

SWS_modelling_stereoplot_THEO_dipping(dips)

dips=40
split_model(incangle,dips)

SWS_modelling_stereoplot_THEO_dipping(dips)

dips=60
split_model(incangle,dips)

SWS_modelling_stereoplot_THEO_dipping(dips)

dips=70
split_model(incangle,dips)

SWS_modelling_stereoplot_THEO_dipping(dips)


dips=90
split_model(15,dips)

SWS_modelling_stereoplot_THEO_dipping(dips)