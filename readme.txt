The zip file contains a Matlab and an Excel implementation of the
opponent-channel model described in:

Briley PM, Kitterick PT, Summerfield AQ (2012) Evidence for
opponent-process analysis of sound-source location in humans. J Assoc
Res Otolaryngol. Epub Oct 23rd, 2012. PMID: 23090057. doi:
10.1007/s10162-012-0356-x

Both implementations reproduce the key elements of Figures 8B-C, 9
and 11 of the associated paper. Specific details of model parameters
and outputs are given at the start of the Matlab file and on the
"Description" sheet of the Excel workbook. In an "opponent-channel"
model, sound-source location is coded by the balance of activity in
two, broadly-tuned channels, one maximally responsive to the left, and
one to the right, auditory hemispace. Our model can be used to
describe the tuning of these two channels, and consequently to predict
psychophysical spatial acuity ("minimum audible angles") for different
azimuthal positions, as well as electroencephalographic (EEG)
responses to abrupt shifts in sound-source location. The default
parameters were obtained by fitting the model to EEG data from a group
of young, normally-hearing adults.

In the Excel implementation of the model, parameters can be set, and
output can be viewed, in the "Figures" sheet. The detailed
computations are performed in the "Model" sheet. The Matlab
implementation is run using the following syntax:

[chan maa ac] = brileyetal2012_oppchanmodel(figs,params)

Inputs (optional):

figs - set to 0 to suppress the production of figures (default is 1)
params - a structure that can be used to set the model parameters
(e.g. the inclusion of params.scale would set the scale parameter to
the specified value), otherwise default parameters are used

Outputs:

chan - contains the channel parameters and tuning functions
maa - contains the predicted minimum audible angles
ac - contains the predicted EEG response sizes

These model files were contributed by Paul Briley.
