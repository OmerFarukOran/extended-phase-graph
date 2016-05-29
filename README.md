# extended-phase-graph

This MATLAB program (phaseDiagram.m) is for the generation and the investigation of the Extended Phase Graph (EPG) for a Magnetic Resonance Imaging (MRI) multi-spin-echo pulse sequence. The program is developed based on the paper "Extended phase graphs: Dephasing, RF pulses, and echoes - pure and simple" by M. Weigel (http://doi.org/10.1002/jmri.24619).

Given the parameters such as flip angle, initial phase and time-position of the RF pulses, the complex magnetization is calculated for each of the dephasing, rephasing and stored pathways. Each pathway is then drawn on the phase diagram.

The program also allows to specify the relative crusher gradient strength for each of the refocusing RF pulse. By this way, the contribution of each echo to the final signal can be tailored.

Here is a snapshot of the figure obtained using the program

![Snapshot](/snapshot.jpg)

Omer Faruk Oran, May 2016
