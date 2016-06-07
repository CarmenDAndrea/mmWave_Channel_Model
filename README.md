On Clustered Statistical MIMO Millimeter Wave Channel Simulation
==========

This code package is related to the following article:

S.Buzzi, C.D'Andrea,"[On Clustered Statistical MIMO Millimeter Wave Channel Simulation](http://arxiv.org/pdf/1604.00648.pdf)". 


## Abstract of the Article

The use of mmWave frequencies is one of the key strategies to achieve the fascinating 1000x increase in the capacity of future 5G wireless systems. While for traditional sub-6 GHz cellular frequencies several well-developed statistical channel models are available for system simulation, similar tools are not available for mmWave frequencies, thus preventing a fair comparison of independently developed transmission and reception schemes. In this paper we provide a simple albeit accurate statistical procedure for the generation of a clustered MIMO channel model operating at mmWaves, for both the cases of slowly and rapidly time-varying channels. Matlab scripts for channel generation are also provided, along with an example of their use.


## Content of Code Package

The package contains 2 Matlab functions: Generate_Channel_frequency_selective_LTI.m and Generate_Channel_frequency_selective_LTV.m for the generation of channel impulse responses according to expression (1) and (6) in article listed above. In the package there are also 3 ausiliary functions Array_response.m, Evaluation_Path_loss.m and Laplace_distribution.m necessary for the generation of channel's impulse responses. In the script Main_Channel.m there is an example of use of these functions.

The paper listed above contains 2 simulation figures, Figure 2 (a) and Figure 2 (b), generated respectively by Matlab script Main_Figure2a.m and Main_Figure2b.m. The auxiliary function for these scripts are Spectral_Efficiency.m and Empirical_CDF.m.
In these figures we report the CDFs of spectral efficiency of MIMO mm-Wave channel using single carrier modulation, in different conditions of number of antennas and for different number of symbols transmitted simultaneously on MIMO channel. 


##License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
