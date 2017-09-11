# Spin-dependent-Interaction

This program provides the sensitivity estimation of this paper https://arxiv.org/abs/1606.01152.

spindependentinteraction [force index] [sample position]
	Calculate the coupling strength as a function of interaction length 
	[force index] : different spin-dependent interactions: (2,3,4,6,8,9,11,12,14,15,16)
	[position] : there are two sample positions: 1. the distance is gap+modulation amplitude.
		     	       	   	  	     2. the distance is gap.

plot [force index] [sample position] [max] [min]
     plot coupling strength as a function of interaction length
     [max] : maximum value of the plot
     [min] : minimum value of the plot

submit.pl :
	This script can calculate all forces.