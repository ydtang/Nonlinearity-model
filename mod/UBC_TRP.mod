TITLE Cerebellum Granule Cell Model leakage

COMMENT
	Reference: Theta-Frequency Bursting and Resonance in Cerebellar Granule Cells:Experimental
	Evidence and Modeling of a Slow K+-Dependent Mechanism
	Egidio D'Angelo,Thierry Nieus,Arianna Maffei,Simona Armano,Paola Rossi,Vanni Taglietti,
	Andrea Fontana and Giovanni Naldi
	
	SAT: its a distributed mechanism with out any dynamic behaviour. Its just modified from Granular cell's leakage channel.
ENDCOMMENT
 
NEURON { 
	SUFFIX Ubc_TRP 
	NONSPECIFIC_CURRENT itrp
	USEION cAMP READ cAMPi VALENCE 0
	RANGE etrp, gtrp, i, itrp, fcAMP, TonicTRP
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	v (mV) 
	gtrp = 4.18e-6 (mho/cm2)
	celsius = 30 (degC)
	etrp = 0 (mV)
	fcAMP = 1
	theta = 0 (1)
	TonicTRP = 0.05 
    } 
    
    ASSIGNED { 
	cAMPi  (mM)
	itrp (mA/cm2) 
	i (mA/cm2)
    }
    
    BREAKPOINT { 
        itrp = gtrp*(v - etrp) * (TonicTRP + cAMPi*fcAMP )
	i = itrp
} 
