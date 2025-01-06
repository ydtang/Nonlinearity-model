
NEURON {
	:NOTE: since this is an interface, there is no POINT_PROCESS name
	RANGE e, i, gmax, g
	RANGE tau_r
	RANGE tau_d
	RANGE dist :distance from the soma
	RANGE ndist :normalized distance from the soma
	RANGE factor
	RANGE tsyn,mg,mggate,gamma
	GLOBAL mg_dep
	
	RANGE lSF,gSF,Vtrg	:local and global synaptic homeostatic plasticity factor between 0-1
	RANGE tau_H	: time constant for the homeostasis, should be orders of magnitude larger than the membrane time constant.

	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	e = 0			(mV)	
	tau_r = 0	 	(ms) 	
	tau_d = 0		(ms) 
	dist
	factor
	gmax
	tsyn
	mg			(mM)
	mggate
	gamma = 0.062
	Vtrg = -60
	:the default is that there is not homeostasis with such a high time constant
	tau_H = 1e300	 	(ms) :assuming membrane potential within the range of 20ms tau_H is three orders of magnitude higher. 	
	mg_dep = 3.57
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
}

STATE {
	A
	B
	lSF
}

FUNCTION get_factor(tau_fast,tau_slow){
	LOCAL tp
	tp = (tau_fast*tau_slow)/(tau_slow - tau_fast) * log(tau_slow/tau_fast)
	get_factor = -exp(-tp/tau_fast) + exp(-tp/tau_slow)
	get_factor = 1/get_factor
}

INITIAL {
	:start fresh
	A = 0
	B = 0
	g = 0
	lSF = 1
	tsyn = -1e80 			:last time the synapse was active is very long ego
	if(tau_r/tau_d>0.999){
		tau_r = tau_d*0.999
	}
	if(!tau_r || !tau_d){
		printf("User must set tau_r and tau_d (zero by default)")
	}
	factor = get_factor(tau_r,tau_d)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	:see Jahr & Stevens, J. Neurosci 10: 1830-1837, 1990;
	:Jahr & Stevens, J. Neurosci 10: 3178-3182, 1990
	: from Jahr & Stevens
	mggate = 1 / (1 + exp( gamma *  -v) * (mg / mg_dep))
	lSF = lSF*(lSF>0)
	g = gmax*(B - A)*mggate*lSF
	i = g*(v - e)
}


DERIVATIVE state {
	A' = -A/tau_r
	B' = -B/tau_d
	lSF' = (Vtrg-v)/tau_H :the global implementation requires a pointer to the somatic voltage.
}


PROCEDURE EPSP(w){
	A = A + w*factor
	B = B + w*factor
	tsyn = t
}
