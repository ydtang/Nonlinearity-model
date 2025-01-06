COMMENT by Johannes Luthman: 
Based on NEURON 6.0's built-in exp2syn.mod.
Changes made to the original: 
* tau1 renamed taurise; tau2, taudecay
* restructuring of NEURON block
* microsiemens changed to siemens for consistency with the other NMODLs.


Original comment: 
Two state kinetic scheme synapse described by rise time taurise,
and decay time constant taudecay. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/taurise and 1/taudecay is
 A = a*exp(-t/taurise) and
 G = a*taudecay/(taudecay-taurise)*(-exp(-t/taurise) + exp(-t/taudecay))
	where taurise < taudecay

If taudecay-taurise -> 0 then we have a alphasynapse.
and if taurise -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS cfAMPA
	NONSPECIFIC_CURRENT i
	RANGE g, i, e, taurise, taudecay,del,v
	RANGE gampamax
}

UNITS {
	(nA) = (nanoamp)
	(nS) = (nanomho)
	(mV) = (millivolt)
}

PARAMETER {
	gampamax=0.8 (ns)
	taurise = 1 (ms)
	taudecay = 1 (ms)
	e = 0 (mV)
	del=0 (ms)
}

ASSIGNED {
	v (mV)
	i (nA)
	g (microsiemens)
	factor
}

STATE {
	A (nanomho)
	B (nanomho)
}

INITIAL {
	LOCAL tp
	if (taurise/taudecay > .9999) {
		taurise = .9999*taudecay
	}
	A = 0
	B = 0
	tp = (taurise*taudecay)/(taudecay - taurise) * log(taudecay/taurise)
	factor = -exp(-tp/taurise) + exp(-tp/taudecay)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
	:printf("%gww\n",i)
}

DERIVATIVE state {
	A' = -A/taurise
	B' = -B/taudecay
}

NET_RECEIVE(weight (nanomho)) {
	
	if(t>del){
		:printf("%g",t)
		state_discontinuity(A, A + weight*factor*gampamax)
		state_discontinuity(B, B + weight*factor*gampamax)
	}
}
