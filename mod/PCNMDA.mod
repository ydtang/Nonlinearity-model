

NEURON {
	POINT_PROCESS PCNMDA
	NONSPECIFIC_CURRENT i
	RANGE  taurise_0, a1, taudecay_0, a2, tauV, e, i, gVI, st_gVD, v0_gVD, Mg, K0, delta, wf,g
	RANGE taurise, taudecay,del,gnmdamax
	GLOBAL inf
	THREADSAFE
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(nS) 	= (nanomho)
	(mM) = (milli/liter)
	(S)  = (siemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(J)  = (joules)
}

PARAMETER {
: Parameters Control Neurotransmitter and Voltage-dependent gating of NMDAR
	: parameters control exponential decay of taurise
	taurise=8 (ms)
	taudecay=30 (ms)
	gnmdamax=0.96  (ns)
	taurise_0 = 2.2340 (ms)
	del=0 (ms)
	a1 = 0.0896		(ms)
	b1 = 0.0324		(1/mV)
	
	: parameters control exponential rise to a maximum
	taudecay_0 = 62.5808 (ms)
	a2 = 10.0374	 (ms)
	b2 = 0.0239		 (1/mV)

: Parameters Control voltage-dependent gating of NMDAR
	tauV = 7		(ms)	<1e-9,1e9>	: Kim11 
							: at 26 degC & [Mg]o = 1 mM, 
							: [Mg]o = 0 reduces value of this parameter
							: Because TauV at room temperature (20) & [Mg]o = 1 mM is 9.12 Clarke08 & Kim11 
							: and because Q10 at 26 degC is 1.52
							: then tauV at 26 degC should be 7 
	st_gVD = 0.007	(1/mV)	: steepness of the gVD-V graph from Clarke08 -> 2 units / 285 mv
	v0_gVD = -100	(mV)	: Membrane potential at which there is no voltage dependent current, from Clarke08 -> -90 or -100
	gVI = 1			(uS)	: Maximum Conductance of Voltage Independent component, This value is used to calculate gVD
	Q10 = 1.52				: Kim11
	T0 = 26			(degC)	: reference temperature 
	celsius 		(degC)	: actual temperature for simulation, defined in Neuron, usually about 35
: Parameters Control Mg block of NMDAR
	Mg = 1			(mM)	: external magnesium concentration from Spruston95
	K0 = 4.1		(mM)	: IC50 at 0 mV from Spruston95
	delta = 0.8 	(1)		: the electrical distance of the Mg2+ binding site from the outside of the membrane from Spruston95
: Parameter Controls Ohm's law in NMDAR
	e = -0.7		(mV)	: in CA1-CA3 region = -0.7 from Spruston95
}

CONSTANT {
	T = 273.16	(degC)
	F = 9.648e4	(coul)	: Faraday's constant (coulombs/mol)
	R = 8.315	(J/degC): universal gas constant (joules/mol/K)
	z = 2		(1)		: valency of Mg2+
}

ASSIGNED {
	v		(mV)
	dt		(ms)
	i		(nA)
	factor
	wf
	g       (us)
	inf		(uS)
	tau		(ms)
	
}

STATE {
	A (ns)
	B (ns)
	C
	gVD (uS)
	
	}

INITIAL {
	LOCAL tp
	rates(v)
	if (taurise/taudecay > .9999) {
		taurise = .9999*taudecay
	}
	A = 0
	B = 0
	tp = (taurise*taudecay)/(taudecay - taurise) * log(taudecay/taurise)
	factor = -exp(-tp/taurise) + exp(-tp/taudecay)
	factor = 1/factor
	
	: temperature-sensitivity of the slow unblock of NMDARs
	tau = tauV * Q10^((T0 - celsius)/10(degC))

	gVD = 0
	wf = 1
	Mgblock(v)
	
}

BREAKPOINT {
	SOLVE state METHOD runge : cnexp
	g=(B - A)*(gVI + gVD)*Mgblock(v)
	i = (B - A)*(gVI + gVD)*Mgblock(v)*(v - e)	
	:printf("%gss\n",i)
}

DERIVATIVE state {
	rates(v)
	A' = -A/taurise
	B' = -B/taudecay
	: Voltage Dapaendent Gating of NMDA needs prior binding to Glutamate Kim11
	gVD' = (B/wf)*(inf-gVD)/tau		
}

NET_RECEIVE(weight(nanomho)) {
	if(t>del){
		wf = weight*factor*gnmdamax
		A = A + wf 	
		B = B + wf
	}	
}

FUNCTION Mgblock(v(mV)) {
	: from Spruston95
	Mgblock = 1 / (1 + (Mg/K0)*exp((0.001)*(-z)*delta*F*v/R/(T+celsius)))
}

PROCEDURE rates(v (mV)) { 
		
	inf = (v - v0_gVD) * st_gVD * gVI
	
}