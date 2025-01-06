


NEURON {
	SUFFIX Linear
        RANGE Nonlinear,ratio
}
PARAMETER {
	Nonlinear
	ratio
}

INITIAL {
Nonlinear=1
}

BREAKPOINT {
Nonlinear=ratio	
}
