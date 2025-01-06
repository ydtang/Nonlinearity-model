: File to track maximum value calcium concentration reached.

NEURON{
	SUFFIX camax
	USEION ca READ cai
	RANGE camx, tmx
}

UNITS {
 	(molar) = (1/liter)
  	(mM)    = (millimolar)
}


ASSIGNED{
	camx (mM)
	tmx (ms)
	cai (mM)
}

INITIAL{
	camx=cai
}

BREAKPOINT{
	if(cai>camx){
		camx=cai
		tmx=t
	}
}
