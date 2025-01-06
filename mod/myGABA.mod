
COMMENT


Use only what is in ISyn and IUtils

ENDCOMMENT

NEURON {
	POINT_PROCESS myGABA
}

INCLUDE "ISyn.inc"
INCLUDE "IUtils.inc"

NET_RECEIVE(w) {
	EPSP(1): finally trigger an epsp
}
