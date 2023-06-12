:Tone to Interneuron Cells AMPA+NMDA with local Ca2+ pool
:NON-Plastic

NEURON {
	POINT_PROCESS tone2vip
	USEION ca READ eca	
	NONSPECIFIC_CURRENT inmda, iampa
	RANGE initW
	RANGE Cdur_nmda, AlphaTmax_nmda, Beta_nmda, Erev_nmda, gbar_nmda, W_nmda, on_nmda, g_nmda
	RANGE Cdur_ampa, AlphaTmax_ampa, Beta_ampa, Erev_ampa, gbar_ampa, W, on_ampa, g_ampa
	RANGE eca, ICan, P0n, fCan, tauCa, Icatotal
	RANGE ICaa, P0a, fCaa
	RANGE Cainf, pooldiam, z
	RANGE lambda1, lambda2, threshold1, threshold2
	RANGE fmax, fmin, Wmax, Wmin, maxChange, normW, scaleW, srcid, destid
	RANGE pregid,postgid
	THREADSAFE
}

UNITS {
	(mV) = (millivolt)
        (nA) = (nanoamp)
	(uS) = (microsiemens)
	FARADAY = 96485 (coul)
	PI = (pi) (1)
}

PARAMETER {

	srcid = -1 (1)
	destid = -1 (1)

	Cdur_nmda = 16.7650 (ms)
	AlphaTmax_nmda = .2659 (/ms)
	Beta_nmda = 0.008 (/ms)
	Erev_nmda = 0 (mV)
	gbar_nmda = .5e-3 (uS)

	Cdur_ampa = 0.713 (ms)
	AlphaTmax_ampa = 10.1571 (/ms)
	Beta_ampa = 0.4167 (/ms)
	Erev_ampa = 0 (mV)
	gbar_ampa = 1e-3 (uS)

	eca = 120

	Cainf = 50e-6 (mM)
	pooldiam =  1.8172 (micrometer)
	z = 2

	tauCa = 50 (ms)
	P0n = .015
	fCan = .024
	
	P0a = .0015 : .001
	fCaa = .024
	
	lambda1 = 3 : 4 : 3 : 4 : 5 : 7 : 10 
	lambda2 = 0.01 
	threshold1 = 0.45 : 0.5 : 0.55 : .4 : 0.35
	threshold2 = 0.50 : 0.55 : 0.6 : 0.45 : 0.4
	
	initW = 4.5 : 6 : 8 : 7.5 : 6 : 5 : 4 : 5 : 4 : 3 : 4	
	fmax = 4.3 : 3.8 : 4 : 3
	fmin = .8		
	
	GAPstart1 = 96000
	GAPstop1 = 196000		
	
}

ASSIGNED {
	v (mV)

	inmda (nA)
	g_nmda (uS)
	on_nmda
	W_nmda

	iampa (nA)
	g_ampa (uS)
	on_ampa
	W

	t0 (ms)

	ICan (mA)
	ICaa (mA)
	Afactor	(mM/ms/nA)
	Icatotal (mA)

	dW_ampa
	Wmax
	Wmin
	maxChange
	normW
	scaleW
	
	pregid
	postgid
}

STATE { r_nmda r_ampa capoolcon }

INITIAL {
	on_nmda = 0
	r_nmda = 0
	W_nmda = initW

	on_ampa = 0
	r_ampa = 0
	W = initW

	t0 = -1

	Wmax = fmax*initW
	Wmin = fmin*initW
	maxChange = (Wmax-Wmin)/10
	dW_ampa = 0

	capoolcon = Cainf
	Afactor	= 1/(z*FARADAY*4/3*PI*(pooldiam/2)^3)*(1e6)
}

BREAKPOINT {
	SOLVE release METHOD cnexp
}

DERIVATIVE release {
	if (t0>0) {
		if (t-t0 < Cdur_nmda) {
			on_nmda = 1
		} else {
			on_nmda = 0
		}
		if (t-t0 < Cdur_ampa) {
			on_ampa = 1
		} else {
			on_ampa = 0
		}
	}
	r_nmda' = AlphaTmax_nmda*on_nmda*(1-r_nmda)-Beta_nmda*r_nmda
	r_ampa' = AlphaTmax_ampa*on_ampa*(1-r_ampa)-Beta_ampa*r_ampa

	dW_ampa = eta(capoolcon)*(lambda1*omega(capoolcon, threshold1, threshold2)-lambda2*GAP1(GAPstart1, GAPstop1)*W)*dt

	: Limit for extreme large weight changes
	if (fabs(dW_ampa) > maxChange) {
		if (dW_ampa < 0) {
			dW_ampa = -1*maxChange
		} else {
			dW_ampa = maxChange
		}
	}

	:Normalize the weight change
	:normW = (W-Wmin)/(Wmax-Wmin)
	:if (dW_ampa < 0) {
	:	scaleW = sqrt(fabs(normW))
	:} else {
	:	scaleW = sqrt(fabs(1.0-normW))
	:}

	:W = W + dW_ampa*scaleW
	
	:Weight value limits
	:if (W > Wmax) { 
	:	W = Wmax
	:} else if (W < Wmin) {
 	:	W = Wmin
	:}

	g_nmda = gbar_nmda*r_nmda
	inmda = W_nmda*g_nmda*(v - Erev_nmda)*sfunc(v)

	g_ampa = gbar_ampa*r_ampa
	iampa = W*g_ampa*(v - Erev_ampa)

	ICan = P0n*g_nmda*(v - eca)*sfunc(v)
	ICaa = P0a*W*g_ampa*(v-eca)/initW	
	Icatotal = ICan + ICaa
	capoolcon'= -fCan*Afactor*Icatotal + (Cainf-capoolcon)/tauCa
}

NET_RECEIVE(dummy_weight) {
	t0 = t
}

:::::::::::: FUNCTIONs and PROCEDUREs ::::::::::::

FUNCTION sfunc (v (mV)) {
	UNITSOFF
	sfunc = 1/(1+0.33*exp(-0.06*v))
	UNITSON
}

FUNCTION eta(Cani (mM)) {
	LOCAL taulearn, P1, P2, P4, Cacon
	P1 = 0.1
	P2 = P1*1e-4
	P4 = 1
	Cacon = Cani*1e3
	taulearn = P1/(P2+Cacon*Cacon*Cacon)+P4
	eta = 1/taulearn*0.001
}

FUNCTION omega(Cani (mM), threshold1 (uM), threshold2 (uM)) {
	LOCAL r, mid, Cacon
	Cacon = Cani*1e3
	r = (threshold2-threshold1)/2
	mid = (threshold1+threshold2)/2
	if (Cacon <= threshold1) { omega = 0}
	else if (Cacon >= threshold2) {	omega = 1/(1+50*exp(-50*(Cacon-threshold2)))}
	else {omega = -sqrt(r*r-(Cacon-mid)*(Cacon-mid))}
}
FUNCTION GAP1(GAPstart1 (ms), GAPstop1 (ms)) {
	LOCAL s
	if (t <= GAPstart1) { GAP1 = 1}
	else if (t >= GAPstart1 && t <= GAPstop1) {GAP1 = 1}					: During the Gap, apply lamda2*2
	else  {	GAP1 = 1}
}





