""" vectfit3.py     => This script uses vector fitting algorithm for fitting a given response of a network.

Author: Simon De Ridder
Date: 02/08/2019

Acknowledgement: This script is a translation of vector fitting algorithm written by
                 Bjorn Gustavsen in MATLAB (http://www.sintef.no/Projectweb/VECTFIT/).
		 

 [1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency
     domain responses by Vector Fitting", IEEE Trans. Power Delivery,
     vol. 14, no. 3, pp. 1052-1061, July 1999.

 [2] B. Gustavsen, "Improving the pole relocating properties of vector
     fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592,
     July 2006.

 [3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter,
     "Macromodeling of Multiport Systems Using a Fast Implementation of
     the Vector Fitting Method", IEEE Microwave and Wireless Components
     Letters, vol. 18, no. 6, pp. 383-385, June 2008.
     
 [4] Simon De Ridder, GitHub. Feb 08, 2019. Accessed on: August 10, 2022, [Online].https://github.com/SimonDeRidder/pyVF
 
"""

## Input: f, and s matrices, and number of pole pairs           Output: poles and residues matrices


import numpy as np
nax = np.newaxis

def vectfit3(f, s, poles, weight, opts=None):
	'''
	PURPOSE: Approximate f(s) with a state-space model
	             f(s)=C*(s*I-A)^(-1)*B +D +s*E
	         where f(s) is a singe element or a vector of elements.
	         When f(s) is a vector, all elements become fitted with a common pole set.

	INPUT:
	f(s): function (vector) to be fitted.
	      dimension : (Nc,Ns)
	          Nc: number of elements in vector
	          Ns: number of frequency samples

	s: vector of frequency points [rad/sec]
	       dimension : (Ns)

	poles: vector of initial poles [rad/sec]
	       dimension : (N)

	weight: the rows in the system matrix are weighted using this array. Can be used  for achieving
	        higher accuracy at desired frequency samples.
	        If no weighting is desired, use unitary weights: weight=np.ones((1,Ns)).
	        Two dimensions are allowed:
	            dimension : (1,Ns) --> Common weighting for all vector elements.
	            dimension : (Nc,Ns)--> Individual weighting for vector elements.

	opts: dict containing fit options:
		opts['relax']: True:  Use relaxed nontriviality constraint
		               False: Use nontriviality constraint of "standard" vector fitting

		opts['stable']: True:  unstable poles are kept unchanged
		                False: unstable poles are made stable by 'flipping' them into the left
		                       half-plane

		opts['asymp']: 0: Fitting with D=0,  E=0
		               1: Fitting with D!=0, E=0
		               2: Fitting with D!=0, E~=0

		opts['skip_pole']: True: The pole identification part is skipped, i.e (C,D,E)  are
		                         identified using the initial poles (A) as final poles.

		opts['skip_res']: True: The residue identification part is skipped, i.e. only the poles (A)
		                        are identified while C,D,E are returned as zero.

		opts['cmplx_ss']: True:  The returned state-space model has real and complex conjugate
		                         parameters. Output variable A is diagonal (and sparse). 
		                  False: The returned state-space model has real parameters only. Output
		                         variable A is square with 2x2 blocks (and sparse).

		opts['spy1']: True: Plotting, after pole identification (A)
		                    magnitude functions
		                    cyan trace:  (sigma*f)fit
		                    red trace:   (sigma)fit
		                    green trace: f*(sigma)fit - (sigma*f)fit

		opts['spy2']: True: Plotting, after residue identification (C,D,E)
		                    1) magnitude functions
		                    2) phase angles

		opts['logx']: True: Plotting using logarithmic absissa axis

		opts['logy']: True: Plotting using logarithmic ordinate axis

		opts['errplot']: True: Include deviation in magnitude plot

		opts['phaseplot']: True: Show plot also for phase angle

		opts['legend']: True: Include legend in plots

		opts['block']: True: Block on plots (must be closed before computation continues)

	OUTPUT :

		fit(s) = C*(s*I-A)^(-1)*B + D + s*E

	SER['A'] (N,N): A-matrix (sparse). If opts['cmplx_ss']==True: Diagonal and complex. Otherwise,
	                square and real with 2x2 blocks.
	SER['B'] (N): B-matrix. If opts['cmplx_ss']==True:  Column of 1's. 
	                        If opts['cmplx_ss']==False: contains 0's, 1's and 2's
	SER['C'] (Nc,N): C-matrix. If opts['cmplx_ss']==True:  complex
	                           If opts['cmplx_ss']==False: real-only
	SER['D'] (Nc): constant term (real). Is not None if opts['asymp']==1 or 2.
	SER['E'] (Nc): proportional term (real). Is not None if opts['asymp']==2.

	poles (N): new poles

	rmserr (1): root-mean-square error of approximation for f(s).
	            (0 is returned if opts['skip_res']==True)

	fit(Nc,Ns): Rational approximation at samples. (None is returned if opts['skip_res']==True).

	APPROACH:
	The identification is done using the pole relocating method known as Vector Fitting [1],
	with relaxed non-triviality constraint for faster convergence and smaller fitting errors [2],
	and utilization of matrix structure for fast solution of the pole identifion step [3].

	***********************************************************************************************
	NOTE: The use of this program is limited to NON-COMMERCIAL usage only.
	If the program code (or a modified version) is used in a scientific work,
	then reference should be made to the following

	[1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain responses by
	    Vector Fitting", IEEE Trans. Power Delivery, vol. 14, no. 3, pp. 1052-1061, July 1999.

	[2] B. Gustavsen, "Improving the pole relocating properties of vector fitting",
	    IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592, July 2006.

	[3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter, "Macromodeling of Multiport
	    Systems Using a Fast Implementation of the Vector Fitting Method",
	    IEEE Microwave and Wireless Components Letters, vol. 18, no. 6, pp. 383-385, June 2008.
	***********************************************************************************************
	Last revised: 08.08.2008.
	Created by:   Bjorn Gustavsen.
	Translated by: Simon De Ridder                                                              '''

	# default options
	options = {}
	options['relax']     = True  # Use vector fitting with relaxed non-triviality constraint
	options['stable']    = True  # Enforce stable poles
	options['asymp']     = 1     # Include only D in fitting (not E)
	options['skip_pole'] = False # Do NOT skip pole identification
	options['skip_res']  = False # Do NOT skip identification of residues (C,D,E)
	options['cmplx_ss']  = True  # Create complex state space model
	options['spy1']      = False # No plotting for first stage of vector fitting
	options['spy2']      = False # Do NOT Create magnitude plot for fitting of f(s)
	options['logx']      = True  # Use logarithmic abscissa axis
	options['logy']      = True  # Use logarithmic ordinate axis
	options['errplot']   = True  # Include deviation in magnitude plot
	options['phaseplot'] = False # exclude plot of phase angle (in addition to magnitiude)
	options['legend']    = True  # Do include legends in plots
	options['block']     = False # Do NOT block on plots
	if not opts is None:
		# Merge default values into opts
		options.update(opts);

	# Tolerances used by relaxed version of vector fitting
	TOLlow  = 1e-18
	TOLhigh = 1e18

	# check options
	if (not np.isscalar(options['relax'])) or (not isinstance(options['relax'], bool)):
		raise TypeError("options['relax'] must be a scalar boolean")
	if (not np.isscalar(options['stable'])) or (not isinstance(options['stable'], bool)):
		raise TypeError("options['stable'] must be a scalar boolean")
	if (not np.isscalar(options['asymp'])) or (not isinstance(options['asymp'], int)):
		raise TypeError("options['asymp'] must be a scalar integer")
	if options['asymp']<0 or options['asymp']>2:
		raise ValueError("options['asymp'] must be 0, 1 or 2, not "+str(options['asymp']))
	if (not np.isscalar(options['skip_pole'])) or (not isinstance(options['skip_pole'], bool)):
		raise TypeError("options['skip_pole'] must be a scalar boolean")
	if (not np.isscalar(options['skip_res'])) or (not isinstance(options['skip_res'], bool)):
		raise TypeError("options['skip_res] must be a scalar boolean")
	if (not np.isscalar(options['cmplx_ss'])) or (not isinstance(options['cmplx_ss'], bool)):
		raise TypeError("options['cmplx_ss'] must be a scalar boolean")
	if (not np.isscalar(options['spy1'])) or (not isinstance(options['spy1'], bool)):
		raise TypeError("options['spy1'] must be a scalar boolean")
	if (not np.isscalar(options['spy2'])) or (not isinstance(options['spy2'], bool)):
		raise TypeError("options['spy2'] must be a scalar boolean")
	if (not np.isscalar(options['logx'])) or (not isinstance(options['logx'], bool)):
		raise TypeError("options['logx'] must be a scalar boolean")
	if (not np.isscalar(options['logy'])) or (not isinstance(options['logy'], bool)):
		raise TypeError("options['logy'] must be a scalar boolean")
	if (not np.isscalar(options['errplot'])) or (not isinstance(options['errplot'], bool)):
		raise TypeError("options['errplot'] must be a scalar boolean")
	if (not np.isscalar(options['phaseplot'])) or (not isinstance(options['phaseplot'], bool)):
		raise TypeError("options['phaseplot'] must be a scalar boolean")
	if (not np.isscalar(options['legend'])) or (not isinstance(options['legend'], bool)):
		raise TypeError("options['legend'] must be a scalar boolean")

	# check s
	if (not isinstance(s, np.ndarray)) or len(s.shape)>1:
		raise TypeError('s must be a one-dimensional numpy ndarray')
	Ns = s.shape[0]

	# check f
	if (not isinstance(f, np.ndarray)) or len(f.shape)!=2:
		raise TypeError('f must be a two-dimensional numpy ndarray')
	if f.shape[1] != Ns:
		raise ValueError('Second dimension of f does not match length of s')
	Nc = f.shape[0]

	# check poles
	if (not isinstance(poles, np.ndarray)) or len(poles.shape)>1:
		raise TypeError('poles must be a one-dimensional numpy ndarray')
	# correct indeterminacy issue
	if s[0]==0:
		if poles[0]==0 and poles[1]!=0:
			poles[0] = -1
		elif poles[1]==0 and poles[0]!=0:
			poles[1] = -1
		elif poles[0]==0 and poles[1]==0:
			poles[0] = -1 + 1j*10
			poles[1] = -1 - 1j*10
	N = poles.shape[0]

	# check weight
	if (not isinstance(weight, np.ndarray)) or len(weight.shape)!=2:
		raise TypeError('weight must be a two-dimensional numpy ndarray')
	if weight.shape[0] != 1 and weight.shape[0] !=Nc:
		raise ValueError('First dimension of weight is neither 1 nor matches first dimension of f.')
	if weight.shape[1] != Ns:
		raise ValueError('Second dimension of weight does not match length of s.')
	common_weight = (weight.shape[0] == 1)

	# import pyplot only if needed
	if ((not options['skip_pole']) and options['spy1']) or\
	   ((not options['skip_res'])  and options['spy2']):
		import matplotlib.pyplot as plt

	# Find out which starting poles are complex conjugate pairs:
	cindex = np.zeros((N), dtype=np.int_)
	m = 0
	while m < N:
		if poles[m].imag != 0:
			if poles[m].real!=poles[m+1].real or poles[m].imag!=-poles[m+1].imag:
				raise ValueError('Initial poles '+str(m)+' and '+str(m+1)\
				                 +' are subsequent but not complex conjugate.')
			cindex[m]   = 1
			cindex[m+1] = 2
			m += 1
		m += 1


	#==========================
	#=    POLE RELOCATION:    =
	#==========================

	if not options['skip_pole']:
		# Build system-matrix:
		Dk = np.empty((Ns,N+max(1, options['asymp'])), dtype=np.complex_)
		for m in range(N):
			if cindex[m]==0: # real pole
				Dk[:,m] = 1 / (s-poles[m])
			elif cindex[m]==1: # complex pole pair, 1st part
				Dk[:,m]   =  1/(s-poles[m]) +  1/(s-np.conj(poles[m]))
				Dk[:,m+1] = 1j/(s-poles[m]) - 1j/(s-np.conj(poles[m]))
		Dk[:,N] = 1
		if options['asymp']==2:
			Dk[:,N+1] = s

		if options['relax']: # relaxed VF
			# Calculate scaling for last row of LS-problem for relaxed VF (pole identification)
			scale=0
			for n in range(Nc):
				if common_weight:
					scale += np.linalg.norm(weight[0,:]*f[n,:])**2
				else:
					scale += np.linalg.norm(weight[n,:]*f[n,:])**2
			scale = np.sqrt(scale) / Ns

			# Construct linear system matrices
			AA = np.empty((Nc*(N+1),N+1))
			bb = np.zeros((Nc*(N+1)))
			Nleft = N + options['asymp']
			Ntot  = Nleft + N + 1
			for n in range(Nc):
				if n==Nc-1:
					A = np.zeros((2*Ns+1,Ntot))
				else:
					A = np.empty((2*Ns,Ntot))
				if common_weight:
					weig = weight[0,:,nax]
				else:
					weig = weight[n,:,nax]
				A[:Ns,    :Nleft] =  weig * Dk[:,:Nleft].real              # left block
				A[Ns:2*Ns,:Nleft] =  weig * Dk[:,:Nleft].imag
				A[:Ns,    Nleft:] = -weig * (Dk[:,:N+1] * f[n,:,nax]).real # right block
				A[Ns:2*Ns,Nleft:] = -weig * (Dk[:,:N+1] * f[n,:,nax]).imag

				# Integral criterion for sigma (nontriviality constraint):
				if n==Nc-1:
					A[2*Ns,Nleft:] = scale * np.sum(Dk[:,:N+1].real, axis=0)

				# partially solve for this element and add to AA and bb
				Q,R = np.linalg.qr(A, mode='reduced')
				AA[n*(N+1):(n+1)*(N+1),:] = R[Nleft:Ntot,Nleft:Ntot]
				if n==Nc-1:
					bb[n*(N+1):(n+1)*(N+1)] = Q[-1,Nleft:] * Ns * scale
			# end for n in range(Nc)

			# Solve the big system
			Escale = np.linalg.norm(AA, axis=0)
			AA /= Escale[nax,:]
			x = np.linalg.lstsq(AA, bb, rcond=-1)[0]
			x /= Escale
		# end if options['relax']

		# Situation: No relaxation, or produced D of sigma extremely small or large.
		# Solve again, without relaxation
		if (not options['relax']) or np.abs(x[-1])<TOLlow or np.abs(x[-1])>TOLhigh:
			# set D of sigma to fixed value
			if not options['relax']:
				Dnew = 1
			else:
				if x[-1]==0:
					Dnew = 1
				elif np.abs(x[-1])<TOLlow:
					Dnew = np.sign(x[-1]) * TOLlow
				elif np.abs(x[-1])>TOLhigh:
					Dnew = np.sign(x[-1]) * TOLhigh

			# Construct linear system matrices
			AA = np.empty((Nc*N,N))
			bb = np.empty((Nc*N))
			Nleft = N + options['asymp']
			Ntot = Nleft + N
			for n in range(Nc):
				A = np.empty((2*Ns,Ntot))
				b = np.empty((2*Ns))
				if common_weight:
					weig = weight[0,:]
				else:
					weig = weight[n,:]
				A[:Ns,:Nleft] =  weig[:,nax] * Dk[:,:Nleft].real        # left block
				A[Ns:,:Nleft] =  weig[:,nax] * Dk[:,:Nleft].imag
				A[:Ns,Nleft:] = -weig[:,nax] * (Dk[:,:N] * f[n,:,nax]).real # right block
				A[Ns:,Nleft:] = -weig[:,nax] * (Dk[:,:N] * f[n,:,nax]).imag
				b[:Ns] = Dnew * weig * f[n,:].real
				b[Ns:] = Dnew * weig * f[n,:].imag

				# partially solve for this element and add to AA and bb
				Q,R = np.linalg.qr(A, mode='reduced')
				AA[n*N:(n+1)*N,:] = R[Nleft:Ntot,Nleft:Ntot]
				bb[n*N:(n+1)*N]   = np.dot(Q[:,Nleft:Ntot].T, b)
			# end for n in range(Nc)

			# Solve the big system
			Escale = np.linalg.norm(AA, axis=0)
			AA /= Escale[nax,:]
			x = np.linalg.lstsq(AA, bb, rcond=-1)[0]
			x /= Escale
			x = np.append(x, Dnew)
		# end if (not options['relax']) or np.abs(x[-1])<TOLlow or np.abs(x[-1])>TOLhigh

		# plot results if spy1
		if options['spy1']:
			# convert x to complex residues and calculate sigma
			C = np.empty((N))
			for m in range(N):
				if cindex[m]==0:
					C[m] = x[m]
				elif cindex[m]==1:
					C[m]   = x[m] + 1j*x[m+1]
					C[m+1] = x[m] - 1j*x[m+1]
			sigma_rat = x[-1] + np.dot(1/(s[:,nax]-poles[nax,:]), C)
			freq = s.imag / (2*np.pi)
			# plot
			plt.figure()
			plt.plot(freq, np.abs(sigma_rat), 'b', label='sigma')
			plt.xlim(freq[0], freq[-1])
			if options['logx']:
				plt.xscale('log')
			if options['logy']:
				plt.yscale('log')
			plt.xlabel('Frequency [Hz]')
			plt.ylabel('Magnitude')
			plt.title('Sigma')
			if options['legend']:
				plt.legend(loc='best')
			if options['block']:
				plt.show()
			else:
				plt.draw()

		# build state space model of sigma
		LAMBD = np.zeros((N,N))
		B = np.ones((N,1))
		for m in range(N):
			if cindex[m]==0:
				LAMBD[m,m] = poles[m].real
			elif cindex[m]==1:
				LAMBD[m+1,m ]  =-poles[m].imag
				LAMBD[m,  m+1] = poles[m].imag
				LAMBD[m,  m]   = poles[m].real
				LAMBD[m+1,m+1] = poles[m].real
				B[m,0]   = 2
				B[m+1,0] = 0
		C = x[nax,:-1]
		D = x[-1]

		# Calculate the zeros of sigma
		newPoles = np.linalg.eigvals(LAMBD-B*C/D)
		if options['stable']:
			# Forcing unstable poles to be stable by flipping
			unstables = (newPoles.real>0)
			newPoles[unstables] -= 2 * newPoles[unstables].real

		# sort the new poles, real ones first, then in conjugate pairs (positive imag. part first)
		sortInd = np.lexsort((-newPoles.imag, np.abs(newPoles.real), np.abs(newPoles.imag)))
		newPoles = newPoles[sortInd]

		poles = newPoles

		# Find out which of the new poles are complex conjugate pairs:
		cindex = np.zeros((N), dtype=np.int_)
		m = 0
		while m < N:
			if poles[m].imag != 0:
				if poles[m].real!=poles[m+1].real or poles[m].imag!=-poles[m+1].imag:
					raise ValueError('Initial poles '+str(m)+' and '+str(m+1)\
					                 +' are subsequent but not complex conjugate.')
				cindex[m]   = 1
				cindex[m+1] = 2
				m += 1
			m += 1
	# end if not options['skip_pole']


	#=================================
	#=    RESIDUE IDENTIFICATION:    =
	#=================================

	SER = {}
	SER['C'] = None
	SER['D'] = None
	SER['E'] = None
	fit = None
	rmserr = 0
	if not options['skip_res']:
		# Build system matrix
		Dk=np.empty((Ns,N+options['asymp']), dtype=np.complex_)
		for m in range(N):
			if cindex[m]==0: # real pole
				Dk[:,m] = 1 / (s-poles[m])
			elif cindex[m]==1: # complex pole pair, 1st part
				Dk[:,m]   =  1/(s-poles[m]) +  1/(s-np.conj(poles[m]))
				Dk[:,m+1] = 1j/(s-poles[m]) - 1j/(s-np.conj(poles[m]))
		if options['asymp']>0:
			Dk[:,N] = 1
			if options['asymp']==2:
				Dk[:,N+1] = s

		A = np.zeros((2*Ns,N+options['asymp']))
		if common_weight:
			# fill linear system matrices
			A[:Ns,:] = weight[0,:,nax] * Dk.real
			A[Ns:,:] = weight[0,:,nax] * Dk.imag
			B = np.empty((2*Ns,Nc))
			B[:Ns,:] = weight[0,:,nax] * f.real.T
			B[Ns:,:] = weight[0,:,nax] * f.imag.T

			# solve the linear system
			Escale = np.linalg.norm(A, axis=0)
			A /= Escale[nax,:]
			X = np.linalg.lstsq(A, B, rcond=-1)[0]
			X /= Escale[:,nax]

			# construct residues from X
			SER['C'] = X[:N,:].T
			if options['asymp']>0:
				SER['D'] = X[N,:]
				if options['asymp']==2:
					SER['E'] = X[N+1,:]
		else: # not common_weight
			SER['C'] = np.empty((Nc,N))
			if options['asymp']>0:
				SER['D'] = np.empty((Nc))
				if options['asymp']==2:
					SER['E'] = np.empty((Nc))
			for n in range(Nc):
				# fill linear system matrices
				A[:Ns,:] = weight[n,:,nax] * Dk.real
				A[Ns:,:] = weight[n,:,nax] * Dk.imag
				b = np.empty((2*Ns))
				b[:Ns] = weight[n,:] * f[n,:].real
				b[Ns:] = weight[n,:] * f[n,:].imag

				# solve the linear system
				Escale = np.linalg.norm(A, axis=0)
				A /= Escale[nax,:]
				x = np.linalg.lstsq(A, b, rcond=-1)[0]
				x /= Escale

				# construct residues from X
				SER['C'][n,:] = x[:N]
				if options['asymp']>0:
					SER['D'][n] = x[N]
					if options['asymp']==2:
						SER['E'][n] = x[N+1]
		# end if common_weight

		# make SER['C'] complex again if wanted (calculate anyway for fit)
		SERC = np.empty((Nc,N), dtype=np.complex_)
		for m in range(N):
			if cindex[m]==0:
				SERC[:,m] = SER['C'][:,m]
			if cindex[m]==1:
				SERC[:,m]   = SER['C'][:,m] + 1j*SER['C'][:,m+1]
				SERC[:,m+1] = SER['C'][:,m] - 1j*SER['C'][:,m+1]
		if options['cmplx_ss']:
			SER['C'] = SERC

		# calculate fit and rms
		fit = np.dot(SERC, (1/(s[nax,:]-poles[:,nax])))
		if options['asymp']>0:
			fit += SER['D'][:,nax]
			if options['asymp']==2:
				fit += s[nax,:] * SER['E'][:,nax]
		rmserr = np.sqrt(np.sum(np.abs((fit-f)**2))) / np.sqrt(Nc*Ns)

		# make second plot(s) if requested
		if options['spy2']:
			# absolute fit plot
			freq = s.imag / (2*np.pi)
			plt.figure()
			handles = []
			labels  = []
			handles.append(plt.plot(freq, np.abs(f.T), 'b-'));    labels.append('Data')
			handles.append(plt.plot(freq, np.abs(fit.T), 'r--')); labels.append('FRVF')
			if options['errplot']:
				handles.append(plt.plot(freq, np.abs(f-fit).T, 'g')); labels.append('Deviation')
			plt.xlim(freq[0], freq[-1])
			if options['logx']:
				plt.xscale('log')
			if options['logy']:
				plt.yscale('log')
			plt.xlabel('Frequency [Hz]')
			plt.ylabel('Magnitude')
			plt.title('Approximation of f')
			if options['legend']:
				plt.legend([h[0] for h in handles], labels, loc='best')
			# phase plot
			if options['phaseplot']:
				plt.figure()
				handles = []
				labels  = []
				handles.append(plt.plot(freq, 180*np.unwrap(np.angle(f)).T/np.pi, 'b'))
				labels.append('Data')
				handles.append(plt.plot(freq, 180*np.unwrap(np.angle(fit)).T/np.pi, 'r--'))
				labels.append('FRVF')
				plt.xlim(freq[0], freq[-1])
				if options['logx']:
					plt.xscale('log')
				plt.xlabel('Frequency [Hz]')
				plt.ylabel('Phase angle [deg]')
				plt.title('Approximation of f')
				if options['legend']:
					plt.legend([h[0] for h in handles], labels, loc='best')
			if options['block']:
				plt.show()
			else:
				plt.draw()
	# end if not options['skip_res']

	# Convert into real state-space model if requested
	SER['B'] = np.ones((N))
	if not options['cmplx_ss']:
		SER['A'] = np.zeros((N,N))
		for m in range(N):
			if cindex[m]==0:
				SER['A'][m,m] = poles[m].real
			elif cindex[m]==1:
				SER['A'][m,  m]   =  poles[m].real
				SER['A'][m,  m+1] =  poles[m].imag
				SER['A'][m+1,m]   = -poles[m].imag
				SER['A'][m+1,m+1] =  poles[m].real
				SER['B'][m]   = 2
				SER['B'][m+1] = 0
	else:
		SER['A'] = np.diag(poles)

	return SER, poles, rmserr, fit


# Auxiliary functions

def ss2pr(SER, tri=False):
	''' Convert state-space model having COMMON pole set into pole-residue model.
	    Input:
	        SER: must have the format produced by vectfit3. Both formats determined by parameter
	               options['cmplx_ss'] are valid.
	               Ouptut from tri2full is also accepted (set tri to True in this case)
	        tri: False: output from vectfit3 (vector-valued) is assumed
	             True:  output from tri2full (matrix-valued) is assumed
	    Output:
	        R(Nc,N)     Residues, tri=False
	         (Nc,Nc,N)  Residues, tri=True
	        a(N)           poles
	        D(Nc)       Constant term (or None if SER['D']==None)
	        E(Nc)       Linear term   (or None if SER['E']==None)
	    This example script is part of the vector fitting package (VFIT3.zip)
	    Last revised: 19.02.2019.
	    Created by:    Bjorn Gustavsen.
	    Translated and modified by: Simon De Ridder                                             '''

	# Converting real-only state-space model into complex model, if necessary:
	A = SER['A']
	B = SER['B']
	C = SER['C']
	if np.amax(np.abs(A-np.diag(np.diag(A))), axis=None)>0.0:
		for m in range(A.shape[0]-1):
			if A[m,m+1]!=0.0:
				A[m,m]     = A[m,m]     + 1j*A[m,m+1] # first use should convert A to complex entirely
				A[m+1,m+1] = A[m+1,m+1] - 1j*A[m,m+1]
				B[m,:]   = (B[m,:]+B[m+1,:]) / 2
				B[m+1,:] = B[m,:]
				C[:,m]   = C[:,m] + 1j*C[:,m+1]
				C[:,m+1] = np.conj(C[:,m])

	# Converting complex state-space model into pole-residue model
	Nc = C.shape[0]
	if tri:
		N = int(A.shape[0]/Nc)
		R = np.empty((Nc,Nc,N), dtype=np.complex_)
		for m in range(N):
			R[:,:,m] = np.sum(C[:,nax,m:Nc*N+m:N] * np.transpose(B[m:Nc*N+m:N,:])[nax,:,:], axis=2)
	else:
		N = A.shape[0]
		R = np.empty((Nc,N), dtype=np.complex_)
		for m in range(N):
			R[:,m] = np.sum(C[:,m:Nc*N+m:N] * B[nax,m:Nc*N+m:N], axis=1)
	a = np.diag(A[:N,:N])
	return R,a,SER['D'],SER['E']

def tri2full(SER):
	''' Convert rational model of lower matrix triangle into state-space model of full matrix.
	    Input:
	        SER: Output structure from vectfit3 when fitting lower triangle of a square matrix
	             using a common pole set.
	             Both formats determined by parameter options['cmplx_ss'] are valid
	    Output:
	        SER: State-space model of full matrix (with common pole set)
	    This example script is part of the vector fitting package (VFIT3.zip) 
	    Last revised: 08.08.2008. 
	    Created by:    Bjorn Gustavsen.
	    Translated by: Simon De Ridder                                                          '''

	Nt = SER['C'].shape[0]
	Nc = -0.5 + np.sqrt(0.25+2.0*Nt)
	if abs(Nc-np.round(Nc))>1e-12:
		raise ValueError('Dimension is not compatible with a vectorized lower triangular matrix')
	Nc = int(np.round(Nc))
	N = SER['A'].shape[0]
	asymp = 2
	if SER['E'] is None:
		asymp = 1
		if SER['D'] is None:
			asymp = 0

	tell = 0
	AA = np.zeros((Nc*N,Nc*N), dtype=np.complex_)
	BB = np.zeros((Nc*N,Nc),   dtype=np.complex_)
	CC = np.empty((Nc,Nc*N),   dtype=np.complex_)
	DD = None
	EE = None
	if asymp >0:
		DD = np.empty((Nc,Nc), dtype=np.complex_)
		if asymp==2:
			EE = np.empty((Nc,Nc), dtype=np.complex_)
	for col in range(Nc):
		AA[col*N:(col+1)*N,col*N:(col+1)*N] = SER['A']
		BB[col*N:(col+1)*N,col] = SER['B']
		CC[col:,col*N:(col+1)*N] = SER['C'][tell:tell+Nc-col,:]
		CC[col,(col+1)*N:]       = SER['C'][tell+1:tell+Nc-col,:].flatten()
		if asymp>0:
			DD[col:Nc,col]   = SER['D'][tell:tell+Nc-col]
			DD[col,col+1:Nc] = SER['D'][tell+1:tell+Nc-col]
			if asymp==2:
				EE[col:Nc,col]   = SER['E'][tell:tell+Nc-col]
				EE[col,col+1:Nc] = SER['E'][tell+1:tell+Nc-col]
		tell += Nc - col

	return {'A':AA, 'B':BB, 'C':CC, 'D':DD, 'E':EE}