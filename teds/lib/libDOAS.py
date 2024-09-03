import sys
import numpy as np
import netCDF4 as nc
from scipy import interpolate
from matplotlib import pyplot
import warnings
from collections import OrderedDict
import math
import datetime

from teds import log
from teds.lib import constants
from teds.lib.libWrite import writevariablefromname

class OptimalEstimationException(Exception):
    pass

class DidNotConvergeWarning(Warning):
    pass

class Model(object):
    """The Model class for optimal estimation

The model class is a helper class for the L{optimal estimation<OE>} class.
It maintains the state vector, prior information, physical limits, and a list of
parameters that are fitted (or kept fixed).

This class must be subclassed and override the C{__call__} method.

Subclasses may define a "diagnostic" method which returns a dictionary with
additional diagnostic information.

@ivar _priorCov: Internal storage of the a priori covariance matrix. Accessed as
                 a property to check the size of the matrix, and handle a pure
                 diagonal matrix.
@type _priorCov: A L{np.matrix} instance.
@ivar _covar:    The current covariance matrix.
@type _covar:    A L{np.matrix} instance.
@ivar _prior:    The a priori state vector.
@type _prior:    A L{np.ndarray} instance.
@ivar _statevec: Internal storage of the state vector.
@type _statevec: A L{np.ndarray} instance.
@ivar modelCalculation: The model calculation for the current state vector.
@type modelCalculation: A L{np.ndarray} instance.
@ivar Jacobian:  The Jacobian matrix, the derivative of the model calculation with
                 respect to the state vector elements.
@type Jacobian:  A L{np.matrix} instance.
@ivar otherModelParameters: Other parameters for the model that are not part of
                 the state-vector, such as solar zenith angle.
@type otherModelParameters: A L{np.ndarray} instance.
@ivar parameterNames: The names of the parameters, for output purposes.
@type parameterNames: a list of strings.
@ivar stateVector: The public facing property of L{_statevec}.
@ivar priorCovariance: The public facing property of L{_covar}.
@ivar prior:     The public facing property of L{_prior}.
@ivar initialStateVector: The starting point of the fit. This may be different
                          from the prior. the default is to use the prior.
@type initialStateVector: A L{np.ndarray} instance.
@ivar observation: The measurement.
@type observation: A L{np.ndarray} instance with length M{n}.
@ivar observationError: The 1-sigma error on the observation.
@type observationError: A L{np.ndarray} instance with length M{n}.
@ivar independentVariable: The independent variable for the model, for our
                           retrievals usually M{\lambda}.
@type independentVariable: A L{np.ndarray} instance with length M{n}.
@ivar verbose:   An integer indicating how chatty (and plotting) the model should be.
                 0: be silent, 1: plot to file, 2: show.
"""
    def __init__(self, prior=None, priorCov=None,
                 otherModelParam=None, parameterNames=None, observation=None,
                 observationError=None, independentVariable=None,
                 stateVector=None, verbose=0, **kwargs):
        """The initializer method.

Here we check that all parameters are consistent (in length, mostly). Most
named arguments map to instance variable in a straightforward way.

@param prior: The value of the a priori state vector (required). This is a vector of length M{m}.
@param priorCov: The a priori covariance error matrix. The matrix can be supplied
              as a full error covariance matrix, including relations between
              parameters. If the parameter is supplied as a 1D-array, it is
              considered to be the diagonal of the covariance matrix. This matrix
              is either $M{m*m} or M{m} elements long.
@param otherModelParam: a vector with the remaining model parameters.
@param parameterNames: The names of all parameters (statevector + parameterNames, in that order).
@param observation: The observation vector (M{B{y}}).
@param observationError: The 1-sigma error on each observation.
@param independentVariable: The wavelength scale of the observation.
@param stateVector: The initial state vector.
@param verbose: flag to indicate chatty-ness.
@param kwargs: extra keyword parameters for subclasses.
@raise ValueError: if the a priori state vecor is missing.
"""
        self._priorCov = None
        self._covar = None
        self.modelCalculation = None
        self.Jacobian = None
        self.verbose = verbose
        if 'logfile' in kwargs and kwargs['logfile'] is not None:
            self.logfile = kwargs['logfile']
        else:
            self.logfile = sys.stdout

        if prior is None:
            raise ValueError("A value for the prior statevector is required.")

        self._prior = prior.copy()
        self.priorCovariance = priorCov

        self._statevec = prior.copy()
        self.otherModelParameters = otherModelParam if otherModelParam is not None else []
        self.parameterNames = parameterNames

        if len(self.parameterNames) != len(self.otherModelParameters) + len(self.prior):
            m = "The number of parameters is not correct, {0} != {1} + {2}"
            raise ValueError(m.format(len(self.parameterNames),
                                      len(self.otherModelParameters),
                                      len(self.prior)))

        if stateVector is not None:
            self.stateVector = stateVector

        self.initialStateVector = self.stateVector.copy()

        self.observation = observation
        self.observationError = observationError
        self.independentVariable = independentVariable

        # ignore floating point errors
        np.seterr(divide='ignore',
                  invalid='ignore',
                  over='ignore',
                  under='ignore')

    def checkObservation(self):
        """Check that the observation vectors are consistent.

The method returns a True value indicating a consistent set. The observation
vectors are consistens if they are all available, and have the same length.

@return: C{True} if C{self.independentVariable}, C{self.observation}, and
         C{self.observationError} are not C{None} and all have the same length.
@rtype: boolean
"""
        if (self.independentVariable is not None
            and self.observation is not None
            and self.observationError is not None):
            l = len(self.independentVariable)
            if (l == len(self.observation) and l == len(self.observationError)):
                return True
        return False

    def setObservationVector(self, independentVariable=None,
                             observation=None, observationError=None):
        """Set the observation vector.

Set the complete observation vector for the Model object. That is, the
independentVariable ('lambda'), the observation ('y') and the observation error
('sigma_y'). These attributes can be accessed directly, but this method also
performs a consistency check, and raises a ValueError if there is a length
mismatch.

@param independentVariable: The value for the new C{self.independentVariable}
                            instance variable
@param observation:         The value for the new C{self.observation} instance
                            variable
@param observationError:    The value for the new C{self.observationError}
                            instance variable
@raise ValueError:          If L{checkObservation} fails, a ValueError is raised.
"""
        if independentVariable is not None:
            self.independentVariable = independentVariable
        if observation is not None:
            self.observation = observation
        if observationError is not None:
            self.observationError = observationError

        if not self.checkObservation():
            raise ValueError("The elements of the observation vector "
                "('independentVariable', 'observation' and 'observationError') "
                "are not consistent.")

    def reset(self):
        self.covariance = None
        self.stateVector = self.initialStateVector

    def getPrior(self):
        return self._prior

    def setPrior(self, prior):
        if len(prior) != len(self._prior):
            m = "The number of state vector elements cannot be changed!"
            raise ValueError(m)
        self._prior = prior

    prior = property(getPrior, setPrior, None,
                     "The a priori state vector for this model")

    def getPriorCov(self):
        return self._priorCov

    def setPriorCov(self, priorCov):
        if priorCov is None:
            m = "A value for the prior covariance matrix is required."
            raise ValueError(m)
        elif (len(priorCov.shape) == 1
              and priorCov.shape[0] == len(self.prior)):
            self._priorCov = np.matrix(np.diag(priorCov**2))
        elif (len(priorCov.shape) == 2
              and priorCov.shape[0] == priorCov.shape[1]
              and priorCov.shape[0] == len(self.prior)):
            self._priorCov = np.matrix(priorCov.copy())
        else:
            m = """Unexpected number of dimensions of the error
covariance matrix, or length of dimensions not
consistent with prior."""
            raise ValueError(m)

    priorCovariance = property(getPriorCov, setPriorCov, None,
                               "The a priori error covariance for this model")

    def getStateVector(self):
        return self._statevec

    def setStateVector(self, stateVector):
        if len(stateVector) != len(self._statevec):
            m = "The number of state vector elements cannot be changed!"
            raise ValueError(m)
        self._statevec = stateVector

    stateVector = property(getStateVector, setStateVector, None,
                  "The state vector of the model. When setting the state vector "
                  "the values are automatically checked against the limits.")

    def getCovariance(self):
        return self._covar

    def setCovariance(self, covar):
        if covar is None:
            self._covar = None
        elif (len(covar.shape) == 1
              and covar.shape[0] == len(self.prior)):
            self._covar = np.matrix(np.diag(covar**2))
        elif (len(covar.shape) == 2
              and covar.shape[0] == covar.shape[1]
              and covar.shape[0] == len(self.prior)):
            self._covar = np.matrix(covar.copy())
        else:
            raise ValueError("Unexpected number of dimensions of the error "
                             "covariance matrix, or length of dimensions not "
                             "consistent with state vector.")

    covariance = property(getCovariance, setCovariance, None,
                          "The current covariance error matrix")

    def __str__(self):
        """informal representation"""
        m = "{name:15s}: {prior} +/- {error}; {current} +/- {curErr}; B F"
        r = [m.format(
             name="Name", prior="Prior", error="S_prior",
             current="Current", curErr="S_curr")]
        r.append("-"*len(r[0]))

        m = "{name:15s}: {prior:.4g} +/- {error:.4g}; {current:.4g} +/- {curErr:.4g}"
        for i in range(len(self.stateVector)):
            r.append(m.format(
                     name=self.parameterNames[i],
                     prior=self.prior[i],
                     error=np.sqrt(self.priorCovariance[i,i]),
                     current=self.stateVector[i],
                     curErr=np.sqrt(self.covariance[i,i])))
        r.append("-"*len(r[0]))

        r.append("Error covariance matrix:")
        for i in range(len(self.stateVector)):
            r.append(", ".join(["{0:.3g}".format(float(v))
                                for v in self.covariance[:,i]]))
        r.append("-"*len(r[0]))

        r.append("Correlation matrix:")
        for i in range(len(self.stateVector)):
            k = []
            covii = self.covariance[i,i]
            for j in range(len(self.stateVector)):
                covjj = self.covariance[j,j]
                v = self.covariance[i,j]/(np.sqrt(covii*covjj))
                k.append("{0:.3g}".format(v))
            r.append(", ".join(k))
        r.append("-"*len(r[0]))

        r.append("Other parameters:")
        len_sv = len(self.stateVector)
        m = "{name:15s}: {value}"
        for i in range(len(self.stateVector), len(self.parameterNames)):
            r.append(m.format(name=self.parameterNames[i],
                              value=self.otherModelParam[i-len_sv]))
        r.append("-"*len(r[0]))

        if self.verbose > 0 and self.modelCalculation is not None:
            r.append("x_indep\ty_measure\ty_error\ty_model\tdelta_y")
            s = "-"*len(r[-1])
            r.append(s)
            m = "{x:.3g}\t{y_measure:.4g}\t{y_error:.4g}\t{y_model:.4g}\t{delta:.4g}"
            for i in range(len(self.modelCalculation)):
                r.append(m.format(
                    x=self.independentVariable[i],
                    y_measure=self.observation[i],
                    y_error=self.observationError[i],
                    y_model=self.modelCalculation[i],
                    delta=((self.observation[i] - self.modelCalculation[i])/
                           self.observationError[i])))
            r.append(s)

        if self.verbose > 0 and self.Jacobian is not None:
            r.append("Jacobian:")
            r.append("x_indep\t" + ("\t".join(self.parameterNames[0:len_sv])))
            J = np.zeros((len(self.modelCalculation), len_sv+1))
            J[:, 0] = self.independentVariable
            J[:, 1:] = self.Jacobian

            for i in range(len(self.modelCalculation)):
                r.append("\t".join(["{0:.4g}".format(v) for v in J[i,:]]))
            del J

        return "\n".join(r)

    def namedModelParameter(self, s):
        """Obtain the value of a named parameter

The named parameter should be found in the L{stateVector} or in the
L{otherModelParameters}, using the names given in the L{parameterNames}
instance variable.

@param s: name of the parameter.
@return: value of named parameter or C{None} if parameter not found.
"""
        try:
            idx = self.parameterNames.index(s)
        except ValueError:
            return None

        if idx >= len(self.stateVector):
            idx -= len(self.stateVector)
            val = self.otherModelParameters[idx]
        else:
            val = self.stateVector[idx]

        return val

    def indexInStateVector(self, s):
        """Obtain the index within the state vector.

This index can be used for building the Jacobian. Other model parameters return
a fill value.

@param s: name of the parameter.
@return: index in state vector, or C{None} if the parameter is not in the state
         vector.
"""
        try:
            idx = self.parameterNames.index(s)
        except ValueError:
            return None

        if idx < len(self.stateVector):
            return idx
        else:
            return None

    def __call__(self, do_jacobian=True):
        """Calculate the model with the state vector that the model-object has.

The model should be calculated at the values given by the independentVariable
parameter. This method raises a NotImplementedError, as it is expected to be
overridden in a subclass.
"""
        raise NotImplementedError("The Model call is not implemented")

    def close(self):
        """Dummy method to allow subclasses to clean up after themselves."""
        pass

    @property
    def rmse(self):
        d = ((self.observation-self.modelCalculation)/self.observationError)**2
        return math.sqrt(np.sum(d)/len(self.observation))

    def plot(self, iteration=None, stateVectorConv=None):
        """A simple plot routine for the model."""
        r = ["{0}".format(self.__class__.__name__)]
        if iteration is not None:
            r.append("i: {0}".format(iteration))
        fmt = lambda a : ", ".join(["{0:.4g}".format(float(v)) for v in a])
        r.append("stateVector: {}".format(fmt(self.stateVector)))
        if stateVectorConv is not None:
            r.append("stateVectorConv: {0:.4g}".format(stateVectorConv))

        s = "; ".join(r)

        if iteration is not None and self.verbose > 0:
            print(s, file=self.logfile)

        if self.verbose > 1:

            r = ["{0}".format(self.__class__.__name__)]
            if iteration is not None:
                r.append(" i: {0}".format(iteration))
            stateVectorstr = ["{:.1E}".format(i) for i in self.stateVector]
            r.append("\nstateVector: {}".format(stateVectorstr))
            if stateVectorConv is not None:
                r.append("\nstateVectorConv: {0:.4g}".format(stateVectorConv))

            s = "".join(r)


            nplot = 2 + len(self.stateVector)
            fig = pyplot.figure(figsize=(9,10))
            fig.subplots_adjust(left=0.17, bottom=0.09, right=0.98,
                                top=0.92, wspace=0.12, hspace=0.2)
            ax = fig.add_subplot(nplot,1,1)
            ax.set_title(s,fontsize=10)
            # ax.set_title('DOAS OE. Iteration {}'.format(iteration))

            # ax.set_ylabel("$R [sr^{-1}]$")
            ax.set_ylabel("$R [sr^{-1}]$",rotation=0, fontsize=12, labelpad=45)

            ax.plot(self.independentVariable, self.observation, 'k',
                    label='measurement')
            ax.plot(self.independentVariable, self.modelCalculation, 'r',
                    label='model')
            ax.legend(loc='lower right',fontsize=8)

            ax.yaxis.tick_right()
            # ax.tick_params(axis='y', labelsize=6)

            ax.xaxis.set_visible(False)

            l = fig.add_subplot(nplot,1,2)
            l.plot(self.independentVariable,
                   (self.observation-self.modelCalculation)/self.observationError,
                   'k', label="err")
            # l.set_ylabel("$\Delta R/\sigma$")
            l.set_ylabel("$\Delta R/\sigma$",rotation=0, fontsize=12, labelpad=45)
            l.yaxis.tick_right()
            # l.tick_params(axis='y', labelsize=6)
            l.xaxis.set_visible(False)


            color = ['k-', 'r-', 'b-', 'g-', 'k--', 'r--', 'b--', 'g--', 'k-.',
                     'r-.', 'b-.', 'g-.', 'k:', 'r:', 'b:', 'g:']
            for i in range(len(self.stateVector)):
                name = self.parameterNames[i]
                k = fig.add_subplot(nplot,1,3+i)
                k.plot(self.independentVariable, self.Jacobian[:, i], 'k')
                # k.set_ylabel(r"$\partial R/\partial ({0})$".format(name.replace("_", " ")))
                k.set_ylabel(r"$\partial R/\partial ({0})$".format(name.replace("_", " ")), rotation=0, fontsize=12, labelpad=45)
                k.yaxis.tick_right()
                # k.tick_params(axis='y', labelsize=12)
                # k.yaxis.offsetText.set_fontsize(6)
                if i < (len(self.stateVector)-1):
                    k.axes.xaxis.set_visible(False)


            with warnings.catch_warnings():
                warnings.simplefilter("ignore")

                if self.verbose <= 1:
                    # fig.savefig("{0}_{1}_{2}.pdf".format(r[0], 
                    #             r[1].split()[1][:-1],
                    #             ("{0:02d}".format(iteration)
                    #              if iteration is not None
                    #              else "final")), transparent=True)
                    None
                else:
                    fig.show()


class IFDOEmodel(Model):
    """
     This subclass of L{Model} that models the reflectance spectrum.
    
     The derivatives of the signal with respect to the state vector 
     parameters are calculated analytically. 
    
     The instance variables L{modelCalculation} and L{Jacobian} are 
     set as well within the L{__call__} routine.
    """
    
    def __init__(self, prior=None, priorCov=None,
                 otherModelParam=None, parameterNames=None, observation=None,
                 observationError=None, independentVariable=None,
                 stateVector=None, verbose=0, **kwargs):
        """
         Init of the IFDOEmodel class.
         The via "kwargs" specified keywords are removed after they are used.
         The rest of the arguments (positional and keywords) is passed on
         to "super()".
        """

        self.IFDOEconfig  = kwargs['IFDOEconfig']
        self.IFDOErefspec = kwargs['IFDOErefspec']
                                                           
        super(IFDOEmodel, self).__init__(prior=prior, priorCov=priorCov,
              otherModelParam=otherModelParam,
              parameterNames=parameterNames, observation=observation,
              observationError=observationError,
              independentVariable=independentVariable,
              stateVector=stateVector, verbose=verbose)


    def __call__(self):
        """
         This part contains the model specific code. 
         The rest of the model-code is used as is.
         
         m  = the model 
         k  = the derivatives (Jacobian)
         sv = the state vesctor
        """

        m  = np.zeros(len(self.observation))
        k  = np.zeros((len(self.observation), len(self.prior)))
        sv = self.stateVector

      # The Rmod equation
      # -----------------
      # Several terms are given their own name, so that these
      # can be re-used easily in the derivaties further down.
        
      # a) polynomial term
        
        npoly = self.IFDOEconfig['polynomial_coefs']

        POLYterm = np.zeros(len(self.observation))

        for ipoly in range (0,npoly,1):
            POLYterm = POLYterm + sv[ipoly]*(self.independentVariable**ipoly)
        
        m = POLYterm
        nterms = npoly
        
      # b) exponential with the summation terms: trace gases, linear Ring
        
        termAdded = 0       # introduced for local book keeping
        
        ngases = self.IFDOEconfig['nr_trace_gases']

        SUMterm = np.zeros(len(self.observation))

        if ( ngases > 0 ):
            termAdded = 1
            for igas in range(0,ngases,1):
                gas = self.IFDOEconfig['trace_gas_list'][igas]
                SUMterm = SUMterm + sv[nterms+igas]*self.IFDOErefspec[gas]
            nterms = nterms+ngases

        if ( self.IFDOEconfig['ring_fit_term'] == 'sigma'):
            termAdded = 1
            SUMterm = SUMterm + sv[nterms]*self.IFDOErefspec['Dring']
            nterms = nterms+1
            
        if ( termAdded == 1 ):
            # Introduced to avoid doing the exp() if nothing was
            # added, i.e. if SUMterm=0.
            EXPterm = np.exp(-SUMterm)
            m = POLYterm * EXPterm
        else:
            EXPterm = np.ones(len(self.observation))
        
      # c) non-linear Ring term
        
        if ( self.IFDOEconfig['ring_fit_term'] == 'Iring'):
            nlRingTerm = 1.0 + sv[nterms]*self.IFDOErefspec['Cring']/self.IFDOErefspec['solar']
            m = m * nlRingTerm
            nterms = nterms + 1
        else:
            nlRingTerm = np.ones(len(self.observation))
       
      # d) intensity offset term
        
        if ( self.IFDOEconfig['intensity_coefs'] > 0 ):
            # C0 term
            OFFSETterm = np.zeros(len(self.observation))
            Soff = self.IFDOEconfig['intensity_scaling_factor']*np.ones(len(self.independentVariable))
            OFFSETterm = sv[nterms] * Soff/self.IFDOErefspec['solar']
            if ( self.IFDOEconfig['intensity_coefs'] == 2 ):
                # C1 term
                OFFSETterm = OFFSETterm + sv[nterms+1]*self.independentVariable * self.IFDOEconfig['intensity_scaling_factor']/self.IFDOErefspec['solar']
            m = m + OFFSETterm
            nterms = nterms + self.IFDOEconfig['intensity_coefs']
        
        # print('>>> At end of model setup: nterms = ',nterms)
        
        
      # The derivatives for the Jacobian matrix
      # ---------------------------------------
        
      # a) polynomial term    
                   
        for ipoly in range (0,npoly,1):
            k[:,ipoly] = (self.independentVariable**ipoly) * EXPterm * nlRingTerm
                   
        nterms = npoly

      # b) exponential with the summation terms: trace gases, linear Ring

        if ( ngases > 0 ):
            for igas in range(0,ngases,1):
                gas = self.IFDOEconfig['trace_gas_list'][igas]
                k[:,nterms+igas] = POLYterm * -self.IFDOErefspec[gas] * EXPterm * nlRingTerm 
                # None
            nterms = nterms+ngases
        
        if ( self.IFDOEconfig['ring_fit_term'] == 'sigma'):
            k[:,nterms] = POLYterm * -self.IFDOErefspec['Dring'] * EXPterm 
            nterms = nterms+1

      # c) non-linear Ring term
        
        if ( self.IFDOEconfig['ring_fit_term'] == 'Iring'):
            k[:,nterms] = POLYterm * EXPterm * self.IFDOErefspec['Cring']/self.IFDOErefspec['solar']
            nterms = nterms+1
                
      # d) intensity offset term
        
        if ( self.IFDOEconfig['intensity_coefs'] > 0 ):
            # C0 derivate
            k[:,nterms] = Soff/self.IFDOErefspec['solar']
            nterms = nterms+1
            if ( self.IFDOEconfig['intensity_coefs'] == 2 ):
                # C1 derivative
                k[:,nterms] = self.independentVariable * Soff/self.IFDOErefspec['solar']
                nterms = nterms+1

        # print('>>> At end of Jacobians setup: nterms = ',nterms)
        

      # Finish up
      # ---------
        
        self.modelCalculation, self.Jacobian = m, k

        return m, k


class OE(object):
    """The optimal estimation engine.

The optimal estimation engine uses thh model-object to obtain the measurements
and the a priori information.

After creating an instance, the actual fit is performed in the L{__call__} method.

@ivar model: The C{Model} instance which performs the forward calculation and
             holds the a priori information.
@ivar priorSinvh: One of the transformation matrices obtained with singular
             value decomposition. This one is derived from the a priori error
             covariance matrix. In pure latex notation: C{\mathbf{S}_a^{-1/2}}
@ivar priorSinv: One of the transformation matrices obtained with singular
             value decomposition. This one is derived from the a priori error
             covariance matrix. In pure latex notation: C{\mathbf{S}_a^{-1}}
@ivar priorSh: One of the transformation matrices obtained with singular value
             decomposition. This one is derived from the a priori error
             covariance matrix. In pure latex notation: C{\mathbf{S}_a^{1/2}}
@ivar errSinvhD: The diagonal of the inverse square root of the mearurement
             error covariance matrix. Since the measurement error covariance
             matrix itself is diagonal, this matrix is diagonal as well.
@ivar errSinvD: The diagonal of the inverse of the mearurement
             error covariance matrix. Since the measurement error covariance
             matrix itself is diagonal, this matrix is diagonal as well.
@ivar errShD: The diagonal of the square root of the mearurement
             error covariance matrix. Since the measurement error covariance
             matrix itself is diagonal, this matrix is diagonal as well.
@ivar U:     Result of the singular value decomposition of the Jacobian.
             See L{w} and L{V} as well.
@ivar w:     Result of the singular value decomposition of the Jacobian.
             See L{U} and L{V} as well.
@ivar V:     Result of the singular value decomposition of the Jacobian.
             See L{U} and L{wV} as well.
@ivar maxiter: The maximum number of iterations.
@ivar stateVectorConvThreshold: Convergence criterium, using the change in the
             state vector bewteen iterations as a criterium.
@ivar answer: Place to collect the diagnostic values.
"""
    def __init__(self, model=None, maxiter=12, stateVectorConvThreshold=1.0):
        """Initializer method.

The initializer method.

@param model: The model object.
@type model: an instance of L{Model}.
@param maxiter: The maximum number of iterations before giving up. Defaults to 12.
@type maxiter: integer
@param stateVectorConvThreshold: Convergence treshold, defaults to 1.0.
"""
        self.model = model

        self.priorSinvh = None
        self.priorSinv = None
        self.priorSh = None

        self.errSinvhD = None
        self.errSinvD = None
        self.errShD = None

        self.U = None
        self.w = None
        self.V = None

        self.maxiter = maxiter
        self.stateVectorConvThreshold = stateVectorConvThreshold

        self.answer = None

    def reset(self):
        self.priorSinvh = None
        self.priorSinv = None
        self.priorSh = None

        self.errSinvhD = None
        self.errSinvD = None
        self.errShD = None

        self.U = None
        self.w = None
        self.V = None

        self.model.reset()

    def start(self):
        """Set up the transformation matrices.

Use the a priori error covariance matrix and measurement errors to initialize
the L{priorSinvh}, L{priorSinv}, L{priorSh}, L{errSinvhD}, L{errSinvD} and
L{errShD} instance variables. These do not change between iterations.

called from the L{__call__} method.
"""
        self.model.stateVector = self.model.initialStateVector
        self.transformPriorErrorCovariance()
        self.transformMeasurementError()

    @property
    def Sinv(self):
        """Calculate inverse of the a posteriori error covariance matrix."""
        Wplus = np.matrix(np.diag(self.w**2 + 1.0))
        return self.priorSinvh * self.V.T * Wplus * self.V * self.priorSinvh

    def DecomposeJacobian(self, K):
        """Decompose the pre-whitened Jacobian using SVD.

This will fill the L{U}, L{w} and L{V} instance variables. Since these change
with each iteration, it is important to call this method during the iteration
sequence.

@param K: The Jacobian.
"""
        tmp = self.errSinvh * K * self.priorSh
        self.U, self.w, VT = np.linalg.linalg.svd(tmp, full_matrices=False)
        self.V = VT.T

    def iterations(self):
        """Iterate to find a solution.

Run the model once for each iteration, and find a new state-vector. Continue
until the maximum number of iterations is reached, or a solution is found
(whichever comes first).

The method itself is described in the Disamar manual and in the support
documentation of this code.
"""
        i = 0
        stateVectorConv = self.stateVectorConvThreshold * 1.0e6
        n = len(self.model.stateVector)
        self.answer = None

        while ((i < self.maxiter)
            and (stateVectorConv > self.stateVectorConvThreshold)):

            F, K = self.model()

            if np.any(np.isnan(F)) or np.any(np.isnan(K)):
                m = "Iteration {0} failure of model."
                raise OptimalEstimationException(m.format(i))

            if self.model.verbose > 0:
                self.model.plot(i+1, stateVectorConv)

            try:
                self.DecomposeJacobian(K)
            except np.linalg.LinAlgError:
                m = "Iteration {0} failure in decomposition."
                raise OptimalEstimationException(m.format(i))

            statevectorOffset = (self.V.T * self.priorSinvh *
                        np.matrix(np.array(self.model.stateVector) - np.array(self.model.prior) ).T)
            measurementOffset = (self.U.T * self.errSinvh *
                                       np.matrix(self.model.observation - F).T)

            newState = np.matrix((self.w *
                              (measurementOffset.A1 +
                              self.w * statevectorOffset.A1))/(self.w**2+1.0)).T
            newState = self.priorSh * self.V * newState
            newState = newState.A1 + self.model.prior

            stateVectorConv = ((np.matrix(newState - self.model.stateVector) *
                self.Sinv * np.matrix(newState - self.model.stateVector).T)/n)[0,0]

            self.model.stateVector = newState

            if i == 0:
                stateVectorConv = self.stateVectorConvThreshold * 1.0e6

            i += 1

        F, K = self.model()
        if self.model.verbose > 0:
            self.model.plot(i+1, stateVectorConv)

        try:
            self.DecomposeJacobian(K)
        except np.linalg.LinAlgError:
            raise OptimalEstimationException("Failure in decomposition.")

        Wplus2 = np.matrix(np.diag(1.0/(self.w**2+1.0)))
        self.model.covariance = (self.priorSh * self.V * Wplus2 *
                                                        self.V.T * self.priorSh)
        return i, stateVectorConv

    def diagnostic(self, numberOfIterations, stateVectorConv):
        """Calculate diagnostic parameters of a fit.

Calculate convergence parameters and return them in a dictionary.
The dictionary has the following keys: 'converged', 'number of iterations',
'state vector convergence criterium', 'state vector boundary violation',
'state vector boundary violation detail', 'fixed parameter list',
'a priori state vector', 'a priori error covariance matrix',
'a posteriori state vector', 'a posteriori error covariance matrix',
'a posteriori noise error covariance matrix', 'cost function',
'degrees of freedom', 'averaging kernel', 'gain matrix', and
'weighted root mean square difference'.
These parameters are decribed in the support documentation of this code.

The L{answer} instance variable is set to this dictionary as well.

@param numberOfIterations: actual number of iterations.
@param stateVectorConv: state vector convergence criterium

@rtype: C{OrderedDict()}
@return: Dictionary with keys given above.
"""
        result = OrderedDict()
        result['number of iterations'] = numberOfIterations
        result['state vector element names'] = self.model.parameterNames
        result['state vector convergence criterium'] = stateVectorConv
        result['a priori state vector'] = self.model.prior
        result['a priori error covariance matrix'] = self.model.priorCovariance
        result['a posteriori state vector'] = self.model.stateVector
        result['a posteriori error covariance matrix'] = self.model.covariance
        Wplus2 = np.matrix(np.diag(self.w**2/(self.w**2+1.0)))
        AK = self.priorSh * self.V * Wplus2 * self.V.T * self.priorSinvh
        result['a posteriori noise error covariance matrix'] = AK * self.model.covariance
        result['cost function'] = self.costFunction
        result['degrees of freedom'] = AK.trace()[0,0]
        result['averaging kernel'] = AK
        Wplus = np.matrix(np.diag(self.w/(self.w**2+1.0)))
        result['gain matrix'] = self.priorSh * self.V * Wplus * self.U.T * self.errSinvh
        result['weighted root mean square difference'] = self.model.rmse
        result['correlation matrix'] = self.correlation_matrix
        result['condition number Jacobian'] = np.linalg.cond(self.model.Jacobian)
        condition_number = np.linalg.cond(result['a posteriori error covariance matrix'])
        result['condition number a posteriori error covariance matrix'] = condition_number
        result['state vector boundary violation'] = False
        result['state vector boundary violation detail'] = None
        result['fixed parameter list'] = []


        result['converged'] = ((stateVectorConv < self.stateVectorConvThreshold)
                                and (numberOfIterations < self.maxiter))

        if hasattr(self.model, "diagnostic"):
            d = self.model.diagnostic()
            for k, v in list(d.items()):
                result[k] = v

        self.answer = result

        return result

    @property
    def correlation_matrix(self):
        """return the correlation matrix"""
        correlation_matrix = self.model.covariance.copy()
        sigmaD = np.sqrt(np.diag(correlation_matrix))
        for ii in range(correlation_matrix.shape[0]):
            for jj in range(correlation_matrix.shape[1]):
                correlation_matrix[ii, jj] /= sigmaD[ii] * sigmaD[jj]
        return correlation_matrix

    def __call__(self):
        """Fit the model parameters to match the data.

After creating the OE object, calling this method will perform the whole fit.

@param verbose: if C{True} create plots of the function, derivatives and
residual for each iteration.
@rtype: C{OrderedDict()}
@return: Dictionary as returned by L{diagnostic}.
"""
        self.start()
        numberOfIterations, stateVectorConv = self.iterations()
        if numberOfIterations <= 1:
            self.answer = None
            raise DidNotConvergeWarning("Number of iterations <= 1.")
        result = self.diagnostic(numberOfIterations, stateVectorConv)
        return result

    def transformPriorErrorCovariance(self):
        """Set up the prior error covariance transformation matrices.

Use the a priori error covariance matrix to initialize
the L{priorSinvh}, L{priorSinv} and L{priorSh} instance variables.
These do not change between iterations.
"""
        U_a, w_a, V_aT = np.linalg.linalg.svd(self.model.priorCovariance,
                                              full_matrices=False)
        V_a = V_aT.T
        self.priorSinvh = V_a * np.matrix(np.diag(np.sqrt(1.0/w_a))) * U_a.T
        self.priorSh = U_a * np.matrix(np.diag(np.sqrt(w_a))) * V_aT
        self.priorSinv = V_a * np.matrix(np.diag(1.0/w_a)) * U_a.T

    def transformMeasurementError(self):
        """Set up the measurement error transformation matrices.

Use the measurement errors to initialize the L{errSinvhD}, L{errSinvD} and
L{errShD} instance variables. These do not change between iterations.
"""
        var = self.model.observationError**2
        self.errShD = self.model.observationError
        self.errSinvD = 1.0/var
        self.errSinvhD = np.sqrt(self.errSinvD)

    @property
    def errSinv(self):
        """Diagonal matrix, inverse measurement variance matrix"""
        return np.matrix(np.diag(self.errSinvD))

    @property
    def errSinvh(self):
        """Diagonal matrix, inverse measurement error matrix"""
        return np.matrix(np.diag(self.errSinvhD))

    @property
    def errSh(self):
        """Diagonal matrix, measurement error matrix"""
        return np.matrix(np.diag(self.errShD))

    @property
    def costFunction(self):
        """Calculate the cost function.

Otherwise know as the chi square value.

@return: chi square
"""
        priorDiff = np.matrix(self.model.stateVector - self.model.prior).T
        measurementDiff = np.matrix(self.model.observation
                                    - self.model.modelCalculation).T
        chisq = measurementDiff.T * self.errSinv * measurementDiff
        chisq += priorDiff.T * self.priorSinv * priorDiff

        return chisq[0,0]

    def _matrixToStr(self, name, mat):
        """Helper method to create a string representation of a matrix"""
        r = []
        r.append("\n" + name)
        for  i in range(len(self.answer['a priori state vector'])):
            r.append(", ".join(["{0:=+10.4g}".format(float(v))
                                for v in mat[:, i]]))
        return "\n".join(r)

    def __str__(self):
        if self.answer is None:
            return "None"
        r = []

        for k in ["converged", "number of iterations"]:
            r.append("{0}: {1}".format(k, self.answer[k]))

        for k in ['state vector convergence criterium', 'cost function',
                  'degrees of freedom', 'weighted root mean square difference']:
            r.append("{0}: {1:.5g}".format(k, self.answer[k]))
        r.append("\na priori state vector:\n{0}".format(
                 ", ".join([str(float(v))
                            for v in self.answer['a priori state vector']])))

        m = self.answer['a priori error covariance matrix']
        r.append(self._matrixToStr("a priori error covariance matrix:", m))

        if self.answer['state vector boundary violation']:
            r.append("\nThe final state has one or more boundary violations")
            r.append("a posteriori state vector: {0}".format(
                     ", ".join([str(float(v))
                           for v in self.answer['a posteriori state vector']])))
            r.append("boundary violation status: {0}".format(
                     ", ".join(["{0}".format(bool(v))
              for v in self.answer['state vector boundary violation detail']])))
        else:
            r.append("\na posteriori state vector: {0}".format(
                     ", ".join([str(float(v))
                           for v in self.answer['a posteriori state vector']])))
        r.append("Fixed parameters list      : {0}".format(
                 ", ".join([str(bool(v))
                                for v in self.answer['fixed parameter list']])))

        m = self.answer['a posteriori error covariance matrix']
        r.append(self._matrixToStr("a posteriori error covariance matrix:", m))

        m = self.answer['a posteriori noise error covariance matrix']
        r.append(self._matrixToStr("a posteriori noise error covariance matrix:", m))

        m = self.answer['averaging kernel']
        r.append(self._matrixToStr("averaging kernel:", m))

        m = self.answer['gain matrix']
        r.append(self._matrixToStr("gain matrix:", m))

        return "\n".join(r)


class IFDOEcalibration(Model):
    """
     This subclass of L{Model} that models the reflectance spectrum.
    
     The derivatives of the signal with respect to the state vector 
     parameters are calculated analytically. 
    
     The instance variables L{modelCalculation} and L{Jacobian} are 
     set as well within the L{__call__} routine.
    """
    
    def __init__(self, prior=None, priorCov=None,
                 otherModelParam=None, parameterNames=None, observation=None,
                 observationError=None, independentVariable=None,
                 stateVector=None, verbose=0, **kwargs):
        """
         Init of the IFDOEmodel class.
         The via "kwargs" specified keywords are removed after they are used.
         The rest of the arguments (positional and keywords) is passed on
         to "super()".
        """

        self.cfg  = kwargs['cfg']
        self.IFDOErefspec = kwargs['IFDOErefspec']
        self.groundPixel = kwargs['groundPixel']

                                                           
        super(IFDOEcalibration, self).__init__(prior=prior, priorCov=priorCov,
              otherModelParam=otherModelParam,
              parameterNames=parameterNames, observation=observation,
              observationError=observationError,
              independentVariable=independentVariable,
              stateVector=stateVector, verbose=verbose)


    def __call__(self):
        """
         This part contains the model specific code. 
         The rest of the model-code is used as is.
         
         m  = the model 
         k  = the derivatives (Jacobian)
         sv = the state vesctor
        """

        m  = np.zeros(len(self.observation))
        k  = np.zeros((len(self.observation), len(self.prior)))
        sv = self.stateVector



        # Wavelength calibration terms
        # ----------------------------

        lambdaCal = self.independentVariable
        lambdaCentre = (self.cfg['fit_window'][1] - self.cfg['fit_window'][0] ) / 2 + self.cfg['fit_window'][0]
        lambdaStar = 2*(self.independentVariable - (lambdaCentre) ) / (self.cfg['fit_window'][1] - self.cfg['fit_window'][0])


        npoly = self.cfg['polynomial_coefs']
        ncal = 0

        if self.cfg['wvl_cal_offset']:
            offsetterm = sv[npoly]
            lambdaCal = lambdaCal + offsetterm
            ncal += 1


        if self.cfg['wvl_cal_stretch']:
            stretchterm = sv[npoly+ncal]* lambdaStar
            lambdaCal = lambdaCal + stretchterm
            ncal += 1


        # Interpolate ref. spec. to lambda cal
        # -----------------------------------

        ngases = self.cfg['nr_trace_gases']

        IFDOErefspecCalGrid = {}
        IFDOErefspecCalGridDeriv = {}

        for igas in range(0,ngases,1):
            gas = self.cfg['trace_gas_list'][igas]

            IFDOErefspecCalGrid[gas], IFDOErefspecCalGridDeriv[gas]  = self.interp_refspec_lambdacal(self.IFDOErefspec['xsWvl'], lambdaCal, self.IFDOErefspec[gas][self.groundPixel,:])

        if ( self.cfg['ring_fit_term'] == 'sigma'):
            IFDOErefspecCalGrid['Dring'], IFDOErefspecCalGridDeriv['Dring']  = self.interp_refspec_lambdacal(self.IFDOErefspec['solarWvl'], lambdaCal, self.IFDOErefspec['Dring'][self.groundPixel,:])

        if 'solarSpline' in self.IFDOErefspec:
            if ( self.cfg['ring_fit_term'] == 'Iring'):
                IFDOErefspecCalGrid['Cring'], IFDOErefspecCalGridDeriv['Cring']  = self.interp_refspec_lambdacal_precalc(self.IFDOErefspec['ringSpline'][self.groundPixel], lambdaCal)
            elif ( self.cfg['ring_fit_term'] == 'sigma'):
                IFDOErefspecCalGrid['Dring'], IFDOErefspecCalGridDeriv['Dring']  = self.interp_refspec_lambdacal_precalc(self.IFDOErefspec['ringSpline'][self.groundPixel], lambdaCal)
        
            IFDOErefspecCalGrid['solar'], IFDOErefspecCalGridDeriv['solar']  = self.interp_refspec_lambdacal_precalc(self.IFDOErefspec['solarSpline'][self.groundPixel], lambdaCal)

        else:
            if ( self.cfg['ring_fit_term'] == 'Iring'):
                IFDOErefspecCalGrid['Cring'], IFDOErefspecCalGridDeriv['Cring']  = self.interp_refspec_lambdacal(self.IFDOErefspec['solarWvl'], lambdaCal, self.IFDOErefspec['Cring'][self.groundPixel,:])
            elif ( self.cfg['ring_fit_term'] == 'sigma'):
                IFDOErefspecCalGrid['Dring'], IFDOErefspecCalGridDeriv['Dring']  = self.interp_refspec_lambdacal(self.IFDOErefspec['solarWvl'], lambdaCal, self.IFDOErefspec['Dring'][self.groundPixel,:])
        
            IFDOErefspecCalGrid['solar'], IFDOErefspecCalGridDeriv['solar']  = self.interp_refspec_lambdacal(self.IFDOErefspec['solarWvl'], lambdaCal, self.IFDOErefspec['solar'][self.groundPixel,:])


      # The Rmod equation
      # -----------------
      # Several terms are given their own name, so that these
      # can be re-used easily in the derivaties further down.
        
      # a) polynomial term
        
        

        POLYterm = np.zeros(len(self.observation))

        for ipoly in range (0,npoly,1):
            POLYterm = POLYterm + sv[ipoly]*(lambdaStar**ipoly)
        
        m = POLYterm
        nterms = npoly + ncal

        
      # b) exponential with the summation terms: trace gases, linear Ring
        
        termAdded = 0       # introduced for local book keeping
        

        SUMterm = np.zeros(len(self.observation))

        if ( ngases > 0 ):
            termAdded = 1
            for igas in range(0,ngases,1):
                gas = self.cfg['trace_gas_list'][igas]
                SUMterm = SUMterm + sv[nterms+igas]*IFDOErefspecCalGrid[gas]
            nterms = nterms+ngases

        if ( self.cfg['ring_fit_term'] == 'sigma'):
            termAdded = 1
            SUMterm = SUMterm + sv[nterms]*IFDOErefspecCalGrid['Dring']
            nterms = nterms+1
            
        if ( termAdded == 1 ):
            # Introduced to avoid doing the exp() if nothing was
            # added, i.e. if SUMterm=0.
            EXPterm = np.exp(-SUMterm)
            m = POLYterm * EXPterm
        else:
            EXPterm = np.ones(len(self.observation))
        
      # c) non-linear Ring term
        
        if ( self.cfg['ring_fit_term'] == 'Iring'):
            nlRingTerm = IFDOErefspecCalGrid['solar'] + sv[nterms]*IFDOErefspecCalGrid['Cring']
            m = m * nlRingTerm
            nterms = nterms + 1

        else: # if linear or no ring is used nlRingTerm contains only E0
            nlRingTerm = IFDOErefspecCalGrid['solar']
            m = m * nlRingTerm
       
        
        # print('>>> At end of model setup: nterms = ',nterms)
        
        
      # The derivatives for the Jacobian matrix
      # ---------------------------------------
        
        # Wavelength calibration terms


        if self.cfg['wvl_cal_offset']:


            SUMderivterm = np.zeros(len(self.observation))

            if ( ngases > 0 ):
                for igas in range(0,ngases,1):
                    gas = self.cfg['trace_gas_list'][igas]
                    SUMderivterm = SUMderivterm + sv[ncal+npoly+igas]*IFDOErefspecCalGridDeriv[gas]

            if ( self.cfg['ring_fit_term'] == 'sigma'):
                SUMderivterm = SUMderivterm + sv[ncal+npoly+ngases]*IFDOErefspecCalGridDeriv['Dring']

            nlRingTermDeriv = np.zeros(len(self.observation))

            if ( self.cfg['ring_fit_term'] == 'Iring'):
                nlRingTermDeriv = nlRingTermDeriv + sv[ncal+npoly+ngases] * IFDOErefspecCalGridDeriv['Cring']

            k[:,npoly] = POLYterm * EXPterm * ( -SUMderivterm*nlRingTerm + (IFDOErefspecCalGridDeriv['solar'] + nlRingTermDeriv )  )


        if self.cfg['wvl_cal_stretch']:
            SUMderivterm2 = np.zeros(len(self.observation))

            if ( ngases > 0 ):
                for igas in range(0,ngases,1):
                    gas = self.cfg['trace_gas_list'][igas]
                    SUMderivterm2 = SUMderivterm2 + sv[ncal+npoly+igas]*lambdaStar*IFDOErefspecCalGridDeriv[gas]

            if ( self.cfg['ring_fit_term'] == 'sigma'):
                SUMderivterm2 = SUMderivterm2 + sv[ncal+npoly+ngases]*lambdaStar*IFDOErefspecCalGridDeriv['Dring']

            nlRingTermDeriv = np.zeros(len(self.observation))

            if ( self.cfg['ring_fit_term'] == 'Iring'):
                nlRingTermDeriv = nlRingTermDeriv + sv[ncal+npoly+ngases] * lambdaStar  * IFDOErefspecCalGridDeriv['Cring']

            k[:,npoly+ncal-1] = POLYterm * EXPterm * ( -SUMderivterm2*nlRingTerm + (lambdaStar*IFDOErefspecCalGridDeriv['solar'] + nlRingTermDeriv )  )


      # a) polynomial term    
                   
        for ipoly in range (0,npoly,1):
            k[:,ipoly] = (lambdaStar**ipoly) * EXPterm * nlRingTerm
                   
        nterms = ncal+npoly

      # b) exponential with the summation terms: trace gases, linear Ring

        if ( ngases > 0 ):
            for igas in range(0,ngases,1):
                gas = self.cfg['trace_gas_list'][igas]
                k[:,nterms+igas] = POLYterm * -IFDOErefspecCalGrid[gas] * EXPterm * nlRingTerm 
                # None
            nterms = nterms+ngases
        
        if ( self.cfg['ring_fit_term'] == 'sigma'):
            # if linear ring is used nlRingTerm contains only E0
            k[:,nterms] = POLYterm * -IFDOErefspecCalGrid['Dring'] * EXPterm * nlRingTerm
            nterms = nterms+1

      # c) non-linear Ring term
        
        if ( self.cfg['ring_fit_term'] == 'Iring'):
            k[:,nterms] = POLYterm * EXPterm * IFDOErefspecCalGrid['Cring']
            nterms = nterms+1
                

        # print('>>> At end of Jacobians setup: nterms = ',nterms)

        
      # Finish up
      # ---------
        
        self.modelCalculation, self.Jacobian = m, k

        # pdb.set_trace()

        return m, k

    def interp_refspec_lambdacal(self,lambdaNom, lambdaCal, refSpec):
        '''
        Interpolate refSpec (given on lambdaNom) to new wavelength grid (lambdaCal).
        Also returns derivative of refSpec to lambdaCal.
        '''

        # create spline for ref.spec. using lambda nom
        spline = interpolate.InterpolatedUnivariateSpline(lambdaNom, refSpec, ext=0, k=4)

        # interpolate to lambda cal 
        refSpecCalGrid = spline(lambdaCal)

        # get derivative to lambda cal
        refSpecCalDeriv = spline.derivative(n=1)(lambdaCal)

        return refSpecCalGrid, refSpecCalDeriv

    def interp_refspec_lambdacal_precalc(self, spline, lambdaCal):
        '''
        Interpolate refSpec (given on lambdaNom) to new wavelength grid (lambdaCal).
        Using precalculated spline
        Also returns derivative of refSpec to lambdaCal.
        '''
        # interpolate to lambda cal 
        refSpecCalGrid = spline(lambdaCal)

        # get derivative to lambda cal
        refSpecCalDeriv = spline.derivative()(lambdaCal)

        return refSpecCalGrid, refSpecCalDeriv

    
    def diagnostic(self):
        """
         A "diagnostic" method that returns a dictionary with
         additional diagnostic information; cf. OE.py
        """

        d = OrderedDict()
        d['modelCalculation'] = self.modelCalculation
        d['Jacobian'] = self.Jacobian
        return d

    
# ==================================================================


def irrCal(irr,irrError,irrWvl,cfg,IFDOERefSpec,groundPixel,verboseLevel):
    ''' 
    Solar spectrum (irradiance) wavelength calibration using OE. Some settings are hard-coded for now.
    '''

    irrCalConfig = {}
   
    irrCalConfig['polynomial_coefs'] = 3
    irrCalConfig['trace_gas_list'] = []

    irrCalConfig['nr_trace_gases'] = len(irrCalConfig['trace_gas_list'])
    irrCalConfig['ring_fit_term'] = None


    irrCalConfig['fit_window'] = np.copy(cfg['wlcal_window'])
    
    # # extend fit window with 1.5nm on both sides to prevent extrapolation
    # if cfg['wlcal_window_extend']:
    #     irrCalConfig['fit_window'][0] -=1.5
    #     irrCalConfig['fit_window'][1] +=1.5


    irrCalConfig['wvl_cal_offset'] = True 
    irrCalConfig['wvl_cal_stretch'] = np.copy(cfg['wlcal_solar_stretch'])

    
    # clip irradiance to wvl cal window using nnearest neighbour

    iLow = (np.abs(irrWvl - irrCalConfig['fit_window'][0])).argmin()
    iHigh = (np.abs(irrWvl - irrCalConfig['fit_window'][1])).argmin()

    irrFit = irr[iLow:iHigh+1]
    irrErrorFit = irrError[iLow:iHigh+1]
    irrWvlFit = irrWvl[iLow:iHigh+1]


    # Assemble the list of parameters
    if cfg['wlcal_solar_stretch']:
        irrCalParameterNames = ['P0','P1','P2','ws','wq']

    else:
        irrCalParameterNames = ['P0','P1','P2','ws']

    irrCalPrior = np.asarray([cfg['prior']['wvl_cal'][key][0] for key in irrCalParameterNames])
    irrCalPriorCov = np.diag(np.asarray([cfg['prior']['wvl_cal'][key][1] for key in irrCalParameterNames])**2)

    # Assemble the model details
    irrCalModel = IFDOEcalibration(prior=irrCalPrior, priorCov=irrCalPriorCov, 
                  otherModelParam=None, 
                  parameterNames=irrCalParameterNames,
                  verbose=verboseLevel-1,
                  observation=irrFit, 
                  observationError=irrErrorFit,
                  independentVariable=irrWvlFit,
                  cfg=irrCalConfig,
                  IFDOErefspec=IFDOERefSpec,
                  groundPixel=groundPixel)
    
    
    # Perform the fit       
    irrCalOE = OE(model=irrCalModel, maxiter=cfg['max_iterations'], stateVectorConvThreshold=cfg['convergence_threshold'])
    irrCalOE()

    irrCalConverged = irrCalOE.answer['converged']

    # Get wavelength shift and apply to solar grid (only fit window)
    irrCalWs = irrCalModel.stateVector[3]
    irrCalWsSigma = np.sqrt(irrCalModel.covariance[3,3])

    irrCalChiSquare = irrCalOE.costFunction

    if cfg['wlcal_solar_stretch']:
        irrCalWq = irrCalModel.stateVector[-1]
        irrCalWqSigma = np.sqrt(irrCalModel.covariance[-1,-1])
        lambdaCentre = irrCalConfig['fit_window'][0] + (irrCalConfig['fit_window'][1]-irrCalConfig['fit_window'][0])/2
        lambdaStar = 2*(irrWvl - (lambdaCentre) ) / (irrCalConfig['fit_window'][1] - irrCalConfig['fit_window'][0])
        irrWvlCal = irrWvl + irrCalWs + irrCalWq*lambdaStar
        log.debug(('irr cal. (ws, ws_sigma, wq, wq_sigma): ',irrCalWs,irrCalWsSigma,irrCalWq,irrCalWqSigma))
        return irrWvlCal, irrCalWs, irrCalWsSigma, irrCalWq, irrCalWqSigma, irrCalChiSquare, irrCalConverged


    else:
        irrWvlCal = irrWvl + irrCalWs
        log.debug(('irr cal. (ws, ws_sigma): ',irrCalWs,irrCalWsSigma))
        return irrWvlCal, irrCalWs, irrCalWsSigma, np.nan, np.nan, irrCalChiSquare, irrCalConverged



def radCal(rad,radError,radWvl,cfg,IFDOERefSpec,groundPixel,verboseLevel):
    ''' 
    Earth spectrum (radiance) wavelength calibration using OE. Some settings are hard-coded for now.
    '''

    radCalConfig = {}
   
    radCalConfig['polynomial_coefs'] = 3
    radCalConfig['trace_gas_list'] = []

    radCalConfig['nr_trace_gases'] = len(radCalConfig['trace_gas_list'])
    radCalConfig['ring_fit_term'] = cfg['ring_fit_term']

    radCalConfig['fit_window'] = np.copy(cfg['wlcal_window'])

    # # extend fit window with 1.5nm on both sides to prevent extrapolation
    # if cfg['wlcal_window_extend']:
    #     radCalConfig['fit_window'][0] -= 1.5
    #     radCalConfig['fit_window'][1] += 1.5

    radCalConfig['wvl_cal_offset'] = True 
    radCalConfig['wvl_cal_stretch'] = np.copy(cfg['wlcal_earth_stretch'])

    # # clip radiance to  wvl cal fit window

    # radFit = rad[(radWvl>=radCalConfig['fit_window'][0]) & (radWvl<=radCalConfig['fit_window'][1])]
    # radErrorFit = radError[(radWvl>=radCalConfig['fit_window'][0]) & (radWvl<=radCalConfig['fit_window'][1])]
    # radWvlFit = radWvl[(radWvl>=radCalConfig['fit_window'][0]) & (radWvl<=radCalConfig['fit_window'][1])]

    # clip radiance to wvl cal window using nnearest neighbour

    iLow = (np.abs(radWvl - radCalConfig['fit_window'][0])).argmin()
    iHigh = (np.abs(radWvl - radCalConfig['fit_window'][1])).argmin()

    radFit = rad[iLow:iHigh+1]
    radErrorFit = radError[iLow:iHigh+1]
    radWvlFit = radWvl[iLow:iHigh+1]

    # radFit = rad
    # radErrorFit = radError
    # radWvlFit = radWvl

    # Assemble the list of parameters
    if cfg['wlcal_earth_stretch']:
        radCalParameterNames = ['P0','P1','P2','ws','wq']
    else:
        radCalParameterNames = ['P0','P1','P2','ws']

    if radCalConfig['ring_fit_term'] == 'Iring':
        radCalParameterNames.append('Cring')
    elif radCalConfig['ring_fit_term'] == 'sigma':
        radCalParameterNames.append('Dring')


    radCalPrior = np.asarray([cfg['prior']['wvl_cal'][key][0] for key in radCalParameterNames])
    radCalPriorCov = np.diag(np.asarray([cfg['prior']['wvl_cal'][key][1] for key in radCalParameterNames])**2)


    # Assemble the model details
    radCalModel = IFDOEcalibration(prior=radCalPrior, priorCov=radCalPriorCov, 
                  otherModelParam=None, 
                  parameterNames=radCalParameterNames,
                  verbose=verboseLevel-1,
                  observation=radFit, 
                  observationError=radErrorFit,
                  independentVariable=radWvlFit,
                  cfg=radCalConfig,
                  IFDOErefspec=IFDOERefSpec,
                  groundPixel=groundPixel)
    
    
    # Perform the fit       
    radCalOE = OE(model=radCalModel, maxiter=cfg['max_iterations'], stateVectorConvThreshold=cfg['convergence_threshold'])
    radCalOE()

    radCalConverged = radCalOE.answer['converged']
    # print(radCalOE.answer['state vector convergence criterium'])


    # Get wavelength stretch and apply to wavelength grid


    radCalWs = radCalModel.stateVector[3]
    radCalWsSigma = np.sqrt(radCalModel.covariance[3,3])

    radCalChiSquare = radCalOE.costFunction



    if cfg['wlcal_earth_stretch']:
        radCalWq = radCalModel.stateVector[4]
        radCalWqSigma = np.sqrt(radCalModel.covariance[4,4])
        lambdaCentre = radCalConfig['fit_window'][0] + (radCalConfig['fit_window'][1]-radCalConfig['fit_window'][0])/2
        lambdaStar = 2*(radWvl - (lambdaCentre) ) / (radCalConfig['fit_window'][1] - radCalConfig['fit_window'][0])
        radWvlCal = radWvl + radCalWs + radCalWq*lambdaStar
        log.debug(('rad cal. (ws, ws_sigma, wq, wq_sigma): ',radCalWs,radCalWsSigma,radCalWq,radCalWqSigma))
        return radWvlCal, radCalWs, radCalWsSigma, radCalWq, radCalWqSigma, radCalChiSquare, radCalConverged


    else:
        radWvlCal = radWvl + radCalWs
        log.debug(('rad cal. (ws, ws_sigma): ',radCalWs,radCalWsSigma))
        return radWvlCal, radCalWs, radCalWsSigma, np.nan, np.nan, radCalChiSquare, radCalConverged


def IFDOEparameters(cfg):
    '''
     Assemble the list of parameters
    '''

    parameterNames = []

    for ipoly in range(0,cfg['polynomial_coefs'],1):
        parameterNames.append('P'+str(ipoly))

    for igas in range(0,cfg['nr_trace_gases'],1):
        # parameterNames.append(TGlist[igas])
        parameterNames.append(cfg['trace_gas_list'][igas])
        
        

    if ( cfg['ring_fit_term'] == 'Iring' ):
        parameterNames.append('Cring')
    elif ( cfg['ring_fit_term'] == 'sigma' ):
        parameterNames.append('Dring')
    else:
        None

    for ipoly in range(0,cfg['intensity_coefs'],1):
        parameterNames.append('C'+str(ipoly))

    log.info('Parameter names: %s' % parameterNames )
    
    return parameterNames


def ReadRefSpecAsc(cfg, groundPixel):
    '''
    Read in ascii ref.spec. for all ground pixels
    '''
    IFDOERefSpec = {}
    IFDOERefSpecWvl = {}

    # read in gas ref.spec.
    for gas in cfg['trace_gas_list']:
        IFDOERefSpecWvl[gas], IFDOERefSpec[gas] = ReadRefSpecAscii(cfg['ref_xs_{}'.format(gas.lower())], groundPixel)

     # read in solar ref.spec.
    IFDOERefSpecWvl['solar'], IFDOERefSpec['solar'] = ReadRefSpecAscii(cfg['ref_solar'], groundPixel)

     # read in differential ring spectrum if linear ring fit is included
    if cfg['ring_fit_term'] == 'sigma':
        IFDOERefSpecWvl['Dring'], IFDOERefSpec['Dring'] = ReadRefSpecAscii(cfg['ref_xs_ring'], groundPixel)

     # read in ring radiance spectrum if non-linear ring fit is included
    if cfg['ring_fit_term'] == 'Iring':
        IFDOERefSpecWvl['Cring'], IFDOERefSpec['Cring'] = ReadRefSpecAscii(cfg['ref_ring_radiance'], groundPixel)


    # check if wavelength grids are identical for ref.spec.
    for gas in IFDOERefSpecWvl:
        if IFDOERefSpecWvl[gas].shape[0] != IFDOERefSpecWvl['solar'].shape[-1]:
            print( ' *** Warning: Ref.spec. for {} not on same wavelength grid as solar ref.spec.'.format(gas) )
            # sys.exit(1)

        # add dimension for groundpixel: simplifies interpolation code

        IFDOERefSpec[gas] = IFDOERefSpec[gas][np.newaxis,...]


    # replace dict with common wvl grid
    IFDOERefSpec['wvl'] = IFDOERefSpecWvl['solar']


    # add wvl keys to dict for wvl calibration
    IFDOERefSpec['solarWvl']=IFDOERefSpecWvl['solar']

    IFDOERefSpec['xsWvl']=IFDOERefSpecWvl['NO2']

    return IFDOERefSpec


def ReadRefSpecAscii(FileName,groundPixel):
    '''
     Read QDOAS reference spectrum ASCII file
     rows are wavelengths, column is groundpixel
     first column in wavelength value
    '''

    data = np.genfromtxt(FileName)

    wvl = data[:,0]

    NPixel = data.shape[1]

    if groundPixel > NPixel:        
        print( ' *** Error: groundPixel {0} out of range. RefSpec {} has max. {} groundpixels'.format(groundPixel,FileName, NPixel) )
        sys.exit(1) 


    RefSpec = data[:,1:].T

    RefSpec = data[:,groundPixel+1].T


    return wvl, RefSpec


def ScaleWvl(Wvl,FitWindow):
    '''
     Rescale wvl to [-1:+1] for fit window
    '''

    ScaledWvl = (Wvl - FitWindow[0]) / (FitWindow[1] - FitWindow[0]) * 2 - 1

    return ScaledWvl

def DeScaleWvl(ScaledWvl,FitWindow):
    '''
     Return Scaled Wvl to normal. Inverse of ScaleWvl() function.
    '''

    Wvl = (ScaledWvl + 1) / 2 * (FitWindow[1] - FitWindow[0]) + FitWindow[0]

    return Wvl


def AddConfig(cfg):

    cfg['nr_trace_gases'] = len(cfg['trace_gas_list'])

    return cfg


def getDimensionsRad(radFile):
    # Get scanline and groundpixel dimensions from radiance file

    with nc.Dataset(radFile) as src:

        scanN = src.dimensions['along_track'].size
        pxlN = src.dimensions['across_track'].size
        spectralN = src.dimensions['wavelength'].size

    return scanN,pxlN, spectralN

def getDimensionsAtm(atmFile):
    # Get scanline and groundpixel dimensions from atm file

    with nc.Dataset(atmFile) as src:

        scanN = src.dimensions['along_track'].size
        pxlN = src.dimensions['across_track'].size

    return scanN,pxlN

def ReadIrrSGM(file, pxlN):
    '''
    Read irradiance from SGM file
    '''

    with nc.Dataset(file) as f:

        irr = f['solar_irradiance'][:] # [spectral_bins] - "photons / (nm m2 s)"

        # photons to mol
        irr /= constants.NA


        irrWvl = f['wavelength'][:] # [spectral_bins] - nm


    irr  = np.repeat(irr[np.newaxis,:], pxlN, axis=0)
    irrWvl  = np.repeat(irrWvl[np.newaxis,:], pxlN, axis=0)


    irrError = irr.copy() * 1e-6 # fixed noise value, otherwise OE does not work
    irrFlag = np.zeros_like(irr)

    return irrWvl, irr, irrError, irrFlag



def ReadRadSGM(f, iscan):
    '''
    Read radiance from SGM file. file is already open for thread safety
    '''
    # with nc.Dataset(file) as f:

    rad = f['radiance'][iscan,:,:] # [alt, act, spectral_bins] - "photons/(sr nm m2 s)

    # photons to mol
    rad /= constants.NA

    radWvl = f['wavelength'][:] # [spectral_bins] - nm

    radWvl  = np.repeat(radWvl[np.newaxis,:], rad.shape[0], axis=0)

    radError = rad.copy() * 1e-6 # fixed noise value, otherwise OE does not work
    radFlag = np.zeros_like(rad)

    return radWvl, rad, radError, radFlag


def ReadRad(f, iscan):
    '''
    Read radiance from L1B file. file is already open for thread safety
    '''

    rad = f['observation_data/radiance'][iscan,:,:] # [alt, act, spectral_bins] - "photons/(sr nm m2 s)
    # radFlag = f['observation_data/radiance_mask'][iscan,:,:] # [alt, act, spectral_bins] -  0 = good, 1 = bad
    radError = f['observation_data/radiance_stdev'][iscan,:,:] # [alt, act, spectral_bins] - 

    radWvl = f['observation_data/wavelength'][:,:] # [act, spectral_bins] - nm

    # photons to mol
    rad /= constants.NA
    radError /= constants.NA

    radFlag = np.zeros(rad.shape) # TODO: placeholder

    return radWvl, rad, radError, radFlag

# ======================

def SolarToEarthWvl_simple(radWvlScaled,irrWvlScaled,irr,irrError,ipxl):
    # Interpolate solar grid to earth grid

    splineIrr = interpolate.InterpolatedUnivariateSpline(irrWvlScaled[ipxl,:],irr[ipxl,:], ext=1, k=4)
    irrRegridded = splineIrr(radWvlScaled)

    splineIrrError = interpolate.InterpolatedUnivariateSpline(irrWvlScaled[ipxl,:],irrError[ipxl,:], ext=1, k=4)
    irrErrorRegridded = splineIrrError(radWvlScaled)

    return irrRegridded, irrErrorRegridded

# ======================

def SolarToEarthWvl(IFDOERefSpec,radWvlScaled,irrWvlScaled,irr,irrError,ipxl):
    # Interpolate solar grid to earth grid using ref_solar


    splineTmp = interpolate.InterpolatedUnivariateSpline(IFDOERefSpec['solarWvlScaled'],IFDOERefSpec['solar'][0,:], ext=1, k=4)
    RefSolarRad = splineTmp(radWvlScaled)
    RefSolarIrr = splineTmp(irrWvlScaled[ipxl,:])

    irrRegridded = RefSolarRad/RefSolarIrr * irr[ipxl,:]
    irrErrorRegridded = RefSolarRad/RefSolarIrr * irrError[ipxl,:]


    return irrRegridded, irrErrorRegridded

# ======================

def InterpRefSpecGrid(IFDOERefSpec,commonWvl,IFDOEconfig,ipxl):
    # Interpolate ref.spec grids to common wvl grid

    IFDOERefSpecRegrid = {}

    for RefSpec in IFDOERefSpec:
        if 'wvl' in RefSpec.lower() or 'spline' in RefSpec.lower():
            continue

        if RefSpec in IFDOEconfig['trace_gas_list']:
            IFDOERefSpecRegrid[RefSpec] = SplineInterp1D(IFDOERefSpec['xsWvlScaled'],IFDOERefSpec[RefSpec][ipxl,:],commonWvl)
        else:
            IFDOERefSpecRegrid[RefSpec] = SplineInterp1D(IFDOERefSpec['solarWvlScaled'],IFDOERefSpec[RefSpec][ipxl,:],commonWvl)


    return IFDOERefSpecRegrid

def SplineInterp1D(x,y,xnew):
    '''
    1D spline interpolation for y on grid x to new grid xnew resulting in ynew.
    Extrapolation disabled, returns zero for out of bounds.
    '''

    splineTmp = interpolate.InterpolatedUnivariateSpline(x, y, ext=1, k=4)

    # splineTmp = interpolate.CubicSpline(x,y, extrapolate=False)

    ynew = splineTmp(xnew)

    return ynew

def initOutput(l2_file,geo):
    # initialise output

    with nc.Dataset(l2_file, 'w', format='NETCDF4') as dst:

        # dst.DOAS_config = str(IFDOEconfig)
        dst.processing_date = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')

        # create dims
        alt_dim = dst.createDimension('scanline', geo['lat'].shape[0])
        act_dim = dst.createDimension('ground_pixel', geo['lat'].shape[1])

        _ = writevariablefromname(dst, 'latitude', ('scanline','ground_pixel'), geo['lat'])
        _ = writevariablefromname(dst, 'longitude', ('scanline','ground_pixel'), geo['lon'])
        _ = writevariablefromname(dst, 'solarzenithangle', ('scanline','ground_pixel'), geo['sza'])
        _ = writevariablefromname(dst, 'viewingzenithangle', ('scanline','ground_pixel'), geo['vza'])
        _ = writevariablefromname(dst, 'solarazimuthangle', ('scanline','ground_pixel'), geo['saa'])
        _ = writevariablefromname(dst, 'viewingazimuthangle', ('scanline','ground_pixel'), geo['vaa'])

        newgroup = dst.createGroup('doas')

    return

def writeOutput(l2_file,IFDOEconfig,parameterNames,results,geo):
    # write results to output file

    dst = nc.Dataset(l2_file, 'a', format='NETCDF4')

    group = 'doas/' 

    if IFDOEconfig['fit_window'][1] < 477.0:
        newgroup3 = dst[group].createGroup('no2')
        group += 'no2/' 
    elif IFDOEconfig['fit_window'][1] > 477.0:
        newgroup3 = dst[group].createGroup('o2o2')
        group += 'o2o2/' 
        
    # convert np.nan to netcdf fill value

    intList = ['nIter','errorFlag']

    for key in results:
        if key in intList:
            results[key][np.isnan(results[key])] = -2147483647
            results[key] = results[key].astype(int)

        elif results[key].dtype == float:
            results[key][np.isnan(results[key])] = 9.96921E36

    # replace quality flag
    ifdoeError = dst[group].createVariable('ifdoe_error', int, ('scanline','ground_pixel'), fill_value=-2147483647)
    ifdoeError[:,:]   =   results['errorFlag']
    ifdoeError.comment = 'Simplified error flag from ifdoe program: 0 is no error, 1 is error'
    # ifdoeError.coordinates = coord_string

    ##### add spectra to netCDF output

    if IFDOEconfig['export_spectra']:
        # add wvl dim

        wvl_dim = dst[group].createDimension('spectral_channel', results['wvl'].shape[-1])

        # write wvl var
        wvl_nc = dst[group].createVariable('wavelength', float,  ('scanline','ground_pixel','spectral_channel'), fill_value=9.96921E36, zlib=True)
        wvl_nc.units = 'nm'
        wvl_nc[:,:,:] = results['wvl']

        # write spectra

        RModel_nc = dst[group].createVariable('R_modeled', float, ('scanline','ground_pixel','spectral_channel'), fill_value=9.96921E36, zlib=True)
        RModel_nc.units = 'mol.m-2.nm-1.sr-1.s-1'
        RModel_nc[:,:,:] = results['R_model']

        RMeas_nc = dst[group].createVariable('R_measured', float, ('scanline','ground_pixel','spectral_channel'), fill_value=9.96921E36, zlib=True)
        RMeas_nc.units = 'mol.m-2.nm-1.sr-1.s-1'
        RMeas_nc[:,:,:] = results['R_meas']

        RRes_nc = dst[group].createVariable('R_residual', float, ('scanline','ground_pixel','spectral_channel'), fill_value=9.96921E36, zlib=True)
        RRes_nc.units = 'mol.m-2.nm-1.sr-1.s-1'
        RRes_nc[:,:,:] = results['R_res']


    def write_field(dictname,fieldname,units=None,i=None):


        if results[dictname].ndim == 1:
            var = dst[group].createVariable(fieldname, float, ('ground_pixel'), fill_value=9.96921E36)
            var[:] =   results[dictname]

        elif results[dictname].ndim == 2 and i==None:
            var = dst[group].createVariable(fieldname, float, ('scanline','ground_pixel'), fill_value=9.96921E36)
            # var.coordinates = coord_string
            var[:,:] =   results[dictname]

        elif results[dictname].ndim == 2 and i!=None:
            if i==0:
                if  'polynomial_coefficients' in fieldname:
                    poly_dim = 'polynomial_exponents'

                elif 'intensity_offset_coefficients' in fieldname:
                    poly_dim = 'intensity_offset_polynomial_exponents'
                
                var = dst[group].createVariable(fieldname, float, ('scanline','ground_pixel',poly_dim), fill_value=9.96921E36)
                # var.coordinates = coord_string

            else:
                var = dst[group+fieldname]
            var[:,:,i] =   results[dictname]

        if units:
            var.units = units
        else:
            var.units = '-'


        return

    # very sophisticated way of writing output:

    if 'NO2' in parameterNames:
        write_field('NO2','nitrogendioxide_slant_column_density',units="mol m-2")
        write_field('NO2Sigma','nitrogendioxide_slant_column_density_precision',units="mol m-2")

        if IFDOEconfig['fit_window'][1] < 475.0:
            write_field('NO2Geo','nitrogendioxide_geometric_column_density',units="mol m-2")
            write_field('NO2GeoSigma','nitrogendioxide_geometric_column_density_precision',units="mol m-2")
        
            write_field('amfGeo','air_mass_factor_geometric')

    if 'O2O2' in parameterNames:
        write_field('O2O2','oxygen_oxygen_dimer_slant_column_density',units="mol2 m-5")
        write_field('O2O2Sigma','oxygen_oxygen_dimer_slant_column_density_precision',units="mol2 m-5")
    
    if 'O3' in parameterNames:
        write_field('O3','ozone_slant_column_density',units="mol m-2")
        write_field('O3Sigma','ozone_slant_column_density_precision',units="mol m-2")

    if 'H2Oliq' in parameterNames:
        write_field('H2Oliq','water_liquid_slant_column_density',units="m")
        write_field('H2OliqSigma','water_liquid_slant_column_density_precision',units="m")

    if 'H2Ovap' in parameterNames:
        write_field('H2Ovap','water_slant_column_density',units="mol m-2")
        write_field('H2OvapSigma','water_slant_column_density_precision',units="mol m-2")

    if 'Cring' in parameterNames:
        write_field('Cring','ring_coefficient')
        write_field('CringSigma','ring_coefficient_precision')

    if 'P0' in parameterNames:
        try:
            poly_dim = dst[group].createDimension('polynomial_exponents',IFDOEconfig['polynomial_coefs'] )
        except:
            pass
        for i in range(IFDOEconfig['polynomial_coefs']):
            write_field('P'+str(i),'polynomial_coefficients',i=i)
            write_field('P'+str(i)+'Sigma','polynomial_coefficients_precision',i=i)

    if 'C0' in parameterNames:
        try:
            intensity_dim = dst[group].createDimension('intensity_offset_polynomial_exponents',IFDOEconfig['intensity_coefs'] )
        except:
            pass
        for i in range(IFDOEconfig['intensity_coefs']):
            write_field('C'+str(i),'intensity_offset_coefficients',i=i)
            write_field('C'+str(i)+'Sigma','intensity_offset_coefficients_precision',i=i)

    write_field('RMS','root_mean_square_error_of_fit')
    write_field('DOF','degrees_of_freedom')
    write_field('nIter','number_of_iterations')
    write_field('chiSquare','chi_square')

    write_field('irrCalWs','wavelength_calibration_irradiance_offset',units="nm")
    write_field('irrCalWsSigma','wavelength_calibration_irradiance_offset_precision',units="nm")
    write_field('irrCalChiSquare','wavelength_calibration_irradiance_chi_square')

    write_field('radCalWs','wavelength_calibration_offset',units="nm")
    write_field('radCalWsSigma','wavelength_calibration_offset_precision',units="nm")
    write_field('radCalWq','wavelength_calibration_stretch')
    write_field('radCalWqSigma','wavelength_calibration_stretch_precision')
    write_field('radCalChiSquare','wavelength_calibration_chi_square')

    if IFDOEconfig['fit_window'][1] < 477.0:
        R_no2 = dst[group].createVariable('NO2_continuum_reflectance_at_reference_wavelength', float, ('scanline','ground_pixel'), fill_value=9.96921E36)
        R_no2.reference_wavelength = 440.0
        # R_no2.coordinates = coord_string
        R_no2[:,:] = results['R_NO2']

        R_no2_precision = dst[group].createVariable('NO2_continuum_reflectance_at_reference_wavelength_precision', float, ('scanline','ground_pixel'), fill_value=9.96921E36)
        R_no2_precision.reference_wavelength = 440.0
        # R_no2_precision.coordinates = coord_string
        R_no2_precision[:,:] = results['R_NO2_precision']


    elif IFDOEconfig['fit_window'][1] > 477.0:
        R_o2o2 = dst[group].createVariable('O2O2_continuum_reflectance_at_reference_wavelength', float, ('scanline','ground_pixel'), fill_value=9.96921E36)
        R_o2o2.reference_wavelength = 477.0
        # R_o2o2.coordinates = coord_string
        R_o2o2[:,:] = results['R_O2O2']

        R_o2o2_precision = dst[group].createVariable('O2O2_continuum_reflectance_at_reference_wavelength_precision', float, ('scanline','ground_pixel'), fill_value=9.96921E36)
        R_o2o2_precision.reference_wavelength = 477.0
        # R_o2o2_precision.coordinates = coord_string
        R_o2o2_precision[:,:] = results['R_O2O2_precision']
        
    dst.close()
    
    return

def readGeometryL1b(rad_file, slice_alt, slice_act):
    # Read geometry from L1B
    geo = {}
    with nc.Dataset(rad_file) as f:
        geo['sza'] = f['geolocation_data/solarzenithangle'][slice_alt,slice_act]
        geo['saa'] = f['geolocation_data/solarazimuthangle'][slice_alt,slice_act]
        geo['vza'] = f['geolocation_data/viewingzenithangle'][slice_alt,slice_act]
        geo['vaa'] = f['geolocation_data/viewingazimuthangle'][slice_alt,slice_act]
        geo['lat'] = f['geolocation_data/latitude'][slice_alt,slice_act]
        geo['lon'] = f['geolocation_data/longitude'][slice_alt,slice_act]
    return geo


def readGeometryGm(gm_file, slice_alt, slice_act):
    # Read geometry from GM file
    geo = {}
    with nc.Dataset(gm_file) as f:
        geo['sza'] = f['solarzenithangle'][slice_alt,slice_act]
        geo['saa'] = f['solarazimuthangle'][slice_alt,slice_act]
        geo['vza'] = f['viewingzenithangle'][slice_alt,slice_act]
        geo['vaa'] = f['viewingazimuthangle'][slice_alt,slice_act]
        geo['lat'] = f['latitude'][slice_alt,slice_act]
        geo['lon'] = f['longitude'][slice_alt,slice_act]
    return geo