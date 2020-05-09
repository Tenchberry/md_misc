import simtk.openmm as mm
import simtk.unit as units
import numpy as np
import time

kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA

class VelocityVerletIntegrator(mm.CustomIntegrator):

    def __init__(self, timestep=1.0 * units.femtoseconds, RemoveComVel=False):
        """Construct a velocity Verlet integrator.
        Parameters
        ----------
        timestep : numpy.unit.Quantity compatible with femtoseconds, default: 1*simtk.unit.femtoseconds
           The integration timestep.
        """

        super(VelocityVerletIntegrator, self).__init__(timestep)

        # Add PerDof Variable
        self.addPerDofVariable("x1", 0)
        self.addPerDofVariable("vint", 0)
        self.addPerDofVariable("positions", 0)
        self.addPerDofVariable("velocities", 0)
        self.addGlobalVariable("nstep", 0)
        self.addGlobalVariable("Epot", 0)
        self.addGlobalVariable("Ekin", 0)

        # Update the state of the Context
        self.addUpdateContextState()

        if RemoveComVel:
            # Remove COM rotation from the first step velocities
            # This prevents the system from starting to rotate
            self.beginIfBlock("nstep = 0")
            self.addRemoveCom("vint")
            self.addComputePerDof("v", "vint")
            self.endBlock()

        self.addComputeGlobal("Epot", "energy")
        self.addComputeSum("Ekin", "0.5*m*v*v")

        # Velocity Verlet Agorithm
        self.addComputePerDof("v", "v+0.5*dt*f/m")
        self.addComputePerDof("x", "x+dt*v")
        self.addComputePerDof("x1", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
        self.addConstrainVelocities()

        self.addComputePerDof("positions", "x") 
        self.addComputePerDof("velocities", "v")
        self.addComputeGlobal("nstep", "nstep + 1")
        

    @property
    def getx(self):
        return self.getPerDofVariableByName("positions")

    @property
    def getv(self):
        return self.getPerDofVariableByName("velocities")

    @property
    def get_vint(self):
        return self.getPerDofVariableByName("vint")

    @property
    def get_Epot(self):
        return self.getGlobalVariableByName("Epot")

    @property
    def get_Ekin(self):
        return self.getGlobalVariableByName("Ekin")


class FilterApplication(mm.CustomIntegrator):

    def __init__(self, coeff, timestep=1.0 * units.femtoseconds, RemoveComVel=False):
        """Construct a velocity Verlet integrator.
        Parameters
        ----------
        timestep : numpy.unit.Quantity compatible with femtoseconds, default: 1*simtk.unit.femtoseconds
           The integration timestep.
        """

        super(FilterApplication, self).__init__(timestep)

        # Add PerDof Variable
        self.addPerDofVariable("x1", 0)
        self.addPerDofVariable("vint", 0)
        # Filter related variables
        self.addPerDofVariable("vfilt", 0)
        self.addPerDofVariable("x0", 0)
        self.addPerDofVariable("v0", 0)
        self.addGlobalVariable("c", 0)

        # Update the state of the Context
        self.addUpdateContextState()

        if RemoveComVel:
            # Remove COM rotation from the first step velocities
            # This prevents the system from starting to rotate
            self.addRemoveCom("vint")
            self.addComputePerDof("v", "vint")


        Nbuffer= len(coeff)
        self.addComputePerDof("vfilt", "0")
        self.addComputePerDof("x0", "x")
        self.addComputePerDof("v0", "v")

        for i in range(Nbuffer):
            self.MD(1)
            if RemoveComVel:
                self.addRemoveCom("vint")
            else:
                self.addComputePerDof("vint","v")
            self.addComputeGlobal("c", str(coeff[i]))
            self.addComputePerDof("vfilt", "vfilt + vint*c")

        self.addComputePerDof("x", "x0")
        self.addComputePerDof("v", "v0")


    def MD(self, nsteps):
        """
            Velocity Verlet integrator
        """
        for step in range(nsteps):
            self.addComputePerDof("v", "v+0.5*dt*f/m")
            self.addComputePerDof("x", "x+dt*v")
            self.addComputePerDof("x1", "x")
            self.addConstrainPositions()
            self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
            self.addConstrainVelocities()


    @property
    def getvfilt(self):
        return self.getPerDofVariableByName('vfilt')


class DFMDIntegrator(mm.CustomIntegrator):

    def __init__(self, coeff, filterFreq, timestep=1*units.femtoseconds):
        """Construct a DFMD Integrator
        Parameters
        ----------
        coeff       : filter coefficients
        filterFreq  : Frequency to apply the Digital Filter
        timestep    : The instegration timestep

        """
        super(DFMDIntegrator, self).__init__(timestep)
        
        # Variables for Velocity Verlet Algorithm
        self.addPerDofVariable("x1", 0)  # for constraints
        
        # Variables for Filter Application
        self.addPerDofVariable("x_mid", 0)
        self.addPerDofVariable("vfilt", 0)
        self.addPerDofVariable("vint", 0)
        self.addGlobalVariable("c", 0)
        
        # Variables for controlling the frequency to apply Filter Application 
        self.addGlobalVariable("nstep", 0)
        self.addGlobalVariable("filterFreq", filterFreq)

        self.addComputeGlobal("nstep", "nstep + 1")
        
        self.addUpdateContextState()

        self.MD() # write NVT integrator 
        
        # If freq of DF apply DF
        self.beginIfBlock("floor(nstep/filterFreq)*filterFreq = nstep")
        self.DFMD(coeff)
        self.endBlock()

    def DFMD(self, coeff):
        """
            Filter application Algorithm for DFMD
            Note: This will reset the position to the middle of buffer.
            and reset the velocities to the filtered velocities.
        """
        Nbuffer = len(coeff)
        self.addComputePerDof("vfilt", "0")
        
        for i in xrange(Nbuffer):
            self.MD()
            self.addRemoveCom("vint")
            self.addComputeGlobal("c", str(coeff[i]))
            self.addComputePerDof("vfilt", "vfilt + vint*c")
            if i == (Nbuffer-1)/2:
                self.addComputePerDof("x_mid", "x")

        # Restore state to the middle of buffer
        self.addComputePerDof("x","x_mid")
        self.addComputePerDof("v","vfilt")

    def MD(self):
        """
            Velocity Verlet algorithm
            with Position and velocity constraints
        """
        self.addComputePerDof("v", "v+0.5*dt*f/m")
        self.addComputePerDof("x", "x+dt*v")
        self.addComputePerDof("x1", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
        self.addConstrainVelocities()

class RDFMDIntegrator(mm.CustomIntegrator):

    def __init__(self, coeff, filterFreq, delay=0, timestep=1*units.femtoseconds):
        """Construct a RDFMD Integrator
        Parameters
        ----------
        coeff       : filter coefficients
        filterFreq  : Frequency to apply the Digital Filter
        delay       : 
        timestep    : The instegration timestep

        """
        super(RDFMDIntegrator, self).__init__(timestep)
        
        # Variables for Velocity Verlet Algorithm
        self.addPerDofVariable("x1", 0)  # for constraints
        
        self.addPerDofVariable("x_mid", 0)
        self.addPerDofVariable("x_0", 0)
        self.addPerDofVariable("v_0", 0)
        self.addPerDofVariable("x_out", 0)
        self.addPerDofVariable("v_out", 0)
        self.addPerDofVariable("vfilt", 0)
        self.addPerDofVariable("vint", 0)

        self.addGlobalVariable("c", 0)
        
        # Variables for controlling the frequency to apply Filter Application 
        self.addGlobalVariable("nstep", 0)
        self.addGlobalVariable("filterFreq", filterFreq)

        self.addComputeGlobal("nstep", "nstep + 1")
        
        self.addUpdateContextState()

        self.MD() # write NVT integrator 
        
        # If freq of DF apply DF
        self.beginIfBlock("floor(nstep/filterFreq)*filterFreq = nstep")
        self.RDFMD(coeff, 0)
        self.endBlock()

    def RDFMD(self, coeff, delay):
        Nbuffer = len(coeff)
        # initialize filter velocities to zero.
        self.addComputePerDof("vfilt", "0")
        mid = (Nbuffer-1)/2
        Nforward = mid + delay
        Nbackward = mid - delay
        coeff_forward = coeff[Nbackward:]
        coeff_backward = coeff[:Nbackward]

        #Save the current state
        self.addComputePerDof("x_0", "x")
        self.addComputePerDof("v_0", "v")

        #Forward MD
        self.addRemoveCom("vint")
        self.addComputeGlobal("c", str(coeff_forward[0]))
        self.addComputePerDof("vfilt", "vfilt + vint*c")
        for i in xrange(Nforward):
            if delay >= 0 and i==delay:
                self.addComputePerDof("x_mid", "x")
            self.MD()
            self.addRemoveCom("vint")
            self.addComputeGlobal("c", str(coeff_forward[i+1]))
            self.addComputePerDof("vfilt", "vfilt + vint*c")

        #Backward MD
        # Negate the velocities to perform reverse MD
        self.addComputePerDof("v", "-v_0")
        self.addComputePerDof("x", "x_0")
        
        for j in xrange(Nbackward):
            if delay < 0 and i==-delay:
                self.addComputePerDof("x_mid", "x")
            self.MD()
            self.addRemoveCom("vint")
            self.addComputeGlobal("c", str(coeff_backward[j]))
            self.addComputePerDof("vfilt", "vfilt + vint*c")

        # Restore state to the middle of buffer
        self.addComputePerDof("x","x_mid")
        self.addComputePerDof("v","vfilt")

    def MD(self):
        """
            Velocity Verlet algorithm
            with Position and velocity constraints
        """
        self.addComputePerDof("v", "v+0.5*dt*f/m")
        self.addComputePerDof("x", "x+dt*v")
        self.addComputePerDof("x1", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
        self.addConstrainVelocities()


class HMCIntegrator(mm.CustomIntegrator):
    # This works can be used

    """
    Hybrid Monte Carlo (HMC) integrator.
    """

    def __init__(self, temperature=298.0 * units.kelvin, MDStep=10, timestep=1 * units.femtoseconds, RemoveComVel=False):
        """
        Create a hybrid Monte Carlo (HMC) integrator.
        Parameters
        ----------
        temperature : numpy.unit.Quantity compatible with kelvin, default: 298*simtk.unit.kelvin
           The temperature.
        MDSteps : int, default: 10
           The number of velocity Verlet steps to take per HMC trial.
        timestep : numpy.unit.Quantity compatible with femtoseconds, default: 1*simtk.unit.femtoseconds
           The integration timestep.
        Warning
        -------
        Because 'nsteps' sets the number of steps taken, a call to integrator.step(1) actually takes 'nsteps' steps.
        Notes
        -----
        The velocity is drawn from a Maxwell-Boltzmann distribution, then 'nsteps' steps are taken,
        and the new configuration is either accepted or rejected.
        Additional global variables 'ntrials' and  'naccept' keep track of how many trials have been attempted and
        accepted, respectively.
        TODO
        ----
        Currently, the simulation timestep is only advanced by 'timestep' each step, rather than timestep*nsteps.  Fix this.
        Examples
        --------
        Create an HMC integrator.
        >>> timestep = 1.0 * simtk.unit.femtoseconds # fictitious timestep
        >>> temperature = 298.0 * simtk.unit.kelvin
        >>> MDSteps = 10 # number of steps per call
        >>> integrator = HMCIntegrator(temperature, MDSteps, timestep)
        """

        super(HMCIntegrator, self).__init__(timestep)

        # Compute the thermal energy.
        kT = kB * temperature

        #
        # Integrator initialization.
        #
        self.addGlobalVariable("naccept", 0)  # number accepted
        self.addGlobalVariable("ntrials", 0)  # number of Metropolization trials

        self.addGlobalVariable("kT", kT)  # thermal energy
        self.addPerDofVariable("sigma", 0)
        self.addGlobalVariable("ke", 0)  # kinetic energy
        self.addPerDofVariable("xold", 0)  # old positions
        self.addGlobalVariable("Eold", 0)  # old energy
        self.addGlobalVariable("Enew", 0)  # new energy
        self.addGlobalVariable("accept", 0)  # accept or reject
        self.addPerDofVariable("x1", 0)  # for constraints
        self.addPerDofVariable("vint", 0)  # to calculate internal velocity

        #
        # Pre-computation.
        # This only needs to be done once, but it needs to be done for each degree of freedom.
        # Could move this to initialization?
        #
        self.addComputePerDof("sigma", "sqrt(kT/m)")

        #
        # Allow Context updating here, outside of inner loop only.
        #
        self.addUpdateContextState()

        #
        # Draw new velocity.
        # from Normal distribution 
        # with std=sigma and mean=0
        #

        self.addComputePerDof("v", "sigma*gaussian")
        self.addConstrainVelocities()

        if RemoveComVel:
            # Remove COM rotation from the drawn velocities
            # This prevents the system to start rotating
            self.addRemoveCom("vint")
            self.addComputePerDof("v", "vint")

        #
        # Store old position and energy.
        #
        self.addComputeSum("ke", "0.5*m*v*v")
        self.addComputeGlobal("Eold", "ke + energy")
        self.addComputePerDof("xold", "x")

        #
        # Inner symplectic steps using velocity Verlet.
        #
        self.MD(MDStep)

        #
        # Store position and energy before the MD steps starts
        #
        self.addComputeSum("ke", "0.5*m*v*v")
        self.addComputeGlobal("Enew", "ke + energy")

        #
        # Accept/reject step.
        #
        self.addComputeGlobal("accept", "step(exp(-(Enew-Eold)/kT) - uniform)")
        self.addComputePerDof("x", "x*accept + xold*(1-accept)")

        #
        # Accumulate statistics.
        #
        self.addComputeGlobal("naccept", "naccept + accept")
        self.addComputeGlobal("ntrials", "ntrials + 1")

    def MD(self, nsteps):
        for step in range(nsteps):
            self.addComputePerDof("v", "v+0.5*dt*f/m")
            self.addComputePerDof("x", "x+dt*v")
            self.addComputePerDof("x1", "x")
            self.addConstrainPositions()
            self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
            self.addConstrainVelocities()

    @property
    def n_accept(self):
        """The number of accepted HMC moves."""
        return self.getGlobalVariableByName("naccept")

    @property
    def n_trials(self):
        """The total number of attempted HMC moves."""
        return self.getGlobalVariableByName("ntrials")

    @property
    def acceptance_rate(self):
        """The acceptance rate: n_accept  / n_trials."""
        return 100*(self.n_accept / float(self.n_trials))


class TestDFHMC(mm.CustomIntegrator):
    # The filter appliction works as expected
    # try the enhancement!!
    def __init__(self, coeff, filterUpdate, atoms=None, scl=1.0, temperature=300.0 * units.kelvin, MDStep=10, timestep=1 * units.femtoseconds):

        super(TestDFHMC, self).__init__(timestep)

        # Compute the thermal energy.
        kT = kB * temperature

        # input variables
        self.addGlobalVariable("filterUpdate", filterUpdate)

        # MC statistics and out put handling
        self.addGlobalVariable("naccept", 0)  # number accepted
        self.addGlobalVariable("ntrials", 0)  # number of Metropolization trials
        self.addPerDofVariable("vfilt_n", 0)    # the last step velocity of the DF application
        self.addPerDofVariable("v_bimodal", 0)

        # variables related to drawing velocities
        self.addGlobalVariable("kT", kT)  # thermal energy
        self.addGlobalVariable("scl", scl)
        self.addPerDofVariable("sigma", 0)
        self.addPerDofVariable("vfilt", 0)
        self.addPerDofVariable("b", 0)
        self.addGlobalVariable("a", 0) # two sided coin ouput 0 or 1 with 0.5 prob

        # variables related to the acceptance criteria
        self.addGlobalVariable("ke", 0)  # kinetic energy
        self.addPerDofVariable("xold", 0)  # old positions
        self.addGlobalVariable("Eold", 0)  # old total energy
        self.addGlobalVariable("Enew", 0)  # new total energy
        self.addGlobalVariable("accept", 0)  # accept or reject
        
        
        self.addGlobalVariable("bias",1)
        self.addGlobalVariable("bvold",0)
        self.addGlobalVariable("bvnew",0)

        # Variables used in applying the filter
        self.addGlobalVariable("c", 0)
        self.addPerDofVariable("x0", 0)
        self.addPerDofVariable("v0", 0)
        self.addPerDofVariable("x_mid", 0)
        self.addPerDofVariable("atoms", 1)
        self.addPerDofVariable("vint",0)
        # self.addGlobalVariable("penalty", 1)
        self.addPerDofVariable("testdof1",0)    
        self.addPerDofVariable("testdof2",0)    

        # MD integrator related variables
        self.addPerDofVariable("x1", 0)  # for constraints

        if atoms is not None:
            self.setPerDofVariableByName("atoms", atoms)
       
        self.addComputePerDof("sigma", "sqrt(kT/m)")

        # Draw velocities from Normal distro with mean 0 and std=sigma
        self.addComputePerDof("v", "sigma*gaussian")
        self.addConstrainVelocities()

        # only update b when step is a mod of filterupdate
        self.beginIfBlock('floor(ntrials/filterUpdate)*filterUpdate = ntrials ')
        self.applyDF(coeff)
        self.addComputePerDof("b", "atoms*vfilt*scl")
        self.endBlock()

        #
        # Draw new velocity
        #

        # Add/Substract the bias b to the velocities to obtain a bimodal distribution
        self.addComputeGlobal("a","step(0.5-uniform)")
        self.addComputePerDof("v","v+b*a + (1-a)*(-b)")
        #
        # Store position and energy before the MD steps starts
        #
        self.addComputeSum("ke", "0.5*m*v*v")
        self.addComputeSum("bvold", "m*b*v")
        self.addComputeGlobal("Eold", "ke + energy")
        self.addComputePerDof("xold", "x")
        #
        # Inner symplectic steps using velocity Verlet.
        #
        self.addComputePerDof("testdof2","x")
        self.MD(MDStep)
        self.addComputePerDof("testdof1","x")    

        #
        # Store position and energy before the MD steps starts
        #
        self.addComputeSum("ke", "0.5*m*v*v")
        self.addComputeSum("bvnew", "m*b*v")
        self.addComputeGlobal("Enew", "ke + energy")

        #
        # Accept/reject step.
        #
        
        self.addComputeGlobal("bias","(exp(-bvnew/kT)+exp(bvnew/kT))/(exp(-bvold/kT)+exp(bvold/kT))")
        self.addComputeGlobal("accept", "step(exp(-(Enew-Eold)/kT)*bias - uniform)")
        # if accept=1 accept the current state
        # if accept=0 reject the state and hence set the current state to xold
        self.addComputePerDof("x", "x*accept + xold*(1-accept)")
        
        #
        # Accumulate statistics.
        #
        self.addComputeGlobal("naccept", "naccept + accept")
        self.addComputeGlobal("ntrials", "ntrials + 1")

    def applyDF(self, coeff):
        """
            Calculate vfilt
            by taking the linear combination
            of the each velocities with the given coefficients
        """

        # TODO: Remove com motion from velocities and position
        # before filter aplication
        # Something is going on here investigate. the forces at the end of the Df changes 
        # when run one after other
        # also the forces before the application fdoes not match after the application

        Nbuffer= len(coeff)
        self.addComputePerDof("vfilt", "0")
        self.addComputePerDof("x0", "x")
        self.addComputePerDof("v0", "v")

        for i in range(Nbuffer):
            self.MD(1)
            self.addRemoveCom("vint")
            self.addComputeGlobal("c", str(coeff[i]))
            self.addComputePerDof("vfilt", "vfilt + vint*c")
            if i == (Nbuffer-1)/2:
                self.addComputePerDof("vfilt_n","v")

        self.addComputePerDof("x", "x0")
        self.addComputePerDof("v", "v0")

    def applyRDF(self, coeff, delay=0):
        Nbuffer = len(coeff)
        mid = (Nbuffer-1)/2
        Nforward = mid + delay
        Nbackward = mid - delay
        coeff_forward = coeff[Nbackward:]
        coeff_backward = coeff[:Nbackward]

        #Save the current state
        self.addComputePerDof("vfilt", "0")
        self.addComputePerDof("x0", "x")
        self.addComputePerDof("v0", "v")

        #Forward MD
        self.addComputeGlobal("c", str(coeff_forward[0]))
        self.addComputePerDof("vfilt", "vfilt + v*c")
        for i in xrange(Nforward):
            if delay >= 0 and i==delay:
                self.addComputePerDof("x_mid", "x")
            self.MD(1)
            self.addComputeGlobal("c", str(coeff_forward[i+1]))
            self.addComputePerDof("vfilt", "vfilt + v*c")

        #Backward MD
        # Negate the velocities to perform reverse MD
        self.addComputePerDof("v", "-v0")
        self.addComputePerDof("x", "x0")

        for j in xrange(Nbackward):
            if delay < 0 and i==-delay:
                self.addComputePerDof("x_mid", "x")
            self.MD(1)
            self.addComputeGlobal("c", str(coeff_backward[j]))
            self.addComputePerDof("vfilt", "vfilt - v*c")

        # Restore state to the middle of buffer
        self.addComputePerDof("x","x_mid")
        self.addComputePerDof("v","vfilt")

    def MD(self, nsteps):
        """
            Velocity Verlet integrator
        """
        for step in range(nsteps):
            self.addComputePerDof("v", "v+0.5*dt*f/m")
            self.addComputePerDof("x", "x+dt*v")
            self.addComputePerDof("x1", "x")
            self.addConstrainPositions()
            self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
            self.addConstrainVelocities()


    @property
    def testdof1(self):
        """Get the bias velocity b= vfilt*scl*atoms"""
        return self.getPerDofVariableByName("testdof1")

    @property
    def testdof2(self):
        """Get the bias velocity b= vfilt*scl*atoms"""
        return self.getPerDofVariableByName("testdof2")

    @property
    def getb(self):
        """Get the bias velocity b= vfilt*scl*atoms"""
        return self.getPerDofVariableByName("b")

    @property
    def getbias(self):
        """The acceptance criteria of the bias part"""
        return self.getGlobalVariableByName("bias")

    @property
    def getbvnew(self):
        """the exponential term of the bias b*v*m*"""
        return self.getGlobalVariableByName("bvnew")

    @property
    def getbvold(self):
        return self.getGlobalVariableByName("bvold")

    @property
    def x_old(self):
        """The positions before the MDSteps"""
        return self.getPerDofVariableByName("xold")

    @property
    def vfilt_n(self):
        """The velocity of the mid point of the buffer"""
        return self.getPerDofVariableByName("vfilt_n")

    @property
    def v_bimodal(self):
        """The velocity drawn from bimodal distro"""
        return self.getPerDofVariableByName("v_bimodal")

    @property
    def getscl(self):
        """The filter scaling factor"""
        return self.getGlobalVariableByName("scl")

    @property
    def vfilt(self):
        """The velocity drawn from bimodal distro"""
        return self.getPerDofVariableByName("vfilt")

    @property
    def n_accept(self):
        """The number of accepted HMC moves."""
        return self.getGlobalVariableByName("naccept")

    @property
    def n_trials(self):
        """The total number of attempted HMC moves."""
        return self.getGlobalVariableByName("ntrials")

    @property
    def acceptance_rate(self):
        """The acceptance rate: n_accept  / n_trials."""
        return 100*self.n_accept / float(self.n_trials)


class DFHMCIntegrator(mm.CustomIntegrator):
    """
    The bias is updated every nstfilter step. But the calculated b value is used till the next update.
    Thus if the molecule rotates during the nstfilter step, the b value will not be inline with the positions.

    """
    def __init__(self, coeff, filterUpdate, atoms=None, scl=1.0, temperature=300.0*units.kelvin,
                 MDStep=10, timestep=1*units.femtoseconds, RemoveComVel=False):

        super(DFHMCIntegrator, self).__init__(timestep)

        # Compute the thermal energy.
        kT = kB * temperature

        # input variables
        self.addGlobalVariable("filterUpdate", filterUpdate)

        # MC statistics and out put handling
        self.addGlobalVariable("naccept", 0)  # number accepted
        self.addGlobalVariable("ntrials", 0)  # number of Metropolization trials

        # variables related to drawing velocities
        self.addGlobalVariable("kT", kT)  # thermal energy
        self.addGlobalVariable("scl", scl)
        self.addPerDofVariable("sigma", 0)
        self.addPerDofVariable("vfilt", 0)
        self.addPerDofVariable("b", 0)
        self.addGlobalVariable("a", 0) # two sided coin ouput 0 or 1 with 0.5 prob

        # variables related to the acceptance criteria
        self.addGlobalVariable("ke", 0)  # kinetic energy
        self.addPerDofVariable("xold", 0)  # old positions
        self.addGlobalVariable("Eold", 0)  # old total energy
        self.addGlobalVariable("Enew", 0)  # new total energy
        self.addGlobalVariable("accept", 0)  # accept or reject
        
        
        self.addGlobalVariable("bias",1)
        self.addGlobalVariable("bvold",0)
        self.addGlobalVariable("bvnew",0)

        # Variables used in applying the filter
        self.addGlobalVariable("c", 0)
        self.addPerDofVariable("x0", 0)
        self.addPerDofVariable("v0", 0)
        self.addPerDofVariable("atoms", 1)
        self.addPerDofVariable("vint",0)

        # MD integrator related variables
        self.addPerDofVariable("x1", 0)  # for constraints

        if atoms is not None:
            self.setPerDofVariableByName("atoms", atoms)
       
        self.addComputePerDof("sigma", "sqrt(kT/m)")

        # Draw velocities from Normal distro with mean 0 and std=sigma
        self.addComputePerDof("v", "sigma*gaussian")

        # only update b when step is a mod of filterupdate
        self.beginIfBlock('floor(ntrials/filterUpdate)*filterUpdate = ntrials ')
        self.applyDF(coeff)
        self.addComputePerDof("b", "atoms*vfilt*scl")
        self.endBlock()

        # Add/Substract the bias b to the velocities to obtain a bimodal distribution
        self.addComputeGlobal("a","step(0.5-uniform)")
        self.addComputePerDof("v","v+b*a + (1-a)*(-b)")
        self.addConstrainVelocities()


        if RemoveComVel:
            # Remove COM rotation from the drawn velocities
            # This prevents the system to start rotating
            self.addRemoveCom("vint")
            self.addComputePerDof("v", "vint")

        #
        # Store position and energy before the MD steps starts
        #
        self.addComputeSum("ke", "0.5*m*v*v")
        self.addComputeSum("bvold", "m*b*v")
        self.addComputeGlobal("Eold", "ke + energy")
        self.addComputePerDof("xold", "x")
        
        #
        # Inner symplectic steps using velocity Verlet.
        #
        self.MD(MDStep)

        #
        # Store energy after the MD steps
        #
        self.addComputeSum("ke", "0.5*m*v*v")
        self.addComputeSum("bvnew", "m*b*v")
        self.addComputeGlobal("Enew", "ke + energy")

        #
        # Accept/reject step.
        #
        
        self.addComputeGlobal("bias","(exp(-bvnew/kT)+exp(bvnew/kT))/(exp(-bvold/kT)+exp(bvold/kT))")
        self.addComputeGlobal("accept", "step(exp(-(Enew-Eold)/kT)*bias - uniform)")
        # if accept=1 accept the current state
        # if accept=0 reject the state and hence set the current state to xold
        self.addComputePerDof("x", "x*accept + xold*(1-accept)")
        
        #
        # Accumulate statistics.
        #
        self.addComputeGlobal("naccept", "naccept + accept")
        self.addComputeGlobal("ntrials", "ntrials + 1")

    def applyDF(self, coeff):
        """
            Calculate vfilt
            by taking the linear combination
            of the each velocities with the given coefficients
        """

        Nbuffer= len(coeff)
        self.addComputePerDof("vfilt", "0")
        self.addComputePerDof("x0", "x")
        self.addComputePerDof("v0", "v")

        for i in range(Nbuffer):
            self.MD(1)
            self.addRemoveCom("vint")
            self.addComputeGlobal("c", str(coeff[i]))
            self.addComputePerDof("vfilt", "vfilt + vint*c")

        self.addComputePerDof("x", "x0")
        self.addComputePerDof("v", "v0")

    def MD(self, nsteps):
        """
            Velocity Verlet integrator
        """
        for step in range(nsteps):
            self.addComputePerDof("v", "v+0.5*dt*f/m")
            self.addComputePerDof("x", "x+dt*v")
            self.addComputePerDof("x1", "x")
            self.addConstrainPositions()
            self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
            self.addConstrainVelocities()

    @property
    def getb(self):
        """Get the bias velocity b= vfilt*scl*atoms"""
        return self.getPerDofVariableByName("b")

    @property
    def n_accept(self):
        """The number of accepted HMC moves."""
        return self.getGlobalVariableByName("naccept")

    @property
    def n_trials(self):
        """The total number of attempted HMC moves."""
        return self.getGlobalVariableByName("ntrials")

    @property
    def acceptance_rate(self):
        """The acceptance rate: n_accept  / n_trials."""
        return 100*self.n_accept / float(self.n_trials)


class DFHMC_DWIntegrator(mm.CustomIntegrator):
    # The filter appliction works as expected
    # try the enhancement!!
    def __init__(self, coeff, filterUpdate, atoms=None, scl=1.0, temperature=300.0 * units.kelvin, MDStep=10, timestep=1 * units.femtoseconds):

        super(DFHMC_DWIntegrator, self).__init__(timestep)

        # Compute the thermal energy.
        kT = kB * temperature

        # input variables
        self.addGlobalVariable("filterUpdate", filterUpdate)

        # MC statistics and out put handling
        self.addGlobalVariable("naccept", 0)  # number accepted
        self.addGlobalVariable("ntrials", 0)  # number of Metropolization trials

        # variables related to drawing velocities
        self.addGlobalVariable("kT", kT)  # thermal energy
        self.addGlobalVariable("scl", scl)
        self.addPerDofVariable("sigma", 0)
        self.addPerDofVariable("vfilt", 0)
        self.addPerDofVariable("b", 0)
        self.addGlobalVariable("a", 0) # two sided coin ouput 0 or 1 with 0.5 prob

        # variables related to the acceptance criteria
        self.addGlobalVariable("ke", 0)  # kinetic energy
        self.addPerDofVariable("xold", 0)  # old positions
        self.addGlobalVariable("Eold", 0)  # old total energy
        self.addGlobalVariable("Enew", 0)  # new total energy
        self.addGlobalVariable("accept", 0)  # accept or reject
        
        
        self.addGlobalVariable("bias",1)
        self.addGlobalVariable("bvold",0)
        self.addGlobalVariable("bvnew",0)

        # Variables used in applying the filter
        self.addGlobalVariable("c", 0)
        self.addPerDofVariable("x0", 0)
        self.addPerDofVariable("v0", 0)
        self.addPerDofVariable("atoms", 1)
        self.addPerDofVariable("vint",0)
        self.addPerDofVariable("vout",0)

        # MD integrator related variables
        self.addPerDofVariable("x1", 0)  # for constraints

        if atoms is not None:
            self.setPerDofVariableByName("atoms", atoms)
       
        self.addComputePerDof("sigma", "sqrt(kT/m)")

        # Draw velocities from Normal distro with mean 0 and std=sigma
        self.addComputePerDof("v", "sigma*gaussian")

        # only update b when step is a mod of filterupdate
        self.beginIfBlock('floor(ntrials/filterUpdate)*filterUpdate = ntrials ')
        self.applyDF(coeff)
        self.addComputePerDof("b", "atoms*vfilt*scl")
        self.endBlock()

        #
        # Draw new velocity
        #

        # Add/Substract the bias b to the velocities to obtain a bimodal distribution
        self.addComputeGlobal("a","step(0.5-uniform)")
        self.addComputePerDof("v","v+b*a + (1-a)*(-b)")
        self.addConstrainVelocities()
        self.addComputePerDof("vout","v")
        #
        # Store position and energy before the MD steps starts
        #
        self.addComputeSum("ke", "0.5*m*v*v")
        self.addComputeSum("bvold", "m*b*v")
        self.addComputeGlobal("Eold", "ke + energy")
        self.addComputePerDof("xold", "x")
        #
        # Inner symplectic steps using velocity Verlet.
        #
        self.MD(MDStep)

        #
        # Store energy after the MD steps
        #
        self.addComputeSum("ke", "0.5*m*v*v")
        self.addComputeSum("bvnew", "m*b*v")
        self.addComputeGlobal("Enew", "ke + energy")

        #
        # Accept/reject step.
        #
        
        self.addComputeGlobal("bias","(exp(-bvnew/kT)+exp(bvnew/kT))/(exp(-bvold/kT)+exp(bvold/kT))")
        self.addComputeGlobal("accept", "step(exp(-(Enew-Eold)/kT)*bias - uniform)")
        # if accept=1 accept the current state
        # if accept=0 reject the state and hence set the current state to xold
        self.addComputePerDof("x", "x*accept + xold*(1-accept)")
        
        #
        # Accumulate statistics.
        #
        self.addComputeGlobal("naccept", "naccept + accept")
        self.addComputeGlobal("ntrials", "ntrials + 1")

    def applyDF(self, coeff):
        """
            Calculate vfilt
            by taking the linear combination
            of the each velocities with the given coefficients
        """

        Nbuffer= len(coeff)
        self.addComputePerDof("vfilt", "0")
        self.addComputePerDof("x0", "x")
        self.addComputePerDof("v0", "v")

        for i in range(Nbuffer):
            self.MD(1)
            self.addComputeGlobal("c", str(coeff[i]))
            self.addComputePerDof("vfilt", "vfilt + v*c")

        self.addComputePerDof("x", "x0")
        self.addComputePerDof("v", "v0")

    def MD(self, nsteps):
        """
            Velocity Verlet integrator
        """
        for step in range(nsteps):
            self.addComputePerDof("v", "v+0.5*dt*f/m")
            self.addComputePerDof("x", "x+dt*v")
            self.addComputePerDof("x1", "x")
            self.addConstrainPositions()
            self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
            self.addConstrainVelocities()


    @property
    def getsigma(self):
        """Get the bias velocity b= vfilt*scl*atoms"""
        return self.getPerDofVariableByName("sigma")

    @property
    def vout(self):
        """Get the bias velocity b= vfilt*scl*atoms"""
        return self.getPerDofVariableByName("vout")

    @property
    def getb(self):
        """Get the bias velocity b= vfilt*scl*atoms"""
        return self.getPerDofVariableByName("b")

    @property
    def getbias(self):
        """The acceptance criteria of the bias part"""
        return self.getGlobalVariableByName("bias")

    @property
    def getbvnew(self):
        """the exponential term of the bias b*v*m*"""
        return self.getGlobalVariableByName("bvnew")

    @property
    def getbvold(self):
        return self.getGlobalVariableByName("bvold")

    @property
    def x_old(self):
        """The positions before the MDSteps"""
        return self.getPerDofVariableByName("xold")

    @property
    def getscl(self):
        """The filter scaling factor"""
        return self.getGlobalVariableByName("scl")

    @property
    def vfilt(self):
        """The velocity drawn from bimodal distro"""
        return self.getPerDofVariableByName("vfilt")

    @property
    def n_accept(self):
        """The number of accepted HMC moves."""
        return self.getGlobalVariableByName("naccept")

    @property
    def n_trials(self):
        """The total number of attempted HMC moves."""
        return self.getGlobalVariableByName("ntrials")

    @property
    def acceptance_rate(self):
        """The acceptance rate: n_accept  / n_trials."""
        return 100*self.n_accept / float(self.n_trials)


class MEHMCIntegrator(mm.CustomIntegrator):
    
    def __init__(self, L=14, temperature=298.0 * units.kelvin, MDSteps=10, timestep=1*units.femtoseconds):
        super(MEHMCIntegrator, self).__init__(timestep)

        # Compute the thermal energy.
        kT = kB * temperature

        #
        # Integrator initialization.
        #
        self.addGlobalVariable("naccept", 0)  # number accepted
        self.addGlobalVariable("ntrials", 0)  # number of Metropolization trials

        # Recursive Filter
        self.addPerDofVariable("vfilt", 0)
        self.addGlobalVariable("L", L)

        self.addGlobalVariable("kT", kT)  # thermal energy
        self.addPerDofVariable("sigma", 0)
        self.addGlobalVariable("ke", 0)  # kinetic energy
        self.addPerDofVariable("xold", 0)  # old positions
        self.addGlobalVariable("Eold", 0)  # old energy
        self.addGlobalVariable("Enew", 0)  # new energy
        self.addGlobalVariable("accept", 0)  # accept or reject
        self.addPerDofVariable("x1", 0)  # for constraints

        #
        # Pre-computation.
        # This only needs to be done once, but it needs to be done for each degree of freedom.
        # Could move this to initialization?
        #
        self.addComputePerDof("sigma", "sqrt(kT/m)")

        #
        # Allow Context updating here, outside of inner loop only.
        #
        self.addUpdateContextState()

        #
        # Draw new velocity.
        #
        self.addComputePerDof("v", "sigma*gaussian")
        self.addConstrainVelocities()

        #
        # Store old position and energy.
        #
        self.addComputeSum("ke", "0.5*m*v*v")
        self.addComputeGlobal("Eold", "ke + energy")
        self.addComputePerDof("xold", "x")

        #
        # Inner symplectic steps using velocity Verlet.
        #
        self.MD(MDSteps)

        #
        # Store position and energy before the MD steps starts
        #
        self.addComputeSum("ke", "0.5*m*v*v")
        self.addComputeGlobal("Enew", "ke + energy")

        #
        # Accept/reject step.
        #
        self.addComputeGlobal("accept", "step(exp(-(Enew-Eold)/kT) - uniform)")
        self.addComputePerDof("x", "x*accept + xold*(1-accept)")

        #
        # Accumulate statistics.
        #
        self.addComputeGlobal("naccept", "naccept + accept")
        self.addComputeGlobal("ntrials", "ntrials + 1")

    def MD(self, nsteps):
        """
            Velocity Verlet integrator
        """
        for step in range(nsteps):
            self.addComputePerDof("v", "v+0.5*dt*f/m")
            self.addComputePerDof("x", "x+dt*v")
            self.addComputePerDof("x1", "x")
            self.addConstrainPositions()
            self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
            self.addConstrainVelocities()
            self.addComputePerDof("vfilt", "(1.0-1.0/L)*vfilt + (1.0/L)*v")

    @property
    def vfilt(self):
        """Recursive filter result"""
        return self.getPerDofVariableByName("vfilt")

    @property
    def n_accept(self):
        """The number of accepted HMC moves."""
        return self.getGlobalVariableByName("naccept")

    @property
    def n_trials(self):
        """The total number of attempted HMC moves."""
        return self.getGlobalVariableByName("ntrials")

    @property
    def acceptance_rate(self):
        """The acceptance rate: n_accept  / n_trials."""
        return 100*self.n_accept / float(self.n_trials)

