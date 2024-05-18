## Written by Abdennour Harche
## Mars 2024

# This is a full package of methods for VQE for fermions in a Pairing Interaction

# Import Packages
from qiskit_aer.primitives import Estimator as AerEstimator
from qiskit.primitives import Estimator
from qiskit.circuit.library import RealAmplitudes, EfficientSU2
from qiskit.circuit import Parameter
from qiskit.quantum_info import SparsePauliOp
import numpy as np
from qiskit import QuantumCircuit
from scipy.optimize import minimize, Bounds
import matplotlib.pyplot as plt

# Define Functions
def Ep(p) -> float:
    counter = 0
    for n in range(100):
        for j in range((n+1)*(n+2)//2):
            if p==counter: return (n+3/2)
            counter += 1

def textmutate(txt:str,index:int,sub:str):
    t = list(txt)
    t[index] = sub
    return "".join(t)

def create_hamiltonian(nstates:int, pairingG:float) -> SparsePauliOp:
    iDEN = 'I'*nstates
    terms = list()
    coeff = list()
    G = pairingG
    for p in range(nstates//2):
        term = [iDEN]*11
        
        terms.append(term[0])
        coeff.append(Ep(p))

        term[1] = textmutate(term[1],2*p,'Z')
        term[1] = term[1][::-1]
        terms.append(term[1])
        coeff.append(-0.5*Ep(p))

        term[2] = textmutate(term[2],2*p+1,'Z')
        term[2] = term[2][::-1]
        terms.append(term[2])
        coeff.append(-0.5*Ep(p))

        for k in range(p+1,nstates//2):
            term[3] = textmutate(term[3],2*p,'Y')
            term[3] = textmutate(term[3],2*p+1,'X')
            term[3] = textmutate(term[3],2*k,'Y')
            term[3] = textmutate(term[3],2*k+1,'X')
            term[3] = term[3][::-1]
            terms.append(term[3])
            coeff.append(-G/8)
            term[3] = iDEN

            term[4] = textmutate(term[4],2*p,'Y')
            term[4] = textmutate(term[4],2*p+1,'X')
            term[4] = textmutate(term[4],2*k,'X')
            term[4] = textmutate(term[4],2*k+1,'Y')
            term[4] = term[4][::-1]
            terms.append(term[4])
            coeff.append(-G/8)
            term[4] = iDEN

            term[5] = textmutate(term[5],2*p,'Y')
            term[5] = textmutate(term[5],2*p+1,'Y')
            term[5] = textmutate(term[5],2*k,'Y')
            term[5] = textmutate(term[5],2*k+1,'Y')
            term[5] = term[5][::-1]
            terms.append(term[5])
            coeff.append(-G/8)
            term[5] = iDEN

            term[6] = textmutate(term[6],2*p,'Y')
            term[6] = textmutate(term[6],2*p+1,'Y')
            term[6] = textmutate(term[6],2*k,'X')
            term[6] = textmutate(term[6],2*k+1,'X')
            term[6] = term[6][::-1]
            terms.append(term[6])
            coeff.append(G/8)
            term[6] = iDEN

            term[7] = textmutate(term[7],2*p,'X')
            term[7] = textmutate(term[7],2*p+1,'X')
            term[7] = textmutate(term[7],2*k,'Y')
            term[7] = textmutate(term[7],2*k+1,'Y')
            term[7] = term[7][::-1]
            terms.append(term[7])
            coeff.append(G/8)
            term[7] = iDEN
            
            term[8] = textmutate(term[8],2*p,'X')
            term[8] = textmutate(term[8],2*p+1,'X')
            term[8] = textmutate(term[8],2*k,'X')
            term[8] = textmutate(term[8],2*k+1,'X')
            term[8] = term[8][::-1]
            terms.append(term[8])
            coeff.append(-G/8)
            term[8] = iDEN

            term[9] = textmutate(term[9],2*p,'X')
            term[9] = textmutate(term[9],2*p+1,'Y')
            term[9] = textmutate(term[9],2*k,'Y')
            term[9] = textmutate(term[9],2*k+1,'X')
            term[9] = term[9][::-1]
            terms.append(term[9])
            coeff.append(-G/8)
            term[9] = iDEN

            term[10] = textmutate(term[10],2*p,'X')
            term[10] = textmutate(term[10],2*p+1,'Y')
            term[10] = textmutate(term[10],2*k,'X')
            term[10] = textmutate(term[10],2*k+1,'Y')
            term[10] = term[10][::-1]
            terms.append(term[10])
            coeff.append(-G/8)
            term[10] = iDEN
    return SparsePauliOp(terms,coeffs=coeff)

def create_auxhamiltonian(nstates:int, pairingG:float, λ: float) -> SparsePauliOp:
    iDEN = 'I'*nstates
    terms = list()
    coeff = list()
    G = pairingG
    for p in range(nstates//2):
        term = [iDEN]*11
        enu = Ep(p) - λ
        terms.append(term[0])
        coeff.append(enu)

        term[1] = textmutate(term[1],2*p,'Z')
        term[1] = term[1][::-1]
        terms.append(term[1])
        coeff.append(-0.5*enu)

        term[2] = textmutate(term[2],2*p+1,'Z')
        term[2] = term[2][::-1]
        terms.append(term[2])
        coeff.append(-0.5*enu)

        for k in range(p+1,nstates//2):
            term[3] = textmutate(term[3],2*p,'Y')
            term[3] = textmutate(term[3],2*p+1,'X')
            term[3] = textmutate(term[3],2*k,'Y')
            term[3] = textmutate(term[3],2*k+1,'X')
            term[3] = term[3][::-1]
            terms.append(term[3])
            coeff.append(-G/8)
            term[3] = iDEN

            term[4] = textmutate(term[4],2*p,'Y')
            term[4] = textmutate(term[4],2*p+1,'X')
            term[4] = textmutate(term[4],2*k,'X')
            term[4] = textmutate(term[4],2*k+1,'Y')
            term[4] = term[4][::-1]
            terms.append(term[4])
            coeff.append(-G/8)
            term[4] = iDEN

            term[5] = textmutate(term[5],2*p,'Y')
            term[5] = textmutate(term[5],2*p+1,'Y')
            term[5] = textmutate(term[5],2*k,'Y')
            term[5] = textmutate(term[5],2*k+1,'Y')
            term[5] = term[5][::-1]
            terms.append(term[5])
            coeff.append(-G/8)
            term[5] = iDEN

            term[6] = textmutate(term[6],2*p,'Y')
            term[6] = textmutate(term[6],2*p+1,'Y')
            term[6] = textmutate(term[6],2*k,'X')
            term[6] = textmutate(term[6],2*k+1,'X')
            term[6] = term[6][::-1]
            terms.append(term[6])
            coeff.append(G/8)
            term[6] = iDEN

            term[7] = textmutate(term[7],2*p,'X')
            term[7] = textmutate(term[7],2*p+1,'X')
            term[7] = textmutate(term[7],2*k,'Y')
            term[7] = textmutate(term[7],2*k+1,'Y')
            term[7] = term[7][::-1]
            terms.append(term[7])
            coeff.append(G/8)
            term[7] = iDEN
            
            term[8] = textmutate(term[8],2*p,'X')
            term[8] = textmutate(term[8],2*p+1,'X')
            term[8] = textmutate(term[8],2*k,'X')
            term[8] = textmutate(term[8],2*k+1,'X')
            term[8] = term[8][::-1]
            terms.append(term[8])
            coeff.append(-G/8)
            term[8] = iDEN

            term[9] = textmutate(term[9],2*p,'X')
            term[9] = textmutate(term[9],2*p+1,'Y')
            term[9] = textmutate(term[9],2*k,'Y')
            term[9] = textmutate(term[9],2*k+1,'X')
            term[9] = term[9][::-1]
            terms.append(term[9])
            coeff.append(-G/8)
            term[9] = iDEN

            term[10] = textmutate(term[10],2*p,'X')
            term[10] = textmutate(term[10],2*p+1,'Y')
            term[10] = textmutate(term[10],2*k,'X')
            term[10] = textmutate(term[10],2*k+1,'Y')
            term[10] = term[10][::-1]
            terms.append(term[10])
            coeff.append(-G/8)
            term[10] = iDEN
    return SparsePauliOp(terms,coeffs=coeff)

def cost_func(params, ansatz, hamiltonian, estimator):
    """Return estimate of energy from estimator

    Parameters:
        params (ndarray): Array of ansatz parameters
        ansatz (QuantumCircuit): Parameterized ansatz circuit
        hamiltonian (SparsePauliOp): Operator representation of Hamiltonian
        estimator (Estimator): Estimator primitive instance

    Returns:
        float: Energy estimate
    """
    energy = estimator.run(ansatz, hamiltonian, parameter_values=params).result().values[0]
    return energy

def build_callback(ansatz, hamiltonian, estimator, callback_dict):
    """Return callback function that uses Estimator instance,
    and stores intermediate values into a dictionary.

    Parameters:
        ansatz (QuantumCircuit): Parameterized ansatz circuit
        hamiltonian (SparsePauliOp): Operator representation of Hamiltonian
        estimator (Estimator): Estimator primitive instance
        callback_dict (dict): Mutable dict for storing values

    Returns:
        Callable: Callback function object
    """

    def callback(current_vector):
        """Callback function storing previous solution vector,
        computing the intermediate cost value, and displaying number
        of completed iterations and average time per iteration.

        Values are stored in pre-defined 'callback_dict' dictionary.

        Parameters:
            current_vector (ndarray): Current vector of parameters
                                      returned by optimizer
        """
        # Keep track of the number of iterations
        callback_dict["iters"] += 1
        # Set the prev_vector to the latest one
        callback_dict["prev_vector"] = current_vector
        # Compute the value of the cost function at the current vector
        # This adds an additional function evaluation
        current_cost = (
            estimator.run(ansatz, hamiltonian, parameter_values=current_vector).result().values[0]
        )
        callback_dict["cost_history"].append(current_cost)
        # Print to screen on single line
        print(
            "Iters. done: {} [Current cost: {}]".format(callback_dict["iters"], current_cost),
            end="\r",
            flush=True,
        )

    return callback

def cost_func2(params, ansatz, nstates, pG, estimator):
    """Return estimate of energy from estimator

    Parameters:
        params (ndarray): Array of ansatz parameters
        ansatz (QuantumCircuit): Parameterized ansatz circuit
        nstates (SparsePauliOp): Number of states
        pg : Pairing Constant G
        estimator (Estimator): Estimator primitive instance

    Returns:
        float: Energy estimate
    """
    hamiltonian = create_auxhamiltonian(nstates, pG, params[0])
    energy = estimator.run(ansatz, hamiltonian, parameter_values=params[1:]).result().values[0]
    return energy

def build_callback2(ansatz, nstates, pG, estimator, callback_dict):
    """Return callback function that uses Estimator instance,
    and stores intermediate values into a dictionary.

    Parameters:
        ansatz (QuantumCircuit): Parameterized ansatz circuit
        nstates (SparsePauliOp): Number of states
        pg : Pairing Constant G
        estimator (Estimator): Estimator primitive instance
        callback_dict (dict): Mutable dict for storing values

    Returns:
        Callable: Callback function object
    """

    def callback(current_vector):
        """Callback function storing previous solution vector,
        computing the intermediate cost value, and displaying number
        of completed iterations and average time per iteration.

        Values are stored in pre-defined 'callback_dict' dictionary.

        Parameters:
            current_vector (ndarray): Current vector of parameters
                                      returned by optimizer
        """
        # Keep track of the number of iterations
        callback_dict["iters"] += 1
        # Set the prev_vector to the latest one
        callback_dict["prev_vector"] = current_vector
        # Compute the value of the cost function at the current vector
        # This adds an additional function evaluation
        hamiltonian = create_auxhamiltonian(nstates, pG, current_vector[0])
        npar = checkN(ansatz, estimator, create_N(nstates), current_vector[1:])
        current_cost = (
            estimator.run(ansatz, hamiltonian, parameter_values=current_vector[1:]).result().values[0]
        )
        callback_dict["cost_history"].append(current_cost)
        # Print to screen on single line
        print(
            "Iters. done: {} [Current cost: {}, number of particles: {}]".format(callback_dict["iters"], current_cost, npar),
            end="\r",
            flush=True,
        )

    return callback

def create_N(nstates:int) -> SparsePauliOp:
    iDEN = 'I'*nstates
    termsN = list()
    coeffN = list()

    for p in range(nstates//2):
        term = [iDEN]*3

        termsN.append(term[0])
        coeffN.append(1)

        term[1] = textmutate(term[1],2*p,'Z')
        #term[1][2*p]='Z'
        term[1] = term[1][::-1]
        termsN.append(term[1])
        coeffN.append(-1/2)

        term[2] = textmutate(term[2],2*p+1,'Z')
        #term[2][2*p+1]='Z'
        term[2] = term[2][::-1]
        termsN.append(term[2])
        coeffN.append(-1/2)

    return SparsePauliOp(termsN,coeffs=coeffN)

def checkN(ansatz ,estimator, numberOperator, params):
    return estimator.run(ansatz, numberOperator, parameter_values=params).result().values[0]

def VQE(hamiltonian: SparsePauliOp
        ,ansatz: QuantumCircuit
        ,numberParticles:int
        ,numberOperator:SparsePauliOp
        ,estimator, x0 = None):
    callback_dict = {
    "prev_vector": None,
    "iters": 0,
    "cost_history": [],
    }
    def cons(params):
        number = estimator.run(ansatz, numberOperator, parameter_values=params).result().values[0]
        return (number - numberParticles)
    const = {'type':'eq', 'fun': cons}
    if x0 is None: x0 = [0.1]*ansatz.num_parameters
    anglebounds = Bounds(0, 4 * np.pi)
    callback = build_callback(ansatz, hamiltonian, estimator, callback_dict)
    opt = {'disp':False,'maxiter':500}
    res = minimize(
        cost_func,
        x0,
        args=(ansatz, hamiltonian, estimator),
        method="SLSQP",
        callback=callback,
        bounds=anglebounds,
        constraints=const,
        options=opt
    )
    return res

def VQE2(nstates, pG ,ansatz: QuantumCircuit
        ,numberParticles:int
        ,numberOperator:SparsePauliOp
        ,estimator, x0 = None):
    callback_dict = {
    "prev_vector": None,
    "iters": 0,
    "cost_history": [],
    }
    def cons(params):
        number = estimator.run(ansatz, numberOperator, parameter_values=params[1:]).result().values[0]
        return (number - numberParticles)
    const = {'type':'eq', 'fun': cons}
    if x0 is None: x0 = [0.1]*(ansatz.num_parameters+1)
    anglebounds = Bounds(0, 4 * np.pi)
    callback = build_callback2(ansatz, nstates, pG, estimator, callback_dict)
    opt = {'disp':False,'maxiter':500}
    res = minimize(
        cost_func2,
        x0,
        args=(ansatz, nstates, pG, estimator),
        method="SLSQP",
        callback=callback,
        bounds=anglebounds,
        constraints=const,
        options=opt
    )
    return res

def fulljob(number_of_states:int,number_of_particles:int, G:float, ansatzDepth:int):
    psitheta = RealAmplitudes(num_qubits=number_of_states, reps=ansatzDepth)
    
    qc = QuantumCircuit(number_of_states)
    #qc.x(range(number_of_particles))
    psitheta = qc.compose(psitheta)
    H = create_hamiltonian(nstates=number_of_states, pairingG=G)
    N = create_N(number_of_states)
    estimator = Estimator()
    print(psitheta.num_parameters)
    res = VQE(hamiltonian=H, ansatz=psitheta, numberParticles=number_of_particles, numberOperator=N, estimator=estimator)
    #print(checkN(psitheta,estimator,N,res.x))
    return res

def fulljob2(number_of_states:int,number_of_particles:int, G:float, ansatzDepth:int):
    psitheta = RealAmplitudes(num_qubits=number_of_states, reps=ansatzDepth)
    
    qc = QuantumCircuit(number_of_states)
    #qc.x(range(number_of_particles))
    psitheta = qc.compose(psitheta)
    H = create_hamiltonian(nstates=number_of_states, pairingG=G)
    N = create_N(number_of_states)
    estimator = AerEstimator(approximation=True)
    print(psitheta.num_parameters)
    res = VQE(hamiltonian=H, ansatz=psitheta, numberParticles=number_of_particles, numberOperator=N, estimator=estimator)
    #print(checkN(psitheta,estimator,N,res.x))
    return res

def fulljob3(psitheta: QuantumCircuit,number_of_particles:int, G:float, x0):
    nstates = psitheta.num_qubits
    H = create_hamiltonian(nstates=nstates, pairingG=G)
    N = create_N(nstates)
    estimator = AerEstimator(approximation=True)
    print(psitheta.num_parameters)
    res = VQE(hamiltonian=H, ansatz=psitheta, numberParticles=number_of_particles, numberOperator=N, estimator=estimator, x0=x0)
    #print(checkN(psitheta,estimator,N,res.x))
    return res
def fulljob4(psitheta: QuantumCircuit,number_of_particles:int, G:float, x0):
    nstates = psitheta.num_qubits
    H = create_hamiltonian(nstates=nstates, pairingG=G).simplify()
    N = create_N(nstates)
    estimator = AerEstimator(approximation=True)
    print(psitheta.num_parameters)
    res = VQE(hamiltonian=H, ansatz=psitheta, numberParticles=number_of_particles, numberOperator=N, estimator=estimator, x0=x0)
    #print(checkN(psitheta,estimator,N,res.x))
    return res

def fulljob5(psitheta: QuantumCircuit,number_of_particles:int, G:float, x0):
    nstates = psitheta.num_qubits
    N = create_N(nstates)
    estimator = AerEstimator(approximation=True)
    print(psitheta.num_parameters)
    res = VQE2(nstates=nstates,pG=G, ansatz=psitheta, numberParticles=number_of_particles, numberOperator=N, estimator=estimator, x0=x0)
    #print(checkN(psitheta,estimator,N,res.x))
    return res

def graphingjob(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    es = list()
    for g in G:
        r = fulljob(number_of_states=number_of_states,number_of_particles=number_of_particles,G=g,ansatzDepth=ansatzDepth)
        es.append(r.fun)
    fig, ax = plt.subplots()
    ax.plot(G, es)
    ax.set_xlabel("Pairing Constant G")
    ax.set_ylabel("Average Energy E")
    plt.draw()

def excitingjob(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    delta = list()
    for g in G:
        r1 = fulljob(number_of_states=number_of_states,number_of_particles=number_of_particles,G=g,ansatzDepth=ansatzDepth)
        r2 = fulljob(number_of_states=number_of_states,number_of_particles=number_of_particles+1,G=g,ansatzDepth=ansatzDepth)
        delta.append(r2.fun-r1.fun)
    fig, ax = plt.subplots()
    ax.plot(G, delta)
    ax.set_xlabel("Pairing Constant G")
    ax.set_ylabel(r"Average Energy Difference $\Delta$")
    plt.draw()

def graphingjob2(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    e0 = 8
    es = list()
    for g in G:
        r = fulljob(number_of_states=number_of_states,number_of_particles=number_of_particles,G=g,ansatzDepth=ansatzDepth)
        es.append(r.fun)
    
    Δsquared = G * (e0 - np.array(es))
    Δ = abs(Δsquared)
    fig, ax = plt.subplots()
    ax.plot(G, Δ)
    ax.set_xlabel("Pairing Constant G")
    ax.set_ylabel("Average Energy E")
    plt.draw()


def graphingjob3(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    e0 = 8
    es = list()
    for g in G:
        r = fulljob(number_of_states=number_of_states,number_of_particles=number_of_particles,G=g,ansatzDepth=ansatzDepth)
        es.append(r.fun)
    
    Δsquared = G * (e0 - np.array(es))
    Δ = np.sqrt(abs(Δsquared))
    fig, ax = plt.subplots()
    ax.plot(G, Δ)
    ax.set_xlabel("Pairing Constant G")
    ax.set_ylabel(r"$\Delta$")
    plt.draw()

def graphingjob4(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    e0 = 8
    es = list()
    for g in G:
        r = fulljob2(number_of_states=number_of_states,number_of_particles=number_of_particles,G=g,ansatzDepth=ansatzDepth)
        es.append(r.fun)
    
    Δsquared = G * (e0 - np.array(es))
    Δ = np.sqrt(abs(Δsquared))
    fig, ax = plt.subplots()
    ax.plot(G, Δ)
    ax.set_xlabel("Pairing Constant G")
    ax.set_ylabel(r"$\Delta$")
    plt.draw()

def graphingjob4(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    e0 = 8
    es = list()
    for g in G:
        r = fulljob2(number_of_states=number_of_states,number_of_particles=number_of_particles,G=g,ansatzDepth=ansatzDepth)
        es.append(r.fun)
    
    Δsquared = G * (e0 - np.array(es))
    Δ = np.sqrt(abs(Δsquared))
    fig, ax = plt.subplots()
    ax.plot(G, Δ)
    ax.set_xlabel("Pairing Constant G")
    ax.set_ylabel(r"$\Delta$")
    plt.draw()


def graphingjob5(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    psitheta = RealAmplitudes(num_qubits=number_of_states, reps=ansatzDepth)
    
    qc = QuantumCircuit(number_of_states)
    #qc.x(range(number_of_particles))
    psitheta = qc.compose(psitheta)

    e0 = 8
    es = list()
    xx = None
    for g in G:
        r = fulljob3(psitheta=psitheta,number_of_particles=number_of_particles,G=g,ansatzDepth=ansatzDepth, x0=xx)
        es.append(r.fun)
        #xx = r.x #+ np.array([0.01]*psitheta.num_parameters)
    
    Δsquared = G * (e0 - np.array(es))
    Δ = np.sqrt(abs(Δsquared))
    fig, ax = plt.subplots()
    ax.plot(G, Δ)
    ax.set_xlabel("Pairing Constant G")
    ax.set_ylabel(r"$\Delta$")
    plt.show()

def graphingjob6(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    psitheta = RealAmplitudes(num_qubits=number_of_states, reps=ansatzDepth)
    
    qc = QuantumCircuit(number_of_states)
    #qc.x(range(number_of_particles))
    psitheta = qc.compose(psitheta)

    e0 = 3
    es = list()
    xx = None
    for g in G:
        r = fulljob3(psitheta=psitheta,number_of_particles=number_of_particles,G=g,ansatzDepth=ansatzDepth, x0=xx)
        es.append(r.fun)
        xx = r.x + np.array([0.25]*psitheta.num_parameters)
    
    Δsquared = G * (e0 - np.array(es))
    Δ = np.sqrt(abs(Δsquared))
    fig, ax = plt.subplots()
    ax.plot(G, Δ)
    ax.set_xlabel("Pairing Constant G")
    ax.set_ylabel(r"$\Delta$")
    plt.show()

def graphingjob7(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    psitheta = RealAmplitudes(num_qubits=number_of_states, reps=ansatzDepth)
    
    qc = QuantumCircuit(number_of_states)
    #qc.x(range(number_of_particles))
    psitheta = qc.compose(psitheta)

    e0 = 1.5 * 2 + 2.5 * 6
    es = list()
    xx = None
    for g in G:
        r = fulljob3(psitheta=psitheta,number_of_particles=number_of_particles,G=g,ansatzDepth=ansatzDepth, x0=xx)
        es.append(r.fun)
        xx = r.x + np.array([0.25]*psitheta.num_parameters)
    
    Δsquared = G * (e0 - np.array(es))
    Δ = np.sqrt(abs(Δsquared))
    fig, ax = plt.subplots()
    ax.plot(G, Δ)
    ax.set_xlabel("Pairing Constant G")
    ax.set_ylabel(r"$\Delta$")
    plt.show()

def graphingjob8(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    psitheta = RealAmplitudes(num_qubits=number_of_states, reps=ansatzDepth)
    
    qc = QuantumCircuit(number_of_states)
    #qc.x(range(number_of_particles))
    psitheta = qc.compose(psitheta)

    e0 = 0
    counter = number_of_particles
    for n in range(20):
        g = (n+1)*(n+2)
        if(g <= counter):
            e0 += g * (n+3/2)
            counter -= g
        else:
            e0 += counter * (n+3/2)
            break
        
        
    es = list()
    xx = None
    for g in G:
        r = fulljob3(psitheta=psitheta,number_of_particles=number_of_particles,G=g,ansatzDepth=ansatzDepth, x0=xx)
        es.append(r.fun)
        #xx = r.x + np.array([0.25]*psitheta.num_parameters)
    
    Δsquared = G * (e0 - np.array(es))
    Δ = np.sqrt(abs(Δsquared))
    fig, ax = plt.subplots()
    ax.plot(G, Δ)
    ax.set_xlabel("Pairing Constant G")
    ax.set_ylabel(r"$\Delta$")
    plt.show()

def graphingjob9(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    psitheta = RealAmplitudes(num_qubits=number_of_states, reps=ansatzDepth)
    
    qc = QuantumCircuit(number_of_states)
    #qc.x(range(number_of_particles))
    psitheta = qc.compose(psitheta)

    e0 = 0
    counter = number_of_particles
    for n in range(20):
        g = (n+1)*(n+2)
        if(g <= counter):
            e0 += g * (n+3/2)
            counter -= g
        else:
            e0 += counter * (n+3/2)
            break
        
        
    es = list()
    xx = None
    for g in G:
        r = fulljob3(psitheta=psitheta,number_of_particles=number_of_particles,G=g, x0=xx)
        es.append(r.fun)
        #xx = r.x + np.array([0.25]*psitheta.num_parameters)
    
    Δsquared = G * (e0 - np.array(es))
    Δ = np.sqrt(abs(Δsquared))
    fig, ax = plt.subplots()
    ax.plot(G, Δ)
    ax.set_xlabel("Pairing Constant G")
    ax.set_ylabel(r"$\Delta$")
    plt.savefig(f'{number_of_states}s{number_of_particles}p{G[0]}to{G[-1]}g{len(G)}s{ansatzDepth}l.png')

def storingjob(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    psitheta = RealAmplitudes(num_qubits=number_of_states, reps=ansatzDepth)
    
    qc = QuantumCircuit(number_of_states)
    #qc.x(range(number_of_particles))
    psitheta = qc.compose(psitheta)

    e0 = 0
    counter = number_of_particles
    for n in range(20):
        g = (n+1)*(n+2)
        if(g <= counter):
            e0 += g * (n+3/2)
            counter -= g
        else:
            e0 += counter * (n+3/2)
            break
        
        
    es = list()
    xx = None
    for g in G:
        r = fulljob3(psitheta=psitheta,number_of_particles=number_of_particles,G=g, x0=xx)
        es.append(r.fun)
        #xx = r.x + np.array([0.25]*psitheta.num_parameters)
    np.savetxt(f'./Energies/{number_of_states}s{number_of_particles}p{G[0]}to{G[-1]}g{len(G)}s{ansatzDepth}l.dat', es)

def storingjob2(number_of_states:int,number_of_particles:int,G,ansatzDepth:int):
    psitheta = EfficientSU2(num_qubits=number_of_states, reps=ansatzDepth)
    
    qc = QuantumCircuit(number_of_states)
    #qc.x(range(number_of_particles))
    psitheta = qc.compose(psitheta)

    e0 = 0
    counter = number_of_particles
    for n in range(20):
        g = (n+1)*(n+2)
        if(g <= counter):
            e0 += g * (n+3/2)
            counter -= g
        else:
            e0 += counter * (n+3/2)
            break
        
        
    es = list()
    xx = None
    for g in G:
        r = fulljob3(psitheta=psitheta,number_of_particles=number_of_particles,G=g, x0=xx)
        es.append(r.fun)
        #xx = r.x + np.array([0.25]*psitheta.num_parameters)
    np.savetxt(f'./SU2Energies/su{number_of_states}s{number_of_particles}p{G[0]}to{G[-1]}g{len(G)}s{ansatzDepth}l.dat', es)

def BCS(number_of_states:int) -> QuantumCircuit:
    qc = QuantumCircuit(number_of_states)
    parameters = [Parameter(fr'$\theta_{i}$') for i in range(number_of_states//2)]
    for i in range(number_of_states//2):
        qc.ry(parameters[i], 2*i)
        qc.cx(2*i,2*i+1)
    return qc

def BCS2(number_of_states:int) -> QuantumCircuit:
    qc = QuantumCircuit(number_of_states)
    parameters = [Parameter(fr'$\theta_{i}$') for i in range(number_of_states*2)]
    for i in range(number_of_states//2):
        qc.ry(parameters[4*i], 2*i)
        qc.ry(parameters[4*i+1], 2*i+1)
        qc.cx(2*i,2*i+1)
        qc.ry(parameters[4*i+2], 2*i)
        qc.ry(parameters[4*i+3], 2*i+1)
    return qc
    
def storingjob3(number_of_states:int,number_of_particles:int,G):
    psitheta = BCS(number_of_states)
    
    #qc = QuantumCircuit(number_of_states)
    # qc.x(1)
    #psitheta = qc.compose(psitheta)

    e0 = 0
    counter = number_of_particles
    for n in range(20):
        g = (n+1)*(n+2)
        if(g <= counter):
            e0 += g * (n+3/2)
            counter -= g
        else:
            e0 += counter * (n+3/2)
            break
        
        
    es = list()
    xx = np.random.rand(psitheta.num_parameters) * np.pi
    for g in G:
        r = fulljob3(psitheta=psitheta,number_of_particles=number_of_particles,G=g, x0=xx)
        es.append(r.fun)
        #xx = r.x + np.array([0.25]*psitheta.num_parameters)
    np.savetxt(f'./BCSEnergies/bcs{number_of_states}s{number_of_particles}p{G[0]}to{G[-1]}g{len(G)}s.dat', es)

def storingjob4(number_of_states:int,number_of_particles:int,G):
    psitheta = BCS(number_of_states)
    
    #qc = QuantumCircuit(number_of_states)
    # qc.x(1)
    #psitheta = qc.compose(psitheta)

    e0 = 0
    counter = number_of_particles
    for n in range(20):
        g = (n+1)*(n+2)
        if(g <= counter):
            e0 += g * (n+3/2)
            counter -= g
        else:
            e0 += counter * (n+3/2)
            break
        
        
    es = list()
    xx = np.random.rand(psitheta.num_parameters) * np.pi
    for g in G:
        r = fulljob4(psitheta=psitheta,number_of_particles=number_of_particles,G=g, x0=xx)
        es.append(r.fun)
        #xx = r.x + np.array([0.25]*psitheta.num_parameters)
    np.savetxt(f'./BCSEnergiesS/Sbcs{number_of_states}s{number_of_particles}p{G[0]}to{G[-1]}g{len(G)}s.dat', es)

def storingjob5(number_of_states:int,number_of_particles:int,G):
    psitheta = BCS(number_of_states)
    
    #qc = QuantumCircuit(number_of_states)
    # qc.x(1)
    #psitheta = qc.compose(psitheta)

    e0 = 0
    counter = number_of_particles
    for n in range(20):
        g = (n+1)*(n+2)
        if(g <= counter):
            e0 += g * (n+3/2)
            counter -= g
        else:
            e0 += counter * (n+3/2)
            break
        
        
    es = list()
    xx = None #np.random.rand(psitheta.num_parameters + 1) * np.pi
    for g in G:
        r = fulljob5(psitheta=psitheta,number_of_particles=number_of_particles,G=g, x0=xx)
        es.append(cost_func(r.x[1:], psitheta, create_hamiltonian(number_of_states, g), AerEstimator(approximation=True)))
        #xx = r.x + np.array([0.25]*psitheta.num_parameters)
    np.savetxt(f'./AuxBCS/auxbcs{number_of_states}s{number_of_particles}p{G[0]}to{G[-1]}g{len(G)}s.dat', es)