from ipie.legacy.estimators.ci import simple_fci_bose_fermi
from ipie.legacy.systems.hubbard_holstein import HubbardHolstein

def ipie_fci(w0: float, g: float, t: float, nsites: int, 
             nboson: int, pbc: bool=False, verbose: bool=False
    ) -> float:
    
    nup = 1
    ndown = 0

    options = {
        "name": "HubbardHolstein",
        "nup": nup,
        "ndown": ndown,
        "nx": nsites,
        "ny": 1,
        "U": 0.,
        "g": g,
        "t": t,
        "w0": w0,
#       "lambda": lmbda,
        "maxiter": 1,
        "xpbc": pbc,
        "ypbc": False,
    }

    system = HubbardHolstein(options, verbose=verbose)
    e, _ = simple_fci_bose_fermi(
        system, system, nboson_max=nboson, hamil=False, verbose=verbose
    )

    return e[0]

energy = ipie_fci(2, 0.1, 1, 4, 1, verbose=True)
print(energy)