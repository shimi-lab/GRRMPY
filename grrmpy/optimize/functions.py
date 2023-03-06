from ase.optimize import FIRE

#USER
from grrmpy.functions import get_fmax

def add_stop_to(optimizer,max_n=12,times=2):
    """stop機能を付けたOptimizerに変換する
    
    max_n回連続,同じEnergy,Forceの繰り返しがtimes回あった場合に計算を終了する.

    Parameters:
    
    optimizer:class
        Optimizer.FIRE等
    max_n: int
        デフォルトは12
    times: int
        デフォルトは2

    Returns:
        class: Optimizer
    """ 
    custun_optimizer = type(f"{optimizer.__name__}withSTOP", (optimizer,),{})
    custun_optimizer.stop = False
    custun_optimizer.count1 = 0
    custun_optimizer.count2 = 0
    custun_optimizer.max_n = max_n
    custun_optimizer.times = times
    custun_optimizer.first=True
    custun_optimizer.stop_func = compare
    custun_optimizer.run = run
    return custun_optimizer

def run(self, fmax=0.05, steps=None):
    """self.stop=Trueになると計算を停止する."""  
    def irun(self):
        """Run dynamics algorithm as generator. This allows, e.g.,
        to easily run two optimizers or MD thermostats at the same time.

        Examples:
        >>> opt1 = BFGS(atoms)
        >>> opt2 = BFGS(StrainFilter(atoms)).irun()
        >>> for _ in opt2:
        >>>     opt1.run()
        """

        # compute initial structure and log the first step
        self.atoms.get_forces()

        # yield the first time to inspect before logging
        yield False

        if self.nsteps == 0:
            self.log()
            self.call_observers()

        # run the algorithm until converged or max_steps reached
        while not self.converged() and self.nsteps < self.max_steps and not self.stop:

            # compute the next step
            self.step()
            self.nsteps += 1

            # let the user inspect the step and change things before logging
            # and predicting the next step
            yield False

            # log the step
            self.log()
            self.call_observers()
            self.stop_func()

        # finally check if algorithm was converged
        yield self.converged()
        
    self.fmax = fmax
    if steps:
        self.max_steps = steps
        
    for converged in irun(self):
        pass
    return converged

def stop_iteration(opt,max_n=15,times=2):
    hef = HoldEnergyAndForce(opt,max_n,times)
    while True:
        yield hef.compare()
        

def compare(self):
    if self.first:
        self.pre_e = round(self.atoms.get_potential_energy(),6)
        self.pre_f = round(get_fmax(self.atoms),4)
        self.first = False
    e = round(self.atoms.get_potential_energy(),6)
    f = round(get_fmax(self.atoms),4)
    if e == self.pre_e and f==self.pre_f:
        self.count1+=1
        if self.count1 >= self.max_n:
            self.count1 = 0
            self.count2 += 1            
            if self.count2 >= self.times:
                self.stop=True
    else:
        self.pre_e = e
        self.pre_f = f
        self.count1 = 0