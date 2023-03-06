def run(self, fmax=0.05, steps=None):
    """ call Dynamics.run and keep track of fmax"""
    
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

        # finally check if algorithm was converged
        yield self.converged()
        
    self.fmax = fmax
    if steps:
        self.max_steps = steps
    for converged in irun(self):
        pass
    return converged

def attach_stop(self, function, interval=1, *args, **kwargs):
    def function1(function=function):
        b = function()
        if b:
            self.stop = True
            
    self.attach(function1, interval=interval, *args, **kwargs) 
    

def OptWithCond(optimizer, atoms, **kwargs):
    """
    
    | OptWithCond.stopがTrueになると計算が停止するOptimizer
    | attach_stop()の引数にTrue,Falseを返す関数をアタッチできる.
    | イタレーション毎にattach_stopで設定した関数が実行され,関数の返り値がTrueになればOptWithCond.stopがTureとなり,計算が停止する.

    Parameters:
    
    optimizer: class
        | 継承するoptimizer
    atoms: Atoms
        | Atomsオブジェクト
    kwargs:
        | optimizerの引数

    Returns:
        class: Stop機能を付けたOptimizer object
    """
    custun_optimizer = type(f"{optimizer.__name__}withCOND", (optimizer,),{})
    custun_optimizer.run = run
    custun_optimizer.attach_stop = attach_stop
    custun_optimizer.stop = False
    return custun_optimizer(atoms, **kwargs)