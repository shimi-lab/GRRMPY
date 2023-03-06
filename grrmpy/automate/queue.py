from heapq import heapify,heappop,heappush
from typing import List
from pathlib import Path
import pickle
from ase.optimize import LBFGS,FIRE
from ase.neb import NEB,interpolate
from ase import Atoms
# USER
from grrmpy.calculator import pfp_calculator
from grrmpy.functions import set_calc2images,set_const2images,copy_images
from grrmpy.path.reaction_path import ReactPath
from grrmpy.neb.repetitive_neb import RepetitiveNEB
from grrmpy.neb.functions import get_energy_list
try:
    from grrmpy.optimize import FIRELBFGS
    OPT_OPTIMIZER = FIRELBFGS
except:
    OPT_OPTIMIZER = FIRE
    
class Queue():
    """    
    executorにクラスを指定して使用する.
    executorのクラスはrun()メソッドを必ず有する必要がある.
    さらにSTOP.txtで処理を停止し,途中再開できるようにするにはtodict()とfromdictのクラスメソッドを有する必要がある.

    Parameters:
    
    executor:
        | 計算実行クラスオブジェクト
    priority: int
        | queueはタプルのリストであるが,階層の情報をリストの何番目に追加するかを指定する.
        | priority=0の場合,階層の浅さが全てに優先されるため,幅優先になる
        | priority>0の時,深さ優先になる.
    maxjob: int
        | maxjob回の計算を実行したら自動で計算を停止する.
    """
    reastart_file = "ReStart.pkl"
    def __init__(self,executor,priority=0,maxjob=float('inf')):
        """    
        executorにクラスを指定して使用する.
        executorのクラスはrun()メソッドを必ず有する必要がある.
        さらにSTOP.txtで処理を停止し,途中再開できるようにするにはtodict()とfromdictのクラスメソッドを有する必要がある.
    
        Parameters:
        
        executor:
            計算実行クラスオブジェクト
        priority: int
            queueはタプルのリストであるが,階層の情報をリストの何番目に追加するかを指定する.
            priority=0の場合,階層の浅さが全てに優先されるため,幅優先になる
            priority>0の時,深さ優先になる.
        maxjob: int
            maxjob回の計算を実行したら自動で計算を停止する.
        """
        self.executor = executor
        self.executor_cls = executor.__class__
        self.queue = []
        self.priority = priority
        self.maxjob = maxjob
        heapify(self.queue)
        self.attach_list = []
        
    def attach(self,functions):
        """1回の計算が終える度に実行する関数をアタッチする"""
        self.attach_list.append(functions)
    
    def topkl(self):
        data_dict_list = []
        for data in self.queue:
            data = list(data)
            layer = data.pop(self.priority)
            data_dict = self.executor.data_to_serializable_obj(data)
            data_dict.insert(self.priority,layer)
            data_dict_list.append(data_dict)
        with open(self.reastart_file,"wb") as f:
            pickle.dump(data_dict_list, f)
        Path("STOP.txt").unlink()
            
    def frompkl(self,file):
        queue = []
        with open(self.reastart_file,"rb") as f:
            data_dict_list = pickle.load(f)
        for data_dict in data_dict_list:
            layer = data_dict.pop(self.priority)
            data_list = list(self.executor_cls.data_from_serializable_obj(data_dict))
            data_list.insert(self.priority, layer)
            queue.append(data_list)
        return queue
    
    def push(self,data,layer):
        for d in data:
            d = list(d)
            d.insert(self.priority, layer)
            heappush(self.queue,d)
        
    def pop(self):
        data = list(heappop(self.queue))
        layer = data.pop(self.priority)
        return data, layer
        
    def run(self,restart=False):
        if restart:
            self.queue = self.frompkl(self.reastart_file)
        else:
            data = self.executor.first_push()
            self.push(data,0)
        
        self.done_job=0
        while len(self.queue) > 0 and self.done_job < self.maxjob:
            data,layer = self.pop()
            pushdata = self.executor.run_and_push(data,layer)
            self.push(pushdata, layer+1)
            for function in self.attach_list:
                function()
            self.done_job += 1
            
            if Path("STOP.txt").exists():
                #保存処理処理
                self.topkl()
                break
            
            
class Executor():
    """Queueクラスのインスタンス引数のexecutorのひな形クラス
    
    完全に自作しても良いが，必ず必要なメソッドをここに示す
    """
    def __init__(self,images,constraints=[],calc_func=pfp_calculator):
        self.images = images
        self.constraints = constraints
        self.calc_func = calc_func
    
    def first_push(self)->list:
        """一番始めにQueueに追加するデータ"""
        return self.images
    
    def run_and_push(self,data,layer)->list:
        """Queueからデータをもらい,新たに計算すべきデータを与える"""
        return []
    
    def data_to_serializable_obj(self,data):
        """dataをpickleできる状態に加工して返す(calculatorを消す等)"""
        for image in data:
            del image.calc
        return data
    
    @classmethod
    def data_from_serializable_obj(cls,data):
        """topyobjで加工した場合,元に戻す(calculatorを再度付け直すなど)"""
        for image in data:
            image.calc = cls.calc_func()
        return data
    
    def to_obj(self):
        return self
    
    @classmethod
    def from_obj(cls,obj):
        return obj