from heapq import heappush,heappop

# heappush(h, (7, 'release product'))
# heappop(h)

class ShortestPath():
    def __init__(self,net_graph,priority=None):
        self.net = net_graph
        self.weight_edge = self.net._weight_edge
        self.group = self.net._group
        self.n_node = len(self.net.node)
        self.n_edge = len(self.net.node)
        if priority is None:
            self.priority = Priority()
        else:
            self.priority = priority
            
    def _search_base(self,source,target,group=True):
        # 無限大コストの要素を作成
        cost_list = [[float('inf') for _ in self.priority]
                     for _ in self.n_node]
        # heapqの作成
        heap_que = []
        
           
    def search(self,source,target,group=True):
        self.costlist = []
        
    def _search_with_group(self,source,target):
        pass
    
class Priority():
    def __init__(self):
        pass
    
    def push(self,weight_edge,group):
        pass
    
    def pop(self):
        pass