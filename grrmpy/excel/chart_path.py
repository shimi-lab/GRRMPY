from xlsxwriter.chart_scatter import ChartScatter

from grrmpy.path.reaction_path import ReactPath

class ChartPath(ChartScatter):
    def __init__(self, options=None):
        super().__init__(options)
        
    def add_path(self, path:ReactPath):
        pass