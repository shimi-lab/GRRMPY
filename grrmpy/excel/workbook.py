import xlsxwriter
from xlsxwriter.chart_area import ChartArea
from xlsxwriter.chart_column import ChartColumn
from xlsxwriter.chart_bar import ChartBar
from xlsxwriter.chart_doughnut import ChartDoughnut
from xlsxwriter.chart_line import ChartLine
from xlsxwriter.chart_pie import ChartPie
from xlsxwriter.chart_radar import ChartRadar
from xlsxwriter.chart_scatter import ChartScatter
from xlsxwriter.chart_stock import ChartStock
from warnings import warn

# Package imports.
from grrmpy.excel.worksheet import Worksheet
from grrmpy.excel.chart_path import ChartPath
  
class Workbook(xlsxwriter.Workbook):
    def __init__(self, filename=None, options=None):
        """_summary_

        Parameters:
        
        filename: str
            作成する新しい Excelファイルの名前
        options dict:
            `オプションのブック パラメーター <https://xlsxwriter.readthedocs.io/workbook.html>`_
        """
        super().__init__(filename, options)
        
    def add_worksheet(self, name=None):
        """ワークブックに新しいワークシートを追加する.

        Parameters:
        
        name : str
            ワークシート名. デフォルトではSheet1 など.

        Returns:
            Worksheet: grrmpy.excell.worksheet.Worksheetオブジェクト
        """
        return super().add_worksheet(name, Worksheet)
    
    def add_chart(self, options):
        """ワークシートに追加可能なグラフ オブジェクトを作成します.
        
        | 基本的にはxlsxwriterのWorkBook.add_chartと同様.
        | optionsのtypeに'path'を設定できるようにした.
        | 'path'を設定した時のみ,xlsxwriter.Chartオブジェクトでなく,grrmpy.excell.ChartPathオブジェクトが作成される.

        Parameters:
        
        options: dict
            チャート タイプ オプションのディクショナリ。

        Returns:
            Chart: xlsxwriterのChartオブジェクト
            
        Example:
        
        >>> chart = book.add_chart({'type': 'path'})
        
        """
        chart_type = options.get('type')
        if chart_type is None:
            warn("Chart type must be defined in add_chart()")
            return

        if chart_type == 'area':
            chart = ChartArea(options)
        elif chart_type == 'bar':
            chart = ChartBar(options)
        elif chart_type == 'column':
            chart = ChartColumn(options)
        elif chart_type == 'doughnut':
            chart = ChartDoughnut(options)
        elif chart_type == 'line':
            chart = ChartLine(options)
        elif chart_type == 'pie':
            chart = ChartPie(options)
        elif chart_type == 'radar':
            chart = ChartRadar(options)
        elif chart_type == 'scatter':
            chart = ChartScatter(options)
        elif chart_type == 'stock':
            chart = ChartStock(options)
        ###########新たに追加##############
        elif chart_type == "path":
            options.update({'type':'scatter','subtype': 'straight'})
            chart = ChartPath(options)
        #################################
        else:
            warn("Unknown chart type '%s' in add_chart()" % chart_type)
            return

        # Set the embedded chart name if present.
        if 'name' in options:
            chart.chart_name = options['name']

        chart.embedded = True
        chart.date_1904 = self.date_1904
        chart.remove_timezone = self.remove_timezone

        self.charts.append(chart)

        return chart