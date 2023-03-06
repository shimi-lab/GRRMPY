from grrmpy.excel.workbook import Workbook
from grrmpy.excel.worksheet import Worksheet
from grrmpy.excel.chart_path import ChartPath

__all__ = ["Workbook","Worksheet","ChartPath"]

    
"""
基本的な操作はxlsxwriterのWorkbook,Worksheetと同様である(継承しているので)

Workbookの変更点
add_chart(options)のoptionsのtypeに'path'を設定できるように変更した.
pathを設定するとgrrmpy.excell.ChartPathが作成される.

Worksheetの変更点
add_tableの第一引数がDataFrameだった場合に'A1'の位置にTableが作成されるように変更した

ChartPathについて
基本的にxlsxwriterのChartScatterと同じである.
add_pathというメソッドを追加している
"""