"""
argperseのヘルプのフォーマットをカスタマイズするクラス

使用方法
parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)

ArgumentParserのformatter_classに設定すれば良い.
"""

import argparse
from argparse import (OPTIONAL, SUPPRESS, ZERO_OR_MORE,
                      ArgumentDefaultsHelpFormatter, ArgumentParser,
                      RawDescriptionHelpFormatter, RawTextHelpFormatter)

class CustomHelpFormatter(RawTextHelpFormatter, RawDescriptionHelpFormatter, ArgumentDefaultsHelpFormatter):
    # 'max_help_position'のデフォルト値を「24」から「30」へ変更
    def __init__(self, prog, indent_increment=2, max_help_position=30, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)
    def _format_action(self, action: argparse.Action) -> str:
        return super()._format_action(action) + "\n"
    def _get_help_string(self, action):
        help = action.help
        if "%(default)" not in action.help:
            if action.default is not SUPPRESS:
                defaulting_nargs = [OPTIONAL, ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    # カスタム設定
                    if action.default is not None and not action.const:
                        help += " (default: %(default)s)"
        return help