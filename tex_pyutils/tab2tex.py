# Copyright (C) 2013-2017 University of Southern California and Liz Ji
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import os
import argparse

################################################################################
#-------------------------------------------------------------------------------
TABLE_NAME="{table}"
TABULAR_NAME="{tabular}"

TABLE_FRAME="""\\begin{type_name}
{content}
\\end{type_name}
"""
TABULAR_FRAME="""\\begin{type_name} {col_format}
{content}
\\end{type_name}
"""

def load_table(fin, col):
    IN = open(fin)
  
    tab = []
    for line in IN:
        items = line.rstrip().split('\t')
        tab.append([items[i] for i in col])
    return tab


def set_tabular_col_format(num_col, align_type='c'):
    col_formats = [align_type] * num_col
    return "{|%s|}" % ('|'.join(col_formats))


def write_rows(tab):
    rows = ["\\hline"]
    for line in tab:
      rows.append(' & '.join(line) + ' \\\\')
      rows.append("\\hline")
    return rows

def fit_in_page(content):
    text = "\\resizebox{\\textwidth}{!}{%s}" % (content)
    return text

def replace_illegal_chars(text):
  text = text.replace('_', '\\_')
  return text


#==============================================================================
def main():
    parser = argparse.ArgumentParser(description = 'BED filter')
    parser.add_argument('-c', '--column', dest='column', nargs='+', \
                        default=[1, 2, 3], \
                        help='the modes you want to run with BED filter')
    parser.add_argument('-i', '--input', dest='fin', required=True, \
                        help="input table")
    parser.add_argument('-o', '--output', dest='fout', required=True, \
                        help="output txt")
   
    args = parser.parse_args()
    tab = load_table(args.fin, args.column)
    
    tabular_col_format = set_tabular_col_format(len(args.column))

    rows = write_rows(tab)

    tabular_content = TABULAR_FRAME.format(type_name=TABULAR_NAME, \
                                           col_format=tabular_col_format, \
                                           content = '\n'.join(rows))
    table_content = TABLE_FRAME.format(type_name=TABLE_NAME, \
                                       content = fit_in_page(tabular_content))

    FOUT = open(args.fout, 'w')
    FOUT.write(replace_illegal_chars(table_content))
    FOUT.close()

#==============================================================================
if __name__ == '__main__':
  main()
