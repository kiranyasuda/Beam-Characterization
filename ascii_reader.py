import sys
import re
from typing import List
from io import StringIO
from matplotlib.figure import Figure


def read_data_section(lines: List[str], s_index: int, e_index: int) -> List[List[float]]:
    data_matrix = []
    for line in lines[s_index:e_index]:
        if not re.search('[a-zA-Z]', line):
            items = line.replace('<', '').replace('>', '').replace('+', '').replace('=', '').split()
            data_matrix.append([float(x.strip()) for x in items])
    return data_matrix


class AsciiProfile(object):

    def __init__(self, ascii_file):
        self.content = {'data': [[]]}
        self.lines = StringIO(ascii_file.getvalue().decode('utf-8')).readlines()
        self.data_sects = []


        #with open(ascii_file) as i_file:
        #lines = i_file.readlines()

        for i, line in enumerate(self.lines):
            if "%VNR 1.0" in line:
                self.read_version1()
                break
            elif "%VERSION 02" in line:
                self.read_version2()
                break

    def read_version1(self):
        found_start = False

        for i, line in enumerate(self.lines):
            if "%" in line:
                items = line.split("\t")
                items[-1] = items[-1].strip()
                self.content[items[0].strip()] = items[1:]
            elif not found_start and "=" in line:
                self.content['start'] = i
                found_start = True
            elif ":EOM" in line:
                self.content['end'] = i
                found_start = False
                self.data_sects.append(self.content)
                self.content = {'data': [[]]}

        for data_sect in self.data_sects:
            field_size = 'x'.join(data_sect['%FSZ'])
            depth = str(data_sect['%STS'][2]).split('#')[0].strip()
            if field_size is not None:
                data_sect['fsz'] = field_size
            if depth is not None:
                data_sect['depth'] = depth

            data_sect['data'] = read_data_section(self.lines, data_sect['start'], data_sect['end'])

            if data_sect['%STS'][0] != data_sect['%EDS'][0]:
                data_sect['axis'] = "X"
            elif data_sect['%STS'][1] != data_sect['%EDS'][1]:
                data_sect['axis'] = 'Y'
            else:
                data_sect['axis'] = 'Z'

    def read_version2(self):
        found_start = False

        for i, line in enumerate(self.lines):
            if "%" in line:
                items = line.split()
                items[-1] = items[-1].strip()
                self.content[items[0].strip()] = items[1:]
            elif not found_start and "<" in line:
                self.content['start'] = i
                found_start = True
            elif "$ENOM" in line:
                self.content['end'] = i
                found_start = False
                self.data_sects.append(self.content)
                self.content = {'data': [[]]}
        for data_sect in self.data_sects:
            field_size = 'x'.join([str(y) for y in [int(x) for x in data_sect['%FLSZ'][0].split('*')]])
            field_axis = data_sect['%AXIS'][0]
            field_depth = int(data_sect['%DPTH'][0])
            if field_size is not None:
                data_sect['fsz'] = field_size
            if field_axis == 'X':
                data_sect['axis'] = 'X'
            else:
                data_sect['axis'] = 'Y'
            if field_depth is not None:
                data_sect['depth'] = field_depth

            data_sect['data'] = read_data_section(self.lines, data_sect['start'], data_sect['end'])


    def graph(self, field_size, field_depth, field_axis):
        for data_sect in self.data_sects:
            if (data_sect['fsz'] == field_size and int(float(data_sect['depth'])) == field_depth and
                    data_sect['axis'] == field_axis):
                if field_axis == "X":
                    axis_array = [x[0] for x in data_sect['data']]
                elif field_axis == "Y":
                    axis_array = [x[1] for x in data_sect['data']]
                else:
                    axis_array = [x[2] for x in data_sect['data']]
                data_array = [x[3] for x in data_sect['data']]
                fig = Figure(figsize=(5, 5), dpi=100)
                ax1 = fig.add_subplot(111)
                if data_sect['axis'] == 'X':
                    ax1.set_xlabel("X (mm)")
                else:
                    ax1.set_xlabel("Y (mm)")
                ax1.set_ylabel("cGray")
                ax1.plot(axis_array, data_array)
                return fig




def main(argv):
    input_file = 'C:/Users/kiran/Profiles/6MV PROFILES.asc'
    AsciiProfile(input_file)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
