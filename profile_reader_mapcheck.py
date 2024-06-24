# profile reader for mapcheck
import sys
import re
from io import StringIO, BytesIO
from scipy.optimize import curve_fit
import pydicom

class MapcheckProfile(object):

    def __init__(self, mapcheck_file):
        content = {}
        lines = StringIO(mapcheck_file.getvalue().decode('utf-8')).readlines()

        found_background = False
        found_calib_factors = False
        found_offset = False
        found_raw_counts = False
        found_corr_counts = False
        found_dose_counts = False
        found_data_flags = False
        found_interpolated = False
        found_dose_interpolated = False
        unkeyed = []
        for i, line in enumerate(lines):
            if ":" in line:
                items = line.split("\t")
                num_params = int(len(items) / 2)
                for j in range(num_params):
                    key = items[2 * j].replace(":", "")
                    value = items[2 * j + 1].replace("\n", "")
                    content[key] = value

            elif 'Background' in line:
                found_background = True
                s_background = i
            elif found_background and 'Calibration Factors' in line:
                e_background = i - 1
                found_calib_factors = True
                s_calib_factors = i
            elif found_calib_factors and 'Offset' in line:
                e_calib_factors = i - 1
                found_offset = True
                s_offset = i
            elif found_offset and 'Raw Counts' in line:
                e_offset = i - 1
                found_raw_counts = True
                s_raw_counts = i
            elif found_raw_counts and 'Corrected Counts' in line:
                e_raw_counts = i - 1
                found_corr_counts = True
                s_corr_counts = i
            elif found_corr_counts and 'Dose Counts' in line:
                e_corr_counts = i - 1
                found_dose_counts = True
                s_dose_counts = i
            elif found_dose_counts and 'Data Flags' in line:
                e_dose_counts = i - 1
                found_data_flags = True
                s_data_flags = i
            elif found_data_flags and 'Interpolated' in line and 'Dose' not in line:
                e_data_flags = i - 1
                found_interpolated = True
                s_interpolated = i
            elif found_interpolated and 'Dose Interpolated' in line:
                e_interpolated = i - 1
                found_dose_interpolated = True
                s_dose_interpolated = i
            elif found_dose_interpolated and 'Xcm' in line:
                line_col_params = i

            elif not found_background:
                if line.strip() != '':
                    unkeyed.append(line.strip())

            e_dose_interpolated = len(lines)

        self.background = read_data_section(lines, s_background, e_background)
        self.calibration_factors = read_data_section(lines, s_calib_factors, e_calib_factors)
        self.offset = read_data_section(lines, s_offset, e_offset)
        self.raw_counts = read_data_section(lines, s_raw_counts, e_raw_counts)
        self.corrected_counts = read_data_section(lines, s_corr_counts, e_corr_counts)
        self.dose_counts = read_data_section(lines, s_dose_counts, e_dose_counts)
        self.data_flags = read_data_section(lines, s_data_flags, e_data_flags)
        self.interpolated = read_data_section(lines, s_interpolated, e_interpolated)
        self.dose_interpolated = read_data_section(lines, s_dose_interpolated, e_dose_interpolated)
        self.col_params = read_col_params(lines[line_col_params])
        self.row_params = read_row_params(self.dose_interpolated)
        self.date = content.get('Date')
        self.time = content.get('Time')

    def select_slice_by_index(self, x_idx, y_idx):

        data_array = []
        if x_idx is None and y_idx is not None:
            data_array = self.dose_interpolated[y_idx][2:]
        elif x_idx is not None and y_idx is None:
            for y_data_array in self.dose_interpolated:
                data_array.append(y_data_array[x_idx])

        return data_array

    def select_slice(self, x_slice, y_slice):

        data_array = self.dose_interpolated
        x_array = self.col_params

        if x_slice is not None:
            x_index = None
            for i in range(len(x_array) - 1):
                if x_array[i] <= x_slice < x_array[i + 1]:
                    x_index = i
                    break
            if x_index is None and x_slice <= x_array[-1]:
                x_index = len(x_array) - 1

            if x_index is None:
                return False, None, None, None

            s_array = []
            d_array = []
            for row in data_array:
                s_array.append(row[0] * 10.0)
                d_array.append(row[x_index + 2])

            return True, x_array[x_index], s_array, d_array

        else:
            y_index = None
            d_array = []
            s_array = [x * 10.0 for x in x_array]
            for i in range(len(data_array) - 1):
                if data_array[i + 1][0] < y_slice <= data_array[i][0]:
                    if y_slice - data_array[i + 1][0] < data_array[i][0] - y_slice:
                        y_index = i + 1
                    else:
                        y_index = i
                    d_array = data_array[y_index][2:]
                    break
            if y_index is None and data_array[-1][0] >= y_slice:
                y_index = len(data_array) - 1
                d_array = data_array[y_index][2:]

            if y_index is None:
                return False, None, None, None

            return True, data_array[y_index][0] * 10., s_array, d_array

    def calc_center_avg(self, d_range):

        cols = self.col_params
        dose_dist = self.dose_interpolated

        count = 0
        total = 0
        center_dose = []
        dose_max = 0
        for slice in dose_dist:
            if abs(slice[0]) <= d_range:
                for col, dose in zip(cols, slice[2:]):
                    if abs(col) <= d_range:
                        count += 1
                        center_dose.append([slice[0], col, dose])
                        total += dose
                        if dose > dose_max:
                            dose_max = dose

        # return total / count
        return dose_max

def read_data_section(lines, s_index, e_index):

    data_matrix = []
    for line in lines[s_index+2:e_index-2]:
        if not re.search('[a-zA-Z]', line):
            data_matrix.append([float(x) for x in line.split('\t')])
    return data_matrix


def read_col_params(line):

    col_params = [float(col) for col in [x for x in line.split('\t')][2:]]
    return col_params


def read_row_params(dose_data):

    row_params = [y[0] for y in dose_data]
    return row_params


def main(argv):

    input_file = argv[0]
    MapcheckProfile(input_file)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
