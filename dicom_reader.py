
import sys
import math
import itertools
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import pydicom
from pydicom.pixel_data_handlers.util import apply_modality_lut
from pydicom.pixel_data_handlers.util import apply_voi_lut
from io import BytesIO

class DicomProfile(object):

    def __init__(self, dicom_file):
        self.dcm_cols = None
        self.dcm_rows = None
        self.data_array = None
        self.dataset = pydicom.dcmread(BytesIO(dicom_file.getvalue()))
        self.dicom_reader()

    def dicom_reader(self):

        conv_to_cgray = 100.
        image_orientation_patient = self.dataset.ImageOrientationPatient
        image_orientation_rows = image_orientation_patient[:3]
        image_orientation_cols = image_orientation_patient[3:]
        image_position = self.dataset.ImagePositionPatient
        pixel_spacing = self.dataset.PixelSpacing
        dose_grid_scaling = self.dataset.DoseGridScaling
        num_cols = self.dataset.Columns
        num_rows = self.dataset.Rows
        pixel_array = self.dataset.pixel_array
        data_array = apply_voi_lut(pixel_array, self.dataset, index=0)

        image_position_x = []
        for i in range(num_cols):
            image_position_x.append(image_orientation_rows[0] * pixel_spacing[0] * float(i) + image_position[0])

        image_position_y = []
        for i in range(num_rows):
            image_position_y.append(image_orientation_cols[2] * pixel_spacing[1] * float(i) + image_position[2])

        image_position_y = [y - (image_position_y[0]+image_position_y[-1])/2. for y in image_position_y]

        dose_array = []
        for i in range(num_cols):
            temp = []
            for j in range(num_rows):
                temp.append(data_array[i, j] * dose_grid_scaling * conv_to_cgray)
            dose_array.append(temp)

        self.dcm_cols = [round(x/10., 2) for x in image_position_x]
        self.dcm_rows = [round(y/10., 2) for y in image_position_y]
        self.data_array = np.array(dose_array)

        return

    def dicom_reader_plt(self):

        # Normal mode:

        print("Storage type.....:", self.dataset.SOPClassUID)
        print()

        pat_name = self.dataset.PatientName
        image_position = self.dataset.ImagePositionPatient
        pixel_spacing = self.dataset.PixelSpacing
        dose_grid_scaling = self.dataset.DoseGridScaling
        display_name = pat_name.family_name + ", " + pat_name.given_name
        print("Patient's name...:", display_name)
        print("Patient id.......:", self.dataset.PatientID)
        print("Modality.........:", self.dataset.Modality)
        print("Study Date.......:", self.dataset.StudyDate)

        if 'PixelData' in self.dataset:
            rows = int(self.dataset.Rows)
            cols = int(self.dataset.Columns)
            print("Image size.......: {rows:d} x {cols:d}, {size:d} bytes".format(
                rows=rows, cols=cols, size=len(self.dataset.PixelData)))
            if 'PixelSpacing' in self.dataset:
                print("Pixel spacing....:", self.dataset.PixelSpacing)

        # use .get() if not sure the item exists, and want a default value if missing
        print("Slice location...:", self.dataset.get('SliceLocation', "(missing)"))

        # plot the image using matplotlib
        plt.imshow(self.dataset.pixel_array, cmap=plt.cm.hot)
        plt.show()

    def get_axes(self):

        image_position = self.dataset.ImagePositionPatient
        pixel_spacing = self.dataset.PixelSpacing
        dose_grid_scaling = self.dataset.DoseGridScaling
        num_cols = self.dataset.Columns
        num_rows = self.dataset.Rows

        x_array = []
        for i in range(num_cols):
            x_array.append(pixel_spacing[0] * float(i) + image_position[0])

        y_array = []
        for j in range(num_rows):
            y_array.append([image_position[2] - float(j) * pixel_spacing[1] for j in range(num_cols)])

        return x_array, y_array

    def get_dataset(self):

        conv_to_cgray = 100.
        dose_grid_scaling = self.dataset.DoseGridScaling
        num_cols = self.dataset.Columns
        num_rows = self.dataset.Rows
        pixel_array = self.dataset.pixel_array
        data_array = apply_voi_lut(pixel_array, self.dataset, index=0)

        for i in range(num_cols):
            for j in range(num_rows):
                data_array[i, j] = data_array[i, j] * dose_grid_scaling * conv_to_cgray

        return data_array

    def calc_center_avg(self, d_range):

        count = 0
        total = 0
        center_dose = []
        dose_max = 0
        for i, x in enumerate(self.dcm_cols):
            if abs(x) <= d_range:
                for j, y in enumerate(self.dcm_rows):
                    if abs(y) <= d_range:
                        count += 1
                        dose = self.data_array[i, j]
                        center_dose.append([x, y, dose])
                        total += dose
                        if dose > dose_max:
                            dose_max = dose


        # return total/count
        return dose_max

    def select_slice_by_index(self, x_idx, y_idx):

        data_array = []
        if x_idx is None and y_idx is not None:
            data_array = self.data_array[y_idx]
        elif x_idx is not None and y_idx is None:
            for y_data_array in self.data_array:
                data_array.append(y_data_array[x_idx])

        return data_array

    def select_slice(self, x_slice, y_slice):

        if x_slice is not None:
            x_slice *= 10.
        if y_slice is not None:
            y_slice *= 10.
        conv_to_cgray = 100.
        image_position = self.dataset.ImagePositionPatient
        pixel_spacing = self.dataset.PixelSpacing
        dose_grid_scaling = self.dataset.DoseGridScaling
        num_cols = self.dataset.Columns
        num_rows = self.dataset.Rows
        pixel_array = self.dataset.pixel_array
        data_array = apply_voi_lut(pixel_array, self.dataset, index=0)

        if x_slice is not None:
            x_index = 0
            for i in range(num_rows):
                x0 = pixel_spacing[0] * float(i) + image_position[0]
                x1 = pixel_spacing[0] * float(i+1) + image_position[0]
                if x0 <= x_slice < x1:
                    if x_slice-x0 <= x1-x_slice:
                        x_index = i
                    else:
                        x_index = i + 1
            x_pos = pixel_spacing[0] * x_index + image_position[0]
            d_array = [x[x_index]*dose_grid_scaling*conv_to_cgray for x in data_array]
            y_array = [round((image_position[2] - float(i)*pixel_spacing[1])/10., 2) for i in range(num_cols)]
            return x_pos, y_array, d_array
        else:
            y_index = 0
            for i in range(num_cols):
                y0 = -pixel_spacing[1]*float(i) + image_position[2]
                y1 = -pixel_spacing[1] * float(i+1) + image_position[2]
                if y1 <= y_slice < y0:
                    if y0-y_slice <= y_slice-y1:
                        y_index = i
                    else:
                        y_index = i + 1
            y_pos = -pixel_spacing[1] * y_index + image_position[2]
            d_array = [x*dose_grid_scaling*conv_to_cgray for x in data_array[y_index]]
            x_array = [round((image_position[0] + float(i)*pixel_spacing[0])/10., 2) for i in range(num_rows)]
            return y_pos, x_array, d_array


def dose_curve(x, eta, t, a, b1, b2, d1, x0):

    return eta * (t + (1 - t) * a * erf(b1 * (x0 - x)) + (1 - a) * erf(b2 * (x0 - x))) + d1


def newton_approx(x, x_data, y_data):

    # find x_data bin for x
    ix = 0
    for i in range(len(x_data)):
        if i < len(x_data)-1 and x_data[i] <= x < x_data[i+1]:
            ix = i
            break

    y_pred = y_data[ix]
    # calculate slope at both ends of the bin
    if x != x_data[ix] and 0 < ix < len(x_data)-1:
        s0 = (y_data[ix] - y_data[ix-1])/(x_data[ix] - x_data[ix-1])
        s1 = (y_data[ix+1] - y_data[ix])/(x_data[ix+1] - x_data[ix])
        y0 = y_data[ix] + (x - x_data[ix]) * s0
        y1 = y_data[ix+1] + (x - x_data[ix+1]) * s1
        w0 = 1./(x-x_data[ix])
        w1 = 1./(x_data[ix+1] - x)
        y_pred = (y0 * w0 + y1 * w1) / (w0 + w1)

    return y_pred

