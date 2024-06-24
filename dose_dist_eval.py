# -*- coding: utf-8 -*-
# Author: Taka
# Created on 5/31/2020 10:04 AM
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy import interpolate
import profile_reader_mapcheck as prm
import dicom_reader as dcr

D_DTA = 0.3
D_DD = 0.03


def find_min_gamma(dcm_cols, dcm_rows, dcm_data_ar, x_pt, y_pt, dose, popt):
    dcm_max = dcm_data_ar.max()
    x_index = None
    for i in range(len(dcm_cols) - 1):
        if dcm_cols[-1] > dcm_cols[0]:
            if dcm_cols[i] <= x_pt < dcm_cols[i + 1]:
                x_index = i
                break
        else:
            if dcm_cols[i + 1] <= x_pt < dcm_cols[i]:
                x_index = i
                break

    y_index = None
    for i in range(len(dcm_rows) - 1):
        if dcm_rows[-1] > dcm_rows[0]:
            if dcm_rows[i] <= y_pt < dcm_rows[i + 1]:
                y_index = i
                break
        else:
            if dcm_rows[i + 1] <= y_pt < dcm_rows[i]:
                y_index = i
                break

    if x_index is None or y_index is None:
        return None, 0., 1.

    # find subset for fitting
    pos_x = []
    pos_y = []
    sub_data = []
    for iy in range(y_index - 2, y_index + 3):
        if iy < len(dcm_data_ar):
            arr = dcm_data_ar[iy]
            if len(sub_data) > 0:
                sub_data = np.concatenate((sub_data, arr[x_index - 2:x_index + 3]), axis=0)
            else:
                sub_data = arr[x_index - 2:x_index + 3]
            for ix in range(x_index - 2, x_index + 3):
                pos_x.append(dcm_cols[ix])
                pos_y.append(dcm_rows[iy])

    try:
        # interpolate.bisplrep: Find a bivariate B-spline representation of a surface.
        # Given a set of data points (x[i], y[i], z[i]) representing a surface z=f(x,y),
        # compute a B-spline representation of the surface.
        tck = interpolate.bisplrep(np.array(pos_x), np.array(pos_y), np.array(sub_data), s=0)

        # find the minimum gamma
        x_min = 0.
        y_min = 0.
        g_min = 1000.
        for x in pos_x:
            for y in pos_y:
                # interpolate.bisplev: Return a rank-2 array of spline function values (or spline derivative values)
                # at points given by the cross-product of the rank-1 arrays x and y.
                z = interpolate.bisplev(x, y, tck)
                g = gamma_func(math.sqrt((x - x_pt) ** 2 + (y - y_pt) ** 2), abs(z - dose) / dose, D_DTA,
                               dcm_max * D_DD)
                if g < g_min:
                    g_min = g
                    x_min = x
                    y_min = y

        x0 = np.array([x_min, y_min])
        # scipy.optimize.minimize: Minimization of scalar function of one or more variables.
        result = minimize(gamma, x0, args=(x_pt, y_pt, tck, dose, D_DTA, dcm_max * D_DD))
        z_fit = interpolate.bisplev(result.x[0], result.x[1], tck)
        dta = math.sqrt((x_pt - result.x[0]) ** 2 + (y_pt - result.x[1]) ** 2)
        dd = (z_fit - dose)

        return z_fit, dta, dd

    except RuntimeError as e:
        return None, 0., 1.


def gamma_func(dta, dd, d_dta, d_dd):
    return math.sqrt((dta / d_dta) ** 2 + (dd / d_dd) ** 2)


def gamma(x, x_pt, y_pt, tck, dose, d_dta, d_dd):
    z = interpolate.bisplev(x[0], x[1], tck)
    dta = math.sqrt((x_pt - x[0]) ** 2 + (y_pt - x[1]) ** 2)
    dd = z - dose

    return gamma_func(dta, dd, d_dta, d_dd)


def heatmap(data, row_labels, col_labels, tick_freq=1, ax=None,
            cbar_kw=None, cbarlabel="", ticks=True, **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    tick_freq
        A frequency of tick
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    ticks
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if cbar_kw is None:
        cbar_kw = {}

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    if ticks:
        # We want to show all ticks...
        ax.set_xticks(np.arange(0, data.shape[1], tick_freq))
        ax.set_yticks(np.arange(0, data.shape[0], tick_freq))
        # ... and label them with the respective list entries.
        ax.set_xticklabels(col_labels[0::tick_freq])
        ax.set_yticklabels(row_labels[0::tick_freq])

        # Let the horizontal axes labeling appear on top.
        ax.tick_params(top=True, bottom=False,
                       labeltop=True, labelbottom=False)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
                 rotation_mode="anchor")

        # Turn spines off and create white grid.
        for edge, spine in ax.spines.items():
            spine.set_visible(False)

        # ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=False)
        # ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=False)
        # ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
        # ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


class DoseProfiler(object):

    def __init__(self, mpc_file, dicom_file):
        self.mpc = prm.MapcheckProfile(mpc_file)
        self.dcm = dcr.DicomProfile(dicom_file)
        self.norm = self.normalize()
        self.mpc_data = None
        self.dd = None
        self.dta = None
        self.gamma = None

    def normalize(self):
        d_range = 1.0
        mpc_avg = self.mpc.calc_center_avg(d_range)
        dic_avg = self.dcm.calc_center_avg(d_range)

        return dic_avg / mpc_avg

    def calculate(self, scrolled_text):

        mpc_rows = [row[0] for row in self.mpc.dose_interpolated]

        mpc_data_ar = []
        dta_ar = []
        dd_ar = []
        gamma_ar = []
        x_slice = None
        x_idx = None
        for y_idx, y_slice in enumerate(mpc_rows):
            print("Calculating gamma at y = " + str(y_slice) + " cm")
            ##scrolled_text.insert(tkinter.END, "Calculating gamma at y = " + str(y_slice) + " cm\n")
            ##scrolled_text.yview(tkinter.END)
            ##scrolled_text.update()

            data_mpc = self.mpc.select_slice_by_index(x_idx, y_idx)

            data_mpc = [x * self.norm for x in data_mpc]

            # minimize gamma
            data_dta = []
            data_dd = []
            data_gamma = []
            for ax_pos, data_pt in zip(self.mpc.col_params, data_mpc):
                if abs(ax_pos) > 10.:
                    data_dd.append(0.)
                    data_dta.append(0.)
                    data_gamma.append(0.)
                    continue
                popt = None
                z_fit = 1.
                if data_pt == 0:
                    dta = 0.
                    dd = 0.
                else:
                    z_fit, dta, dd = find_min_gamma(self.dcm.dcm_cols, self.dcm.dcm_rows,
                                                    self.dcm.data_array, ax_pos, y_slice, data_pt, popt)
                if dta is not None and dd is not None:
                    if z_fit is None:
                        z_fit = 1.
                    g_value = gamma_func(dta, dd / z_fit, D_DTA, self.dcm.data_array.max() * D_DD)
                    data_dta.append(dta)
                    data_dd.append(abs(dd))
                    data_gamma.append(g_value)

            mpc_data_ar.append(data_mpc)
            dta_ar.append(data_dta)
            dd_ar.append(data_dd)
            gamma_ar.append(data_gamma)

        self.mpc_data = mpc_data_ar
        self.dd = np.array(dd_ar)
        self.dta = np.array(dta_ar)
        self.gamma = np.array(gamma_ar)

    def select_dd_slice_by_index(self, x_idx, y_idx):

        data_array = []
        if x_idx is None and y_idx is not None:
            data_array = self.dd[y_idx]
        elif x_idx is not None and y_idx is None:
            for y_data_array in self.dd:
                data_array.append(y_data_array[x_idx])

        return data_array

    def select_dta_slice_by_index(self, x_idx, y_idx):

        data_array = []
        if x_idx is None and y_idx is not None:
            data_array = self.dta[y_idx]
        elif x_idx is not None and y_idx is None:
            for y_data_array in self.dta:
                data_array.append(y_data_array[x_idx])

        return data_array

    def select_gamma_slice_by_index(self, x_idx, y_idx):

        data_array = []
        if x_idx is None and y_idx is not None:
            data_array = self.gamma[y_idx]
        elif x_idx is not None and y_idx is None:
            for y_data_array in self.gamma:
                data_array.append(y_data_array[x_idx])

        return data_array

    def do_it(self):

        self.calculate()

        # dcm x slice at y = 0
        dcm_data_ar_x0 = self.dcm.data_array[int(len(self.dcm.data_array) / 2) + 1]
        # dcm y slice at x = 0
        x_ind = int(len(self.dcm.data_array[0]) / 2 + 1)
        dcm_data_ar_y0 = []
        for i in range(len(self.dcm.data_array)):
            dcm_data_ar_y0.append(self.dcm.data_array[i][x_ind])

        mpc_cols = self.mpc.col_params
        mpc_rows = [row[0] for row in self.mpc.dose_interpolated]

        dta_hist_pars = {'min': 0., 'max': 2., 'n_bins': 20}
        dd_hist_pars = {'min': 0., 'max': 1., 'n_bins': 20}
        gamma_hist_pars = {'min': 0., 'max': 2., 'n_bins': 20}

        dta_hist = []
        dd_hist = []
        gamma_hist = []

        fig1, ax1 = plt.subplots()
        ax1.plot(self.dcm.dcm_cols, dcm_data_ar_x0)
        fig2, ax2 = plt.subplots()
        ax2.plot(self.dcm.dcm_rows, dcm_data_ar_y0)
        fig3, ax3 = plt.subplots()
        heatmap(np.array(self.mpc_data), mpc_rows, mpc_cols, ax=ax3)
        fig4, ax4 = plt.subplots()
        heatmap(np.array(self.dcm.data_array), self.dcm.dcm_rows, self.dcm.dcm_cols, ticks=False, ax=ax4)
        fig5, ax5 = plt.subplots()
        heatmap(np.array(self.dta), mpc_rows, mpc_cols, ax=ax5)
        fig6, ax6 = plt.subplots()
        heatmap(np.array(self.dd), mpc_rows, mpc_cols, ax=ax6)
        fig7, ax7 = plt.subplots()
        heatmap(np.array(self.gamma), mpc_rows, mpc_cols, ax=ax7)
        fig8, ax8 = plt.subplots()
        ax8.hist(dta_hist, range=(dta_hist_pars['min'], dta_hist_pars['max']), bins=dta_hist_pars['n_bins'])
        fig9, ax9 = plt.subplots()
        ax9.hist(dd_hist, range=(dd_hist_pars['min'], dd_hist_pars['max']), bins=dd_hist_pars['n_bins'])
        fig10, ax10 = plt.subplots()
        ax10.hist(gamma_hist, range=(gamma_hist_pars['min'], gamma_hist_pars['max']), bins=gamma_hist_pars['n_bins'])
        x_ind = int(len(mpc_rows) / 2) + 1
        y_ind = int(len(mpc_cols) / 2) + 1
        fig11, ax11 = plt.subplots()
        ax11.set_ylim([0., 2.])
        ax11.plot(mpc_cols, self.gamma[x_ind])
        gamma_ar_y = []
        for i in range(len(mpc_rows)):
            gamma_ar_y.append(self.gamma[i][y_ind])
        fig12, ax12 = plt.subplots()
        ax12.set_ylim([0., 2.])
        ax12.plot(mpc_rows, gamma_ar_y)
        plt.show()

        """    
        s_data_mpc = [x*0.01 for x in data_mpc]
        s_data_dc = [x*0.01 for x in data_dcm]
        ax.plot(axis_mpc, s_data_mpc)
        ax.plot(axis, s_data_dc)
        ax.plot(axis_mpc, data_dta)
        ax.plot(axis_mpc, data_gamma)
        """

        plt.title("Test")


def main(argv):
    input_file_mpc = argv[0]
    input_file_dcm = argv[1]

    dp = DoseProfiler(input_file_mpc, input_file_dcm)
    dp.do_it()


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
