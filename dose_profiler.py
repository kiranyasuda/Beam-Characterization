import streamlit as st
import numpy as np
import dose_dist_eval as dde
import ascii_reader as ar
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


class MainDoseProfiler:

    def __init__(self):
        st.set_page_config(layout="wide")
        self.mpc_file = ""
        self.dcm_file = ""
        self.ascii_file = ""
        self.dp = None
        self.ar = None
        self.profile_dir = 'C:/Users/kiran/Profiles/'

    def setup(self):
        container = st.container()
        uploaded_mapcheck = container.file_uploader("Choose a Mapcheck File ", "txt")
        uploaded_dicom = container.file_uploader("Choose a Dicom File ", "dcm")
        uploaded_ascii = container.file_uploader("Choose an ASCII File", "asc")
        if uploaded_mapcheck is not None and uploaded_dicom is not None and uploaded_ascii is not None:
            st.session_state['setup1'] = uploaded_mapcheck
            st.session_state['setup2'] = uploaded_dicom
            st.session_state['setup3'] = uploaded_ascii
            self.dp = dde.DoseProfiler(uploaded_mapcheck,
                                       uploaded_dicom)
            self.ar = ar.AsciiProfile(uploaded_ascii)
            container.success("Data has been uploaded")
        if 'setup1' in st.session_state:
            if st.button("Calculate"):
                with st.spinner('Calculating, please wait.'):
                    self.dp.calculate(None)
                container.success('Calculations Complete')

    def mapcheck(self):
        col1, col2 = st.columns([0.5, 0.5], gap="small")
        mapcheck_slice_axis = col1.selectbox("Mapcheck Slice Axis: ",
                                             ["X", "Y"])
        if self.dp:
            slice_axis = None
            if mapcheck_slice_axis == "X":
                slice_axis = sorted(self.dp.mpc.col_params)
            elif mapcheck_slice_axis == "Y":
                slice_axis = sorted(self.dp.mpc.row_params)
            mapcheck_slice_coordinate = col2.selectbox("Mapcheck Slice Coordinate: ",
                                                       slice_axis)
            if col1.button("Show Mapcheck Heat Map"):
                figure1 = self.create_mpc_heatmap()
                st.session_state['mapcheck1'] = figure1
            if col2.button("Show Mapcheck Slice"):
                figure2 = self.mpc_on_select(slice_axis, mapcheck_slice_coordinate)
                st.session_state['mapcheck2'] = figure2

            if 'mapcheck1' in st.session_state:
                col1.write(st.session_state['mapcheck1'])

            if 'mapcheck2' in st.session_state:
                col2.write(st.session_state['mapcheck2'])

    def dcm(self):
        col1, col2 = st.columns([0.5, 0.5], gap="small")
        dicom_slice_axis = col1.selectbox("Dicom Slice Axis: ",
                                          ["X", "Y"])
        if self.dp:
            if dicom_slice_axis == "X":
                dcm_slice_axis = sorted(self.dp.dcm.dcm_cols)
            else:
                dcm_slice_axis = sorted(self.dp.dcm.dcm_rows)
            dicom_slice_coordinate = col2.selectbox("Dicom Slice Coordinate: ",
                                                    dcm_slice_axis)
            if col1.button("Show DCM Heat Map"):
                figure1 = self.create_dcm_heatmap()
                st.session_state['dcm1'] = figure1
            if col2.button("Show DCM Slice"):
                figure2 = self.dcm_on_select(dicom_slice_axis, dicom_slice_coordinate)
                st.session_state['dcm2'] = figure2

            if 'dcm1' in st.session_state:
                col1.write(st.session_state['dcm1'])

            if 'dcm2' in st.session_state:
                col2.write(st.session_state['dcm2'])

    def overlay(self):
        col1, col2 = st.columns([0.5, 0.5], gap="small")
        overlay_slice_axis = col1.selectbox("Overlay Slice Axis: ",
                                            ["X", "Y"])
        if self.dp:
            if max(self.dp.mpc.col_params) > max(self.dp.mpc.row_params):
                overlay_slice_coord = col1.selectbox("Overlay Slice Coordinate", sorted(self.dp.mpc.row_params))
            else:
                overlay_slice_coord = col1.selectbox("Overlay Slice Coordinate", sorted(self.dp.mpc.col_params))
            if col1.button("View Overlayed Slices"):
                fig1 = self.overlay_on_select(overlay_slice_axis, overlay_slice_coord)
                col2.write(fig1)

    def dose_difference(self):
        col1, col2 = st.columns([0.5, 0.5], gap="small")
        difference_slice_axis = col1.selectbox("Dose Difference Slice Axis: ",
                                               ["X", "Y"])
        if self.dp:
            if difference_slice_axis == "X":
                slice_axis = sorted(self.dp.mpc.col_params)
            else:
                slice_axis = sorted(self.dp.mpc.row_params)
            difference_slice_coordinate = col2.selectbox("Dose Difference Slice Coordinate: ",
                                                         slice_axis)
            if col1.button("Show DD Heat Map"):
                figure1 = self.create_dd_heatmap()
                st.session_state['dd1'] = figure1
            if col2.button("Show DD Slice"):
                figure2 = self.dd_on_select(difference_slice_axis, difference_slice_coordinate)
                st.session_state['dd2'] = figure2

            if 'dd1' in st.session_state:
                col1.write(st.session_state['dd1'])

            if 'dd2' in st.session_state:
                col2.write(st.session_state['dd2'])

    def dta(self):
        col1, col2 = st.columns([0.5, 0.5], gap="small")
        dta_slice_axis = col1.selectbox("DTA Slice Axis: ",
                                        ["X", "Y"])
        if self.dp:
            if dta_slice_axis == "X":
                slice_axis = sorted(self.dp.mpc.col_params)
            else:
                slice_axis = sorted(self.dp.mpc.row_params)
            dta_slice_coordinate = col2.selectbox("DTA Slice Coordinate: ",
                                                  slice_axis)
            if col1.button("Show DTA Heat Map"):
                figure1 = self.create_dta_heatmap()
                st.session_state['dta1'] = figure1
            if col2.button("Show DTA Slice"):
                figure2 = self.dta_on_select(dta_slice_axis, dta_slice_coordinate)
                st.session_state['dta2'] = figure2

            if 'dta1' in st.session_state:
                col1.write(st.session_state['dta1'])

            if 'dta2' in st.session_state:
                col2.write(st.session_state['dta2'])

    def gamma(self):
        col1, col2 = st.columns([0.5, 0.5], gap="small")
        gamma_slice_axis = col1.selectbox("Gamma Slice Axis: ",
                                          ["X", "Y"])
        if self.dp:
            if gamma_slice_axis == "X":
                slice_axis = sorted(self.dp.mpc.col_params)
            else:
                slice_axis = sorted(self.dp.mpc.row_params)
            gamma_slice_coordinate = col2.selectbox("Gamma Slice Coordinate: ",
                                                    slice_axis)

            if col1.button("Show Gamma Heatmap"):
                figure1 = self.create_gamma_heatmap()
                st.session_state['gamma1'] = figure1
            if col2.button("Show Gamma Slice"):
                figure2 = self.gamma_on_select(gamma_slice_axis, gamma_slice_coordinate)
                st.session_state['gamma2'] = figure2

            if 'gamma1' in st.session_state:
                col1.write(st.session_state['gamma1'])

            if 'gamma2' in st.session_state:
                col2.write(st.session_state['gamma2'])

    def depth_dose(self):
        col1, col2 = st.columns(2, gap="small")

        if self.dp:
            field_sizes = set([data_sect['fsz'] for data_sect in self.ar.data_sects])
            field_sizes = sorted([[int(x[0]), int(x[1])] for x in [x.split('x') for x in field_sizes]])
            field_sizes = ['x'.join([str(x[0]), str(x[1])]) for x in field_sizes]
            field_size = col1.selectbox("Choose Field Size", field_sizes)
            field_depths = set([data_sect['depth'] for data_sect in self.ar.data_sects])
            field_depths = sorted([int(float(x)) for x in field_depths])
            field_depth = col1.selectbox("Choose Field Depth", field_depths)
            field_axis = col1.selectbox("Choose Field Axis", ['X', 'Y', 'Z'])
            figure = None
            if col1.button("Generate Graph"):
                figure = self.ar.graph(field_size, field_depth, field_axis)
                st.session_state['depthdose'] = figure
            if 'depthdose' in st.session_state:
                if figure:
                    col2.write(st.session_state['depthdose'])
                else:
                    col2.write('Choose New Axis')

    def create_mpc_heatmap(self):
        figure = plt.Figure(figsize=(6, 5), dpi=100)
        ax = figure.add_subplot(111)
        if self.dp is not None:
            mpc = self.dp.mpc
            mpc_cols = mpc.col_params
            mpc_rows = [row[0] for row in mpc.dose_interpolated]
            dde.heatmap(np.array(mpc.dose_interpolated), mpc_rows, mpc_cols, 10, ax)
            return figure

    def mpc_on_select(self, axis, coord):
        if axis == 'X':
            x_idx = self.dp.mpc.col_params.index(float(coord)) + 2
            y_idx = None
            axis_array = self.dp.mpc.row_params
        else:
            x_idx = None
            y_idx = self.dp.mpc.row_params.index(float(coord))
            axis_array = self.dp.mpc.col_params

        data_array = [x * self.dp.norm for x in self.dp.mpc.select_slice_by_index(x_idx, y_idx)]

        fig = Figure(figsize=(5, 5), dpi=100)
        ax1 = fig.add_subplot(111)
        if axis == 'X':
            ax1.set_xlabel("Y (cm)")
        else:
            ax1.set_xlabel("X (cm)")
        ax1.set_ylabel("cGray")
        ax1.plot(axis_array, data_array)
        return fig

    def create_dcm_heatmap(self):
        figure = plt.Figure(figsize=(6, 5), dpi=100)
        ax = figure.add_subplot(111)
        if self.dp is not None:
            dcm = self.dp.dcm
            dcm_cols = dcm.dcm_cols
            dcm_rows = dcm.dcm_rows
            dde.heatmap(dcm.data_array, dcm_rows, dcm_cols, 40, ax)
            return figure

    def dcm_on_select(self, axis, coord):
        if axis == 'X':
            x_idx = self.dp.dcm.dcm_cols.index(float(coord)) + 2
            y_idx = None
            axis_array = self.dp.dcm.dcm_rows
        else:
            x_idx = None
            y_idx = self.dp.dcm.dcm_rows.index(float(coord))
            axis_array = self.dp.dcm.dcm_cols

        data_array = self.dp.dcm.select_slice_by_index(x_idx, y_idx)

        fig = Figure(figsize=(5, 5), dpi=100)
        ax1 = fig.add_subplot(111)
        if axis == 'X':
            ax1.set_xlabel("Y (cm)")
        else:
            ax1.set_xlabel("X (cm)")
        ax1.set_ylabel("cGray")
        ax1.plot(axis_array, data_array, color="C1")
        return fig

    def overlay_on_select(self, axis, coord):

        if axis == 'X':
            x_slice = float(coord)
            y_slice = None
            x_idx = self.dp.mpc.col_params.index(float(coord)) + 2
            y_idx = None
            axis_array = self.dp.mpc.row_params
        else:
            x_slice = None
            y_slice = float(coord)
            x_idx = None
            y_idx = self.dp.mpc.row_params.index(float(coord))
            axis_array = self.dp.mpc.col_params

        data_array = [x * self.dp.norm for x in self.dp.mpc.select_slice_by_index(x_idx, y_idx)]

        fig1 = Figure(figsize=(6, 6), dpi=100)
        ax1 = fig1.add_subplot(111, label="1")
        ax1.set_xlim(-15., 15.)
        y_max = max(data_array) * 1.2
        ax1.set_ylim(0., y_max)
        if axis == 'X':
            ax1.set_xlabel("Y (cm)")
        else:
            ax1.set_xlabel("X (cm)")
        ax1.set_ylabel("cGray")
        ax1.scatter(axis_array, data_array, s=2, color="C0")

        coord, c_array, d_array = self.dp.dcm.select_slice(x_slice, y_slice)
        ax2 = fig1.add_subplot(111, label="2", frame_on=False)
        ax2.set_xlim(-15., 15.)
        ax2.set_ylim(0., y_max)
        ax2.scatter(c_array, d_array, s=2, color="C1")
        return fig1

    def create_dd_heatmap(self):
        figure = plt.Figure(figsize=(6, 5), dpi=100)
        ax = figure.add_subplot(111)
        if self.dp is not None:
            dd = self.dp.dd
            mpc = self.dp.mpc
            mpc_cols = mpc.col_params
            mpc_rows = [row[0] for row in mpc.dose_interpolated]
            dde.heatmap(dd, mpc_rows, mpc_cols, 10, ax)
            return figure

    def dd_on_select(self, axis, coord):
        if axis == 'X':
            x_idx = self.dp.mpc.col_params.index(float(coord)) + 2
            y_idx = None
            axis_array = self.dp.mpc.row_params
        else:
            x_idx = None
            y_idx = self.dp.mpc.row_params.index(float(coord))
            axis_array = self.dp.mpc.col_params

        data_array = self.dp.select_dd_slice_by_index(x_idx, y_idx)

        fig = Figure(figsize=(5, 5), dpi=100)
        ax1 = fig.add_subplot(111)
        # ax1.set_ylim([-0.1, 1.2])
        if axis == 'X':
            ax1.set_xlabel("Y (cm)")
        else:
            ax1.set_xlabel("X (cm)")
        ax1.set_ylabel("Dose Difference")
        ax1.plot(axis_array, data_array)
        return fig

    def create_dta_heatmap(self):
        figure = plt.Figure(figsize=(6, 5), dpi=100)
        ax = figure.add_subplot(111)
        if self.dp is not None:
            dta = self.dp.dta
            mpc = self.dp.mpc
            mpc_cols = mpc.col_params
            mpc_rows = [row[0] for row in mpc.dose_interpolated]
            dde.heatmap(dta, mpc_rows, mpc_cols, 10, ax)
            return figure

    def dta_on_select(self, axis, coord):
        if axis == 'X':
            x_idx = self.dp.mpc.col_params.index(float(coord)) + 2
            y_idx = None
            axis_array = self.dp.mpc.row_params
        else:
            x_idx = None
            y_idx = self.dp.mpc.row_params.index(float(coord))
            axis_array = self.dp.mpc.col_params

        data_array = self.dp.select_dta_slice_by_index(x_idx, y_idx)

        fig = Figure(figsize=(5, 5), dpi=100)
        ax1 = fig.add_subplot(111)
        ax1.set_ylim([0., 1.2])
        if axis == 'X':
            ax1.set_xlabel("Y (cm)")
        else:
            ax1.set_xlabel("X (cm)")
        ax1.set_ylabel("Distance-to-Agreement (cm)")
        ax1.plot(axis_array, data_array)
        return fig

    def create_gamma_heatmap(self):
        figure = plt.Figure(figsize=(6, 5), dpi=100)
        ax = figure.add_subplot(111)
        if self.dp is not None:
            gamma = self.dp.gamma
            mpc = self.dp.mpc
            mpc_cols = mpc.col_params
            mpc_rows = [row[0] for row in mpc.dose_interpolated]
            dde.heatmap(gamma, mpc_rows, mpc_cols, 10, ax)
            return figure

    def gamma_on_select(self, axis, coord):
        if axis == 'X':
            x_idx = self.dp.mpc.col_params.index(float(coord)) + 2
            y_idx = None
            axis_array = self.dp.mpc.row_params
        else:
            x_idx = None
            y_idx = self.dp.mpc.row_params.index(float(coord))
            axis_array = self.dp.mpc.col_params

        data_array = self.dp.select_gamma_slice_by_index(x_idx, y_idx)

        fig = Figure(figsize=(5, 5), dpi=100)
        ax1 = fig.add_subplot(111)
        ax1.set_ylim([0., 1.2])
        if axis == 'X':
            ax1.set_xlabel("Y (cm)")
        else:
            ax1.set_xlabel("X (cm)")
        ax1.set_ylabel("Gamma")
        ax1.plot(axis_array, data_array)
        return fig
