import streamlit as st
import numpy as np
import dose_dist_eval as dde
import matplotlib.pyplot as plt
from streamlit_option_menu import option_menu
from dose_profiler2 import MainDoseProfiler


if 'dp' not in st.session_state:
    dp = MainDoseProfiler()
    st.session_state['dp'] = dp

class Sidebar:
    def __init__(self):
        with st.sidebar:
            self.app = option_menu(
                menu_title="Main Menu",
                options=["Setup", "Mapcheck", "Dicom", "Overlay", "Dose Difference", "DTA", "Gamma", "Depth Dose"],

            )

    def run(self):
        dp = st.session_state['dp']
        if self.app == "Setup":
            dp.setup()
        if self.app == "Mapcheck":
            dp.mapcheck()
        if self.app == "Dicom":
            dp.dcm()
        if self.app == "Overlay":
            dp.overlay()
        if self.app == "Dose Difference":
            dp.dose_difference()
        if self.app == "DTA":
            dp.dta()
        if self.app == "Gamma":
            dp.gamma()
        if self.app == "Depth Dose":
            dp.depth_dose()
        st.session_state['dp'] = dp


sidebar = Sidebar()
sidebar.run()

