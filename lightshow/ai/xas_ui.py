import os
import numpy as np
import dash
from dash import dcc
import plotly.express as px

from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from mp_api.client import MPRester

import crystal_toolkit.components as ctc
from crystal_toolkit.helpers.layouts import (
    Box,
    Column,
    Columns,
    Loading
)

from lightshow.ai.models import predict


app = dash.Dash(prevent_initial_callbacks=True, title="OmniXAS@Lightshow.ai")

struct_component = ctc.StructureMoleculeComponent(id="st_vis")
search_component = ctc.SearchComponent(id='mpid_search')
upload_component = ctc.StructureMoleculeUploadComponent(id='file_loader')
xas_plot = dcc.Graph(id='xas_plot')

onmixas_layout = Columns([
        Column(Box([
                    search_component.layout(),
                    upload_component.layout()],
                style={"width": "350px"}), narrow=True),
        Column(Loading(struct_component.layout(size="100%"))),
        Column(xas_plot)
    ],
    desktop_only=False,
    centered=False
)


@app.callback(
    Output(struct_component.id(), "data", allow_duplicate=True),
    Input(search_component.id(), "data")
)
def update_structure_by_mpid(search_mpid: str) -> Structure:
    if not search_mpid:
        raise PreventUpdate
    
    with MPRester() as mpr:
        st = mpr.get_structure_by_material_id(search_mpid)
        print("Struct from material.")

    return st


@app.callback(
    Output(struct_component.id(), "data", allow_duplicate=True),
    Input(upload_component.id(), "data")
)
def update_structure_by_file(upload_data: str) -> Structure:
    if not upload_data:
        raise PreventUpdate
    st = Structure.from_dict(upload_data['data'])
    return st


@app.callback(
    Output("xas_plot", "figure"),
    Input(struct_component.id(), "data")
)
def predict_xas(st_data: str) -> Structure:
    if not st_data:
        raise PreventUpdate
    st = Structure.from_dict(st_data)
    absorbing_site = 'Cu'
    spectroscopy_type = 'FEFF'
    specs = predict(st, absorbing_site, spectroscopy_type)
    spectrum = specs.mean(axis=0)
    ene = np.arange(spectrum.shape[0])
    fig = px.scatter(x=ene, y=ene)
    return fig



ctc.register_crystal_toolkit(app=app, layout=onmixas_layout)

if __name__ == "__main__":
    if "MP_API_KEY" not in os.environ:
        print("Environment variable MP_API_KEY not found, "
              "please set your materials project API key to "
              "this environment variable before running this app")
        exit()
    app.run_server(debug=True, port=8050, host='0.0.0.0')
