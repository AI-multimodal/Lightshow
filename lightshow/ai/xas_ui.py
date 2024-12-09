import os
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

app = dash.Dash(prevent_initial_callbacks=True, title="OmniXAS Prediction @ Lightshow.ai")

struct_component = ctc.StructureMoleculeComponent(id="st_vis")
search_component = ctc.SearchComponent(id='mpid_search')
upload_component = ctc.StructureMoleculeUploadComponent(id='file_loader')
xas_plot = dcc.Graph()

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
    Output(struct_component.id(), "data"),
    Input(search_component.id(), "data")
)
def update_structure_by_mpid(search_mpid: str) -> Structure:
    if not search_mpid:
        raise PreventUpdate
    
    with MPRester() as mpr:
        struct = mpr.get_structure_by_material_id(search_mpid)
        print("Struct from material.")

    return struct


ctc.register_crystal_toolkit(app=app, layout=onmixas_layout)

if __name__ == "__main__":
    if "MP_API_KEY" not in os.environ:
        print("Environment variable MP_API_KEY not found, "
              "please set your materials project API key to "
              "this environment variable before running this app")
        exit()
    app.run_server(debug=True, port=8050, host='0.0.0.0')
