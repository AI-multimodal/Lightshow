import dash
from dash import dcc
import plotly.express as px

from dash.dependencies import Input, Output
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

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
        Column([search_component.layout(),
                upload_component.layout()]),
        Column(Loading(struct_component.layout(size="100%"))),
        Column(xas_plot)
    ],
    desktop_only=False,
    centered=False
)

ctc.register_crystal_toolkit(app=app, layout=onmixas_layout)

if __name__ == "__main__":
    app.run_server(debug=True, port=8050, host='0.0.0.0')
