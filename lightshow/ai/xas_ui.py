import os
import numpy as np
import dash
from dash import dcc, html
import plotly.express as px

from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
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

all_elements = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']
ene_start = {'Ti': 4964.504, 'V': 5464.097, 'Cr': 5989.168, 'Mn': 6537.886, 
             'Fe': 7111.23, 'Co': 7709.282, 'Ni': 8332.181, 'Cu': 8983.173}
ene_grid = {el: np.linspace(start, start + 35, 141) for el, start in ene_start.items()}
xas_model_names = [f'{el} FEFF' for el in all_elements] + ['Ti VASP', 'Cu VASP']
absorber_dropdown = dcc.Dropdown(xas_model_names, value='Ti VASP', id='absorber')

onmixas_layout = Columns([
        Column(Box([Loading(search_component.layout()),
                    Loading(upload_component.layout()),
                    html.Br(), html.Br(),
                    html.Div("Element and Theory:"),
                    Loading(absorber_dropdown)],
                style={"width": "350px"}), narrow=True),
        Column(Loading(struct_component.layout(size="100%"))),
        Column(xas_plot)
    ],
    desktop_only=False,
    centered=False
)


@app.callback(
    Output(struct_component.id(), "data", allow_duplicate=True),
    Input(search_component.id(), "data"),
    State('absorber', 'value')
)
def update_structure_by_mpid(search_mpid: str, el_type) -> Structure:
    if not search_mpid:
        raise PreventUpdate
    
    with MPRester() as mpr:
        st = mpr.get_structure_by_material_id(search_mpid)
        print("Struct from material.")
    st_dict = decorate_structure_with_xas(st, el_type)
    return st_dict


def decorate_structure_with_xas(st: Structure, el_type):
    absorbing_site, spectroscopy_type = el_type.split(' ')
    st_dict = st.as_dict()
    if absorbing_site in st.composition:
        specs = predict(st, absorbing_site, spectroscopy_type)
        st_dict['xas'] = specs
    else:
        st_dict['xas'] = {}
    return st_dict


@app.callback(
    Output(struct_component.id(), "data", allow_duplicate=True),
    Input(upload_component.id(), "data"),
    State('absorber', 'value')
)
def update_structure_by_file(upload_data: dict, el_type) -> Structure:
    if not upload_data:
        raise PreventUpdate
    st = Structure.from_dict(upload_data['data'])
    st_dict = decorate_structure_with_xas(st, el_type)
    return st_dict


@app.callback(
    Output("xas_plot", "figure", allow_duplicate=True),
    Input(struct_component.id(), "data"),
    State('absorber', 'value')
)
def predict_average_xas(st_data: dict, el_type) -> Structure:
    if not st_data:
        raise PreventUpdate
    specs = st_data['xas']
    if len(specs) == 0:
        fig = build_figure(None, el_type, is_average=True, no_element=True, sel_mismatch=False)
    else:
        specs = np.array(list(specs.values()))
        spectrum = specs.mean(axis=0)
        fig = build_figure(spectrum, el_type, is_average=True, no_element=False, sel_mismatch=False)
    return fig


def build_figure(spectrum, el_type, is_average, no_element, sel_mismatch):
    element = el_type.split(" ")[0]
    if spectrum is None:
        ene = None
    else:
        ene = ene_grid[element]
    if no_element:
        title = f"This structure doesn't contain {element}"
    elif sel_mismatch:
        title = f"The selected atom is not a {element} atom"
    elif is_average:
        title = f'Average K-edge XANES Spectrum of {el_type}'
    else:
        title = f'K-edge XANES Spectrum for the selected {element} atom'
    fig = px.scatter(x=ene, y=spectrum, title=title, 
                     labels={'x': "Energy (eV)", "y": "Absorption"})
    return fig


@app.callback(
    Output("xas_plot", "figure", allow_duplicate=True),
    Input(struct_component.id('scene'), "selectedObject"),
    State(struct_component.id(), 'data'),
    State('absorber', 'value')
)
def predict_site_specific_xas(sel, st_data, el_type) -> Structure:
    specs = st_data['xas']
    element = el_type.split(' ')[0]
    if len(specs) == 0:
        fig = build_figure(None, el_type, is_average=False, no_element=True, sel_mismatch=False)
    elif len(sel) == 0:
        specs = np.array(list(specs.values()))
        spectrum = specs.mean(axis=0)
        fig = build_figure(spectrum, el_type, is_average=True, no_element=False, sel_mismatch=False)
    else:
        st = Structure.from_dict(st_data)
        spheres = st._get_sites_to_draw()
        spheres = list(spheres)
        i_sphere = int(sel[0]['id'].split('--')[-1])
        cur_sphere = spheres[i_sphere]
        i_site = cur_sphere[0]
        if st[i_site].specie.symbol != element:
            fig = build_figure(None, el_type, is_average=False, no_element=False, sel_mismatch=True)
        else:
            spectrum = np.array(specs[str(i_site)])
            fig = build_figure(spectrum, el_type, is_average=False, no_element=False, sel_mismatch=False)
    return fig


@app.callback(
    Output(struct_component.id(), "data", allow_duplicate=True),
    Input('absorber', 'value'),
    State(struct_component.id(), "data")
)
def update_structure_by_mpid(el_type, st_data) -> Structure:
    st = Structure.from_dict(st_data)
    st_dict = decorate_structure_with_xas(st, el_type)
    return st_dict
    

ctc.register_crystal_toolkit(app=app, layout=onmixas_layout)

if __name__ == "__main__":
    if "MP_API_KEY" not in os.environ:
        print("Environment variable MP_API_KEY not found, "
              "please set your materials project API key to "
              "this environment variable before running this app")
        exit()
    app.run_server(debug=False, port=8050, host='0.0.0.0')
