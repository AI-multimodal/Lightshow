import pathlib
from functools import cache
from typing import List

import numpy as np
import torch
from lightning import LightningModule
from matgl import load_model
from matgl.ext.pymatgen import Structure2Graph
from matgl.graph.compute import (
    compute_pair_vector_and_distance,
    compute_theta_and_phi,
    create_line_graph,
)
from matgl.utils.cutoff import polynomial_cutoff
from pymatgen.core import Structure as PymatgenStructure
from torch import nn

PARENT_DIRECTORY = pathlib.Path(__file__).parent.resolve()
XASBLOCKS_PATH = PARENT_DIRECTORY / "models" / "xasblock" / "v1.1.1"
AVAILABLE_COMBINATIONS = [f.stem for f in XASBLOCKS_PATH.glob("*.ckpt")]


class XASBlock(nn.Sequential):
    DROPOUT = 0.5

    def __init__(self, input_dim: int, hidden_dims: List[int], output_dim: int):
        dims = [input_dim] + hidden_dims + [output_dim]
        layers = []
        for i, (w1, w2) in enumerate(zip(dims[:-1], dims[1:])):
            layers.append(nn.Linear(w1, w2))
            if i < len(dims) - 2:  # not last layer
                layers.append(nn.BatchNorm1d(w2))
                layers.append(nn.SiLU())
                layers.append(nn.Dropout(self.DROPOUT))
            else:
                layers.append(nn.Softplus())  # last layer
        super().__init__(*layers)


class XASBlockModule(LightningModule):
    def __init__(self, model: nn.Module):
        super().__init__()
        self.model = model

    def forward(self, x):
        return self.model(x)

    @classmethod
    def load(
        cls,
        element: str,
        spectroscopy_type: str,
        pattern=XASBLOCKS_PATH / "{element}_{type}.ckpt",
    ):
        pattern = str(pattern)
        path = pattern.format(element=element, type=spectroscopy_type)
        model = XASBlock(
            input_dim=64,
            hidden_dims=[
                500,
                500,
                550,
            ],  # TODO: hardcoded dims for version v1.1.1
            output_dim=141,
        )
        module = cls.load_from_checkpoint(checkpoint_path=path, model=model)
        return module


class M3GNetFeaturizer:
    def __init__(self, model=None, n_blocks=None):
        self.model = model or M3GNetFeaturizer._load_m3gnet()
        self.model.eval()
        self.n_blocks = n_blocks or self.model.n_blocks

    def featurize(
        self,
        structure: PymatgenStructure,
    ):
        graph_converter = Structure2Graph(
            self.model.element_types, self.model.cutoff
        )
        g, state_attr = graph_converter.get_graph(structure)

        node_types = g.ndata["node_type"]
        bond_vec, bond_dist = compute_pair_vector_and_distance(g)

        g.edata["bond_vec"] = bond_vec.to(g.device)
        g.edata["bond_dist"] = bond_dist.to(g.device)

        with torch.no_grad():
            expanded_dists = self.model.bond_expansion(g.edata["bond_dist"])

            l_g = create_line_graph(g, self.model.threebody_cutoff)

            l_g.apply_edges(compute_theta_and_phi)
            g.edata["rbf"] = expanded_dists
            three_body_basis = self.model.basis_expansion(l_g)
            three_body_cutoff = polynomial_cutoff(
                g.edata["bond_dist"], self.model.threebody_cutoff
            )
            node_feat, edge_feat, state_feat = self.model.embedding(
                node_types, g.edata["rbf"], state_attr
            )

            for i in range(self.n_blocks):
                edge_feat = self.model.three_body_interactions[i](
                    g,
                    l_g,
                    three_body_basis,
                    three_body_cutoff,
                    node_feat,
                    edge_feat,
                )
                edge_feat, node_feat, state_feat = self.model.graph_layers[i](
                    g, edge_feat, node_feat, state_feat
                )
        res = np.array(node_feat.detach().numpy())
        return res

    @cache
    @staticmethod
    def _load_m3gnet(path=PARENT_DIRECTORY / "models/M3GNet-MP-2021.2.8-PES"):
        model = load_model(path).model
        model.eval()
        return model


class XASModel:
    featurizer = M3GNetFeaturizer()

    def __init__(self, element: str, spectroscopy_type: str):
        self.element = element
        self.spectroscopy_type = spectroscopy_type
        self.model = XASBlockModule.load(
            element=element, spectroscopy_type=spectroscopy_type
        )
        self.model.eval()

    def _get_feature(self, structure: PymatgenStructure):
        return self.featurizer.featurize(structure)

    def predict(
        self,
        structure: PymatgenStructure,
    ):
        with torch.no_grad():
            feature = self._get_feature(structure)
            # Normalize by the proper amount
            feature = feature * 1000.0
            device = (
                torch.device("cuda") if torch.cuda.is_available() else "cpu"
            )
            feature = torch.tensor(feature, device=device)
            spectrum = self.model(feature)
        spectrum = spectrum.detach().cpu().numpy().squeeze()
        return spectrum


def predict(structure, absorbing_site, spectroscopy_type):
    site_idxs = [
        ii
        for ii, site in enumerate(structure.sites)
        if site.specie.symbol == absorbing_site
    ]
    if len(site_idxs) == 0:
        raise ValueError(
            f"element {absorbing_site} not found in provided structure"
        )
    spec = XASModel(
        element=absorbing_site, spectroscopy_type=spectroscopy_type
    ).predict(structure)
    result = {ii: spec[ii] for ii in site_idxs}
    return result
