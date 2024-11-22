# Copyright (c) 2022, Materials Virtual Lab
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import pathlib
import sys
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
        type: str,
        pattern=XASBLOCKS_PATH / "{element}_{type}.ckpt",
    ):
        pattern = str(pattern)
        path = pattern.format(element=element, type=type)
        print(f"Loading XASBlock model from {path}")
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
        module.eval()
        module.freeze()
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
        return np.array(node_feat.detach().numpy())

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
            element=element, type=spectroscopy_type
        )

    def _get_feature(self, structure: PymatgenStructure):
        return self.featurizer.featurize(structure)

    def predict(
        self,
        structure: PymatgenStructure,
    ):
        feature = self._get_feature(structure)
        print(feature.shape)
        spectrum = self.model(torch.tensor(feature))
        return spectrum.detach().numpy().squeeze()


def entrypoint():
    # TODO: proper CLI
    spectroscopy_type = sys.argv[1]
    el = sys.argv[2]
    path = sys.argv[3]
    struct = PymatgenStructure.from_file(path)
    site_idxs = [
        ii for ii, site in enumerate(struct.sites) if site.specie.symbol == el
    ]
    if len(site_idxs) == 0:
        raise ValueError(f"element {el} not found in provided structure")
    spec = XASModel(element=el, spectroscopy_type=spectroscopy_type).predict(
        struct
    )
    result = {ii: spec[ii] for ii in site_idxs}
    print(result)


# if __name__ == "__main__":
#     print(AVAILABLE_COMBINATIONS)
# material_structure_file = (
#     "examples/material/mp-1005792/POSCAR"  # TODO: hardcoded path
# )
# strucutre = PymatgenStructure.from_file(material_structure_file)
# spectrum = XASModel(element="Cu", type="FEFF").predict(strucutre)

# %%
