import json

from flask import Flask, jsonify, request
from pymatgen.core.structure import Structure

from lightshow.ai.models import predict

app = Flask(__name__)


# Example: Echo back the received data with a status
def predict_xas(payload):
    structure_as_json = payload["structure_as_json"]
    absorbing_site = payload["absorbing_site"]
    spectroscopy_type = payload["spectroscopy_type"]
    d = json.loads(structure_as_json)
    structure = Structure.from_dict(d)
    return predict(structure, absorbing_site, spectroscopy_type)


@app.route("/predict_xas", methods=["POST"])
def process_predict_xas():
    payload = request.get_json()
    if payload is None:
        return jsonify({"error": "Invalid JSON payload"}), 400
    result = predict_xas(payload)
    result = {k: v.tolist() for k, v in result.items()}
    return jsonify(result)


def serve():
    app.run(host="0.0.0.0", port=8095)
