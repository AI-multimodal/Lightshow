import sys

import requests
from pymatgen.core.structure import Structure

# Thanks,Chat GPT!


def send_json_payload(server_url, payload):
    """
    Sends a JSON payload to the server and retrieves the response.

    Args:
        server_url (str): The URL of the server endpoint.
        payload (dict): The JSON payload to send.

    Returns:
        dict: The JSON response from the server.
    """
    try:
        # Send the POST request with the JSON payload
        response = requests.post(server_url, json=payload)

        # Raise an exception if the request failed
        response.raise_for_status()

        # Return the JSON response
        return response.json()
    except requests.exceptions.RequestException as e:
        # Handle any errors
        return {"error": str(e)}


if __name__ == "__main__":
    # Server URL (update if the server runs on a different host/port)
    server_url = "http://127.0.0.1:8095/predict_xas"

    path = sys.argv[1]

    payload = {
        "structure_as_json": Structure.from_file(path).to_json(),
        "absorbing_site": "Cu",
        "spectroscopy_type": "FEFF",
    }

    # Send the payload and print the response
    response = send_json_payload(server_url, payload)
    print("Server Response:", response)
