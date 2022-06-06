import json
import bz2
import base64
import hashlib
import re


pspLibrary = "SSSP_precision"
# pspJsonOut = 'SSSP_precision_pseudos.json'


pspJsonIn = pspLibrary + ".json"

with open(pspJsonIn, "r") as f:
    inJSON = json.load(f)

pspJsonIn = pspLibrary + "_pseudos.json"

with open(pspJsonIn, "r") as f:
    pspJSON = json.load(f)

for element in inJSON:
    print(element)

    print(inJSON[element]["filename"])
    fileName = inJSON[element]["filename"]

    if fileName not in pspJSON:
        print("Incomplete psp database")
        exit()

    pspString = bz2.decompress(base64.b64decode(pspJSON[fileName]))

    print("Expected hash:  " + inJSON[element]["md5"])
    print("Resultant hash: " + hashlib.md5(pspString).hexdigest())

    m = re.search(
        'z_valence="\s*(\d+\.?\d*([Ed][+-]?\d+)?)\s*"',
        pspString.decode(),
        re.IGNORECASE,
    )
    if m:
        print(int(float(m.group(1))), " :    ", m.group(0))
        zval = int(float(m.group(1)))
    else:
        m = re.search(
            "\s*(\d+\.?\d*)\s+z valence", pspString.decode(), re.IGNORECASE
        )
        print(int(float(m.group(1))), " :    ", m.group(0))
        #        print( int(float(m.group(1))))

        zval = int(float(m.group(1)))

    inJSON[element]["Z_val"] = zval


with open("zval.json", "w") as f:
    f.write(
        json.dumps(inJSON, sort_keys=True, indent=3, separators=(",", ": "))
    )
