import bz2
import gzip
import lzma
import os

import requests

BASE_URL = "https://nomad-lab.eu/prod/v1/api/v1"


def request(element: str, author: str) -> bool:
    json_body = {
        "query": {
            "results.method.simulation.program_name:any": ["VASP"],
            "authors.name:any": [author],
            "results.properties.available_properties:all": ["geometry_optimization"],
            "results.material.elements_exclusive": element,
        },
        "pagination": {"page_size": 10, "page": 1},
        "required": {"exclude": ["quantities", "sections"]},
    }
    response = requests.post(f"{BASE_URL}/entries/query", json=json_body)
    data = response.json()
    for entry in data["data"]:
        vasp_input_files = get_file_names(entry["files"])
        if vasp_input_files["INCAR"] == "" or vasp_input_files["KPOINTS"] == "":
            continue
        get_vasp_inputs(vasp_input_files, entry["entry_id"])
        print(entry["entry_id"])
        return True
    return False


def get_file_names(file_list: list) -> dict:
    vasp_input_files = {"INCAR": "", "KPOINTS": ""}
    for file_path in file_list:
        file_name = file_path.split("/")[-1]
        for key in vasp_input_files:
            if key in file_name:
                vasp_input_files[key] = file_name
    return vasp_input_files


def get_vasp_inputs(vasp_input_files: dict, entry_id: str) -> None:
    output_dir = os.path.dirname(os.path.abspath(__file__))
    url = f"{BASE_URL}/entries/{entry_id}/raw/"

    for file_name in vasp_input_files:
        file_url = url + vasp_input_files[file_name]
        response = requests.get(file_url)
        text = decompress_text(response, vasp_input_files[file_name])
        if "INCAR" in file_name:
            text += "\nGGA = PE"

        if response.status_code == 200:
            file_path = os.path.join(output_dir, file_name)
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(text)
            print(f"Downloaded {file_name} to {file_path}")
        else:
            print(f"Failed to download {file_name}: {response.status_code}")

def decompress_text(response, file_name: str) -> str:
    text = (
        gzip.decompress(response.content).decode("utf-8")
        if file_name.endswith(".gz")
        else (
            bz2.decompress(response.content).decode("utf-8")
            if file_name.endswith(".bz2")
            else (
                lzma.decompress(response.content).decode("utf-8")
                if file_name.endswith(".xz")
                else response.text
            )
        )
    )
    return text
