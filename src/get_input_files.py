import bz2
import gzip
import lzma
import os
import re

import requests

BASE_URL = "https://nomad-lab.eu/prod/v1/api/v1"


def request(element: str, author: str, atom_count: int) -> bool:
    try:
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
        response = requests.post(f"{BASE_URL}/entries/query", json=json_body, timeout = 10)
        data = response.json()
        for entry in data["data"]:
            vasp_input_files = get_file_names(entry["files"])
            if vasp_input_files["INCAR"] == "" or vasp_input_files["KPOINTS"] == "":
                continue
            if not get_vasp_inputs(vasp_input_files, entry["entry_id"], atom_count):
                continue
            print(entry["entry_id"])
            return True
        print(f"Failed to find valid files from {author}")
        return False
    except Exception as e:
        print(f"An error occurred: {e}")
        return False

def get_file_names(file_list: list) -> dict:
    vasp_input_files = {"INCAR": "", "KPOINTS": ""}
    for file_path in file_list:
        file_name = file_path.split("/")[-1]
        for key in vasp_input_files:
            if key in file_name:
                vasp_input_files[key] = file_name
    return vasp_input_files

def verify_INCAR(incar: str) -> bool:
    if re.search(r"^\s*MAGNETISM\s*=.*$", incar, re.IGNORECASE | re.MULTILINE):
        return False
    ncore_match = re.search(r"^\s*NCORE\s*=\s*(\d+)", incar, re.IGNORECASE | re.MULTILINE)
    kpar_match = re.search(r"^\s*KPAR\s*=\s*(\d+)", incar, re.IGNORECASE | re.MULTILINE)
    if ncore_match and kpar_match and int(ncore_match.group(1)) % int(kpar_match.group(1)) != 0:
        return False
    tags = re.findall(r"^\s*(XC|GGA|METAGGA)\s*=.*$", incar, re.IGNORECASE | re.MULTILINE)
    if len(tags) > 1:
        return False
    return True



def get_vasp_inputs(vasp_input_files: dict, entry_id: str, atom_count: int) -> bool:
    output_dir = os.path.dirname(os.path.abspath(__file__))
    url = f"{BASE_URL}/entries/{entry_id}/raw/"

    for file_name in vasp_input_files:
        file_url = url + vasp_input_files[file_name]
        response = requests.get(file_url, timeout = 10)
        text = decompress_text(response, vasp_input_files[file_name])
        if "INCAR" in file_name:
            if not verify_INCAR(text):
                return False
            text = re.sub(r'^\s*(NBANDS|ICHARG|METAGGA)\s*=.*$', '', text, flags=re.MULTILINE | re.IGNORECASE)
            text = re.sub(
                r'^(?P<prefix>\s*MAGMOM\s*=\s*)(?P<count>\d+)(?P<suffix>\*.*)$',
                lambda m: f"{m.group('prefix')}{atom_count}{m.group('suffix')}",
                text,
                flags=re.MULTILINE | re.IGNORECASE
            )
            if not re.search(r'^\s*(XC|GGA|METAGGA)\s*=.*$', text, re.IGNORECASE | re.MULTILINE):
                text += "\nGGA = PE"

        if response.status_code == 200:
            file_path = os.path.join(output_dir, file_name)
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(text)
            print(f"Downloaded {file_name} to {file_path}")
        else:
            print(f"Failed to download {file_name}: {response.status_code}")
            return False
    return True


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
