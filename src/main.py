import csv
import os
import shutil

from ase.build import bulk
from ase.calculators.calculator import kptdensity2monkhorstpack
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS

from get_input_files import request

FCC = ["Al"]
DIAMOND = ["Ge"]

AUTHORS_LIST = ["Stefano Curtarolo", "Chris Wolverton", "Materials Project", "Miguel Marques", "Bjoern Bieniek", "Oliver Hoffmann", "Kurt Lejaeghere"]

def get_molecule_name(molecule) -> str:
    atoms = {}
    for atom in molecule:
        if atom.symbol in atoms:
            atoms[atom.symbol] = atoms[atom.symbol] + 1
        else:
            atoms[atom.symbol] = 1
    result = ""
    for symbol, count in atoms.items():
        result += symbol + str(count)
    print(result)
    return result

def run(element: str, defect: str) -> list:
    if defect != "":
        fcc = bulk(element, crystalstructure="fcc", a=4, cubic=True)
        super_cell = fcc.repeat((2, 2, 2))
        super_cell[0].symbol = defect
    else:
        if element in FCC:
            super_cell = bulk(element, 'fcc', a=4, cubic=False)
        elif element in DIAMOND:
            super_cell = bulk(element, 'diamond', a=5.658)
        else:
            super_cell = bulk(element, 'bcc', a=4, cubic=False)
    super_cell.calc = Vasp()
    super_cell.calc.read_incar("INCAR")
    super_cell.calc.kpts = tuple(kptdensity2monkhorstpack(super_cell, 5, False))
    optimizer = BFGS(super_cell)
    optimizer.run(fmax=0.02)
    return [get_molecule_name(super_cell), super_cell.get_potential_energy()]

def move_all_files(author_name: str, folder_name: str) -> None:
    author_name = author_name.lower().replace(" ", "_")
    current_dir = os.getcwd()
    destination_folder = os.path.join(current_dir, author_name)
    os.makedirs(destination_folder, exist_ok=True)

    new_folder = os.path.join(destination_folder, folder_name)

    os.makedirs(new_folder, exist_ok=True)
    for file in os.listdir(current_dir):
        file_path = os.path.join(current_dir, file)
        if ".py" in file or not os.path.isfile(file_path):
            continue
        if file == "run.sh":
            continue
        shutil.move(file_path, os.path.join(new_folder, file))


def main():
    result = [["Settings",  "Molecule", "Energy"]]
    element = "Al"
    defect = "Ge"
    for name in AUTHORS_LIST:
        row = [name]
        if request(element, name):
            row.extend(run(element, defect))
            move_all_files(name, element + defect)

            request(element, name)
            row.extend(run(element, ""))
            move_all_files(name, element)
        else:
            row.extend(["", "", "", ""])
        if request(defect, name):
            row.extend(run(defect, ""))
            move_all_files(name, defect)

        result.append(row)
        result.append([])
    with open('energies.csv', mode='w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerows(result)

if __name__ == '__main__':
    main()
