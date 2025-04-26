import os
import shutil
import csv
from ase.build import bulk
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS

from get_input_files import request

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
    return result

def run(element: str, defect: str) -> list:
    fcc = bulk(element, crystalstructure="fcc", a=4, cubic=True)
    super_cell = fcc.repeat((2, 2, 2))
    if defect != "":
        super_cell[0].symbol = defect
    super_cell.calc = Vasp()
    super_cell.calc.read_incar("INCAR")

    optimizer = BFGS(super_cell)
    optimizer.run(fmax=0.02)
    return [get_molecule_name(super_cell), super_cell.get_potential_energy()]

def move_all_files(author_name: str, has_defects: bool):
    author_name = author_name.lower().replace(" ", "_")
    current_dir = os.getcwd()
    destination_folder = os.path.join(current_dir, author_name)
    os.makedirs(destination_folder, exist_ok=True)

    if has_defects:
        new_folder = os.path.join(destination_folder, "defects")
    else:
        new_folder = os.path.join(destination_folder, "clean")
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
        if request(element, name):
            row = [name, *run(element, "")]
            result.append(row)
            move_all_files(name, False)

            request(element, name)
            row = [name, *run(element, defect)]
            result.append(row)
            move_all_files(name, True)
    with open('energies.csv', mode='w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerows(result)

if __name__ == '__main__':
    main()

