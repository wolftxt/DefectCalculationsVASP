import csv
import os
import shutil

from ase import Atoms
from ase.build import bulk
from ase.calculators.calculator import kptdensity2monkhorstpack
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS

from get_input_files import request, set_incar_tag

FCC = ["Al"]
DIAMOND = ["Ge"]

AUTHORS_LIST = ["Materials Project"]

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

def set_calc(cell: Atoms, kpts: int, ismear: float) -> None:
    set_incar_tag("ISMEAR", str(ismear))
    set_incar_tag("LWAVE", ".FALSE.")
    set_incar_tag("LAECHG", ".FALSE.")
    set_incar_tag("LCHARG", ".FALSE.")
    cell.calc = Vasp()
    cell.calc.read_incar("INCAR")
    cell.calc.kpts = tuple(kptdensity2monkhorstpack(cell, kpts, False))

def run(element: str, defect: str, kpts: int) -> list:
    if defect != "":
        ismear = 1  # Bad temporary measure to ensure correct ismear
        fcc = bulk(element, crystalstructure="fcc", a=4, cubic=True)
        super_cell = fcc.repeat((2, 2, 2))
        super_cell[0].symbol = defect
    else:
        if element in FCC:
            ismear = 1 # Bad temporary measure to ensure correct ismear
            super_cell = bulk(element, 'fcc', a=4, cubic=False)
        elif element in DIAMOND:
            ismear = 0  # Bad temporary measure to ensure correct ismear
            super_cell = bulk(element, 'diamond', a=5.658)
        else:
            ismear = 1  # Bad temporary measure to ensure correct ismear
            super_cell = bulk(element, 'bcc', a=4, cubic=False)
    set_calc(super_cell, kpts, ismear)
    optimizer = BFGS(super_cell)
    optimizer.run(fmax=0.05)

    set_calc(super_cell, kpts, ismear)
    optimizer = BFGS(super_cell)
    optimizer.run(fmax=0.01)

    set_calc(super_cell, kpts, -5)
    return [get_molecule_name(super_cell), super_cell.get_potential_energy()]

def move_all_files(author_name: str, folder_name: str, kpts: int) -> None:
    author_name += str(kpts)
    author_name = author_name.lower().replace(" ", "_")
    current_dir = os.getcwd()
    destination_folder = os.path.join(current_dir, author_name)
    os.makedirs(destination_folder, exist_ok=True)

    new_folder = os.path.join(destination_folder, folder_name)

    os.makedirs(new_folder, exist_ok=True)
    for file in os.listdir(current_dir):
        file_path = os.path.join(current_dir, file)
        if not os.path.isfile(file_path):
            continue
        if ".py" in file or ".sh" in file or ".csv" in file:
            continue
        shutil.move(file_path, os.path.join(new_folder, file))

def calculate(name: str, element: str, defect: str, kpts: int) -> list:
    if not request(element, name):
        return ["", ""]
    result = run(element, defect, kpts)
    move_all_files(name, element + defect, kpts)
    return result

def main():
    result = [["Settings", "KPTS in points per Ang^-1",  "Molecule", "Energy"]]
    element = "Al"
    defect = "Ge"
    for name in AUTHORS_LIST:
        for kpts in range(8, 21, 3):
            row = [name, kpts]
            row.extend(calculate(name, element, defect, kpts))
            row.extend(calculate(name, element, "", kpts))
            row.extend(calculate(name, defect, "", kpts))

            result.append(row)
            result.append([])
            with open('energies.csv', mode='w', newline='') as file: # Save more often so results aren't deleted in the case of a crash
                csv_writer = csv.writer(file)
                csv_writer.writerows(result)

if __name__ == '__main__':
    main()