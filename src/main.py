import os
import shutil

from ase import Atoms
from ase.build import bulk
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS

from get_input_files import request

#AUTHORS_LIST = ["Stefano Curtarolo", "Chris Wolverton", "Materials Project", "Miguel Marques", "Bjoern Bieniek", "Oliver Hoffmann", "Kurt Lejaeghere"]
AUTHORS_LIST = ["Stefano Curtarolo"]

def run(element: str, defect: str):
    fcc = bulk(element, crystalstructure="fcc", a=4)
    super_cell: Atoms = fcc.repeat((2, 2, 2))
    if defect != "":
        super_cell[0].symbol = defect
    super_cell.calc = Vasp(xc="PBE")

    optimizer = BFGS(super_cell)
    optimizer.run(fmax=0.02)

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
        if ".py" not in file and os.path.isfile(file_path):
            shutil.move(file_path, os.path.join(new_folder, file))


def main():
    element = "Al"
    defect = "Ge"
    for name in AUTHORS_LIST:
        if request(element, name):
            run(name, "")
            move_all_files(name, False)
            run(name, defect)
            move_all_files(name, True)

if __name__ == '__main__':
    main()

