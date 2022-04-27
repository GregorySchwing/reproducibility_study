"""Create atomistic representation of ethanol."""
import os

import mbuild as mb

from reproducibility_project.src import molecules


class opc3_rigid(mb.Compound):
    """Create a single particle water compound."""

    def __init__(self):
        super(opc3_rigid, self).__init__()
        abs_path = os.path.dirname(os.path.abspath(molecules.__file__))
        self.add(mb.load(f"{abs_path}/opc3_rigid.mol2"), label="WAT")


def main():
    """Create a opc3 compound and print basic properties."""
    water = opc3_rigid()
    print(water)
    print(water.name)
    print(water.labels)
    print(water["WAT"])


if __name__ == "__main__":
    main()
