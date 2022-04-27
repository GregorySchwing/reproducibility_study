"""Create atomistic representation of ethanol."""
import os

import mbuild as mb

from reproducibility_project.src import molecules


class tip4p_2005_rigid(mb.Compound):
    """Create a single particle water compound."""

    def __init__(self):
        super(tip4p_2005_rigid, self).__init__()
        abs_path = os.path.dirname(os.path.abspath(molecules.__file__))
        self.add(mb.load(f"{abs_path}/tip4p_2005_rigid.mol2"), label="WAT")


def main():
    """Create a tip4p_2005 compound and print basic properties."""
    water = tip4p_2005_rigid()
    print(water)
    print(water.name)
    print(water.labels)
    print(water["WAT"])


if __name__ == "__main__":
    main()