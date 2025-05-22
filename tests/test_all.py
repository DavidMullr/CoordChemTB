import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from CoordChenTB.utils.complexorbitalssplitting import get_d_electron_count, get_pairing_energy
from CoordChenTB.utils.metals_db import METALS

def test_get_d_electron_count():
    assert get_d_electron_count(26, 2) == 6
    assert get_d_electron_count(25, 3) == 4

def test_get_pairing_energy():
    assert get_pairing_energy("3d") == 17000
    assert get_pairing_energy("4d") == 13000
    assert get_pairing_energy("5d") == 12000
    assert get_pairing_energy("f") == 11000
    # Default fallback
    assert get_pairing_energy("unknown") == 15000
