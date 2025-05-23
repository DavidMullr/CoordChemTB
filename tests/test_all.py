import pytest
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from CoordChenTB.utils.metals_db import METALS
from CoordChenTB.utils.complexorbitalssplitting import get_d_electron_count, get_pairing_energy

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

from CoordChenTB.utils.complexorbitalssplitting import get_metal_data

def test_get_metal_data_valid():
    fe = get_metal_data("Fe")
    assert fe is not None
    assert fe["symbol"] == "Fe"
    assert "oxidation_states" in fe

def test_get_metal_data_invalid():
    assert get_metal_data("Xx") is None

from CoordChenTB.utils.complexorbitalssplitting import predict_spin_state

def test_predict_spin_state_low_spin():
    assert predict_spin_state(20000, 15000) == "Low-spin"

def test_predict_spin_state_high_spin():
    assert predict_spin_state(10000, 15000) == "High-spin"

from CoordChenTB.utils.complexorbitalssplitting import load_ligand_data

def test_load_ligand_data_structure():
    ligand_data = load_ligand_data("all_ligands.sdf")
    assert isinstance(ligand_data, dict)
    if ligand_data:
        sample = next(iter(ligand_data.values()))
        assert "Delta" in sample
        assert "Charge" in sample

from CoordChenTB.utils.complexorbitalssplitting import find_closest_ligand

def test_find_closest_ligand_exact_match(monkeypatch):
    data = {"H2O": {"Delta": 10000, "Charge": 0}}

    # Simulate user entering "1" to accept first suggestion
    monkeypatch.setattr("builtins.input", lambda _: "1")

    match = find_closest_ligand("H2O", data)
    assert match == "H2O"

def test_find_closest_ligand_close_match(monkeypatch):
    data = {"ammonia": {"Delta": 15000, "Charge": 0}}

    monkeypatch.setattr("builtins.input", lambda _: "1")
    match = find_closest_ligand("ammnia", data)
    assert match == "ammonia"
