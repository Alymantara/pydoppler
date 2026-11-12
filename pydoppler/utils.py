"""Utility helpers for PyDoppler native."""
import os
import shutil

def test_data():
    """Copy U Gem test dataset into the current directory."""
    src = os.path.join(os.path.dirname(__file__), 'test_data')
    if not os.path.isdir(src):
        raise FileNotFoundError(f"Test data not found at {src}")
    for item in os.listdir(src):
        if item == '__pycache__':     # created by pip when installing
            continue
        s = os.path.join(src, item)
        d = os.path.join('.', item)
        if os.path.isdir(s):
            if os.path.exists(d):
                shutil.rmtree(d)
            shutil.copytree(s, d)
        else:
            shutil.copy2(s, d)
    print("  ✓  Test data copied to current directory.")
