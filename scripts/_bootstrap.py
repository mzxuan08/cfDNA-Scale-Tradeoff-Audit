"""Direct script execution uses the repository package without installation."""
import argparse
import os
from pathlib import Path
import sys

os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def parser(description):
    p = argparse.ArgumentParser(description=description)
    p.add_argument('--config', default=None, help='Config YAML; default config/config.yaml')
    return p
