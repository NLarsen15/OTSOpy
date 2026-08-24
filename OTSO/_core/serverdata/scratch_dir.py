import os
import tempfile

SCRATCH_DIR = os.path.join(tempfile.gettempdir(), "OTSO_serverdata_scratch")
os.makedirs(SCRATCH_DIR, exist_ok=True)
