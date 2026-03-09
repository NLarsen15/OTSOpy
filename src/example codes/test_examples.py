"""
Pytest-based tests for OTSO example codes.

This module contains tests that run each example script and verify
they execute without errors.
"""

import pytest
import subprocess
import sys
from pathlib import Path


# Get the directory containing the example files
EXAMPLES_DIR = Path(__file__).parent

# Add project root to Python path so we can import OTSO without installation
PROJECT_ROOT = EXAMPLES_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def get_example_files():
    """Get all Python example files except test files."""
    return [f for f in EXAMPLES_DIR.glob("*.py") 
            if not f.name.startswith("test_") and f.name != "conftest.py"]


class TestOTSOExamples:
    """Test class for OTSO example scripts."""
    
    @pytest.mark.parametrize("example_file", get_example_files(), ids=lambda x: x.stem)
    def test_example_runs_successfully(self, example_file):
        """Test that an example script runs without errors and produces valid output."""
        # Set up environment to include project root in Python path
        env = {**subprocess.os.environ, 'PYTHONPATH': str(PROJECT_ROOT)}
        
        result = subprocess.run(
            [sys.executable, example_file.name],
            cwd=EXAMPLES_DIR,
            capture_output=True,
            text=True,
            timeout=180,  # 3 minute timeout
            env=env
        )
        
        # Check that the script ran successfully
        if result.returncode != 0:
            pytest.fail(
                f"Example {example_file.name} failed with return code {result.returncode}\n"
                f"STDOUT: {result.stdout}\n"
                f"STDERR: {result.stderr}"
            )
        
        # Validate that there's actual output (not just empty)
        if not result.stdout.strip():
            pytest.fail(f"Example {example_file.name} produced no output")
        
        # Check for common error patterns in output
        stdout_lower = result.stdout.lower()
        error_patterns = [
            'traceback', 'exception occurred', 'importerror', 
            'modulenotfounderror', 'attributeerror'
        ]
        
        for pattern in error_patterns:
            if pattern in stdout_lower:
                pytest.fail(
                    f"Example {example_file.name} output contains '{pattern}': {result.stdout[:200]}"
                )
        
        # Validate specific outputs based on example type
        self._validate_example_output(example_file.stem, result.stdout)
    
    def _validate_example_output(self, example_name, output):
        """Validate specific outputs based on the example type."""
        output_lines = output.strip().split('\n')
        
        if example_name == 'trajectory':
            # trajectory should return dictionary and text info
            assert len(output_lines) >= 2, f"trajectory output too short: {len(output_lines)} lines"
            # Look for dictionary-like output or structured data
            has_structure = any('{' in line or 'OULU' in line or 'ROME' in line for line in output_lines)
            assert has_structure, "trajectory should output dictionary or station data"
        
        elif example_name == 'cone':
            # cone should return dataframes and text info  
            assert len(output_lines) >= 2, f"cone output too short: {len(output_lines)} lines"
            # Look for meaningful content
            assert any(line.strip() for line in output_lines), "cone should produce non-empty output"
        
        elif example_name == 'cutoff':
            # cutoff should return dataframe and text info
            assert len(output_lines) >= 2, f"cutoff output too short: {len(output_lines)} lines"
        
        elif example_name == 'magfield':
            # magfield should return magnetic field data table and info
            assert len(output_lines) >= 5, f"magfield output too short: {len(output_lines)} lines"
            # Check for coordinate/magnetic field column headers or data
            has_field_data = any('GSM_B' in line or 'GEO' in line or 'nT' in line for line in output_lines)
            assert has_field_data, "magfield should output magnetic field data with coordinates"
        
        elif example_name == 'planet':
            # planet computations should have output
            assert len(output_lines) >= 1, "planet should produce some output"
        
        elif example_name == 'coordtrans':
            # coordinate transformations should show numeric results
            has_numbers = any(any(char.isdigit() for char in line) for line in output_lines)
            assert has_numbers, "coordtrans should output numeric coordinate values"
        
        elif example_name == 'flight':
            # flight trajectory should have output
            assert len(output_lines) >= 1, "flight should produce some output"
        
        elif example_name == 'maglinetrace':
            # magnetic line tracing should have output
            assert len(output_lines) >= 1, "maglinetrace should produce some output"
    
    def test_otso_imports(self):
        """Test that OTSO package imports work correctly."""
        try:
            import OTSO
            from OTSO import trajectory, cone, cutoff, magfield, planet, coordtrans, flight
        except ImportError as e:
            pytest.fail(f"Failed to import OTSO package: {e}")