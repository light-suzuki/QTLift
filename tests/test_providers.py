import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from qtlift.providers import validate_provider


class ProviderValidationTests(unittest.TestCase):
    def test_configured_native_executable_avoids_false_unavailable_warning(self):
        with patch("qtlift.providers.shutil.which", return_value=None):
            warnings = validate_provider("windows", native_blastn_path="/custom/blastn")
        self.assertEqual(warnings, [])

    def test_missing_native_executable_still_warns(self):
        with patch("qtlift.providers.shutil.which", return_value=None):
            warnings = validate_provider("windows")
        self.assertTrue(any("unavailable" in warning for warning in warnings))

    def test_native_executable_does_not_suppress_missing_wsl_warning(self):
        with patch("qtlift.providers.shutil.which", return_value=None):
            warnings = validate_provider("wsl", native_blastn_path="/custom/blastn")
        self.assertTrue(any("wsl.exe is unavailable" in warning for warning in warnings))


if __name__ == "__main__":
    unittest.main()
