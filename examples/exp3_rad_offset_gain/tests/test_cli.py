import sys
import unittest
from unittest import mock

from rad_offset_gain import cli


class TestCLI(unittest.TestCase):

    @mock.patch("rad_offset_gain.cli.add_radiance_offset")
    def test_cli_success(self, mock_add):
        test_args = [
            "prog",
            "--input",
            "input.nc",
            "--output",
            "output.nc",
            "--offset",
            "0.05",
        ]
        with (
            mock.patch.object(sys, "argv", test_args),
            mock.patch("os.path.exists", return_value=True),
        ):
            cli.main()
            mock_add.assert_called_once_with("input.nc", "output.nc", 0.05)

    def test_missing_arguments(self):
        test_args = ["prog", "--input", "input.nc"]
        with mock.patch.object(sys, "argv", test_args), self.assertRaises(SystemExit):
            cli.main()

    def test_input_file_not_found(self):
        test_args = ["prog", "--input", "missing.nc", "--output", "output.nc"]
        with (
            mock.patch.object(sys, "argv", test_args),
            mock.patch("os.path.exists", return_value=False),
            self.assertRaises(FileNotFoundError),
        ):
            cli.main()

    @mock.patch("rad_offset_gain.cli.add_radiance_offset", side_effect=Exception("Mock error"))
    def test_processing_error(self, mock_add):
        test_args = ["prog", "--input", "input.nc", "--output", "output.nc"]
        with (
            mock.patch.object(sys, "argv", test_args),
            mock.patch("os.path.exists", return_value=True),
            self.assertRaises(RuntimeError),
        ):
            cli.main()


if __name__ == "__main__":
    unittest.main()
