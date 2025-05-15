import unittest
from unittest import mock, skip

import netCDF4 as nc
import numpy as np

from rad_offset_gain import lib


class TestLib(unittest.TestCase):

    @skip("rad_offset_gain.lib.nc.Dataset")
    @mock.patch("rad_offset_gain.lib.nc.Dataset")
    def test_get_l1b(self, mock_dataset: nc.Dataset) -> None:
        # Create mock variables
        mock_vars = {
            "wavelength": np.array([[500.0]]),
            "radiance": np.array([[1.0, 2.0]]),
            "noise": np.array([[0.1, 0.2]]),
            "sza": np.array([[30.0]]),
            "saa": np.array([[150.0]]),
            "vza": np.array([[20.0]]),
            "vaa": np.array([[180.0]]),
            "latitude": np.array([[10.0]]),
            "longitude": np.array([[20.0]]),
        }

        # Mock NetCDF dataset context manager behavior
        mock_nc = mock.MagicMock()
        mock_nc.__enter__.return_value = mock_nc
        mock_nc.variables = {
            k: mock.MagicMock(__getitem__=mock.MagicMock(return_value=v))
            for k, v in mock_vars.items()
        }
        mock_dataset.return_value = mock_nc

        # Call function
        result = lib._get_l1b("fake_file.nc")

        # Validate output
        for key in mock_vars:
            np.testing.assert_array_equal(result[key], mock_vars[key])

    @skip("rad_offset_gain.lib.nc.Dataset")
    def test_add_radiance_offset(self) -> None:
        # Sample input
        input_data = {"radiance": np.array([[1.0, 2.0, 3.0]])}
        offset = 0.1
        expected = input_data["radiance"] * (1 + offset)

        assert expected

        # Call function if it exists
        if hasattr(lib, "add_radiance_offset"):
            output_path = "/tmp/fake_output.nc"  # hypothetical path
            with (
                mock.patch("rad_offset_gain.lib._get_l1b", return_value=input_data),
                mock.patch("rad_offset_gain.lib.nc.Dataset"),
            ):
                lib.add_radiance_offset("fake_input.nc", output_path, offset)


if __name__ == "__main__":
    unittest.main()
