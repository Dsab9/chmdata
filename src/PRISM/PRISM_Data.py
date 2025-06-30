import os

import pyPRISMClimate as ppc


class PRISM:
    """
        PRISM data downloader

         PRISM (Parameter-elevation Regressions on Independent Slopes Model) data downloader
         Utilizes the pyPRISMClimate package.

         Note: variables ppt & tmean required for Seasonality index and %precip as snow calc

         Args:
             start_year (int): daily data 1981 to present
             end_year (int): see above
             variable (str): options include 'tmean', 'tmin', 'tmax', 'ppt', 'vpdmon', 'vpdmax'
             root_directory (str): absolute filepath of directory where data is downloaded to
    """

    def __init__(self, start_year: int, end_year: int, variable: str, root_directory: str):
        self.years = list(range(start_year, end_year + 1))
        self.months = list(range(1, 12 + 1))
        self.variable = variable
        self.root_directory = root_directory
        # Ensure the output directory exists
        if not os.path.exists(root_directory):
            os.makedirs(root_directory)
        self.monthly_data()

    def monthly_data(self):
        """ gets averaged monthly 4km rasterized data for given time interval, saves in 'root_directory' """

        print(f"Downloading monthly {self.variable} data ...")
        v_directory = os.path.join(self.root_directory, self.variable)
        if not os.path.exists(v_directory):
            os.makedirs(v_directory)
        try:
            ppc.get_prism_monthlys(
                variable=self.variable,
                years=self.years,
                months=self.months,
                dest_path=v_directory,
                # keep_zip=False  # Crashing when set to True, need to delete zip folders manually
            )
            print("Download complete!")
            print(f"Files saved in: {v_directory}")

        except Exception as e:
            print(f"An error occurred: {e}")
            print(
                "Please ensure you have the pyPRISMClimate library installed and that your internet connection is stable.")
        print("Script finished.")
