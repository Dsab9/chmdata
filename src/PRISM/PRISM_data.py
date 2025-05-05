import pyPRISMClimate as ppc
import os


class PRISM:
    def __init__(self, start_year, end_year, variable, root_directory):
        """
        PRISM (Parameter-elevation Regressions on Independent Slopes Model) data downloader
        Utilizes the pyPRISMClimate package
        Available variables include: 'tmean', 'tmin', 'tmax', 'ppt', 'vpdmon', 'vpdmax'
        ppt & tmean required for Seasonality index and %precip as snow calc

        -monthly_data(): gets averaged monthly 4km rasterized data for given time interval

        """

        self.years = list(range(start_year, end_year + 1))
        self.months = list(range(1, 12 + 1))
        self.variable = variable
        self.root_directory = root_directory
        # Ensure the output directory exists
        if not os.path.exists(root_directory):
            os.makedirs(root_directory)

    def monthly_data(self):
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

