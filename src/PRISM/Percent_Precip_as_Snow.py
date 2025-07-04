import os

import matplotlib.pyplot as plt
import numpy as np
import rasterio as rio
from rasterio.plot import show


class PPS:
    """
        Calculates percent precip as snow from downloaded monthly PRISM data.

        Calculation for percent precipitation as snow from McGabe and Wolock 2009:
        https://journals.ametsoc.org/view/journals/eint/13/12/2009ei283.1.xml

        Percent snowfall for values between Train and Tsnow follow a linear gradient,
        default values are provided with the ability to change if desired.

        Note: Seasonality index module must be run prior to running PPS to aggregate precip data

        Args:
            start_year (int): daily data 1981 to present
            end_year (int): see above
            root_directory (str): absolute filepath of root directory where data is downloaded to
            train (int): Temp, in Celsius, at which % precip is 0% snow (default value is 3)
            tsnow (int): Temp, in Celsius, at which % precip is 100% snow (default value is -1)
    """

    def __init__(self, start_year: int, end_year: int, root_directory: str, train: int = 3, tsnow: int = -1):
        self.years = list(range(start_year, end_year + 1))
        self.months = list(range(1, 12 + 1))
        self.root_directory = root_directory
        self.train = train
        self.tsnow = tsnow
        # Ensure directories exist
        self.output_directory = os.path.join(root_directory, 'PPS')
        self.precip_directory = os.path.join(root_directory, 'ppt')
        self.tmean_directory = os.path.join(root_directory, 'tmean')
        if not os.path.exists(self.output_directory):
            os.makedirs(root_directory)
        if not os.path.exists(self.precip_directory):
            print("No ppt directory")
        if not os.path.exists(self.tmean_directory):
            print("No tmean directory")

    def create_file_lists(self) -> list:
        """Compiles list of daily averages (mean temp) from directory"""
        file_list = []
        for m in self.months:
            month_list = []
            for y in self.years:
                file_name = f'PRISM_tmean_stable_4kmM3_{str(y)}{str(m)}_bil.bil'
                month_list.append(file_name)
            file_list.append(month_list)
        print(file_list)
        return file_list

    def get_monthly_averages(self, raster_file_lists: list):
        """Uses output from def create_file_lists to generate monthly average temp rasters."""
        for month_list in raster_file_lists:
            raster_arrays = []
            for raster_file in month_list:
                path = os.path.join(self.tmean_directory, raster_file)
                with rio.open(path) as src:
                    raster_arrays.append(src.read(1).astype(np.float32))  # Read the first band, convert to float
                    profile = src.profile  # Store metadata for output
                summed_array = np.nansum(raster_arrays,
                                         axis=0)  # axis zero sums across arrays, axis 1 sums within array
                ave_array = summed_array / len(raster_arrays)
            print(f"Doing Maths for Month {self.months[raster_file_lists.index(month_list)]}...")

            filename = f'tmean_month{self.months[raster_file_lists.index(month_list)]}.tif'
            out_path = os.path.join(self.output_directory, filename)
            profile.update(dtype=rio.float32, compress='lzw')
            with rio.open(out_path, 'w', **profile) as dst:
                dst.write(ave_array, 1)

    def calc_percent_snow(self):
        """Calculates pps from tmean and ppt raster data, outputs data as raster."""
        month_temp_arrays = []
        for m in self.months:
            file = f'tmean_month{m}.tif'
            path = os.path.join(self.output_directory, file)
            with rio.open(path) as src:
                month_temp_arrays.append(src.read(1).astype(np.float32))

        month_ppt_arrays = []
        for m in self.months:
            file = f'avePPT_month{m}.tif'
            path = os.path.join(self.precip_directory, file)
            with rio.open(path) as src:
                month_ppt_arrays.append(src.read(1).astype(np.float32))

        file = 'annual_ppt.tif'
        path = os.path.join(self.precip_directory, file)
        with rio.open(path) as src:
            annual_ppt = src.read(1).astype(np.float32)
            profile = src.profile

        snow_ratio_arrays = []
        for t in month_temp_arrays:
            ratio = np.copy(t).astype(float)
            ratio = np.where((t > self.tsnow) & (t < self.train), ((self.train - t) / (self.train - self.tsnow)), ratio)
            ratio = np.where(t >= 3, 0, ratio)
            ratio = np.where(t <= -1, 1, ratio)
            snow_ratio_arrays.append(ratio)

        precip_as_snow_arrays = [a * b for a, b in zip(snow_ratio_arrays, month_ppt_arrays)]
        total_snow = np.nansum(precip_as_snow_arrays, axis=0)
        percent_snow = (total_snow / annual_ppt) * 100

        outpath = os.path.join(self.output_directory, 'percent_snow.tif')
        with rio.open(outpath, 'w', **profile) as dst:
            dst.write(percent_snow, 1)

        return percent_snow

    def calc_PPS(self):
        """""Calculates percent precip as snow skipping intermediary steps, outputs data as raster."""

        raster_file_lists = self.create_file_lists()
        self.get_monthly_averages(raster_file_lists)
        self.calc_percent_snow()

        # # Open raster and plot
        raster_path = os.path.join(self.output_directory, 'percent_snow.tif')
        with rio.open(raster_path) as src:
            raster_data = src.read(1)  # Read the first band
            raster_meta = src.meta
            omit_value = 100
            masked_data = np.where(raster_data == omit_value, np.nan, raster_data)
        fig, ax = plt.subplots()
        show(masked_data, transform=raster_meta['transform'], ax=ax, cmap='viridis')
        ax.set_title('% Precip as Snow')

        plt.show()
