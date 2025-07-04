import os

import matplotlib.pyplot as plt
import numpy as np
import rasterio as rio
from rasterio.plot import show


class SI:
    """
    Calculates Seasonality Index from downloaded monthly PRISM data.

    Calculation for precipitation seasonality index from Imteaz & Hossain 2022:
    https://link.springer.com/article/10.1007/s11269-022-03320-z

    Note: Requires monthly ppt rasters in root directory (see PRISM_data module)

    Args:
        start_year (int): daily data 1981 to present
        end_year (int): see above
        root_directory (str): absolute filepath of root directory where data is downloaded to
    """

    def __init__(self, start_year: int, end_year: int, root_directory: str):
        self.years = list(range(start_year, end_year + 1))
        self.months = list(range(1, 12 + 1))
        self.root_directory = root_directory
        # Ensure directories exist
        self.output_directory = os.path.join(root_directory, 'SI')
        self.precip_directory = os.path.join(root_directory, 'ppt')
        if not os.path.exists(self.output_directory):
            os.makedirs(root_directory)
        if not os.path.exists(self.precip_directory):
            print("No ppt directory")

    def create_file_lists(self) -> list:
        """Compiles list of daily averages from directory"""
        file_list = []
        for m in self.months:
            month_list = []
            for y in self.years:
                file_name = f'PRISM_ppt_stable_4kmM3_{str(y)}{str(m)}_bil.bil'
                month_list.append(file_name)
            file_list.append(month_list)
        print(file_list)
        return file_list

    def get_monthly_averages(self, raster_file_lists: list):
        """Uses output from def create_file_lists to generate monthly average ppt rasters."""
        for month_list in raster_file_lists:
            raster_arrays = []
            for raster_file in month_list:
                path = os.path.join(self.precip_directory, raster_file)
                with rio.open(path) as src:
                    raster_arrays.append(src.read(1).astype(np.float32))  # Read the first band, convert to float
                    profile = src.profile  # Store metadata for output
                summed_array = np.nansum(raster_arrays,
                                         axis=0)  # axis zero sums across arrays, axis 1 sums within array
                ave_array = summed_array / len(raster_arrays)
            print(f"Doing Maths for Month {self.months[raster_file_lists.index(month_list)]}...")

            filename = f'avePPT_month{self.months[raster_file_lists.index(month_list)]}.tif'
            out_path = os.path.join(self.output_directory, filename)
            profile.update(dtype=rio.float32, compress='lzw')
            with rio.open(out_path, 'w', **profile) as dst:
                dst.write(ave_array, 1)

    def get_annual_average(self):
        """Averages monthly ave ppt raster to produce average annual ppt raster."""
        raster_arrays = []
        for m in self.months:
            file = f'avePPT_month{m}.tif'
            path = os.path.join(self.output_directory, file)
            with rio.open(path) as src:
                raster_arrays.append(src.read(1).astype(np.float32))
                profile = src.profile
        summed_array = np.nansum(raster_arrays, axis=0)

        profile.update(dtype=rio.float32, compress='lzw')
        outpath = os.path.join(self.output_directory, 'annual_ppt.tif')
        with rio.open(outpath, 'w', **profile) as dst:
            dst.write(summed_array, 1)

        ave_array = summed_array / len(raster_arrays)
        profile.update(dtype=rio.float32, compress='lzw')
        outpath = os.path.join(self.output_directory, 'ave_month_ppt.tif')
        with rio.open(outpath, 'w', **profile) as dst:
            dst.write(ave_array, 1)

    def calc_seasonality_index(self):
        """Calculates SI from monthly average and annual average rasters, outputs data as raster."""
        month_arrays = []
        for m in self.months:
            file = f'avePPT_month{m}.tif'
            path = os.path.join(self.output_directory, file)
            with rio.open(path) as src:
                month_arrays.append(src.read(1).astype(np.float32))

        file = 'annual_ppt.tif'
        path = os.path.join(self.output_directory, file)
        with rio.open(path) as src:
            annual_array = src.read(1).astype(np.float32)
            profile = src.profile

        file = 'ave_month_ppt.tif'
        path = os.path.join(self.output_directory, file)
        with rio.open(path) as src:
            ave_month_array = src.read(1).astype(np.float32)

        abs_dev_arrays = []
        for a in month_arrays:
            dif = np.subtract(a, ave_month_array)
            abs_ = np.absolute(dif)
            abs_dev_arrays.append(abs_)
        sum_abs_dev_array = np.nansum(abs_dev_arrays, axis=0)
        si = sum_abs_dev_array / annual_array

        outpath = os.path.join(self.output_directory, 'SI.tif')
        with rio.open(outpath, 'w', **profile) as dst:
            dst.write(si, 1)

        return si

    def calc_SI(self):
        """Calculates precip seasonality index skipping intermediary steps, outputs data as raster."""
        raster_file_lists = self.create_file_lists()
        self.get_monthly_averages(raster_file_lists)
        self.get_annual_average()
        self.calc_seasonality_index()

        # # Open raster and plot
        raster_path = os.path.join(self.output_directory, 'SI.tif')
        with rio.open(raster_path) as src:
            raster_data = src.read(1)  # Read the first band
            raster_meta = src.meta
            omit_value = 0
            masked_data = np.where(raster_data == omit_value, np.nan, raster_data)
        fig, ax = plt.subplots()
        show(masked_data, transform=raster_meta['transform'], ax=ax, cmap='viridis')  # You can change the cmap
        ax.set_title('Seasonality Index')
        plt.show()
