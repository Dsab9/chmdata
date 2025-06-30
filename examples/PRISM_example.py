from src.PRISM import PRISM_Data, Percent_Precip_as_Snow, Seasonality_Index

# get mean temp and pecip monthly averaged data, will create relative file path named "Downloads"
PRISM_Data.PRISM(2013, 2023, 'tmean', 'Downloads')
PRISM_Data.PRISM(2013, 2023, 'ppt', 'Downloads')

# First we need to calculate precip seasonality index, this will aggregate and save our precip data
SI = Seasonality_Index.SI(2013, 2023, 'Downloads')
SI.calc_SI()

# We can now calculate percent annual precip fallen as snow
PPS = Percent_Precip_as_Snow.PPS(2013, 2023, 'Downloads')
PPS.calc_PPS()
