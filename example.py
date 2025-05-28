from src.chmdata import Agrimet, Mesonet, met_utils, thredds
from src.PRISM import PRISM_data, Percent_Precip_as_Snow, Seasonality_Index


prism = PRISM_data.PRISM(2022, 2023, 'tmean', 'Downloads')
prism.monthly_data()