from netCDF4 import Dataset

dataset = Dataset("wghc_params.nc")

x = dataset['LON']
print(x)
