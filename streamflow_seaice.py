import pandas as pd
import numpy as np
import scipy.signal as ss

# verify New Castle, VA

extent = pd.read_excel('Sea_Ice_Index_Daily_Extent_G02135_v3.0.xlsx', sheet_name=1)
extent = extent.T
extent = extent.drop(['Unnamed: 0', 'Unnamed: 1', 1978, 1979, 2021, 2022, '1981-2010 mean', '1981-2010 median'])
extent = extent.astype('float')
extent = extent.interpolate(method='linear', axis=0)
extent = extent.values.flatten()
seaice_detrended = ss.detrend(extent)
seaice = np.array(seaice_detrended).reshape(41, 366)
seaice = pd.DataFrame(seaice)

# get sea ice data for DJF
seaice_jf = seaice.loc[:, 0:59]
seaice_d = seaice.loc[:, 335:]
seaice_djf = pd.concat([seaice_d, seaice_jf], axis=1)
seaice_mean_djf = seaice_djf.mean(axis=1).to_list()

# get sea ice data for JJA
seaice_jja = seaice.loc[:, 152:243]
seaice_jja = pd.DataFrame(seaice_jja)
seaice_mean_jja = seaice_jja.mean(axis=1).to_list()

oni = pd.read_excel('oni.xlsx')
oni_djf = oni['DJF'].to_list()
oni_jja = oni['JJA'].to_list()

num_stations = 15
seaice_corr = []
oni_corr = []

for sheet in range(num_stations):
    streamflow_df = pd.read_excel('Streamflow Data.xlsx', sheet_name=sheet)
    # extract column with streamflow data
    streamflow_list = streamflow_df['ft3/s'].to_list()
    print(streamflow_list)
    # convert data from ft^3/s to m^3/s
    streamflow_list = [value / 35.3147 for value in streamflow_list]
    # compute mean and SD
    mean = np.mean(streamflow_list)
    sd = np.std(streamflow_list)
    # standardize the data
    standardized_values = []
    for item in range(len(streamflow_list)):
        standardized_values.append(((streamflow_list[item] - mean) / sd))

    june_streamflow = standardized_values[5::12]
    july_streamflow = standardized_values[6::12]
    august_streamflow = standardized_values[7::12]

    sf_df = pd.DataFrame([june_streamflow, july_streamflow, august_streamflow])
    sf_df = sf_df.T
    sf_df = sf_df.mean(axis=1).to_list()

    r1 = np.corrcoef(seaice_mean_jja, sf_df)
    r2 = np.corrcoef(oni_djf, sf_df)
    seaice_corr.append(r1[0, 1])
    oni_corr.append(r2[0,1])

avg_corr = np.mean(seaice_corr)
sd_corr = np.std(seaice_corr)

oni_avg = np.mean(oni_corr)

print(avg_corr)
print(sd_corr)

print(oni_avg)
