import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from io import StringIO

file_path = 'Persistance Hot Fire 1'

# Reading the file content first to handle initial metadata lines
with open(file_path, 'r') as file:
    file_content = file.readlines()

# Combine the file content into a single string for reading into pandas
file_content_str = ''.join(file_content)

# Read the data into a DataFrame, skipping the first 3 metadata lines
data = pd.read_csv(StringIO(file_content_str), skiprows=3)

# Display the first few rows of the DataFrame to verify
print(data.head())

# Smoothing the data for TCNitrousSupply and TCNitrousRun
data['TCNitrousSupply_Smooth'] = data['TCNitrousSupply'].rolling(window=50).mean()
data['TCNitrousRun_Smooth'] = data['TCNitrousRun'].rolling(window=50).mean()

# Define the list of columns to plot
columns_to_plot = ['TCNitrousSupply', 'TCNitrousSupply_Smooth', 'TCNitrousRun', 'TCNitrousRun_Smooth', 
                   'PTN2OSupply', 'PTRun', 'PTPreInjector', 'PTEngine', 'PTN2Supply', 'LCNitrousFill', 'LCThrust']

# Create a subplot figure with shared x-axis
fig = make_subplots(rows=len(columns_to_plot), cols=1, shared_xaxes=True, 
                    subplot_titles=columns_to_plot, vertical_spacing=0.03)

# Add traces for each column
for i, column in enumerate(columns_to_plot, 1):
    fig.add_trace(go.Scatter(x=data['Time'], y=data[column], mode='lines', name=column), row=i, col=1)

# Update layout
fig.update_layout(height=3000, width=1000, title_text='Interactive Plots of Rocket Engine Test Data', showlegend=False)

# Save plot to an HTML file and open in a web browser
fig.write_html("rocket_engine_test_data.html")

# Display the time windows where LCThrust spikes
print("Time windows where LCThrust spikes:")
print(data[data['LCThrust'] > (data['LCThrust'].mean() + 2 * data['LCThrust'].std())]['Time'])
