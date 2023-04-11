"""
IMPRoV 0.5.0


This is the first program of Interactive Proper motions for VLBI (IMProV).

The aim of the program is to interactively identify maser features on single epoch spot maps. You can upload a .csv file with the "upload" button on the top. Note that the CSV file needs the columns "['RA','RAERR','DEC','DECERR','VLSR','FLUX','DFLUX']" and must be comma separated to work. It has a scatter plot of the masers and shows the spectral distribution of the selected masers. You can use the Zoom and Select tools to navigate the data. Use the rectangle tool to select maser spots, and if they show a Gaussian spectral distribution (for data with spectral resolution < 0.5 km/s) you can press the "Save features" button to save a rectangle (x0,x1,y0,y1) containing the specific data points. Once you have selected all the maser features you want to select, you can press the "Export" button to save a csv file containing the rectangles for each feature. Do this for each epoch which you have maser spot maps. The second and third scripts of the program will be used to calculate spot statistics and proper motions.

Author: Job Vorster
Date: April 11, 2023

Usage: python one_feature_identification.py. The script will give an IP which you copy paste into your browser to run the program.

Requirements:
- Python 3, dash, numpy, pandas, matplotlib, scipy, datetime

Usage:
Interaction with the program is via GUI as described above.

Any questions, comments or suggestions can be sent to:
jobvorster8@gmail.com
"""




#########################################################################################
#                                                                                       #
#                               Importing Libraries                                     #
#                                                                                       #
#########################################################################################


from dash import Dash, dcc, html, callback_context, dash_table
import numpy as np
import pandas as pd
from dash.dependencies import Input, Output, State
import plotly.express as px
import scipy.optimize as optim
import plotly.graph_objs as go
from dash.exceptions import PreventUpdate
from matplotlib import cm
from matplotlib import pyplot as plt
import math

import base64
import datetime
import io



#########################################################################################
#                                                                                       #
#                               Global Variables                                        #
#                                                                                       #
#########################################################################################


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
plotted = False
app = Dash(__name__, external_stylesheets=external_stylesheets)
dfcolumns = ['RA min','RA max','Dec min','Dec max']
# make a sample data frame with 6 columns
df = pd.DataFrame(columns = ['RA','RAERR','DEC','DECERR','VLSR','FLUX','DFLUX','Epoch'])



#########################################################################################
#                                                                                       #
#                               GUI Setup                                               #
#                                                                                       #
#########################################################################################


app.layout = html.Div([
    html.Div([html.H1(children='1. Feature Identification', style={'textAlign': 'center'}),
              dcc.Upload(
        id='upload-data',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select Files')
        ]),
        style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        # Allow multiple files to be uploaded
        multiple=True
    ),
        html.Div(id='output-data-upload'),
        dcc.Store(id = 'main-data-store'),
        dcc.Store(id ='current_selected_data_store'),
        dcc.Store(id ='all_selected_data_store'),
        html.Button("GROUP FEATURES", id='btn_group', n_clicks=0),
        html.Div(id ='group_features')],

        ),
    html.Div([

        html.Div(
            dcc.Graph(id='g1', config={'displayModeBar': True, 'scrollZoom': True}, style={'height': '80vh'}),
            className='four columns', style={"border": "2px black solid"}),
            dcc.Store(id='g1-relayout-data'),
        html.Div
            ([ dcc.Store(id='Gauss_Fit_results_data'),
                    html.Button("SAVE FEATURE", id='btn_curvefit', n_clicks=0),
            html.Div(
                dcc.Graph(id='g2', config={'displayModeBar': False}),
                className='row', style={"border": "2px black solid"}),
            html.Div(
                [
                    dcc.Checklist(
                        id='Gauss_Fit_checkbox',
                        options=[{'label': "Do Gaussian Fit", 'value': 'bool_Gaussian'}],
                        value=[]
                    ),
                    dcc.Store(id='Gauss_Fit_checkbox_data'),
                    html.Div
                        ([
                        html.Div
                            ([
                            html.Div(['RA min'], className='two columns'),
                            html.Div(['RA max'], className='two columns'),
                            html.Div(['Dec min'], className='two columns'),
                            html.Div(['Dec max'], className='two columns'), html.Div(['Peak'], className='two columns')
                        ], className='row'),
                        html.Div
                            ([
                            html.Div([' '], id='cont_RA min', className='two columns'),
                            html.Div([' '], id='cont_RA max', className='two columns'),
                            html.Div([' '], id='cont_Dec min', className='two columns'),
                            html.Div([' '], id='cont_Dec max', className='two columns'),
                            html.Div([' '], id='cont_Peak', className='two columns')
                        ], className='row'),
                        html.Div
                            ([
                            html.Div(['Peak err'], className='two columns'),
                            html.Div(['Center'], className='two columns'),
                            html.Div(['Center err'], className='two columns'),
                            html.Div(['FWHM'], className='two columns'), html.Div(['FWHM err'], className='two columns')
                        ], className='row'),
                        html.Div
                            ([
                            html.Div([' '], id='cont_Peak err', className='two columns'),
                            html.Div([' '], id='cont_Center', className='two columns'),
                            html.Div([' '], id='cont_Center err', className='two columns'),
                            html.Div([' '], id='cont_FWHM', className='two columns'),
                            html.Div([' '], id='cont_FWHM err', className='two columns')
                        ], className='row'),
                    ], className='row', style={'height': '13vh'}),
                    # dcc.Store(id='Gauss_Fit_results_data'),
                    # html.Button("SAVE CURRENT GAUSS FIT", id='btn_curvefit', n_clicks=0),
                ], className='row', style={"border": "2px black solid", 'height': '31.2vh'})
        ], className='four columns', style={"border": "2px black solid"}),
        html.Div([
            dash_table.DataTable(id="table_output",columns = [{"name": i,'id':i} for i in dfcolumns],export_format="csv"),
            dcc.Store(id='store_Gauss_Fit_all_results'),
            dcc.Store(id = 'store_Table')
        ], className='four columns', style={"border": "2px black solid", 'overflowY': 'auto'})
    ], className='row', style={"border": "2px black solid"})
])

#########################################################################################
#                                                                                       #
#                               Function Definitions                                    #
#                                                                                       #
#########################################################################################





def parse_contents(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    df1 = []
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df1 = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df1 = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])
    return df1

def monte_carlo_weighted_average(x,dx,weight,dweight,N,weight_method = 'normal'):
    #x = 1/sum(weights)* sum(weight*x)
    weighted_average_arr = []
    for i in range(int(N)):
        if weight_method =='uniform':
            rand_x = np.random.uniform(x-dx,x+dx)
            rand_weight = np.random.uniform(weight-dweight,weight+dweight)
        else:
            rand_x = np.random.normal(x,dx)
            rand_weight = np.random.normal(weight,dweight)
        weighted_average = sum(rand_x*rand_weight)/(sum(rand_weight))
        weighted_average_arr.append(weighted_average)
    return np.mean(weighted_average_arr),np.std(weighted_average_arr)*3

#There are two ways to group the features by VLSR window:
    #First, start at the spot with the lowest VLSR, group all the spots
    #with its VLSR+VLSR window, and just run through all the spots.

#How many monte carlo iterations give a good uncertainty?
    #Which method of the monte carlo method gives the most reliable results?

def group_features_simple(RA,RAERR,DEC,DECERR,VLSR,FLUX,FLUXERR,vlsrrange):
    data_vlsr_range = max(VLSR) - min(VLSR)
    N_window = math.ceil(data_vlsr_range/vlsrrange)
    vlsr_windows = np.arange(min(VLSR),min(VLSR)+(N_window+1)*vlsrrange,vlsrrange)
    RAf = []
    RAf_error = []
    DECf = []
    DECf_error = []
    VLSRf = []
    N = 1000
    for i in range(N_window):

        inds_vlsr = np.where(np.logical_and(VLSR>=vlsr_windows[i],VLSR<=vlsr_windows[i+1]))[0]
        if len(inds_vlsr)!=0:
            app_RA,app_RAerr = monte_carlo_weighted_average(RA[inds_vlsr],RAERR[inds_vlsr],
                                                            FLUX[inds_vlsr],FLUXERR[inds_vlsr],N)
            RAf.append(app_RA)
            RAf_error.append(app_RAerr)
            app_DEC,app_DECerr = monte_carlo_weighted_average(DEC[inds_vlsr],DECERR[inds_vlsr],
                                                            FLUX[inds_vlsr],FLUXERR[inds_vlsr],N)
            DECf.append(app_DEC)
            DECf_error.append(app_DECerr)
            app_VLSR = np.mean([vlsr_windows[i],vlsr_windows[i+1]])
            VLSRf.append(app_VLSR)
    return RAf,RAf_error,DECf,DECf_error,VLSRf
#Secondly, you can take the brightest spot, take +- half the vlsr window,
    #take all the spots, go to next brightest spot, until all the spots are gone. 
# def group_features_brightest(RA,RAERR,DEC,DECERR,VLSR,FLUX,FLUXERR,vlsrrange):

#     return RAf,RAf_error,DECf,DECf_error,VLSRf
def show_scatter(RA,RAERR,DEC,DECERR,VLSR,FLUX,FLUXERR,epoch,selections):
    if selections:
        selection_ids = []
        fitted_lims = []
        for xlims,ylims,ids in selections:
            selection_ids.append(ids)
            fitted_lims.append([xlims,ylims])
    else:
        fitted_lims = []
    dfplot = pd.DataFrame(columns = ['RA','RAERR','DEC','DECERR','VLSR','Epoch'])
    dfplot["RA"]= RA
    dfplot['RAERR'] = RAERR
    dfplot['DEC'] = DEC
    dfplot['DECERR'] = DECERR
    dfplot['VLSR']= VLSR
    dfplot['Epoch']= np.array(epoch,dtype='str')
    cmap1 = plt.cm.ScalarMappable(plt.Normalize(),cmap='jet')
    colors = cmap1.to_rgba(VLSR)
    colors_string = ['rgba(%d,%d,%d,%d)'%(round(r*256),round(g*256),round(b*256),round(a*256)) for r,g,b,a in colors]
    fig = px.scatter(dfplot,x='RA',y='DEC',color = 'VLSR',color_continuous_scale='Jet',text='Epoch',opacity = 0)
    # for i, bar in enumerate(DECERR):
    #     fig.add_trace(go.Scatter(x=[RA[i]],y=[DEC[i]],
    #         error_y=dict( type='data',color = colors_string[i],array=[bar],visible=True),
    #         error_x=dict( type='data',color = colors_string[i],array=[RAERR[i]],visible=True),
    #         marker=dict(color='rgba(0,0,0,0)', size=12),showlegend=False
    #                 ))
        
    fig.update_traces(textfont_size=25,textfont_color = colors_string)
    if fitted_lims:
        for xlim,ylim in fitted_lims:
            fig.add_shape(type='rect',x0=xlim[0],x1=xlim[1],y0=ylim[0],y1=ylim[1],fillcolor='green',opacity=0.1)
    fig.update_shapes(dict(xref='x', yref='y'))
    return fig



#########################################################################################
#                                                                                       #
#                               GUI Callbacks                                           #
#                                                                                       #
#########################################################################################




@app.callback(
    Output('g1','figure'),
    Input('main-data-store','data'),
    Input("all_selected_data_store",'data'))
def update_figure(df,selections):
    if dict(df):
        df = pd.DataFrame.from_dict(df)
        fig = show_scatter(df['RA'],df['RAERR'],df['DEC'],df['DECERR'],
            df['VLSR'],df['FLUX'],df['DFLUX'],df['Epoch'],selections)
        return fig
    else:
        raise PreventUpdate

@app.callback(Output('output-data-upload', 'children'),
              Output('main-data-store','data'),
              Input('upload-data', 'contents'),
              State('upload-data', 'filename'))
def update_output(list_of_contents, list_of_names):
    global df
    if list_of_contents is not None:
        df = pd.DataFrame(columns = ['RA','RAERR','DEC','DECERR','VLSR','FLUX','DFLUX','Epoch'])
        list_of_names,list_of_contents = zip(*sorted(zip(list_of_names, list_of_contents)))
        print(list_of_names)
        for i,(content,name) in enumerate(zip(list_of_contents,list_of_names)):
            df_app = parse_contents(content,name)
            df_app['Epoch'] = i+1
            df = pd.concat([df,df_app])
        return ['Successfully uploaded the files: '+str(list_of_names),dict(df)]
    else:
        return ['Please upload a file.',[]]

#Store selection into the STORE object.
@app.callback(
    Output('current_selected_data_store', 'data'),
    Input('g1', 'selectedData'))
def store_selections(selectedData):
    # print(selectedData)
    if selectedData:
        ids = []
        points = selectedData['points']
        for point in points:
            ids.append(point['pointIndex'])
        selections = [selectedData['range']['x'],selectedData['range']['y'],ids]
        return selections
    else:
        raise PreventUpdate

@app.callback(
    Output('g2','figure'),
    Input('current_selected_data_store', 'data'),
    State('main-data-store','data')
    )
def show_spectral_profile(selectedData,df):
    #selectedData is the xlims and ylims of the data.
    dfplot = pd.DataFrame.from_dict(df)
    ids = selectedData[2]

    fig = px.scatter(dfplot.iloc[ids],x='VLSR',y='FLUX',color = 'VLSR',color_continuous_scale='Jet')
    fig.update_traces(marker=dict(size=12,line=dict(width=2,color='DarkSlateGrey')))
    return fig

@app.callback(
    Output('all_selected_data_store','data'),
    Input('btn_curvefit','n_clicks'),
    State('current_selected_data_store','data'),
    State('all_selected_data_store', 'data')
             )
def store_all_gauss_fits(n_clicks,data_now,data_all):
    changed_id = [p['prop_id'] for p in callback_context.triggered][0]
    if changed_id=='btn_curvefit.n_clicks':
        if not data_all:
            data_all = []
        data_all.append(data_now)
        return data_all
    else:
        raise PreventUpdate

@app.callback(
        Output('table_output','data'),
        Input('all_selected_data_store','data')
              )
def update_table(data):
    if data:
        df_table = pd.DataFrame(columns = dfcolumns)
        RA_min = []
        RA_max= []
        DEC_min = []
        DEC_max = []
        for xlims,ylims,ids in data:
            RA_min.append(xlims[0])
            RA_max.append(xlims[1])
            DEC_min.append(ylims[0])
            DEC_max.append(ylims[1])
        df_table['RA min'] = RA_min
        df_table['RA max'] = RA_max
        df_table['Dec min'] = DEC_min
        df_table['Dec max'] = DEC_max
        return df_table.to_dict('records')
    else:
        raise PreventUpdate
@app.callback(
    Output('group_features','children'),
    Input('btn_group','n_clicks'),
    State('all_selected_data_store','data')
    )
def group_all_features(n_clicks,limits):
    changed_id = [p['prop_id'] for p in callback_context.triggered][0]
    if changed_id=='btn_group.n_clicks':
        if limits:
            vlsrrange = 0.1
            RAg = []
            RAg_error = []
            DECg = []
            DECg_error = []
            VLSRg = []
            for xlim,ylim,inds in limits:
                dfarr = [df['RA'],df['RAERR'],
                    df['DEC'],df['DECERR'],df['VLSR'],df['FLUX'],df['DFLUX']]
                
                ra1,raerr1,dec1,decerr1,vlsr1,flux1,fluxerr1 = [arr1.values[inds] for arr1 in dfarr]
                RAf,RAf_error,DECf,DECf_error,VLSRf = group_features_simple(ra1,raerr1,dec1,
                    decerr1,vlsr1,flux1,fluxerr1 ,vlsrrange)
                RAg.append(RAf)
                RAg_error.append(RAf_error)
                DECg.append(DECf)
                DECg_error.append(DECf_error)
                VLSRg.append(VLSRf)
            # print('RA')
            # print(RAg)
            # print('DRA')
            # print(RAg_error)
            # print('DEC')
            # print(DECg)
            # print('DDEC')
            # print(DECg_error)    
            # print('VLSR')
            # print(VLSRg)
            # return 
        else:
            raise PreventUpdate
    else:
        raise PreventUpdate

if __name__ == '__main__':
    app.run_server(debug=True)
