from dash import Dash, dcc, html, callback_context, dash_table
import numpy as np
import pandas as pd
from dash.dependencies import Input, Output, State
import plotly.express as px
import scipy.optimize as optim
import plotly.figure_factory as ff
import plotly.graph_objs as go
from dash.exceptions import PreventUpdate
from matplotlib import cm
from matplotlib import pyplot as plt
import matplotlib
import base64
import datetime
import io

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
plotted = False
app = Dash(__name__, external_stylesheets=external_stylesheets)
df = []
cols = ['ID','L value','Epochs','RA','RAerr','Dec','Decerr','mu_alpha','dmu_alpha','mu_delta','dmu_delta','VLSR','EpID']
all_dates = []
#Epdate_float: G9.62 [2011+330/365.25,2012+256/365.25,2013+263/365.25] Source_dec: -20-31/60.-35/3600
#Epdate_float: NGC6334I  [2014.717808219178, 2014.9013698630138, 2015.0849315068492, 2015.2849315068493, 2015.8821917808218, 2016.109589041096, 2016.194520547945]
epdate_float = [2008 + 327./365., 2009 + 32./365., 2009 + 138./365, 2009 + 138./365., 2009 + 240./365.,
                2009 + 258./365., 2009 + 270./365., 2009 + 297./365., 2009 + 347./365., 2010 + 41./365., 2010 + 94./365., 2010 + 233./365. ]
#NGC6334I source_Dec= -35.78375
source_Dec = 17 + 59/60 + 22.95890/3600
def straight_line(x,a,b):
    return a*x+b
def calc_mu(x,dx,t):
    m = (x[1]-x[0])/(t[1]-t[0])
    p0 = np.array([m,x[0]+m*t[0]],dtype='float64')
    sigma = np.array(dx,dtype='float64')/np.array(x,dtype='float64')
    popt,pcov = optim.curve_fit(straight_line,xdata=t,ydata=x,sigma=sigma,p0=p0)
    return popt[0],np.sqrt(np.diag(pcov))[0]


def parse_contents(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    df1 = 'Error'
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

def calculate_likeness(df,inds,a=1,b=3,c=2):
    peaks = df['Peak'].values[inds]
    epochs = df['Epoch'].values[inds]
    fwhms = df['FWHM'].values[inds]
    vlsrs = df['VLSR'].values[inds]
    L = 0
    if epochs.any() and len(epochs) > 1:
        for i in range(len(epochs)):
            for j in range(i+1,len(epochs)):
                L += a*abs(np.log10(peaks[i])-np.log10(peaks[j])) + b*abs(vlsrs[i]-vlsrs[j]) + c*abs(fwhms[i]-fwhms[j])
                # L+= a*abs(np.log10(peaks[i])/np.log10(peaks[j])) + b*abs(2*(vlsrs[i]-vlsrs[j])/(vlsrs[i]+vlsrs[j]))+c*abs(fwhms[i]/fwhms[j])
        L = L/((a+b+c)*sum(np.arange(1,len(epochs))))
        return L
    else:
        return 'N/A'


app.layout = html.Div([
html.Div([html.H1(children='IMProV', style={'textAlign': 'center'}),
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
        dcc.Store(id = 'main-data-store')
          ]),
        html.Div([
        html.Div(
            dcc.Graph(id='g1', config={'displayModeBar': True, 'scrollZoom': True}, style={'height': '80vh'}),
            className='five columns', style={"border": "2px black solid"}),
            dcc.Store(id='g1-selection-data'),
        html.Div(
            dash_table.DataTable(id="table_Gauss_pars",columns = [{"name": i,'id':i} for i in cols],
                                 style_table={'overflowX': 'scroll'},export_format="csv"),
                 className= 'four columns', style={"border": "2px black solid"}),
            dcc.Store(id='all_feature_info'),
            dcc.Store(id='proper_motions'),
        html.Div([
                html.Div([
                        html.Button("SAVE FEATURES", id='btn_features', n_clicks=0),
                ],className='rows'),
                html.Div([
                    html.Button("CALCULATE PROPER MOTIONS", id='btn_proper_motions', n_clicks=0)
                ],className='rows'),

                html.Div([dcc.Markdown(r"The likeness parameter $L$ is calculated by:" + "\n"
                    r'$L_{ij} = a|log_{10}A_i-log_{10}A_j| + b|\mu_i-\mu_j| + '+
                    r'c |FWHM_i-FWHM_j|$' + '\n \n' +
                    r'with $L = \frac{1}{(a+b+c)\Sigma_n^N (n-1)}\Sigma_i^N\Sigma_{j>i}^N L_{ij}$,' + '\n' +
                    r'$a = 1, b = 1, c = 1$',mathjax=True,
                                       style={"white-space": "pre","overflow-x": "scroll"})],className='rows'),
                html.Div([
                    dcc.Markdown('L = ',style={"white-space": "pre","overflow-x": "scroll"}, mathjax=True),
                    html.Div(children='',id = 'Likeness Parameter')
                ],className='rows'),
                html.Div(id = 'Feature_output',children='')
        ],className='three columns', style={"border": "2px black solid"})


])
])
@app.callback(
    Output('g1-selection-data','data'),
    Output('Likeness Parameter','children'),
    Input('g1','selectedData')
)
def save_selection(selection):
    inds = []
    if selection:
        for i in range(len(selection['points'])):
            inds.append(selection['points'][i]['pointIndex'])
        if inds:
            if type(calculate_likeness(df,inds))== float:
                l_string = str(round(calculate_likeness(df,inds),7))
            else:
                l_string = calculate_likeness(df,inds)
            return [inds,l_string]
    else:
        raise PreventUpdate



@app.callback(
    Output('g1','figure'),
    Input('main-data-store','data'),
    Input('proper_motions','data'),
    State('all_feature_info','data')
             )
def show_figure(data,proper_motions,feature_data):

    if data:
        print(feature_data)
        dfdata = pd.DataFrame.from_dict(data)
        fig = px.scatter(dfdata, x="RA", y='DEC', color='Epoch', hover_data={
            'Epoch': True,
            'RA': False,
            'DEC': False,
            'DRA': False,
            'DDEC': False,
            'VLSR': True,
            'Peak': True,
            'Peak err': False,
            'Center': False,
            'Center err': False,
            'FWHM': True,
            'FWHM err': False
        })
        dot_opacity = np.ones(len(dfdata['RA']))
        if feature_data:
            for i in range(len(feature_data)):
                inds = feature_data[i][0]
                dot_opacity[inds] = 0
        fig.update_traces(marker=dict(opacity=dot_opacity))
        fig.update_layout(clickmode='event+select')
        if proper_motions:
            ra = np.array(proper_motions['RA'],dtype='float')
            dec = np.array(proper_motions['Dec'],dtype='float')
            mu_ra = np.array(proper_motions['mu_alpha'],dtype='float')
            mu_dec= np.array(proper_motions['mu_delta'],dtype='float')
            fig2 = ff.create_quiver(x = ra,y = dec,u = mu_ra,v= mu_dec,scale = 100)
            fig.add_traces(data=fig2.data)

        return fig
    else:
        raise PreventUpdate
@app.callback(Output('output-data-upload', 'children'),
              Output('main-data-store','data'),
              Input('upload-data', 'contents'),
              State('upload-data', 'filename'))
def update_output(list_of_contents, list_of_names):
    global df
    global all_dates
    dfcolumns = ['Epoch','RA', 'DEC', 'DRA', 'DDEC', 'VLSR', 'Peak', 'Peak err', 'Center', 'Center err', 'FWHM', 'FWHM err']
    df = pd.DataFrame(columns= dfcolumns)
    all_dates = []
    if list_of_contents is not None:
        tuples = zip(*sorted(zip(list_of_names,list_of_contents)))
        list_of_names,list_of_contents = [list(tuple1) for tuple1 in tuples]
        for i,(content,name) in enumerate(zip(list_of_contents,list_of_names)):
            df1 = parse_contents(content , name)
            ind1 = name.find('EP')
            all_dates.append(epdate_float[int(name[ind1+2])-1])
            df1['Epoch'] = np.array(np.zeros(len(df1['RA']))+i+1,dtype='int')
            df = pd.concat([df,df1])
        return ['Successfully uploaded the files: '+str(list_of_names),dict(df)]
    else:
        return ['Please upload a file.',[]]

@app.callback(
        Output("Feature_output","children"),
        Output('all_feature_info','data'),
        Input('btn_features','n_clicks'),
        State('g1-selection-data','data'),
        State('Likeness Parameter','children'),
        State('all_feature_info','data')
             )
def save_feature(n_clicks,data,likenes_par,feature_data):
    changed_id = [p['prop_id'] for p in callback_context.triggered][0]
    if changed_id == 'btn_features.n_clicks':
        all_epochs = df['Epoch'].values[data]
        all_epochs_set = set(all_epochs)
        if len(all_epochs_set) != len(all_epochs):
            return ['Please choose only one datapoint per epoch. Data not saved.',feature_data]
        else:
            if feature_data:
                fets = feature_data
                fets.append([data,likenes_par])
            else:
                fets =[[data, likenes_par]]
            return ['Data points saved. ',fets]
    else:
        raise PreventUpdate


# @app.callback(
#     Output('')
# )


@app.callback(
    Output('table_Gauss_pars','data'),
    Output('proper_motions','data'),
    Input('all_feature_info','data')
)
def update_table(data):
    if data:
        l_pars = []
        epochs = []
        ids = list(range(1, len(data) + 1))
        ra = []
        raerr = []
        dec = []
        decerr = []
        vlsr = []
        mu_alpha = []
        dmu_alpha = []
        mu_delta = []
        dmu_delta = []
        ep_id = []
        for i in range(len(data)):
            l_pars.append(round(data[i][1],5))
            inds = data[i][0]
            t = np.array(all_dates)[list(df['Epoch'].values[inds]-1)]
            ep_id.append(str(list(df.index.values[inds])))
            ra.append(str(round(df['RA'].values[inds][np.argmin(df['Epoch'].values[inds])],5)))
            raerr.append(str(round(df['DRA'].values[inds][np.argmin(df['Epoch'].values[inds])], 5)))
            dec.append(str(round(df['DEC'].values[inds][np.argmin(df['Epoch'].values[inds])], 5)))
            decerr.append(str(round(df['DDEC'].values[inds][np.argmin(df['Epoch'].values[inds])], 5)))
            epochs.append(str(df['Epoch'].values[inds]))
            vlsr.append(str(round(np.mean(df['VLSR'].values[inds]),5)))
            mu,dmu = calc_mu(df['RA'].values[inds],df['DRA'].values[inds],t)
            mu_alpha.append(str(round(mu*np.cos(np.deg2rad(source_Dec)),5)))
            dmu_alpha.append(str(round(dmu*np.cos(np.deg2rad(source_Dec)),5)))
            mu,dmu = calc_mu(df['DEC'].values[inds],df['DDEC'].values[inds],t)
            mu_delta.append(str(round(mu,5)))
            dmu_delta.append(str(round(dmu,5)))

        df_dict = pd.DataFrame.from_dict(
            {'ID': ids,"L value":l_pars, 'Epochs': epochs,"RA":ra,"RAerr":raerr,'Dec':dec,'Decerr':decerr,
             'mu_alpha':mu_alpha,'dmu_alpha':dmu_alpha,'mu_delta':mu_delta,'dmu_delta':dmu_delta,'VLSR':vlsr,'EpID':ep_id})
        return [df_dict.to_dict('records'),dict(df_dict)]
    else:
        raise PreventUpdate
if __name__ == '__main__':
    app.run_server(debug=True)
