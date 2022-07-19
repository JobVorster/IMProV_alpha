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

import base64
import datetime
import io

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
plotted = False
app = Dash(__name__, external_stylesheets=external_stylesheets)
dfcolumns = ['RA min','RA max','Dec min','Dec max','Peak','Peak err','Center','Center err','FWHM','FWHM err']
# make a sample data frame with 6 columns
df = []#pd.read_csv("Ep1_VERA_Shifted.txt", sep=',')
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
        dcc.Store(id = 'main-data-store')]),
    html.Div([

        html.Div(
            dcc.Graph(id='g1', config={'displayModeBar': True, 'scrollZoom': True}, style={'height': '80vh'}),
            className='four columns', style={"border": "2px black solid"}),
            dcc.Store(id='g1-relayout-data'),
        html.Div
            ([
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
                    dcc.Store(id='Gauss_Fit_results_data'),
                    html.Button("SAVE CURRENT GAUSS FIT", id='btn_curvefit', n_clicks=0),
                ], className='row', style={"border": "2px black solid", 'height': '31.2vh'})
        ], className='four columns', style={"border": "2px black solid"}),
        html.Div([
            dash_table.DataTable(id="table_Gauss_pars",columns = [{"name": i,'id':i} for i in dfcolumns],export_format="csv"),
            html.Div(id = 'test_div'),
            dcc.Store(id='store_Gauss_Fit_all_results'),
            dcc.Store(id = 'store_Table')
        ], className='four columns', style={"border": "2px black solid", 'overflowY': 'auto'})
    ], className='row', style={"border": "2px black solid"})
])

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


def isolate(x, y, xlim, ylim):
    '''Gives indices for elements in x,y range.
    Parameters :
    ------------
    x : array-like
    x coordinates of element set.
    y : array-like
    y coordinates of element set.
    xlim : array-like of form :[xlimmin,xlimmax]
    specify maximum and minimum x values for elements to isolate.
    ylim : array-like of form :[ylimmin,ylimmax]
    specify maximum and minimum y values for elements to isolate.

    Returns :
    --------
    ind : array-like
    Array containing indices of isolated elements. Can then be easily called, e.g. x-coords : x[ind]'''
    if len(xlim) == 0 or len(ylim) == 0:
        return []
    else:
        ind = np.where(
            np.logical_and(np.logical_and(y >= ylim[0], y <= ylim[1]), np.logical_and(x <= xlim[1], x >= xlim[0])))
        return ind[0]


def gauss(x, A, mu, sigma):
    return A * np.exp(-1 * (x - mu) ** 2 / (2 * sigma ** 2))


def get_figure(df, x_col, y_col, selectedpoints, selectedpoints_local,fitted_lims,relayout_data):
    # global plotted
    # if not plotted:
    if dict(df):
        fig = px.scatter(df, x=df[x_col], y=df[y_col],error_x = df['RAERR'].values,error_y=df['DECERR'].values, color_continuous_scale='Jet')
        if fitted_lims:
            for xlim,ylim in fitted_lims:
                fig.add_shape(type='rect',x0=xlim[0],x1=xlim[1],y0=ylim[0],y1=ylim[1],fillcolor='green',opacity=0.1)

        fig.update_coloraxes(autocolorscale=True, colorscale='Jet', showscale=True)
        fig.update_traces(selectedpoints=selectedpoints,
                          customdata=df.index,
                          mode='markers', marker={'color': df['VLSR'], 'size': 12},
                          unselected={'marker': {'opacity': 1}})
        # fig.update_yaxes(scaleanchor='x',scaleratio=1) #This options makes an equal aspect ratio.
        fig.update_layout(margin={'l': 20, 'r': 0, 'b': 15, 't': 5}, dragmode='select', colorscale_diverging='Jet',hovermode=False)
        if relayout_data:
            if 'xaxis.range[0]' in relayout_data.keys():
                fig.update_xaxes(range=[relayout_data['xaxis.range[0]'], relayout_data['xaxis.range[1]']])
                fig.update_yaxes(range=[relayout_data['yaxis.range[0]'], relayout_data['yaxis.range[1]']])
        fig.update_shapes(dict(xref='x', yref='y'))
        # if all_inds:
        #     fig.add_trace(go.Scatter(x=df[x_col][all_inds],y=df[y_col][all_inds],mode='markers'))
        plotted = True
        return fig
    else:
        raise PreventUpdate


def get_figure2(df, x_col, y_col, selectedpoints, selectedpoints_local, do_Gaussian,all_inds, Gauss_pars=[]):
    if dict(df):
        fig = px.scatter(df, x=df['VLSR'][selectedpoints], y=df['FLUX'][selectedpoints],error_y = df['DFLUX'][selectedpoints])
        vlsr_ranges = [min(df['VLSR'].values[selectedpoints]), max(df['VLSR'].values[selectedpoints])]
        flux_ranges = [min(df['FLUX'][selectedpoints]), max(df['FLUX'][selectedpoints])]
        fig.update_traces(selectedpoints=selectedpoints,
                          customdata=df.index,
                          mode='markers', marker={'color': df['VLSR'], 'size': 12}, marker_colorscale='Jet',
                          marker_cmax=vlsr_ranges[1], marker_cmin=vlsr_ranges[0],
                          unselected={'marker': {'opacity': 1}, 'textfont': {'color': 'rgba(0, 0, 0, 0)'}})
        if do_Gaussian:
            popt = Gauss_pars[0]
            pcov = Gauss_pars[1]
            x = np.linspace(vlsr_ranges[0], vlsr_ranges[1], 10000)
            fig.add_trace(go.Scatter(x=x, y=gauss(x, *popt)))
        fig.update_xaxes(
            range=[vlsr_ranges[0] - 0.1 * np.diff(vlsr_ranges)[0], vlsr_ranges[1] + 0.1 * np.diff(vlsr_ranges)[0]])
        fig.update_yaxes(range=[0 - 0.07 * flux_ranges[1], flux_ranges[1] + 0.1 * flux_ranges[1]])

        fig.update_layout(margin={'l': 20, 'r': 0, 'b': 15, 't': 5}, dragmode=False,
                          hovermode=False,showlegend=False)
        return fig
    else:
        raise PreventUpdate


# this callback defines 3 figures
# as a function of the intersection of their 2 selections
@app.callback(
    [Output('g1', 'figure'),
     Output('g2', 'figure'),
     Output('Gauss_Fit_results_data', 'data')],
    Input('g1', 'selectedData'),
    Input('g2', 'selectedData'),
    Input('main-data-store','data'),
    Input("Gauss_Fit_checkbox", 'value'),
    State('g1-relayout-data','data'),
    State('store_Gauss_Fit_all_results','data')
)
def callback(selection1, selection2, df, bool_Gaussian, relayout_data ,already_fitted):
    if dict(df):
        df = pd.DataFrame.from_dict(df)
    if 'bool_Gaussian' in bool_Gaussian:
        do_Gaussian = True
        Gauss_pars = []
    else:
        do_Gaussian = False
        Gauss_pars = []
        points_range = []
    if dict(df):
        selectedpoints = df.index
    else:
        selectedpoints = []
    for selected_data in [selection1, selection2]:
        if selected_data and selected_data['points']:
            selectedpoints = np.intersect1d(selectedpoints,
                                            [p['customdata'] for p in selected_data['points']])
    fitted_lims = []
    all_inds = []

    if already_fitted:
        for i in range(len(already_fitted)):
            xlim = [already_fitted[i]['RA min'],already_fitted[i]['RA max']]
            ylim = [already_fitted[i]['Dec min'],already_fitted[i]['Dec max']]
            fitted_lims.append([xlim,ylim])

        for lim in fitted_lims:
            inds = isolate(df['RA'].values,df['DEC'].values,lim[0],lim[1])
            all_inds.extend(inds)
    points_range = []
    if do_Gaussian :
        fluxes = np.array(df['FLUX'][selectedpoints])
        vlsrs = np.array(df['VLSR'][selectedpoints])
        dfluxes = np.array(df['DFLUX'][selectedpoints])
        p0 = [max(fluxes), vlsrs[np.argmax(fluxes)], 1]
        Gauss_pars = optim.curve_fit(gauss, xdata=vlsrs, ydata=fluxes,absolute_sigma=True,sigma=dfluxes,
                                     p0=p0)  # I still have to add the calculation with the errors.
        points_range = selection1['range']

    return [get_figure(df, "RA", "DEC", selectedpoints, selection1,fitted_lims,relayout_data),
            get_figure2(df, "VLSR", "FLUX", selectedpoints, selection2, do_Gaussian, all_inds, Gauss_pars),
            [points_range, Gauss_pars]
            ]


@app.callback(
    [Output('cont_RA min', 'children'),
     Output('cont_RA max', 'children'),
     Output('cont_Dec min', 'children'),
     Output('cont_Dec max', 'children'),
     Output('cont_Peak', 'children'),
     Output('cont_Peak err', 'children'),
     Output('cont_Center', 'children'),
     Output('cont_Center err', 'children'),
     Output('cont_FWHM', 'children'),
     Output('cont_FWHM err', 'children')],
    [Input("Gauss_Fit_checkbox", 'value'),
     Input('Gauss_Fit_results_data', 'data')]
)
def show_Gaussian_results(value, data):
    if "bool_Gaussian" in value:
        if len(data) == 0:
            return ['a\t'] * 10
        else:
            lims = data[0]
            xlims = list(np.around(np.array(lims['x'], dtype=float), decimals=6))
            ylims = list(np.around(np.array(lims['y'], dtype=float), decimals=6))
            gauss_pars = data[1]
            popt = list(np.around(gauss_pars[0], decimals=6))
            pcov = list(np.around(np.sqrt(np.diag(gauss_pars[1])), decimals=6))
            Peaks = [popt[0], pcov[0]]
            Centers = [popt[1], pcov[1]]
            FWHMs = list(np.around(np.array([popt[2] * 2 * np.sqrt(2 * np.log(2)), pcov[2] * 2 * np.sqrt(2 * np.log(2))]),decimals=6))
            return_list = []
            for lists in [xlims, ylims, Peaks, Centers, FWHMs]:
                return_list.extend(lists)
            return return_list
    else:
        return ['0.000000 \n'] * 10


#Callback to store the gauss fit parameters into a store dataframe.
@app.callback(
    Output('store_Gauss_Fit_all_results','data'),
    Output("Gauss_Fit_checkbox", 'value'),
    Input('btn_curvefit','n_clicks'),
    State('store_Gauss_Fit_all_results','data'),
    State('Gauss_Fit_results_data', 'data'),
    State("Gauss_Fit_checkbox", 'value')
             )
def store_all_gauss_fits(n_clicks,data_all,data_now,value):
    changed_id = [p['prop_id'] for p in callback_context.triggered][0]
    if changed_id=='btn_curvefit.n_clicks':
        if "bool_Gaussian" in value:
            if not data_all:
                data_all = pd.DataFrame(columns = ['RA min','RA max','Dec min','Dec max','Peak','Peak err','Center','Center err','FWHM','FWHM err'])
            else:
                data_all = pd.DataFrame.from_dict(data_all)
            app_dict = pd.DataFrame(columns=data_all.columns)
            lims = data_now[0]
            xlims = list(np.around(np.array(lims['x'], dtype=float), decimals=6))
            ylims = list(np.around(np.array(lims['y'], dtype=float), decimals=6))
            gauss_pars = data_now[1]
            popt = list(np.around(gauss_pars[0], decimals=6))
            pcov = list(np.around(np.sqrt(np.diag(gauss_pars[1])), decimals=6))
            Peaks = [popt[0], pcov[0]]
            Centers = [popt[1], pcov[1]]
            FWHMs = list(
                np.around(np.array([popt[2] * 2 * np.sqrt(2 * np.log(2)), pcov[2] * 2 * np.sqrt(2 * np.log(2))]),
                          decimals=6))
            return_list = []
            for lists in [xlims, ylims, Peaks, Centers, FWHMs]:
                return_list.extend(lists)
            for i,column in enumerate(data_all.columns):
                app_dict[column] = [return_list[i]]

            data_all1 = pd.concat([data_all,app_dict])
            return [data_all1.to_dict('records'),[]]
        else:
            raise PreventUpdate
    else:
        raise PreventUpdate

# @app.callback(
#     [Output('table_Gauss_pars','data'),
#      Output('table_Gauss_pars','columns')],
#     Input('store_Table','data')
# )
# def show_pars_data(data):
#     return [data,['1']]
# @app.callback(
#         [Output('store_Gauss_Fit_all_results','data')],
#         [Input('store_Table','data')]
#              )
# def move_data(data):
#     changed_id = [p['prop_id'] for p in callback_context.triggered][0]
#     if "store_Table.data" in changed_id:
#         return [data]
#     else:
#     else:
#         raise PreventUpdate


@app.callback(
        Output('table_Gauss_pars','data'),
        Input('store_Gauss_Fit_all_results','data')
              )
def update_table(data):
        return data

@app.callback(Output('output-data-upload', 'children'),
              Output('main-data-store','data'),
              Input('upload-data', 'contents'),
              State('upload-data', 'filename'))
def update_output(list_of_contents, list_of_names):
    global df
    if list_of_contents is not None:
        df = parse_contents(list_of_contents[0],list_of_names[0])
        return ['Successfully uploaded the file: '+str(list_of_names[0]),dict(df)]
    else:
        return ['Please upload a file.',[]]

@app.callback(
        Output('g1-relayout-data','data'),
        Input('g1','relayoutData')
             )
def update_relayoutdata(relayoutdata):
    if relayoutdata:
        if "xaxis.range[0]" in relayoutdata.keys():
            return relayoutdata
        else:
            raise PreventUpdate
    else:
        raise PreventUpdate


if __name__ == '__main__':
    app.run_server(debug=True)
