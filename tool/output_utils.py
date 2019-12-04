import pandas as pd
import json
import copy
import numpy as np
# <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/jquery.dataTables.min.css">
# <script type="text/javascript" charset="utf8" src="https://code.jquery.com/jquery-3.3.1.js"></script>
# <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"></script>

# <script>

#  $(document).ready(function() {{
#     $('table.dataframe').DataTable();
# }} );

# </script>

pd.set_option('display.max_colwidth', -1)

def df_float_formatting(df):
    """
    Casts all float formats in a df to precision 2
    """

    def cast_and_format(val):

        if isinstance(val, np.floating) or isinstance(val, float):
            #print('{0:.2f}'.format(val))
            return('{0:.2f}'.format(val))
        else:
            return(val)

    df = df.applymap(lambda x: cast_and_format(x)) 

    return(df)


def Results_To_Json(results):

    #make a copy of the dict and its items
    #otherwise the .to_dict() will change the df and 
    #the Results_To_Html does not work
    results_copy = copy.deepcopy(results)
    
    if results["error"]:
        return(json.dumps(results))

    results_copy['cluster_info']['Ontology'] = results_copy['cluster_info']['Ontology'].to_dict()
    results_copy['classify_info_b']['fp_imp'] = results_copy['classify_info_b']['fp_imp'].to_dict()
    results_copy['classify_info_t']['fp_imp'] = results_copy['classify_info_t']['fp_imp'].to_dict()

    json_results = json.dumps(results_copy)

    return(json_results)


def Results_To_Django(results):

    #make a copy of the dict and its items
    #otherwise the .to_dict() will change the df and 
    #the Results_To_Html does not work
    results_copy = copy.deepcopy(results)
    
    if results["error"]:
        return(results)

    # results_copy['cluster_info']['Ontology'] = results_copy['cluster_info']['Ontology'].T.to_dict()
    # results_copy['classify_info_b']['fp_imp'] = results_copy['classify_info_b']['fp_imp'].T.to_dict()
    # results_copy['classify_info_t']['fp_imp'] = results_copy['classify_info_t']['fp_imp'].T.to_dict()

    results_copy['cluster_info']['Ontology'] = \
    results_copy['cluster_info']['Ontology'].rename(columns={
    "ChEBI_ID": "ChEBI ID",
    "ChEBI_Name": "ChEBI Name",
    "Corr-PValue": "Corr. p",
    "PValue": "p",
    "Fold": "Fold",
    "SamplePercentage": "Sample Percentage",
     })


    results_copy['cluster_info']['Ontology'] = df_float_formatting(results_copy['cluster_info']['Ontology']).to_html(escape=False, index = False)
    results_copy['classify_info_b']['fp_imp'] = df_float_formatting(results_copy['classify_info_b']['fp_imp']).to_html(escape=False,)
    results_copy['classify_info_t']['fp_imp'] = df_float_formatting(results_copy['classify_info_t']['fp_imp']).to_html(escape=False,)

    # results_copy = json.dumps(results_copy)

    return(results_copy)

# <link rel="stylesheet" href="https://cdn.datatables.net/1.10.19/css/jquery.dataTables.min.css">

# </head>
# <body>

# <script src="https://code.jquery.com/jquery-3.3.1.js"></script>
# <script src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"></script>

# <script>

# $(document).ready(function() {
#     $('table.dataframe').DataTable();
# } );

# </script>

def Results_To_Html(results):

    top = """
    <!DOCTYPE html>
    <html>
    <head>
    <title>Example</title>

    </head>
    <body>
    """

    end = """
    </body>
    </html>
    """

    if results["error"]:

        error ="""
        <h3>Parsing Error: </h3>
        {0}
        """.format(results["error"].replace("\n","<br>"))

        html_text = top + error + end

        return(html_text)


    c_inf = results['cluster_info']
    c_inf['Ontology'] = df_float_formatting(c_inf['Ontology'])
    c_inf['Ontology'] = c_inf['Ontology'].to_html(escape=False)

    b_cell = results['classify_info_b']
    b_cell['fp_imp'] = df_float_formatting(b_cell['fp_imp'])
    b_cell['fp_imp'] = b_cell['fp_imp'].to_html(escape=False)

    t_cell = results['classify_info_t']
    t_cell['fp_imp'] = df_float_formatting(t_cell['fp_imp'])
    t_cell['fp_imp'] = t_cell['fp_imp'].to_html(escape=False)

    #print(c_inf["Cluster"]

    input_svg = """
    <h3>Input: </h3>
    {0}
    """.format(results['input_svg'])

    cluster = """
    <h3>Cluster: {Cluster}</h3>
    <h3>Members: {Members}</h3>
    <h3>Warning: {Warning}</h3>
    <h3>Name: {Name}</h3>
    <h3>Ontology description:</h3>
    {Ontology}

    """.format(**c_inf)

    b_cell = """
    <h3>B cell probability: {proba}</h3>
    <h3>Info: {fitting_info}</h3>
    <h3>Known epitopes: {p_support}</h3>
    <h3>Important features: </h3>

    {fp_imp}

    """.format(**b_cell)

    t_cell = """
    <h3>T cell probability: {proba}</h3>
    <h3>Info: {fitting_info}</h3>
    <h3>Known epitopes: {p_support}</h3>
    <h3>Important features: </h3>

    {fp_imp}
    """.format(**t_cell)

    end = """
    </body>
    </html>
    """

    html_text = top + input_svg + cluster +  b_cell + t_cell

    #html_text = html_text.replace('border="1"','border="0"')

    return(html_text)