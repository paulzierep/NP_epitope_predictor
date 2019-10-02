import pandas as pd
import json

# <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/jquery.dataTables.min.css">
# <script type="text/javascript" charset="utf8" src="https://code.jquery.com/jquery-3.3.1.js"></script>
# <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"></script>

# <script>

#  $(document).ready(function() {{
#     $('table.dataframe').DataTable();
# }} );

# </script>

pd.set_option('display.max_colwidth', -1)

def Results_To_Json(results):

    # for k,v in results.items():
    # print(k)
    # for k2,v2 in v.items():
    #     print(type(v2))
    

    results['cluster_info']['Ontology'] = results['cluster_info']['Ontology'].to_dict()
    results['classify_info_b']['fp_imp'] = results['classify_info_b']['fp_imp'].to_dict()
    results['classify_info_t']['fp_imp'] = results['classify_info_t']['fp_imp'].to_dict()

    json_results = json.dumps(results)

    return(json_results)

def Results_To_Html(results):

    #.to_html(escape=False)

    c_inf = results['cluster_info']
    c_inf['Ontology'] = c_inf['Ontology'].to_html(escape=False)

    b_cell = results['classify_info_b']
    b_cell['fp_imp'] = b_cell['fp_imp'].to_html(escape=False)

    t_cell = results['classify_info_t']
    t_cell['fp_imp'] = t_cell['fp_imp'].to_html(escape=False)

    #print(c_inf["Cluster"]

    top = """

    <!DOCTYPE html>
    <html>
    <head>
    <title>Example</title>
    </head>
    <body>

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

    </body>
    </html>

    """.format(**t_cell)

    html_text = top + b_cell + t_cell

    #html_text = html_text.replace('border="1"','border="0"')

    return(html_text)