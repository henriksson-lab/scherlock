from IPython.display import HTML, display
from random import randint
import plotly.express as px
import matplotlib
import matplotlib.cm
import plotly.express as px
import matplotlib





###################################################################
###################################################################
###################################################################

# Pretty Differential Expression table

###################################################################
###################################################################
###################################################################




def tabulate_data_frame_in_html(df):
    """Make an HTML table out of a pandas data frame, along with clever formatting.

    Args:
        df (dataframe): The dataframe to turn into HTML

    Returns:
        string: The table as HTML
    """

    ## Format a table row
    def format_row(df, ii, m):

        import validators
        
        # Loop over columns
        row = df.iloc[ii,:]
        row_elements = []
        for jj in range(m):
            if validators.url(str(row[jj])):
                row_elements.append('<td><a target="_blank" href="' + str(row[jj]) + '">link</a></td>\n')
            else:
                row_elements.append('<td>' + str(row[jj]) + '</td>\n')
        html_row = ''.join(row_elements)
        
        return(html_row)

    
    ## Format the table
    n = df.shape[0]
    m = df.shape[1]
    table_elements = []

    # Header
    header_elements = []
    header_elements.append('<table><tr>\n')
    for jj in range(m):
        header_elements.append('<th>' + df.columns[jj] + '</th>')
    header_elements.append('</tr>\n')
    html_header = ''.join(header_elements)
    table_elements.append(html_header)

    # Body
    for ii in range(n):
        table_elements.append('<tr>\n' + format_row(df, ii, m) + '</tr>\n')

    # Table end
    table_elements.append('</table>\n')

    # Return table
    html_table = ''.join(table_elements)
    return(html_table)
#     display(HTML(html_table))
    
    
    
    
def multipage_html_table(panel_dict):
    """Makes a multipanel (multipage) HTML table out of a dictionary of dataframes

    Args:
        panel_dict (Dictionary): Keys are names of panels, 
                    and values are pandas dataframes with the data to be
                    displayed in panels
    Returns:
        string: The multipage HTML table
    """
    
    ## The buttons above the panel
    all_tab_buttons=''.join(
        ["<button class=\"INTANCEID_tablinks\" onclick=\"INTANCEID_open_genepanel(event, 'INTANCEID_"+
         str(panel_name)+"')\">"+str(panel_name)+"</button>" for panel_name in panel_dict.keys()])
    all_tab_buttons = '<div class="INTANCEID_tab">{}</div>'.format(all_tab_buttons)

    ## The panels themselves
    all_panels=''.join(
        ['<div id="INTANCEID_' + str(panel_name) + 
         '" class="INTANCEID_tabcontent">{}</div>'.format(tabulate_data_frame_in_html(panel_dict[panel_name])) 
         for panel_name in panel_dict.keys()])
    
    ## Put all the HTML together
    js="""
    <style>
        /* Tab as a whole */
        .INTANCEID_tab {
          overflow: hidden;
          border: 1px solid #ccc;
          background-color: #f1f1f1;
        }

        /* Tab buttons */
        .INTANCEID_tab button {
          background-color: inherit;
          float: left;
          border: none;
          outline: none;
          cursor: pointer;
          padding: 14px 16px;
          transition: 0.3s;
        }

        /* Change background color of buttons on hover */
        .INTANCEID_tab button:hover {
          background-color: #ddd;
        }

        /* Create an active/current tablink class */
        .INTANCEID_tab button.active {
          background-color: #ccc;
        }

        .INTANCEID_tabcontent {
          /*display: none;*/
          padding: 6px 12px;
          border: 1px solid #ccc;
          border-top: none;
        }
    </style>

    <script>
        function INTANCEID_open_genepanel(evt, panelName) {
          var i;

          // Hide all tabs
          var tabcontent = document.getElementsByClassName("INTANCEID_tabcontent");
          for (i = 0; i < tabcontent.length; i++) {
            tabcontent[i].style.display = "none";
          }

          // Set all tab buttons as inactive
          var tablinks = document.getElementsByClassName("INTANCEID_tablinks");
          for (i = 0; i < tablinks.length; i++) {
            tablinks[i].className = tablinks[i].className.replace(" active", "");
          }

          // Show selected tab and set tab button to active
          document.getElementById(panelName).style.display = "block";
          evt.currentTarget.className += " active";
        } 

    </script>
    """
    js+=all_tab_buttons
    js+=all_panels


    ## Make one panel selected by default
    js+="""
    <script>
    tablinks = document.getElementsByClassName("INTANCEID_tablinks");
    tablinks[0].click();
    </script>
    """

    ## Add big random number to the js/html elements to avoid naming conflicts.
    ## Probably not industry standard ;) but will hopefully work often enough.
    js = js.replace("INTANCEID_","foo"+str(randint(100000, 200000))+"_")
    
    
    ## Hide our horrible formatting ;)
    from bs4 import BeautifulSoup as bs
    aa = bs(js)  # Read html from string
    pretty_js = aa.prettify()   # prettify
    
    return(pretty_js)









def build_pretty_DE_table(adata, n_genes = 25, 
                          adata_var_columns_to_include = None,
                          rank_genes_groups_keys_to_include = ['names', 'logfoldchanges', 'pvals'],
                          gene_ID_column = None, gene_name_column = None,
                          show_table=True):
    """Build and render an improved Differential Expression (DE) table.
    
    Builds an improved DE table including
    - Additional useful gene statistics
    - Direct links to common gene and protein databases
    - Users choice of additional data from adata.var and adata.uns['rank_genes_groups']
    and renders it in html.
    It is a multipage table where every page corresponds to
    one leiden cluster in the data.

    Args:
        adata (AnnData object): The anndata/count table
        n_genes (int):  Number of top DE genes to include in table
        adata_var_columns_to_include ([strings]): 
             Column names of columns from adata.var to include in DE table
        rank_genes_groups_keys_to_include ([strings]):
             Keys from adata.uns['rank_genes_groups']. Listed keys will be included in DE table
        gene_ID_column (string):
             Column in adata.var with Ensemble gene_IDs
             If both this and gene_name_column is provid
             the database links will attempt to be a bit smarter
        gene_name_column (string):
             Column in adata.var with HGNC gene symbols
             If both this and gene_ID_column is provided,
             the database links will attempt to be a bit smarter
        show_table (bool):
             If True, the table will be rendered.
             If False, the table will be returned as a dictionary.

    Returns:
        Only if show_table is False.
        dictionary:
            Each key:value pair corresponds to one page in table.
            Keys are the names of the leiden categories in
            adata.obs.leiden.cat.categories, or single key = '0' if it doesnt exist.
            Values are [n_genes x total number of columns] data frames 
            with data for the corresponding pages.
    """

    # To do:
    # - Ordering of columns in table. How is best?
    # - Test if it works when adata.obs.leiden.cat.categories doesnt exist


    # Check that user havent overwritten adata.var.index with something
    # after running scanpy.tl.rank_genes_groups()
    # That is hard without knowing the code that calls this
    # function, so doing the next best and checking that all the
    # genes in adata.uns['rank_genes_groups']['names'] are still present
    # in adata.var.index
    d = {gene: 1 for gene in adata.var.index}  # For O(1) time lookup
    l = []
    for sublist in adata.uns['rank_genes_groups']['names']:
        l = l + list(sublist)
    if not all([gene in d for gene in l]):
        raise Exception("Error: \n" +
              "Failed to find at least one gene from adata.uns['rank_genes_groups']['names'] \n" +
              "in adata.var.index. \n" +
              "One possible cause is if adata.var.index was changed after \n" +
              "scanpy.tl.rank_genes_groups() was run. That would break the connection \n" +
              "between adata.uns['rank_genes_groups']['names'] and adata.var.index \n" + 
              "If that is the case, run scanpy.tl.rank_genes_groups() and this function \n" +
              "without changing adata.var.index inbetween")

    # Get cluster/group names
    try:
        groups = adata.obs.leiden.cat.categories
    except:
        groups = ['0']
    
    # Calculate some default stats
    percent_cells = adata.var.n_cells/adata.shape[0]
    global_means = pd.Series(np.sum(adata.X, axis=0)/adata.shape[0], index=adata.var.index)

    # Loop over clusters/groups/pages. Build dict of pages
    page_dict = {}
    for ii_group in range(len(groups)):

        # Get identifiers for the top DE genes in current page/cluster
        page_genes = [adata.uns['rank_genes_groups']['names'][ii][ii_group] for ii in range(n_genes)]

        # Data from adata.uns['rank_genes_groups']
        column_dict = {}
        for col in rank_genes_groups_keys_to_include:
            column_dict[col] = [adata.uns['rank_genes_groups'][col][ii][ii_group] for ii in range(n_genes)]
        df = pd.DataFrame(column_dict)
        
        # Some stats
        df['percent_cells'] = np.array(percent_cells[page_genes])
        df['mean_expression'] = np.array(global_means[page_genes])

        ### Links to databases
        if gene_ID_column is not None and gene_name_column is not None:
            gene_ids = adata.var[gene_ID_column][page_genes]
            gene_names = adata.var[gene_name_column][page_genes]


            for gene_name, gene_id in zip(gene_names, gene_ids):

                # For link format, see: https://www.genecards.org/Guide/AboutGeneCards
                df['gene_cards_link'] = "https://www.genecards.org/cgi-bin/carddisp.pl?id=" + gene_id + "&id_type=ensembl"

#                 # With the current really shitty implementation, retrieving uniProtIDs on the
#                 # fly takes way too long. Easily fixed, but waiting for response from uniProt
#                 # on how to get the most relevant entry before doing anything.
#                 # If impossible to get relevance for the different protein IDs 
#                 # its better to not use this approach at all.
#                 import datetime
#                 t0 = datetime.datetime.now()
#                 print('Retrieving of uniProt ID started: ' + str(t0))
#                 uniProt_id = get_uniProt_id_from_ensemble_id('ENSG00000187608')
#                 print('Retrieving of uniProt ID finished:' + str(datetime.datetime.now()))
#                 print('It took ' + str(datetime.datetime.now() - t0))
#                 link_uniprot = 'https://www.uniprot.org/uniprot/' + uniProt_id
                df['uniprot_link'] = "https://www.uniprot.org/uniprot/?query="+gene_name+"&sort=score"

                # Not fixed to go straight to gene page yet, so just links to database internal search results
                df['wikipedia_link'] = "https://en.wikipedia.org/w/index.php?search="+gene_name+"&go=Go"
                df['ensembl_link'] = "http://www.ensembl.org/Multi/Search/Results?q="+gene_name+";site=ensembl_all"
        else:
            for gene in page_genes:
                # Search databases for the gene name or whatever was in adata.var.index upon running
                # scanpy.tl.rank_genes_groups
                df['gene_cards_link']  = "https://www.genecards.org/cgi-bin/carddisp.pl?gene="+gene+"&keywords="+gene
                df['uniprot_link']     = "https://www.uniprot.org/uniprot/?query="+gene+"&sort=score"
                df['wikipedia_link']   = "https://en.wikipedia.org/w/index.php?search="+gene+"&go=Go"
                df['ensembl_link']     = "http://www.ensembl.org/Multi/Search/Results?q="+gene+";site=ensembl_all"


        # User requested columns from adata.var
        if adata_var_columns_to_include is not None:
            for column in adata_var_columns_to_include:
                df[column] = np.array(adata.var[column][page_genes])

        # Add current page to dictionary
        page_dict[groups[ii_group]] = df
    
    
    # Turn page dictionary into html
    js = multipage_html_table(page_dict)
    
    # Render or return page dictionary
    if show_table: 
        display(HTML(js))
    else:
        return(page_dict)


    
# # Testing
# build_pretty_DE_table(adata, n_genes = 10, 
#                           adata_var_columns_to_include = ['gene_names', 'from', 'feature_types', 'highly_variable'],
#                           rank_genes_groups_keys_to_include = ['logfoldchanges', 'pvals'],
#                           gene_ID_column = 'gene_ids_isocount', gene_name_column = 'gene_names',
#                           show_table=True)



###################################################################
###################################################################
###################################################################
# Interactive UMAPS ###############################################
###################################################################
###################################################################
###################################################################




def plot_umaps_sidebyside(adata1, adata2, obsname="leiden", palette="Set3"):
    """Plot two UMAPs and allow the user to see which points correspond using the mouse.

    Args:
        adata1 (anndata): The first count data object
        adata2 (anndata): The second count data object
        obsname (string, optional): Name of the .obs column to show
        palette (string, optional): Colors to use for plotting categorical annotation groups.
            Either a name of a matplotlib predefined palette to pick from
            :class:`~matplotlib.colors.ListedColormap`,
            or a list of matplotlib predefined named colors.
            (see :func:`~matplotlib.colors.is_color_like`).

    Returns:
        Nothing; displays interactive HTML
    """

    ### Figure out the color mapping
    if not isinstance(palette, list):
        ppp = matplotlib.cm.get_cmap(palette).colors
    else:
        ppp = palette
    #pp = [matplotlib.colors.rgb2hex(color).replace("#","%23") for color in ppp]
    pp = [matplotlib.colors.rgb2hex(color) for color in ppp]

    ### Turn points from one umap into a dataframe
    def one_adata_tojs(adata):
        cats = adata.obs[obs_name].cat.categories
        color_dict = {cats[ii]: pp[ii] for ii in range(len(cats))}
        x=adata.obsm["X_umap"][:,0].tolist()
        y=adata.obsm["X_umap"][:,1].tolist()
        cols = {'x': x, 'y': y, 
                'bc': adata.obs_names,
                'color':[color_dict[i] for i in adata.obs[obs_name].tolist()]}
        df = pd.DataFrame(data=cols)
        return(df)

    df_ptA=one_adata_tojs(adata1)
    df_ptB=one_adata_tojs(adata2)

    ### Generate the point-to-point correspondence table
    df1=pd.DataFrame(data={'indA':range(0,len(adata1.obs_names)),'bc':adata1.obs_names})
    df2=pd.DataFrame(data={'indB':range(0,len(adata2.obs_names)),'bc':adata2.obs_names})
    df_corresp=df1.set_index('bc').join(df2.set_index('bc'), how="outer").dropna()

    ### Return a value as quoted if it is a string
    def quoteIfNeeded(x):
        if(isinstance(x,str)):
            return('"'+x+'"')
        else:
            return(str(x))

    ### Turn a pandas dataframe into a javascript variable declaration
    def df2js(df,varname):
        outcols=[];
        for (columnName, columnData) in df.iteritems():
            outcols.append(columnName+":[{}]".format(",".join([quoteIfNeeded(x) for x in columnData.values])))
        return("var "+varname+"={"+','.join(outcols)+"};");


    ### Construct all the HTML
    js="""

    <style>
        .INSTANCEID_svg {
          margin: 0px;
        }

        .INSTANCEID_correspline {
            stroke: black;
            stroke-width: 2
        }

        .INSTANCEID_borderrect {
            fill:none;
            stroke:black;
            stroke-width:1
        }
    </style>

    <svg width="0" height="0" class="INSTANCEID_svg">
        <g id="INSTANCEID_svg_points"/>
        <line x1="0" y1="0" x2="0" y2="0" class="INSTANCEID_correspline"/>
        <rect x="0" y="0" width="0" height="0" class="INSTANCEID_borderrect"/>
        <rect x="0" y="0" width="0" height="0" class="INSTANCEID_borderrect"/>
    </svg>


    <script>

        ///////////// Size settings
        var totalw=800;
        var totalh=400;
        var panelw=totalw/2-2;
        var panelh=totalh-2;

        ///////////// Set up sizes of the panel
        var svgpanel = document.getElementsByClassName('INSTANCEID_svg')[0];
        svgpanel.setAttribute("width", totalw);
        svgpanel.setAttribute("height", totalh);
        var svgpanel_points = document.getElementById('INSTANCEID_svg_points');
        var cline = document.getElementsByClassName('INSTANCEID_correspline')[0];
        var borderpanel = document.getElementsByClassName('INSTANCEID_borderrect');
        borderpanel[0].setAttribute("width", panelw);
        borderpanel[1].setAttribute("width", panelw);
        borderpanel[0].setAttribute("height", panelh);
        borderpanel[1].setAttribute("height", panelh);
        borderpanel[1].setAttribute("x", totalw/2);

    """

    ### Add the adata declarations
    js += df2js(df_ptA,"scatterA")+"\n"
    js += df2js(df_ptB,"scatterB")+"\n"
    js += df2js(df_corresp,"corresp")+"\n"
    
    js += """

        ///////////// Function to turn an abstract point into a circle
        function putPoints(pset, pointdata, shiftx, shifty) {
            var data_x=pointdata["x"];
            var data_y=pointdata["y"];
            var data_color=pointdata["color"];
             for(var i=0;i<data_x.length;i++){
                var pt = document.createElementNS("http://www.w3.org/2000/svg","circle");

                var sx=data_x[i]+shiftx;
                var sy=data_y[i]+shifty;

                pt.setAttribute("class", "INSTANCEID_circ");
                pt.setAttribute("cx", sx);
                pt.setAttribute("cy", sy);
                pt.setAttribute("r", 2);
                pt.setAttribute("style", "fill:"+data_color[i]);

                pt.setAttribute("id", "INSTANCEID_pt"+pset+"_"+i);

                pt.addEventListener('mouseover', mouseOverEffect);
                pt.addEventListener('mouseout', mouseOutEffect);

                svgpanel_points.appendChild(pt);
            }
        }

        /////////////// Make quick lookups for point-point correspondence
        var dictCorresp = {};
        var indA=corresp["indA"];
        var indB=corresp["indB"];
        for(var i=0;i<indA.length;i++){
            dictCorresp["INSTANCEID_ptA_"+indA[i]]="INSTANCEID_ptB_"+indB[i];
            dictCorresp["INSTANCEID_ptB_"+indB[i]]="INSTANCEID_ptA_"+indA[i];
        }

        ///////////// Scale point positions onto the panel size
        function scalepoint(arr, theoffset, thescale) {
            var min = Math.min.apply(null, arr);
            var max = Math.max.apply(null, arr);
            for(var i=0;i<arr.length;i++) {
                arr[i] = (arr[i]-min)/(max-min)*thescale+theoffset;
            }
        }
        function transformXY(pointdata) {
            scalepoint(pointdata["x"], 10, panelw-20);
            scalepoint(pointdata["y"], 10, panelh-20);
        }
        transformXY(scatterA);
        transformXY(scatterB);

        ///////////// Put the points into the SVG
        putPoints("A",scatterA, 0, 0);
        putPoints("B",scatterB, totalw/2, 0);


        ///////////// Callback whenever the mouse hovers a point
        function mouseOverEffect() {
            var id=this.getAttribute("id");
            console.log(id);

            var otherId=dictCorresp[id];
            if(otherId==null){
                cline.setAttribute("x1", 0);
                cline.setAttribute("y1", 0);
                cline.setAttribute("x2", 0);
                cline.setAttribute("y2", 0);
            } else {
                var ptA = this;
                var ptB = document.getElementById(otherId);

                cline.setAttribute("x1", ptA.getAttribute("cx"));
                cline.setAttribute("y1", ptA.getAttribute("cy"));
                cline.setAttribute("x2", ptB.getAttribute("cx"));
                cline.setAttribute("y2", ptB.getAttribute("cy"));
            }
        }
        function mouseOutEffect() {
        }

    </script>

    """
    display(HTML(
       js.replace("INTANCEID_","foo"+str(randint(0, 100000))+"_")
    ))










###################################################################
###################################################################
###################################################################
# 3D UMAPS ########################################################
###################################################################
###################################################################
###################################################################



# Functions for continuous 3D umaps

# User exposed function for 3D umap. 
# Handles arguments, loops if many indentifiers given or if identifier is not unique, and calls inner_plot_3d_umap
def plot_3d_umap_continuous(adata, identifiers, column='index', palette='viridis', marker_size=1):
    """Plot an interactive 3D UMAP colored by a continuous .obs (using plotly).

    Args:
        adata (anndata): The count data object
        identifiers (string / [string]): 
            One or several elements in chosen_column: Identifier (name, id, ...) of thing to color UMAP on.
            Will make one umap for each identifier.
        column (string): Column in adata.var to find identifier in. Defaults to adata.var.index
        palette (string / [string]):
            Either a name of a matplotlib predefined palette to pick from
            :class:`~matplotlib.colors.ListedColormap`,
            or a list of matplotlib predefined named colors.
            (see :func:`~matplotlib.colors.is_color_like`).
            The colors in the list will be interpolated to get a continuous color scale

    Returns:
        Nothing; displays interactive HTML
    """

    import matplotlib
    import matplotlib.cm as cm
    
    # adata.var['index'] doesnt work
    if column == 'index':
        chosen_column = adata.var.index
    else:
        chosen_column = adata.var[column]
        
    # Make identifiers a list if not already
    if isinstance(identifiers, str):
        identifiers = [identifiers]
        
    # Handle options for palettes
    if isinstance(palette, str):
        pal = cm.get_cmap(palette)
        # Reversed() is to match the scanpy default using the viridis palette
        pal = list(reversed([matplotlib.colors.to_hex(color) for color in pal.colors]))
    else:
        pal = [matplotlib.colors.cnames[color] for color in palette]
        
    for identifier in identifiers:
        identifier_index = np.where(chosen_column == identifier)[0]
        if len(identifier_index) != 1:
            print('Warning: ' + identifier + ' appears several times in ' + 
                  column + '\nMaking one UMAP for each occurrence..')
            for ii in identifier_index:
                bb = adata.X[:, identifier_index[ii]]
                __plot_3d_umap_continuous_inner(adata, color_values=bb, title=identifier, 
                                              palette=pal, marker_size=marker_size)
        else:
            bb = adata.X[:, identifier_index[0]]
            __plot_3d_umap_continuous_inner(adata, color_values=bb, title=identifier, palette=pal, marker_size=marker_size)


# (Inner) Function for 3D UMAP. Takes array (or series or list) of values for color values
def __plot_3d_umap_continuous_inner(adata, color_values, title, palette='viridis', marker_size=1):
    # adata:          anndata object
    # color_values:   Values used for color of markers. Same length as the number of observations in adata (adata.shape[0])
    # palette:        List of strings. Colors in hexadecimal. 
    #                 Colors in list will be interpolated to get a continuous color scale

    import plotly.express as px

    # Format data to pandas data frame
    x = adata.obsm["X_umap"][:,0].tolist()
    y = adata.obsm["X_umap"][:,1].tolist()
    z = adata.obsm["X_umap"][:,2].tolist()
    cols = {'UMAP1': x, 'UMAP2': y, 'UMAP3': z, title: color_values}
    df = pd.DataFrame(data=cols)

    # Currently going without titles since this leaves more space for data,
    # and since the name of the identifier (same as title currently) is displayed on legend anyways.
    # If titles wanted, enable title arg. Also disable margin=dict(l=0, r=0, b=0, t=0) below,
    # otherwise the title will be outside canvas.
    fig = px.scatter_3d(df, x='UMAP1', y='UMAP2', z='UMAP3', color=title, color_continuous_scale=palette,
#                         title=title
                       )
    fig.update_traces(marker=dict(size=marker_size))
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0), legend= {'itemsizing': 'constant'})
    fig.show()




# # Test continuous 3D umap

# print('Test: One unique identifier..')
# plot_3d_umap_continuous(adata, identifiers = 'ENSG00000187608_6_e')
# print('Test: One identifier from custom column that has duplicates..')
# plot_3d_umap_continuous(adata, identifiers = 'ISG15', column='gene_names')
# print('Test: List of two unique identifiers')
# plot_3d_umap_continuous(adata, identifiers = ['ENSG00000198692_3_e', 'ENSG00000187608_6_e'])
# print('Test: Custom colors')
# plot_3d_umap_continuous(adata, 'ENSG00000111732_1_e', palette=['cornflowerblue', 'lemonchiffon', 'brown'])
# plot_3d_umap_continuous(adata, 'ENSG00000111732_1_e', palette=['cornflowerblue', 'brown', 'lemonchiffon'])

# # All pass 2020-11-24 :)
# # Commenting out to make notebook less heavy








## Function for categorical 3D UMAPs
def plot_3d_umap_categorical(adata, column, palette=None, marker_size=1):
    """Plot an interactive 3D UMAP colored by a discrete/categorial .obs (using plotly).

    Args:
        adata (anndata): The count data object
        identifiers (string / [string]): 
            One or several elements in chosen_column: Identifier (name, id, ...) of thing to color UMAP on.
            Will make one umap for each identifier. (TODO support multiple)
        column (string): Column in adata.var to find identifier in. Defaults to adata.var.index
        palette (string / [string]):
            Either a name of a matplotlib predefined palette to pick from,
            or a list of colors with length matching the number of colors needed.
            Default: None --> matplotlib.rcParams["axes.prop_cycle"] used
            This is the default used by scanpy. See docs for scanpy.pl.umap.
            Not sure I like it, because scanpy also overwrites axes.prop_cycle when user plots
            with a custom palette, which has been more of an annoying side effect than a
            convenience for me so far in data analysis. However, it is easily worked around
            by the user, so going for compatibility/consistency here
            The colors in the list will be interpolated to get a continuous color scale

        marker_size (number): 
            Marker size in scatter.

    Returns:
        Nothing; displays interactive HTML
    """
    
    
    if palette is None:
        ppp = [d['color'] for d in matplotlib.rcParams["axes.prop_cycle"]]
    elif not isinstance(palette, list):
        ppp = matplotlib.cm.get_cmap(palette).colors
    else:
        ppp = palette
    
    pp = [matplotlib.colors.rgb2hex(color) for color in ppp]
    
    cats = adata.obs[column].cat.categories
    color_dict = {cats[ii]: pp[ii] for ii in range(len(cats))}
    
    x=adata.obsm["X_umap"][:,0].tolist()
    y=adata.obsm["X_umap"][:,1].tolist()
    z=adata.obsm["X_umap"][:,2].tolist()
    cats = adata.obs[column].tolist()
    cols = {'x': x, 'y': y, 'z': z, column: cats}
    df = pd.DataFrame(data=cols)
    
    fig = px.scatter_3d(df, x='x', y='y', z='z', color=column, color_discrete_map=color_dict)
    fig.update_traces(marker=dict(size=marker_size))
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0), legend= {'itemsizing': 'constant'})
    fig.show()


# # Test categorical 3D umap
# plot_3d_umap_categorical(adata, "gc_zone")
# plot_3d_umap_categorical(adata, "leiden", "tab10")
# plot_3d_umap_categorical(adata, "gc_zone", palette=['cornflowerblue', 'lemonchiffon', 'brown'])












###################################################################
###################################################################
###################################################################

# Helper functions

###################################################################
###################################################################
###################################################################




##### Currently not used (commented out), but might be updated
##### to be used by build_pretty_DE_table()
##### after response from uniProt on how to best
##### rank relevance of results from translating
##### ensembl gene IDs to uniProt protein IDs

# def get_uniProt_id_from_ensemble_id(ensemble_id):
#     # Translates Ensemble gene IDs to uniProt IDs
#     # Currently only takes a single gene ID,
#     # and picks the first in the list of returned 
#     # uniProt IDs. Both not optimal.
#     # Waiting for uniProt email answer on how to best
#     # rank the uniProt IDs on relevance before fixing
#     # anything further.

#     import urllib.parse
#     import urllib.request

#     url = 'https://www.uniprot.org/uploadlists/'
#     params = {
#     'from': 'ENSEMBL_ID',
#     'to': 'ACC',
#     'format': 'tab',
#     'query': ensemble_id
#     }

#     data = urllib.parse.urlencode(params)
#     data = data.encode('utf-8')
#     req = urllib.request.Request(url, data)
#     with urllib.request.urlopen(req) as f:
#        response = f.read()
# #     print(response.decode('utf-8'))
#     return(response.decode('utf-8').split('\n')[1].split('\t')[1])

# uniProt_id = get_uniProt_id_from_ensemble_id('ENSG00000187608')
# print(uniProt_id)








