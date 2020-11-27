from IPython.display import HTML, display
from random import randint
import plotly.express as px
import matplotlib
import matplotlib.cm





def pretty_rank_gene_groups(adata, maxgenes=20):
    maxgenes=20

    ### How many comparisons and genes?
    numgrp=len(adata.uns['rank_genes_groups']["names"][0])
    numgene=len(adata.uns['rank_genes_groups']["names"])
    numgene=min(maxgenes,numgene)


    ## Format one gene as a table row
    def format_onegene(grp, i):
        gname=(adata.uns['rank_genes_groups']["names"][i])[grp]
        score=(adata.uns['rank_genes_groups']["scores"][i])[grp]
        lfc = (adata.uns['rank_genes_groups']["logfoldchanges"][i])[grp]
        pval = (adata.uns['rank_genes_groups']["pvals"][i])[grp]
        link_genecards = "https://www.genecards.org/cgi-bin/carddisp.pl?gene="+gname+"&keywords="+gname
        link_uniprot = "https://www.uniprot.org/uniprot/?query="+gname+"&sort=score"
        link_wp = "https://en.wikipedia.org/w/index.php?search="+gname+"&go=Go"

        #if ENSG/ENSMUS given, can go straight instead. 
        link_ensembl = "http://www.ensembl.org/Multi/Search/Results?q="+gname+";site=ensembl_all"

        #If species known, can make ucsc a bit clever too
        #link_ucsc = "http://www.ensembl.org/Multi/Search/Results?q="+gname+";site=ensembl_all"

        return('<td>'+gname+'</td>'+
               '<td>'+str(score)+'</td>'+
               '<td>'+str(pval)+'</td>'+
               '<td>'+str(lfc)+'</td>'+
               '<td><a target="_blank" href="'+link_genecards+'">Genecards</a></td>' +
               '<td><a target="_blank" href="'+link_wp+'">Wikipedia</a></td>' +
               '<td><a target="_blank" href="'+link_ensembl+'">Ensembl</a></td>' +
               '<td><a target="_blank" href="'+link_uniprot+'">Uniprot</a></td>' 
              )

    ## Format the table for one gene
    def format_onetable(grp):
        tab_format = """<table><tr>
            <th>Gene</th>
            <th>Score</th>
            <th>p-value</th>
            <th>LogFC</th>
            <th></th>
            <th></th>
            <th></th>
            <th></th>
            </tr><tr>{}</tr></table>
            """
        return(tab_format.format('</tr><tr>'.join([format_onegene(grp,i) for i in range(0,numgene)])))

    ## The buttons above the panel
    all_tab_buttons=''.join(
        ["<button class=\"INTANCEID_tablinks\" onclick=\"INTANCEID_open_genepanel(event, 'INTANCEID_"+
         str(grp)+"')\">"+str(grp)+"</button>" for grp in range(0,numgrp)])
    all_tab_buttons = '<div class="INTANCEID_tab">{}</div>'.format(all_tab_buttons)

    ## The panels themselves
    all_panels=''.join(
        ['<div id="INTANCEID_'+str(grp)+'" class="INTANCEID_tabcontent">{}</div>'.format(format_onetable(grp)) 
         for grp in range(0,numgrp)])

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
        function INTANCEID_open_genepanel(evt, cityName) {
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
          document.getElementById(cityName).style.display = "block";
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

    display(HTML(
       js.replace("INTANCEID_","foo"+str(randint(0, 100000))+"_")
    ))






###################################################################
###################################################################
###################################################################




def plot_umaps_sidebyside(adata1, adata2, obsname="leiden", palette="Set3"):

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
