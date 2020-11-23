from IPython.display import HTML, display
from random import randint


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


